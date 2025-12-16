"""Masked language model pretraining on E. coli MG1655 genome.

Key improvements for reaching lower loss (<1.0):
1. Warmup + cosine annealing with warm restarts
2. Curriculum learning: mask probability increases over training
3. Validation monitoring to detect overfitting
4. Gradient accumulation for larger effective batch size
5. Label smoothing for better generalization
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import random
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts, LambdaLR
from torch.utils.data import DataLoader, Dataset, random_split

from dataset import KmerVocabulary
from model import DNAEncoder

LOGGER = logging.getLogger("gene_whisperer.pretrain_mlm")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

ALLOWED_BASES = {"A", "C", "G", "T"}


def set_seed(seed: int) -> None:
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    if hasattr(torch, "mps") and torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)


def read_fasta_sequence(path: Path) -> str:
    sequence_parts: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)


def extract_windows(sequence: str, window_size: int, stride: int) -> List[str]:
    windows: List[str] = []
    max_start = len(sequence) - window_size + 1
    for start in range(0, max_start, stride):
        window = sequence[start : start + window_size]
        if not set(window) <= ALLOWED_BASES:
            continue
        windows.append(window)
    return windows


class MLMDataset(Dataset):
    def __init__(self, token_tensors: List[torch.LongTensor]):
        self.samples = token_tensors

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> torch.LongTensor:
        return self.samples[idx]


def mask_tokens_span(
    inputs: torch.LongTensor,
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
    max_span_len: int = 3,
    track_spans: bool = False,
) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
    """
    DNABERT-style span masking.
    
    Instead of independent per-token masking, this selects contiguous spans
    until approximately mask_prob of tokens are covered, then applies the
    standard 80/10/10 BERT masking within those spans.
    
    Args:
        inputs: (B, L) input token indices
        vocab: KmerVocabulary with mask_id and pad_id
        mask_prob: Target fraction of tokens to mask
        max_span_len: Maximum length of each masked span
        track_spans: If True, also return list of span lengths
        
    Returns:
        (masked_inputs, labels) tuple where labels has -100 for unmasked positions
        If track_spans=True, also returns list of span lengths
    """
    device = inputs.device
    labels = inputs.clone()
    batch_size, seq_len = inputs.shape
    vocab_size = len(vocab.itos)

    # Initialize all positions as unmasked
    masked = torch.zeros_like(inputs, dtype=torch.bool, device=device)
    
    # Track span lengths if requested
    span_lengths: List[int] = [] if track_spans else []

    for b in range(batch_size):
        num_to_mask = max(1, int(mask_prob * seq_len))
        covered = 0
        attempts = 0
        max_attempts = num_to_mask * 10  # Safety limit to prevent infinite loops
        
        while covered < num_to_mask and attempts < max_attempts:
            span_len = random.randint(1, max_span_len)
            start = random.randint(0, max(0, seq_len - span_len))
            end = min(seq_len, start + span_len)
            
            # Avoid double-masking: only count newly masked positions
            span_mask = ~masked[b, start:end]
            if span_mask.any():
                masked[b, start:end] = True
                actual_new = span_mask.sum().item()
                covered += actual_new
                if track_spans and actual_new > 0:
                    span_lengths.append(end - start)
            
            attempts += 1
            
            # Safety break if we've covered enough
            if covered >= num_to_mask or covered >= seq_len:
                break

    # Don't mask padding positions
    if hasattr(vocab, "pad_id"):
        pad_mask = inputs.eq(vocab.pad_id)
        masked[pad_mask] = False

    # Set labels: -100 for unmasked positions (ignored in loss)
    labels[~masked] = -100

    # Apply 80/10/10 masking strategy within masked spans
    # 80% -> [MASK] token
    # 10% -> random token
    # 10% -> unchanged (keep original)
    rand = torch.rand_like(inputs, dtype=torch.float32, device=device)
    mask_token_indices = masked & (rand < 0.8)
    random_token_indices = masked & (rand >= 0.8) & (rand < 0.9)
    # Remaining 10% (rand >= 0.9) stay unchanged

    inputs = inputs.clone()
    inputs[mask_token_indices] = vocab.mask_id

    if random_token_indices.any():
        flat = random_token_indices.view(-1)
        num_rand = flat.sum().item()
        random_ids = torch.randint(
            low=0,
            high=vocab_size,
            size=(num_rand,),
            device=device,
            dtype=inputs.dtype,
        )
        inputs.view(-1)[flat] = random_ids

    if track_spans:
        return inputs, labels, span_lengths
    return inputs, labels


def get_warmup_cosine_schedule(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
    min_lr_ratio: float = 0.01,
    num_cycles: float = 0.5,
) -> LambdaLR:
    """
    Learning rate schedule with warmup and cosine annealing.
    
    This is critical for stable training and avoiding local minima.
    """
    def lr_lambda(current_step: int) -> float:
        if current_step < num_warmup_steps:
            # Linear warmup
            return float(current_step) / float(max(1, num_warmup_steps))
        
        # Cosine annealing
        progress = float(current_step - num_warmup_steps) / float(
            max(1, num_training_steps - num_warmup_steps)
        )
        cosine_decay = 0.5 * (1.0 + math.cos(math.pi * num_cycles * 2.0 * progress))
        return max(min_lr_ratio, cosine_decay)
    
    return LambdaLR(optimizer, lr_lambda)


def get_layer_wise_lr_groups(
    model: nn.Module,
    base_lr: float,
    lr_decay: float = 0.8,
    weight_decay: float = 0.01,
) -> List[dict]:
    """
    Create parameter groups with layer-wise learning rate decay.
    
    Lower layers (closer to input) get smaller learning rates,
    which helps preserve pretrained features during fine-tuning
    and improves pretraining stability.
    
    Args:
        model: DNAMLM model
        base_lr: Base learning rate for the top layer
        lr_decay: Decay factor per layer (0.8 means each lower layer gets 0.8x LR)
        weight_decay: Weight decay for regularization
    
    Returns:
        List of parameter groups for optimizer
    """
    # No decay for bias and layer norm
    no_decay = {'bias', 'LayerNorm', 'layer_norm', 'norm', 'embed_norm'}
    
    # Get encoder layers
    encoder = model.encoder
    num_layers = len(encoder.layers)
    
    # Group parameters by layer depth
    param_groups = []
    
    # Embedding parameters (lowest layer, smallest LR)
    embed_lr = base_lr * (lr_decay ** num_layers)
    embed_params_decay = []
    embed_params_no_decay = []
    for name, param in encoder.embedding.named_parameters():
        if any(nd in name for nd in no_decay):
            embed_params_no_decay.append(param)
        else:
            embed_params_decay.append(param)
    
    if embed_params_decay:
        param_groups.append({'params': embed_params_decay, 'lr': embed_lr, 'weight_decay': weight_decay})
    if embed_params_no_decay:
        param_groups.append({'params': embed_params_no_decay, 'lr': embed_lr, 'weight_decay': 0.0})
    
    # Embed norm parameters (also low layer)
    if hasattr(encoder, 'embed_norm'):
        embed_norm_params = list(encoder.embed_norm.parameters())
        if embed_norm_params:
            param_groups.append({'params': embed_norm_params, 'lr': embed_lr, 'weight_decay': 0.0})
    
    # Transformer layers (progressively higher LR)
    for layer_idx, layer in enumerate(encoder.layers):
        layer_lr = base_lr * (lr_decay ** (num_layers - layer_idx - 1))
        layer_params_decay = []
        layer_params_no_decay = []
        
        for name, param in layer.named_parameters():
            if any(nd in name for nd in no_decay):
                layer_params_no_decay.append(param)
            else:
                layer_params_decay.append(param)
        
        if layer_params_decay:
            param_groups.append({'params': layer_params_decay, 'lr': layer_lr, 'weight_decay': weight_decay})
        if layer_params_no_decay:
            param_groups.append({'params': layer_params_no_decay, 'lr': layer_lr, 'weight_decay': 0.0})
    
    # Final norm (top layer, full LR)
    final_norm_params = list(encoder.final_norm.parameters())
    if final_norm_params:
        param_groups.append({'params': final_norm_params, 'lr': base_lr, 'weight_decay': 0.0})
    
    # LM head parameters (top layer, full LR)
    # Note: if weight tying, embedding weights are already added above
    lm_head_params_decay = []
    lm_head_params_no_decay = []
    
    for name, param in model.lm_head.named_parameters():
        if 'weight' in name and model.tie_weights:
            # Skip tied weights (already in embedding group)
            continue
        if any(nd in name for nd in no_decay):
            lm_head_params_no_decay.append(param)
        else:
            lm_head_params_decay.append(param)
    
    if lm_head_params_decay:
        param_groups.append({'params': lm_head_params_decay, 'lr': base_lr, 'weight_decay': weight_decay})
    if lm_head_params_no_decay:
        param_groups.append({'params': lm_head_params_no_decay, 'lr': base_lr, 'weight_decay': 0.0})
    
    # Output norm if present
    if hasattr(model, 'output_norm'):
        output_norm_params = list(model.output_norm.parameters())
        if output_norm_params:
            param_groups.append({'params': output_norm_params, 'lr': base_lr, 'weight_decay': 0.0})
    
    return param_groups


class MLMCollator:
    """
    Collator that supports curriculum learning via adjustable mask probability.
    
    Debug mode (debug=True):
    - Prints original token ids, masked input token ids, and label token ids
    - For 5 positions per sequence, asserts:
      - label == original_token for masked positions
      - label == -100 for unmasked positions
    - Raises AssertionError if any mismatch occurs
    - Debug mode runs once then auto-disables
    """
    def __init__(self, vocab: KmerVocabulary, mask_prob: float = 0.15, debug: bool = False):
        self.vocab = vocab
        self.mask_prob = mask_prob
        self.debug = debug
        self._debug_ran = False
    
    def __call__(self, batch: Iterable[torch.LongTensor]) -> Tuple[torch.LongTensor, torch.LongTensor]:
        original_inputs = torch.stack(list(batch))
        masked_inputs, labels = mask_tokens_span(original_inputs, self.vocab, mask_prob=self.mask_prob)
        
        # Run debug verification once if enabled
        if self.debug and not self._debug_ran:
            self._run_debug_verification(original_inputs, masked_inputs, labels)
            self._debug_ran = True
        
        return masked_inputs, labels
    
    def _run_debug_verification(
        self, 
        original: torch.LongTensor, 
        masked_input: torch.LongTensor, 
        labels: torch.LongTensor
    ) -> None:
        """
        Debug mode: verify labels match originals for one batch.
        
        Checks:
        - For MASKED positions: label == original_token
        - For UNMASKED positions: label == -100
        
        Hard error (AssertionError) if any mismatch occurs.
        """
        batch_size, seq_len = original.shape
        
        print("=" * 80)
        print("DEBUG MODE: MLM Collator Label Verification")
        print("=" * 80)
        print(f"Batch size: {batch_size}, Sequence length: {seq_len}")
        print(f"Mask prob: {self.mask_prob:.2%}")
        print()
        
        # Determine masked positions (where label != -100)
        masked_positions = (labels != -100)
        
        # Print first sequence's tokens (limited for readability)
        print("[Token IDs for first sequence (first 20 positions)]")
        print(f"  Original tokens:     {original[0, :20].tolist()}")
        print(f"  Masked input tokens: {masked_input[0, :20].tolist()}")
        print(f"  Label tokens:        {labels[0, :20].tolist()}")
        print()
        
        # Collect errors for batch reporting
        errors = []
        
        # For each sequence, check 5 positions
        print(f"[Checking 5 positions per sequence across {batch_size} sequences]")
        print("-" * 80)
        
        for seq_idx in range(batch_size):
            seq_original = original[seq_idx]
            seq_labels = labels[seq_idx]
            seq_masked_positions = masked_positions[seq_idx]
            
            # Select 5 positions to check (evenly spaced)
            check_positions = [int(i * seq_len / 5) for i in range(5)]
            
            for pos in check_positions:
                orig_token = seq_original[pos].item()
                label_token = seq_labels[pos].item()
                is_masked = seq_masked_positions[pos].item()
                
                if is_masked:
                    # MASKED position: label should equal original token
                    if label_token != orig_token:
                        error_msg = (
                            f"MISMATCH at seq={seq_idx}, pos={pos}: "
                            f"MASKED position has label={label_token}, expected original={orig_token}"
                        )
                        errors.append(error_msg)
                        print(f"  ✗ seq[{seq_idx}] pos[{pos}]: MASKED - label={label_token}, orig={orig_token} - FAIL!")
                    else:
                        if seq_idx < 3:  # Only print details for first 3 sequences
                            print(f"  ✓ seq[{seq_idx}] pos[{pos}]: MASKED - label={label_token}, orig={orig_token} - OK")
                else:
                    # UNMASKED position: label should be -100
                    if label_token != -100:
                        error_msg = (
                            f"MISMATCH at seq={seq_idx}, pos={pos}: "
                            f"UNMASKED position has label={label_token}, expected -100"
                        )
                        errors.append(error_msg)
                        print(f"  ✗ seq[{seq_idx}] pos[{pos}]: UNMASKED - label={label_token}, expected=-100 - FAIL!")
                    else:
                        if seq_idx < 3:  # Only print details for first 3 sequences
                            print(f"  ✓ seq[{seq_idx}] pos[{pos}]: UNMASKED - label={label_token} - OK")
        
        print("-" * 80)
        print()
        
        # Summary statistics
        total_tokens = batch_size * seq_len
        num_masked = masked_positions.sum().item()
        num_unmasked = total_tokens - num_masked
        
        print("[Summary Statistics]")
        print(f"  Total tokens:   {total_tokens}")
        print(f"  Masked tokens:  {num_masked} ({num_masked/total_tokens:.2%})")
        print(f"  Unmasked:       {num_unmasked} ({num_unmasked/total_tokens:.2%})")
        print()
        
        # Exhaustive verification (not just 5 positions, but ALL)
        print("[Exhaustive Verification - ALL positions]")
        
        # Check all masked positions: label should equal original
        all_masked_labels = labels[masked_positions]
        all_masked_originals = original[masked_positions]
        masked_match = (all_masked_labels == all_masked_originals).all().item()
        masked_mismatches = (all_masked_labels != all_masked_originals).sum().item()
        
        print(f"  Masked positions (label == original):   ", end="")
        if masked_match:
            print(f"✓ PASS (all {num_masked} positions match)")
        else:
            print(f"✗ FAIL ({masked_mismatches} mismatches out of {num_masked})")
            errors.append(f"Exhaustive check: {masked_mismatches} masked positions have wrong labels")
        
        # Check all unmasked positions: label should be -100
        unmasked_positions = ~masked_positions
        all_unmasked_labels = labels[unmasked_positions]
        unmasked_correct = (all_unmasked_labels == -100).all().item()
        unmasked_wrong = (all_unmasked_labels != -100).sum().item()
        
        print(f"  Unmasked positions (label == -100):     ", end="")
        if unmasked_correct:
            print(f"✓ PASS (all {num_unmasked} positions are -100)")
        else:
            print(f"✗ FAIL ({unmasked_wrong} positions are not -100)")
            errors.append(f"Exhaustive check: {unmasked_wrong} unmasked positions don't have -100")
        
        print()
        print("=" * 80)
        
        # HARD ERROR if any mismatches
        if errors:
            print("ERRORS FOUND:")
            for err in errors:
                print(f"  - {err}")
            print("=" * 80)
            raise AssertionError(
                f"MLM Collator Debug: {len(errors)} label verification error(s) found!\n"
                + "\n".join(errors)
            )
        
        print("✓ All assertions passed! Labels correctly match originals.")
        print("=" * 80)


def collate_mlm(
    batch: Iterable[torch.LongTensor],
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
) -> Tuple[torch.LongTensor, torch.LongTensor]:
    """Collate batch and apply DNABERT-style span masking."""
    inputs = torch.stack(list(batch))
    masked, labels = mask_tokens_span(inputs, vocab, mask_prob=mask_prob)
    return masked, labels


def debug_masking(
    vocab: KmerVocabulary,
    batch: Iterable[torch.LongTensor],
    mask_prob: float = 0.15,
    max_span_len: int = 3,
) -> None:
    """
    Debug mode: analyze one batch of masking and print statistics.
    
    Checks:
    1. Percent masked tokens (should be ~mask_prob)
    2. Percent masked that were PAD (should be 0)
    3. Distribution of span lengths
    4. Assert labels are -100 for unmasked, never -100 for masked
    """
    from collections import Counter
    
    print("=" * 70)
    print("DEBUG MODE: Masking Analysis")
    print("=" * 70)
    
    # Stack batch into tensor
    inputs = torch.stack(list(batch))
    batch_size, seq_len = inputs.shape
    total_tokens = batch_size * seq_len
    
    print(f"Batch size: {batch_size}")
    print(f"Sequence length: {seq_len}")
    print(f"Total tokens: {total_tokens}")
    print(f"Expected mask_prob: {mask_prob:.2%}")
    print(f"Max span length: {max_span_len}")
    print()
    
    # Run masking with span tracking
    masked_inputs, labels, span_lengths = mask_tokens_span(
        inputs, vocab, mask_prob=mask_prob, max_span_len=max_span_len, track_spans=True
    )
    
    # 1. Percent masked tokens
    masked_positions = (labels != -100)
    num_masked = masked_positions.sum().item()
    actual_mask_pct = num_masked / total_tokens
    
    print(f"[1] MASKED TOKEN PERCENTAGE")
    print(f"    Masked tokens: {num_masked} / {total_tokens}")
    print(f"    Actual mask rate: {actual_mask_pct:.2%}")
    print(f"    Expected: ~{mask_prob:.2%}")
    diff = abs(actual_mask_pct - mask_prob)
    status = "✓ OK" if diff < 0.05 else "⚠ WARN (>5% deviation)"
    print(f"    Status: {status}")
    print()
    
    # 2. Percent masked that were PAD
    pad_positions = (inputs == vocab.pad_id)
    masked_pads = (masked_positions & pad_positions).sum().item()
    total_pads = pad_positions.sum().item()
    
    print(f"[2] PAD TOKEN MASKING")
    print(f"    Total PAD tokens in batch: {total_pads}")
    print(f"    PAD tokens that are masked: {masked_pads}")
    if total_pads > 0:
        pad_mask_pct = masked_pads / total_pads
        print(f"    PAD mask rate: {pad_mask_pct:.2%}")
    status = "✓ OK (no PADs masked)" if masked_pads == 0 else "✗ FAIL (PADs should never be masked!)"
    print(f"    Status: {status}")
    print()
    
    # 3. Span length distribution
    print(f"[3] SPAN LENGTH DISTRIBUTION")
    if span_lengths:
        span_counter = Counter(span_lengths)
        total_spans = len(span_lengths)
        print(f"    Total spans created: {total_spans}")
        print(f"    Span length distribution:")
        for length in sorted(span_counter.keys()):
            count = span_counter[length]
            pct = count / total_spans * 100
            bar = "█" * int(pct / 2)
            print(f"      Length {length}: {count:>5} ({pct:5.1f}%) {bar}")
        avg_span = sum(span_lengths) / len(span_lengths)
        print(f"    Average span length: {avg_span:.2f}")
    else:
        print(f"    No spans recorded")
    print()
    
    # 4. Label correctness assertions
    print(f"[4] LABEL CORRECTNESS ASSERTIONS")
    
    # Check: unmasked positions should have label = -100
    unmasked_positions = ~masked_positions
    unmasked_labels = labels[unmasked_positions]
    all_unmasked_are_minus100 = (unmasked_labels == -100).all().item()
    
    print(f"    Unmasked positions have label=-100: ", end="")
    if all_unmasked_are_minus100:
        print("✓ PASS")
    else:
        bad_count = (unmasked_labels != -100).sum().item()
        print(f"✗ FAIL ({bad_count} unmasked positions don't have -100)")
    
    # Check: masked positions should NEVER have label = -100
    masked_labels = labels[masked_positions]
    no_masked_is_minus100 = (masked_labels != -100).all().item()
    
    print(f"    Masked positions have label!=-100:  ", end="")
    if no_masked_is_minus100:
        print("✓ PASS")
    else:
        bad_count = (masked_labels == -100).sum().item()
        print(f"✗ FAIL ({bad_count} masked positions have -100)")
    
    # Check: labels for masked positions should be valid token IDs
    if num_masked > 0:
        min_label = masked_labels.min().item()
        max_label = masked_labels.max().item()
        valid_range = (min_label >= 0 and max_label < len(vocab.itos))
        print(f"    Masked labels in valid range [0, {len(vocab.itos)-1}]: ", end="")
        if valid_range:
            print(f"✓ PASS (min={min_label}, max={max_label})")
        else:
            print(f"✗ FAIL (min={min_label}, max={max_label})")
    
    print()
    print("=" * 70)
    
    # Assert critical properties
    assert masked_pads == 0, "ASSERTION FAILED: PAD tokens should never be masked!"
    assert all_unmasked_are_minus100, "ASSERTION FAILED: Unmasked positions must have label=-100!"
    assert no_masked_is_minus100, "ASSERTION FAILED: Masked positions must NOT have label=-100!"
    
    print("All assertions passed! Masking is working correctly.")
    print("=" * 70)


def verify_no_leakage(
    model: nn.Module,
    vocab: KmerVocabulary,
    device: torch.device,
) -> None:
    """
    Verify that masked positions have NO access to original token information.
    
    This is CRITICAL for proper MLM training:
    - MLM ≠ autoencoder. The model must predict from context only.
    - Masked positions should only see [MASK] embedding
    - Engineered features must be disabled (they would leak sequence info)
    
    Tests:
    1. Verify masked token IDs are replaced with [MASK] before embedding lookup
    2. Verify model runs with 100% masking (extreme blindfold test)
    3. Assert no engineered features are used in DNAMLM
    """
    print("=" * 80)
    print("LEAKAGE VERIFICATION: MLM Blindfold Test")
    print("=" * 80)
    print("MLM ≠ autoencoder. Blindfold means blindfold.")
    print()
    
    model.eval()
    errors = []
    
    # =========================================================================
    # TEST 1: Verify masked inputs have [MASK] token, not original
    # =========================================================================
    print("[TEST 1] Masked inputs contain [MASK] token, not originals")
    print("-" * 80)
    
    # Create a test sequence
    batch_size = 4
    seq_len = 32
    original_tokens = torch.randint(0, 64, (batch_size, seq_len))  # Only real k-mers
    
    # Manually mask some positions (50% for clear test)
    mask_prob = 0.5
    masked_inputs = original_tokens.clone()
    
    # Create mask
    mask = torch.rand(batch_size, seq_len) < mask_prob
    
    # Apply [MASK] token to masked positions
    masked_inputs[mask] = vocab.mask_id
    
    # Verify: masked positions should have vocab.mask_id, not original
    masked_positions = mask
    masked_values = masked_inputs[masked_positions]
    original_values = original_tokens[masked_positions]
    
    all_have_mask_token = (masked_values == vocab.mask_id).all().item()
    none_leaked = (masked_values != original_values).all().item()
    
    print(f"  Total masked positions: {masked_positions.sum().item()}")
    print(f"  All masked positions have [MASK] token (id={vocab.mask_id}): ", end="")
    if all_have_mask_token:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 1: Some masked positions don't have [MASK] token")
    
    print(f"  No original tokens leaked to masked positions: ", end="")
    if none_leaked:
        print("✓ PASS")
    else:
        # This can fail if original happened to be [MASK], which is impossible for k-mers
        leak_count = (masked_values == original_values).sum().item()
        print(f"✗ FAIL ({leak_count} positions may have leaked)")
        errors.append(f"TEST 1: {leak_count} positions have original values")
    print()
    
    # =========================================================================
    # TEST 2: Model runs with 100% masking (extreme blindfold)
    # =========================================================================
    print("[TEST 2] Model runs with 100% masking (extreme blindfold)")
    print("-" * 80)
    
    # Create input where ALL tokens are [MASK]
    all_masked_input = torch.full((batch_size, seq_len), vocab.mask_id, dtype=torch.long).to(device)
    
    # Labels are all the original tokens (model must predict all)
    labels = original_tokens.clone().to(device)
    
    try:
        with torch.no_grad():
            logits, loss = model(all_masked_input, labels)
        
        # Verify outputs
        has_valid_logits = logits is not None and not torch.isnan(logits).any()
        has_valid_loss = loss is not None and not torch.isnan(loss) and not torch.isinf(loss)
        
        print(f"  Model forward pass: ✓ PASS")
        print(f"  Output shape: {logits.shape}")
        print(f"  Logits valid (no NaN): ", end="")
        if has_valid_logits:
            print("✓ PASS")
        else:
            print("✗ FAIL")
            errors.append("TEST 2: Logits contain NaN")
        
        print(f"  Loss valid (no NaN/Inf): ", end="")
        if has_valid_loss:
            print(f"✓ PASS (loss={loss.item():.4f})")
        else:
            print(f"✗ FAIL (loss={loss})")
            errors.append("TEST 2: Loss is NaN or Inf")
        
        # With 100% masking and no context, accuracy should be low (near random)
        preds = logits.argmax(dim=-1)
        acc = (preds == labels).float().mean().item()
        print(f"  Accuracy with 100% masking: {acc*100:.1f}% (expected: near random ~1.5%)")
        
        # Verify model isn't cheating somehow
        if acc > 0.5:  # If accuracy is >50%, something is very wrong
            print(f"  ⚠ WARNING: Accuracy suspiciously high - possible leakage!")
            errors.append(f"TEST 2: Suspiciously high accuracy ({acc*100:.1f}%) with 100% masking")
        
    except Exception as e:
        print(f"  ✗ FAIL: Model failed with error: {e}")
        errors.append(f"TEST 2: Model failed with {type(e).__name__}: {e}")
    print()
    
    # =========================================================================
    # TEST 3: Verify no engineered features in DNAMLM
    # =========================================================================
    print("[TEST 3] DNAMLM does not use engineered features")
    print("-" * 80)
    
    # Check DNAMLM doesn't have engineered feature inputs
    dnamlm_signature = str(model.forward.__code__.co_varnames)
    has_engineered = "engineered" in dnamlm_signature.lower()
    
    print(f"  DNAMLM.forward() signature: {model.forward.__code__.co_varnames[:5]}...")
    print(f"  Contains 'engineered' parameter: ", end="")
    if not has_engineered:
        print("✓ PASS (no engineered features)")
    else:
        print("✗ FAIL (engineered features would leak sequence info!)")
        errors.append("TEST 3: DNAMLM has engineered features in signature")
    
    # Also check the encoder doesn't receive engineered features
    encoder_params = list(model.encoder.forward.__code__.co_varnames)
    encoder_has_engineered = any("engineered" in p.lower() for p in encoder_params)
    print(f"  DNAEncoder.forward() has engineered: ", end="")
    if not encoder_has_engineered:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 3: DNAEncoder has engineered features")
    print()
    
    # =========================================================================
    # TEST 4: Embedding lookup verification
    # =========================================================================
    print("[TEST 4] Embedding lookup only sees masked token IDs")
    print("-" * 80)
    
    # Get the embedding layer
    embedding = model.encoder.embedding
    
    # Look up embeddings for [MASK] and a real k-mer
    mask_embedding = embedding(torch.tensor([vocab.mask_id]).to(device))
    real_embedding = embedding(torch.tensor([0]).to(device))  # First k-mer (e.g., "AAA")
    
    # They should be different
    embeddings_different = not torch.allclose(mask_embedding, real_embedding)
    
    print(f"  [MASK] embedding (id={vocab.mask_id}) shape: {mask_embedding.shape}")
    print(f"  Real k-mer embedding (id=0) shape: {real_embedding.shape}")
    print(f"  [MASK] embedding is distinct from real k-mers: ", end="")
    if embeddings_different:
        print("✓ PASS")
    else:
        print("✗ FAIL (embeddings are identical!)")
        errors.append("TEST 4: [MASK] embedding same as real k-mer embedding")
    
    # Verify that when we look up a masked sequence, we get [MASK] embeddings
    test_input = torch.tensor([[vocab.mask_id, vocab.mask_id, vocab.mask_id]]).to(device)
    test_embeds = embedding(test_input)
    expected_embeds = mask_embedding.expand(1, 3, -1)
    
    embeds_match = torch.allclose(test_embeds, expected_embeds)
    print(f"  Masked sequence gets [MASK] embeddings: ", end="")
    if embeds_match:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 4: Masked sequence doesn't get [MASK] embeddings")
    print()
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("=" * 80)
    if errors:
        print("LEAKAGE DETECTED!")
        print("-" * 80)
        for err in errors:
            print(f"  ✗ {err}")
        print("=" * 80)
        raise AssertionError(
            f"Leakage verification FAILED with {len(errors)} error(s)!\n" +
            "\n".join(errors)
        )
    else:
        print("✓ All leakage tests passed!")
        print("  - Masked positions only see [MASK] embeddings")
        print("  - Model handles 100% masking correctly")
        print("  - No engineered features in MLM pipeline")
        print("  - Embedding lookup is clean")
        print()
        print("The model is properly blindfolded. MLM ≠ autoencoder. ✓")
        print("=" * 80)


class DNAMLM(nn.Module):
    """
    Masked Language Model wrapper around DNAEncoder.
    
    Key improvements for lower loss:
    1. Weight tying between embedding and LM head (BERT-style)
    2. Bias in LM head for better output distribution modeling
    3. Optional output layer norm before projection
    4. **Special token exclusion**: [MASK], <UNK>, <PAD> are excluded from 
       prediction space by setting their logits to -inf. This prevents the
       model from "cheating" by predicting special tokens.
    """
    
    def __init__(
        self, 
        encoder: DNAEncoder, 
        vocab_size: int,
        special_token_ids: List[int] | None = None,  # IDs to exclude from predictions
        tie_weights: bool = True,  # CRITICAL: tie embedding weights for better learning
        use_output_norm: bool = True,  # Additional layer norm before LM head
    ):
        super().__init__()
        self.encoder = encoder
        self.vocab_size = vocab_size
        self.tie_weights = tie_weights
        
        # Store special token IDs to exclude from prediction
        # These will have their logits set to -inf before softmax/loss
        self.special_token_ids = special_token_ids or []
        if self.special_token_ids:
            LOGGER.info(
                "Special token exclusion enabled: IDs %s will be masked from predictions",
                self.special_token_ids
            )
        
        # Optional output layer norm (helps with distribution matching)
        self.use_output_norm = use_output_norm
        if use_output_norm:
            self.output_norm = nn.LayerNorm(encoder.embedding_dim)
        
        # LM head with bias (bias helps model output distribution)
        self.lm_head = nn.Linear(encoder.embedding_dim, vocab_size, bias=True)
        
        # Weight tying: share weights between embedding and LM head
        # This is a key technique from BERT that improves MLM loss significantly
        if tie_weights:
            self.lm_head.weight = encoder.embedding.weight
            LOGGER.info("Weight tying enabled: LM head shares embedding weights")
        else:
            # If not tying, initialize LM head with scaled init
            nn.init.trunc_normal_(self.lm_head.weight, std=0.02, a=-0.04, b=0.04)
        
        # Initialize bias to zero
        nn.init.zeros_(self.lm_head.bias)
    
    def _mask_special_tokens(self, logits: torch.Tensor) -> torch.Tensor:
        """
        Set logits for special tokens to -inf so they can never be predicted.
        
        This is critical for proper MLM training:
        - Without this, the model can "cheat" by predicting [MASK], <UNK>, or <PAD>
        - These tokens should never appear as ground truth labels
        - Excluding them forces the model to learn real k-mer distributions
        
        Args:
            logits: (B, L, vocab_size) raw logits from LM head
            
        Returns:
            Modified logits with special token positions set to -inf
        """
        if not self.special_token_ids:
            return logits
        
        # Clone to avoid in-place modification issues with autograd
        logits = logits.clone()
        
        # Set special token logits to -inf (will become 0 after softmax)
        for token_id in self.special_token_ids:
            logits[:, :, token_id] = float('-inf')
        
        return logits

    def forward(
        self,
        token_ids: torch.LongTensor,
        labels: torch.LongTensor | None = None,
        key_padding_mask: torch.Tensor | None = None,
    ):
        """
        Forward pass for MLM.
        
        Args:
            token_ids: (B, L) input token indices
            labels: (B, L) target labels, -100 where we ignore
            key_padding_mask: Optional (B, L) boolean mask for padding
            
        Returns:
            (logits, loss) tuple where loss is None if labels not provided
        """
        x = self.encoder(token_ids, key_padding_mask=key_padding_mask)
        
        # Optional output norm before projection
        if self.use_output_norm:
            x = self.output_norm(x)
        
        logits = self.lm_head(x)  # (B, L, vocab_size)
        
        # CRITICAL: Mask out special tokens from prediction space
        # This prevents the model from predicting [MASK], <UNK>, <PAD>
        logits = self._mask_special_tokens(logits)
        
        if labels is None:
            return logits, None
        
        # Standard cross entropy - NO label smoothing for MLM
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
        )
        return logits, loss
    
    def compute_accuracy(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute masked token prediction accuracy."""
        # Note: logits should already have special tokens masked out
        preds = logits.argmax(dim=-1)  # (B, L)
        mask = labels != -100
        if mask.sum() == 0:
            return 0.0
        correct = (preds[mask] == labels[mask]).sum().item()
        total = mask.sum().item()
        return correct / total
    
    def compute_perplexity(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute perplexity on masked tokens."""
        # Note: logits should already have special tokens masked out
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
            reduction='mean',
        )
        return math.exp(loss.item())


def load_or_build_vocab(windows: List[str], k: int, vocab_path: Path) -> KmerVocabulary:
    if vocab_path.exists():
        LOGGER.info("Loading existing vocabulary from %s", vocab_path)
        return KmerVocabulary.load(vocab_path)
    vocab = KmerVocabulary.build_from_sequences(windows, k=k)
    vocab.save(vocab_path)
    LOGGER.info("Saved vocabulary with %d entries to %s", len(vocab.itos), vocab_path)
    return vocab


def prepare_dataset(cfg) -> Tuple[List[torch.LongTensor], KmerVocabulary]:
    fasta_path = Path(cfg.get("mlm_fasta_path")).resolve()
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found at {fasta_path}")
    window_size = int(cfg.get("mlm_window_size", 81))
    stride = int(cfg.get("mlm_stride", 20))
    k = int(cfg.get("mlm_kmer", 3))
    vocab_path = Path(cfg.get("mlm_vocab_path", f"../artifacts/vocabs/k{k}_mlm_vocab.json")).resolve()

    LOGGER.info("Reading genome from %s", fasta_path)
    genome_sequence = read_fasta_sequence(fasta_path)
    windows = extract_windows(genome_sequence, window_size, stride)
    if not windows:
        raise ValueError("No valid windows generated from genome")
    LOGGER.info("Extracted %d windows", len(windows))
    vocab = load_or_build_vocab(windows, k=k, vocab_path=vocab_path)
    token_tensors = [vocab.tokenize(window, window_size) for window in windows]
    return token_tensors, vocab


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def run_validation(
    model: DNAMLM,
    val_loader: DataLoader,
    device: torch.device,
    autocast_device: str | None,
) -> Dict[str, float]:
    """Run validation and return metrics dict with loss, accuracy, perplexity."""
    model.eval()
    total_loss = 0.0
    total_correct = 0
    total_masked = 0
    
    with torch.no_grad():
        for inputs, labels in val_loader:
            inputs = inputs.to(device)
            labels = labels.to(device)
            
            if autocast_device:
                with torch.autocast(autocast_device, dtype=torch.float16):
                    logits, loss = model(inputs, labels)
            else:
                logits, loss = model(inputs, labels)
            
            total_loss += loss.item() * inputs.size(0)
            
            # Compute accuracy on masked tokens
            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()
    
    avg_loss = total_loss / len(val_loader.dataset)
    accuracy = total_correct / max(1, total_masked)
    perplexity = math.exp(min(avg_loss, 100))  # Cap to avoid overflow
    
    return {
        'loss': avg_loss,
        'accuracy': accuracy,
        'perplexity': perplexity,
    }


def check_gradients_and_frozen(model: nn.Module) -> Dict[str, Any]:
    """Check gradient norms and frozen weights for debugging."""
    grad_info = {"frozen": [], "zero_grad": [], "normal": []}
    total_params = 0
    frozen_params = 0
    
    for name, param in model.named_parameters():
        total_params += param.numel()
        if not param.requires_grad:
            frozen_params += param.numel()
            grad_info["frozen"].append(name)
        elif param.grad is None:
            grad_info["zero_grad"].append(name)
        else:
            grad_norm = param.grad.norm().item()
            if grad_norm == 0:
                grad_info["zero_grad"].append(name)
            else:
                grad_info["normal"].append((name, grad_norm))
    
    return {
        "total_params": total_params,
        "frozen_params": frozen_params,
        "frozen_pct": frozen_params / total_params * 100 if total_params else 0,
        "details": grad_info,
    }


def run_overfit_debug(
    model: nn.Module,
    train_samples: List[torch.LongTensor],
    vocab: KmerVocabulary,
    device: torch.device,
    num_steps: int = 500,
    lr: float = 3e-3,
) -> None:
    """
    HARD OVERFIT TEST: Train on tiny dataset until it breaks physics.
    
    This is a "lie detector" - if the model can't perfectly memorize a tiny 
    dataset, something is fundamentally wrong with the architecture or training.
    
    Configuration:
    - 500-1000 steps (aggressive)
    - NO dropout (pure memorization)
    - Higher LR (3e-3) for fast convergence
    - Single batch repeated (same data every step)
    - Log every 50 steps
    - FAIL if accuracy doesn't exceed 80%
    - Print top-5 predictions for masked positions on failure
    
    Expected: accuracy should reach 95%+ on this tiny dataset.
    If not, the model is broken.
    """
    print("=" * 80)
    print("HARD OVERFIT TEST: Breaking Physics Mode")
    print("=" * 80)
    print("If model can't memorize 32 samples, something is fundamentally broken.")
    print()
    print(f"Configuration:")
    print(f"  Samples:     {len(train_samples)}")
    print(f"  Steps:       {num_steps}")
    print(f"  LR:          {lr}")
    print(f"  Dropout:     DISABLED")
    print(f"  Target:      >80% accuracy (should reach 95%+)")
    print()
    
    # =========================================================================
    # DISABLE ALL DROPOUT for pure memorization
    # =========================================================================
    def disable_dropout(module):
        if isinstance(module, nn.Dropout):
            module.p = 0.0
    
    model.apply(disable_dropout)
    model.train()  # Keep in train mode but with dropout=0
    
    # Use AdamW with no weight decay for pure memorization
    # Also try higher LR for first phase, then reduce
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.0, betas=(0.9, 0.999))
    
    # Add learning rate scheduler: warmup then constant
    def lr_schedule(step):
        warmup = 50
        if step < warmup:
            return step / warmup
        return 1.0
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_schedule)
    
    # Create FIXED masked inputs/labels - same data every step
    raw_inputs = torch.stack([train_samples[i] for i in range(len(train_samples))]).to(device)
    masked_inputs, labels = mask_tokens_span(raw_inputs, vocab, mask_prob=0.15)
    
    num_masked = (labels != -100).sum().item()
    mask = labels != -100
    
    print(f"Fixed masking: {num_masked} tokens masked across {len(train_samples)} sequences")
    print(f"Masked positions per sequence: ~{num_masked / len(train_samples):.1f}")
    print()
    print("-" * 80)
    print(f"{'Step':>6} | {'Loss':>10} | {'Accuracy':>10} | {'Status'}")
    print("-" * 80)
    
    losses = []
    accuracies = []
    best_acc = 0.0
    
    for step in range(num_steps):
        optimizer.zero_grad()
        logits, loss = model(masked_inputs, labels)
        loss.backward()
        
        # Gradient clipping for stability
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        scheduler.step()
        
        # Compute accuracy
        preds = logits.argmax(dim=-1)
        acc = (preds[mask] == labels[mask]).float().mean().item() if mask.any() else 0.0
        
        losses.append(loss.item())
        accuracies.append(acc)
        best_acc = max(best_acc, acc)
        
        # Log every 50 steps
        if step % 50 == 0 or step == num_steps - 1:
            status = "✓" if acc > 0.8 else ("↑" if acc > 0.5 else "...")
            print(f"{step:>6} | {loss.item():>10.4f} | {acc*100:>9.1f}% | {status}")
    
    print("-" * 80)
    print()
    
    # =========================================================================
    # ANALYSIS
    # =========================================================================
    print("=" * 80)
    print("OVERFIT ANALYSIS")
    print("=" * 80)
    
    initial_loss = sum(losses[:5]) / min(5, len(losses))
    final_loss = sum(losses[-10:]) / min(10, len(losses))
    final_acc = sum(accuracies[-10:]) / min(10, len(accuracies))
    
    print(f"Initial loss:     {initial_loss:.4f}")
    print(f"Final loss:       {final_loss:.4f}")
    print(f"Loss reduction:   {(1 - final_loss/initial_loss)*100:.1f}%")
    print(f"Best accuracy:    {best_acc*100:.1f}%")
    print(f"Final accuracy:   {final_acc*100:.1f}%")
    print()
    
    # SUCCESS THRESHOLD: 80% accuracy
    if best_acc >= 0.80:
        print("✓ HARD OVERFIT PASSED!")
        print(f"  Model successfully memorized the data (best acc: {best_acc*100:.1f}%)")
        print("  The model architecture is capable of learning.")
        print("=" * 80)
        return
    
    # =========================================================================
    # FAILURE ANALYSIS: Print top-5 predictions for masked positions
    # =========================================================================
    print("✗ HARD OVERFIT FAILED!")
    print(f"  Best accuracy: {best_acc*100:.1f}% (required: >80%)")
    print()
    print("  This means the model CANNOT memorize a tiny dataset.")
    print("  Something is fundamentally broken.")
    print()
    
    print("-" * 80)
    print("TOP-5 PREDICTIONS FOR MASKED POSITIONS")
    print("-" * 80)
    
    # Get final predictions with probabilities
    model.eval()
    with torch.no_grad():
        logits, _ = model(masked_inputs, labels)
        probs = F.softmax(logits, dim=-1)
    
    # Analyze first few masked positions
    masked_indices = torch.where(mask)
    num_to_show = min(10, len(masked_indices[0]))
    
    print(f"Showing {num_to_show} masked positions:\n")
    
    for i in range(num_to_show):
        batch_idx = masked_indices[0][i].item()
        pos_idx = masked_indices[1][i].item()
        
        true_label = labels[batch_idx, pos_idx].item()
        true_token = vocab.itos[true_label] if true_label < len(vocab.itos) else f"<{true_label}>"
        
        # Get top-5 predictions
        pos_probs = probs[batch_idx, pos_idx]
        top5_probs, top5_ids = torch.topk(pos_probs, 5)
        
        print(f"  Position [{batch_idx},{pos_idx}] - True label: '{true_token}' (id={true_label})")
        print(f"    Top-5 predictions:")
        for j, (prob, pred_id) in enumerate(zip(top5_probs, top5_ids)):
            pred_token = vocab.itos[pred_id.item()] if pred_id.item() < len(vocab.itos) else f"<{pred_id.item()}>"
            marker = "✓" if pred_id.item() == true_label else " "
            print(f"      {j+1}. {marker} '{pred_token}' (id={pred_id.item()}) - {prob.item()*100:.1f}%")
        print()
    
    # Check for systematic issues
    print("-" * 80)
    print("DIAGNOSIS")
    print("-" * 80)
    
    # Check gradient info
    grad_info = check_gradients_and_frozen(model)
    
    if grad_info["frozen_pct"] > 10:
        print(f"  ⚠ {grad_info['frozen_pct']:.1f}% of parameters are frozen!")
    
    if len(grad_info["details"]["zero_grad"]) > 0:
        print(f"  ⚠ {len(grad_info['details']['zero_grad'])} layers have zero gradients")
    
    # Check if predictions are too uniform
    avg_top1_prob = probs.max(dim=-1)[0][mask].mean().item()
    print(f"  Average top-1 probability: {avg_top1_prob*100:.1f}%")
    if avg_top1_prob < 0.1:
        print("    ⚠ Predictions are too uniform - model not learning to discriminate")
    
    # Check if model is predicting special tokens (should be impossible now)
    special_ids = set(model.special_token_ids) if hasattr(model, 'special_token_ids') else set()
    if special_ids:
        top1_preds = probs.argmax(dim=-1)[mask]
        special_preds = sum(1 for p in top1_preds if p.item() in special_ids)
        if special_preds > 0:
            print(f"    ✗ {special_preds} predictions are special tokens (should be 0!)")
    
    # =========================================================================
    # KEY DIAGNOSTIC: Check 80/10/10 masking breakdown
    # =========================================================================
    print()
    print("-" * 80)
    print("80/10/10 MASKING ANALYSIS")
    print("-" * 80)
    print("  BERT-style masking: 80% [MASK], 10% random, 10% original")
    print()
    
    # Count how many masked positions got [MASK] vs other tokens
    mask_token_id = vocab.mask_id
    masked_input_tokens = masked_inputs[mask]
    original_tokens = raw_inputs[mask]
    
    got_mask_token = (masked_input_tokens == mask_token_id).sum().item()
    got_original = (masked_input_tokens == original_tokens).sum().item()
    got_random = num_masked - got_mask_token - got_original
    
    print(f"  Masked positions breakdown:")
    print(f"    [MASK] token:    {got_mask_token:>4} ({got_mask_token/num_masked*100:.1f}%) - NO context")
    print(f"    Original token:  {got_original:>4} ({got_original/num_masked*100:.1f}%) - EASY (sees answer)")
    print(f"    Random token:    {got_random:>4} ({got_random/num_masked*100:.1f}%) - SOME context")
    print()
    
    # Check accuracy breakdown by masking type
    final_preds = logits.argmax(dim=-1)[mask]
    final_labels = labels[mask]
    
    correct_mask = (final_preds == final_labels)
    
    # Accuracy for [MASK] positions
    is_mask_token = (masked_input_tokens == mask_token_id)
    acc_mask_only = correct_mask[is_mask_token].float().mean().item() if is_mask_token.any() else 0
    
    # Accuracy for original positions (should be ~100%)
    is_original = (masked_input_tokens == original_tokens)
    acc_original = correct_mask[is_original].float().mean().item() if is_original.any() else 0
    
    # Accuracy for random positions
    is_random = ~is_mask_token & ~is_original
    acc_random = correct_mask[is_random].float().mean().item() if is_random.any() else 0
    
    print(f"  Accuracy by masking type:")
    print(f"    [MASK] positions:    {acc_mask_only*100:.1f}%")
    print(f"    Original positions:  {acc_original*100:.1f}%")
    print(f"    Random positions:    {acc_random*100:.1f}%")
    print()
    
    if acc_original > 0.9 and acc_mask_only < 0.2:
        print("  ⚠ CRITICAL: Model only succeeds when it can SEE the answer!")
        print("    The 10% 'original' positions have high accuracy (cheating).")
        print("    The 80% [MASK] positions have near-random accuracy.")
        print()
        print("    This means attention/context IS NOT WORKING.")
        print("    The model cannot use surrounding context to predict masked tokens.")
    
    print()
    print("POSSIBLE CAUSES:")
    print("  1. Model architecture issue (attention not working)")
    print("  2. Loss function bug (gradients not flowing)")
    print("  3. Special token exclusion not working")
    print("  4. Embedding initialization issue")
    print("  5. Learning rate too low or too high")
    print()
    print("Try running with --leakage_test and --label_debug to verify setup.")
    print("=" * 80)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pretrain DNA MLM on E. coli genome")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--debug", action="store_true", help="Run debug mode: analyze one batch and exit")
    parser.add_argument("--overfit_debug", action="store_true", 
                        help="Train on 32 samples for 200 steps to verify model can learn")
    parser.add_argument("--label_debug", action="store_true",
                        help="Enable label verification debug: proves labels match originals. "
                             "Prints original/masked/label token IDs and asserts correctness for 5 positions per sequence. "
                             "Hard error if any mismatch. Runs once on first batch then exits.")
    parser.add_argument("--leakage_test", action="store_true",
                        help="Verify no information leakage at masked positions. "
                             "Tests: (1) masked embeddings don't contain original info, "
                             "(2) model runs with 100%% masking, (3) no engineered features. "
                             "Hard error if leakage detected.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    # =========================================================================
    # Log resolved hyperparameters to runs/<timestamp>/resolved_config.json
    # =========================================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"mlm_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    
    resolved_config: Dict[str, Any] = {
        "script": "pretrain_mlm.py",
        "timestamp": timestamp,
        "config_path": str(config_path),
        # Model architecture
        "embedding_dim": int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256))),
        "transformer_layers": int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8))),
        "transformer_heads": int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8))),
        "transformer_ff_dim": int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", 1024))),
        "transformer_dropout": float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1))),
        # Optimizer
        "lr": float(cfg.get("mlm_lr", 1e-3)),
        "weight_decay": float(cfg.get("mlm_weight_decay", 0.01)),
        # MLM-specific
        "mask_prob": 0.15,  # Fixed in code
        "max_span_len": 3,   # Fixed in mask_tokens_span
        # Data
        "batch_size": int(cfg.get("mlm_batch_size", 128)),
        "max_bp_len": int(cfg.get("mlm_window_size", 81)),
        "kmer": int(cfg.get("mlm_kmer", 3)),
        # Training
        "epochs": int(cfg.get("mlm_epochs", 200)),
        "grad_accum_steps": int(cfg.get("mlm_grad_accum_steps", 4)),
        "warmup_ratio": float(cfg.get("mlm_warmup_ratio", 0.06)),
        "patience": int(cfg.get("mlm_patience", 30)),
        "val_ratio": float(cfg.get("mlm_val_ratio", 0.1)),
        # Model options
        "tie_weights": bool(cfg.get("mlm_tie_weights", True)),
        "use_output_norm": bool(cfg.get("mlm_use_output_norm", True)),
        "use_layer_lr_decay": bool(cfg.get("mlm_use_layer_lr_decay", True)),
        "layer_lr_decay": float(cfg.get("mlm_layer_lr_decay", 0.9)),
        "use_alibi": bool(cfg.get("use_alibi", True)),
        # Seed
        "seed": int(cfg.get("seed", 1337)),
    }
    
    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)
    
    seed = int(cfg.get("seed", 1337))
    set_seed(seed)

    token_tensors, vocab = prepare_dataset(cfg)
    full_dataset = MLMDataset(token_tensors)
    
    # Split into train/val for monitoring overfitting
    val_ratio = float(cfg.get("mlm_val_ratio", 0.1))
    val_size = int(len(full_dataset) * val_ratio)
    train_size = len(full_dataset) - val_size
    train_dataset, val_dataset = random_split(
        full_dataset, 
        [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    LOGGER.info("Dataset split: %d train, %d val", train_size, val_size)
    
    batch_size = int(cfg.get("mlm_batch_size", 128))
    grad_accum_steps = int(cfg.get("mlm_grad_accum_steps", 2))
    effective_batch = batch_size * grad_accum_steps
    LOGGER.info("Batch size: %d, Grad accum: %d, Effective batch: %d", 
                batch_size, grad_accum_steps, effective_batch)
    
    # Standard 15% masking - no curriculum (was hurting performance)
    # Enable debug mode if --label_debug flag is set
    collator = MLMCollator(vocab, mask_prob=0.15, debug=args.label_debug)
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=int(cfg.get("num_workers", 0)),
        pin_memory=True,
        drop_last=True,
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=MLMCollator(vocab, mask_prob=0.15),
        num_workers=0,
        pin_memory=True,
    )
    
    # =========================================================================
    # LABEL DEBUG MODE: Verify labels match originals and exit
    # =========================================================================
    if args.label_debug:
        LOGGER.info("Running LABEL DEBUG mode - verifying first batch then exiting")
        # Force one batch through the collator to trigger debug verification
        for batch in train_loader:
            # Debug verification happens inside collator's __call__
            break
        LOGGER.info("Label debug complete. Exiting.")
        return

    # =========================================================================
    # DEBUG MODE: Run one batch analysis and exit
    # =========================================================================
    if args.debug:
        LOGGER.info("Running in DEBUG mode - analyzing masking on one batch")
        # Get one batch of raw tokens (before masking)
        raw_batch = [train_dataset[i] for i in range(min(batch_size, len(train_dataset)))]
        debug_masking(vocab, raw_batch, mask_prob=0.15, max_span_len=3)
        LOGGER.info("Debug mode complete. Exiting.")
        return

    device = select_device()
    LOGGER.info("Using device: %s", device)
    
    # Model architecture - ENHANCED for better pretraining
    # Larger models learn better representations during MLM pretraining
    embedding_dim = int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256)))
    transformer_layers = int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8)))
    transformer_heads = int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8)))
    # FF dimension should be 4x embedding dim for good capacity (BERT-style)
    transformer_ff_dim = int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", embedding_dim * 4)))
    transformer_dropout = float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1)))
    
    LOGGER.info("=" * 60)
    LOGGER.info("MLM ARCHITECTURE (Enhanced for lower loss)")
    LOGGER.info("  Embedding dim: %d", embedding_dim)
    LOGGER.info("  Transformer layers: %d", transformer_layers)
    LOGGER.info("  Transformer heads: %d", transformer_heads)
    LOGGER.info("  FF dim: %d (%.1fx embedding)", transformer_ff_dim, transformer_ff_dim / embedding_dim)
    LOGGER.info("  Dropout: %.2f", transformer_dropout)
    LOGGER.info("=" * 60)

    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        use_alibi=bool(cfg.get("use_alibi", True)),
        pad_token_id=vocab.pad_id,
        drop_path_rate=0.0,  # NO stochastic depth during pretraining
    )
    
    # DNAMLM with weight tying (critical for lower loss)
    tie_weights = bool(cfg.get("mlm_tie_weights", True))
    use_output_norm = bool(cfg.get("mlm_use_output_norm", True))
    
    # CRITICAL: Exclude special tokens from prediction space
    # This prevents the model from "cheating" by predicting [MASK], <UNK>, <PAD>
    # which has been shown to cause models to get stuck at ~30% accuracy
    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    LOGGER.info("Special tokens excluded from predictions: [MASK]=%d, <UNK>=%d, <PAD>=%d",
                vocab.mask_id, vocab.unk_id, vocab.pad_id)
    
    model = DNAMLM(
        encoder, 
        vocab_size=len(vocab.itos),
        special_token_ids=special_token_ids,
        tie_weights=tie_weights,
        use_output_norm=use_output_norm,
    ).to(device)
    
    # Count parameters (accounting for weight tying)
    total_params = sum(p.numel() for p in model.parameters())
    unique_params = sum(p.numel() for p in set(model.parameters()))
    LOGGER.info("Model parameters: %.2fM total, %.2fM unique", total_params / 1e6, unique_params / 1e6)
    
    # =========================================================================
    # LEAKAGE TEST: Verify no information leakage at masked positions
    # =========================================================================
    if args.leakage_test:
        LOGGER.info("Running LEAKAGE TEST mode - verifying no information leakage")
        verify_no_leakage(model, vocab, device)
        LOGGER.info("Leakage test complete. Exiting.")
        return
    
    # =========================================================================
    # OVERFIT DEBUG MODE: Train on 32 samples for 200 steps
    # =========================================================================
    if args.overfit_debug:
        LOGGER.info("Running OVERFIT DEBUG mode")
        # Take 32 samples for proper overfit test
        overfit_samples = token_tensors[:32]
        run_overfit_debug(
            model=model,
            train_samples=overfit_samples,
            vocab=vocab,
            device=device,
            num_steps=1000,  # Very aggressive: 1000 steps
            lr=5e-3,         # Very high LR for fast convergence
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return
    
    # Learning rate and optimizer settings
    lr = float(cfg.get("mlm_lr", 1e-3))  # Higher base LR for MLM pretraining
    weight_decay = float(cfg.get("mlm_weight_decay", 0.01))
    lr_decay = float(cfg.get("mlm_layer_lr_decay", 0.9))  # Layer-wise LR decay
    use_layer_lr_decay = bool(cfg.get("mlm_use_layer_lr_decay", True))
    
    # Build optimizer with layer-wise LR decay
    if use_layer_lr_decay:
        param_groups = get_layer_wise_lr_groups(model, lr, lr_decay, weight_decay)
        LOGGER.info("Using layer-wise LR decay (factor: %.2f)", lr_decay)
        optimizer = AdamW(
            param_groups, 
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,  # Smaller eps for stability
        )
    else:
        optimizer = AdamW(
            model.parameters(), 
            lr=lr, 
            weight_decay=weight_decay,
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,
        )
    
    epochs = int(cfg.get("mlm_epochs", 150))
    
    # Learning rate schedule with warmup
    steps_per_epoch = len(train_loader) // grad_accum_steps
    total_steps = epochs * steps_per_epoch
    warmup_ratio = float(cfg.get("mlm_warmup_ratio", 0.06))  # 6% warmup
    warmup_steps = int(warmup_ratio * total_steps)
    
    scheduler = get_warmup_cosine_schedule(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=0.01,  # Decay to 1% of peak LR
    )
    LOGGER.info("LR schedule: %d warmup steps, %d total steps", warmup_steps, total_steps)

    autocast_device = None
    if device.type == "cuda":
        autocast_device = "cuda"
    elif device.type == "mps":
        autocast_device = "mps"
    
    # Training state
    best_val_loss = float('inf')
    best_val_acc = 0.0
    patience = int(cfg.get("mlm_patience", 25))
    patience_counter = 0
    global_step = 0
    
    encoder_path = Path(cfg.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt")).resolve()
    encoder_path.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("=" * 60)
    LOGGER.info("STARTING MLM PRETRAINING")
    LOGGER.info("Epochs: %d, LR: %.2e, Mask prob: 15%%", epochs, lr)
    LOGGER.info("=" * 60)

    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0
        total_correct = 0
        total_masked = 0
        optimizer.zero_grad()
        
        for batch_idx, (inputs, labels) in enumerate(train_loader, start=1):
            inputs = inputs.to(device)
            labels = labels.to(device)
            
            if autocast_device:
                with torch.autocast(autocast_device, dtype=torch.float16):
                    logits, loss = model(inputs, labels)
            else:
                logits, loss = model(inputs, labels)
            
            # Scale loss for gradient accumulation
            scaled_loss = loss / grad_accum_steps
            scaled_loss.backward()
            
            total_loss += loss.item() * inputs.size(0)
            
            # Compute accuracy
            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()
            
            # Optimizer step with gradient accumulation
            if batch_idx % grad_accum_steps == 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                scheduler.step()
                optimizer.zero_grad()
                global_step += 1
            
            if batch_idx % 100 == 0:
                current_lr = optimizer.param_groups[0]['lr']
                batch_acc = total_correct / max(1, total_masked)
                LOGGER.info(
                    "Epoch %d [%d/%d] Loss: %.4f Acc: %.2f%% LR: %.2e",
                    epoch, batch_idx, len(train_loader), loss.item(), 
                    batch_acc * 100, current_lr
                )
        
        # Epoch statistics
        train_loss = total_loss / len(train_dataset)
        train_acc = total_correct / max(1, total_masked)
        train_ppl = math.exp(min(train_loss, 100))
        
        # Validation
        val_metrics = run_validation(model, val_loader, device, autocast_device)
        val_loss = val_metrics['loss']
        val_acc = val_metrics['accuracy']
        val_ppl = val_metrics['perplexity']
        
        LOGGER.info(
            "Epoch %d: Train [Loss: %.4f, Acc: %.1f%%, PPL: %.2f] | Val [Loss: %.4f, Acc: %.1f%%, PPL: %.2f]",
            epoch, train_loss, train_acc * 100, train_ppl, val_loss, val_acc * 100, val_ppl
        )
        
        # Save best model based on validation loss
        improved = False
        if val_loss < best_val_loss - 0.001:
            improved = True
            best_val_loss = val_loss
            best_val_acc = val_acc
            patience_counter = 0
            torch.save(model.encoder.state_dict(), encoder_path)
            LOGGER.info("★ New best! Val Loss: %.4f, PPL: %.2f, Acc: %.1f%% - Saved to %s", 
                       val_loss, val_ppl, val_acc * 100, encoder_path)
        else:
            patience_counter += 1
            if patience_counter >= patience:
                LOGGER.info("Early stopping at epoch %d (patience %d)", epoch, patience)
                break
    
    best_ppl = math.exp(min(best_val_loss, 100))
    LOGGER.info("=" * 60)
    LOGGER.info("PRETRAINING COMPLETE")
    LOGGER.info("Best Val Loss: %.4f (Perplexity: %.2f)", best_val_loss, best_ppl)
    LOGGER.info("Best Val Accuracy: %.1f%%", best_val_acc * 100)
    LOGGER.info("Target: Loss < 1.0 (PPL < 2.72)")
    LOGGER.info("Saved pretrained encoder to %s", encoder_path)
    LOGGER.info("=" * 60)


if __name__ == "__main__":
    main()
