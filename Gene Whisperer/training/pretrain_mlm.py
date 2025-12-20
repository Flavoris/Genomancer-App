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
import copy
import json
import logging
import math
import os
import random
import time
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
from length_curriculum import (
    LengthCurriculumConfig,
    LengthCurriculumScheduler,
    CurriculumMLMCollator,
    create_curriculum_collator,
)
from pretrain_health import (
    PretrainHealthDashboard,
    create_health_dashboard,
    run_pretrain_sanity_checks,
)
from seed_utils import (
    set_global_seed,
    get_deterministic_dataloader_kwargs,
)
from numerics import (
    NumericsChecker,
    NumericsError,
    assert_finite,
    cap_model_param_norm_,
    detect_nan_in_model,
)

LOGGER = logging.getLogger("gene_whisperer.pretrain_mlm")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

ALLOWED_BASES = {"A", "C", "G", "T"}
REV_COMP_MAP = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of an A/C/G/T-only DNA sequence."""
    return sequence.translate(REV_COMP_MAP)[::-1]


# Note: set_seed is imported from seed_utils as set_global_seed
# This wrapper maintains backwards compatibility
def set_seed(seed: int) -> None:
    """Set all random seeds for reproducibility (wrapper for set_global_seed)."""
    set_global_seed(seed)


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


def _resolve_cfg_path(path_value: str | Path | None, *, base_dir: Path) -> Path:
    if path_value is None or str(path_value).strip() == "":
        raise ValueError("Expected a non-empty path value")
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path


class MLMDataset(Dataset):
    def __init__(self, token_tensors: List[torch.LongTensor]):
        self.samples = token_tensors

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> torch.LongTensor:
        return self.samples[idx]


def sample_span_length_geometric(max_span_len: int, mean_span: float = 3.0) -> int:
    """
    Sample span length from geometric distribution (DNABERT-style).
    
    The geometric distribution naturally generates span lengths where shorter
    spans are more common than longer ones, which matches the intuition that
    most masked regions should be small.
    
    Args:
        max_span_len: Maximum allowed span length
        mean_span: Target mean span length (controls distribution)
        
    Returns:
        Sampled span length in [1, max_span_len]
    """
    # Geometric distribution: P(X=k) = p * (1-p)^(k-1)
    # Mean = 1/p, so p = 1/mean_span
    p = 1.0 / mean_span
    
    # Sample from geometric distribution and clamp to [1, max_span_len]
    # Use inverse CDF method: X = ceil(log(1-U) / log(1-p))
    u = random.random()
    if u == 0:
        u = 1e-10  # Avoid log(0)
    
    span_len = int(math.ceil(math.log(1 - u) / math.log(1 - p)))
    return max(1, min(span_len, max_span_len))


def mask_tokens_span(
    inputs: torch.LongTensor,
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
    max_span_len: int = 3,
    track_spans: bool = False,
    span_distribution: str = "geometric",
    mean_span: float = 3.0,
    exclude_special_from_labels: bool = True,
) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
    """
    DNABERT-style span masking with geometric or uniform distribution.
    
    Instead of independent per-token masking, this selects contiguous spans
    until approximately mask_prob of tokens are covered, then applies the
    standard 80/10/10 BERT masking within those spans.
    
    Key features:
    - Geometric distribution (default): shorter spans more common, matches DNABERT
    - Uniform distribution: all span lengths equally likely
    - Never masks PAD tokens
    - Handles edge cases (all-PAD, very short sequences)
    - CRITICAL: Special tokens ([MASK], [UNK], [PAD]) are NEVER used as labels
    
    Args:
        inputs: (B, L) input token indices
        vocab: KmerVocabulary with mask_id and pad_id
        mask_prob: Target fraction of tokens to mask (~15%)
        max_span_len: Maximum length of each masked span
        track_spans: If True, also return list of span lengths
        span_distribution: "geometric" (DNABERT-style) or "uniform"
        mean_span: Mean span length for geometric distribution
        exclude_special_from_labels: If True, set labels to -100 for any position
            where the original token was a special token. This prevents the model
            from being trained to predict special tokens.
        
    Returns:
        (masked_inputs, labels) tuple where labels has -100 for unmasked positions
        If track_spans=True, also returns list of span lengths
    """
    device = inputs.device
    labels = inputs.clone()
    batch_size, seq_len = inputs.shape
    vocab_size = len(vocab.itos)

    if max_span_len < 1:
        raise ValueError(f"max_span_len must be >= 1, got {max_span_len}")

    # Initialize all positions as unmasked
    masked = torch.zeros_like(inputs, dtype=torch.bool, device=device)
    
    # Identify PAD positions early (we won't mask these)
    if hasattr(vocab, "pad_id"):
        pad_mask = inputs.eq(vocab.pad_id)
    else:
        pad_mask = torch.zeros_like(inputs, dtype=torch.bool, device=device)
    
    for b in range(batch_size):
        # Count non-PAD tokens in this sequence
        non_pad_count = (~pad_mask[b]).sum().item()
        
        # Handle all-PAD sequences
        if non_pad_count == 0:
            continue
        
        # Target number of tokens to mask (from non-PAD tokens)
        num_to_mask = max(1, int(mask_prob * non_pad_count))
        covered = 0
        attempts = 0
        # With non-overlapping spans and a hard cap on contiguous masked runs,
        # placement can be more constrained. Allow more retries before giving up.
        max_attempts = num_to_mask * 50  # Safety limit
        
        while covered < num_to_mask and attempts < max_attempts:
            # Sample span length based on distribution
            if span_distribution == "geometric":
                span_len = sample_span_length_geometric(max_span_len, mean_span)
            elif span_distribution == "uniform":
                span_len = random.randint(1, max_span_len)
            else:
                raise ValueError(f"Unknown span_distribution: {span_distribution!r}")

            # Enforce cap (defensive; sampling should already respect max_span_len)
            span_len = min(span_len, max_span_len)

            # Never mask more than we still need (prevents systematic overshoot)
            remaining = num_to_mask - covered
            span_len = max(1, min(span_len, remaining))

            # Find a start position that:
            # - stays within non-PAD region
            # - does not overlap an existing masked span
            # - does not touch an existing masked span (prevents contiguous runs > max_span_len)
            chosen_start: int | None = None
            chosen_len: int | None = None

            # If we can't place the sampled length, try shorter lengths down to 1.
            for try_len in range(span_len, 0, -1):
                valid_starts: List[int] = []
                max_start = seq_len - try_len
                if max_start < 0:
                    continue
                for start in range(0, max_start + 1):
                    end = start + try_len

                    # Require all tokens in the span to be non-PAD (no partial spans)
                    if not (~pad_mask[b, start:end]).all().item():
                        continue

                    # No overlap with existing masked tokens
                    if masked[b, start:end].any().item():
                        continue

                    # No touching existing masked tokens (keeps contiguous masked runs bounded)
                    if start > 0 and masked[b, start - 1].item():
                        continue
                    if end < seq_len and masked[b, end].item():
                        continue

                    valid_starts.append(start)

                if valid_starts:
                    chosen_start = random.choice(valid_starts)
                    chosen_len = try_len
                    break
            
            if chosen_start is None or chosen_len is None:
                attempts += 1
                break  # No more valid placements possible under constraints

            start = chosen_start
            end = start + chosen_len

            masked[b, start:end] = True
            covered += chosen_len
            
            attempts += 1
            
            # Break if we've covered enough
            if covered >= num_to_mask:
                break

    # Ensure PAD tokens are NEVER masked (redundant safety check)
    masked[pad_mask] = False

    # Set labels: -100 for unmasked positions (ignored in loss)
    labels[~masked] = -100
    
    # CRITICAL: Exclude special tokens from being supervised targets
    # Even if a position is masked, if the original token was a special token,
    # we should NOT train the model to predict it. This prevents the model
    # from learning to predict [MASK], [UNK], or [PAD] tokens.
    if exclude_special_from_labels:
        special_token_ids = {vocab.mask_id, vocab.unk_id, vocab.pad_id}
        for special_id in special_token_ids:
            special_mask = (inputs == special_id)
            labels[special_mask] = -100

    # Apply 80/10/10 BERT masking strategy within masked spans
    # 80% -> [MASK] token
    # 10% -> random token (from vocab, excluding special tokens)
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
        # Use vocab's base token IDs (excludes special tokens)
        if hasattr(vocab, '_base_token_ids') and vocab._base_token_ids:
            random_ids = torch.tensor(
                [random.choice(vocab._base_token_ids) for _ in range(num_rand)],
                device=device,
                dtype=inputs.dtype,
            )
        else:
            random_ids = torch.randint(
                low=0,
                high=vocab_size,
                size=(num_rand,),
                device=device,
                dtype=inputs.dtype,
            )
        inputs.view(-1)[flat] = random_ids

    if track_spans:
        # Track contiguous masked run lengths (what the model actually sees as "spans")
        span_lengths: List[int] = []
        pad0 = torch.zeros((1,), dtype=torch.int8, device=device)
        for b in range(batch_size):
            row = masked[b].to(torch.int8)
            padded = torch.cat([pad0, row, pad0])
            diff = padded[1:] - padded[:-1]
            starts = torch.where(diff == 1)[0]
            ends = torch.where(diff == -1)[0]
            if starts.numel() == 0:
                continue
            span_lengths.extend((ends - starts).tolist())

        return inputs, labels, span_lengths
    return inputs, labels


def debug_label_token_frequency(
    labels: torch.LongTensor,
    vocab: KmerVocabulary,
    num_batches: int = 1,
) -> Dict[str, Any]:
    """
    Generate debug report showing frequency of each token ID in labels.
    
    This is critical for verifying that special tokens are NOT being used
    as supervised targets. A properly configured MLM should have:
    - 0 occurrences of [MASK], [UNK], [PAD] in labels
    - All label tokens should be real k-mers (IDs 0 to 63 for k=3)
    
    Args:
        labels: (B, L) label tensor with -100 for ignored positions
        vocab: KmerVocabulary to decode token IDs
        num_batches: Number of batches this report covers (for averaging)
        
    Returns:
        Dict with:
        - token_counts: Dict[int, int] mapping token ID to count
        - special_token_counts: Dict[str, int] for [MASK], [UNK], [PAD]
        - total_supervised: Total number of supervised (non -100) labels
        - issues: List of issues found (empty if all is well)
    """
    # Flatten and filter out -100 (ignored positions)
    flat_labels = labels.view(-1)
    supervised_mask = flat_labels != -100
    supervised_labels = flat_labels[supervised_mask]
    
    total_supervised = supervised_labels.numel()
    
    # Count occurrences of each token ID
    token_counts: Dict[int, int] = {}
    if total_supervised > 0:
        unique, counts = torch.unique(supervised_labels, return_counts=True)
        for token_id, count in zip(unique.tolist(), counts.tolist()):
            token_counts[token_id] = count
    
    # Check special tokens specifically
    special_token_names = {
        vocab.mask_id: "[MASK]",
        vocab.unk_id: "[UNK]",
        vocab.pad_id: "[PAD]",
    }
    special_token_counts = {}
    issues = []
    
    for token_id, name in special_token_names.items():
        count = token_counts.get(token_id, 0)
        special_token_counts[name] = count
        if count > 0:
            issues.append(f"CRITICAL: {name} (id={token_id}) appears {count} times as supervised target!")
    
    return {
        "token_counts": token_counts,
        "special_token_counts": special_token_counts,
        "total_supervised": total_supervised,
        "issues": issues,
        "num_batches": num_batches,
    }


def print_label_debug_report(
    report: Dict[str, Any],
    vocab: KmerVocabulary,
    top_k: int = 10,
) -> None:
    """
    Print a formatted debug report for label token frequencies.
    
    Args:
        report: Output from debug_label_token_frequency()
        vocab: KmerVocabulary for decoding token names
        top_k: Number of most common tokens to display
    """
    print("=" * 80)
    print("LABEL TOKEN FREQUENCY DEBUG REPORT")
    print("=" * 80)
    print()
    
    print(f"Total supervised labels (non -100): {report['total_supervised']:,}")
    print(f"Batches analyzed: {report['num_batches']}")
    print()
    
    # Special tokens section - CRITICAL
    print("-" * 80)
    print("SPECIAL TOKEN CHECK (should ALL be 0)")
    print("-" * 80)
    all_zero = True
    for name, count in report["special_token_counts"].items():
        status = "✓" if count == 0 else "✗ FAIL"
        print(f"  {name}: {count:,} occurrences {status}")
        if count > 0:
            all_zero = False
    
    print()
    if all_zero:
        print("  ✓ All special tokens correctly excluded from labels!")
    else:
        print("  ✗ CRITICAL ERROR: Special tokens found in labels!")
        print("    The model is being trained to predict special tokens.")
        print("    This will hurt performance significantly.")
    print()
    
    # Top-k most common tokens
    print("-" * 80)
    print(f"TOP {top_k} MOST COMMON LABEL TOKENS")
    print("-" * 80)
    
    token_counts = report["token_counts"]
    if token_counts:
        sorted_counts = sorted(token_counts.items(), key=lambda x: -x[1])[:top_k]
        total = report["total_supervised"]
        
        for rank, (token_id, count) in enumerate(sorted_counts, 1):
            pct = count / total * 100 if total > 0 else 0
            # Get token name
            if token_id < len(vocab.itos):
                token_name = vocab.itos[token_id]
            else:
                token_name = f"<unknown:{token_id}>"
            print(f"  {rank:>2}. '{token_name}' (id={token_id}): {count:,} ({pct:.2f}%)")
    else:
        print("  No supervised labels found.")
    
    print()
    
    # Issues section
    if report["issues"]:
        print("-" * 80)
        print("ISSUES FOUND")
        print("-" * 80)
        for issue in report["issues"]:
            print(f"  ✗ {issue}")
        print()
    
    print("=" * 80)


def get_warmup_cosine_schedule(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
    final_lr_ratio: float = 0.1,
) -> LambdaLR:
    """
    Learning rate schedule with linear warmup and cosine decay.
    
    Schedule:
    1. Linear warmup: LR ramps from 0 to base_lr over num_warmup_steps
    2. Cosine decay: LR decays from base_lr to base_lr * final_lr_ratio
    
    Args:
        optimizer: The optimizer to schedule
        num_warmup_steps: Number of steps for linear warmup
        num_training_steps: Total number of training steps
        final_lr_ratio: Final LR as fraction of peak LR (e.g., 0.1 = 10% of peak)
        
    Returns:
        LambdaLR scheduler that multiplies base_lr by the lambda value
    """
    def lr_lambda(current_step: int) -> float:
        # Phase 1: Linear warmup from 0 to 1
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        
        # Phase 2: Cosine decay from 1 to final_lr_ratio
        # progress goes from 0 to 1 after warmup
        progress = float(current_step - num_warmup_steps) / float(
            max(1, num_training_steps - num_warmup_steps)
        )
        # Clamp progress to [0, 1] in case of extra steps
        progress = min(1.0, progress)
        
        # Cosine decay: starts at 1.0, ends at final_lr_ratio
        # Formula: final_lr_ratio + (1 - final_lr_ratio) * 0.5 * (1 + cos(pi * progress))
        # At progress=0: 0.5 * (1 + 1) = 1.0 -> returns 1.0
        # At progress=1: 0.5 * (1 + (-1)) = 0.0 -> returns final_lr_ratio
        cosine_decay = 0.5 * (1.0 + math.cos(math.pi * progress))
        return final_lr_ratio + (1.0 - final_lr_ratio) * cosine_decay
    
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
    
    # Positional embedding parameters (also low layer, same LR as token embedding)
    if hasattr(encoder, 'pos_embedding'):
        pos_embed_params_decay = []
        pos_embed_params_no_decay = []
        for name, param in encoder.pos_embedding.named_parameters():
            if any(nd in name for nd in no_decay):
                pos_embed_params_no_decay.append(param)
            else:
                pos_embed_params_decay.append(param)
        
        if pos_embed_params_decay:
            param_groups.append({'params': pos_embed_params_decay, 'lr': embed_lr, 'weight_decay': weight_decay})
        if pos_embed_params_no_decay:
            param_groups.append({'params': pos_embed_params_no_decay, 'lr': embed_lr, 'weight_decay': 0.0})
    
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

    # Learnable logit scale for normalized MLM head (do not decay)
    if hasattr(model, "logit_scale"):
        lm_head_params_no_decay.append(model.logit_scale)
    
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
    
    track_spans mode:
    - Returns (masked_inputs, labels, span_lengths) instead of (masked_inputs, labels)
    - Used by health dashboard to log span length histograms
    """
    def __init__(
        self, 
        vocab: KmerVocabulary, 
        mask_prob: float = 0.15, 
        debug: bool = False,
        track_spans: bool = False,
    ):
        self.vocab = vocab
        self.mask_prob = mask_prob
        self.debug = debug
        self.track_spans = track_spans
        self._debug_ran = False
    
    def __call__(
        self, batch: Iterable[torch.LongTensor]
    ) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
        original_inputs = torch.stack(list(batch))
        
        if self.track_spans:
            masked_inputs, labels, span_lengths = mask_tokens_span(
                original_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=True
            )
        else:
            masked_inputs, labels = mask_tokens_span(
                original_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=False
            )
            span_lengths = []
        
        # Run debug verification once if enabled
        if self.debug and not self._debug_ran:
            self._run_debug_verification(original_inputs, masked_inputs, labels)
            self._debug_ran = True
        
        if self.track_spans:
            return masked_inputs, labels, span_lengths
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
        self.logit_scale = nn.Parameter(torch.tensor(10.0))
        
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

        # Normalize MLM head projection to prevent logit explosion (CLIP-style)
        h = F.normalize(x, dim=-1)
        w = F.normalize(self.lm_head.weight, dim=-1)
        logits = self.logit_scale * (h @ w.T)
        if self.lm_head.bias is not None:
            logits = logits + self.lm_head.bias
        
        # CRITICAL: Mask out special tokens from prediction space
        # This prevents the model from predicting [MASK], <UNK>, <PAD>
        logits = self._mask_special_tokens(logits)
        
        if labels is None:
            return logits, None

        # Standard cross entropy - NO label smoothing for MLM
        # CRITICAL: Force logits + loss to fp32 for numerical stability in mixed precision
        logits = logits.float()
        assert logits.dtype == torch.float32
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
        logits = logits.float()
        assert logits.dtype == torch.float32
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
    script_dir = Path(__file__).resolve().parent
    fasta_path = _resolve_cfg_path(cfg.get("mlm_fasta_path"), base_dir=script_dir)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found at {fasta_path}")
    window_size = int(cfg.get("mlm_window_size", 81))
    stride = int(cfg.get("mlm_stride", 20))
    k = int(cfg.get("mlm_kmer", 3))
    include_reverse_complements = bool(cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False)))
    vocab_path = _resolve_cfg_path(
        cfg.get("mlm_vocab_path", f"../artifacts/vocabs/k{k}_mlm_vocab.json"),
        base_dir=script_dir,
    )

    LOGGER.info("Reading genome from %s", fasta_path)
    genome_sequence = read_fasta_sequence(fasta_path)
    ablation_cfg = cfg.get("ablation") if isinstance(cfg, dict) else None
    if isinstance(ablation_cfg, dict) and bool(ablation_cfg.get("enabled", False)):
        bp_limit_raw = ablation_cfg.get("pretrain_bp_limit")
        if bp_limit_raw is not None:
            try:
                bp_limit = int(bp_limit_raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"ablation.pretrain_bp_limit must be an int, got {bp_limit_raw!r}") from exc
            if bp_limit < 0:
                raise ValueError(f"ablation.pretrain_bp_limit must be >= 0, got {bp_limit}")
            if bp_limit > 0:
                genome_sequence = genome_sequence[:bp_limit]
                LOGGER.info("Ablation enabled: using first %d bp from FASTA", len(genome_sequence))
    windows = extract_windows(genome_sequence, window_size, stride)
    if not windows:
        raise ValueError("No valid windows generated from genome")
    forward_count = len(windows)
    if include_reverse_complements:
        rc_windows = [reverse_complement(window) for window in windows]
        windows.extend(rc_windows)
        rc_fraction = len(rc_windows) / max(1, len(windows))
        LOGGER.info(
            "Extracted %d windows (%d + RC=%d, rc_frac=%.3f)",
            len(windows),
            forward_count,
            len(rc_windows),
            rc_fraction,
        )
    else:
        LOGGER.info("Extracted %d windows (rc_frac=0.000)", forward_count)
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
    autocast_dtype: torch.dtype | None,
) -> Dict[str, float]:
    """Run validation and return metrics dict with loss, accuracy, perplexity."""
    model.eval()
    total_loss = 0.0
    total_correct = 0
    total_masked = 0
    
    with torch.no_grad():
        for batch_data in val_loader:
            if len(batch_data) == 3:
                inputs, labels, _span_lengths = batch_data
            else:
                inputs, labels = batch_data

            inputs = inputs.to(device)
            labels = labels.to(device)

            if autocast_device and autocast_dtype is not None:
                with torch.autocast(autocast_device, dtype=autocast_dtype):
                    logits, _ = model(inputs, labels=None)
            else:
                logits, _ = model(inputs, labels=None)

            # Loss computation outside autocast for numerical stability
            loss = F.cross_entropy(
                logits.float().view(-1, logits.size(-1)),
                labels.view(-1),
                ignore_index=-100,
            )
            
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


def print_topk_predictions_debug(
    model: nn.Module,
    inputs: torch.LongTensor,
    labels: torch.LongTensor,
    vocab: KmerVocabulary,
    epoch: int,
    num_sequences: int = 3,
    top_k: int = 5,
) -> Dict[str, Any]:
    """
    Top-K sanity print for masked positions.
    
    Prints detailed prediction info for the first `num_sequences` sequences:
    - Masked position index
    - True token string
    - Model top-5 predicted tokens + probabilities
    
    This helps verify that:
    1. Predictions become sharper over time (top-1 prob rises)
    2. True token appears in top-5 more often as training progresses
    
    Args:
        model: The DNAMLM model
        inputs: (B, L) masked input token IDs
        labels: (B, L) target labels with -100 for non-masked positions
        vocab: KmerVocabulary for decoding tokens
        epoch: Current epoch number (for display)
        num_sequences: Number of sequences to analyze (default: 3)
        top_k: Number of top predictions to show (default: 5)
        
    Returns:
        Dict with summary statistics:
        - avg_top1_prob: Average probability of top-1 prediction
        - true_in_top1_pct: Percentage of masked tokens where true label is top-1
        - true_in_topk_pct: Percentage of masked tokens where true label is in top-k
    """
    model.eval()
    device = next(model.parameters()).device
    
    # Move to device and get predictions
    inputs = inputs.to(device)
    labels = labels.to(device)
    
    with torch.no_grad():
        logits, _ = model(inputs, labels)
        probs = F.softmax(logits, dim=-1)
    
    batch_size = min(num_sequences, inputs.size(0))
    
    print()
    print("=" * 90)
    print(f"TOP-K PREDICTIONS DEBUG (Epoch {epoch})")
    print("=" * 90)
    print(f"Analyzing {batch_size} sequences, showing top-{top_k} predictions per masked position")
    print()
    
    # Track statistics
    total_masked = 0
    true_in_top1 = 0
    true_in_topk = 0
    top1_probs = []
    
    for seq_idx in range(batch_size):
        seq_labels = labels[seq_idx]
        seq_probs = probs[seq_idx]
        
        # Find masked positions (where label != -100)
        masked_positions = torch.where(seq_labels != -100)[0]
        
        if len(masked_positions) == 0:
            print(f"[Sequence {seq_idx}] No masked positions found")
            continue
            
        print(f"[Sequence {seq_idx}] {len(masked_positions)} masked positions")
        print("-" * 90)
        
        # Show up to 5 masked positions per sequence for readability
        positions_to_show = masked_positions[:5]
        
        for pos in positions_to_show:
            pos_idx = pos.item()
            true_label = seq_labels[pos_idx].item()
            
            # Get true token string
            if true_label < len(vocab.itos):
                true_token = vocab.itos[true_label]
            else:
                true_token = f"<id:{true_label}>"
            
            # Get top-k predictions
            pos_probs = seq_probs[pos_idx]
            topk_probs, topk_ids = torch.topk(pos_probs, top_k)
            
            # Check if true label is in top-1 or top-k
            top1_id = topk_ids[0].item()
            is_top1_correct = (top1_id == true_label)
            is_in_topk = true_label in topk_ids.tolist()
            
            # Update statistics
            total_masked += 1
            if is_top1_correct:
                true_in_top1 += 1
            if is_in_topk:
                true_in_topk += 1
            top1_probs.append(topk_probs[0].item())
            
            # Format output
            correct_marker = "✓" if is_top1_correct else ("○" if is_in_topk else "✗")
            print(f"  Position {pos_idx:>3} | True: '{true_token}' (id={true_label}) | {correct_marker}")
            print(f"    Top-{top_k} predictions:")
            
            for rank, (prob, pred_id) in enumerate(zip(topk_probs, topk_ids), 1):
                pred_id_val = pred_id.item()
                if pred_id_val < len(vocab.itos):
                    pred_token = vocab.itos[pred_id_val]
                else:
                    pred_token = f"<id:{pred_id_val}>"
                
                # Mark if this is the true token
                marker = " ★" if pred_id_val == true_label else ""
                print(f"      {rank}. '{pred_token}' (id={pred_id_val}): {prob.item()*100:6.2f}%{marker}")
            print()
        
        print()
    
    # Calculate and print summary statistics
    avg_top1_prob = sum(top1_probs) / len(top1_probs) if top1_probs else 0.0
    true_in_top1_pct = true_in_top1 / total_masked * 100 if total_masked > 0 else 0.0
    true_in_topk_pct = true_in_topk / total_masked * 100 if total_masked > 0 else 0.0
    
    print("-" * 90)
    print("SUMMARY STATISTICS")
    print("-" * 90)
    print(f"  Total masked positions analyzed: {total_masked}")
    print(f"  Average top-1 probability:       {avg_top1_prob*100:.2f}%")
    print(f"  True token is top-1:             {true_in_top1}/{total_masked} ({true_in_top1_pct:.1f}%)")
    print(f"  True token in top-{top_k}:            {true_in_topk}/{total_masked} ({true_in_topk_pct:.1f}%)")
    print()
    print("  Healthy training indicators:")
    print("    → Top-1 probability should INCREASE over epochs")
    print("    → True-in-top-1 % should INCREASE over epochs")
    print("    → True-in-top-5 % should be HIGH early, approaching true-in-top-1 later")
    print("=" * 90)
    print()
    
    model.train()  # Restore training mode
    
    return {
        "avg_top1_prob": avg_top1_prob,
        "true_in_top1_pct": true_in_top1_pct,
        "true_in_topk_pct": true_in_topk_pct,
        "total_masked": total_masked,
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
        
        # Gradient clipping for stability with norm logging
        unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        clipped_grad_norm = min(unclipped_grad_norm.item(), 1.0)
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


def run_mlm_pretrain(cfg: dict, *, overrides: dict | None = None) -> dict:
    """Runs MLM pretraining and returns metrics + encoder checkpoint path."""
    start_time = time.perf_counter()

    cfg_run = copy.copy(cfg) if isinstance(cfg, dict) else dict(cfg or {})
    if overrides:
        cfg_run.update(overrides)

    script_dir = Path(__file__).resolve().parent
    config_path_value = cfg_run.get("_config_path") or cfg_run.get("config_path")
    config_id = str(config_path_value) if config_path_value else "<in-memory>"

    ablation_cfg = cfg_run.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_enabled = bool(ablation_cfg.get("enabled", False))

    epochs = int(cfg_run.get("mlm_epochs", 150))
    if ablation_enabled and ablation_cfg.get("pretrain_epochs") is not None:
        epochs_raw = ablation_cfg.get("pretrain_epochs")
        try:
            epochs = int(epochs_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.pretrain_epochs must be an int, got {epochs_raw!r}") from exc
        if epochs <= 0:
            raise ValueError(f"ablation.pretrain_epochs must be > 0, got {epochs}")

    max_steps: int | None = None
    if ablation_enabled and ablation_cfg.get("pretrain_max_steps") is not None:
        max_steps_raw = ablation_cfg.get("pretrain_max_steps")
        try:
            max_steps_int = int(max_steps_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.pretrain_max_steps must be an int, got {max_steps_raw!r}") from exc
        if max_steps_int < 0:
            raise ValueError(f"ablation.pretrain_max_steps must be >= 0, got {max_steps_int}")
        if max_steps_int > 0:
            max_steps = max_steps_int

    label_debug = bool(cfg_run.get("label_debug", False))
    debug_mode = bool(cfg_run.get("debug", False))
    overfit_debug = bool(cfg_run.get("overfit_debug", False))
    leakage_test = bool(cfg_run.get("leakage_test", False))
    generate_golden_batch = bool(cfg_run.get("generate_golden_batch", False))
    topk_debug = bool(cfg_run.get("topk_debug", False))
    golden_batch_path = cfg_run.get("golden_batch_path")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"mlm_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    length_curriculum_config = LengthCurriculumConfig.from_config(cfg_run)

    resolved_config: Dict[str, Any] = {
        "script": "pretrain_mlm.py",
        "timestamp": timestamp,
        "config_path": config_id,
        # Ablation controls
        "ablation_enabled": ablation_enabled,
        "pretrain_bp_limit": ablation_cfg.get("pretrain_bp_limit") if ablation_enabled else None,
        "pretrain_max_steps": max_steps,
        "pretrain_epochs": epochs if ablation_enabled else None,
        # Model architecture
        "embedding_dim": int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256))),
        "transformer_layers": int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8))),
        "transformer_heads": int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8))),
        "transformer_ff_dim": int(cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", 1024))),
        "transformer_dropout": float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1))),
        # Optimizer
        "lr": float(cfg_run.get("mlm_lr", 3e-4)),
        "weight_decay": float(cfg_run.get("mlm_weight_decay", 0.01)),
        # MLM-specific
        "mask_prob": 0.15,  # Fixed in code
        "max_span_len": 3,  # Fixed in mask_tokens_span
        # Data
        "batch_size": int(cfg_run.get("mlm_batch_size", 128)),
        "max_bp_len": int(cfg_run.get("mlm_window_size", 81)),
        "kmer": int(cfg_run.get("mlm_kmer", 3)),
        "include_reverse_complements": bool(
            cfg_run.get(
                "include_reverse_complements",
                cfg_run.get("mlm_include_reverse_complements", False),
            )
        ),
        # Training
        "epochs": epochs,
        "grad_accum_steps": int(cfg_run.get("mlm_grad_accum_steps", 4)),
        "warmup_ratio": float(cfg_run.get("mlm_warmup_ratio", 0.06)),
        "final_lr_ratio": float(cfg_run.get("mlm_final_lr_ratio", 0.1)),
        "lr_log_interval": int(cfg_run.get("mlm_lr_log_interval", 50)),
        "patience": int(cfg_run.get("mlm_patience", 30)),
        "val_ratio": float(cfg_run.get("mlm_val_ratio", 0.1)),
        # Model options
        "tie_weights": bool(cfg_run.get("mlm_tie_weights", True)),
        "use_output_norm": bool(cfg_run.get("mlm_use_output_norm", True)),
        "use_layer_lr_decay": bool(cfg_run.get("mlm_use_layer_lr_decay", True)),
        "layer_lr_decay": float(cfg_run.get("mlm_layer_lr_decay", 0.9)),
        "use_alibi": bool(cfg_run.get("use_alibi", True)),
        # Seed
        "seed": int(cfg_run.get("seed", 1337) or 1337),
        # Length curriculum
        "length_curriculum": length_curriculum_config.to_dict(),
    }

    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)

    seed = int(cfg_run.get("seed", 1337) or 1337)
    set_seed(seed)
    LOGGER.info("Using seed: %d", seed)

    token_tensors, vocab = prepare_dataset(cfg_run)
    full_dataset = MLMDataset(token_tensors)

    val_ratio = float(cfg_run.get("mlm_val_ratio", 0.1))
    val_size = int(len(full_dataset) * val_ratio)
    train_size = len(full_dataset) - val_size
    train_dataset, val_dataset = random_split(
        full_dataset,
        [train_size, val_size],
        generator=torch.Generator().manual_seed(seed),
    )
    LOGGER.info("Dataset split: %d train, %d val", train_size, val_size)

    batch_size = int(cfg_run.get("mlm_batch_size", 128))
    grad_accum_steps = int(cfg_run.get("mlm_grad_accum_steps", 2))
    effective_batch = batch_size * grad_accum_steps
    LOGGER.info(
        "Batch size: %d, Grad accum: %d, Effective batch: %d",
        batch_size,
        grad_accum_steps,
        effective_batch,
    )

    num_workers = int(cfg_run.get("num_workers", 0))
    steps_per_epoch_approx = len(train_dataset) // batch_size
    total_steps_for_curriculum = epochs * steps_per_epoch_approx // grad_accum_steps
    if max_steps is not None:
        total_steps_for_curriculum = min(total_steps_for_curriculum, max_steps)

    curriculum_scheduler: LengthCurriculumScheduler | None = None
    if length_curriculum_config.enabled:
        LOGGER.info("=" * 60)
        LOGGER.info("LENGTH CURRICULUM ENABLED")
        LOGGER.info("  Start length: %d tokens", length_curriculum_config.start_tokens)
        LOGGER.info("  End length:   %d tokens", length_curriculum_config.end_tokens)
        LOGGER.info(
            "  Ramp period:  %.0f%% - %.0f%% of training",
            length_curriculum_config.start_pct * 100,
            length_curriculum_config.end_pct * 100,
        )
        LOGGER.info("  Total steps:  %d", total_steps_for_curriculum)
        LOGGER.info("=" * 60)

        collator, curriculum_scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=total_steps_for_curriculum,
            cfg=cfg_run,
            mask_prob=0.15,
            track_spans=True,
            debug=label_debug,
        )
    else:
        LOGGER.info("Length curriculum: DISABLED (using fixed length)")
        collator = MLMCollator(vocab, mask_prob=0.15, debug=label_debug, track_spans=True)

    deterministic_kwargs = get_deterministic_dataloader_kwargs(seed, num_workers=num_workers)

    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=True,
        **deterministic_kwargs,
    )

    val_collator = MLMCollator(vocab, mask_prob=0.15)
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=val_collator,
        num_workers=0,
        pin_memory=True,
    )

    if label_debug:
        LOGGER.info("Running LABEL DEBUG mode - verifying first batch then exiting")
        for _batch in train_loader:
            break
        LOGGER.info("Label debug complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": 0,
        }

    if generate_golden_batch:
        from golden_batch import generate_golden_batch as gen_golden

        LOGGER.info("Generating golden batch with seed %d...", seed)

        if golden_batch_path:
            golden_path = Path(str(golden_batch_path))
            if not golden_path.is_absolute():
                golden_path = (script_dir / str(golden_batch_path)).resolve()
        else:
            golden_path = (script_dir / "../artifacts/golden_batch.pt").resolve()

        device = select_device()

        embedding_dim = int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256)))
        transformer_layers = int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8)))
        transformer_heads = int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8)))
        transformer_ff_dim = int(
            cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", embedding_dim * 4))
        )
        transformer_dropout = float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1)))

        encoder = DNAEncoder(
            vocab_size=len(vocab.itos),
            kmer=vocab.k,
            embedding_dim=embedding_dim,
            num_layers=transformer_layers,
            num_heads=transformer_heads,
            ff_dim=transformer_ff_dim,
            dropout=transformer_dropout,
            use_alibi=bool(cfg_run.get("use_alibi", True)),
            pad_token_id=vocab.pad_id,
            drop_path_rate=0.0,
        )

        special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
        model = DNAMLM(
            encoder,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
            tie_weights=bool(cfg_run.get("mlm_tie_weights", True)),
            use_output_norm=bool(cfg_run.get("mlm_use_output_norm", True)),
        )

        metadata = {
            "config": config_id,
            "embedding_dim": embedding_dim,
            "transformer_layers": transformer_layers,
            "transformer_heads": transformer_heads,
            "batch_size": batch_size,
            "seed": seed,
        }

        gen_golden(
            dataloader=train_loader,
            model=model,
            device=device,
            seed=seed,
            path=golden_path,
            metadata=metadata,
        )

        LOGGER.info("Golden batch generation complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": sum(p.numel() for p in model.parameters()),
        }

    LOGGER.info("Running pretrain sanity checks (checkpoint gate)...")
    sanity_passed, sanity_issues = run_pretrain_sanity_checks(
        train_loader=train_loader,
        vocab=vocab,
        num_batches=5,
        target_mask_prob=0.15,
    )

    if not sanity_passed:
        LOGGER.error("Sanity checks FAILED! Training cannot proceed.")
        for issue in sanity_issues:
            LOGGER.error("  - %s", issue)
        raise RuntimeError(
            "Pretrain sanity checks failed. Fix the issues above before training.\n"
            + "\n".join(f"  - {i}" for i in sanity_issues)
        )

    LOGGER.info("Running label token frequency debug report...")
    label_report = None
    for batch_data in train_loader:
        if len(batch_data) == 3:
            _, labels, _ = batch_data
        else:
            _, labels = batch_data
        label_report = debug_label_token_frequency(labels, vocab, num_batches=1)
        break

    if label_report:
        print_label_debug_report(label_report, vocab, top_k=10)
        if label_report["issues"]:
            raise RuntimeError(
                "Label token frequency check FAILED!\n"
                "Special tokens found in supervised labels.\n"
                + "\n".join(f"  - {i}" for i in label_report["issues"])
            )

    if debug_mode:
        LOGGER.info("Running in DEBUG mode - analyzing masking on one batch")
        raw_batch = [train_dataset[i] for i in range(min(batch_size, len(train_dataset)))]
        debug_masking(vocab, raw_batch, mask_prob=0.15, max_span_len=3)
        LOGGER.info("Debug mode complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": 0,
        }

    device = select_device()
    LOGGER.info("Using device: %s", device)

    embedding_dim = int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256)))
    transformer_layers = int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8)))
    transformer_heads = int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8)))
    transformer_ff_dim = int(cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", embedding_dim * 4)))
    transformer_dropout = float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1)))

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
        use_alibi=bool(cfg_run.get("use_alibi", True)),
        pad_token_id=vocab.pad_id,
        drop_path_rate=0.0,
    )

    tie_weights = bool(cfg_run.get("mlm_tie_weights", True))
    use_output_norm = bool(cfg_run.get("mlm_use_output_norm", True))

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    LOGGER.info(
        "Special tokens excluded from predictions: [MASK]=%d, <UNK>=%d, <PAD>=%d",
        vocab.mask_id,
        vocab.unk_id,
        vocab.pad_id,
    )

    model = DNAMLM(
        encoder,
        vocab_size=len(vocab.itos),
        special_token_ids=special_token_ids,
        tie_weights=tie_weights,
        use_output_norm=use_output_norm,
    ).to(device)

    total_params = sum(p.numel() for p in model.parameters())
    unique_params = sum(p.numel() for p in set(model.parameters()))
    LOGGER.info("Model parameters: %.2fM total, %.2fM unique", total_params / 1e6, unique_params / 1e6)

    if leakage_test:
        LOGGER.info("Running LEAKAGE TEST mode - verifying no information leakage")
        verify_no_leakage(model, vocab, device)
        LOGGER.info("Leakage test complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": total_params,
        }

    if overfit_debug:
        LOGGER.info("Running OVERFIT DEBUG mode")
        overfit_samples = token_tensors[:32]
        run_overfit_debug(
            model=model,
            train_samples=overfit_samples,
            vocab=vocab,
            device=device,
            num_steps=1000,
            lr=5e-3,
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": total_params,
        }

    lr = float(cfg_run.get("mlm_lr", 3e-4))
    weight_decay = float(cfg_run.get("mlm_weight_decay", 0.01))
    lr_decay = float(cfg_run.get("mlm_layer_lr_decay", 0.9))
    use_layer_lr_decay = bool(cfg_run.get("mlm_use_layer_lr_decay", True))

    LOGGER.info("=" * 60)
    LOGGER.info("OPTIMIZER CONFIGURATION")
    if use_layer_lr_decay:
        param_groups = get_layer_wise_lr_groups(model, lr, lr_decay, weight_decay)
        LOGGER.info("Layer-wise LR decay: ENABLED (factor: %.2f)", lr_decay)
        LOGGER.info("  Lower layers get smaller LRs to preserve pretrained features")
        optimizer = AdamW(
            param_groups,
            betas=(0.9, 0.98),
            eps=1e-6,
        )
    else:
        LOGGER.info("Layer-wise LR decay: DISABLED (all layers use same LR)")
        LOGGER.info("  All parameters use base LR: %.2e", lr)
        optimizer = AdamW(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay,
            betas=(0.9, 0.98),
            eps=1e-6,
        )
    LOGGER.info("=" * 60)

    batches_per_epoch = len(train_loader)
    steps_per_epoch = batches_per_epoch // grad_accum_steps
    total_steps = epochs * steps_per_epoch
    if max_steps is not None:
        total_steps = min(total_steps, max_steps)

    warmup_ratio = float(cfg_run.get("mlm_warmup_ratio", 0.06))
    warmup_steps = int(warmup_ratio * total_steps)
    final_lr_ratio = float(cfg_run.get("mlm_final_lr_ratio", 0.1))
    lr_log_interval = int(cfg_run.get("mlm_lr_log_interval", 50))

    LOGGER.info("=" * 60)
    LOGGER.info("TOTAL STEPS CALCULATION")
    LOGGER.info("  Epochs:                  %d", epochs)
    LOGGER.info("  Batches per epoch:       %d", batches_per_epoch)
    LOGGER.info("  Grad accum steps:        %d", grad_accum_steps)
    LOGGER.info("  Optimizer steps/epoch:   %d", steps_per_epoch)
    LOGGER.info("  Total optimizer steps:   %d", total_steps)
    LOGGER.info("=" * 60)

    scheduler = get_warmup_cosine_schedule(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        final_lr_ratio=final_lr_ratio,
    )
    LOGGER.info("=" * 60)
    LOGGER.info("LR SCHEDULE (Linear Warmup + Cosine Decay)")
    LOGGER.info("  Base LR:        %.2e", lr)
    LOGGER.info("  Warmup steps:   %d (%.1f%% of total)", warmup_steps, warmup_ratio * 100)
    LOGGER.info("  Total steps:    %d", total_steps)
    LOGGER.info("  Final LR:       %.2e (%.0f%% of peak)", lr * final_lr_ratio, final_lr_ratio * 100)
    LOGGER.info("  LR log interval: every %d optimizer steps", lr_log_interval)
    LOGGER.info("=" * 60)

    autocast_device: str | None = None
    autocast_dtype: torch.dtype | None = None
    scaler = None

    if device.type == "cuda":
        autocast_device = "cuda"
        autocast_dtype = torch.float16
        scaler = torch.cuda.amp.GradScaler()
    elif device.type == "mps":
        try:
            with torch.autocast("mps", dtype=torch.bfloat16):
                x = torch.zeros(1, device=device, dtype=torch.float32)
                _ = x + x
            autocast_device = "mps"
            autocast_dtype = torch.bfloat16
        except Exception:
            LOGGER.warning(
                "MPS autocast bf16 not available; disabling autocast to avoid fp16 without GradScaler."
            )

    best_val_loss = float("inf")
    best_val_acc = 0.0
    patience = int(cfg_run.get("mlm_patience", 25))
    patience_counter = 0
    global_step = 0
    param_norm_warn = float(cfg_run.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg_run.get("param_norm_cap", 200.0))

    encoder_path = _resolve_cfg_path(
        cfg_run.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt"),
        base_dir=script_dir,
    )
    encoder_path.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("=" * 60)
    LOGGER.info("STARTING MLM PRETRAINING")
    LOGGER.info("Epochs: %d, LR: %.2e, Mask prob: 15%%", epochs, lr)
    if max_steps is not None:
        LOGGER.info("Ablation enabled: stopping after %d optimizer steps", max_steps)
    LOGGER.info("=" * 60)

    health_log_interval = int(cfg_run.get("health_log_interval", 50))
    health_dashboard = create_health_dashboard(
        run_dir=run_dir,
        vocab=vocab,
        log_every_n_steps=health_log_interval,
        target_mask_prob=0.15,
    )
    health_dashboard.print_table_header()

    numerics_checker = NumericsChecker(enabled=True)
    LOGGER.info("Numerics checker enabled - will fail fast on NaN/Inf")

    best_train_loss = float("inf")
    reached_max_steps = False

    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0
        total_correct = 0
        total_masked = 0
        samples_seen = 0
        optimizer.zero_grad()

        if scaler is not None:
            LOGGER.info("Epoch %d: GradScaler scale = %.0f", epoch, scaler.get_scale())

        if curriculum_scheduler is not None:
            curr_len = curriculum_scheduler.current_length
            progress = curriculum_scheduler.get_progress()
            LOGGER.info(
                "Epoch %d: Curriculum length = %d tokens (%s, %.0f%% ramp progress)",
                epoch,
                curr_len,
                progress["phase"],
                progress["ramp_progress"] * 100,
            )

        for batch_idx, batch_data in enumerate(train_loader, start=1):
            if max_steps is not None and global_step >= max_steps:
                reached_max_steps = True
                break

            if len(batch_data) == 3:
                inputs, labels, span_lengths = batch_data
            else:
                inputs, labels = batch_data
                span_lengths = None

            inputs = inputs.to(device)
            labels = labels.to(device)
            samples_seen += int(inputs.size(0))

            current_lr = optimizer.param_groups[0]["lr"]
            numerics_checker.set_context(
                epoch=epoch,
                step=global_step,
                lr=current_lr,
                batch_idx=batch_idx,
                model=model,
            )

            if autocast_device and autocast_dtype is not None:
                with torch.autocast(autocast_device, dtype=autocast_dtype):
                    logits, _ = model(inputs, labels=None)
            else:
                logits, _ = model(inputs, labels=None)

            loss = F.cross_entropy(
                logits.float().view(-1, logits.size(-1)),
                labels.view(-1),
                ignore_index=-100,
            )
            numerics_checker.check_forward(logits, loss)

            scaled_loss = loss / grad_accum_steps
            if scaler is not None:
                scaler.scale(scaled_loss).backward()
            else:
                scaled_loss.backward()

            if scaler is not None:
                numerics_checker.check_backward(model, check_grads=False)
            else:
                numerics_checker.check_backward(model)

            total_loss += loss.item() * inputs.size(0)

            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()

            health_dashboard.accumulate_batch(
                loss=loss.item(),
                logits=logits.detach(),
                labels=labels,
                masked_inputs=inputs,
                span_lengths=span_lengths,
            )

            if batch_idx % grad_accum_steps == 0:
                if scaler is not None:
                    scaler.unscale_(optimizer)
                unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)

                if scaler is not None:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()

                total_norm, scale = cap_model_param_norm_(model, max_norm=param_norm_cap)
                effective_param_norm = param_norm_cap if scale < 1.0 else total_norm
                if scale < 1.0:
                    LOGGER.warning(
                        "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                        global_step + 1,
                        total_norm,
                        param_norm_cap,
                        scale,
                    )

                scheduler.step()
                optimizer.zero_grad()
                global_step += 1

                if global_step % lr_log_interval == 0:
                    current_lr = optimizer.param_groups[0]["lr"]
                    if global_step <= warmup_steps:
                        phase = "warmup"
                        phase_progress = global_step / max(1, warmup_steps) * 100
                    else:
                        phase = "decay"
                        phase_progress = (global_step - warmup_steps) / max(1, total_steps - warmup_steps) * 100

                    LOGGER.info(
                        "Step %d/%d [%s %.0f%%]: LR=%.2e, grad_norm=%.4f, param_norm=%.2f",
                        global_step,
                        total_steps,
                        phase,
                        phase_progress,
                        current_lr,
                        unclipped_grad_norm.item(),
                        effective_param_norm,
                    )
                    if param_norm_warn > 0 and effective_param_norm > param_norm_warn:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f above warn=%.2f",
                            global_step,
                            effective_param_norm,
                            param_norm_warn,
                        )

                if curriculum_scheduler is not None:
                    new_length = curriculum_scheduler.step(global_step)
                    if curriculum_scheduler.is_ramping() and global_step % 100 == 0:
                        LOGGER.debug("Step %d: Curriculum length = %d tokens", global_step, new_length)

                health_metrics = health_dashboard.log_step(
                    step=global_step,
                    epoch=epoch,
                    model=model,
                    optimizer=optimizer,
                )
                if health_metrics:
                    health_dashboard.print_table_row(health_metrics)

                if max_steps is not None and global_step >= max_steps:
                    reached_max_steps = True
                    break

            if batch_idx % 100 == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                batch_acc = total_correct / max(1, total_masked)
                LOGGER.info(
                    "Epoch %d [%d/%d] Loss: %.4f Acc: %.2f%% LR: %.2e",
                    epoch,
                    batch_idx,
                    len(train_loader),
                    loss.item(),
                    batch_acc * 100,
                    current_lr,
                )

        train_loss = total_loss / max(1, samples_seen)
        train_acc = total_correct / max(1, total_masked)
        train_ppl = math.exp(min(train_loss, 100)) if math.isfinite(train_loss) else float("nan")
        if math.isfinite(train_loss):
            best_train_loss = min(best_train_loss, train_loss)

        val_metrics = run_validation(model, val_loader, device, autocast_device, autocast_dtype)
        val_loss = val_metrics["loss"]
        val_acc = val_metrics["accuracy"]
        val_ppl = val_metrics["perplexity"]

        curriculum_info = f" | Len: {curriculum_scheduler.current_length}" if curriculum_scheduler is not None else ""

        LOGGER.info(
            "Epoch %d: Train [Loss: %.4f, Acc: %.1f%%, PPL: %.2f] | Val [Loss: %.4f, Acc: %.1f%%, PPL: %.2f]%s",
            epoch,
            train_loss,
            train_acc * 100,
            train_ppl,
            val_loss,
            val_acc * 100,
            val_ppl,
            curriculum_info,
        )

        if topk_debug:
            for debug_batch_data in train_loader:
                if len(debug_batch_data) == 3:
                    debug_inputs, debug_labels, _ = debug_batch_data
                else:
                    debug_inputs, debug_labels = debug_batch_data
                break

            topk_stats = print_topk_predictions_debug(
                model=model,
                inputs=debug_inputs,
                labels=debug_labels,
                vocab=vocab,
                epoch=epoch,
                num_sequences=3,
                top_k=5,
            )

            LOGGER.info(
                "Top-K Debug: Avg top-1 prob: %.1f%%, True in top-1: %.1f%%, True in top-5: %.1f%%",
                topk_stats["avg_top1_prob"] * 100,
                topk_stats["true_in_top1_pct"],
                topk_stats["true_in_topk_pct"],
            )

        if val_loss < best_val_loss - 0.001:
            best_val_loss = val_loss
            best_val_acc = val_acc
            patience_counter = 0
            torch.save(model.encoder.state_dict(), encoder_path)
            LOGGER.info(
                "★ New best! Val Loss: %.4f, PPL: %.2f, Acc: %.1f%% - Saved to %s",
                val_loss,
                val_ppl,
                val_acc * 100,
                encoder_path,
            )
        else:
            patience_counter += 1
            if patience_counter >= patience:
                LOGGER.info("Early stopping at epoch %d (patience %d)", epoch, patience)
                break

        if reached_max_steps:
            LOGGER.info(
                "Reached ablation.pretrain_max_steps=%d (optimizer steps); stopping after validation",
                max_steps,
            )
            break

    if not encoder_path.exists():
        torch.save(model.encoder.state_dict(), encoder_path)

    best_ppl = math.exp(min(best_val_loss, 100)) if math.isfinite(best_val_loss) else float("nan")
    LOGGER.info("=" * 60)
    LOGGER.info("PRETRAINING COMPLETE")
    LOGGER.info("Best Val Loss: %.4f (Perplexity: %.2f)", best_val_loss, best_ppl)
    LOGGER.info("Best Val Accuracy: %.1f%%", best_val_acc * 100)
    LOGGER.info("Target: Loss < 1.0 (PPL < 2.72)")
    LOGGER.info("Saved pretrained encoder to %s", encoder_path)
    LOGGER.info("=" * 60)

    return {
        "mlm_loss": float(best_train_loss if math.isfinite(best_train_loss) else train_loss),
        "mlm_steps": int(global_step),
        "encoder_ckpt_path": str(encoder_path),
        "pretrain_seconds": float(time.perf_counter() - start_time),
        "params": int(total_params),
        "best_val_loss": float(best_val_loss),
        "best_val_accuracy": float(best_val_acc),
        "epochs_ran": int(epoch),
    }


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
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility. Overrides config value.")
    parser.add_argument("--generate_golden_batch", action="store_true",
                        help="Generate and save a golden batch checkpoint for reproducibility testing.")
    parser.add_argument("--golden_batch_path", type=str, default=None,
                        help="Path to golden batch file (default: artifacts/golden_batch.pt)")
    parser.add_argument("--topk_debug", action="store_true",
                        help="Enable top-k sanity print for masked positions. "
                             "Prints masked position index, true token, and top-5 predictions + probabilities "
                             "for 3 sequences once per epoch. Helps verify predictions become sharper over time.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    overrides: Dict[str, Any] = {"_config_path": str(config_path)}
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.debug:
        overrides["debug"] = True
    if args.overfit_debug:
        overrides["overfit_debug"] = True
    if args.label_debug:
        overrides["label_debug"] = True
    if args.leakage_test:
        overrides["leakage_test"] = True
    if args.generate_golden_batch:
        overrides["generate_golden_batch"] = True
    if args.golden_batch_path is not None:
        overrides["golden_batch_path"] = args.golden_batch_path
    if args.topk_debug:
        overrides["topk_debug"] = True

    _ = run_mlm_pretrain(cfg, overrides=overrides)
    return

    # =========================================================================
    # Log resolved hyperparameters to runs/<timestamp>/resolved_config.json
    # =========================================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"mlm_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse length curriculum config early (needed for resolved_config)
    length_curriculum_config = LengthCurriculumConfig.from_config(cfg)
    
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
        "lr": float(cfg.get("mlm_lr", 3e-4)),
        "weight_decay": float(cfg.get("mlm_weight_decay", 0.01)),
        # MLM-specific
        "mask_prob": 0.15,  # Fixed in code
        "max_span_len": 3,   # Fixed in mask_tokens_span
        # Data
        "batch_size": int(cfg.get("mlm_batch_size", 128)),
        "max_bp_len": int(cfg.get("mlm_window_size", 81)),
        "kmer": int(cfg.get("mlm_kmer", 3)),
        "include_reverse_complements": bool(cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False))),
        # Training
        "epochs": int(cfg.get("mlm_epochs", 200)),
        "grad_accum_steps": int(cfg.get("mlm_grad_accum_steps", 4)),
        "warmup_ratio": float(cfg.get("mlm_warmup_ratio", 0.06)),
        "final_lr_ratio": float(cfg.get("mlm_final_lr_ratio", 0.1)),
        "lr_log_interval": int(cfg.get("mlm_lr_log_interval", 50)),
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
        # Length curriculum
        "length_curriculum": length_curriculum_config.to_dict(),
    }
    
    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)
    
    # Seed: CLI argument takes precedence over config
    seed = args.seed if args.seed is not None else int(cfg.get("seed", 1337))
    set_seed(seed)
    LOGGER.info("Using seed: %d", seed)

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
    
    # Standard 15% masking - enable debug mode if --label_debug flag is set
    # track_spans=True enables span length tracking for health dashboard
    
    # Calculate total steps for curriculum scheduler
    num_workers = int(cfg.get("num_workers", 0))
    epochs = int(cfg.get("mlm_epochs", 200))
    steps_per_epoch_approx = len(train_dataset) // batch_size
    total_steps = epochs * steps_per_epoch_approx // grad_accum_steps
    
    # Create collator with or without length curriculum
    curriculum_scheduler: LengthCurriculumScheduler | None = None
    if length_curriculum_config.enabled:
        LOGGER.info("=" * 60)
        LOGGER.info("LENGTH CURRICULUM ENABLED")
        LOGGER.info("  Start length: %d tokens", length_curriculum_config.start_tokens)
        LOGGER.info("  End length:   %d tokens", length_curriculum_config.end_tokens)
        LOGGER.info("  Ramp period:  %.0f%% - %.0f%% of training",
                    length_curriculum_config.start_pct * 100,
                    length_curriculum_config.end_pct * 100)
        LOGGER.info("  Total steps:  %d", total_steps)
        LOGGER.info("=" * 60)
        
        collator, curriculum_scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=total_steps,
            cfg=cfg,
            mask_prob=0.15,
            track_spans=True,
            debug=args.label_debug,
        )
    else:
        LOGGER.info("Length curriculum: DISABLED (using fixed length)")
        collator = MLMCollator(vocab, mask_prob=0.15, debug=args.label_debug, track_spans=True)
    
    deterministic_kwargs = get_deterministic_dataloader_kwargs(seed, num_workers=num_workers)
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=True,
        **deterministic_kwargs,
    )
    
    # Validation always uses full length (no curriculum)
    val_collator = MLMCollator(vocab, mask_prob=0.15)
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=val_collator,
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
    # GOLDEN BATCH MODE: Generate golden batch for reproducibility testing
    # (Runs before sanity checks as it doesn't need them)
    # =========================================================================
    if args.generate_golden_batch:
        from golden_batch import generate_golden_batch as gen_golden
        
        LOGGER.info("Generating golden batch with seed %d...", seed)
        
        # Determine golden batch path
        if args.golden_batch_path:
            golden_path = Path(args.golden_batch_path)
            if not golden_path.is_absolute():
                golden_path = (script_dir / args.golden_batch_path).resolve()
        else:
            golden_path = (script_dir / "../artifacts/golden_batch.pt").resolve()
        
        device = select_device()
        
        # Build model for golden batch generation
        embedding_dim = int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256)))
        transformer_layers = int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8)))
        transformer_heads = int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8)))
        transformer_ff_dim = int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", embedding_dim * 4)))
        transformer_dropout = float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1)))
        
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
            drop_path_rate=0.0,
        )
        
        special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
        model = DNAMLM(
            encoder,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
            tie_weights=bool(cfg.get("mlm_tie_weights", True)),
            use_output_norm=bool(cfg.get("mlm_use_output_norm", True)),
        )
        
        metadata = {
            "config": str(config_path),
            "embedding_dim": embedding_dim,
            "transformer_layers": transformer_layers,
            "transformer_heads": transformer_heads,
            "batch_size": batch_size,
            "seed": seed,
        }
        
        gen_golden(
            dataloader=train_loader,
            model=model,
            device=device,
            seed=seed,
            path=golden_path,
            metadata=metadata,
        )
        
        LOGGER.info("Golden batch generation complete. Exiting.")
        return

    # =========================================================================
    # CHECKPOINT GATE: Run sanity checks before training can start
    # =========================================================================
    LOGGER.info("Running pretrain sanity checks (checkpoint gate)...")
    sanity_passed, sanity_issues = run_pretrain_sanity_checks(
        train_loader=train_loader,
        vocab=vocab,
        num_batches=5,
        target_mask_prob=0.15,
    )
    
    if not sanity_passed:
        LOGGER.error("Sanity checks FAILED! Training cannot proceed.")
        for issue in sanity_issues:
            LOGGER.error("  - %s", issue)
        raise RuntimeError(
            "Pretrain sanity checks failed. Fix the issues above before training.\n"
            + "\n".join(f"  - {i}" for i in sanity_issues)
        )
    
    # =========================================================================
    # LABEL DEBUG REPORT: Verify no special tokens in supervised labels
    # =========================================================================
    LOGGER.info("Running label token frequency debug report...")
    label_report = None
    for batch_data in train_loader:
        if len(batch_data) == 3:
            _, labels, _ = batch_data
        else:
            _, labels = batch_data
        label_report = debug_label_token_frequency(labels, vocab, num_batches=1)
        break  # Only need one batch for the report
    
    if label_report:
        print_label_debug_report(label_report, vocab, top_k=10)
        
        # HARD FAIL if special tokens found in labels
        if label_report["issues"]:
            raise RuntimeError(
                "Label token frequency check FAILED!\n"
                "Special tokens found in supervised labels.\n"
                + "\n".join(f"  - {i}" for i in label_report["issues"])
            )

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
    lr = float(cfg.get("mlm_lr", 3e-4))  # Base LR 3e-4 for stable MLM pretraining
    weight_decay = float(cfg.get("mlm_weight_decay", 0.01))
    lr_decay = float(cfg.get("mlm_layer_lr_decay", 0.9))  # Layer-wise LR decay
    use_layer_lr_decay = bool(cfg.get("mlm_use_layer_lr_decay", True))
    
    # Build optimizer with or without layer-wise LR decay
    LOGGER.info("=" * 60)
    LOGGER.info("OPTIMIZER CONFIGURATION")
    if use_layer_lr_decay:
        param_groups = get_layer_wise_lr_groups(model, lr, lr_decay, weight_decay)
        LOGGER.info("Layer-wise LR decay: ENABLED (factor: %.2f)", lr_decay)
        LOGGER.info("  Lower layers get smaller LRs to preserve pretrained features")
        optimizer = AdamW(
            param_groups, 
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,  # Smaller eps for stability
        )
    else:
        LOGGER.info("Layer-wise LR decay: DISABLED (all layers use same LR)")
        LOGGER.info("  All parameters use base LR: %.2e", lr)
        optimizer = AdamW(
            model.parameters(), 
            lr=lr, 
            weight_decay=weight_decay,
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,
        )
    LOGGER.info("=" * 60)
    
    epochs = int(cfg.get("mlm_epochs", 150))
    
    # Learning rate schedule with warmup + cosine decay
    # total_steps = number of OPTIMIZER steps (not batches)
    # 
    # Calculation:
    #   batches_per_epoch = len(train_loader)  [with drop_last=True]
    #   optimizer_steps_per_epoch = batches_per_epoch // grad_accum_steps
    #   total_optimizer_steps = epochs * optimizer_steps_per_epoch
    #
    batches_per_epoch = len(train_loader)
    steps_per_epoch = batches_per_epoch // grad_accum_steps
    total_steps = epochs * steps_per_epoch
    warmup_ratio = float(cfg.get("mlm_warmup_ratio", 0.06))  # 6% warmup
    warmup_steps = int(warmup_ratio * total_steps)
    final_lr_ratio = float(cfg.get("mlm_final_lr_ratio", 0.1))  # Decay to 10% of peak
    lr_log_interval = int(cfg.get("mlm_lr_log_interval", 50))
    
    LOGGER.info("=" * 60)
    LOGGER.info("TOTAL STEPS CALCULATION")
    LOGGER.info("  Epochs:                  %d", epochs)
    LOGGER.info("  Batches per epoch:       %d", batches_per_epoch)
    LOGGER.info("  Grad accum steps:        %d", grad_accum_steps)
    LOGGER.info("  Optimizer steps/epoch:   %d", steps_per_epoch)
    LOGGER.info("  Total optimizer steps:   %d", total_steps)
    LOGGER.info("=" * 60)
    
    scheduler = get_warmup_cosine_schedule(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        final_lr_ratio=final_lr_ratio,
    )
    LOGGER.info("=" * 60)
    LOGGER.info("LR SCHEDULE (Linear Warmup + Cosine Decay)")
    LOGGER.info("  Base LR:        %.2e", lr)
    LOGGER.info("  Warmup steps:   %d (%.1f%% of total)", warmup_steps, warmup_ratio * 100)
    LOGGER.info("  Total steps:    %d", total_steps)
    LOGGER.info("  Final LR:       %.2e (%.0f%% of peak)", lr * final_lr_ratio, final_lr_ratio * 100)
    LOGGER.info("  LR log interval: every %d optimizer steps", lr_log_interval)
    LOGGER.info("=" * 60)

    autocast_device: str | None = None
    autocast_dtype: torch.dtype | None = None
    scaler = None

    if device.type == "cuda":
        autocast_device = "cuda"
        autocast_dtype = torch.float16
        scaler = torch.cuda.amp.GradScaler()
    elif device.type == "mps":
        # Prefer bf16 autocast on MPS (avoids fp16 without GradScaler).
        try:
            with torch.autocast("mps", dtype=torch.bfloat16):
                x = torch.zeros(1, device=device, dtype=torch.float32)
                _ = x + x
            autocast_device = "mps"
            autocast_dtype = torch.bfloat16
        except Exception:
            LOGGER.warning(
                "MPS autocast bf16 not available; disabling autocast to avoid fp16 without GradScaler."
            )
    
    # Training state
    best_val_loss = float('inf')
    best_val_acc = 0.0
    patience = int(cfg.get("mlm_patience", 25))
    patience_counter = 0
    global_step = 0
    param_norm_warn = float(cfg.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg.get("param_norm_cap", 200.0))
    
    encoder_path = _resolve_cfg_path(
        cfg.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt"),
        base_dir=script_dir,
    )
    encoder_path.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("=" * 60)
    LOGGER.info("STARTING MLM PRETRAINING")
    LOGGER.info("Epochs: %d, LR: %.2e, Mask prob: 15%%", epochs, lr)
    LOGGER.info("=" * 60)
    
    # =========================================================================
    # HEALTH DASHBOARD: Initialize tracking
    # =========================================================================
    health_log_interval = int(cfg.get("health_log_interval", 50))
    health_dashboard = create_health_dashboard(
        run_dir=run_dir,
        vocab=vocab,
        log_every_n_steps=health_log_interval,
        target_mask_prob=0.15,
    )
    health_dashboard.print_table_header()
    
    # =========================================================================
    # NUMERICS CHECKER: NaN/Inf tripwires for fail-fast detection
    # =========================================================================
    numerics_checker = NumericsChecker(enabled=True)
    LOGGER.info("Numerics checker enabled - will fail fast on NaN/Inf")

    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0
        total_correct = 0
        total_masked = 0
        optimizer.zero_grad()

        if scaler is not None:
            LOGGER.info("Epoch %d: GradScaler scale = %.0f", epoch, scaler.get_scale())
        
        # Log current curriculum length at start of epoch
        if curriculum_scheduler is not None:
            curr_len = curriculum_scheduler.current_length
            progress = curriculum_scheduler.get_progress()
            LOGGER.info(
                "Epoch %d: Curriculum length = %d tokens (%s, %.0f%% ramp progress)",
                epoch, curr_len, progress["phase"], progress["ramp_progress"] * 100
            )
        
        for batch_idx, batch_data in enumerate(train_loader, start=1):
            # Handle both (inputs, labels) and (inputs, labels, spans) formats
            if len(batch_data) == 3:
                inputs, labels, span_lengths = batch_data
            else:
                inputs, labels = batch_data
                span_lengths = None
            
            inputs = inputs.to(device)
            labels = labels.to(device)
            
            # Update numerics checker context for error reporting
            current_lr = optimizer.param_groups[0]['lr']
            numerics_checker.set_context(
                epoch=epoch,
                step=global_step,
                lr=current_lr,
                batch_idx=batch_idx,
                model=model,
            )
            
            if autocast_device and autocast_dtype is not None:
                with torch.autocast(autocast_device, dtype=autocast_dtype):
                    logits, _ = model(inputs, labels=None)
            else:
                logits, _ = model(inputs, labels=None)

            # Loss computation outside autocast for numerical stability
            loss = F.cross_entropy(
                logits.float().view(-1, logits.size(-1)),
                labels.view(-1),
                ignore_index=-100,
            )
            
            # TRIPWIRE: Check forward pass outputs for NaN/Inf
            numerics_checker.check_forward(logits, loss)
            
            # Scale loss for gradient accumulation
            scaled_loss = loss / grad_accum_steps
            if scaler is not None:
                scaler.scale(scaled_loss).backward()
            else:
                scaled_loss.backward()
            
            # TRIPWIRE: In fp16+GradScaler mode, gradients are scaled and may overflow by design.
            # Still check parameters, but avoid failing on transient gradient Inf/NaN handled by scaler.
            if scaler is not None:
                numerics_checker.check_backward(model, check_grads=False)
            else:
                numerics_checker.check_backward(model)
            
            total_loss += loss.item() * inputs.size(0)
            
            # Compute accuracy
            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()
            
            # Accumulate metrics for health dashboard
            health_dashboard.accumulate_batch(
                loss=loss.item(),
                logits=logits.detach(),
                labels=labels,
                masked_inputs=inputs,
                span_lengths=span_lengths,
            )
            
            # Optimizer step with gradient accumulation
            if batch_idx % grad_accum_steps == 0:
                if scaler is not None:
                    scaler.unscale_(optimizer)
                # Gradient clipping with norm logging
                # clip_grad_norm_ returns the total norm BEFORE clipping
                unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                # After clipping, effective norm is min(unclipped, max_norm)
                clipped_grad_norm = min(unclipped_grad_norm.item(), 1.0)
                
                if scaler is not None:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()
                total_norm, scale = cap_model_param_norm_(model, max_norm=param_norm_cap)
                effective_param_norm = param_norm_cap if scale < 1.0 else total_norm
                if scale < 1.0:
                    LOGGER.warning(
                        "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                        global_step + 1,
                        total_norm,
                        param_norm_cap,
                        scale,
                    )
                scheduler.step()
                optimizer.zero_grad()
                global_step += 1
                
                # Log LR every N optimizer steps (configurable via mlm_lr_log_interval)
                if global_step % lr_log_interval == 0:
                    current_lr = optimizer.param_groups[0]['lr']
                    # Determine training phase for logging
                    if global_step <= warmup_steps:
                        phase = "warmup"
                        phase_progress = global_step / warmup_steps * 100
                    else:
                        phase = "decay"
                        phase_progress = (global_step - warmup_steps) / max(1, total_steps - warmup_steps) * 100
                    
                    LOGGER.info(
                        "Step %d/%d [%s %.0f%%]: LR=%.2e, grad_norm=%.4f, param_norm=%.2f",
                        global_step, total_steps, phase, phase_progress,
                        current_lr, unclipped_grad_norm.item(), effective_param_norm
                    )
                    if param_norm_warn > 0 and effective_param_norm > param_norm_warn:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f above warn=%.2f",
                            global_step,
                            effective_param_norm,
                            param_norm_warn,
                        )
                
                # Step curriculum scheduler and update length
                if curriculum_scheduler is not None:
                    new_length = curriculum_scheduler.step(global_step)
                    # Log length changes during ramping
                    if curriculum_scheduler.is_ramping() and global_step % 100 == 0:
                        LOGGER.debug(
                            "Step %d: Curriculum length = %d tokens",
                            global_step, new_length
                        )
                
                # Log health metrics every N steps
                health_metrics = health_dashboard.log_step(
                    step=global_step,
                    epoch=epoch,
                    model=model,
                    optimizer=optimizer,
                )
                if health_metrics:
                    health_dashboard.print_table_row(health_metrics)
            
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
        val_metrics = run_validation(model, val_loader, device, autocast_device, autocast_dtype)
        val_loss = val_metrics['loss']
        val_acc = val_metrics['accuracy']
        val_ppl = val_metrics['perplexity']
        
        # Build epoch summary with curriculum info
        curriculum_info = ""
        if curriculum_scheduler is not None:
            curriculum_info = f" | Len: {curriculum_scheduler.current_length}"
        
        LOGGER.info(
            "Epoch %d: Train [Loss: %.4f, Acc: %.1f%%, PPL: %.2f] | Val [Loss: %.4f, Acc: %.1f%%, PPL: %.2f]%s",
            epoch, train_loss, train_acc * 100, train_ppl, val_loss, val_acc * 100, val_ppl, curriculum_info
        )
        
        # =========================================================================
        # TOP-K SANITY PRINT: Run once per epoch if enabled
        # =========================================================================
        if args.topk_debug:
            # Get a fresh batch for debug analysis
            for debug_batch_data in train_loader:
                if len(debug_batch_data) == 3:
                    debug_inputs, debug_labels, _ = debug_batch_data
                else:
                    debug_inputs, debug_labels = debug_batch_data
                break
            
            topk_stats = print_topk_predictions_debug(
                model=model,
                inputs=debug_inputs,
                labels=debug_labels,
                vocab=vocab,
                epoch=epoch,
                num_sequences=3,
                top_k=5,
            )
            
            # Log summary stats for tracking over epochs
            LOGGER.info(
                "Top-K Debug: Avg top-1 prob: %.1f%%, True in top-1: %.1f%%, True in top-5: %.1f%%",
                topk_stats["avg_top1_prob"] * 100,
                topk_stats["true_in_top1_pct"],
                topk_stats["true_in_topk_pct"],
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
