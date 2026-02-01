#!/usr/bin/env python3
"""
Quick test with REAL genome data to see learning dynamics.
"""

import json
import logging
import math
import random
import sys
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR

# Setup path
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from bpe_tokenizer import DNABPETokenizer
from model import DNAEncoder
from pretrain_components.mlm_model import DNAMLM
from pretrain_mlm import mask_tokens_span
from genome_io import read_fasta_records

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
LOGGER = logging.getLogger(__name__)


def load_vocab(vocab_path: Path) -> DNABPETokenizer:
    """Load BPE vocabulary."""
    vocab = DNABPETokenizer()
    with open(vocab_path, "r") as f:
        vocab_data = json.load(f)
    vocab.vocab = vocab_data["vocab"]
    vocab.id_to_token = {v: k for k, v in vocab_data["vocab"].items()}
    vocab.merges = [tuple(m) for m in vocab_data["merges"]]
    return vocab


def extract_windows(sequence: str, window_size: int = 234, stride: int = 50) -> list:
    """Extract sliding windows from a sequence."""
    windows = []
    seq = sequence.upper()
    # Only keep windows with valid DNA characters
    for i in range(0, len(seq) - window_size + 1, stride):
        window = seq[i:i + window_size]
        if all(c in 'ACGT' for c in window):
            windows.append(window)
    return windows


def main():
    """Run quick training test with real genome data."""
    LOGGER.info("=" * 60)
    LOGGER.info("REAL GENOME TRAINING TEST")
    LOGGER.info("=" * 60)

    # Setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LOGGER.info("Device: %s", device)

    # Load vocabulary
    vocab_path = SCRIPT_DIR.parent / "artifacts" / "vocabs" / "bpe_vocab.json"
    vocab = load_vocab(vocab_path)
    LOGGER.info("Loaded vocabulary with %d tokens", len(vocab.vocab))

    # Load E. coli genome
    fasta_path = SCRIPT_DIR.parent / "data" / "GCF_000005845.2_ASM584v2_genomic.fna"
    if not fasta_path.exists():
        LOGGER.error("Genome file not found: %s", fasta_path)
        return

    LOGGER.info("Loading genome from %s", fasta_path)
    records = read_fasta_records(fasta_path, keep_header=True)
    # Records are (header, sequence) tuples
    genome = "".join(seq for _, seq in records)
    LOGGER.info("Genome length: %d bp", len(genome))

    # Extract windows
    windows = extract_windows(genome, window_size=234, stride=50)
    random.shuffle(windows)
    LOGGER.info("Extracted %d windows", len(windows))

    # Tokenize windows
    max_len = 53
    tokenized = []
    for window in windows[:10000]:  # Limit for speed
        tokens = vocab.tokenize(window)
        if len(tokens) < max_len:
            tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
        else:
            tokens = tokens[:max_len]
        tokenized.append(tokens)

    LOGGER.info("Tokenized %d sequences", len(tokenized))

    # Create model - use same settings as actual training
    embedding_dim = 384
    transformer_layers = 6  # Smaller for speed
    transformer_heads = 12
    vocab_size = len(vocab.vocab)

    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=embedding_dim * 4,
        dropout=0.1,
        pad_token_id=vocab.pad_id,
        use_rope=False,  # BERT-style
        ffn_type="gelu",
        norm_type="layernorm",
    )

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=vocab_size,
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,
    ).to(device)

    total_params = sum(p.numel() for p in model.parameters())
    LOGGER.info("Model parameters: %.2fM", total_params / 1e6)

    # Training settings
    batch_size = 64
    n_epochs = 10
    n_batches = len(tokenized) // batch_size

    # Test different learning rates
    learning_rates = [5e-4, 1e-3, 2e-3]

    for lr in learning_rates:
        LOGGER.info("\n" + "=" * 60)
        LOGGER.info("Testing LR = %.0e", lr)
        LOGGER.info("=" * 60)

        # Reset model weights
        def reset_weights(m):
            if hasattr(m, 'reset_parameters'):
                m.reset_parameters()
        model.apply(reset_weights)

        optimizer = AdamW(model.parameters(), lr=lr, weight_decay=0.01, betas=(0.9, 0.999))

        # Warmup + cosine schedule
        warmup_steps = n_batches  # 1 epoch warmup
        total_steps = n_epochs * n_batches

        def lr_lambda(step):
            if step < warmup_steps:
                return step / warmup_steps
            progress = (step - warmup_steps) / (total_steps - warmup_steps)
            return 0.1 + 0.9 * (1 + math.cos(math.pi * progress)) / 2

        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

        for epoch in range(n_epochs):
            model.train()
            total_loss = 0.0
            total_correct = 0
            total_masked = 0

            # Shuffle data
            indices = list(range(len(tokenized)))
            random.shuffle(indices)

            for batch_idx in range(n_batches):
                # Get batch
                batch_indices = indices[batch_idx * batch_size:(batch_idx + 1) * batch_size]
                batch = [tokenized[i] for i in batch_indices]
                input_ids = torch.tensor(batch, device=device)

                # Apply masking
                masked_ids, labels = mask_tokens_span(
                    input_ids,
                    vocab,
                    mask_prob=0.15,
                    max_span_len=1,
                    exclude_special_from_labels=True,
                )

                optimizer.zero_grad()
                logits, loss = model(masked_ids, labels)
                loss.backward()

                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)

                optimizer.step()
                scheduler.step()

                total_loss += loss.item()

                # Accuracy
                with torch.no_grad():
                    preds = logits.argmax(dim=-1)
                    mask = labels != -100
                    total_correct += (preds[mask] == labels[mask]).sum().item()
                    total_masked += mask.sum().item()

            avg_loss = total_loss / n_batches
            acc = total_correct / max(1, total_masked)
            current_lr = scheduler.get_last_lr()[0]

            LOGGER.info("  Epoch %2d: loss=%.4f, acc=%.4f (%.2f%%), lr=%.6f",
                       epoch + 1, avg_loss, acc, acc * 100, current_lr)

        LOGGER.info("  Final accuracy with LR=%.0e: %.2f%%", lr, acc * 100)


if __name__ == "__main__":
    main()
