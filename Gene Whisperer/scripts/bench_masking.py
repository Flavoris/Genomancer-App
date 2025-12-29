#!/usr/bin/env python3
"""Micro-benchmark for mask_tokens_span to identify MLM masking performance bottlenecks."""

import argparse
import random
import sys
import time
from pathlib import Path

import torch

# Add training directory to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = SCRIPT_DIR.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

from dataset import KmerVocabulary
from pretrain_mlm import mask_tokens_span


def build_fake_vocab(k: int = 6) -> KmerVocabulary:
    """Build a small vocabulary for benchmarking."""
    bases = ["A", "C", "G", "T"]
    # Generate all k-mers
    from itertools import product
    all_kmers = ["".join(kmer) for kmer in product(bases, repeat=k)]
    return KmerVocabulary(
        k=k,
        tokens=all_kmers,
        unk_token="<UNK>",
        pad_token="<PAD>",
        mask_token="[MASK]",
    )


def create_fake_batch(
    batch_size: int,
    seq_len: int,
    vocab: KmerVocabulary,
    device: str = "cpu",
) -> torch.LongTensor:
    """Create a fake token batch with non-PAD token ids."""
    # Get valid token IDs (excluding special tokens)
    valid_ids = vocab._base_token_ids
    if not valid_ids:
        valid_ids = list(range(len(vocab.itos) - 3))  # Exclude last 3 (special tokens)

    # Generate random token IDs
    tokens = torch.tensor(
        [[random.choice(valid_ids) for _ in range(seq_len)] for _ in range(batch_size)],
        dtype=torch.long,
        device=device,
    )
    return tokens


def benchmark_masking(
    seq_len: int,
    batch_size: int,
    iters: int,
    vocab: KmerVocabulary,
    device: str = "cpu",
    warmup_iters: int = 10,
) -> dict:
    """Run the masking benchmark and return timing results."""
    # Create fake batch
    inputs = create_fake_batch(batch_size, seq_len, vocab, device)

    # Warmup
    print(f"Warming up ({warmup_iters} iterations)...")
    for _ in range(warmup_iters):
        mask_tokens_span(inputs, vocab)

    # Benchmark
    print(f"Running benchmark ({iters} iterations)...")
    times = []
    for i in range(iters):
        start = time.perf_counter()
        mask_tokens_span(inputs, vocab)
        end = time.perf_counter()
        times.append(end - start)

        if (i + 1) % 50 == 0:
            print(f"  Progress: {i + 1}/{iters}")

    # Calculate stats
    total_time = sum(times)
    avg_ms = (total_time / iters) * 1000
    total_tokens = batch_size * seq_len * iters
    tokens_per_sec = total_tokens / total_time

    return {
        "seq_len": seq_len,
        "batch_size": batch_size,
        "iters": iters,
        "total_time_s": total_time,
        "avg_ms_per_iter": avg_ms,
        "tokens_per_sec": tokens_per_sec,
        "min_ms": min(times) * 1000,
        "max_ms": max(times) * 1000,
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark mask_tokens_span performance")
    parser.add_argument("--seq_len", type=int, default=81, help="Sequence length (default: 81)")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size (default: 64)")
    parser.add_argument("--iters", type=int, default=200, help="Number of iterations (default: 200)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--k", type=int, default=6, help="K-mer size for vocab (default: 6)")
    parser.add_argument("--device", type=str, default="cpu", help="Device: cpu or cuda (default: cpu)")
    parser.add_argument("--warmup", type=int, default=10, help="Warmup iterations (default: 10)")
    args = parser.parse_args()

    # Set seeds
    random.seed(args.seed)
    torch.manual_seed(args.seed)

    print("=" * 60)
    print("mask_tokens_span Benchmark")
    print("=" * 60)
    print(f"Configuration:")
    print(f"  seq_len:    {args.seq_len}")
    print(f"  batch_size: {args.batch_size}")
    print(f"  iters:      {args.iters}")
    print(f"  seed:       {args.seed}")
    print(f"  k:          {args.k}")
    print(f"  device:     {args.device}")
    print("=" * 60)

    # Build vocab
    print("\nBuilding vocabulary...")
    vocab = build_fake_vocab(k=args.k)
    print(f"Vocabulary size: {len(vocab.itos)} tokens")

    # Run benchmark
    print("\n" + "-" * 60)
    results = benchmark_masking(
        seq_len=args.seq_len,
        batch_size=args.batch_size,
        iters=args.iters,
        vocab=vocab,
        device=args.device,
        warmup_iters=args.warmup,
    )

    # Print results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"  Sequence length:     {results['seq_len']}")
    print(f"  Batch size:          {results['batch_size']}")
    print(f"  Iterations:          {results['iters']}")
    print(f"  Total time:          {results['total_time_s']:.3f} s")
    print(f"  Avg time per iter:   {results['avg_ms_per_iter']:.3f} ms")
    print(f"  Min time per iter:   {results['min_ms']:.3f} ms")
    print(f"  Max time per iter:   {results['max_ms']:.3f} ms")
    print(f"  Tokens/sec:          {results['tokens_per_sec']:,.0f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
