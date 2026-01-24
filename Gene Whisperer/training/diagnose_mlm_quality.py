#!/usr/bin/env python3
"""Diagnose MLM pre-training quality across k-mer values.

Analyzes whether MLM pre-training is producing useful representations by:
1. Measuring masked token prediction accuracy (top-1, top-5)
2. Computing prediction entropy (confidence)
3. Evaluating embedding clustering quality (silhouette score)

Usage:
    cd "Gene Whisperer/training"
    python diagnose_mlm_quality.py
    python diagnose_mlm_quality.py --kmers 3 4
    python diagnose_mlm_quality.py --max-samples 200 --mask-prob 0.20
    python diagnose_mlm_quality.py --save-embeddings
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import torch

from dataset import KmerVocabulary
from mlm_quality_metrics import (
    EncoderConfig,
    MLMMetrics,
    build_mlm_model,
    compute_mlm_predictions,
    compute_silhouette,
    extract_embeddings,
    interpret_silhouette,
    load_encoder_config,
)

LOGGER = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Default paths relative to training/ directory
ARTIFACTS_DIR = Path("../artifacts")
VOCABS_DIR = ARTIFACTS_DIR / "vocabs"
DATA_DIR = Path("../data")
VAL_DATA_PATH = DATA_DIR / "stage1_val.tsv"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Diagnose MLM pre-training quality across k-mer values",
    )
    parser.add_argument(
        "--kmers", type=int, nargs="+", default=[3, 4, 5, 6],
        help="K-mer values to analyze (default: 3 4 5 6)",
    )
    parser.add_argument(
        "--max-samples", type=int, default=None,
        help="Maximum validation samples to use (default: all)",
    )
    parser.add_argument(
        "--mask-prob", type=float, default=0.15,
        help="Fraction of tokens to mask (default: 0.15)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=32,
        help="Inference batch size (default: 32)",
    )
    parser.add_argument(
        "--save-embeddings", action="store_true",
        help="Save t-SNE embedding plots to artifacts/",
    )
    parser.add_argument(
        "--val-data", type=str, default=None,
        help="Path to validation TSV (default: ../data/stage1_val.tsv)",
    )
    parser.add_argument(
        "--artifacts-dir", type=str, default=None,
        help="Path to artifacts directory (default: ../artifacts)",
    )
    parser.add_argument(
        "--device", type=str, default=None,
        help="Force device: 'cpu', 'cuda', or 'mps' (default: auto-detect)",
    )
    parser.add_argument(
        "--seed", type=int, default=1337,
        help="Random seed for reproducibility (default: 1337)",
    )
    return parser.parse_args()


def load_validation_data(
    val_path: Path, max_samples: Optional[int] = None,
) -> Tuple[List[str], np.ndarray]:
    """Load validation sequences and labels from TSV.

    Returns:
        Tuple of (sequences, labels) where labels are 0/1 numpy array
    """
    if not val_path.exists():
        raise FileNotFoundError(f"Validation data not found: {val_path}")

    df = pd.read_csv(val_path, sep="\t")
    if "sequence" not in df.columns or "is_promoter" not in df.columns:
        raise ValueError(
            f"Expected columns 'sequence' and 'is_promoter' in {val_path}, "
            f"got: {list(df.columns)}"
        )

    sequences = df["sequence"].tolist()
    labels = df["is_promoter"].values.astype(int)

    if max_samples and len(sequences) > max_samples:
        indices = np.random.choice(len(sequences), max_samples, replace=False)
        sequences = [sequences[i] for i in indices]
        labels = labels[indices]

    return sequences, labels


def tokenize_sequences(
    sequences: List[str], vocab: KmerVocabulary, max_bp_len: int = 81
) -> torch.LongTensor:
    """Tokenize all sequences using the given vocabulary."""
    tokens = [vocab.tokenize(seq, max_bp_len) for seq in sequences]
    return torch.stack(tokens)


def get_device() -> torch.device:
    """Select best available device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def save_tsne_plot(
    embeddings: np.ndarray,
    labels: np.ndarray,
    k: int,
    output_dir: Path,
) -> None:
    """Generate and save t-SNE visualization of embeddings."""
    try:
        from sklearn.manifold import TSNE
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("sklearn/matplotlib not available, skipping t-SNE plot")
        return

    # Subsample for t-SNE speed
    max_points = 1000
    if len(embeddings) > max_points:
        idx = np.random.choice(len(embeddings), max_points, replace=False)
        embeddings = embeddings[idx]
        labels = labels[idx]

    perplexity = min(30, len(embeddings) - 1)
    # Use method='exact' to avoid Barnes-Hut segfault on macOS ARM;
    # acceptable since we subsample to <=1000 points above.
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity, method="exact")
    coords = tsne.fit_transform(embeddings.astype(np.float64))

    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(
        coords[:, 0], coords[:, 1],
        c=labels, cmap="RdYlBu", alpha=0.6, s=15,
    )
    ax.set_title(f"MLM Embeddings k={k} (colored by promoter label)")
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    plt.colorbar(scatter, ax=ax, label="is_promoter")

    output_path = output_dir / f"mlm_tsne_k{k}.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info("Saved t-SNE plot to %s", output_path)


def analyze_kmer(
    k: int,
    sequences: List[str],
    labels: np.ndarray,
    artifacts_dir: Path,
    device: torch.device,
    batch_size: int,
    mask_prob: float,
    save_embeddings: bool,
) -> Optional[MLMMetrics]:
    """Run full MLM quality analysis for a single k-mer value."""
    vocab_path = artifacts_dir / "vocabs" / f"mlm_k{k}_vocab.json"
    checkpoint_path = artifacts_dir / f"mlm_encoder_k{k}.pt"
    metadata_path = artifacts_dir / f"mlm_encoder_k{k}.metadata.json"

    # Load vocabulary
    if not vocab_path.exists():
        LOGGER.error("Vocab not found for k=%d: %s", k, vocab_path)
        return None
    vocab = KmerVocabulary.load(vocab_path)
    LOGGER.info("k=%d: vocab_size=%d, checkpoint=%s", k, len(vocab.itos), checkpoint_path.name)

    # Load architecture config from metadata
    encoder_cfg = load_encoder_config(metadata_path)
    LOGGER.info(
        "  Architecture: dim=%d, layers=%d, heads=%d",
        encoder_cfg.embedding_dim, encoder_cfg.num_layers, encoder_cfg.num_heads,
    )

    # Build model and load weights
    encoder, weights_loaded = build_mlm_model(vocab, encoder_cfg, checkpoint_path)
    if not weights_loaded:
        LOGGER.warning("  Skipping k=%d: checkpoint not available", k)
        return MLMMetrics(
            k=k, top1_accuracy=0.0, top5_accuracy=0.0,
            avg_entropy=0.0, checkpoint_exists=False,
        )

    # Tokenize sequences
    max_bp_len = 81
    token_ids = tokenize_sequences(sequences, vocab, max_bp_len)
    LOGGER.info("  Tokenized %d sequences → shape %s", len(sequences), tuple(token_ids.shape))

    # Compute MLM prediction metrics
    metrics = compute_mlm_predictions(
        encoder, vocab, token_ids, device, batch_size, mask_prob
    )

    # Compute embedding quality
    embeddings = extract_embeddings(encoder, token_ids, device, batch_size)
    sil_score = compute_silhouette(embeddings, labels)
    metrics.silhouette_score = sil_score

    # Save t-SNE plot if requested
    if save_embeddings:
        save_tsne_plot(embeddings, labels, k, artifacts_dir)

    return metrics


def print_results(all_metrics: List[MLMMetrics]) -> None:
    """Print formatted analysis results."""
    print("\n" + "=" * 70)
    print("MLM QUALITY ANALYSIS")
    print("=" * 70)

    # Find best performing k-mer
    valid = [m for m in all_metrics if m.checkpoint_exists]
    best_k = max(valid, key=lambda m: m.top1_accuracy).k if valid else None

    for m in all_metrics:
        if not m.checkpoint_exists:
            print(f"k={m.k}: CHECKPOINT NOT FOUND - run pretrain_mlm.py first")
            continue

        annotation = ""
        if m.k == best_k:
            annotation = " ← Best"
        elif m.top1_accuracy < 0.45:
            annotation = " ← Undertrained"

        print(
            f"k={m.k}: Top-1 Acc: {m.top1_accuracy*100:.1f}%, "
            f"Top-5 Acc: {m.top5_accuracy*100:.1f}%, "
            f"Avg Entropy: {m.avg_entropy:.2f}"
            f"{annotation}"
        )

    print()
    print("Embedding Clustering (Silhouette Score):")
    for m in all_metrics:
        if not m.checkpoint_exists:
            print(f"k={m.k}: N/A (no checkpoint)")
            continue
        score = m.silhouette_score if m.silhouette_score is not None else float("nan")
        interpretation = interpret_silhouette(score)
        if np.isnan(score):
            print(f"k={m.k}: N/A")
        else:
            print(f"k={m.k}: {score:.3f} ({interpretation})")

    print("=" * 70)

    # Print summary interpretation
    if valid:
        print("\nINTERPRETATION:")
        for m in valid:
            if m.top1_accuracy > 0.60:
                print(f"  k={m.k}: MLM learned strong k-mer patterns ✓")
            elif m.top1_accuracy > 0.40:
                print(f"  k={m.k}: MLM learned some patterns, could benefit from more training")
            else:
                print(f"  k={m.k}: MLM may be undertrained or vocab too large")


def main() -> None:
    """Main entry point for MLM quality diagnosis."""
    args = parse_args()

    # Set seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # Resolve paths
    artifacts_dir = Path(args.artifacts_dir) if args.artifacts_dir else ARTIFACTS_DIR
    val_path = Path(args.val_data) if args.val_data else VAL_DATA_PATH

    device = torch.device(args.device) if args.device else get_device()
    LOGGER.info("Device: %s", device)

    # Load validation data
    sequences, labels = load_validation_data(val_path, args.max_samples)
    LOGGER.info(
        "Loaded %d validation sequences (%d promoters, %d non-promoters)",
        len(sequences), labels.sum(), len(labels) - labels.sum(),
    )

    # Analyze each k-mer value
    all_metrics: List[MLMMetrics] = []
    for k in args.kmers:
        LOGGER.info("\n--- Analyzing k=%d ---", k)
        metrics = analyze_kmer(
            k=k,
            sequences=sequences,
            labels=labels,
            artifacts_dir=artifacts_dir,
            device=device,
            batch_size=args.batch_size,
            mask_prob=args.mask_prob,
            save_embeddings=args.save_embeddings,
        )
        if metrics:
            all_metrics.append(metrics)

    # Print results
    print_results(all_metrics)


if __name__ == "__main__":
    main()
