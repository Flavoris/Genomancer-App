#!/usr/bin/env python3
"""
Reality check: Sweep decision thresholds on Stage 1 validation set.

Determines whether the 0.5 threshold is optimal or if a different threshold
would yield higher accuracy.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import torch

# Add training directory to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = SCRIPT_DIR.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

import yaml

from dataset import build_dataloaders
from ensemble_infer import build_model, load_checkpoint


def select_device() -> torch.device:
    """Select best available device: CUDA > MPS > CPU."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def compute_mcc(tp: int, tn: int, fp: int, fn: int) -> float:
    """Compute Matthews Correlation Coefficient."""
    numerator = tp * tn - fp * fn
    denominator = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if denominator == 0:
        return 0.0
    return numerator / np.sqrt(denominator)


def compute_metrics_at_threshold(
    probs: np.ndarray, labels: np.ndarray, threshold: float
) -> dict:
    """Compute accuracy, F1, and MCC at a given threshold."""
    preds = (probs >= threshold).astype(int)
    correct = (preds == labels).sum()
    accuracy = correct / len(labels)

    tp = ((preds == 1) & (labels == 1)).sum()
    tn = ((preds == 0) & (labels == 0)).sum()
    fp = ((preds == 1) & (labels == 0)).sum()
    fn = ((preds == 0) & (labels == 1)).sum()

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    mcc = compute_mcc(int(tp), int(tn), int(fp), int(fn))

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "mcc": mcc,
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Sweep thresholds on Stage 1 val set to find optimal decision boundary."
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (default: config.yaml)",
    )
    parser.add_argument(
        "--ckpt",
        default=None,
        help="Path to Stage 1 checkpoint. If omitted, auto-picks from stage1_ckpt_by_k[kmer].",
    )
    parser.add_argument(
        "--min_thr",
        type=float,
        default=0.01,
        help="Minimum threshold (default: 0.01)",
    )
    parser.add_argument(
        "--max_thr",
        type=float,
        default=0.99,
        help="Maximum threshold (default: 0.99)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=99,
        help="Number of threshold steps (default: 99)",
    )
    args = parser.parse_args()

    # Resolve config path
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (TRAINING_DIR / config_path).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    # Load config
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    # Change to training directory for relative paths in config
    os.chdir(TRAINING_DIR)

    # Determine checkpoint path
    ckpt_path = args.ckpt
    if ckpt_path is None:
        kmer = int(cfg.get("kmer", 6))
        ckpt_by_k = cfg.get("stage1_ckpt_by_k", {})
        if kmer in ckpt_by_k:
            ckpt_path = ckpt_by_k[kmer]
        elif str(kmer) in ckpt_by_k:
            ckpt_path = ckpt_by_k[str(kmer)]
        else:
            print(f"ERROR: No checkpoint found for kmer={kmer} in stage1_ckpt_by_k.")
            print(f"Available keys: {list(ckpt_by_k.keys())}")
            print("Use --ckpt to specify a checkpoint path manually.")
            sys.exit(1)

    ckpt_path = Path(ckpt_path)
    if not ckpt_path.is_absolute():
        ckpt_path = (TRAINING_DIR / ckpt_path).resolve()
    if not ckpt_path.exists():
        print(f"ERROR: Checkpoint not found: {ckpt_path}")
        sys.exit(1)

    print(f"Config: {config_path}")
    print(f"Checkpoint: {ckpt_path}")
    print(f"Threshold range: [{args.min_thr:.2f}, {args.max_thr:.2f}], steps={args.steps}")

    # Select device
    device = select_device()
    print(f"Device: {device}")

    # Build dataloaders and get val loader
    print("\nLoading validation data...")
    loaders = build_dataloaders(cfg)
    val_loader = loaders["stage1"]["val"]
    vocab = val_loader.dataset.vocab

    print(f"Validation samples: {len(val_loader.dataset)}")
    print(f"Vocab size: {len(vocab.itos)}, k={vocab.k}")

    # Build and load model
    print("\nBuilding model...")
    model = build_model(cfg, vocab, device)
    load_checkpoint(model, ckpt_path, device)
    model.eval()

    # Collect predictions with logits vs probs sanity check
    print("\nRunning inference with logits/probs sanity check...")
    all_probs_model = []
    all_probs2 = []
    all_logits = []
    all_labels = []

    with torch.no_grad():
        for batch_idx, batch in enumerate(val_loader):
            tokens, engineered, labels = batch
            tokens = tokens.to(device)
            engineered = engineered.to(device)

            # Get model default output (should be probs after sigmoid)
            probs_model = model(tokens, engineered_features=engineered)

            # Get both probs and logits via return_logits=True
            probs2, logits = model(tokens, engineered_features=engineered, return_logits=True)

            all_probs_model.append(probs_model.squeeze(-1).cpu())
            all_probs2.append(probs2.squeeze(-1).cpu())
            all_logits.append(logits.squeeze(-1).cpu())
            all_labels.append(labels)

    # Concatenate all batches
    probs_model = torch.cat(all_probs_model).numpy()
    probs2 = torch.cat(all_probs2).numpy()
    logits = torch.cat(all_logits).numpy()
    labels = np.concatenate(all_labels).astype(int)

    # =========================================================================
    # LOGITS VS PROBS SANITY CHECK
    # =========================================================================
    print("\n" + "=" * 70)
    print("LOGITS VS PROBS SANITY CHECK")
    print("=" * 70)

    # Check 1: probs_model should equal sigmoid(logits)
    sigmoid_of_logits = 1.0 / (1.0 + np.exp(-logits))  # numpy sigmoid
    diff1 = np.abs(probs_model - sigmoid_of_logits)
    max_diff1 = diff1.max()

    # Check 2: probs_model should equal probs2 (both from same forward pass logic)
    diff2 = np.abs(probs_model - probs2)
    max_diff2 = diff2.max()

    print(f"\nAssertion 1: probs_model ≈ sigmoid(logits)")
    print(f"  Max absolute difference: {max_diff1:.2e}")
    if max_diff1 < 1e-6:
        print("  ✓ PASSED: probs_model equals sigmoid(logits)")
    else:
        print("  ✗ FAILED: Mismatch! probs not equal to sigmoid(logits).")
        print("           Check for double sigmoid or inconsistent forward path.")

    print(f"\nAssertion 2: probs_model ≈ probs2 (from return_logits path)")
    print(f"  Max absolute difference: {max_diff2:.2e}")
    if max_diff2 < 1e-6:
        print("  ✓ PASSED: probs_model equals probs2")
    else:
        print("  ✗ FAILED: Mismatch! Default output differs from return_logits output.")
        print("           Check for inconsistent forward path.")

    # Print logits and probs statistics
    print(f"\nLogits statistics:")
    print(f"  min:  {logits.min():.4f}")
    print(f"  max:  {logits.max():.4f}")
    print(f"  mean: {logits.mean():.4f}")

    print(f"\nprobs_model statistics:")
    print(f"  min:  {probs_model.min():.4f}")
    print(f"  max:  {probs_model.max():.4f}")
    print(f"  mean: {probs_model.mean():.4f}")

    # Proportion of probs < 0.5 and >= 0.5
    below_half = (probs_model < 0.5).sum()
    above_half = (probs_model >= 0.5).sum()
    print(f"\nProbs distribution:")
    print(f"  probs < 0.5:  {below_half} ({100*below_half/len(probs_model):.1f}%)")
    print(f"  probs >= 0.5: {above_half} ({100*above_half/len(probs_model):.1f}%)")

    # Interpretation
    print(f"\nInterpretation:")
    if max_diff1 < 1e-6:
        print("  → Model forward path is consistent (no double sigmoid).")
    else:
        print("  → WARNING: Forward path may have double sigmoid or other issue!")

    if logits.mean() > 0.5:
        print(f"  → Logits are shifted positive (mean={logits.mean():.4f} > 0).")
        print("    Model may need calibration or bias adjustment.")
    elif logits.mean() < -0.5:
        print(f"  → Logits are shifted negative (mean={logits.mean():.4f} < 0).")
        print("    Model may need calibration or bias adjustment.")
    else:
        print(f"  → Logits are reasonably centered (mean={logits.mean():.4f}).")

    print("=" * 70)

    # Use probs_model for threshold sweep (the correct probabilities)
    probs = probs_model

    print(f"\nTotal predictions: {len(probs)}")
    print(f"Label distribution: {(labels == 1).sum()} positive, {(labels == 0).sum()} negative")
    print(f"Prob range: [{probs.min():.4f}, {probs.max():.4f}], mean={probs.mean():.4f}")

    # Sweep thresholds
    thresholds = np.linspace(args.min_thr, args.max_thr, args.steps)
    results = []

    for thr in thresholds:
        metrics = compute_metrics_at_threshold(probs, labels, thr)
        metrics["threshold"] = thr
        results.append(metrics)

    # Sort by accuracy descending
    results_sorted = sorted(results, key=lambda x: x["accuracy"], reverse=True)

    # Find baseline (0.5) metrics
    baseline_metrics = compute_metrics_at_threshold(probs, labels, 0.5)

    # Print results
    print("\n" + "=" * 70)
    print("THRESHOLD SWEEP RESULTS")
    print("=" * 70)

    print(f"\nBaseline @ threshold=0.5:")
    print(f"  Accuracy: {baseline_metrics['accuracy']:.4f}")
    print(f"  F1:       {baseline_metrics['f1']:.4f}")
    print(f"  MCC:      {baseline_metrics['mcc']:.4f}")
    print(f"  TP={baseline_metrics['tp']}, TN={baseline_metrics['tn']}, FP={baseline_metrics['fp']}, FN={baseline_metrics['fn']}")

    best = results_sorted[0]
    print(f"\nBest accuracy @ threshold={best['threshold']:.4f}:")
    print(f"  Accuracy: {best['accuracy']:.4f} ({(best['accuracy'] - baseline_metrics['accuracy'])*100:+.2f}% vs 0.5)")
    print(f"  F1:       {best['f1']:.4f}")
    print(f"  MCC:      {best['mcc']:.4f}")
    print(f"  TP={best['tp']}, TN={best['tn']}, FP={best['fp']}, FN={best['fn']}")

    print(f"\nTop 10 thresholds by accuracy:")
    print("-" * 70)
    print(f"{'Rank':<6}{'Threshold':<12}{'Accuracy':<12}{'F1':<10}{'MCC':<10}")
    print("-" * 70)
    for rank, r in enumerate(results_sorted[:10], 1):
        print(f"{rank:<6}{r['threshold']:<12.4f}{r['accuracy']:<12.4f}{r['f1']:<10.4f}{r['mcc']:<10.4f}")

    # Summary
    print("\n" + "=" * 70)
    if best["accuracy"] > baseline_metrics["accuracy"] + 0.001:
        improvement = (best["accuracy"] - baseline_metrics["accuracy"]) * 100
        print(f"FINDING: Threshold {best['threshold']:.4f} improves accuracy by {improvement:.2f}% over 0.5")
    else:
        print("FINDING: Threshold 0.5 is near-optimal for this model/dataset.")
    print("=" * 70)


if __name__ == "__main__":
    main()
