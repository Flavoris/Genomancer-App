"""Evaluate soft voting ensemble on validation set.

This module provides ensemble evaluation for multi-kmer Stage 1 models.
Supports both equal-weight soft voting and learned ensemble weights.

Usage:
    python evaluate_ensemble.py --config config.yaml
    python evaluate_ensemble.py --config config.yaml --checkpoints k3:path/k3.pt k4:path/k4.pt
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
import yaml
from tqdm.auto import tqdm

from dataset import (
    KmerVocabulary,
    compute_cksnap,
    compute_pseeiip,
    compute_pstnp,
    compute_tnc,
    _sanitize_sequence,
)
from ensemble_infer import build_model, load_checkpoint, infer_model_architecture
from ensemble_weights import optimize_ensemble_weights, save_ensemble_weights

LOGGER = logging.getLogger("gene_whisperer.evaluate_ensemble")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def select_device() -> torch.device:
    """Select the best available device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def compute_engineered_features(sequence: str, cfg: Optional[Dict] = None) -> torch.Tensor:
    """Compute engineered features: TNC(64) + PseEIIP(64) + CKSNAP(96) + PSTNP(64) = 288."""
    tnc = compute_tnc(sequence)
    pseeiip = compute_pseeiip(sequence)
    cksnap = compute_cksnap(sequence)
    pstnp = compute_pstnp(sequence)

    if cfg is not None:
        if not bool(cfg.get("stage1_feature_enable_tnc", True)):
            tnc = torch.zeros_like(tnc)
        if not bool(cfg.get("stage1_feature_enable_pseeiip", True)):
            pseeiip = torch.zeros_like(pseeiip)
        if not bool(cfg.get("stage1_feature_enable_cksnap", True)):
            cksnap = torch.zeros_like(cksnap)
        if not bool(cfg.get("stage1_feature_enable_pstnp", True)):
            pstnp = torch.zeros_like(pstnp)

    return torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)


def prepare_inputs(
    sequence: str, vocab: KmerVocabulary, max_bp_len: int, cfg: Optional[Dict] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Prepare tokenized inputs and engineered features."""
    tokens = vocab.tokenize(sequence, max_bp_len)
    engineered = compute_engineered_features(sequence, cfg)
    return tokens.unsqueeze(0), engineered.unsqueeze(0)


def _binary_roc_auc(probabilities: np.ndarray, targets: np.ndarray) -> float:
    """Compute ROC AUC for binary targets without external dependencies."""
    probs = np.asarray(probabilities, dtype=np.float64).reshape(-1)
    y = np.asarray(targets, dtype=np.float64).reshape(-1)
    y = (y >= 0.5).astype(np.int32)

    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return 0.5

    order = np.argsort(probs, kind="mergesort")
    sorted_probs = probs[order]
    sorted_y = y[order]

    n = len(sorted_probs)
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_probs[j + 1] == sorted_probs[i]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        ranks[i : j + 1] = avg_rank
        i = j + 1

    rank_sum_pos = float(ranks[sorted_y == 1].sum())
    auc = (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def calculate_metrics(probabilities: np.ndarray, labels: np.ndarray) -> Dict[str, float]:
    """Calculate classification metrics from probabilities and labels."""
    probs = probabilities.reshape(-1)
    targets = labels.reshape(-1)
    preds = (probs >= 0.5).astype(np.float32)

    correct = (preds == targets).sum()
    total = float(len(targets))
    accuracy = correct / total if total else 0.0

    tp = ((preds == 1) & (targets == 1)).sum()
    tn = ((preds == 0) & (targets == 0)).sum()
    fp = ((preds == 1) & (targets == 0)).sum()
    fn = ((preds == 0) & (targets == 1)).sum()

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0

    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / float(np.sqrt(denom)) if denom > 0 else 0.0

    roc_auc = _binary_roc_auc(probs, targets)

    return {
        "accuracy": float(accuracy),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "mcc": float(mcc),
        "roc_auc": float(roc_auc),
        "auc": float(roc_auc),  # Alias
        "specificity": float(specificity),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }


def load_validation_data(cfg: Dict, script_dir: Path) -> pd.DataFrame:
    """Load validation data from config-specified TSV file."""
    val_path = cfg.get("stage1_val", "../data/stage1_val.tsv")
    if not Path(val_path).is_absolute():
        val_path = (script_dir / val_path).resolve()
    else:
        val_path = Path(val_path)

    delimiter = cfg.get("delimiter", "\t")
    df = pd.read_csv(val_path, delimiter=delimiter)
    LOGGER.info("Loaded %d validation examples from %s", len(df), val_path)
    return df


def load_model_for_kmer(
    kmer: int,
    checkpoint_path: Path,
    cfg: Dict,
    device: torch.device,
) -> Tuple[torch.nn.Module, KmerVocabulary]:
    """Load a model and its vocabulary for a specific k-mer size."""
    script_dir = Path(__file__).resolve().parent

    # Find vocabulary
    vocab_cache_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
    if not vocab_cache_dir.is_absolute():
        vocab_cache_dir = (script_dir / vocab_cache_dir).resolve()

    # Try different vocab naming conventions
    vocab_paths = [
        vocab_cache_dir / f"k{kmer}_vocab.json",
        vocab_cache_dir / f"mlm_k{kmer}_vocab.json",
    ]

    vocab_path = None
    for vp in vocab_paths:
        if vp.exists():
            vocab_path = vp
            break

    if vocab_path is None:
        raise FileNotFoundError(f"No vocabulary found for k={kmer} in {vocab_cache_dir}")

    vocab = KmerVocabulary.load(vocab_path)
    LOGGER.info("Loaded k=%d vocabulary with %d tokens from %s", kmer, len(vocab.itos), vocab_path)

    # Infer architecture from checkpoint
    arch_overrides = infer_model_architecture(checkpoint_path)
    model_cfg = {**cfg, **arch_overrides}

    # Build and load model
    model = build_model(model_cfg, vocab, device)
    load_checkpoint(model, checkpoint_path, device)

    return model, vocab


def evaluate_soft_voting_ensemble(
    config_path: str,
    checkpoint_paths: Dict[int, str],
    weights: Optional[Dict[int, float]] = None,
) -> Dict:
    """
    Evaluate soft voting ensemble.

    Args:
        config_path: Path to config.yaml
        checkpoint_paths: Dict mapping k-mer to checkpoint path
        weights: Optional dict mapping k-mer to ensemble weight. Defaults to equal weights.

    Returns:
        Dict with ensemble metrics (accuracy, mcc, auc, etc.)
    """
    script_dir = Path(__file__).resolve().parent
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    device = select_device()
    LOGGER.info("Using device: %s", device)

    # Load models
    models: Dict[int, Tuple[torch.nn.Module, KmerVocabulary]] = {}
    kmers = sorted(checkpoint_paths.keys())

    for kmer in kmers:
        ckpt_path = Path(checkpoint_paths[kmer])
        if not ckpt_path.is_absolute():
            ckpt_path = (script_dir / ckpt_path).resolve()

        if not ckpt_path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

        LOGGER.info("Loading model for k=%d from %s", kmer, ckpt_path)
        model, vocab = load_model_for_kmer(kmer, ckpt_path, cfg, device)
        models[kmer] = (model, vocab)

    # Load validation data
    df = load_validation_data(cfg, script_dir)

    # Compute baseline weights
    if weights is None:
        ensemble_cfg = cfg.get("ensemble") if isinstance(cfg.get("ensemble"), dict) else {}
        config_weights = ensemble_cfg.get("weights")
        if not isinstance(config_weights, dict):
            config_weights = cfg.get("ensemble_weights", {})
        if not isinstance(config_weights, dict):
            config_weights = {}
        weights = {}
        for kmer in kmers:
            weights[kmer] = float(config_weights.get(str(kmer), config_weights.get(kmer, 1.0 / len(kmers))))
        total_weight = sum(weights.values())
        weights = {k: v / total_weight for k, v in weights.items()} if total_weight else {
            k: 1.0 / len(kmers) for k in kmers
        }

    LOGGER.info("Baseline ensemble weights: %s", weights)

    # Evaluate each model
    max_bp_len = int(cfg.get("max_bp_len", 81))
    all_probs: Dict[int, List[float]] = {k: [] for k in kmers}
    all_labels: List[float] = []

    LOGGER.info("Evaluating ensemble on %d samples...", len(df))

    with torch.no_grad():
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Evaluating ensemble"):
            sequence = str(row.get("sequence", row.get("Sequence", "")))
            label = float(row.get("label", row.get("Label", row.get("is_promoter", 0))))
            all_labels.append(label)

            for kmer in kmers:
                model, vocab = models[kmer]
                tokens, engineered = prepare_inputs(sequence, vocab, max_bp_len, cfg)
                tokens = tokens.to(device)
                engineered = engineered.to(device)

                logits = model(tokens, engineered)
                prob = torch.sigmoid(logits).squeeze().cpu().item()
                all_probs[kmer].append(prob)

    # Convert to numpy
    labels_np = np.array(all_labels)

    # Compute individual model metrics
    individual_metrics: Dict[int, Dict] = {}
    for kmer in kmers:
        probs_np = np.array(all_probs[kmer])
        metrics = calculate_metrics(probs_np, labels_np)
        individual_metrics[kmer] = metrics
        LOGGER.info(
            "k=%d: Acc=%.4f, MCC=%.4f, AUC=%.4f",
            kmer, metrics["accuracy"], metrics["mcc"], metrics["auc"]
        )

    def compute_ensemble_probs(weight_map: Dict[int, float]) -> np.ndarray:
        """Compute weighted ensemble probabilities."""
        probs = np.zeros(len(labels_np), dtype=np.float64)
        for kmer in kmers:
            probs += weight_map[kmer] * np.asarray(all_probs[kmer], dtype=np.float64)
        return probs

    baseline_probs = compute_ensemble_probs(weights)
    baseline_metrics = calculate_metrics(baseline_probs, labels_np)
    LOGGER.info(
        "Baseline ensemble: Acc=%.4f (%.2f%%), MCC=%.4f, AUC=%.4f",
        baseline_metrics["accuracy"],
        baseline_metrics["accuracy"] * 100,
        baseline_metrics["mcc"],
        baseline_metrics["auc"],
    )

    ensemble_cfg = cfg.get("ensemble") if isinstance(cfg.get("ensemble"), dict) else {}
    optimize_weights = bool(ensemble_cfg.get("optimize_weights", False))
    optimization_metric = str(ensemble_cfg.get("optimization_metric", "mcc"))
    min_weight = float(ensemble_cfg.get("min_weight", 0.01))
    grid_step = float(ensemble_cfg.get("grid_step", 0.02))
    max_candidates = int(ensemble_cfg.get("max_candidates", 50000))
    random_seed = int(ensemble_cfg.get("random_seed", 13))

    final_weights = dict(weights)
    final_metrics = baseline_metrics
    optimization_metadata = None

    if optimize_weights:
        try:
            optimized_weights, optimization_metadata = optimize_ensemble_weights(
                val_probs={k: np.asarray(all_probs[k], dtype=np.float64) for k in kmers},
                val_labels=labels_np,
                metric=optimization_metric,
                min_weight=min_weight,
                grid_step=grid_step,
                max_candidates=max_candidates,
                random_seed=random_seed,
                return_metadata=True,
            )
        except ValueError as exc:
            LOGGER.warning("Weight optimization failed: %s", exc)
        else:
            optimized_probs = compute_ensemble_probs(optimized_weights)
            optimized_metrics = calculate_metrics(optimized_probs, labels_np)
            final_weights = optimized_weights
            final_metrics = optimized_metrics

            LOGGER.info(
                "Optimized ensemble: Acc=%.4f (%.2f%%), MCC=%.4f, AUC=%.4f",
                optimized_metrics["accuracy"],
                optimized_metrics["accuracy"] * 100,
                optimized_metrics["mcc"],
                optimized_metrics["auc"],
            )

            weights_path_value = ensemble_cfg.get("weights_path")
            if weights_path_value:
                weights_path = Path(weights_path_value)
                if not weights_path.is_absolute():
                    weights_path = (script_dir / weights_path).resolve()
            else:
                weights_path = (script_dir / "../artifacts/ensemble_weights.json").resolve()

            metadata = dict(optimization_metadata or {})
            metadata.update(
                {
                    "kmers": kmers,
                    "num_samples": int(len(labels_np)),
                }
            )
            save_ensemble_weights(weights_path, optimized_weights, metadata)
            LOGGER.info("Saved optimized weights to %s", weights_path)

    results = {
        **final_metrics,
        "individual_metrics": individual_metrics,
        "weights": final_weights,
        "kmers": kmers,
        "num_samples": len(labels_np),
    }

    if optimize_weights and optimization_metadata:
        results.update(
            {
                "baseline_weights": weights,
                "baseline_metrics": baseline_metrics,
                "optimized_weights": final_weights,
                "optimized_metrics": final_metrics,
                "optimization": optimization_metadata,
            }
        )

    return results


def main() -> None:
    """CLI entry point for ensemble evaluation."""
    parser = argparse.ArgumentParser(
        description="Evaluate soft voting ensemble on validation set"
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file",
    )
    parser.add_argument(
        "--checkpoints",
        nargs="+",
        default=None,
        help="Checkpoint paths as k:path pairs (e.g., 3:path/k3.pt 4:path/k4.pt)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save results JSON",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    # Parse checkpoint paths
    if args.checkpoints:
        checkpoint_paths = {}
        for item in args.checkpoints:
            kmer_str, path = item.split(":", 1)
            checkpoint_paths[int(kmer_str)] = path
    else:
        # Use default paths from config
        stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k", {})
        if not stage1_ckpt_by_k:
            # Fall back to artifacts directory
            artifacts_dir = (script_dir / "../artifacts/checkpoints").resolve()
            kmers = cfg.get("multi_scale_kmers", [3, 4, 5, 6])
            checkpoint_paths = {}
            for k in kmers:
                ckpt_path = artifacts_dir / f"stage1_k{k}.pt"
                if ckpt_path.exists():
                    checkpoint_paths[k] = str(ckpt_path)
        else:
            checkpoint_paths = {int(k): v for k, v in stage1_ckpt_by_k.items()}

    if not checkpoint_paths:
        LOGGER.error("No checkpoint paths found. Use --checkpoints or configure stage1_ckpt_by_k in config.")
        return

    LOGGER.info("Evaluating ensemble with k-mers: %s", sorted(checkpoint_paths.keys()))

    results = evaluate_soft_voting_ensemble(
        config_path=str(config_path),
        checkpoint_paths=checkpoint_paths,
    )

    # Print results
    print("\n" + "=" * 70)
    print("ENSEMBLE EVALUATION RESULTS")
    print("=" * 70)

    print("\nIndividual Model Results:")
    for kmer in sorted(results["individual_metrics"].keys()):
        m = results["individual_metrics"][kmer]
        print(f"  k={kmer}: Acc={m['accuracy']*100:.2f}%, MCC={m['mcc']:.4f}, AUC={m['auc']:.4f}")

    if "baseline_metrics" in results and "optimized_metrics" in results:
        baseline = results["baseline_metrics"]
        optimized = results["optimized_metrics"]
        print(f"\nBaseline Ensemble (weights={results['baseline_weights']}):")
        print(f"  Accuracy:    {baseline['accuracy']:.4f} ({baseline['accuracy']*100:.2f}%)")
        print(f"  Precision:   {baseline['precision']:.4f}")
        print(f"  Recall:      {baseline['recall']:.4f}")
        print(f"  F1:          {baseline['f1']:.4f}")
        print(f"  MCC:         {baseline['mcc']:.4f}")
        print(f"  AUC:         {baseline['auc']:.4f}")
        print(f"  Specificity: {baseline['specificity']:.4f}")

        print(f"\nOptimized Ensemble (weights={results['optimized_weights']}):")
        print(f"  Accuracy:    {optimized['accuracy']:.4f} ({optimized['accuracy']*100:.2f}%)")
        print(f"  Precision:   {optimized['precision']:.4f}")
        print(f"  Recall:      {optimized['recall']:.4f}")
        print(f"  F1:          {optimized['f1']:.4f}")
        print(f"  MCC:         {optimized['mcc']:.4f}")
        print(f"  AUC:         {optimized['auc']:.4f}")
        print(f"  Specificity: {optimized['specificity']:.4f}")
    else:
        print(f"\nEnsemble Results (weights={results['weights']}):")
        print(f"  Accuracy:    {results['accuracy']:.4f} ({results['accuracy']*100:.2f}%)")
        print(f"  Precision:   {results['precision']:.4f}")
        print(f"  Recall:      {results['recall']:.4f}")
        print(f"  F1:          {results['f1']:.4f}")
        print(f"  MCC:         {results['mcc']:.4f}")
        print(f"  AUC:         {results['auc']:.4f}")
        print(f"  Specificity: {results['specificity']:.4f}")
    print("=" * 70)

    # Save results if output path specified
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)

        print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
