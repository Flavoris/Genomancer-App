"""Evaluate Stage 1 ensemble over k=4 and k=6 models on validation data.

This script evaluates:
1. k=4 model alone
2. k=6 model alone
3. Ensemble average of probabilities (0.5*(p4+p6))

Reports accuracy, precision, recall, F1, MCC, and AUROC for each.
"""
from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
import yaml
from tqdm.auto import tqdm

from dataset import KmerVocabulary, compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc, _sanitize_sequence
from ensemble_infer import build_model, load_checkpoint, select_device

LOGGER = logging.getLogger("gene_whisperer.eval_ensemble")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def load_vocab(vocab_path: Path) -> KmerVocabulary:
    """Load vocabulary from JSON file."""
    return KmerVocabulary.load(vocab_path)


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
        "specificity": float(specificity),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }


# Small model architecture parameters (must match train_stage1_k*_small.sh)
SMALL_MODEL_CONFIG = {
    "embedding_dim": 192,
    "transformer_layers": 4,
    "transformer_heads": 4,
    "transformer_ff_dim": 768,
    "transformer_dropout": 0.12,
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


def evaluate_model(
    model: torch.nn.Module,
    vocab: KmerVocabulary,
    df: pd.DataFrame,
    cfg: Dict,
    device: torch.device,
) -> Tuple[np.ndarray, np.ndarray]:
    """Evaluate a single model on the validation data.

    Returns:
        Tuple of (probabilities, labels) as numpy arrays
    """
    max_bp_len = int(cfg.get("max_bp_len", 81))
    model.eval()

    all_probs = []
    all_labels = []

    with torch.no_grad():
        for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"Evaluating k={vocab.k}"):
            sequence = str(row.get("sequence", row.get("Sequence", "")))
            # Support multiple label column names: label, Label, is_promoter
            label = float(row.get("label", row.get("Label", row.get("is_promoter", 0))))

            tokens, engineered = prepare_inputs(sequence, vocab, max_bp_len, cfg)
            tokens = tokens.to(device)
            engineered = engineered.to(device)

            prob = model(tokens, engineered)
            prob = prob.squeeze().cpu().item()

            all_probs.append(prob)
            all_labels.append(label)

    return np.array(all_probs), np.array(all_labels)


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate Stage 1 ensemble (k=4, k=6)")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()

    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    device = select_device()
    LOGGER.info("Using device: %s", device)

    # Define paths for vocabs and checkpoints
    artifacts_dir = (script_dir / "../artifacts").resolve()
    vocab_k4_path = artifacts_dir / "vocabs" / "k4_vocab.json"
    vocab_k6_path = artifacts_dir / "vocabs" / "k6_vocab.json"
    ckpt_k4_path = artifacts_dir / "checkpoints" / "stage1_k4.pt"
    ckpt_k6_path = artifacts_dir / "checkpoints" / "stage1_k6.pt"

    # Check that files exist
    missing = []
    for p in [vocab_k4_path, vocab_k6_path, ckpt_k4_path, ckpt_k6_path]:
        if not p.exists():
            missing.append(str(p))
    if missing:
        LOGGER.error("Missing required files:\n  %s", "\n  ".join(missing))
        raise FileNotFoundError(f"Missing: {missing}")

    # Load vocabularies
    LOGGER.info("Loading vocabularies...")
    vocab_k4 = load_vocab(vocab_k4_path)
    vocab_k6 = load_vocab(vocab_k6_path)
    LOGGER.info("k=4 vocab size: %d", len(vocab_k4.itos))
    LOGGER.info("k=6 vocab size: %d", len(vocab_k6.itos))

    # Merge small model config with base config
    # The small model training scripts use reduced architecture parameters
    model_cfg = {**cfg, **SMALL_MODEL_CONFIG}
    LOGGER.info(
        "Using small model config: embedding_dim=%d, transformer_layers=%d, transformer_heads=%d",
        model_cfg["embedding_dim"],
        model_cfg["transformer_layers"],
        model_cfg["transformer_heads"],
    )

    # Build and load models
    LOGGER.info("Building and loading k=4 model...")
    model_k4 = build_model(model_cfg, vocab_k4, device)
    load_checkpoint(model_k4, ckpt_k4_path, device)

    LOGGER.info("Building and loading k=6 model...")
    model_k6 = build_model(model_cfg, vocab_k6, device)
    load_checkpoint(model_k6, ckpt_k6_path, device)

    # Load validation data
    df = load_validation_data(cfg, script_dir)

    # Evaluate k=4 model
    LOGGER.info("=" * 60)
    LOGGER.info("Evaluating k=4 model...")
    probs_k4, labels = evaluate_model(model_k4, vocab_k4, df, cfg, device)
    metrics_k4 = calculate_metrics(probs_k4, labels)

    # Evaluate k=6 model
    LOGGER.info("=" * 60)
    LOGGER.info("Evaluating k=6 model...")
    probs_k6, _ = evaluate_model(model_k6, vocab_k6, df, cfg, device)
    metrics_k6 = calculate_metrics(probs_k6, labels)

    # Compute ensemble (simple average)
    LOGGER.info("=" * 60)
    LOGGER.info("Computing ensemble (0.5 * p4 + 0.5 * p6)...")
    probs_ensemble = 0.5 * probs_k4 + 0.5 * probs_k6
    metrics_ensemble = calculate_metrics(probs_ensemble, labels)

    # Print results
    print("\n" + "=" * 70)
    print("STAGE 1 ENSEMBLE EVALUATION RESULTS")
    print("=" * 70)

    def print_metrics(name: str, metrics: Dict[str, float]) -> None:
        print(f"\n{name}:")
        print(f"  Accuracy:    {metrics['accuracy']:.4f} ({metrics['accuracy']*100:.2f}%)")
        print(f"  Precision:   {metrics['precision']:.4f}")
        print(f"  Recall:      {metrics['recall']:.4f}")
        print(f"  F1:          {metrics['f1']:.4f}")
        print(f"  MCC:         {metrics['mcc']:.4f}")
        print(f"  ROC-AUC:     {metrics['roc_auc']:.4f}")
        print(f"  Specificity: {metrics['specificity']:.4f}")

    print_metrics("k=4 Model", metrics_k4)
    print_metrics("k=6 Model", metrics_k6)
    print_metrics("Ensemble (k=4 + k=6)", metrics_ensemble)

    # Check if ensemble is better
    best_single = max(metrics_k4["mcc"], metrics_k6["mcc"])
    ensemble_mcc = metrics_ensemble["mcc"]
    if ensemble_mcc >= best_single:
        print(f"\n* Ensemble MCC ({ensemble_mcc:.4f}) >= best single model ({best_single:.4f})")
    else:
        print(f"\n* Ensemble MCC ({ensemble_mcc:.4f}) < best single model ({best_single:.4f})")

    # Build report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report = {
        "timestamp": timestamp,
        "config_path": str(config_path),
        "device": str(device),
        "validation_samples": len(df),
        "models": {
            "k4": {
                "vocab_path": str(vocab_k4_path),
                "checkpoint_path": str(ckpt_k4_path),
                "vocab_size": len(vocab_k4.itos),
                "metrics": metrics_k4,
            },
            "k6": {
                "vocab_path": str(vocab_k6_path),
                "checkpoint_path": str(ckpt_k6_path),
                "vocab_size": len(vocab_k6.itos),
                "metrics": metrics_k6,
            },
        },
        "ensemble": {
            "method": "average",
            "weights": {"k4": 0.5, "k6": 0.5},
            "metrics": metrics_ensemble,
        },
        "comparison": {
            "best_single_mcc": best_single,
            "ensemble_mcc": ensemble_mcc,
            "ensemble_wins": ensemble_mcc >= best_single,
        },
    }

    # Save report to runs directory
    runs_dir = (script_dir / "../runs" / f"stage1_{timestamp}").resolve()
    runs_dir.mkdir(parents=True, exist_ok=True)
    report_path = runs_dir / "stage1_ensemble_report.json"

    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)

    print(f"\nReport saved to: {report_path}")
    print("=" * 70)

    # Also print JSON summary to stdout
    print("\nJSON Summary:")
    summary = {
        "k4_accuracy": metrics_k4["accuracy"],
        "k4_mcc": metrics_k4["mcc"],
        "k4_auroc": metrics_k4["roc_auc"],
        "k6_accuracy": metrics_k6["accuracy"],
        "k6_mcc": metrics_k6["mcc"],
        "k6_auroc": metrics_k6["roc_auc"],
        "ensemble_accuracy": metrics_ensemble["accuracy"],
        "ensemble_mcc": metrics_ensemble["mcc"],
        "ensemble_auroc": metrics_ensemble["roc_auc"],
        "ensemble_wins": ensemble_mcc >= best_single,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
