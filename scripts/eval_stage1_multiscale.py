#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import numpy as np
import pandas as pd
import torch
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
TRAINING_DIR = REPO_ROOT / "gene_whisperer" / "training"
sys.path.insert(0, str(TRAINING_DIR))

from dataset import KmerVocabulary  # noqa: E402
from ensemble_infer import (  # noqa: E402
    build_model,
    compute_engineered_features,
    infer_model_architecture,
    load_checkpoint,
    load_stage1_report_weight,
    select_device,
)

LOGGER = logging.getLogger("gene_whisperer.eval_multiscale")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


@dataclass
class ModelEntry:
    k: int
    vocab: KmerVocabulary
    model: torch.nn.Module


@dataclass
class ProbStats:
    count: int = 0
    total: float = 0.0
    min: float = float("inf")
    max: float = float("-inf")
    below: int = 0
    above: int = 0

    def update(self, probs: torch.Tensor) -> None:
        if not torch.isfinite(probs).all().item():
            raise AssertionError("Non-finite probabilities detected.")
        batch_min = float(probs.min().item())
        batch_max = float(probs.max().item())
        self.min = min(self.min, batch_min)
        self.max = max(self.max, batch_max)
        self.total += float(probs.sum().item())
        self.count += int(probs.numel())
        self.below += int((probs < 0.5).sum().item())
        self.above += int((probs >= 0.5).sum().item())


def _load_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"config.yaml not found at {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"Expected config to be a mapping, got {type(cfg).__name__}")
    return cfg


def _resolve_config_path(config_arg: str, training_dir: Path) -> Path:
    config_path = Path(config_arg).expanduser()
    if not config_path.is_absolute():
        config_path = (Path.cwd() / config_path).resolve()
    if config_path.exists():
        return config_path
    fallback = (training_dir / config_arg).resolve()
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"config.yaml not found at {config_arg}")


def _resolve_path(base_dir: Path, value: str | Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _parse_kmers_arg(value: str) -> list[int]:
    tokens = [item.strip() for item in value.split(",")]
    kmers: list[int] = []
    for token in tokens:
        if not token:
            continue
        try:
            kmers.append(int(token))
        except ValueError as exc:
            raise ValueError(f"Invalid k-mer value: {token!r}") from exc
    if not kmers:
        raise ValueError("Expected at least one k-mer value in --kmers.")
    if len(set(kmers)) != len(kmers):
        raise ValueError(f"Duplicate k-mer values in --kmers: {kmers}")
    return kmers


def _resolve_artifacts_dir(cfg: Dict[str, Any], base_dir: Path) -> Path:
    vocab_cache_dir = cfg.get("vocab_cache_dir")
    if isinstance(vocab_cache_dir, str) and vocab_cache_dir:
        vocab_dir = _resolve_path(base_dir, vocab_cache_dir)
        return vocab_dir.parent
    return (base_dir / "../artifacts").resolve()


def _get_kmers(cfg: Dict[str, Any]) -> list[int]:
    kmers = cfg.get("multi_scale_kmers")
    if not isinstance(kmers, list) or not kmers:
        raise ValueError(f"multi_scale_kmers must be a non-empty list, got {kmers!r}")
    if not all(isinstance(k, int) for k in kmers):
        raise ValueError(f"multi_scale_kmers must be list[int], got {kmers!r}")
    return kmers


def _get_ckpt_path(stage1_ckpt_by_k: Dict[Any, Any], k: int) -> Optional[str]:
    value = stage1_ckpt_by_k.get(k)
    if value is None:
        value = stage1_ckpt_by_k.get(str(k))
    return value


def _find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str:
    for name in candidates:
        if name in df.columns:
            return name
    raise ValueError(f"Missing expected columns. Tried: {', '.join(candidates)}")


def _resolve_ensemble_weights(kmers: list[int], artifacts_dir: Path) -> Optional[list[float]]:
    weights: list[float] = []
    missing: list[int] = []

    for k in kmers:
        report_path = artifacts_dir / f"stage1_report_k{k}.json"
        weight = load_stage1_report_weight(report_path)
        if weight is None:
            missing.append(k)
            weight = 1.0
        weights.append(float(weight))

    if missing:
        LOGGER.warning(
            "Missing validation metrics for k=%s; defaulting weights to 1.0 for those models",
            ",".join(str(k) for k in missing),
        )

    weights = [max(0.0, weight) for weight in weights]
    weight_sum = float(sum(weights))
    if not math.isfinite(weight_sum) or weight_sum <= 0.0:
        LOGGER.warning("Weighted ensemble requested but weights sum to 0; using unweighted mean instead")
        return None

    normalized = [weight / weight_sum for weight in weights]
    if not all(math.isfinite(weight) for weight in normalized):
        LOGGER.warning("Weighted ensemble requested but weights are non-finite; using unweighted mean instead")
        return None

    return normalized


def _build_models(
    cfg: Dict[str, Any],
    base_dir: Path,
    device: torch.device,
) -> list[ModelEntry]:
    stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k")
    if not isinstance(stage1_ckpt_by_k, dict):
        raise ValueError("stage1_ckpt_by_k missing or invalid in config")

    vocab_cache_dir = cfg.get("vocab_cache_dir")
    if not vocab_cache_dir:
        raise ValueError("vocab_cache_dir missing from config")

    vocab_dir = _resolve_path(base_dir, vocab_cache_dir)
    kmers = _get_kmers(cfg)

    models: list[ModelEntry] = []
    missing: list[str] = []

    for k in kmers:
        ckpt_value = _get_ckpt_path(stage1_ckpt_by_k, k)
        if not ckpt_value:
            missing.append(f"k={k} missing stage1_ckpt_by_k entry")
            continue
        ckpt_path = _resolve_path(base_dir, ckpt_value)
        vocab_path = vocab_dir / f"k{k}_vocab.json"

        if not ckpt_path.exists():
            missing.append(f"k={k} missing checkpoint: {ckpt_path}")
            continue
        if not vocab_path.exists():
            missing.append(f"k={k} missing vocab: {vocab_path}")
            continue

        vocab = KmerVocabulary.load(vocab_path)
        arch_overrides = infer_model_architecture(ckpt_path)
        model_cfg = {**cfg, **arch_overrides}
        model = build_model(model_cfg, vocab, device)
        load_checkpoint(model, ckpt_path, device)
        models.append(ModelEntry(k=k, vocab=vocab, model=model))

    if missing:
        raise FileNotFoundError("Missing required artifacts:\n  " + "\n  ".join(missing))

    return models


def _precompute_engineered(sequences: list[str], cfg: Dict[str, Any]) -> torch.Tensor:
    try:
        features = [compute_engineered_features(seq, cfg) for seq in sequences]
        return torch.stack(features, dim=0)
    except MemoryError as exc:
        raise MemoryError("Failed to precompute engineered features; dataset too large.") from exc


def _precompute_tokens(
    sequences: list[str],
    vocab: KmerVocabulary,
    max_bp_len: int,
) -> Optional[torch.Tensor]:
    try:
        tokens = [vocab.tokenize(seq, max_bp_len) for seq in sequences]
        return torch.stack(tokens, dim=0)
    except MemoryError:
        LOGGER.warning("Token cache too large for k=%d; falling back to streaming tokenization.", vocab.k)
        return None


def _compute_metrics_at_threshold(
    probs: np.ndarray,
    labels: np.ndarray,
    threshold: float,
) -> Dict[str, float]:
    preds = (probs >= threshold).astype(int)
    targets = labels.astype(int)

    correct = (preds == targets).sum()
    accuracy = float(correct / len(targets)) if len(targets) else 0.0

    tp = int(((preds == 1) & (targets == 1)).sum())
    tn = int(((preds == 0) & (targets == 0)).sum())
    fp = int(((preds == 1) & (targets == 0)).sum())
    fn = int(((preds == 0) & (targets == 1)).sum())

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / float(np.sqrt(denom)) if denom > 0 else 0.0

    return {
        "accuracy": float(accuracy),
        "f1": float(f1),
        "mcc": float(mcc),
    }


def _best_threshold_accuracy(probs: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    thresholds = np.linspace(0.01, 0.99, 99)
    best_thr = 0.5
    best_acc = -1.0
    targets = labels.astype(int)

    for thr in thresholds:
        preds = (probs >= thr).astype(int)
        accuracy = float((preds == targets).mean()) if len(targets) else 0.0
        if accuracy > best_acc:
            best_acc = accuracy
            best_thr = float(thr)

    return best_thr, best_acc


def _iter_batches(total: int, batch_size: int) -> Iterable[tuple[int, int]]:
    for start in range(0, total, batch_size):
        end = min(total, start + batch_size)
        yield start, end


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate Stage 1 multiscale ensemble on val set")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (defaults to training/config.yaml).",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=None,
        help="Override batch size for evaluation (default: config batch_size).",
    )
    parser.add_argument(
        "--kmers",
        default=None,
        help="Comma-separated list of k-mers to evaluate (overrides config).",
    )
    parser.add_argument(
        "--debug_first_batch",
        action="store_true",
        help="Check ensemble output matches mean of per-k outputs on the first batch.",
    )
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--weighted",
        action="store_true",
        help="Weight ensemble outputs using validation metrics from artifacts/stage1_report_k{k}.json (default).",
    )
    mode_group.add_argument(
        "--unweighted",
        action="store_true",
        help="Use unweighted mean ensemble.",
    )
    args = parser.parse_args()

    if not TRAINING_DIR.exists():
        raise SystemExit(f"Training dir not found: {TRAINING_DIR}")

    config_path = _resolve_config_path(args.config, TRAINING_DIR)
    cfg = _load_config(config_path)
    base_dir = config_path.parent
    if args.kmers:
        cfg["multi_scale_kmers"] = _parse_kmers_arg(args.kmers)

    val_path = cfg.get("stage1_val", "../data/stage1_val.tsv")
    val_path = _resolve_path(base_dir, val_path)
    delimiter = cfg.get("delimiter", "\t")
    if not val_path.exists():
        raise FileNotFoundError(f"Stage1 val set not found: {val_path}")

    LOGGER.info("Loading Stage1 val set: %s", val_path)
    df = pd.read_csv(val_path, delimiter=delimiter)
    if df.empty:
        raise ValueError(f"Stage1 val set is empty: {val_path}")

    seq_col = _find_column(df, ["sequence", "Sequence", "seq", "Seq"])
    label_col = _find_column(df, ["label", "Label", "is_promoter", "isPromoter"])
    sequences = df[seq_col].astype(str).tolist()
    labels = df[label_col].astype(float).to_numpy()

    device = select_device()
    LOGGER.info("Using device: %s", device)

    models = _build_models(cfg, base_dir, device)
    kmers = [entry.k for entry in models]
    LOGGER.info("Loaded %d models: %s", len(models), kmers)

    max_bp_len = int(cfg.get("max_bp_len", 81))
    batch_size = int(args.batch_size or cfg.get("batch_size", 32))

    use_weighted = True
    ensemble_weights: Optional[list[float]] = None
    if args.unweighted:
        use_weighted = False
    if use_weighted:
        artifacts_dir = _resolve_artifacts_dir(cfg, base_dir)
        ensemble_weights = _resolve_ensemble_weights(kmers, artifacts_dir)
        if ensemble_weights is None:
            LOGGER.warning("Weighted ensemble requested but weights unavailable; using unweighted mean instead")
            use_weighted = False
        else:
            weight_summary = ", ".join(
                f"k{k}={weight:.4f}" for k, weight in zip(kmers, ensemble_weights)
            )
            LOGGER.info("Using weighted ensemble (%s)", weight_summary)

    LOGGER.info("Precomputing engineered features...")
    engineered_all = _precompute_engineered(sequences, cfg)

    token_cache: dict[int, Optional[torch.Tensor]] = {}
    for entry in models:
        LOGGER.info("Precomputing tokens for k=%d...", entry.k)
        token_cache[entry.k] = _precompute_tokens(sequences, entry.vocab, max_bp_len)

    total = len(sequences)
    all_probs = np.zeros(total, dtype=np.float32)
    prob_stats: dict[int, ProbStats] = {entry.k: ProbStats() for entry in models}
    debug_checked = False

    LOGGER.info("Running ensemble inference...")
    with torch.no_grad():
        for start, end in _iter_batches(total, batch_size):
            batch_eng = engineered_all[start:end].to(device)
            batch_probs = []
            sum_probs = None
            weighted_sum = None
            batch_sequences = sequences[start:end]

            for idx, entry in enumerate(models):
                cached = token_cache.get(entry.k)
                if cached is None:
                    tokens = [entry.vocab.tokenize(seq, max_bp_len) for seq in batch_sequences]
                    tokens = torch.stack(tokens, dim=0)
                else:
                    tokens = cached[start:end]
                tokens = tokens.to(device)

                logits = entry.model(tokens, engineered_features=batch_eng)
                probs = torch.sigmoid(logits).squeeze(-1)
                prob_stats[entry.k].update(probs)
                batch_probs.append(probs)
                sum_probs = probs if sum_probs is None else (sum_probs + probs)
                if use_weighted and ensemble_weights is not None:
                    weight = ensemble_weights[idx]
                    weighted_sum = probs * weight if weighted_sum is None else (weighted_sum + probs * weight)

            if use_weighted and weighted_sum is not None:
                ensemble_probs = weighted_sum
            else:
                ensemble_probs = sum_probs / len(batch_probs)
            if args.debug_first_batch and not debug_checked:
                manual = torch.stack(batch_probs, dim=0)
                if use_weighted and ensemble_weights is not None:
                    weight_tensor = torch.tensor(
                        ensemble_weights,
                        device=manual.device,
                        dtype=manual.dtype,
                    ).view(-1, 1)
                    manual_mean = (manual * weight_tensor).sum(dim=0)
                else:
                    manual_mean = manual.mean(dim=0)
                max_abs_diff = (manual_mean - ensemble_probs).abs().max().item()
                assert max_abs_diff < 1e-7, (
                    f"Ensemble mismatch on first batch: max_abs_diff={max_abs_diff:.3e}"
                )
                LOGGER.info(
                    "Debug first batch check passed (max_abs_diff=%.3e).",
                    max_abs_diff,
                )
                debug_checked = True
            all_probs[start:end] = ensemble_probs.detach().cpu().numpy()

    metrics_05 = _compute_metrics_at_threshold(all_probs, labels, 0.5)
    best_thr, best_acc = _best_threshold_accuracy(all_probs, labels)

    print("\n" + "=" * 70)
    print("STAGE 1 MULTISCALE ENSEMBLE EVALUATION")
    print("=" * 70)
    print(f"Samples:    {total}")
    print(f"K-mers:     {kmers}")
    print(f"Mode:       {'weighted' if use_weighted else 'unweighted'}")
    if use_weighted and ensemble_weights is not None:
        weight_summary = ", ".join(
            f"k{k}={weight:.4f}" for k, weight in zip(kmers, ensemble_weights)
        )
        print(f"Weights:    {weight_summary}")
    print("-" * 70)
    print("Per-k probability sanity:")
    for k in kmers:
        stats = prob_stats[k]
        if stats.count == 0:
            raise AssertionError(f"No probabilities collected for k={k}")
        if stats.min < 0.0:
            raise AssertionError(f"Probabilities below 0 for k={k}: min={stats.min:.6f}")
        if stats.max > 1.0:
            raise AssertionError(f"Probabilities above 1 for k={k}: max={stats.max:.6f}")
        mean = stats.total / stats.count
        frac_below = stats.below / stats.count
        frac_above = stats.above / stats.count
        print(
            f"k={k}: min={stats.min:.6f} mean={mean:.6f} max={stats.max:.6f} "
            f"frac<0.5={frac_below:.4f} frac>=0.5={frac_above:.4f}"
        )
    print("-" * 70)
    print(f"Accuracy@0.5: {metrics_05['accuracy']:.4f} ({metrics_05['accuracy']*100:.2f}%)")
    print(f"F1@0.5:       {metrics_05['f1']:.4f}")
    print(f"MCC@0.5:      {metrics_05['mcc']:.4f}")
    print(f"Best thr:     {best_thr:.2f}")
    print(f"Best acc:     {best_acc:.4f} ({best_acc*100:.2f}%)")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
