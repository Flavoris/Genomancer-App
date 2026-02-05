#!/usr/bin/env python3
from __future__ import annotations

__test__ = False

import argparse
import json
import math
from pathlib import Path
from typing import Any

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
TRAINING_DIR = REPO_ROOT / "gene_whisperer" / "training"


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


def _resolve_artifacts_dir(cfg: dict[str, Any], base_dir: Path) -> Path:
    vocab_cache_dir = cfg.get("vocab_cache_dir")
    if isinstance(vocab_cache_dir, str) and vocab_cache_dir:
        vocab_dir = _resolve_path(base_dir, vocab_cache_dir)
        return vocab_dir.parent
    return (base_dir / "../artifacts").resolve()


def _get_kmers(cfg: dict[str, Any]) -> list[int]:
    kmers = cfg.get("multi_scale_kmers")
    if not isinstance(kmers, list) or not kmers:
        raise ValueError(f"multi_scale_kmers must be a non-empty list, got {kmers!r}")
    if not all(isinstance(k, int) for k in kmers):
        raise ValueError(f"multi_scale_kmers must be list[int], got {kmers!r}")
    return kmers


def _load_report(report_path: Path) -> dict[str, Any]:
    if not report_path.exists():
        raise FileNotFoundError(f"Stage 1 report not found: {report_path}")
    try:
        with report_path.open("r", encoding="utf-8") as handle:
            report = json.load(handle)
    except (json.JSONDecodeError, OSError) as exc:
        raise ValueError(f"Failed to read report: {report_path}") from exc
    if not isinstance(report, dict):
        raise ValueError(f"Expected report mapping at {report_path}, got {type(report).__name__}")
    return report


def _extract_weight(report: dict[str, Any]) -> tuple[float, str]:
    val_metrics = report.get("val_metrics") or {}
    for key in ("mcc", "mcc_best", "accuracy", "f1"):
        if key in val_metrics:
            return float(val_metrics[key]), f"val_metrics.{key}"

    for key in (
        "val_mcc",
        "val_mcc_best",
        "best_val_mcc",
        "stage1_mcc",
        "val_acc@best",
        "val_accuracy",
        "val_acc_best",
    ):
        if key in report:
            return float(report[key]), key

    raise KeyError("No validation metric found for weighting")


def main() -> int:
    parser = argparse.ArgumentParser(description="Test weight normalization for multiscale ensemble")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (defaults to training/config.yaml).",
    )
    parser.add_argument(
        "--require-reports",
        action="store_true",
        help="Fail if any stage1_report_k{k}.json files are missing.",
    )
    args = parser.parse_args()

    if not TRAINING_DIR.exists():
        raise SystemExit(f"Training dir not found: {TRAINING_DIR}")

    config_path = _resolve_config_path(args.config, TRAINING_DIR)
    cfg = _load_config(config_path)
    base_dir = config_path.parent

    kmers = _get_kmers(cfg)
    artifacts_dir = _resolve_artifacts_dir(cfg, base_dir)
    if not artifacts_dir.exists():
        raise FileNotFoundError(f"Artifacts dir not found: {artifacts_dir}")

    raw_weights: list[float] = []
    sources: list[str] = []

    missing: list[int] = []
    for k in kmers:
        report_path = artifacts_dir / f"stage1_report_k{k}.json"
        if not report_path.exists():
            if args.require_reports:
                raise FileNotFoundError(f"Stage 1 report not found: {report_path}")
            missing.append(k)
            raw_weights.append(1.0)
            sources.append("default_missing")
            continue
        report = _load_report(report_path)
        weight, source = _extract_weight(report)
        raw_weights.append(weight)
        sources.append(source)

    clipped_weights = [max(0.0, weight) for weight in raw_weights]
    if not all(math.isfinite(weight) for weight in clipped_weights):
        raise AssertionError(f"Non-finite weight detected: {clipped_weights}")

    weight_sum = float(sum(clipped_weights))
    if not math.isfinite(weight_sum) or weight_sum <= 0.0:
        raise AssertionError(f"Weights sum to {weight_sum}; cannot normalize")

    weights = [weight / weight_sum for weight in clipped_weights]
    if any(weight < 0.0 for weight in weights):
        raise AssertionError(f"Negative weight after clipping: {weights}")

    final_sum = float(sum(weights))
    if abs(final_sum - 1.0) > 1e-6:
        raise AssertionError(f"sum(weights)={final_sum:.10f} (expected 1.0)")

    if not all(math.isfinite(weight) for weight in weights):
        raise AssertionError(f"Non-finite normalized weight detected: {weights}")

    print("Weight normalization check")
    if missing:
        missing_str = ",".join(str(k) for k in missing)
        print(
            "WARNING: missing stage1_report_k{{k}}.json for k={missing}; defaulting to 1.0".format(
                missing=missing_str
            )
        )
    for k, raw, clipped, normalized, source in zip(kmers, raw_weights, clipped_weights, weights, sources):
        print(
            "k={k} raw={raw:.6f} clipped={clipped:.6f} weight={normalized:.6f} source={source}".format(
                k=k,
                raw=raw,
                clipped=clipped,
                normalized=normalized,
                source=source,
            )
        )
    print(f"sum(weights)={final_sum:.10f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
