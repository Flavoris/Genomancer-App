#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


def _load_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"config.yaml not found at {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"Expected config to be a mapping, got {type(cfg).__name__}")
    return cfg


def _validate_multiscale(cfg: dict[str, Any]) -> list[int]:
    use_multi_scale = cfg.get("use_multi_scale")
    if use_multi_scale is not True:
        raise ValueError(f"use_multi_scale must be true, got {use_multi_scale!r}")

    multi_scale_kmers = cfg.get("multi_scale_kmers")
    if not isinstance(multi_scale_kmers, list) or not multi_scale_kmers:
        raise ValueError(
            f"multi_scale_kmers must be a non-empty list, got {multi_scale_kmers!r}"
        )
    if not all(isinstance(k, int) for k in multi_scale_kmers):
        raise ValueError(f"multi_scale_kmers must be list[int], got {multi_scale_kmers!r}")
    return multi_scale_kmers


def _read_report_metrics(report_path: Path) -> tuple[float, float]:
    if not report_path.exists():
        raise FileNotFoundError(f"Expected report JSON at {report_path}")
    with report_path.open("r", encoding="utf-8") as handle:
        report = json.load(handle)
    val_metrics = report.get("val_metrics")
    if not isinstance(val_metrics, dict):
        raise ValueError(f"val_metrics missing or invalid in report: {report_path}")

    acc = val_metrics.get("acc_best", val_metrics.get("accuracy"))
    mcc = val_metrics.get("mcc_best", val_metrics.get("mcc"))
    if acc is None or mcc is None:
        raise ValueError(f"Missing val accuracy/MCC in report: {report_path}")
    return float(acc), float(mcc)


def _parse_kmers_arg(kmers_arg: str) -> list[int]:
    items = [item.strip() for item in kmers_arg.split(",") if item.strip()]
    if not items:
        raise ValueError("Expected --kmers to be a comma-separated list of integers")
    kmers: list[int] = []
    for item in items:
        try:
            kmers.append(int(item))
        except ValueError as exc:
            raise ValueError(f"Invalid k-mer value: {item!r}") from exc
    return kmers


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


def _run_stage1(training_dir: Path, config_path: Path, kmer: int) -> None:
    cmd = [sys.executable, "train_stage1.py", "--config", str(config_path)]
    result = subprocess.run(cmd, cwd=training_dir)
    if result.returncode != 0:
        raise RuntimeError(f"Stage1 training failed for k={kmer} (exit {result.returncode})")


def main() -> int:
    parser = argparse.ArgumentParser(description="Train Stage 1 across multiple k-mers")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (defaults to training/config.yaml).",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Only generate per-k configs and print their settings.",
    )
    parser.add_argument(
        "--kmers",
        default=None,
        help="Comma-separated list of k-mers to run (overrides multi_scale_kmers).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    training_dir = repo_root / "gene_whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    config_path = _resolve_config_path(args.config, training_dir)
    base_cfg = _load_config(config_path)
    config_kmers = _validate_multiscale(base_cfg)
    kmers = config_kmers
    if args.kmers is not None:
        kmers = _parse_kmers_arg(args.kmers)
        extra = [k for k in kmers if k not in config_kmers]
        if extra:
            print(
                f"Warning: --kmers includes values not in multi_scale_kmers: {extra}",
                file=sys.stderr,
            )

    artifacts_root = (training_dir / "../artifacts").resolve()
    runs_dir = artifacts_root / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    summaries: list[tuple[int, Path, float, float]] = []
    dry_run = bool(args.dry_run)

    for k in kmers:
        temp_cfg = copy.deepcopy(base_cfg)
        temp_cfg["kmer"] = k
        temp_cfg["stage1_checkpoint_name"] = f"stage1_k{k}.pt"
        mlm_kmer = temp_cfg.get("mlm_kmer")
        if mlm_kmer is None:
            raise ValueError("mlm_kmer missing from config")
        if k != mlm_kmer:
            temp_cfg["stage1_load_mlm_weights"] = False

        temp_cfg_path = runs_dir / f"stage1_multiscale_k{k}.yaml"
        with temp_cfg_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(temp_cfg, handle, sort_keys=False, allow_unicode=False)

        if dry_run:
            print(
                "config={config} k={k} stage1_load_mlm_weights={load} "
                "stage1_checkpoint_name={ckpt}".format(
                    config=temp_cfg_path,
                    k=k,
                    load=temp_cfg.get("stage1_load_mlm_weights"),
                    ckpt=temp_cfg.get("stage1_checkpoint_name"),
                )
            )
            continue

        _run_stage1(training_dir, temp_cfg_path, k)

        report_path = artifacts_root / f"stage1_report_k{k}.json"
        val_acc, val_mcc = _read_report_metrics(report_path)
        checkpoint_path = artifacts_root / "checkpoints" / f"stage1_k{k}.pt"
        summaries.append((k, checkpoint_path, val_acc, val_mcc))

    if dry_run:
        return 0

    print("\nMultiscale Stage 1 Summary")
    for k, checkpoint_path, val_acc, val_mcc in summaries:
        print(
            f"k={k} ckpt={checkpoint_path} val_acc={val_acc:.4f} val_mcc={val_mcc:.4f}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        raise SystemExit(str(exc)) from None
