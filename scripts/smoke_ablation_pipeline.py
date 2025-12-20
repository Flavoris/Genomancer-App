#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any


def _resolve_config_path(config_arg: str, *, training_dir: Path, initial_cwd: Path) -> Path:
    candidate = Path(config_arg)
    if candidate.is_absolute() and candidate.exists():
        return candidate

    by_cwd = (initial_cwd / candidate).resolve()
    if by_cwd.exists():
        return by_cwd

    by_training = (training_dir / candidate).resolve()
    if by_training.exists():
        return by_training

    raise FileNotFoundError(
        f"Config not found: tried {by_cwd} and {by_training} (arg={config_arg!r})"
    )


def _apply_smoke_ablation(cfg: dict[str, Any]) -> None:
    ablation_cfg = cfg.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_cfg["enabled"] = True
    ablation_cfg["pretrain_max_steps"] = 5
    ablation_cfg["stage1_max_steps"] = 5
    ablation_cfg["stage2_max_steps"] = 5
    ablation_cfg["stage1_max_train_samples"] = 64
    ablation_cfg["stage1_max_val_samples"] = 64
    ablation_cfg["stage2_max_train_samples"] = 64
    ablation_cfg["stage2_max_val_samples"] = 64
    cfg["ablation"] = ablation_cfg


def _assert_finite(name: str, value: float) -> None:
    if not math.isfinite(value):
        raise AssertionError(f"{name} must be finite, got {value!r}")


def _get_float(report: dict[str, Any], key: str) -> float:
    value = report.get(key)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _get_required_path(report: dict[str, Any], key: str) -> Path:
    value = str(report.get(key, "")).strip()
    if not value:
        raise AssertionError(f"Missing {key} in report")
    return Path(value)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a tiny end-to-end ablation smoke test.")
    parser.add_argument(
        "--config",
        default="config_smoke.yaml",
        help="Path to YAML config (relative to repo root or training/).",
    )
    args = parser.parse_args()

    initial_cwd = Path.cwd()
    repo_root = Path(__file__).resolve().parents[1]
    training_dir = repo_root / "Gene Whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    sys.path.insert(0, str(training_dir))
    from pretrain_mlm import run_mlm_pretrain  # type: ignore
    from train_stage1 import load_config, run_stage1_training  # type: ignore
    from train_stage2 import run_stage2_training  # type: ignore

    config_path = _resolve_config_path(args.config, training_dir=training_dir, initial_cwd=initial_cwd)
    cfg: dict[str, Any] = load_config(config_path)

    _apply_smoke_ablation(cfg)
    cfg["stage1_load_mlm_weights"] = True

    artifacts_root = (training_dir / "../artifacts").resolve()
    checkpoints_dir = artifacts_root / "checkpoints"
    encoder_ckpt_path = (checkpoints_dir / "mlm_encoder_smoke_pipeline.pt").resolve()
    cfg["mlm_encoder_path"] = str(encoder_ckpt_path)
    cfg["mlm_encoder_checkpoint"] = str(encoder_ckpt_path)
    cfg["stage1_checkpoint_name"] = "stage1_smoke_pipeline.pt"
    cfg["stage2_checkpoint_name"] = "stage2_smoke_pipeline.pt"

    os.chdir(config_path.parent)

    pretrain_report = run_mlm_pretrain(cfg, overrides={"_config_path": str(config_path)})
    encoder_path = _get_required_path(pretrain_report, "encoder_ckpt_path")
    if not encoder_path.exists():
        raise AssertionError(f"Encoder checkpoint missing at {encoder_path}")

    stage1_report = run_stage1_training(
        cfg,
        overrides={
            "_config_path": str(config_path),
            "mlm_encoder_path": str(encoder_path),
        },
    )
    stage1_auc = _get_float(stage1_report, "stage1_auc")
    stage1_mcc = _get_float(stage1_report, "stage1_mcc")
    _assert_finite("stage1_auc", stage1_auc)
    _assert_finite("stage1_mcc", stage1_mcc)

    stage1_ckpt = _get_required_path(stage1_report, "best_ckpt_path")
    if not stage1_ckpt.exists():
        raise AssertionError(f"Stage1 checkpoint missing at {stage1_ckpt}")

    stage2_report = run_stage2_training(
        cfg,
        overrides={
            "_config_path": str(config_path),
            "stage1_ckpt": str(stage1_ckpt),
        },
    )
    stage2_auc = _get_float(stage2_report, "stage2_auc")
    stage2_mcc = _get_float(stage2_report, "stage2_mcc")
    _assert_finite("stage2_auc", stage2_auc)
    _assert_finite("stage2_mcc", stage2_mcc)

    print("SMOKE OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
