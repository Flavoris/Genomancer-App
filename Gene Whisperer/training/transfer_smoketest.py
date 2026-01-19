#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import math
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import torch

from dataset import build_dataloaders
from model import GeneWhispererStage1
from pretrain_mlm import run_mlm_pretrain
from train_stage1 import load_config, run_stage1_training
from train_stage2 import run_stage2_training


SMALL_DELTA_THRESHOLD = 1e-3
# Smaller batches keep the smoketest fast without changing model shapes.
SMOKETEST_MLM_BATCH_SIZE = 32
SMOKETEST_MLM_GRAD_ACCUM = 1
SMOKETEST_STAGE_BATCH_SIZE = 16
SMOKETEST_STAGE_GRAD_ACCUM = 1


def _resolve_config_path(config_arg: str, *, script_dir: Path, initial_cwd: Path) -> Path:
    candidate = Path(config_arg)
    if candidate.is_absolute() and candidate.exists():
        return candidate

    by_cwd = (initial_cwd / candidate).resolve()
    if by_cwd.exists():
        return by_cwd

    by_script = (script_dir / candidate).resolve()
    if by_script.exists():
        return by_script

    fallback = script_dir / "config.yaml"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        f"Config not found: tried {by_cwd} and {by_script} (arg={config_arg!r})"
    )


def _coerce_positive_int(value: Any, *, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be an int, got {value!r}") from exc
    if parsed <= 0:
        raise ValueError(f"{label} must be > 0, got {parsed}")
    return parsed


def _get_smoketest_steps(cfg: dict, key: str, fallback_ablation_key: str, default: int) -> int:
    raw = cfg.get(key)
    if raw is None:
        ablation_cfg = cfg.get("ablation")
        if isinstance(ablation_cfg, dict):
            raw = ablation_cfg.get(fallback_ablation_key)
    if raw is None:
        raw = default
    return _coerce_positive_int(raw, label=key)


def _apply_pretrain_ablation(cfg: dict, *, steps: int) -> None:
    ablation_cfg = cfg.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_cfg["enabled"] = True
    ablation_cfg["pretrain_max_steps"] = steps
    existing_epochs = ablation_cfg.get("pretrain_epochs", 1)
    try:
        existing_epochs_int = int(existing_epochs)
    except (TypeError, ValueError):
        existing_epochs_int = 1
    ablation_cfg["pretrain_epochs"] = max(existing_epochs_int, steps)
    cfg["ablation"] = ablation_cfg


def _apply_stage_ablation(cfg: dict, *, stage: str, steps: int) -> None:
    ablation_cfg = cfg.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_cfg["enabled"] = True
    ablation_cfg[f"{stage}_max_steps"] = steps
    existing_epochs = ablation_cfg.get(f"{stage}_epochs", 1)
    try:
        existing_epochs_int = int(existing_epochs)
    except (TypeError, ValueError):
        existing_epochs_int = 1
    ablation_cfg[f"{stage}_epochs"] = max(existing_epochs_int, steps)
    cfg["ablation"] = ablation_cfg


def _load_checkpoint_state(path: Path) -> Dict[str, Any]:
    state = torch.load(path, map_location="cpu", weights_only=False)
    if isinstance(state, dict) and "model_state" in state:
        return state["model_state"]
    if isinstance(state, dict):
        return state
    raise ValueError(f"Unsupported checkpoint format at {path} (type={type(state)})")


def _build_stage1_model(
    cfg: dict, *, engineered_dim: int, use_engineered: bool
) -> GeneWhispererStage1:
    pad_token_id = cfg.get("pad_token_id")
    pad_token_id = None if pad_token_id is None else int(pad_token_id)
    return GeneWhispererStage1(
        vocab_size=int(cfg.get("vocab_size", 67)),
        kmer=int(cfg.get("kmer", 3)),
        embedding_dim=int(cfg.get("embedding_dim", 256)),
        num_layers=int(cfg.get("transformer_layers", 6)),
        num_heads=int(cfg.get("transformer_heads", 8)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1024)),
        dropout=float(cfg.get("transformer_dropout", 0.15)),
        pad_token_id=pad_token_id,
        engineered_dim=engineered_dim,
        use_engineered_features=use_engineered,
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
    )


def _first_batch_from_loader(loader) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
    batch = next(iter(loader))
    if isinstance(batch, (tuple, list)) and len(batch) == 3:
        tokens, engineered, _ = batch
        return tokens, engineered
    if isinstance(batch, (tuple, list)) and len(batch) == 4:
        tokens, tnc, pseeiip, _ = batch
        engineered = torch.cat([tnc, pseeiip], dim=1)
        return tokens, engineered
    raise ValueError(f"Unexpected batch format: {type(batch)} len={getattr(batch, '__len__', None)}")


def _compare_first_batch_logits(
    cfg: dict,
    *,
    stage: str,
    pretrained_ckpt: Path,
    scratch_ckpt: Path,
    use_engineered: bool,
    engineered_dim: int,
) -> None:
    loaders = build_dataloaders(cfg)
    stage_loaders = loaders.get(stage) if isinstance(loaders, dict) else None
    if not stage_loaders:
        print(f"{stage}_logit_check: missing {stage} dataloaders")
        return
    loader = stage_loaders.get("val") or stage_loaders.get("train")
    if loader is None:
        print(f"{stage}_logit_check: no loader available")
        return

    tokens, engineered = _first_batch_from_loader(loader)
    tokens = tokens.to(dtype=torch.long)
    if not use_engineered:
        engineered = None

    model_pre = _build_stage1_model(cfg, engineered_dim=engineered_dim, use_engineered=use_engineered)
    model_scratch = _build_stage1_model(cfg, engineered_dim=engineered_dim, use_engineered=use_engineered)
    model_pre.load_state_dict(_load_checkpoint_state(pretrained_ckpt), strict=False)
    model_scratch.load_state_dict(_load_checkpoint_state(scratch_ckpt), strict=False)
    model_pre.eval()
    model_scratch.eval()

    with torch.no_grad():
        if engineered is None:
            out_pre = model_pre(tokens)
            out_scratch = model_scratch(tokens)
        else:
            engineered = engineered.to(dtype=torch.float32)
            out_pre = model_pre(tokens, engineered)
            out_scratch = model_scratch(tokens, engineered)

    diff = (out_pre - out_scratch).abs()
    mean_diff = float(diff.mean().item())
    max_diff = float(diff.max().item())
    equal_flag = mean_diff < 1e-6
    print(
        f"{stage}_first_batch_logit_mean_abs_diff={mean_diff:.6f} "
        f"max_abs_diff={max_diff:.6f} equal={equal_flag}"
    )


def _print_delta(label: str, pretrained: float, scratch: float) -> float:
    delta = pretrained - scratch
    print(
        f"{label}_pretrained={pretrained:.6f} "
        f"{label}_scratch={scratch:.6f} "
        f"delta={delta:.6f}"
    )
    return delta


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a tiny transfer smoketest (MLM -> Stage1/Stage2).")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (relative to training/ or current working directory).",
    )
    args = parser.parse_args()

    initial_cwd = Path.cwd()
    script_dir = Path(__file__).resolve().parent
    config_path = _resolve_config_path(args.config, script_dir=script_dir, initial_cwd=initial_cwd)
    cfg: dict[str, Any] = load_config(config_path)

    pretrain_steps = _get_smoketest_steps(cfg, "smoketest_steps_pretrain", "pretrain_max_steps", 300)
    stage1_steps = _get_smoketest_steps(cfg, "smoketest_steps_stage1", "stage1_max_steps", 300)
    stage2_steps = _get_smoketest_steps(cfg, "smoketest_steps_stage2", "stage2_max_steps", 300)

    artifacts_root = (script_dir / "../artifacts").resolve()
    checkpoints_dir = artifacts_root / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    encoder_ckpt = (checkpoints_dir / "mlm_encoder_transfer_smoketest.pt").resolve()

    os.chdir(config_path.parent)

    cfg_pretrain = copy.deepcopy(cfg)
    cfg_pretrain["ablation_disable_wandb"] = True
    cfg_pretrain["mlm_encoder_path"] = str(encoder_ckpt)
    cfg_pretrain["mlm_encoder_checkpoint"] = str(encoder_ckpt)
    cfg_pretrain["mlm_batch_size"] = min(
        int(cfg_pretrain.get("mlm_batch_size", SMOKETEST_MLM_BATCH_SIZE)),
        SMOKETEST_MLM_BATCH_SIZE,
    )
    cfg_pretrain["mlm_grad_accum_steps"] = SMOKETEST_MLM_GRAD_ACCUM
    _apply_pretrain_ablation(cfg_pretrain, steps=pretrain_steps)
    pretrain_report = run_mlm_pretrain(cfg_pretrain, overrides={"_config_path": str(config_path)})
    encoder_path = Path(str(pretrain_report.get("encoder_ckpt_path", ""))).resolve()
    if not encoder_path.exists():
        raise FileNotFoundError(f"Encoder checkpoint missing at {encoder_path}")

    cfg_stage1_pre = copy.deepcopy(cfg)
    cfg_stage1_pre["ablation_disable_wandb"] = True
    cfg_stage1_pre["stage1_load_mlm_weights"] = True
    cfg_stage1_pre["mlm_encoder_path"] = str(encoder_path)
    cfg_stage1_pre["mlm_encoder_checkpoint"] = str(encoder_path)
    cfg_stage1_pre["stage1_checkpoint_name"] = "stage1_transfer_smoketest_pretrained.pt"
    cfg_stage1_pre["batch_size"] = min(
        int(cfg_stage1_pre.get("batch_size", SMOKETEST_STAGE_BATCH_SIZE)),
        SMOKETEST_STAGE_BATCH_SIZE,
    )
    cfg_stage1_pre["grad_accum_steps"] = SMOKETEST_STAGE_GRAD_ACCUM
    _apply_stage_ablation(cfg_stage1_pre, stage="stage1", steps=stage1_steps)

    stage1_pre_report = run_stage1_training(cfg_stage1_pre, overrides={"_config_path": str(config_path)})
    stage1_mcc_pre = float(stage1_pre_report.get("stage1_mcc", float("nan")))
    stage1_pre_ckpt = Path(str(stage1_pre_report.get("best_ckpt_path", ""))).resolve()
    if not stage1_pre_ckpt.exists():
        raise FileNotFoundError(f"Stage1 pretrained checkpoint missing at {stage1_pre_ckpt}")

    cfg_stage1_scratch = copy.deepcopy(cfg)
    cfg_stage1_scratch["ablation_disable_wandb"] = True
    cfg_stage1_scratch["stage1_load_mlm_weights"] = False
    cfg_stage1_scratch["stage1_mlm_transfer_mode"] = "none"
    cfg_stage1_scratch["stage1_checkpoint_name"] = "stage1_transfer_smoketest_scratch.pt"
    cfg_stage1_scratch["batch_size"] = min(
        int(cfg_stage1_scratch.get("batch_size", SMOKETEST_STAGE_BATCH_SIZE)),
        SMOKETEST_STAGE_BATCH_SIZE,
    )
    cfg_stage1_scratch["grad_accum_steps"] = SMOKETEST_STAGE_GRAD_ACCUM
    _apply_stage_ablation(cfg_stage1_scratch, stage="stage1", steps=stage1_steps)

    stage1_scratch_report = run_stage1_training(cfg_stage1_scratch, overrides={"_config_path": str(config_path)})
    stage1_mcc_scratch = float(stage1_scratch_report.get("stage1_mcc", float("nan")))
    stage1_scratch_ckpt = Path(str(stage1_scratch_report.get("best_ckpt_path", ""))).resolve()
    if not stage1_scratch_ckpt.exists():
        raise FileNotFoundError(f"Stage1 scratch checkpoint missing at {stage1_scratch_ckpt}")

    stage1_delta = _print_delta("stage1_mcc", stage1_mcc_pre, stage1_mcc_scratch)
    if math.isfinite(stage1_delta) and abs(stage1_delta) < SMALL_DELTA_THRESHOLD:
        use_engineered_stage1 = bool(cfg_stage1_pre.get("stage1_use_engineered_features", True))
        _compare_first_batch_logits(
            cfg_stage1_pre,
            stage="stage1",
            pretrained_ckpt=stage1_pre_ckpt,
            scratch_ckpt=stage1_scratch_ckpt,
            use_engineered=use_engineered_stage1,
            engineered_dim=int(cfg_stage1_pre.get("engineered_dim", 288)),
        )

    cfg_stage2_pre = copy.deepcopy(cfg)
    cfg_stage2_pre["ablation_disable_wandb"] = True
    cfg_stage2_pre["stage2_checkpoint_name"] = "stage2_transfer_smoketest_pretrained.pt"
    cfg_stage2_pre["batch_size"] = min(
        int(cfg_stage2_pre.get("batch_size", SMOKETEST_STAGE_BATCH_SIZE)),
        SMOKETEST_STAGE_BATCH_SIZE,
    )
    cfg_stage2_pre["grad_accum_steps"] = SMOKETEST_STAGE_GRAD_ACCUM
    _apply_stage_ablation(cfg_stage2_pre, stage="stage2", steps=stage2_steps)

    stage2_pre_report = run_stage2_training(
        cfg_stage2_pre,
        overrides={
            "_config_path": str(config_path),
            "stage1_ckpt": str(stage1_pre_ckpt),
        },
    )
    stage2_mcc_pre = float(stage2_pre_report.get("stage2_mcc", float("nan")))
    stage2_pre_ckpt = Path(str(stage2_pre_report.get("best_ckpt_path", ""))).resolve()
    if not stage2_pre_ckpt.exists():
        raise FileNotFoundError(f"Stage2 pretrained checkpoint missing at {stage2_pre_ckpt}")

    cfg_stage2_scratch = copy.deepcopy(cfg)
    cfg_stage2_scratch["ablation_disable_wandb"] = True
    cfg_stage2_scratch["stage2_checkpoint_name"] = "stage2_transfer_smoketest_scratch.pt"
    cfg_stage2_scratch["batch_size"] = min(
        int(cfg_stage2_scratch.get("batch_size", SMOKETEST_STAGE_BATCH_SIZE)),
        SMOKETEST_STAGE_BATCH_SIZE,
    )
    cfg_stage2_scratch["grad_accum_steps"] = SMOKETEST_STAGE_GRAD_ACCUM
    _apply_stage_ablation(cfg_stage2_scratch, stage="stage2", steps=stage2_steps)

    stage2_scratch_report = run_stage2_training(
        cfg_stage2_scratch,
        overrides={
            "_config_path": str(config_path),
            "stage1_ckpt": str(stage1_scratch_ckpt),
        },
    )
    stage2_mcc_scratch = float(stage2_scratch_report.get("stage2_mcc", float("nan")))
    stage2_scratch_ckpt = Path(str(stage2_scratch_report.get("best_ckpt_path", ""))).resolve()
    if not stage2_scratch_ckpt.exists():
        raise FileNotFoundError(f"Stage2 scratch checkpoint missing at {stage2_scratch_ckpt}")

    stage2_delta = _print_delta("stage2_mcc", stage2_mcc_pre, stage2_mcc_scratch)
    if math.isfinite(stage2_delta) and abs(stage2_delta) < SMALL_DELTA_THRESHOLD:
        use_engineered_stage2 = bool(cfg_stage2_pre.get("stage2_use_engineered_features", True))
        _compare_first_batch_logits(
            cfg_stage2_pre,
            stage="stage2",
            pretrained_ckpt=stage2_pre_ckpt,
            scratch_ckpt=stage2_scratch_ckpt,
            use_engineered=use_engineered_stage2,
            engineered_dim=int(cfg_stage2_pre.get("engineered_dim", 288)),
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
