"""Training script for Gene Whisperer Stage 2 classifier (strong vs weak promoters).

Stage 2 reuses the Stage 1 backbone + binary head, training on the Stage 2 TSV
(`sequence`, `strength`) via the stage2 dataloaders.

Usage:
    python train_stage2.py --config config.yaml
    python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import torch
import torch.nn as nn
from torch.optim.swa_utils import AveragedModel, SWALR

from amp_utils import get_amp_context, log_amp_status
from dataset import build_dataloaders
from model import GeneWhispererStage1, GeneWhispererStage2

import train_stage1 as stage1

LOGGER = logging.getLogger("gene_whisperer.stage2")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Ensure shared helpers (run_epoch, save_checkpoint, etc.) log under stage2.
stage1.LOGGER = LOGGER


class _Stage2LoaderAdapter:
    """Adapter that converts Stage 2 batches to the (tokens, engineered, labels) triple."""

    def __init__(self, loader):
        self._loader = loader
        self.dataset = getattr(loader, "dataset", None)

    def __len__(self) -> int:  # pragma: no cover - delegated behavior
        return len(self._loader)

    def __iter__(self):
        for batch in self._loader:
            # Accept either Stage1-like (tokens, engineered, label) or Stage2-like batches.
            if isinstance(batch, (tuple, list)) and len(batch) == 3:
                tokens, engineered, labels = batch
            elif isinstance(batch, (tuple, list)) and len(batch) == 4:
                tokens, tnc, pseeiip, labels = batch
                engineered = torch.cat([tnc, pseeiip], dim=1)
            else:
                raise ValueError(f"Unexpected Stage 2 batch format: {type(batch)} len={getattr(batch, '__len__', None)}")
            yield tokens, engineered, labels


def load_config(path: Path) -> Dict:
    return stage1.load_config(path)


def _resolve_optional_path(value: Any, *, base_dir: Path) -> Optional[Path]:
    if value is None:
        return None
    if isinstance(value, Path):
        candidate = value
    else:
        text = str(value).strip()
        if not text:
            return None
        candidate = Path(text)
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    return candidate


def _normalize_state_dict(state: dict[str, Any]) -> dict[str, Any]:
    if not state:
        return state
    if any(k.startswith("module.") for k in state):
        return {k.removeprefix("module."): v for k, v in state.items()}
    return state


def _load_checkpoint_compatible(model: nn.Module, checkpoint_path: Path, *, device: torch.device) -> int:
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    if isinstance(checkpoint, dict) and "model_state" in checkpoint:
        state = checkpoint["model_state"]
    elif isinstance(checkpoint, dict):
        state = checkpoint
    else:
        raise ValueError(f"Unsupported checkpoint format: {type(checkpoint)}")

    state = _normalize_state_dict(state)
    model_state = model.state_dict()

    compatible: dict[str, Any] = {}
    skipped: list[str] = []
    for key, value in state.items():
        target = model_state.get(key)
        if target is None:
            skipped.append(key)
            continue
        if getattr(target, "shape", None) != getattr(value, "shape", None):
            skipped.append(key)
            continue
        compatible[key] = value

    missing = [k for k in model_state.keys() if k not in compatible]
    model.load_state_dict(compatible, strict=False)

    loaded = len(compatible)
    LOGGER.info(
        "Loaded %d compatible tensors from %s (skipped=%d, missing=%d)",
        loaded,
        checkpoint_path,
        len(skipped),
        len(missing),
    )
    if skipped:
        preview = ", ".join(skipped[:5])
        suffix = "…" if len(skipped) > 5 else ""
        LOGGER.info("Skipped keys (preview): %s%s", preview, suffix)
    if loaded == 0:
        LOGGER.warning("No compatible tensors loaded from %s", checkpoint_path)

    return loaded


def main() -> None:
    parser = argparse.ArgumentParser(description="Train Stage 2 (strong vs weak) classifier")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument(
        "--stage1_ckpt",
        default=None,
        help="Optional Stage 1 checkpoint to initialize from (relative to training dir unless absolute).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility. Overrides config value.",
    )
    parser.add_argument(
        "--kmer",
        type=int,
        default=None,
        help="K-mer size (e.g. 3, 4, 6). Overrides stage2_kmer or kmer config value.",
    )
    parser.add_argument(
        "--overfit_debug",
        action="store_true",
        help="Train on 32 samples for 200 steps to verify model can learn",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Enable torch.profiler for ~50 steps and export chrome trace",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()

    if not config_path.exists():
        fallback = script_dir / "config.yaml"
        if fallback.exists():
            config_path = fallback
            LOGGER.warning("Using fallback config: %s", config_path)

    cfg = load_config(config_path)

    overrides: Dict[str, Any] = {"_config_path": str(config_path)}
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.stage1_ckpt:
        overrides["stage1_ckpt"] = args.stage1_ckpt
    if args.kmer is not None:
        overrides["stage2_kmer"] = args.kmer
    if args.overfit_debug:
        overrides["overfit_debug"] = True
    if args.profile:
        overrides["profile"] = True

    _ = run_stage2_training(cfg, overrides=overrides)


def run_stage2_training(cfg: dict, *, overrides: dict | None = None) -> dict:
    """Trains Stage2 and returns best val metrics + timing + params."""
    start_time = time.perf_counter()

    cfg_run = copy.copy(cfg) if isinstance(cfg, dict) else dict(cfg or {})
    if overrides:
        cfg_run.update(overrides)

    if bool(cfg_run.get("ablation_disable_wandb", False)):
        os.environ["WANDB_MODE"] = "disabled"
        os.environ["WANDB_DISABLED"] = "true"

    script_dir = Path(__file__).resolve().parent
    config_path_value = cfg_run.get("_config_path") or cfg_run.get("config_path")
    config_id = str(config_path_value) if config_path_value else "<in-memory>"

    ablation_cfg = cfg_run.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_enabled = bool(ablation_cfg.get("enabled", False))

    epochs = int(cfg_run.get("epochs", 120))
    if ablation_enabled and ablation_cfg.get("stage2_epochs") is not None:
        epochs_raw = ablation_cfg.get("stage2_epochs")
        try:
            epochs = int(epochs_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.stage2_epochs must be an int, got {epochs_raw!r}") from exc
        if epochs <= 0:
            raise ValueError(f"ablation.stage2_epochs must be > 0, got {epochs}")

    ablation_max_steps: Optional[int] = None
    if ablation_enabled and ablation_cfg.get("stage2_max_steps") is not None:
        max_steps_raw = ablation_cfg.get("stage2_max_steps")
        try:
            max_steps_int = int(max_steps_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.stage2_max_steps must be an int, got {max_steps_raw!r}") from exc
        if max_steps_int < 0:
            raise ValueError(f"ablation.stage2_max_steps must be >= 0, got {max_steps_int}")
        if max_steps_int > 0:
            ablation_max_steps = max_steps_int
    else:
        ablation_max_steps_raw = cfg_run.get("ablation_max_steps_stage2")
        if ablation_max_steps_raw is not None:
            try:
                max_steps_int = int(ablation_max_steps_raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"ablation_max_steps_stage2 must be an int, got {ablation_max_steps_raw!r}"
                ) from exc
            if max_steps_int < 0:
                raise ValueError(f"ablation_max_steps_stage2 must be >= 0, got {max_steps_int}")
            if max_steps_int > 0:
                ablation_max_steps = max_steps_int

    stage2_max_train_samples = ablation_cfg.get("stage2_max_train_samples") if ablation_enabled else None
    stage2_max_val_samples = ablation_cfg.get("stage2_max_val_samples") if ablation_enabled else None

    # Resolve kmer for Stage 2: use stage2_kmer if present, else fall back to kmer
    stage2_kmer = cfg_run.get("stage2_kmer")
    if stage2_kmer is not None:
        stage2_kmer = int(stage2_kmer)
    else:
        stage2_kmer = int(cfg_run.get("kmer", 3))
    # Update cfg_run so build_dataloaders uses the correct kmer for Stage 2
    cfg_run["kmer"] = stage2_kmer
    LOGGER.info("Stage 2 using k-mer size: %d", stage2_kmer)

    # Helper to get stage2-specific config with fallback to global value when null
    def get_with_fallback(stage_key: str, global_key: str, default: Any) -> Any:
        """Get stage2-specific value, falling back to global if null/missing."""
        val = cfg_run.get(stage_key)
        if val is not None:
            return val
        return cfg_run.get(global_key, default)

    # =========================================================================
    # Compute "effective" Stage 2 hyperparameters with fallback to global values
    # =========================================================================
    stage2_label_smoothing = float(get_with_fallback("stage2_label_smoothing", "label_smoothing", 0.05))
    stage2_use_mixup = bool(get_with_fallback("stage2_use_mixup", "use_mixup", True))
    stage2_mixup_alpha = float(get_with_fallback("stage2_mixup_alpha", "mixup_alpha", 0.2))
    stage2_use_swa = bool(get_with_fallback("stage2_use_swa", "use_swa", True))
    stage2_swa_start_epoch = int(get_with_fallback("stage2_swa_start_epoch", "swa_start_epoch", 60))
    stage2_swa_lr = float(get_with_fallback("stage2_swa_lr", "swa_lr", 5e-5))
    stage2_drop_path_rate = float(get_with_fallback("stage2_drop_path_rate", "drop_path_rate", 0.1))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"stage2_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    # Resolve Stage 1 checkpoint path with auto-selection by k-mer
    # Priority: 1) --stage1_ckpt CLI arg  2) stage1_ckpt_by_k[stage2_kmer]  3) None (fall back to MLM weights)
    stage1_ckpt_path = _resolve_optional_path(cfg_run.get("stage1_ckpt"), base_dir=script_dir)
    stage1_ckpt_auto_selected = False
    if stage1_ckpt_path is None:
        # Auto-select from stage1_ckpt_by_k mapping
        stage1_ckpt_by_k = cfg_run.get("stage1_ckpt_by_k")
        if isinstance(stage1_ckpt_by_k, dict):
            # Keys may be int or str depending on YAML parsing
            ckpt_for_k = stage1_ckpt_by_k.get(stage2_kmer) or stage1_ckpt_by_k.get(str(stage2_kmer))
            if ckpt_for_k:
                stage1_ckpt_path = _resolve_optional_path(ckpt_for_k, base_dir=script_dir)
                if stage1_ckpt_path and stage1_ckpt_path.exists():
                    stage1_ckpt_auto_selected = True
                    LOGGER.info("Auto-selected Stage 1 checkpoint for k=%d: %s", stage2_kmer, stage1_ckpt_path)
                elif stage1_ckpt_path:
                    LOGGER.warning(
                        "stage1_ckpt_by_k[%d] points to non-existent file: %s; will fall back to MLM weights",
                        stage2_kmer,
                        stage1_ckpt_path,
                    )
                    stage1_ckpt_path = None
        else:
            LOGGER.debug("No stage1_ckpt_by_k mapping configured; will use MLM weights if available")

    resolved_config: Dict[str, Any] = {
        "script": "train_stage2.py",
        "timestamp": timestamp,
        "config_path": config_id,
        "stage1_ckpt": str(stage1_ckpt_path) if stage1_ckpt_path else None,
        "stage1_ckpt_auto_selected": stage1_ckpt_auto_selected,
        # System info for reproducibility
        "system_info": {
            "python_version": f"{os.sys.version_info.major}.{os.sys.version_info.minor}.{os.sys.version_info.micro}",
            "torch_version": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
            "mps_available": hasattr(torch.backends, "mps") and torch.backends.mps.is_available(),
        },
        # Ablations / debugging
        "ablation_enabled": ablation_enabled,
        "ablation_max_steps_stage2": ablation_max_steps,
        "ablation_stage2_epochs": epochs if ablation_enabled else None,
        "ablation_stage2_max_train_samples": stage2_max_train_samples if ablation_enabled else None,
        "ablation_stage2_max_val_samples": stage2_max_val_samples if ablation_enabled else None,
        "ablation_disable_wandb": bool(cfg_run.get("ablation_disable_wandb", False)),
        "stage2_use_engineered_features": bool(cfg_run.get("stage2_use_engineered_features", True)),
        # Model architecture
        "embedding_dim": int(cfg_run.get("embedding_dim", 256)),
        "transformer_layers": int(cfg_run.get("transformer_layers", 6)),
        "transformer_heads": int(cfg_run.get("transformer_heads", 8)),
        "transformer_ff_dim": int(cfg_run.get("transformer_ff_dim", 1024)),
        "transformer_dropout": float(cfg_run.get("transformer_dropout", 0.15)),
        # Optimizer
        "lr": float(cfg_run.get("lr", 2e-4)),
        "weight_decay": float(cfg_run.get("weight_decay", 0.02)),
        # Data
        "batch_size": int(cfg_run.get("batch_size", 32)),
        "max_bp_len": int(cfg_run.get("max_bp_len", 81)),
        "kmer": stage2_kmer,  # Resolved stage2_kmer (from --kmer, stage2_kmer, or kmer fallback)
        # Engineered features
        "engineered_dim": int(cfg_run.get("engineered_dim", 288)),
        # Architecture toggles
        "use_attention_pool": bool(cfg_run.get("use_attention_pool", True)),
        "use_tcn": bool(cfg_run.get("use_tcn", True)),
        "post_cnn_transformer_layers": int(cfg_run.get("post_cnn_transformer_layers", 3)),
        "stage2_mlm_transfer_mode": str(cfg_run.get("stage2_mlm_transfer_mode", "embed_only")),
        # Stage 2 head configuration
        "stage2_head_type": str(cfg_run.get("stage2_head_type", "transformer")),
        "stage2_transformer_layers": int(cfg_run.get("stage2_transformer_layers", 2)),
        "stage2_transformer_heads": int(cfg_run.get("stage2_transformer_heads", 4)),
        "stage2_bilstm_hidden": int(cfg_run.get("stage2_bilstm_hidden", 128)),
        "stage2_bilstm_layers": int(cfg_run.get("stage2_bilstm_layers", 1)),
        "stage2_combine_mode": str(cfg_run.get("stage2_combine_mode", "concat")),
        # Seed
        "seed": int(cfg_run.get("seed", 1337) or 1337),
        # AMP (Mixed Precision)
        "amp_enabled": bool(cfg_run.get("amp_enabled", True)),
        "amp_dtype": str(cfg_run.get("amp_dtype", "bfloat16")),
        "amp_force_fp32_loss": bool(cfg_run.get("amp_force_fp32_loss", True)),
    }

    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as handle:
        json.dump(resolved_config, handle, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)

    seed = int(cfg_run.get("seed", 1337) or 1337)
    stage1.set_seed(seed)
    LOGGER.info("Using seed: %d", seed)

    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Using device: %s", device)

    # AMP (Mixed Precision) setup
    amp_ctx = get_amp_context(device, cfg_run, logger=LOGGER)
    log_amp_status(amp_ctx, logger=LOGGER)

    # GradScaler for CUDA fp16 (not needed for MPS bfloat16)
    scaler = None
    if amp_ctx.enabled and device.type == "cuda" and amp_ctx.dtype == torch.float16:
        scaler = torch.cuda.amp.GradScaler()
        LOGGER.info("Using GradScaler for CUDA fp16")

    dataloaders = build_dataloaders(cfg_run)
    train_loader_raw = dataloaders["stage2"]["train"]
    val_loader_raw = dataloaders["stage2"]["val"]

    train_loader = _Stage2LoaderAdapter(train_loader_raw)
    val_loader = _Stage2LoaderAdapter(val_loader_raw)

    train_dataset = getattr(train_loader_raw, "dataset", None)
    stage_vocab = getattr(train_dataset, "vocab", None)
    if stage_vocab is None:
        raise ValueError("Stage 2 dataset is missing vocabulary")

    vocab_size = len(stage_vocab.itos)
    pad_token_id = stage_vocab.pad_id
    kmer = stage_vocab.k

    embedding_dim = int(cfg_run.get("embedding_dim", 192))
    transformer_layers = int(cfg_run.get("transformer_layers", 6))
    transformer_heads = int(cfg_run.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg_run.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg_run.get("transformer_dropout", 0.15))

    use_attention_pool = bool(cfg_run.get("use_attention_pool", True))
    use_tcn = bool(cfg_run.get("use_tcn", True))
    tcn_hidden = int(cfg_run.get("tcn_hidden", 256))
    tcn_levels = int(cfg_run.get("tcn_levels", 4))
    tcn_kernel = int(cfg_run.get("tcn_kernel", 3))
    multiscale_channels = int(cfg_run.get("multiscale_channels", 64))
    multiscale_kernels_cfg = cfg_run.get("multiscale_kernels", [3, 5, 7, 9, 15])
    multiscale_kernels = tuple(multiscale_kernels_cfg)
    lstm_hidden = int(cfg_run.get("lstm_hidden", 192))

    post_cnn_transformer_layers = int(cfg_run.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg_run.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg_run.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg_run.get("fusion_hidden", 256))

    stage2_use_engineered = bool(cfg_run.get("stage2_use_engineered_features", True))
    engineered_dim = int(cfg_run.get("engineered_dim", 288))

    # Stage 2 head configuration
    stage2_head_type = str(cfg_run.get("stage2_head_type", "transformer"))
    if stage2_head_type not in {"transformer", "bilstm", "both"}:
        raise ValueError(f"Invalid stage2_head_type: {stage2_head_type}. Must be 'transformer', 'bilstm', or 'both'")
    stage2_transformer_layers = int(cfg_run.get("stage2_transformer_layers", 2))
    stage2_transformer_heads = int(cfg_run.get("stage2_transformer_heads", 4))
    stage2_transformer_ff_dim = cfg_run.get("stage2_transformer_ff_dim")  # None uses default (2*dim)
    if stage2_transformer_ff_dim is not None:
        stage2_transformer_ff_dim = int(stage2_transformer_ff_dim)
    stage2_bilstm_hidden = int(cfg_run.get("stage2_bilstm_hidden", 128))
    stage2_bilstm_layers = int(cfg_run.get("stage2_bilstm_layers", 1))
    stage2_combine_mode = str(cfg_run.get("stage2_combine_mode", "concat"))
    stage2_combine_order = str(cfg_run.get("stage2_combine_order", "transformer_first"))

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=pad_token_id,
        engineered_dim=engineered_dim,
        use_engineered_features=stage2_use_engineered,
        use_attention_pool=use_attention_pool,
        use_tcn=use_tcn,
        tcn_hidden=tcn_hidden,
        tcn_levels=tcn_levels,
        tcn_kernel=tcn_kernel,
        multiscale_channels=multiscale_channels,
        multiscale_kernels=multiscale_kernels,
        lstm_hidden=lstm_hidden,
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        # Stage 2 head configuration
        stage2_head_type=stage2_head_type,
        stage2_transformer_layers=stage2_transformer_layers,
        stage2_transformer_heads=stage2_transformer_heads,
        stage2_transformer_ff_dim=stage2_transformer_ff_dim,
        stage2_bilstm_hidden=stage2_bilstm_hidden,
        stage2_bilstm_layers=stage2_bilstm_layers,
        stage2_combine_mode=stage2_combine_mode,
        stage2_combine_order=stage2_combine_order,
        # Engineered features MLP
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        fusion_hidden=fusion_hidden,
        # Stochastic depth (drop path)
        drop_path_rate=stage2_drop_path_rate,
    ).to(device)

    LOGGER.info("=" * 60)
    LOGGER.info("MODEL ARCHITECTURE (Stage 2)")
    LOGGER.info("=" * 60)
    LOGGER.info("Model: GeneWhispererStage2 (dedicated strength head)")
    LOGGER.info("Task: strong vs weak promoters (binary)")
    LOGGER.info("Architecture: Embedding → CNN/TCN → Transformer → StrengthHead → Fusion → Classifier")
    if use_tcn:
        LOGGER.info(
            "CNN/TCN: Multi-scale kernels %s → TCN(%d levels, %d hidden, k=%d)",
            multiscale_kernels,
            tcn_levels,
            tcn_hidden,
            tcn_kernel,
        )
    LOGGER.info("Post-CNN Transformer: %d layers, %d heads", post_cnn_transformer_layers, transformer_heads)
    LOGGER.info("Strength Head Type: %s", stage2_head_type)
    if stage2_head_type in {"transformer", "both"}:
        LOGGER.info(
            "  Transformer Head: %d layers, %d heads, ff_dim=%s",
            stage2_transformer_layers,
            stage2_transformer_heads,
            stage2_transformer_ff_dim or "auto",
        )
    if stage2_head_type in {"bilstm", "both"}:
        LOGGER.info(
            "  BiLSTM Head: hidden=%d, layers=%d",
            stage2_bilstm_hidden,
            stage2_bilstm_layers,
        )
    if stage2_head_type == "both":
        LOGGER.info(
            "  Combine Mode: %s (%s)",
            stage2_combine_mode,
            stage2_combine_order,
        )
    if stage2_use_engineered:
        LOGGER.info(
            "Engineered Features: %d dims (TNC+PseEIIP) via MLP %d→%d",
            engineered_dim,
            engineered_mlp_hidden,
            engineered_mlp_output,
        )
    LOGGER.info("=" * 60)

    # Log effective Stage 2 hyperparameters
    use_focal_loss_preview = bool(cfg_run.get("stage2_use_focal_loss", cfg_run.get("stage1_use_focal_loss", True)))
    LOGGER.info("=" * 60)
    LOGGER.info("EFFECTIVE STAGE 2 HYPERPARAMETERS")
    LOGGER.info("=" * 60)
    LOGGER.info(
        "  label_smoothing=%.3f | use_mixup=%s (alpha=%.2f)",
        stage2_label_smoothing,
        stage2_use_mixup,
        stage2_mixup_alpha,
    )
    LOGGER.info(
        "  use_swa=%s (start_epoch=%d, lr=%.2e)",
        stage2_use_swa,
        stage2_swa_start_epoch,
        stage2_swa_lr,
    )
    LOGGER.info(
        "  drop_path_rate=%.3f | focal_loss=%s",
        stage2_drop_path_rate,
        use_focal_loss_preview,
    )
    LOGGER.info("=" * 60)

    load_mlm_weights = bool(cfg_run.get("stage2_load_mlm_weights", False))
    transfer_mode = str(cfg_run.get("stage2_mlm_transfer_mode", "embed_only"))
    if transfer_mode not in {"none", "embed_only", "embed_plus_adapter"}:
        raise ValueError(f"Invalid stage2_mlm_transfer_mode: {transfer_mode}")
    mlm_loaded = False
    mlm_checkpoint = cfg_run.get("mlm_encoder_checkpoint") or cfg_run.get("mlm_encoder_path")
    # Stage 1 checkpoint takes priority over MLM weights - skip MLM loading if stage1_ckpt is provided
    should_load_mlm = load_mlm_weights and transfer_mode != "none" and not stage1_ckpt_path
    with torch.no_grad():
        emb_before = model.embedding.weight.detach().float().clone()
        emb_before_norm = emb_before.norm().item()
    load_attempted = False
    if should_load_mlm and mlm_checkpoint:
        checkpoint_path = _resolve_optional_path(mlm_checkpoint, base_dir=script_dir)
        if checkpoint_path and checkpoint_path.exists():
            LOGGER.info(
                "Loading MLM encoder weights from %s (transfer_mode=%s)",
                checkpoint_path,
                transfer_mode,
            )
            load_attempted = True
            model.load_pretrained_weights(
                checkpoint_path, strict=False, transfer_mode=transfer_mode
            )
            if transfer_mode == "embed_plus_adapter":
                LOGGER.info(
                    "Pretrained MLM weights loaded (embedding + %d transformer layers)",
                    post_cnn_transformer_layers,
                )
            else:
                LOGGER.info("Pretrained MLM weights loaded (embedding only)")
            mlm_loaded = True
            freeze_layers = int(cfg_run.get("stage2_freeze_layers", cfg_run.get("stage1_freeze_layers", 0)) or 0)
            if freeze_layers > 0:
                LOGGER.info("Freezing bottom %d encoder layers", freeze_layers)
                model.encoder.freeze_bottom_layers(freeze_layers)
        else:
            LOGGER.warning("MLM checkpoint not found; training without MLM init (path=%s)", checkpoint_path)
    elif load_mlm_weights and transfer_mode == "none":
        LOGGER.info("stage2_mlm_transfer_mode=none; skipping MLM weight load")
    elif load_mlm_weights and stage1_ckpt_path:
        LOGGER.info("Skipping MLM weight load (stage1_ckpt takes priority)")
    with torch.no_grad():
        emb_after = model.embedding.weight.detach().float()
        emb_after_norm = emb_after.norm().item()
        diff_norm = (emb_after - emb_before).norm().item()
    LOGGER.info(
        "Load sanity: embedding norm before=%.6f after=%.6f diff_norm=%.6f (load_attempted=%s)",
        emb_before_norm,
        emb_after_norm,
        diff_norm,
        load_attempted,
    )
    if load_attempted and diff_norm < 1e-6:
        LOGGER.warning("Pretrain load had no effect (likely key mismatch or already loaded).")

    stage1_ckpt_loaded = False
    if stage1_ckpt_path and stage1_ckpt_path.exists():
        LOGGER.info("Initializing from Stage 1 checkpoint: %s", stage1_ckpt_path)
        # Use GeneWhispererStage2's dedicated method for loading Stage 1 weights
        loaded = model.load_stage1_weights(stage1_ckpt_path, device=device)
        stage1_ckpt_loaded = loaded > 0
        if not stage1_ckpt_loaded:
            suggestion = (script_dir / "../artifacts/checkpoints" / f"stage1_k{kmer}.pt").resolve()
            if suggestion.exists():
                LOGGER.warning("Stage 1 checkpoint appears incompatible; try %s", suggestion)
    elif stage1_ckpt_path:
        LOGGER.warning("Stage 1 checkpoint not found at %s; training without Stage 1 init", stage1_ckpt_path)

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    LOGGER.info("Model: %.2fM total params, %.2fM trainable", total_params / 1e6, trainable_params / 1e6)

    # Update resolved_config with device and model info, then re-save
    resolved_config["system_info"]["device"] = str(device)
    resolved_config["system_info"]["stage2_total_params"] = total_params
    resolved_config["system_info"]["stage2_trainable_params"] = trainable_params
    with resolved_config_path.open("w", encoding="utf-8") as handle:
        json.dump(resolved_config, handle, indent=2)

    lr = float(cfg_run.get("lr", 2e-4))
    weight_decay = float(cfg_run.get("weight_decay", 0.02))
    layer_lr_decay = float(cfg_run.get("layer_lr_decay", 0.8))
    optimizer_foreach = bool(cfg_run.get("optimizer_foreach", True))

    LOGGER.info("Optimizer: AdamW foreach=%s", optimizer_foreach)

    use_llrd = bool(mlm_loaded or stage1_ckpt_loaded)
    if use_llrd:
        LOGGER.info("Using layer-wise LR decay (factor=%.2f) for pretrained model", layer_lr_decay)
        param_groups = stage1.get_layer_wise_lr_groups(
            model,
            base_lr=lr,
            lr_decay=layer_lr_decay,
            weight_decay=weight_decay,
        )
        optimizer = torch.optim.AdamW(param_groups, betas=(0.9, 0.999), foreach=optimizer_foreach)
        LOGGER.info("Created %d parameter groups with LLRD", len(param_groups))
    else:
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay,
            betas=(0.9, 0.999),
            foreach=optimizer_foreach,
        )

    warmup_ratio = float(cfg_run.get("warmup_ratio", 0.15))
    steps_per_epoch = len(train_loader)
    total_steps = epochs * steps_per_epoch
    warmup_steps = int(warmup_ratio * total_steps)

    scheduler = stage1.get_cosine_schedule_with_warmup(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=0.02,
    )
    LOGGER.info("Using cosine schedule with %d warmup steps", warmup_steps)

    # Use effective Stage 2 SWA parameters (computed earlier with fallbacks)
    swa_model = None
    swa_scheduler = None

    if stage2_use_swa:
        swa_model = AveragedModel(model)
        swa_scheduler = SWALR(optimizer, swa_lr=stage2_swa_lr, anneal_epochs=10)
        LOGGER.info("SWA enabled: starts at epoch %d with LR %.2e", stage2_swa_start_epoch, stage2_swa_lr)

    # Use effective Stage 2 label smoothing (computed earlier with fallbacks)
    use_focal_loss = bool(cfg_run.get("stage2_use_focal_loss", cfg_run.get("stage1_use_focal_loss", True)))
    focal_alpha = float(cfg_run.get("stage2_focal_alpha", cfg_run.get("stage1_focal_alpha", 0.5)))
    focal_gamma = float(cfg_run.get("stage2_focal_gamma", cfg_run.get("stage1_focal_gamma", 2.0)))

    if use_focal_loss:
        criterion = stage1.FocalLossWithSmoothing(
            alpha=focal_alpha,
            gamma=focal_gamma,
            smoothing=stage2_label_smoothing,
        )
        LOGGER.info("Using focal loss with smoothing=%.2f", stage2_label_smoothing)
    else:
        criterion = stage1.LabelSmoothingBCE(smoothing=stage2_label_smoothing)
        LOGGER.info("Using BCE with smoothing=%.2f", stage2_label_smoothing)

    # =========================================================================
    # OVERFIT DEBUG MODE: Train on 32 samples for 200 steps
    # =========================================================================
    if bool(cfg_run.get("overfit_debug", False)):
        LOGGER.info("Running OVERFIT DEBUG mode for Stage 2")
        # Use simple BCEWithLogitsLoss for overfit test (no label smoothing or focal loss)
        simple_criterion = nn.BCEWithLogitsLoss()
        stage1.run_overfit_debug_stage1(
            model=model,
            train_loader=train_loader,
            criterion=simple_criterion,
            device=device,
            num_steps=200,
            lr=1e-3,
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return {
            "overfit_debug": True,
            "stage2_acc": float("nan"),
            "stage2_auc": float("nan"),
            "stage2_mcc": float("nan"),
            "best_ckpt_path": None,
            "stage2_seconds": float(time.perf_counter() - start_time),
            "params": int(total_params),
            "best_val_acc": float("nan"),
            "best_val_auc": float("nan"),
            "best_val_mcc": float("nan"),
            "best_epoch": 0,
        }

    # =========================================================================
    # PROFILING MODE: Run profiler for ~50 steps
    # =========================================================================
    if bool(cfg_run.get("profile", False)):
        LOGGER.info("Running PROFILING mode for Stage 2")
        grad_accum_steps_prof = int(cfg_run.get("grad_accum_steps", 4))
        grad_clip_norm_prof = float(cfg_run.get("grad_clip_norm", 1.0))

        profile_metrics = stage1.run_profiling_epoch(
            loader=train_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            run_dir=run_dir,
            num_profile_steps=50,
            grad_accum_steps=grad_accum_steps_prof,
            grad_clip_norm=grad_clip_norm_prof,
        )

        # Save profile metrics to run_dir
        import json as json_module
        profile_report_path = run_dir / "profile_report.json"
        with profile_report_path.open("w", encoding="utf-8") as f:
            json_module.dump(profile_metrics, f, indent=2)
        LOGGER.info("Saved profile report to %s", profile_report_path)
        LOGGER.info("Profiling complete. Exiting.")

        return {
            "profile": True,
            "profile_metrics": profile_metrics,
            "stage2_acc": float("nan"),
            "stage2_auc": float("nan"),
            "stage2_mcc": float("nan"),
            "best_ckpt_path": None,
            "stage2_seconds": float(time.perf_counter() - start_time),
            "params": int(total_params),
            "best_val_acc": float("nan"),
            "best_val_auc": float("nan"),
            "best_val_mcc": float("nan"),
            "best_epoch": 0,
        }

    grad_accum_steps = int(cfg_run.get("grad_accum_steps", 2))
    grad_clip_norm = float(cfg_run.get("grad_clip_norm", 1.0))
    param_norm_warn = float(cfg_run.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg_run.get("param_norm_cap", 200.0))
    # Use effective Stage 2 mixup parameters (computed earlier with fallbacks)
    patience = int(cfg_run.get("early_stopping_patience", 15))

    global_optimizer_steps = 0

    class_balance = stage1.compute_class_balance(train_loader_raw)
    LOGGER.info(
        "Class balance: %.2f%% positive (%d/%d)",
        class_balance["positive_ratio"] * 100.0,
        class_balance["positive"],
        class_balance["total"],
    )

    artifacts_root = (script_dir / "../artifacts").resolve()
    checkpoint_name = cfg_run.get("stage2_checkpoint_name") or f"stage2_k{kmer}.pt"
    checkpoint_path = artifacts_root / "checkpoints" / checkpoint_name
    report_path = artifacts_root / f"stage2_report_k{kmer}.json"

    best_val_mcc = -float("inf")
    best_val_acc = 0.0
    best_epoch = 0
    best_train_metrics = None
    best_val_metrics = None
    patience_counter = 0
    swa_started = False

    for epoch in range(1, epochs + 1):
        in_swa_phase = stage2_use_swa and epoch >= stage2_swa_start_epoch
        if in_swa_phase and not swa_started:
            LOGGER.info("=" * 40)
            LOGGER.info("Starting SWA phase at epoch %d", epoch)
            LOGGER.info("=" * 40)
            swa_started = True

        current_lr = optimizer.param_groups[0]["lr"]
        phase_str = "[SWA]" if in_swa_phase else ""
        LOGGER.info("Epoch %d/%d %s (LR: %.2e)", epoch, epochs, phase_str, current_lr)

        remaining_steps = None if ablation_max_steps is None else max(ablation_max_steps - global_optimizer_steps, 0)
        train_metrics = stage1.run_epoch(
            train_loader,
            model,
            criterion,
            device,
            optimizer=optimizer,
            scheduler=None if in_swa_phase else scheduler,
            desc="Train",
            grad_accum_steps=grad_accum_steps,
            max_optimizer_steps=remaining_steps,
            grad_clip_norm=grad_clip_norm,
            param_norm_warn=param_norm_warn,
            param_norm_cap=param_norm_cap,
            use_mixup=stage2_use_mixup and not in_swa_phase,
            mixup_alpha=stage2_mixup_alpha,
            amp_ctx=amp_ctx,
            scaler=scaler,
        )
        global_optimizer_steps += int(train_metrics.get("optimizer_steps", 0))
        reached_max_steps = ablation_max_steps is not None and global_optimizer_steps >= ablation_max_steps

        if in_swa_phase and swa_model is not None and swa_scheduler is not None:
            swa_model.update_parameters(model)
            swa_scheduler.step()

        val_metrics = stage1.run_epoch(
            val_loader,
            model,
            criterion,
            device,
            optimizer=None,
            desc="Val",
            amp_ctx=amp_ctx,
        )

        swa_val_metrics = None
        if in_swa_phase and swa_model is not None and epoch % 5 == 0:
            try:
                swa_model.train()
                with torch.no_grad():
                    for i, (tokens, eng, _) in enumerate(train_loader):
                        if i >= 50:
                            break
                        tokens = tokens.to(device)
                        eng = eng.to(device)
                        _ = swa_model(tokens, eng)

                swa_model.eval()
                swa_val_metrics = stage1.run_epoch(
                    val_loader,
                    swa_model,
                    criterion,
                    device,
                    optimizer=None,
                    desc="Val-SWA",
                    amp_ctx=amp_ctx,
                )
            except Exception as exc:
                LOGGER.warning("SWA evaluation failed: %s", exc)
                swa_val_metrics = None

        LOGGER.info(
            "Train - Loss: %.4f | Acc: %.4f | MCC: %.4f | AUC: %.4f",
            train_metrics["loss"],
            train_metrics["accuracy"],
            train_metrics["mcc"],
            train_metrics["roc_auc"],
        )
        LOGGER.info(
            "Val   - Loss: %.4f | Acc: %.4f | MCC: %.4f | AUC: %.4f",
            val_metrics["loss"],
            val_metrics["accuracy"],
            val_metrics["mcc"],
            val_metrics["roc_auc"],
        )

        if swa_val_metrics:
            LOGGER.info(
                "SWA   - Loss: %.4f | Acc: %.4f | MCC: %.4f | AUC: %.4f",
                swa_val_metrics["loss"],
                swa_val_metrics["accuracy"],
                swa_val_metrics["mcc"],
                swa_val_metrics["roc_auc"],
            )

        eval_metrics = val_metrics
        if swa_val_metrics and swa_val_metrics["mcc"] > val_metrics["mcc"]:
            eval_metrics = swa_val_metrics
            LOGGER.info("(Using SWA metrics for checkpointing)")

        improved = False
        if eval_metrics["mcc"] > best_val_mcc + 1e-4:
            improved = True
        elif abs(eval_metrics["mcc"] - best_val_mcc) < 1e-4 and eval_metrics["accuracy"] > best_val_acc:
            improved = True

        if improved:
            patience_counter = 0
            best_val_mcc = eval_metrics["mcc"]
            best_val_acc = eval_metrics["accuracy"]
            best_epoch = epoch
            best_train_metrics = train_metrics
            best_val_metrics = eval_metrics

            if swa_val_metrics and swa_val_metrics["mcc"] > val_metrics["mcc"] and swa_model is not None:
                stage1.save_checkpoint(checkpoint_path, swa_model.module, optimizer, scheduler, epoch, eval_metrics)
            else:
                stage1.save_checkpoint(checkpoint_path, model, optimizer, scheduler, epoch, val_metrics)
            LOGGER.info(
                "★ New best! MCC: %.4f, Acc: %.4f (%.2f%%)",
                best_val_mcc,
                best_val_acc,
                best_val_acc * 100,
            )
        else:
            patience_counter += 1
            effective_patience = patience + 10 if in_swa_phase else patience
            LOGGER.info("No improvement. Patience %d/%d", patience_counter, effective_patience)
            if patience_counter >= effective_patience:
                LOGGER.info("Early stopping triggered at epoch %d", epoch)
                break

        if reached_max_steps:
            LOGGER.info(
                "Reached ablation_max_steps_stage2=%d (optimizer steps); stopping after validation",
                ablation_max_steps,
            )
            break

    if best_val_metrics is None:
        best_train_metrics = train_metrics
        best_val_metrics = val_metrics
        best_epoch = epochs

    report = {
        "config": config_id,
        "best_epoch": best_epoch,
        "train_metrics": best_train_metrics,
        "val_metrics": best_val_metrics,
        "class_balance": class_balance,
        "device": str(device),
        "seed": seed,
        "model_params": {
            "total": total_params,
            "trainable": trainable_params,
            "embedding_dim": embedding_dim,
            "layers": transformer_layers,
            "heads": transformer_heads,
            "kmer": kmer,
            "engineered_dim": engineered_dim,
            "stage2_use_engineered_features": stage2_use_engineered,
            "stage2_head_type": stage2_head_type,
            "stage2_transformer_layers": stage2_transformer_layers if stage2_head_type in {"transformer", "both"} else None,
            "stage2_transformer_heads": stage2_transformer_heads if stage2_head_type in {"transformer", "both"} else None,
            "stage2_bilstm_hidden": stage2_bilstm_hidden if stage2_head_type in {"bilstm", "both"} else None,
            "stage2_combine_mode": stage2_combine_mode if stage2_head_type == "both" else None,
        },
        "training_enhancements": {
            "layer_wise_lr_decay": layer_lr_decay,
            "label_smoothing": stage2_label_smoothing,
            "mixup_alpha": stage2_mixup_alpha,
            "use_mixup": stage2_use_mixup,
            "use_swa": stage2_use_swa,
            "swa_start_epoch": stage2_swa_start_epoch if stage2_use_swa else None,
            "swa_lr": stage2_swa_lr if stage2_use_swa else None,
            "drop_path_rate": stage2_drop_path_rate,
            "ablation_max_steps_stage2": ablation_max_steps,
            "ablation_disable_wandb": bool(cfg_run.get("ablation_disable_wandb", False)),
            "stage1_ckpt": str(stage1_ckpt_path) if stage1_ckpt_path else None,
        },
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    LOGGER.info("Saved training report to %s", report_path)

    stage2_seconds = float(time.perf_counter() - start_time)
    best_val_auc = float(best_val_metrics.get("roc_auc", 0.5)) if best_val_metrics else 0.5
    return {
        "stage2_acc": float(best_val_acc),
        "stage2_auc": best_val_auc,
        "stage2_mcc": float(best_val_mcc),
        "stage2_seconds": stage2_seconds,
        "params": int(total_params),
        "best_val_acc": float(best_val_acc),
        "best_val_auc": best_val_auc,
        "best_val_mcc": float(best_val_mcc),
        "best_epoch": int(best_epoch),
        "best_ckpt_path": str(checkpoint_path),
    }


if __name__ == "__main__":
    main()
