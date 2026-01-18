#!/usr/bin/env python3
"""
Final Training Script for Best Model Configuration.

Trains the model with the best hyperparameters from tuning for full 100 epochs.
Includes comprehensive evaluation and comparison with the original model.

Features:
    - Learning rate warmup (10% of steps)
    - Cosine annealing schedule
    - Early stopping (patience=15)
    - Model checkpointing (save best and last)
    - Detailed logging
    - Final evaluation with all metrics
    - Visualization outputs

Usage:
    python train_final_model.py --config config.yaml
    python train_final_model.py --best-config results/hyperparameter_tuning/best_config.yaml
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import torch
import torch.nn as nn
import yaml
from torch.utils.data import DataLoader

LOGGER = logging.getLogger("gene_whisperer.final_training")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Import MLM weight loading function
from train_stage1 import get_mlm_checkpoint_for_kmer


# Best hyperparameters from tuning (from screenshot)
BEST_CONFIG = {
    "variant": "combined",
    "stage1_lr": 0.0001,
    "transformer_dropout": 0.1,
    "classifier_hidden": 256,
    "pooling_type": "attention",
    "weight_decay": 0.02,
}

# Original model baseline for comparison
ORIGINAL_MODEL_BASELINE = {
    "val_acc": 0.8293,  # 82.93%
    "val_mcc": 0.6587,
    "val_auc": 0.8954,
    "param_count": 39_934_436,  # ~38.4M
    "training_time_minutes": 21.1,
}


@dataclass
class TrainingConfig:
    """Configuration for final training."""

    variant: str = "combined"
    lr: float = 0.0001
    dropout: float = 0.1
    classifier_hidden: int = 256
    pooling_type: str = "attention"
    weight_decay: float = 0.02
    epochs: int = 100
    patience: int = 15
    warmup_ratio: float = 0.10
    grad_accum_steps: int = 4
    grad_clip_norm: float = 1.0
    seed: int = 42

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "TrainingConfig":
        """Create config from dictionary."""
        return cls(
            variant=d.get("variant", "combined"),
            lr=d.get("stage1_lr", d.get("lr", 0.0001)),
            dropout=d.get("transformer_dropout", d.get("dropout", 0.1)),
            classifier_hidden=d.get("classifier_hidden", 256),
            pooling_type=d.get("pooling_type", "attention"),
            weight_decay=d.get("weight_decay", 0.02),
            epochs=d.get("epochs", 100),
            patience=d.get("patience", 15),
            warmup_ratio=d.get("warmup_ratio", 0.10),
            grad_accum_steps=d.get("grad_accum_steps", 4),
            grad_clip_norm=d.get("grad_clip_norm", 1.0),
            seed=d.get("seed", 42),
        )


@dataclass
class TrainingResult:
    """Results from final training."""

    best_val_acc: float = 0.0
    best_val_mcc: float = 0.0
    best_val_auc: float = 0.0
    best_val_f1: float = 0.0
    best_epoch: int = 0
    final_train_loss: float = 0.0
    param_count: int = 0
    training_time_seconds: float = 0.0
    stopped_early: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class CosineAnnealingWithWarmup:
    """Learning rate scheduler with linear warmup and cosine annealing."""

    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        warmup_steps: int,
        total_steps: int,
        min_lr_ratio: float = 0.1,
    ):
        self.optimizer = optimizer
        self.warmup_steps = warmup_steps
        self.total_steps = total_steps
        self.min_lr_ratio = min_lr_ratio
        self.base_lrs = [group["lr"] for group in optimizer.param_groups]
        self.current_step = 0

    def step(self) -> None:
        """Update learning rate based on current step."""
        self.current_step += 1
        lr_multiplier = self._get_lr_multiplier()

        for param_group, base_lr in zip(self.optimizer.param_groups, self.base_lrs):
            param_group["lr"] = base_lr * lr_multiplier

    def _get_lr_multiplier(self) -> float:
        """Calculate LR multiplier for current step."""
        if self.current_step < self.warmup_steps:
            # Linear warmup
            return self.current_step / max(1, self.warmup_steps)
        else:
            # Cosine annealing
            progress = (self.current_step - self.warmup_steps) / max(
                1, self.total_steps - self.warmup_steps
            )
            return max(
                self.min_lr_ratio,
                0.5 * (1 + math.cos(math.pi * progress)),
            )

    def get_last_lr(self) -> List[float]:
        """Get current learning rates."""
        return [group["lr"] for group in self.optimizer.param_groups]


def load_base_config(config_path: Path) -> dict:
    """Load base configuration from YAML file."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def create_tunable_model_wrapper(
    variant: str,
    cfg: dict,
    training_config: TrainingConfig,
    device: torch.device,
) -> nn.Module:
    """Create tunable model with hyperparameters (for experimentation)."""
    from tune_best_variant import create_tunable_model

    hp_config = {
        "lr": training_config.lr,
        "dropout": training_config.dropout,
        "classifier_hidden": training_config.classifier_hidden,
        "pooling_type": training_config.pooling_type,
        "weight_decay": training_config.weight_decay,
    }

    return create_tunable_model(variant, cfg, hp_config, device)


def create_standard_model(
    cfg: dict,
    training_config: TrainingConfig,
    device: torch.device,
    vocab: "KmerVocabulary",
) -> nn.Module:
    """
    Create legacy GeneWhispererStage1 model compatible with ensemble scripts.

    This model can be used with ensemble_infer.py and eval_stage1_ensemble.py.
    Applies tuned hyperparameters (dropout) to the standard architecture.
    """
    from model import GeneWhispererStage1Legacy

    kmer = int(cfg.get("kmer", 6))
    max_bp_len = int(cfg.get("max_bp_len", 81))

    model = GeneWhispererStage1Legacy(
        vocab_size=len(vocab.itos),
        kmer=kmer,
        embedding_dim=int(cfg.get("embedding_dim", 384)),
        num_layers=int(cfg.get("transformer_layers", 12)),
        num_heads=int(cfg.get("transformer_heads", 12)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1536)),
        dropout=training_config.dropout,  # Use tuned dropout
        pad_token_id=vocab.pad_id,
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 384)),
        tcn_levels=int(cfg.get("tcn_levels", 5)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 512)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 384)),
        fusion_hidden=int(cfg.get("fusion_hidden", 512)),
        max_seq_len=max_bp_len - kmer + 1,
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", True)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", True)),
    ).to(device)

    return model


def create_model(
    variant: str,
    cfg: dict,
    training_config: TrainingConfig,
    device: torch.device,
    use_standard_model: bool = False,
    vocab: Optional["KmerVocabulary"] = None,
) -> nn.Module:
    """
    Create model with tuned hyperparameters.

    Args:
        variant: Model variant (combined, transformer_only, features_only)
        cfg: Base configuration dictionary
        training_config: Training configuration with tuned hyperparameters
        device: Target device
        use_standard_model: If True, use GeneWhispererStage1Legacy (compatible with ensemble scripts)
        vocab: Required if use_standard_model=True

    Returns:
        Created model on the specified device
    """
    if use_standard_model:
        if vocab is None:
            raise ValueError("vocab is required when use_standard_model=True")
        LOGGER.info("Using legacy GeneWhispererStage1 model (ensemble-compatible)")
        return create_standard_model(cfg, training_config, device, vocab)
    else:
        LOGGER.info("Using tunable model variant: %s", variant)
        return create_tunable_model_wrapper(variant, cfg, training_config, device)


def count_parameters(model: nn.Module) -> int:
    """Count trainable parameters."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def load_pretrained_mlm_weights(
    model: nn.Module,
    cfg: dict,
    transfer_mode: str = "embed_only",
) -> bool:
    """
    Load pretrained MLM weights into the model.

    Args:
        model: Model with load_pretrained_weights method
        cfg: Configuration dictionary
        transfer_mode: One of "embed_only", "embed_plus_adapter", or "none"

    Returns:
        True if weights were successfully loaded, False otherwise
    """
    if transfer_mode == "none":
        LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
        return False

    if not hasattr(model, "load_pretrained_weights"):
        LOGGER.warning("Model does not support pretrained weight loading")
        return False

    kmer = int(cfg.get("kmer", 6))
    checkpoint_path = get_mlm_checkpoint_for_kmer(cfg, kmer)

    if checkpoint_path is None:
        LOGGER.warning("No MLM checkpoint found for k=%d - training from scratch", kmer)
        return False

    LOGGER.info(
        "Loading MLM encoder weights for k=%d from %s (transfer_mode=%s)",
        kmer,
        checkpoint_path,
        transfer_mode,
    )

    # Capture embedding norm before loading for sanity check
    embedding_layer = getattr(model, "embedding", None)
    has_embedding = isinstance(embedding_layer, nn.Embedding)
    if has_embedding:
        with torch.no_grad():
            emb_before = embedding_layer.weight.detach().float().clone()
            emb_before_norm = emb_before.norm().item()
    else:
        emb_before = None
        emb_before_norm = 0.0

    # Load the weights
    model.load_pretrained_weights(
        checkpoint_path=checkpoint_path,
        strict=False,
        transfer_mode=transfer_mode,
    )
    LOGGER.info("Pretrained MLM weights loaded successfully")

    # Verify weights actually changed
    if has_embedding:
        with torch.no_grad():
            emb_after = embedding_layer.weight.detach().float()
            emb_after_norm = emb_after.norm().item()
            diff_norm = (emb_after - emb_before).norm().item()
        LOGGER.info(
            "Load sanity: embedding norm before=%.6f after=%.6f diff_norm=%.6f",
            emb_before_norm,
            emb_after_norm,
            diff_norm,
        )
        if diff_norm < 1e-6:
            LOGGER.warning(
                "Pretrain load had no effect (likely key mismatch or already loaded)"
            )
            return False

    return True


def train_final_model(
    cfg: dict,
    training_config: TrainingConfig,
    train_loader: DataLoader,
    val_loader: DataLoader,
    device: torch.device,
    output_dir: Path,
    use_standard_model: bool = False,
    vocab: Optional["KmerVocabulary"] = None,
    load_mlm_weights: bool = True,
    mlm_transfer_mode: str = "embed_only",
) -> Tuple[nn.Module, TrainingResult, Dict[str, List[float]]]:
    """
    Train the final model with best configuration.

    Args:
        cfg: Base configuration dictionary
        training_config: Training configuration with tuned hyperparameters
        train_loader: Training data loader
        val_loader: Validation data loader
        device: Target device
        output_dir: Directory to save checkpoints and results
        use_standard_model: If True, use GeneWhispererStage1Legacy (compatible with ensemble scripts)
        vocab: Vocabulary (required if use_standard_model=True)
        load_mlm_weights: Whether to load pretrained MLM weights (default: True)
        mlm_transfer_mode: Transfer mode for MLM weights ("embed_only", "embed_plus_adapter", "none")

    Returns:
        Tuple of (trained model, training results, training history)
    """
    from amp_utils import get_amp_context, create_grad_scaler
    from early_stopping import EarlyStopping
    from train_stage1 import run_epoch

    result = TrainingResult()
    history = {
        "train_loss": [],
        "train_acc": [],
        "val_loss": [],
        "val_acc": [],
        "val_mcc": [],
        "val_auc": [],
        "val_f1": [],
        "learning_rate": [],
    }

    # Create model
    model = create_model(
        training_config.variant, cfg, training_config, device,
        use_standard_model=use_standard_model, vocab=vocab
    )

    # Determine model variant for run_epoch (standard model is always "combined"-like)
    model_variant = "combined" if use_standard_model else training_config.variant
    result.param_count = count_parameters(model)
    LOGGER.info("Model created with %d parameters (%.2fM)",
                result.param_count, result.param_count / 1e6)

    # Load pretrained MLM weights if enabled
    if load_mlm_weights:
        mlm_loaded = load_pretrained_mlm_weights(
            model=model,
            cfg=cfg,
            transfer_mode=mlm_transfer_mode,
        )
        if not mlm_loaded:
            LOGGER.warning("MLM weight loading was requested but failed or had no effect")

    # Optimizer
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=training_config.lr,
        weight_decay=training_config.weight_decay,
    )

    # Scheduler with warmup + cosine annealing
    total_steps = training_config.epochs * len(train_loader)
    warmup_steps = int(total_steps * training_config.warmup_ratio)
    scheduler = CosineAnnealingWithWarmup(
        optimizer,
        warmup_steps=warmup_steps,
        total_steps=total_steps,
    )

    LOGGER.info("Training for %d epochs (%d total steps, %d warmup steps)",
                training_config.epochs, total_steps, warmup_steps)

    # Loss and early stopping
    criterion = nn.BCEWithLogitsLoss()
    early_stopping = EarlyStopping(
        patience=training_config.patience,
        mode="max",
        min_delta=0.001,
        min_epochs=10,
        restore_best=True,
        verbose=True,
    )

    # AMP setup
    amp_ctx = get_amp_context(device, cfg)
    amp_enabled = bool(cfg.get("amp_enabled", True))
    scaler = create_grad_scaler(device) if amp_enabled else None

    # Checkpoint directory
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    best_model_state = None
    start_time = time.time()

    for epoch in range(training_config.epochs):
        epoch_start = time.time()

        # Training epoch
        model.train()
        train_metrics = run_epoch(
            loader=train_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            scheduler=scheduler,
            desc=f"Epoch {epoch+1}/{training_config.epochs} [Train]",
            grad_accum_steps=training_config.grad_accum_steps,
            grad_clip_norm=training_config.grad_clip_norm,
            use_mixup=False,
            amp_ctx=amp_ctx,
            scaler=scaler,
            model_variant=model_variant,
        )

        # Validation epoch
        model.eval()
        with torch.no_grad():
            val_metrics = run_epoch(
                loader=val_loader,
                model=model,
                criterion=criterion,
                device=device,
                optimizer=None,
                scheduler=None,
                desc=f"Epoch {epoch+1}/{training_config.epochs} [Val]",
                grad_accum_steps=1,
                grad_clip_norm=0.0,
                use_mixup=False,
                amp_ctx=amp_ctx,
                scaler=None,
                model_variant=model_variant,
            )

        # Record history
        history["train_loss"].append(train_metrics.get("loss", 0.0))
        history["train_acc"].append(train_metrics.get("accuracy", 0.0))
        history["val_loss"].append(val_metrics.get("loss", 0.0))
        history["val_acc"].append(val_metrics.get("accuracy", 0.0))
        history["val_mcc"].append(val_metrics.get("mcc", 0.0))
        history["val_auc"].append(val_metrics.get("roc_auc", 0.0))
        history["val_f1"].append(val_metrics.get("f1", 0.0))
        history["learning_rate"].append(scheduler.get_last_lr()[0])

        val_acc = val_metrics.get("accuracy", 0.0)
        val_mcc = val_metrics.get("mcc", 0.0)
        val_auc = val_metrics.get("roc_auc", 0.0)
        val_f1 = val_metrics.get("f1", 0.0)

        # Update best results
        if val_acc > result.best_val_acc:
            result.best_val_acc = val_acc
            result.best_val_mcc = val_mcc
            result.best_val_auc = val_auc
            result.best_val_f1 = val_f1
            result.best_epoch = epoch + 1
            best_model_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

            # Save best checkpoint (use 'model_state' key for ensemble script compatibility)
            torch.save({
                "epoch": epoch + 1,
                "model_state": model.state_dict(),  # Compatible with ensemble_infer.py
                "model_state_dict": model.state_dict(),  # Backward compatibility
                "optimizer_state_dict": optimizer.state_dict(),
                "val_acc": val_acc,
                "val_mcc": val_mcc,
                "val_auc": val_auc,
                "best_threshold": 0.5,  # Default threshold for ensemble scripts
                "config": asdict(training_config),
            }, checkpoint_dir / "best_model.pt")

        epoch_time = time.time() - epoch_start
        LOGGER.info(
            "Epoch %d/%d: train_loss=%.4f, val_acc=%.4f, val_mcc=%.4f, "
            "val_auc=%.4f, lr=%.2e (%.1fs)",
            epoch + 1, training_config.epochs,
            train_metrics.get("loss", 0.0),
            val_acc, val_mcc, val_auc,
            scheduler.get_last_lr()[0],
            epoch_time,
        )

        # Early stopping check
        early_stopping(val_acc, model, epoch, checkpoint_dir)
        if early_stopping.should_stop:
            LOGGER.info("Early stopping triggered at epoch %d", epoch + 1)
            result.stopped_early = True
            break

        result.final_train_loss = train_metrics.get("loss", 0.0)

    # Save last checkpoint (use 'model_state' key for ensemble script compatibility)
    torch.save({
        "epoch": epoch + 1,
        "model_state": model.state_dict(),  # Compatible with ensemble_infer.py
        "model_state_dict": model.state_dict(),  # Backward compatibility
        "optimizer_state_dict": optimizer.state_dict(),
        "val_acc": val_acc,
        "val_mcc": val_mcc,
        "best_threshold": 0.5,  # Default threshold for ensemble scripts
        "config": asdict(training_config),
    }, checkpoint_dir / "last_model.pt")

    result.training_time_seconds = time.time() - start_time

    # Restore best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
        LOGGER.info("Restored best model from epoch %d", result.best_epoch)

    return model, result, history


def run_final_evaluation(
    model: nn.Module,
    val_loader: DataLoader,
    device: torch.device,
    cfg: dict,
    variant: str,
) -> Dict[str, Any]:
    """Run comprehensive final evaluation."""
    from amp_utils import get_amp_context
    from sklearn.metrics import (
        confusion_matrix,
        roc_curve,
        auc,
        precision_recall_curve,
        classification_report,
    )
    from sklearn.calibration import calibration_curve

    model.eval()
    amp_ctx = get_amp_context(device, cfg)

    all_probs = []
    all_labels = []
    all_logits = []

    use_tokens = variant not in {"features_only"}
    use_engineered = variant not in {"transformer_only"}

    with torch.no_grad():
        for tokens, engineered, labels in val_loader:
            tokens = tokens.to(device) if use_tokens else None
            engineered = engineered.to(device) if use_engineered else None
            labels = labels.to(device)

            if amp_ctx is not None and amp_ctx.enabled:
                with amp_ctx.autocast():
                    logits = model(tokens, engineered)
            else:
                logits = model(tokens, engineered)

            probs = torch.sigmoid(logits.float()).squeeze(-1)
            all_probs.append(probs.cpu())
            all_labels.append(labels.cpu())
            all_logits.append(logits.float().squeeze(-1).cpu())

    probs = torch.cat(all_probs).numpy()
    labels = torch.cat(all_labels).numpy()
    logits = torch.cat(all_logits).numpy()
    preds = (probs >= 0.5).astype(int)

    # Confusion matrix
    cm = confusion_matrix(labels, preds)
    tn, fp, fn, tp = cm.ravel()

    # ROC curve
    fpr, tpr, roc_thresholds = roc_curve(labels, probs)
    roc_auc = auc(fpr, tpr)

    # Precision-recall curve
    precision, recall, pr_thresholds = precision_recall_curve(labels, probs)
    pr_auc = auc(recall, precision)

    # Calibration curve
    prob_true, prob_pred = calibration_curve(labels, probs, n_bins=10)

    # Classification report
    report = classification_report(labels, preds, output_dict=True)

    # Calculate MCC manually
    mcc_denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / mcc_denom if mcc_denom > 0 else 0.0

    return {
        "accuracy": (tp + tn) / (tp + tn + fp + fn),
        "precision": tp / (tp + fp) if (tp + fp) > 0 else 0.0,
        "recall": tp / (tp + fn) if (tp + fn) > 0 else 0.0,
        "f1": report["1"]["f1-score"] if "1" in report else 0.0,
        "mcc": mcc,
        "specificity": tn / (tn + fp) if (tn + fp) > 0 else 0.0,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "confusion_matrix": cm.tolist(),
        "roc_curve": {"fpr": fpr.tolist(), "tpr": tpr.tolist()},
        "pr_curve": {"precision": precision.tolist(), "recall": recall.tolist()},
        "calibration_curve": {
            "prob_true": prob_true.tolist(),
            "prob_pred": prob_pred.tolist(),
        },
        "classification_report": report,
    }


def save_visualizations(
    history: Dict[str, List[float]],
    eval_results: Dict[str, Any],
    output_dir: Path,
) -> None:
    """Save training curves and evaluation plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("matplotlib not available, skipping visualizations")
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Training curves
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Loss curves
    axes[0, 0].plot(history["train_loss"], label="Train Loss")
    axes[0, 0].plot(history["val_loss"], label="Val Loss")
    axes[0, 0].set_xlabel("Epoch")
    axes[0, 0].set_ylabel("Loss")
    axes[0, 0].set_title("Training and Validation Loss")
    axes[0, 0].legend()
    axes[0, 0].grid(True)

    # Accuracy curves
    axes[0, 1].plot(history["train_acc"], label="Train Acc")
    axes[0, 1].plot(history["val_acc"], label="Val Acc")
    axes[0, 1].set_xlabel("Epoch")
    axes[0, 1].set_ylabel("Accuracy")
    axes[0, 1].set_title("Training and Validation Accuracy")
    axes[0, 1].legend()
    axes[0, 1].grid(True)

    # MCC and AUC
    axes[1, 0].plot(history["val_mcc"], label="Val MCC")
    axes[1, 0].plot(history["val_auc"], label="Val AUC")
    axes[1, 0].set_xlabel("Epoch")
    axes[1, 0].set_ylabel("Score")
    axes[1, 0].set_title("Validation MCC and AUC")
    axes[1, 0].legend()
    axes[1, 0].grid(True)

    # Learning rate
    axes[1, 1].plot(history["learning_rate"])
    axes[1, 1].set_xlabel("Epoch")
    axes[1, 1].set_ylabel("Learning Rate")
    axes[1, 1].set_title("Learning Rate Schedule")
    axes[1, 1].set_yscale("log")
    axes[1, 1].grid(True)

    plt.tight_layout()
    plt.savefig(plots_dir / "training_curves.png", dpi=150)
    plt.close()

    # ROC curve
    fig, ax = plt.subplots(figsize=(8, 8))
    roc = eval_results["roc_curve"]
    ax.plot(roc["fpr"], roc["tpr"], lw=2,
            label=f'ROC (AUC = {eval_results["roc_auc"]:.4f})')
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="Random")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend(loc="lower right")
    ax.grid(True)
    plt.savefig(plots_dir / "roc_curve.png", dpi=150)
    plt.close()

    # Calibration plot
    fig, ax = plt.subplots(figsize=(8, 8))
    cal = eval_results["calibration_curve"]
    ax.plot(cal["prob_pred"], cal["prob_true"], "s-", label="Model")
    ax.plot([0, 1], [0, 1], "k--", label="Perfectly calibrated")
    ax.set_xlabel("Mean Predicted Probability")
    ax.set_ylabel("Fraction of Positives")
    ax.set_title("Calibration Plot")
    ax.legend(loc="lower right")
    ax.grid(True)
    plt.savefig(plots_dir / "calibration_plot.png", dpi=150)
    plt.close()

    # Confusion matrix
    fig, ax = plt.subplots(figsize=(8, 6))
    cm = np.array(eval_results["confusion_matrix"])
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Negative", "Positive"])
    ax.set_yticklabels(["Negative", "Positive"])
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Actual")
    ax.set_title("Confusion Matrix")

    # Add text annotations
    for i in range(2):
        for j in range(2):
            text = ax.text(j, i, cm[i, j], ha="center", va="center",
                          color="white" if cm[i, j] > cm.max() / 2 else "black",
                          fontsize=14)

    plt.colorbar(im)
    plt.tight_layout()
    plt.savefig(plots_dir / "confusion_matrix.png", dpi=150)
    plt.close()

    LOGGER.info("Saved visualizations to %s", plots_dir)


def print_comparison_table(result: TrainingResult) -> str:
    """Print comparison table with original model."""
    original = ORIGINAL_MODEL_BASELINE

    # Calculate changes
    acc_change = (result.best_val_acc - original["val_acc"]) * 100
    mcc_change = result.best_val_mcc - original["val_mcc"]
    auc_change = result.best_val_auc - original["val_auc"]
    param_change = (result.param_count - original["param_count"]) / original["param_count"] * 100
    time_change = (
        (result.training_time_seconds / 60 - original["training_time_minutes"])
        / original["training_time_minutes"] * 100
    )

    table = """
================================================================================
                         MODEL COMPARISON
================================================================================

| Metric          | Original      | Simplified    | Change        |
|-----------------|---------------|---------------|---------------|
| Val Accuracy    | {:.2%}       | {:.2%}       | {:+.2f}%      |
| Val MCC         | {:.4f}        | {:.4f}        | {:+.4f}       |
| Val AUC         | {:.4f}        | {:.4f}        | {:+.4f}       |
| Parameters      | {:.2f}M       | {:.2f}M       | {:.1f}%       |
| Train Time      | {:.1f}min     | {:.1f}min     | {:.1f}%       |

================================================================================
""".format(
        original["val_acc"],
        result.best_val_acc,
        acc_change,
        original["val_mcc"],
        result.best_val_mcc,
        mcc_change,
        original["val_auc"],
        result.best_val_auc,
        auc_change,
        original["param_count"] / 1e6,
        result.param_count / 1e6,
        param_change,
        original["training_time_minutes"],
        result.training_time_seconds / 60,
        time_change,
    )

    return table


def get_model_architecture_summary(model: nn.Module) -> str:
    """Generate model architecture summary."""
    lines = ["Model Architecture Summary", "=" * 50]

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    lines.append(f"Total parameters: {total_params:,}")
    lines.append(f"Trainable parameters: {trainable_params:,}")
    lines.append(f"Non-trainable parameters: {total_params - trainable_params:,}")
    lines.append("")

    # Module-wise breakdown
    lines.append("Module-wise parameter count:")
    lines.append("-" * 50)

    for name, module in model.named_children():
        params = sum(p.numel() for p in module.parameters())
        lines.append(f"  {name}: {params:,} ({params/total_params*100:.1f}%)")

    return "\n".join(lines)


def save_results(
    result: TrainingResult,
    eval_results: Dict[str, Any],
    history: Dict[str, List[float]],
    training_config: TrainingConfig,
    model: nn.Module,
    output_dir: Path,
) -> None:
    """Save all training results and artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Final metrics JSON
    final_metrics = {
        "timestamp": datetime.now().isoformat(),
        "training_config": asdict(training_config),
        "training_result": result.to_dict(),
        "evaluation_metrics": {
            "accuracy": eval_results["accuracy"],
            "precision": eval_results["precision"],
            "recall": eval_results["recall"],
            "f1": eval_results["f1"],
            "mcc": eval_results["mcc"],
            "specificity": eval_results["specificity"],
            "roc_auc": eval_results["roc_auc"],
            "pr_auc": eval_results["pr_auc"],
        },
        "confusion_matrix": eval_results["confusion_matrix"],
        "training_history": history,
    }

    with open(output_dir / "final_metrics.json", "w") as f:
        json.dump(final_metrics, f, indent=2)

    # Model architecture summary
    arch_summary = get_model_architecture_summary(model)
    with open(output_dir / "model_architecture.txt", "w") as f:
        f.write(arch_summary)

    # Comparison table
    comparison = print_comparison_table(result)
    with open(output_dir / "comparison.txt", "w") as f:
        f.write(comparison)

    LOGGER.info("Saved results to %s", output_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Train final model with best hyperparameters"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to base config file",
    )
    parser.add_argument(
        "--best-config",
        type=Path,
        default=None,
        help="Path to best hyperparameter config (YAML)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/final_training"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Number of training epochs (default: 100)",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=15,
        help="Early stopping patience (default: 15)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )
    parser.add_argument(
        "--use-standard-model",
        action="store_true",
        help="Use legacy GeneWhispererStage1 model (compatible with ensemble scripts)",
    )
    parser.add_argument(
        "--load-mlm-weights",
        action="store_true",
        default=True,
        help="Load pretrained MLM weights (default: True)",
    )
    parser.add_argument(
        "--no-load-mlm-weights",
        action="store_false",
        dest="load_mlm_weights",
        help="Disable loading pretrained MLM weights",
    )
    parser.add_argument(
        "--mlm-transfer-mode",
        type=str,
        default="embed_only",
        choices=["embed_only", "embed_plus_adapter", "none"],
        help="MLM weight transfer mode (default: embed_only)",
    )
    args = parser.parse_args()

    # Load configurations
    cfg = load_base_config(args.config)

    # Load best hyperparameters
    if args.best_config and args.best_config.exists():
        with open(args.best_config) as f:
            best_hp = yaml.safe_load(f)
        LOGGER.info("Loaded best config from %s", args.best_config)
    else:
        best_hp = BEST_CONFIG.copy()
        LOGGER.info("Using default best config (from tuning screenshot)")

    # Create training config
    training_config = TrainingConfig.from_dict({
        **best_hp,
        "epochs": args.epochs,
        "patience": args.patience,
        "seed": args.seed,
    })

    LOGGER.info("=" * 60)
    LOGGER.info("FINAL MODEL TRAINING")
    LOGGER.info("=" * 60)
    LOGGER.info("Variant: %s", training_config.variant)
    LOGGER.info("Learning Rate: %s", training_config.lr)
    LOGGER.info("Dropout: %s", training_config.dropout)
    LOGGER.info("Classifier Hidden: %s", training_config.classifier_hidden)
    LOGGER.info("Pooling Type: %s", training_config.pooling_type)
    LOGGER.info("Weight Decay: %s", training_config.weight_decay)
    LOGGER.info("Epochs: %s", training_config.epochs)
    LOGGER.info("Patience: %s", training_config.patience)
    LOGGER.info("Use Standard Model: %s", args.use_standard_model)
    LOGGER.info("Load MLM Weights: %s", args.load_mlm_weights)
    LOGGER.info("MLM Transfer Mode: %s", args.mlm_transfer_mode)
    LOGGER.info("=" * 60)

    # Set seed
    from seed_utils import set_global_seed
    set_global_seed(training_config.seed)

    # Select device
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Device: %s", device)

    # Build data loaders
    from dataset import build_dataloaders, KmerVocabulary
    loaders = build_dataloaders(cfg)
    train_loader = loaders["stage1"]["train"]
    val_loader = loaders["stage1"]["val"]
    LOGGER.info("Train samples: %d, Val samples: %d",
                len(train_loader.dataset), len(val_loader.dataset))

    # Load vocabulary if using standard model
    vocab = None
    if args.use_standard_model:
        kmer = int(cfg.get("kmer", 6))
        vocab_path = cfg.get("mlm_vocab_path") or cfg.get("mlm_vocab_by_k", {}).get(kmer)
        if vocab_path and Path(vocab_path).exists():
            vocab = KmerVocabulary.load(Path(vocab_path))
            LOGGER.info("Loaded vocabulary from %s (k=%d, vocab_size=%d)",
                        vocab_path, kmer, len(vocab.itos))
        else:
            # Build vocabulary from config
            vocab_cache_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
            vocab_cache_path = vocab_cache_dir / f"mlm_k{kmer}_vocab.json"
            if vocab_cache_path.exists():
                vocab = KmerVocabulary.load(vocab_cache_path)
                LOGGER.info("Loaded vocabulary from cache %s", vocab_cache_path)
            else:
                raise FileNotFoundError(
                    f"Vocabulary not found. Please run MLM pretraining first or "
                    f"provide vocab at {vocab_cache_path}"
                )

    # Create output directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output / timestamp
    output_dir.mkdir(parents=True, exist_ok=True)

    # Train model
    model, result, history = train_final_model(
        cfg=cfg,
        training_config=training_config,
        train_loader=train_loader,
        val_loader=val_loader,
        device=device,
        output_dir=output_dir,
        use_standard_model=args.use_standard_model,
        vocab=vocab,
        load_mlm_weights=args.load_mlm_weights,
        mlm_transfer_mode=args.mlm_transfer_mode,
    )

    # Run final evaluation
    LOGGER.info("Running final evaluation...")
    eval_variant = "combined" if args.use_standard_model else training_config.variant
    eval_results = run_final_evaluation(
        model=model,
        val_loader=val_loader,
        device=device,
        cfg=cfg,
        variant=eval_variant,
    )

    # Save visualizations
    save_visualizations(history, eval_results, output_dir)

    # Save all results
    save_results(result, eval_results, history, training_config, model, output_dir)

    # Print summary
    print("\n" + "=" * 60)
    print("FINAL TRAINING COMPLETE")
    print("=" * 60)
    print(f"\nBest Epoch: {result.best_epoch}")
    print(f"Best Val Accuracy: {result.best_val_acc:.4f} ({result.best_val_acc*100:.2f}%)")
    print(f"Best Val MCC: {result.best_val_mcc:.4f}")
    print(f"Best Val AUC: {result.best_val_auc:.4f}")
    print(f"Best Val F1: {result.best_val_f1:.4f}")
    print(f"Parameters: {result.param_count:,} ({result.param_count/1e6:.2f}M)")
    print(f"Training Time: {result.training_time_seconds/60:.1f} minutes")
    if result.stopped_early:
        print("(Stopped early due to no improvement)")

    # Print comparison table
    comparison_table = print_comparison_table(result)
    print(comparison_table)

    print(f"\nResults saved to: {output_dir}")
    print(f"  - Best checkpoint: {output_dir}/checkpoints/best_model.pt")
    print(f"  - Training curves: {output_dir}/plots/training_curves.png")
    print(f"  - ROC curve: {output_dir}/plots/roc_curve.png")
    print(f"  - Final metrics: {output_dir}/final_metrics.json")

    return result


if __name__ == "__main__":
    main()
