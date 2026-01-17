#!/usr/bin/env python3
"""
Variant Comparison Experiment Runner for Gene Whisperer Stage 1.

Runs fair comparison of all model variants with:
- Same random seed, train/val split, and hyperparameters
- Consistent metrics tracking across variants
- Results saved to CSV, JSON, and checkpoints

Usage:
    python run_variant_comparison.py --config config.yaml --epochs 50 --output results/variant_comparison/
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import logging
import os
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import torch
import torch.nn as nn
import yaml

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

from amp_utils import get_amp_context, create_grad_scaler
from dataset import build_dataloaders
from early_stopping import EarlyStopping
from model import GeneWhispererStage1, create_model_variant
from seed_utils import set_global_seed
from train_stage1 import (
    run_epoch,
    save_checkpoint,
    calculate_metrics,
    load_config,
)

LOGGER = logging.getLogger("gene_whisperer.variant_comparison")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Variants to compare
VARIANTS = ["transformer_only", "features_only", "combined", "original"]


@dataclass
class VariantResult:
    """Results for a single variant training run."""

    variant: str
    best_val_acc: float = 0.0
    best_val_mcc: float = 0.0
    best_val_auc: float = 0.0
    best_epoch: int = 0
    param_count: int = 0
    training_time_seconds: float = 0.0
    final_train_loss: float = 0.0
    best_threshold: float = 0.5
    training_curves: Dict[str, List[float]] = field(default_factory=dict)


def count_parameters(model: nn.Module) -> int:
    """Count trainable parameters in a model."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def create_original_model(cfg: dict, device: torch.device) -> nn.Module:
    """Create the original full GeneWhispererStage1 model."""
    kmer = int(cfg.get("kmer", 6))
    max_bp_len = int(cfg.get("max_bp_len", 81))
    max_seq_len = max_bp_len - kmer + 1

    return GeneWhispererStage1(
        vocab_size=int(cfg.get("vocab_size", 4099)),
        embedding_dim=int(cfg.get("embedding_dim", 384)),
        kmer=kmer,
        num_layers=int(cfg.get("transformer_layers", 12)),
        num_heads=int(cfg.get("transformer_heads", 12)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1536)),
        dropout=float(cfg.get("transformer_dropout", 0.12)),
        max_seq_len=max_seq_len,
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 384)),
        tcn_levels=int(cfg.get("tcn_levels", 5)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 512)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 384)),
        fusion_hidden=int(cfg.get("fusion_hidden", 512)),
        pad_token_id=int(cfg.get("pad_token_id", 4098)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", True)),
        relative_position_num_buckets=int(cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", True)),
        glu_activation=str(cfg.get("glu_activation", "gelu")),
    ).to(device)


def create_model_for_variant(
    variant: str, cfg: dict, device: torch.device
) -> nn.Module:
    """Create model for a specific variant."""
    if variant == "original":
        return create_original_model(cfg, device)
    return create_model_variant(variant, cfg).to(device)


def train_single_variant(
    variant: str,
    cfg: dict,
    train_loader: torch.utils.data.DataLoader,
    val_loader: torch.utils.data.DataLoader,
    device: torch.device,
    output_dir: Path,
    epochs: int,
    patience: int = 10,
) -> VariantResult:
    """Train a single variant and track metrics."""
    LOGGER.info("=" * 60)
    LOGGER.info("Training variant: %s", variant)
    LOGGER.info("=" * 60)

    result = VariantResult(variant=variant)
    result.training_curves = {
        "train_loss": [],
        "train_acc": [],
        "val_acc": [],
        "val_mcc": [],
        "val_auc": [],
    }

    # Create model
    model = create_model_for_variant(variant, cfg, device)
    result.param_count = count_parameters(model)
    LOGGER.info("Variant %s: %d parameters", variant, result.param_count)

    # Setup optimizer and scheduler
    lr = float(cfg.get("stage1_lr", cfg.get("lr", 2e-5)))
    weight_decay = float(cfg.get("weight_decay", 0.01))
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=lr,
        weight_decay=weight_decay,
    )

    # Warmup + cosine annealing
    warmup_ratio = float(cfg.get("warmup_ratio", 0.1))
    total_steps = epochs * len(train_loader)
    warmup_steps = int(total_steps * warmup_ratio)

    def lr_lambda(step: int) -> float:
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return max(0.1, 0.5 * (1 + torch.cos(torch.tensor(progress * 3.14159)).item()))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    # Loss function
    label_smoothing = float(cfg.get("stage1_label_smoothing", cfg.get("label_smoothing", 0.0)))
    if label_smoothing > 0:
        criterion = nn.BCEWithLogitsLoss(
            pos_weight=torch.ones(1).to(device),
            reduction="mean",
        )
    else:
        criterion = nn.BCEWithLogitsLoss()

    # Early stopping
    early_stopping = EarlyStopping(
        patience=patience,
        mode="max",
        min_delta=0.001,
        min_epochs=5,
        restore_best=False,  # We handle checkpointing ourselves
        verbose=False,
    )

    # AMP context
    amp_ctx = get_amp_context(device, cfg)
    amp_enabled = bool(cfg.get("amp_enabled", True))
    scaler = create_grad_scaler(device) if amp_enabled else None

    # Training config
    grad_accum_steps = int(cfg.get("grad_accum_steps", 4))
    grad_clip_norm = float(cfg.get("grad_clip_norm", 1.0))
    use_mixup = bool(cfg.get("use_mixup", False))
    mixup_alpha = float(cfg.get("mixup_alpha", 0.1))
    param_norm_warn = float(cfg.get("param_norm_warn", 400))
    param_norm_cap = float(cfg.get("param_norm_cap", 500))

    # Map variant name for run_epoch compatibility
    model_variant_key = variant if variant != "original" else "combined"

    start_time = time.time()
    best_checkpoint_path = output_dir / f"{variant}_best.pt"

    for epoch in range(epochs):
        # Training
        train_metrics = run_epoch(
            loader=train_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            scheduler=scheduler,
            desc=f"{variant} Train E{epoch+1}",
            grad_accum_steps=grad_accum_steps,
            grad_clip_norm=grad_clip_norm,
            use_mixup=use_mixup,
            mixup_alpha=mixup_alpha,
            amp_ctx=amp_ctx,
            scaler=scaler,
            param_norm_warn=param_norm_warn,
            param_norm_cap=param_norm_cap,
            model_variant=model_variant_key,
        )

        # Validation
        val_metrics = run_epoch(
            loader=val_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=None,
            desc=f"{variant} Val E{epoch+1}",
            amp_ctx=amp_ctx,
            model_variant=model_variant_key,
        )

        # Track curves
        result.training_curves["train_loss"].append(train_metrics.get("loss", 0.0))
        result.training_curves["train_acc"].append(train_metrics.get("accuracy", 0.0))
        result.training_curves["val_acc"].append(val_metrics.get("accuracy", 0.0))
        result.training_curves["val_mcc"].append(val_metrics.get("mcc", 0.0))
        result.training_curves["val_auc"].append(val_metrics.get("roc_auc", 0.0))

        val_acc = val_metrics.get("accuracy", 0.0)
        val_mcc = val_metrics.get("mcc", 0.0)
        val_auc = val_metrics.get("roc_auc", 0.0)

        LOGGER.info(
            "Epoch %d/%d - Val Acc: %.4f, Val MCC: %.4f, Val AUC: %.4f",
            epoch + 1, epochs, val_acc, val_mcc, val_auc,
        )

        # Track best
        if val_acc > result.best_val_acc:
            result.best_val_acc = val_acc
            result.best_val_mcc = val_mcc
            result.best_val_auc = val_auc
            result.best_epoch = epoch + 1
            result.best_threshold = val_metrics.get("best_threshold", 0.5)
            result.final_train_loss = train_metrics.get("loss", 0.0)

            save_checkpoint(
                path=best_checkpoint_path,
                model=model,
                optimizer=optimizer,
                scheduler=scheduler,
                epoch=epoch + 1,
                metrics=val_metrics,
                model_variant=variant,
                use_simplified_architecture=(variant != "original"),
            )

        # Early stopping check
        if early_stopping(val_acc, model, epoch + 1, checkpoint_dir=None):
            LOGGER.info("Early stopping triggered at epoch %d", epoch + 1)
            break

    result.training_time_seconds = time.time() - start_time
    LOGGER.info(
        "Variant %s complete: Best Val Acc=%.4f at epoch %d (%.1fs)",
        variant, result.best_val_acc, result.best_epoch, result.training_time_seconds,
    )

    return result


def save_results_csv(results: List[VariantResult], path: Path) -> None:
    """Save results summary to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Variant", "Val Acc", "Val MCC", "Val AUC",
            "Params", "Time (s)", "Best Epoch", "Best Threshold"
        ])
        for r in results:
            writer.writerow([
                r.variant,
                f"{r.best_val_acc:.4f}",
                f"{r.best_val_mcc:.4f}",
                f"{r.best_val_auc:.4f}",
                r.param_count,
                f"{r.training_time_seconds:.1f}",
                r.best_epoch,
                f"{r.best_threshold:.4f}",
            ])
    LOGGER.info("Saved CSV results to %s", path)


def save_results_json(results: List[VariantResult], path: Path) -> None:
    """Save full results including training curves to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "timestamp": datetime.now().isoformat(),
        "variants": [asdict(r) for r in results],
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    LOGGER.info("Saved JSON results to %s", path)


def print_summary_table(results: List[VariantResult]) -> None:
    """Print formatted summary table."""
    print("\n" + "=" * 90)
    print("VARIANT COMPARISON RESULTS")
    print("=" * 90)
    print(f"{'Variant':<18} {'Val Acc':>10} {'Val MCC':>10} {'Val AUC':>10} "
          f"{'Params':>12} {'Time (s)':>10} {'Best Epoch':>12}")
    print("-" * 90)
    for r in results:
        print(f"{r.variant:<18} {r.best_val_acc:>10.4f} {r.best_val_mcc:>10.4f} "
              f"{r.best_val_auc:>10.4f} {r.param_count:>12,} "
              f"{r.training_time_seconds:>10.1f} {r.best_epoch:>12}")
    print("=" * 90)


def print_statistical_comparison(results: List[VariantResult]) -> None:
    """Print statistical analysis of variant comparison."""
    print("\n" + "=" * 70)
    print("STATISTICAL ANALYSIS")
    print("=" * 70)

    results_dict = {r.variant: r for r in results}
    threshold = 0.01  # 1% improvement threshold

    combined = results_dict.get("combined")
    transformer = results_dict.get("transformer_only")
    features = results_dict.get("features_only")
    original = results_dict.get("original")

    if combined and transformer:
        diff = combined.best_val_acc - transformer.best_val_acc
        if diff > threshold:
            print(f"✓ Engineered features help: combined > transformer_only by {diff:.2%}")
        elif diff < -threshold:
            print(f"✗ Engineered features hurt: combined < transformer_only by {-diff:.2%}")
        else:
            print(f"≈ Engineered features neutral: combined ≈ transformer_only (diff: {diff:.2%})")

    if combined and features:
        diff = combined.best_val_acc - features.best_val_acc
        if diff > threshold:
            print(f"✓ Transformer helps: combined > features_only by {diff:.2%}")
        elif diff < -threshold:
            print(f"✗ Transformer hurts: combined < features_only by {-diff:.2%}")
        else:
            print(f"≈ Transformer neutral: combined ≈ features_only (diff: {diff:.2%})")

    if combined and transformer and features:
        max_simple = max(transformer.best_val_acc, features.best_val_acc)
        synergy = combined.best_val_acc - max_simple
        if synergy > threshold:
            print(f"✓ Synergy detected: combined > max(transformer, features) by {synergy:.2%}")
        elif synergy < -threshold:
            print(f"✗ Interference: combined < max(transformer, features) by {-synergy:.2%}")
        else:
            print(f"≈ No synergy: combined ≈ max(transformer, features) (diff: {synergy:.2%})")

    if original and combined:
        diff = original.best_val_acc - combined.best_val_acc
        param_ratio = original.param_count / max(combined.param_count, 1)
        print(f"\nOriginal vs Combined:")
        print(f"  Accuracy difference: {diff:+.2%}")
        print(f"  Parameter ratio: {param_ratio:.2f}x ({original.param_count:,} vs {combined.param_count:,})")

    # Find best variant
    if results:
        best = max(results, key=lambda r: r.best_val_acc)
        print(f"\n★ Best variant: {best.variant} (Val Acc: {best.best_val_acc:.4f})")

    print("=" * 70)


def run_experiment(
    config_path: Path,
    output_dir: Path,
    epochs: int = 50,
    patience: int = 10,
    seed: int = 1337,
    variants: Optional[List[str]] = None,
) -> List[VariantResult]:
    """Run the full variant comparison experiment."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load and set up config
    cfg = load_config(config_path)
    cfg["seed"] = seed
    set_global_seed(seed)

    LOGGER.info("Starting variant comparison experiment")
    LOGGER.info("Config: %s", config_path)
    LOGGER.info("Output: %s", output_dir)
    LOGGER.info("Epochs: %d, Patience: %d, Seed: %d", epochs, patience, seed)

    # Device setup
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Device: %s", device)

    # Build dataloaders (shared across all variants)
    LOGGER.info("Building dataloaders...")
    loaders = build_dataloaders(cfg)
    train_loader = loaders["stage1"]["train"]
    val_loader = loaders["stage1"]["val"]

    if train_loader is None or val_loader is None:
        raise RuntimeError("Failed to create dataloaders")

    LOGGER.info("Train batches: %d, Val batches: %d", len(train_loader), len(val_loader))

    # Select variants to run
    variants_to_run = variants or VARIANTS

    # Run each variant
    results: List[VariantResult] = []
    for variant in variants_to_run:
        if variant not in VARIANTS:
            LOGGER.warning("Skipping unknown variant: %s", variant)
            continue

        # Reset seed before each variant for fair comparison
        set_global_seed(seed)

        result = train_single_variant(
            variant=variant,
            cfg=cfg,
            train_loader=train_loader,
            val_loader=val_loader,
            device=device,
            output_dir=output_dir,
            epochs=epochs,
            patience=patience,
        )
        results.append(result)

    # Save results
    save_results_csv(results, output_dir / "results.csv")
    save_results_json(results, output_dir / "results.json")

    # Print summary
    print_summary_table(results)
    print_statistical_comparison(results)

    return results


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run variant comparison experiment for Gene Whisperer Stage 1"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to config file (default: config.yaml)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/variant_comparison"),
        help="Output directory (default: results/variant_comparison)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Number of training epochs (default: 50)",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=10,
        help="Early stopping patience (default: 10)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1337,
        help="Random seed (default: 1337)",
    )
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=VARIANTS,
        default=None,
        help="Specific variants to run (default: all)",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for variant comparison experiment."""
    args = parse_args()

    # Resolve paths relative to script location
    script_dir = Path(__file__).parent
    config_path = args.config
    if not config_path.is_absolute():
        config_path = script_dir / config_path

    output_path = args.output
    if not output_path.is_absolute():
        output_path = script_dir / output_path

    run_experiment(
        config_path=config_path,
        output_dir=output_path,
        epochs=args.epochs,
        patience=args.patience,
        seed=args.seed,
        variants=args.variants,
    )


if __name__ == "__main__":
    main()
