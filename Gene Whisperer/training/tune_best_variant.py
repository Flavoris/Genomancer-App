#!/usr/bin/env python3
"""
Hyperparameter Tuning for Best Model Variant.

Runs a hyperparameter search over:
- Learning rate: [1e-5, 2e-5, 5e-5, 1e-4]
- Dropout: [0.1, 0.15, 0.2]
- Classifier hidden dim: [128, 256, 512]
- Pooling type: ["mean", "cls", "attention"]
- Weight decay: [0.01, 0.02, 0.05]

Uses validation accuracy as optimization target.
Outputs best config as YAML snippet.

Usage:
    python tune_best_variant.py --results results/variant_comparison/results.json
    python tune_best_variant.py --variant combined --trials 20
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import logging
import os
import random
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn
import yaml

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

from amp_utils import get_amp_context, create_grad_scaler
from dataset import build_dataloaders
from early_stopping import EarlyStopping
from seed_utils import set_global_seed
from train_stage1 import run_epoch, load_config, get_mlm_checkpoint

LOGGER = logging.getLogger("gene_whisperer.hyperparameter_tuning")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


# Hyperparameter search space as defined in requirements
HYPERPARAM_SPACE = {
    "lr": [1e-5, 2e-5, 5e-5, 1e-4],
    "dropout": [0.1, 0.15, 0.2],
    "classifier_hidden": [128, 256, 512],
    "pooling_type": ["mean", "cls", "attention"],
    "weight_decay": [0.01, 0.02, 0.05],
}


@dataclass
class TuningResult:
    """Results for a single hyperparameter configuration."""

    config_id: int
    lr: float
    dropout: float
    classifier_hidden: int
    pooling_type: str
    weight_decay: float
    best_val_acc: float = 0.0
    best_val_mcc: float = 0.0
    best_epoch: int = 0
    training_time_seconds: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def load_comparison_results(results_path: Path) -> Dict[str, Any]:
    """Load variant comparison results to find best variant."""
    with open(results_path) as f:
        data = json.load(f)
    return data


def get_best_variant(results: Dict[str, Any]) -> str:
    """Determine the best variant from comparison results."""
    variants = results.get("variants", [])
    if not variants:
        raise ValueError("No variant results found in comparison file")

    # Sort by validation accuracy descending
    sorted_variants = sorted(variants, key=lambda x: x.get("best_val_acc", 0), reverse=True)
    best = sorted_variants[0]

    LOGGER.info("Best variant from comparison: %s (val_acc=%.4f)",
                best["variant"], best["best_val_acc"])
    return best["variant"]


def generate_search_configs(
    search_type: str = "random",
    num_trials: int = 20,
    seed: int = 42,
) -> List[Dict[str, Any]]:
    """Generate hyperparameter configurations for search."""

    if search_type == "grid":
        # Full grid search
        keys = list(HYPERPARAM_SPACE.keys())
        values = [HYPERPARAM_SPACE[k] for k in keys]
        configs = []
        for combo in itertools.product(*values):
            configs.append(dict(zip(keys, combo)))
        LOGGER.info("Grid search: %d total configurations", len(configs))
        return configs

    # Random search
    random.seed(seed)
    configs = []
    for _ in range(num_trials):
        config = {
            key: random.choice(vals) for key, vals in HYPERPARAM_SPACE.items()
        }
        configs.append(config)

    LOGGER.info("Random search: %d configurations sampled", len(configs))
    return configs


def create_tunable_model(
    variant: str,
    cfg: dict,
    hp_config: Dict[str, Any],
    device: torch.device,
) -> nn.Module:
    """
    Create model with tunable hyperparameters.

    Supports tuning:
    - dropout: Applied to transformer and classifier
    - classifier_hidden: Hidden dimension of classifier head
    - pooling_type: mean, cls, or attention pooling
    """
    from model import (
        DNAEncoder,
        AttentionPooling,
        GeneWhispererFeaturesOnly,
    )

    dropout = hp_config["dropout"]
    classifier_hidden = hp_config["classifier_hidden"]
    pooling_type = hp_config["pooling_type"]

    if variant == "features_only":
        # Features-only model doesn't need pooling
        return GeneWhispererFeaturesOnly(
            engineered_dim=int(cfg.get("engineered_dim", 288)),
            hidden_dim=classifier_hidden,
            second_hidden_dim=classifier_hidden // 2,
            dropout=dropout,
        ).to(device)

    # For transformer-based variants, create custom tunable model
    return _create_tunable_transformer_model(
        variant, cfg, dropout, classifier_hidden, pooling_type, device
    )


def _create_tunable_transformer_model(
    variant: str,
    cfg: dict,
    dropout: float,
    classifier_hidden: int,
    pooling_type: str,
    device: torch.device,
) -> nn.Module:
    """Create a transformer-based model with tunable parameters."""
    from model import DNAEncoder, AttentionPooling

    # Common encoder parameters
    max_token_len = int(cfg.get("max_token_len", 24))
    vocab_size = int(cfg.get("vocab_size", 4096))
    embedding_dim = int(cfg.get("embedding_dim", 384))
    num_layers = int(cfg.get("transformer_layers", 12))
    num_heads = int(cfg.get("transformer_heads", 12))
    ff_dim = int(cfg.get("transformer_ff_dim", 1536))
    pad_token_id = int(cfg.get("pad_token_id", 0))

    # Create encoder
    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=dropout,
        pad_token_id=pad_token_id,
        max_seq_len=max_token_len,
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", True)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", True)),
    )

    # Build tunable model
    if variant == "transformer_only":
        model = TunableTransformerOnly(
            encoder=encoder,
            embedding_dim=embedding_dim,
            classifier_hidden=classifier_hidden,
            classifier_dropout=dropout,
            pooling_type=pooling_type,
            pad_token_id=pad_token_id,
        )
    elif variant == "combined":
        engineered_dim = int(cfg.get("engineered_dim", 288))
        engineered_mlp_hidden = int(cfg.get("engineered_mlp_hidden", 512))
        engineered_mlp_output = int(cfg.get("engineered_mlp_output", 384))

        model = TunableCombined(
            encoder=encoder,
            embedding_dim=embedding_dim,
            engineered_dim=engineered_dim,
            engineered_mlp_hidden=engineered_mlp_hidden,
            engineered_mlp_output=engineered_mlp_output,
            classifier_hidden=classifier_hidden,
            classifier_dropout=dropout,
            pooling_type=pooling_type,
            pad_token_id=pad_token_id,
        )
    else:
        raise ValueError(f"Unsupported variant for tuning: {variant}")

    return model.to(device)


class TunableTransformerOnly(nn.Module):
    """Transformer-only model with configurable pooling and classifier."""

    def __init__(
        self,
        encoder: nn.Module,
        embedding_dim: int,
        classifier_hidden: int,
        classifier_dropout: float,
        pooling_type: str,
        pad_token_id: Optional[int] = None,
    ):
        super().__init__()
        self.pad_token_id = pad_token_id
        self._full_encoder = encoder
        self.embedding = encoder.embedding
        self.embedding_dim = embedding_dim
        self.pooling_type = pooling_type

        # Pooling layer
        if pooling_type == "attention":
            from model import AttentionPooling
            self.pool = AttentionPooling(embedding_dim, num_heads=4)
        else:
            self.pool = None

        # Classifier
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim, classifier_hidden),
            nn.LayerNorm(classifier_hidden),
            nn.GELU(),
            nn.Dropout(classifier_dropout),
            nn.Linear(classifier_hidden, classifier_hidden // 2),
            nn.LayerNorm(classifier_hidden // 2),
            nn.GELU(),
            nn.Dropout(classifier_dropout * 0.8),
            nn.Linear(classifier_hidden // 2, 1),
        )

    @property
    def encoder(self):
        return self._full_encoder

    def load_pretrained_weights(self, checkpoint_path, strict=False, transfer_mode="embed_only"):
        if transfer_mode != "none":
            self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is None:
            return None
        return tokens.eq(self.pad_token_id)

    def _pool(self, features: torch.Tensor, mask: Optional[torch.Tensor]) -> torch.Tensor:
        """Apply pooling based on pooling_type."""
        if self.pooling_type == "attention":
            return self.pool(features, mask)
        elif self.pooling_type == "cls":
            # CLS token pooling - first token
            return features[:, 0, :]
        else:  # mean
            if mask is not None:
                # Mask out padding
                mask_expanded = mask.unsqueeze(-1).float()
                features = features * (1 - mask_expanded)
                lengths = (~mask).sum(dim=1, keepdim=True).float().clamp(min=1)
                return features.sum(dim=1) / lengths
            return features.mean(dim=1)

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ):
        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        pooled = self._pool(token_features, padding_mask)
        logits = self.classifier(pooled)
        if return_logits:
            return torch.sigmoid(logits), logits
        return logits


class TunableCombined(nn.Module):
    """Combined model with configurable pooling and classifier."""

    def __init__(
        self,
        encoder: nn.Module,
        embedding_dim: int,
        engineered_dim: int,
        engineered_mlp_hidden: int,
        engineered_mlp_output: int,
        classifier_hidden: int,
        classifier_dropout: float,
        pooling_type: str,
        pad_token_id: Optional[int] = None,
    ):
        super().__init__()
        self.pad_token_id = pad_token_id
        self._full_encoder = encoder
        self.embedding = encoder.embedding
        self.embedding_dim = embedding_dim
        self.engineered_dim = engineered_dim
        self.pooling_type = pooling_type

        # Pooling layer
        if pooling_type == "attention":
            from model import AttentionPooling
            self.pool = AttentionPooling(embedding_dim, num_heads=4)
        else:
            self.pool = None

        # Engineered features MLP
        self.engineered_mlp = nn.Sequential(
            nn.Linear(engineered_dim, engineered_mlp_hidden),
            nn.GELU(),
            nn.Dropout(classifier_dropout),
            nn.Linear(engineered_mlp_hidden, engineered_mlp_output),
            nn.GELU(),
        )

        # Classifier (takes concatenated transformer + engineered features)
        classifier_input_dim = embedding_dim + engineered_mlp_output
        self.classifier = nn.Sequential(
            nn.Linear(classifier_input_dim, classifier_hidden),
            nn.LayerNorm(classifier_hidden),
            nn.GELU(),
            nn.Dropout(classifier_dropout),
            nn.Linear(classifier_hidden, classifier_hidden // 2),
            nn.LayerNorm(classifier_hidden // 2),
            nn.GELU(),
            nn.Dropout(classifier_dropout * 0.8),
            nn.Linear(classifier_hidden // 2, 1),
        )

    @property
    def encoder(self):
        return self._full_encoder

    def load_pretrained_weights(self, checkpoint_path, strict=False, transfer_mode="embed_only"):
        if transfer_mode != "none":
            self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is None:
            return None
        return tokens.eq(self.pad_token_id)

    def _pool(self, features: torch.Tensor, mask: Optional[torch.Tensor]) -> torch.Tensor:
        """Apply pooling based on pooling_type."""
        if self.pooling_type == "attention":
            return self.pool(features, mask)
        elif self.pooling_type == "cls":
            return features[:, 0, :]
        else:  # mean
            if mask is not None:
                mask_expanded = mask.unsqueeze(-1).float()
                features = features * (1 - mask_expanded)
                lengths = (~mask).sum(dim=1, keepdim=True).float().clamp(min=1)
                return features.sum(dim=1) / lengths
            return features.mean(dim=1)

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ):
        if engineered_features is None:
            raise ValueError("Engineered features required for combined variant")

        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        pooled = self._pool(token_features, padding_mask)

        eng_proj = self.engineered_mlp(engineered_features)
        fused = torch.cat([pooled, eng_proj], dim=-1)

        logits = self.classifier(fused)
        if return_logits:
            return torch.sigmoid(logits), logits
        return logits


def count_parameters(model: nn.Module) -> int:
    """Count trainable parameters."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def train_with_config(
    variant: str,
    cfg: dict,
    hp_config: Dict[str, Any],
    train_loader: torch.utils.data.DataLoader,
    val_loader: torch.utils.data.DataLoader,
    device: torch.device,
    epochs: int = 30,
    patience: int = 10,
    config_id: int = 0,
    load_mlm_weights: bool = True,
) -> TuningResult:
    """Train model with a specific hyperparameter configuration."""

    result = TuningResult(
        config_id=config_id,
        lr=hp_config["lr"],
        dropout=hp_config["dropout"],
        classifier_hidden=hp_config["classifier_hidden"],
        pooling_type=hp_config["pooling_type"],
        weight_decay=hp_config["weight_decay"],
    )

    LOGGER.info("-" * 50)
    LOGGER.info("Config %d: lr=%.0e, dropout=%.2f, hidden=%d, pool=%s, wd=%.2f",
                config_id, hp_config["lr"], hp_config["dropout"],
                hp_config["classifier_hidden"], hp_config["pooling_type"],
                hp_config["weight_decay"])

    # Create model
    model = create_tunable_model(variant, cfg, hp_config, device)
    param_count = count_parameters(model)
    LOGGER.info("Parameters: %d", param_count)

    # Load pretrained MLM weights if enabled
    if load_mlm_weights and hasattr(model, "load_pretrained_weights"):
        checkpoint_path = get_mlm_checkpoint(cfg)
        if checkpoint_path is not None:
            LOGGER.info("Loading BPE MLM weights from %s", checkpoint_path)
            model.load_pretrained_weights(
                checkpoint_path=checkpoint_path,
                strict=False,
                transfer_mode="embed_only",
            )
        else:
            LOGGER.warning("No MLM checkpoint found")

    # Optimizer with tuned parameters
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=hp_config["lr"],
        weight_decay=hp_config["weight_decay"],
    )

    # Scheduler: warmup + cosine
    warmup_ratio = float(cfg.get("warmup_ratio", 0.1))
    total_steps = epochs * len(train_loader)
    warmup_steps = int(total_steps * warmup_ratio)

    def lr_lambda(step: int) -> float:
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return max(0.1, 0.5 * (1 + torch.cos(torch.tensor(progress * 3.14159)).item()))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    # Loss and early stopping
    criterion = nn.BCEWithLogitsLoss()
    early_stopping = EarlyStopping(
        patience=patience,
        mode="max",
        min_delta=0.001,
        min_epochs=5,
        restore_best=False,
        verbose=False,
    )

    # AMP
    amp_ctx = get_amp_context(device, cfg)
    amp_enabled = bool(cfg.get("amp_enabled", True))
    scaler = create_grad_scaler(device) if amp_enabled else None

    # Training config
    grad_accum_steps = int(cfg.get("grad_accum_steps", 4))
    grad_clip_norm = float(cfg.get("grad_clip_norm", 1.0))

    start_time = time.time()
    best_val_acc = 0.0
    best_val_mcc = 0.0
    best_epoch = 0

    for epoch in range(epochs):
        # Training epoch
        model.train()
        train_metrics = run_epoch(
            loader=train_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            scheduler=scheduler,
            grad_accum_steps=grad_accum_steps,
            grad_clip_norm=grad_clip_norm,
            use_mixup=False,
            amp_ctx=amp_ctx,
            scaler=scaler,
            model_variant=variant,
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
                grad_accum_steps=1,
                grad_clip_norm=0.0,
                use_mixup=False,
                amp_ctx=amp_ctx,
                scaler=None,
                model_variant=variant,
            )

        val_acc = val_metrics.get("accuracy", 0.0)
        val_mcc = val_metrics.get("mcc", 0.0)

        if val_acc > best_val_acc:
            best_val_acc = val_acc
            best_val_mcc = val_mcc
            best_epoch = epoch + 1

        # Check early stopping (uses __call__, returns True if should stop)
        early_stopping(val_acc, model, epoch)
        if early_stopping.should_stop:
            LOGGER.info("Early stopping at epoch %d", epoch + 1)
            break

    training_time = time.time() - start_time

    result.best_val_acc = best_val_acc
    result.best_val_mcc = best_val_mcc
    result.best_epoch = best_epoch
    result.training_time_seconds = training_time

    LOGGER.info("Result: val_acc=%.4f, val_mcc=%.4f @ epoch %d (%.1fs)",
                best_val_acc, best_val_mcc, best_epoch, training_time)

    return result


def save_results_csv(results: List[TuningResult], output_path: Path) -> None:
    """Save all tuning results to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "config_id", "lr", "dropout", "classifier_hidden",
        "pooling_type", "weight_decay", "best_val_acc",
        "best_val_mcc", "best_epoch", "training_time_seconds"
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r.to_dict())

    LOGGER.info("Results saved to %s", output_path)


def save_results_json(
    results: List[TuningResult],
    best_result: TuningResult,
    variant: str,
    output_path: Path,
) -> None:
    """Save complete tuning results to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = {
        "timestamp": datetime.now().isoformat(),
        "variant": variant,
        "search_space": HYPERPARAM_SPACE,
        "best_config": {
            "lr": best_result.lr,
            "dropout": best_result.dropout,
            "classifier_hidden": best_result.classifier_hidden,
            "pooling_type": best_result.pooling_type,
            "weight_decay": best_result.weight_decay,
            "best_val_acc": best_result.best_val_acc,
            "best_val_mcc": best_result.best_val_mcc,
            "best_epoch": best_result.best_epoch,
        },
        "all_results": [r.to_dict() for r in results],
    }

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

    LOGGER.info("Full results saved to %s", output_path)


def print_best_config_yaml(best_result: TuningResult, variant: str) -> str:
    """Generate YAML snippet for best configuration."""
    yaml_config = {
        "# Best hyperparameters for variant": variant,
        "# Validation accuracy": f"{best_result.best_val_acc:.4f}",
        "# Validation MCC": f"{best_result.best_val_mcc:.4f}",
        "stage1_lr": best_result.lr,
        "transformer_dropout": best_result.dropout,
        "classifier_hidden": best_result.classifier_hidden,
        "pooling_type": best_result.pooling_type,
        "weight_decay": best_result.weight_decay,
    }

    yaml_str = "\n# ===== Best Hyperparameters =====\n"
    yaml_str += f"# Variant: {variant}\n"
    yaml_str += f"# Val Accuracy: {best_result.best_val_acc:.4f}\n"
    yaml_str += f"# Val MCC: {best_result.best_val_mcc:.4f}\n"
    yaml_str += f"# Best Epoch: {best_result.best_epoch}\n"
    yaml_str += "# ================================\n\n"
    yaml_str += f"stage1_lr: {best_result.lr}\n"
    yaml_str += f"transformer_dropout: {best_result.dropout}\n"
    yaml_str += f"classifier_hidden: {best_result.classifier_hidden}\n"
    yaml_str += f"pooling_type: \"{best_result.pooling_type}\"\n"
    yaml_str += f"weight_decay: {best_result.weight_decay}\n"

    return yaml_str


def main():
    parser = argparse.ArgumentParser(
        description="Hyperparameter tuning for best model variant"
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("results/variant_comparison/results.json"),
        help="Path to variant comparison results JSON",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Override variant to tune (default: best from results)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to config file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/hyperparameter_tuning/"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--trials",
        type=int,
        default=20,
        help="Number of random search trials (default: 20)",
    )
    parser.add_argument(
        "--search",
        type=str,
        choices=["random", "grid"],
        default="random",
        help="Search type: random or grid",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=30,
        help="Max epochs per trial (default: 30)",
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
        default=42,
        help="Random seed",
    )
    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    # Determine variant to tune
    if args.variant:
        variant = args.variant
        LOGGER.info("Using specified variant: %s", variant)
    else:
        if not args.results.exists():
            LOGGER.warning(
                "Comparison results not found at %s, defaulting to 'combined'",
                args.results
            )
            variant = "combined"
        else:
            comparison_results = load_comparison_results(args.results)
            variant = get_best_variant(comparison_results)

    LOGGER.info("=" * 60)
    LOGGER.info("Hyperparameter Tuning for Variant: %s", variant)
    LOGGER.info("=" * 60)

    # Set seed
    set_global_seed(args.seed)

    # Select device
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Device: %s", device)

    # Build data loaders
    loaders = build_dataloaders(cfg)
    train_loader = loaders["stage1"]["train"]
    val_loader = loaders["stage1"]["val"]
    LOGGER.info("Train samples: %d, Val samples: %d",
                len(train_loader.dataset), len(val_loader.dataset))

    # Generate search configs
    search_configs = generate_search_configs(
        search_type=args.search,
        num_trials=args.trials,
        seed=args.seed,
    )

    # Run hyperparameter search
    results: List[TuningResult] = []

    for i, hp_config in enumerate(search_configs):
        # Reset seed for fair comparison
        set_global_seed(args.seed)

        result = train_with_config(
            variant=variant,
            cfg=cfg,
            hp_config=hp_config,
            train_loader=train_loader,
            val_loader=val_loader,
            device=device,
            epochs=args.epochs,
            patience=args.patience,
            config_id=i + 1,
        )
        results.append(result)

    # Find best result
    best_result = max(results, key=lambda r: r.best_val_acc)

    # Save results
    args.output.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    save_results_csv(results, args.output / f"tuning_results_{timestamp}.csv")
    save_results_json(results, best_result, variant, args.output / f"tuning_results_{timestamp}.json")

    # Print summary
    print("\n" + "=" * 60)
    print("HYPERPARAMETER TUNING COMPLETE")
    print("=" * 60)
    print(f"\nBest Configuration (Config #{best_result.config_id}):")
    print(f"  - Learning Rate: {best_result.lr}")
    print(f"  - Dropout: {best_result.dropout}")
    print(f"  - Classifier Hidden: {best_result.classifier_hidden}")
    print(f"  - Pooling Type: {best_result.pooling_type}")
    print(f"  - Weight Decay: {best_result.weight_decay}")
    print(f"\nPerformance:")
    print(f"  - Best Val Accuracy: {best_result.best_val_acc:.4f}")
    print(f"  - Best Val MCC: {best_result.best_val_mcc:.4f}")
    print(f"  - Best Epoch: {best_result.best_epoch}")
    print(f"  - Training Time: {best_result.training_time_seconds:.1f}s")

    # Print YAML snippet
    yaml_snippet = print_best_config_yaml(best_result, variant)
    print("\n" + "=" * 60)
    print("YAML CONFIG SNIPPET (copy to config.yaml):")
    print("=" * 60)
    print(yaml_snippet)

    # Also save YAML snippet to file
    yaml_path = args.output / f"best_config_{timestamp}.yaml"
    with open(yaml_path, "w") as f:
        f.write(yaml_snippet)
    print(f"\nYAML snippet saved to: {yaml_path}")

    return best_result


if __name__ == "__main__":
    main()
