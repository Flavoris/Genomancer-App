"""Diagnostic script: Test whether the transformer learns useful patterns.

Runs three experiments on the validation set:
  A) Full model (transformer + engineered features)
  B) Transformer-only (zero out engineered features)
  C) Features-only (randomize transformer output)

If Experiment B shows ~50-55% accuracy, the transformer is NOT learning
and pre-training needs improvement.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Tuple

import torch
import torch.nn as nn
import yaml

from experiment_runner import (
    build_model_from_config,
    load_checkpoint,
    build_val_dataloader,
    run_experiment_a,
    run_experiment_b,
    run_experiment_c,
    print_analysis,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)
LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Diagnose transformer contribution vs engineered features",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to config YAML (default: config.yaml)",
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        default=None,
        help="Path to Stage 1 checkpoint (overrides config)",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Device to run on (default: auto-detect)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Limit validation samples for faster testing",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Decision threshold (default: use checkpoint best_threshold)",
    )
    return parser.parse_args()


def main() -> None:
    """Run transformer contribution diagnosis."""
    args = parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        LOGGER.error("Config file not found: %s", config_path)
        sys.exit(1)

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    device = _resolve_device(args.device)
    LOGGER.info("Using device: %s", device)

    # Build model
    model = build_model_from_config(config)
    model = model.to(device)

    # Load checkpoint
    checkpoint_path = _resolve_checkpoint_path(args.checkpoint, config)
    threshold = load_checkpoint(model, checkpoint_path, device, args.threshold)
    LOGGER.info("Decision threshold: %.4f", threshold)

    model.eval()

    # Build validation dataloader
    val_loader = build_val_dataloader(config, args.max_samples)
    LOGGER.info("Validation samples: %d", len(val_loader.dataset))

    # Run experiments
    print("\n" + "=" * 60)
    print("TRANSFORMER CONTRIBUTION ANALYSIS")
    print("=" * 60 + "\n")

    acc_a = run_experiment_a(model, val_loader, device, threshold)
    acc_b = run_experiment_b(model, val_loader, device, threshold)
    acc_c = run_experiment_c(model, val_loader, device, threshold)

    print_analysis(acc_a, acc_b, acc_c)


def _resolve_device(device_arg: str | None) -> torch.device:
    """Determine compute device."""
    if device_arg:
        return torch.device(device_arg)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _resolve_checkpoint_path(
    cli_path: str | None,
    config: dict,
) -> Path:
    """Find checkpoint path from CLI arg or config."""
    if cli_path:
        path = Path(cli_path)
    else:
        kmer = int(config.get("kmer", 6))
        ckpt_map = config.get("stage1_ckpt_by_k", {})
        if kmer not in ckpt_map:
            ckpt_map = config.get("checkpoints", {}).get("stage1_ckpt_by_k", {})
        ckpt_str = ckpt_map.get(kmer, f"../artifacts/checkpoints/stage1_k{kmer}.pt")
        path = Path(ckpt_str)

    if not path.exists():
        LOGGER.error("Checkpoint not found: %s", path)
        sys.exit(1)

    LOGGER.info("Loading checkpoint: %s", path)
    return path


if __name__ == "__main__":
    main()
