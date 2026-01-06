"""Helpers for resolving and logging Stage 1 training settings."""
from __future__ import annotations

import logging
from typing import Any, Mapping


def resolve_stage1_lr(cfg: Mapping[str, Any], default_lr: float = 1e-4) -> float:
    """Resolve Stage 1 learning rate, falling back to global lr."""
    raw_value = cfg.get("stage1_lr")
    if raw_value is None:
        raw_value = cfg.get("lr", default_lr)
    try:
        return float(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"stage1_lr must be numeric, got {raw_value!r}") from exc


def log_stage1_training_configuration(
    *,
    max_bp_len: int,
    lr: float,
    batch_size: int,
    logger: logging.Logger,
) -> None:
    """Log the resolved Stage 1 training configuration."""
    logger.info("=" * 60)
    logger.info("STAGE 1 TRAINING CONFIGURATION")
    logger.info("=" * 60)
    logger.info("Max BP length: %d", max_bp_len)
    logger.info("Learning rate: %.6f", lr)
    logger.info("Batch size: %d", batch_size)
    logger.info("=" * 60)
