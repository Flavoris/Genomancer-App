"""
Tests for Stage 1 training settings helpers.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stage1_settings import log_stage1_training_configuration, resolve_stage1_lr


class _ListHandler(logging.Handler):
    """Collect log messages for assertions."""

    def __init__(self) -> None:
        super().__init__()
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def test_resolve_stage1_lr_prefers_stage1_lr() -> None:
    """Stage 1 LR should override global LR when present."""
    cfg = {"stage1_lr": "0.0002", "lr": "0.0001"}
    lr = resolve_stage1_lr(cfg)
    assert abs(lr - 0.0002) < 1e-12, f"Expected 0.0002, got {lr}"


def test_resolve_stage1_lr_falls_back_to_lr() -> None:
    """Stage 1 LR should fall back to global LR when missing."""
    cfg = {"lr": 0.00015}
    lr = resolve_stage1_lr(cfg)
    assert abs(lr - 0.00015) < 1e-12, f"Expected 0.00015, got {lr}"


def test_resolve_stage1_lr_default() -> None:
    """Stage 1 LR should fall back to default when missing."""
    lr = resolve_stage1_lr({})
    assert abs(lr - 0.0001) < 1e-12, f"Expected 0.0001, got {lr}"


def test_log_stage1_training_configuration() -> None:
    """Log output should include the configuration header and values."""
    logger = logging.getLogger("stage1_settings_test")
    logger.setLevel(logging.INFO)
    handler = _ListHandler()
    logger.addHandler(handler)

    try:
        log_stage1_training_configuration(
            max_bp_len=81,
            lr=0.0001,
            batch_size=32,
            logger=logger,
        )
    finally:
        logger.removeHandler(handler)

    joined = "\n".join(handler.messages)
    assert "STAGE 1 TRAINING CONFIGURATION" in joined, "Missing configuration header"
    assert "Max BP length: 81" in joined, "Missing max_bp_len log"
    assert "Learning rate: 0.000100" in joined, "Missing learning rate log"
    assert "Batch size: 32" in joined, "Missing batch size log"


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT VALIDATION: Stage 1 Settings Helpers")
    print("=" * 60)

    test_resolve_stage1_lr_prefers_stage1_lr()
    print("OK: test_resolve_stage1_lr_prefers_stage1_lr")
    test_resolve_stage1_lr_falls_back_to_lr()
    print("OK: test_resolve_stage1_lr_falls_back_to_lr")
    test_resolve_stage1_lr_default()
    print("OK: test_resolve_stage1_lr_default")
    test_log_stage1_training_configuration()
    print("OK: test_log_stage1_training_configuration")

    print("\n" + "=" * 60)
    print("ALL STAGE 1 SETTINGS TESTS PASSED ✓")
    print("=" * 60)
