"""Tests for ensemble weight optimization utilities."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ensemble_weights import (
    load_ensemble_weights,
    optimize_ensemble_weights,
    save_ensemble_weights,
)


def _accuracy_from_probs(probs: np.ndarray, labels: np.ndarray) -> float:
    preds = (probs >= 0.5).astype(np.int32)
    targets = (labels >= 0.5).astype(np.int32)
    return float((preds == targets).mean())


def test_optimize_ensemble_weights_improves_accuracy() -> None:
    labels = np.array([1, 1, 1, 0, 0, 0], dtype=np.float64)

    val_probs = {
        3: np.array([0.9, 0.9, 0.8, 0.2, 0.3, 0.1], dtype=np.float64),
        4: np.array([0.4, 0.3, 0.2, 0.7, 0.6, 0.8], dtype=np.float64),
        5: np.array([0.5, 0.4, 0.3, 0.6, 0.7, 0.6], dtype=np.float64),
    }

    baseline_weights = {k: 1.0 / len(val_probs) for k in val_probs}
    baseline_probs = sum(
        baseline_weights[k] * val_probs[k] for k in val_probs
    )
    baseline_acc = _accuracy_from_probs(baseline_probs, labels)

    weights, metadata = optimize_ensemble_weights(
        val_probs=val_probs,
        val_labels=labels,
        metric="accuracy",
        min_weight=0.05,
        grid_step=0.1,
        max_candidates=5000,
        return_metadata=True,
    )

    assert metadata["metric"] == "accuracy"
    assert pytest.approx(sum(weights.values()), rel=1e-6, abs=1e-6) == 1.0
    assert all(weight >= 0.05 for weight in weights.values())

    optimized_probs = sum(weights[k] * val_probs[k] for k in weights)
    optimized_acc = _accuracy_from_probs(optimized_probs, labels)

    assert optimized_acc >= baseline_acc


def test_save_and_load_weights_roundtrip(tmp_path: Path) -> None:
    weights = {3: 0.15, 4: 0.35, 5: 0.5}
    metadata = {"metric": "mcc", "min_weight": 0.05}
    weights_path = tmp_path / "ensemble_weights.json"

    save_ensemble_weights(weights_path, weights, metadata)
    loaded_weights, loaded_metadata = load_ensemble_weights(weights_path)

    assert loaded_metadata["metric"] == "mcc"
    assert pytest.approx(loaded_weights[3]) == weights[3]
    assert pytest.approx(loaded_weights[4]) == weights[4]
    assert pytest.approx(loaded_weights[5]) == weights[5]
