"""Metrics helpers for binary classification."""
from __future__ import annotations

from typing import Dict

import numpy as np


def compute_binary_metrics(logits: np.ndarray, labels: np.ndarray) -> Dict[str, float]:
    probs = 1 / (1 + np.exp(-logits))
    preds = (probs >= 0.5).astype(np.int64)
    labels = labels.astype(np.int64)

    tp = int(((preds == 1) & (labels == 1)).sum())
    tn = int(((preds == 0) & (labels == 0)).sum())
    fp = int(((preds == 1) & (labels == 0)).sum())
    fn = int(((preds == 0) & (labels == 1)).sum())

    accuracy = (tp + tn) / max(1, tp + tn + fp + fn)
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2 * precision * recall / max(1e-8, precision + recall)

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }
