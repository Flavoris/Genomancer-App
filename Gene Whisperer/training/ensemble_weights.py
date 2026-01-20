"""Ensemble weight optimization and persistence helpers."""
from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

LOGGER = logging.getLogger("gene_whisperer.ensemble_weights")


def _sorted_kmers(kmers: Iterable) -> List:
    """Sort k-mer keys consistently, preferring numeric order when possible."""
    def _sort_key(value: object) -> Tuple[int, object]:
        try:
            return (0, int(value))
        except (TypeError, ValueError):
            return (1, str(value))

    return sorted(kmers, key=_sort_key)


def _normalize_labels(labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.float64).reshape(-1)
    return (labels >= 0.5).astype(np.int32)


def _binary_roc_auc(probabilities: np.ndarray, targets: np.ndarray) -> float:
    """Compute ROC AUC for binary targets without external dependencies."""
    probs = np.asarray(probabilities, dtype=np.float64).reshape(-1)
    y = np.asarray(targets, dtype=np.float64).reshape(-1)
    y = (y >= 0.5).astype(np.int32)

    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return 0.5

    order = np.argsort(probs, kind="mergesort")
    sorted_probs = probs[order]
    sorted_y = y[order]

    n = len(sorted_probs)
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_probs[j + 1] == sorted_probs[i]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        ranks[i : j + 1] = avg_rank
        i = j + 1

    rank_sum_pos = float(ranks[sorted_y == 1].sum())
    auc = (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def _score_metric(probabilities: np.ndarray, labels: np.ndarray, metric: str) -> float:
    probs = np.asarray(probabilities, dtype=np.float64).reshape(-1)
    targets = _normalize_labels(labels)

    metric_key = metric.strip().lower()
    if metric_key in {"auc", "roc_auc"}:
        return _binary_roc_auc(probs, targets)

    preds = (probs >= 0.5).astype(np.int32)
    if metric_key in {"accuracy", "acc"}:
        return float((preds == targets).mean())

    if metric_key == "mcc":
        tp = int(((preds == 1) & (targets == 1)).sum())
        tn = int(((preds == 0) & (targets == 0)).sum())
        fp = int(((preds == 1) & (targets == 0)).sum())
        fn = int(((preds == 0) & (targets == 1)).sum())
        denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        return float(((tp * tn) - (fp * fn)) / math.sqrt(denom)) if denom > 0 else 0.0

    raise ValueError(f"Unsupported metric: {metric}")


def _estimate_grid_candidates(grid_units: int, num_models: int) -> int:
    """Estimate grid search candidate count for simplex partitioning."""
    if grid_units < 0 or num_models < 1:
        return 0
    if grid_units == 0:
        return 1
    return math.comb(grid_units + num_models - 1, num_models - 1)


def _iter_weight_compositions(units: int, slots: int) -> Iterable[Tuple[int, ...]]:
    """Yield integer compositions of `units` into `slots` non-negative parts."""
    if slots == 1:
        yield (units,)
        return
    for i in range(units + 1):
        for remainder in _iter_weight_compositions(units - i, slots - 1):
            yield (i,) + remainder


def optimize_ensemble_weights(
    val_probs: Dict,
    val_labels: np.ndarray,
    metric: str = "accuracy",
    min_weight: float = 0.01,
    grid_step: float = 0.02,
    max_candidates: int = 50000,
    random_seed: int = 13,
    return_metadata: bool = False,
) -> Dict | Tuple[Dict, Dict]:
    """
    Find optimal ensemble weights for a validation set.

    Args:
        val_probs: Mapping of k-mer to probability array.
        val_labels: Ground-truth labels (0/1).
        metric: Metric to maximize ("accuracy", "mcc", or "auc").
        min_weight: Minimum weight per model.
        grid_step: Approximate grid step for simplex search.
        max_candidates: Max number of weight candidates to evaluate.
        random_seed: RNG seed for random search fallback.
        return_metadata: When True, return (weights, metadata).

    Returns:
        Dict of k-mer to optimized weight (and metadata if requested).
    """
    if not val_probs:
        raise ValueError("val_probs must be a non-empty mapping")

    kmers = _sorted_kmers(val_probs.keys())
    probs_matrix = np.stack(
        [np.asarray(val_probs[k], dtype=np.float64).reshape(-1) for k in kmers],
        axis=1,
    )
    labels = np.asarray(val_labels, dtype=np.float64).reshape(-1)
    if probs_matrix.shape[0] != labels.shape[0]:
        raise ValueError("val_probs and val_labels must have the same number of samples")

    num_models = len(kmers)
    if num_models == 1:
        weights = {kmers[0]: 1.0}
        metadata = {
            "metric": metric,
            "min_weight": min_weight,
            "search_method": "single",
            "candidates_evaluated": 1,
            "best_score": _score_metric(probs_matrix[:, 0], labels, metric),
        }
        return (weights, metadata) if return_metadata else weights

    if min_weight < 0.0:
        raise ValueError("min_weight must be non-negative")

    remaining = 1.0 - (min_weight * num_models)
    if remaining < 0.0:
        raise ValueError("min_weight is too large for the number of models")

    grid_units = 0
    grid_step = float(grid_step)
    if grid_step > 0.0 and remaining > 0.0:
        grid_units = max(1, int(round(remaining / grid_step)))
    candidate_estimate = _estimate_grid_candidates(grid_units, num_models)

    def evaluate_weights(weights_vec: np.ndarray) -> float:
        ensemble_probs = probs_matrix @ weights_vec
        return _score_metric(ensemble_probs, labels, metric)

    best_score = -float("inf")
    best_weights = None
    candidates_evaluated = 0
    search_method = "grid"

    # Always consider equal weights + single-model heavy weights.
    equal_weights = np.full(num_models, 1.0 / num_models, dtype=np.float64)
    best_score = evaluate_weights(equal_weights)
    best_weights = equal_weights
    candidates_evaluated += 1

    for idx in range(num_models):
        weights = np.full(num_models, min_weight, dtype=np.float64)
        weights[idx] += remaining
        score = evaluate_weights(weights)
        candidates_evaluated += 1
        if score > best_score:
            best_score = score
            best_weights = weights

    if candidate_estimate <= max_candidates and remaining > 0.0:
        step = remaining / grid_units if grid_units > 0 else remaining
        for units in _iter_weight_compositions(grid_units, num_models):
            weights = np.array(units, dtype=np.float64) * step + min_weight
            score = evaluate_weights(weights)
            candidates_evaluated += 1
            if score > best_score:
                best_score = score
                best_weights = weights
    else:
        search_method = "random"
        rng = np.random.default_rng(random_seed)
        num_samples = max(1, int(max_candidates))
        for _ in range(num_samples):
            draw = rng.dirichlet(np.ones(num_models))
            weights = min_weight + remaining * draw
            score = evaluate_weights(weights)
            candidates_evaluated += 1
            if score > best_score:
                best_score = score
                best_weights = weights

    if best_weights is None:
        raise RuntimeError("Failed to optimize ensemble weights")

    weight_sum = float(best_weights.sum())
    if weight_sum <= 0.0:
        raise RuntimeError("Optimized weights sum to zero")
    best_weights = best_weights / weight_sum

    weights_dict = {k: float(w) for k, w in zip(kmers, best_weights)}
    metadata = {
        "metric": metric,
        "min_weight": min_weight,
        "grid_step": grid_step,
        "search_method": search_method,
        "candidates_evaluated": int(candidates_evaluated),
        "best_score": float(best_score),
    }
    return (weights_dict, metadata) if return_metadata else weights_dict


def save_ensemble_weights(
    weights_path: Path,
    weights: Dict,
    metadata: Dict | None = None,
) -> None:
    """Persist optimized ensemble weights to JSON."""
    payload = {"weights": {str(k): float(v) for k, v in weights.items()}}
    if metadata:
        payload["metadata"] = metadata

    weights_path.parent.mkdir(parents=True, exist_ok=True)
    with weights_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def load_ensemble_weights(weights_path: Path) -> Tuple[Dict[int, float], Dict]:
    """Load ensemble weights from JSON."""
    with weights_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if isinstance(payload, dict) and "weights" in payload:
        weights_payload = payload["weights"]
        metadata = payload.get("metadata", {})
    else:
        weights_payload = payload
        metadata = {}

    if not isinstance(weights_payload, dict):
        raise ValueError("weights payload must be a mapping")

    weights: Dict[int, float] = {}
    for key, value in weights_payload.items():
        try:
            parsed_key = int(key)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid k-mer key in weights file: {key!r}") from None
        weights[parsed_key] = float(value)

    return weights, metadata


def align_weights_to_kmers(weights: Dict[int, float], kmers: List[int]) -> List[float] | None:
    """Align weight dict to an ordered list of kmers, returning normalized weights."""
    missing = [k for k in kmers if k not in weights]
    if missing:
        LOGGER.warning("Weights missing for k-mers: %s", ",".join(str(k) for k in missing))
        return None

    values = [float(weights[k]) for k in kmers]
    weight_sum = float(sum(values))
    if weight_sum <= 0.0:
        LOGGER.warning("Weight sum is non-positive: %s", weight_sum)
        return None

    return [value / weight_sum for value in values]
