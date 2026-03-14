"""Structured metrics and manifests for MLM pretraining."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def summarize_series(values: Sequence[float] | Sequence[int], prefix: str) -> dict[str, float]:
    if not values:
        return {
            f"{prefix}_min": 0.0,
            f"{prefix}_mean": 0.0,
            f"{prefix}_p50": 0.0,
            f"{prefix}_p90": 0.0,
            f"{prefix}_max": 0.0,
        }

    array = np.asarray(values, dtype=np.float64)
    return {
        f"{prefix}_min": float(array.min()),
        f"{prefix}_mean": float(array.mean()),
        f"{prefix}_p50": float(np.percentile(array, 50)),
        f"{prefix}_p90": float(np.percentile(array, 90)),
        f"{prefix}_max": float(array.max()),
    }


def build_run_manifest(
    config: Any,
    tokenizer_vocab_size: int,
    tokenizer_metadata: Mapping[str, Any],
    sequence_count: int,
    total_bases: int,
    batches_per_epoch: int,
    source_summary: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "sequence_count": sequence_count,
        "total_bases": total_bases,
        "batches_per_epoch": batches_per_epoch,
        "tokenizer_vocab_size": tokenizer_vocab_size,
        "tokenizer_metadata": _json_safe(dict(tokenizer_metadata)),
        "source_summary": _json_safe(dict(source_summary)),
        "config": _json_safe(vars(config)),
    }


def build_epoch_metrics(
    *,
    epoch: int,
    global_step: int,
    avg_loss: float,
    best_loss: float,
    best_epoch: int,
    learning_rate: float,
    elapsed_seconds: float,
    trained_batches: int,
    skipped_batches: int,
    supervised_tokens: int,
    tokenized_lengths: Sequence[int],
    maskable_counts: Sequence[int],
    masked_counts: Sequence[int],
    batch_times: Sequence[float],
    grad_norms: Sequence[float],
) -> dict[str, Any]:
    total_batches = trained_batches + skipped_batches
    total_tokenized_tokens = int(sum(tokenized_lengths))
    total_maskable_tokens = int(sum(maskable_counts))
    total_masked_tokens = int(sum(masked_counts))
    skipped_batch_rate = (skipped_batches / total_batches) if total_batches else 0.0
    maskable_token_rate = (
        total_maskable_tokens / total_tokenized_tokens
        if total_tokenized_tokens
        else 0.0
    )
    masked_token_rate = (
        total_masked_tokens / total_tokenized_tokens
        if total_tokenized_tokens
        else 0.0
    )

    metrics: dict[str, Any] = {
        "epoch": epoch,
        "global_step": global_step,
        "avg_loss_token_weighted": avg_loss,
        "best_loss": best_loss,
        "best_epoch": best_epoch,
        "learning_rate": learning_rate,
        "elapsed_seconds": elapsed_seconds,
        "trained_batches": trained_batches,
        "skipped_batches": skipped_batches,
        "skipped_batch_rate": skipped_batch_rate,
        "supervised_tokens": supervised_tokens,
        "supervised_tokens_per_batch": (
            supervised_tokens / trained_batches if trained_batches else 0.0
        ),
        "supervised_tokens_per_second": (
            supervised_tokens / elapsed_seconds if elapsed_seconds > 0 else 0.0
        ),
        "tokenized_tokens_total": total_tokenized_tokens,
        "maskable_tokens_total": total_maskable_tokens,
        "masked_tokens_total": total_masked_tokens,
        "maskable_token_rate": maskable_token_rate,
        "masked_token_rate": masked_token_rate,
    }
    metrics.update(summarize_series(tokenized_lengths, "tokenized_length"))
    metrics.update(summarize_series(maskable_counts, "maskable_count"))
    metrics.update(summarize_series(masked_counts, "masked_count"))
    metrics.update(summarize_series(batch_times, "batch_time"))
    metrics.update(summarize_series(grad_norms, "grad_norm"))
    return metrics


class EpochMetricsWriter:
    """Persist run artifacts in JSONL, CSV, and manifest form."""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.metrics_jsonl_path = output_dir / "metrics.jsonl"
        self.epoch_csv_path = output_dir / "epoch_metrics.csv"
        self.run_manifest_path = output_dir / "run_manifest.json"
        self.metrics_jsonl_path.unlink(missing_ok=True)
        self.epoch_csv_path.unlink(missing_ok=True)
        self._csv_fieldnames: list[str] | None = None

    def write_run_manifest(self, manifest: Mapping[str, Any]) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with self.run_manifest_path.open("w", encoding="utf-8") as handle:
            json.dump(_json_safe(dict(manifest)), handle, indent=2, sort_keys=True)

    def append_epoch(self, metrics: Mapping[str, Any]) -> None:
        safe_metrics = _json_safe(dict(metrics))
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with self.metrics_jsonl_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(safe_metrics, sort_keys=True))
            handle.write("\n")

        if self._csv_fieldnames is None:
            self._csv_fieldnames = list(safe_metrics.keys())
            with self.epoch_csv_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=self._csv_fieldnames)
                writer.writeheader()
                writer.writerow(safe_metrics)
            return

        with self.epoch_csv_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self._csv_fieldnames)
            writer.writerow(safe_metrics)
