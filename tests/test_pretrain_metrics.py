from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.training.pretrain_metrics import (
    EpochMetricsWriter,
    build_epoch_metrics,
)


def test_build_epoch_metrics_computes_token_rates() -> None:
    metrics = build_epoch_metrics(
        epoch=3,
        global_step=12,
        avg_loss=2.5,
        best_loss=2.4,
        best_epoch=2,
        learning_rate=1e-4,
        elapsed_seconds=4.0,
        trained_batches=4,
        skipped_batches=1,
        supervised_tokens=40,
        tokenized_lengths=[10, 12, 8, 10],
        maskable_counts=[6, 7, 5, 6],
        masked_counts=[2, 2, 2, 2],
        batch_times=[0.3, 0.4, 0.5, 0.6, 0.2],
        grad_norms=[1.1, 0.9],
    )

    assert metrics["skipped_batch_rate"] == 0.2
    assert metrics["maskable_token_rate"] == 24 / 40
    assert metrics["masked_token_rate"] == 8 / 40
    assert metrics["supervised_tokens_per_batch"] == 10.0
    assert metrics["tokenized_length_p50"] == 10.0


def test_epoch_metrics_writer_creates_manifest_jsonl_and_csv(tmp_path: Path) -> None:
    writer = EpochMetricsWriter(tmp_path)
    writer.write_run_manifest({"config": {"seed": 7}})
    writer.append_epoch(
        build_epoch_metrics(
            epoch=1,
            global_step=2,
            avg_loss=1.5,
            best_loss=1.5,
            best_epoch=1,
            learning_rate=2e-4,
            elapsed_seconds=2.0,
            trained_batches=2,
            skipped_batches=0,
            supervised_tokens=20,
            tokenized_lengths=[10, 10],
            maskable_counts=[7, 7],
            masked_counts=[2, 2],
            batch_times=[0.2, 0.2],
            grad_norms=[0.8],
        )
    )

    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    metrics_lines = (tmp_path / "metrics.jsonl").read_text().strip().splitlines()
    csv_lines = (tmp_path / "epoch_metrics.csv").read_text().strip().splitlines()

    assert manifest["config"]["seed"] == 7
    assert len(metrics_lines) == 1
    assert len(csv_lines) == 2
