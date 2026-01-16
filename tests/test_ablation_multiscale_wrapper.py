from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import run_ablation_multiscale_ensemble as wrapper  # noqa: E402


def test_slugify_handles_symbols_and_empty() -> None:
    assert wrapper._slugify("no-cksnap") == "no_cksnap"
    assert wrapper._slugify("  ") == "variant"


def test_parse_metrics_line_extracts_values() -> None:
    metrics: dict[str, float] = {}
    wrapper._parse_metrics_line(metrics, "Accuracy@0.5: 0.8123 (81.23%)")
    wrapper._parse_metrics_line(metrics, "F1@0.5:       0.7012")
    wrapper._parse_metrics_line(metrics, "MCC@0.5:      0.6221")
    wrapper._parse_metrics_line(metrics, "Best thr:     0.42")
    wrapper._parse_metrics_line(metrics, "Best acc:     0.8555 (85.55%)")

    assert metrics["accuracy"] == 0.8123
    assert metrics["f1"] == 0.7012
    assert metrics["mcc"] == 0.6221
    assert metrics["best_thr"] == 0.42
    assert metrics["best_acc"] == 0.8555
