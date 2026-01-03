#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional


def _run_command(cmd: list[str], *, cwd: Path, label: str) -> int:
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"{label} failed with exit code {result.returncode}.", file=sys.stderr)
    return result.returncode


def _parse_metrics_line(metrics: Dict[str, float], line: str) -> None:
    patterns = {
        "accuracy": r"Accuracy@0\.5:\s+([0-9.]+)",
        "f1": r"F1@0\.5:\s+([0-9.]+)",
        "mcc": r"MCC@0\.5:\s+([0-9.]+)",
        "best_thr": r"Best thr:\s+([0-9.]+)",
        "best_acc": r"Best acc:\s+([0-9.]+)",
    }
    for key, pattern in patterns.items():
        if key in metrics:
            continue
        match = re.search(pattern, line)
        if match:
            metrics[key] = float(match.group(1))


def _run_eval_and_capture(
    cmd: list[str],
    *,
    cwd: Path,
) -> tuple[int, Dict[str, float]]:
    metrics: Dict[str, float] = {}
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="")
        _parse_metrics_line(metrics, line)
    return_code = process.wait()
    if return_code != 0:
        print(f"Stage1 multiscale evaluation failed with exit code {return_code}.", file=sys.stderr)
    return return_code, metrics


def _format_metrics(metrics: Dict[str, float]) -> Optional[str]:
    required = ["accuracy", "f1", "mcc", "best_thr", "best_acc"]
    if not all(key in metrics for key in required):
        return None
    return (
        "Final ensemble metrics\n"
        f"Accuracy@0.5: {metrics['accuracy']:.4f}\n"
        f"F1@0.5:       {metrics['f1']:.4f}\n"
        f"MCC@0.5:      {metrics['mcc']:.4f}\n"
        f"Best thr:     {metrics['best_thr']:.2f}\n"
        f"Best acc:     {metrics['best_acc']:.4f}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Stage1 multiscale training + eval pipeline.")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (relative to repo root or training/).",
    )
    parser.add_argument(
        "--skip_train",
        action="store_true",
        help="Skip training and only run the Stage1 multiscale evaluation.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]

    if not args.skip_train:
        train_cmd = [
            sys.executable,
            "scripts/train_stage1_multiscale.py",
            "--config",
            args.config,
        ]
        train_code = _run_command(train_cmd, cwd=repo_root, label="Stage1 multiscale training")
        if train_code != 0:
            return train_code

    eval_cmd = [
        sys.executable,
        "scripts/eval_stage1_multiscale.py",
        "--config",
        args.config,
    ]
    eval_code, metrics = _run_eval_and_capture(eval_cmd, cwd=repo_root)
    if eval_code != 0:
        return eval_code

    formatted = _format_metrics(metrics)
    if formatted is None:
        print("Warning: unable to parse final ensemble metrics from eval output.", file=sys.stderr)
    else:
        print("\n" + formatted)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
