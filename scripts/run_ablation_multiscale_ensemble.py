#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from ablation_variants import STAGE1_VARIANTS  # noqa: E402


@dataclass
class VariantResult:
    name: str
    kmers: list[int]
    accuracy: float
    f1: float
    mcc: float
    best_thr: float
    best_acc: float
    config_path: Path
    output_dir: Path


def _parse_variants(value: str) -> list[str]:
    parts = [chunk.strip() for chunk in value.replace(",", " ").split()]
    parts = [chunk for chunk in parts if chunk]
    if not parts:
        raise argparse.ArgumentTypeError("Expected at least one variant name")
    return parts


def _parse_kmers_arg(value: str) -> list[int]:
    tokens = [item.strip() for item in value.split(",")]
    kmers: list[int] = []
    for token in tokens:
        if not token:
            continue
        try:
            kmers.append(int(token))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid k-mer value: {token!r}") from exc
    if not kmers:
        raise argparse.ArgumentTypeError("Expected at least one k-mer value in --kmers.")
    if len(set(kmers)) != len(kmers):
        raise argparse.ArgumentTypeError(f"Duplicate k-mer values in --kmers: {kmers}")
    return kmers


def _resolve_config_path(config_arg: str, *, training_dir: Path, initial_cwd: Path) -> Path:
    candidate = Path(config_arg)
    if candidate.is_absolute() and candidate.exists():
        return candidate

    by_cwd = (initial_cwd / candidate).resolve()
    if by_cwd.exists():
        return by_cwd

    by_training = (training_dir / candidate).resolve()
    if by_training.exists():
        return by_training

    raise FileNotFoundError(
        f"Config not found: tried {by_cwd} and {by_training} (arg={config_arg!r})"
    )


def _resolve_path_like(value: Any, *, base_dir: Path) -> Any:
    if value is None:
        return None
    if isinstance(value, Path):
        candidate = value
    elif isinstance(value, str):
        text = value.strip()
        if not text:
            return value
        candidate = Path(text)
    else:
        return value
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    return str(candidate)


def _materialize_config_paths(cfg: dict[str, Any], *, base_dir: Path) -> dict[str, Any]:
    export_cfg = dict(cfg)
    for key in (
        "stage1_train",
        "stage1_val",
        "stage2_train",
        "stage2_val",
        "stage2_test",
        "vocab_cache_dir",
        "mlm_encoder_checkpoint",
        "mlm_encoder_path",
        "mlm_fasta_path",
        "mlm_vocab_path",
    ):
        if key in export_cfg:
            export_cfg[key] = _resolve_path_like(export_cfg.get(key), base_dir=base_dir)
    return export_cfg


def _resolve_stage1_ckpt_map(
    cfg: dict[str, Any],
    *,
    base_dir: Path,
    artifacts_root: Path,
    kmers: Iterable[int],
) -> dict[str, str]:
    stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k")
    if not isinstance(stage1_ckpt_by_k, dict):
        stage1_ckpt_by_k = {}

    resolved: dict[str, str] = {}
    for k in kmers:
        raw = stage1_ckpt_by_k.get(k) or stage1_ckpt_by_k.get(str(k))
        if not raw:
            raw = artifacts_root / "checkpoints" / f"stage1_k{k}.pt"
        resolved[str(k)] = str(_resolve_path_like(raw, base_dir=base_dir))
    return resolved


def _slugify(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return cleaned or "variant"


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


def _run_command(cmd: list[str], *, cwd: Path, label: str) -> int:
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"{label} failed with exit code {result.returncode}.", file=sys.stderr)
    return result.returncode


def _run_eval_and_capture(cmd: list[str], *, cwd: Path) -> tuple[int, Dict[str, float], list[str]]:
    metrics: Dict[str, float] = {}
    captured: list[str] = []
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
        captured.append(line)
        _parse_metrics_line(metrics, line)
    return_code = process.wait()
    if return_code != 0:
        print(f"Stage1 multiscale evaluation failed with exit code {return_code}.", file=sys.stderr)
    return return_code, metrics, captured


def _write_summary_csv(path: Path, results: list[VariantResult]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "variant",
                "kmers",
                "accuracy",
                "f1",
                "mcc",
                "best_thr",
                "best_acc",
                "config_path",
                "output_dir",
            ],
        )
        writer.writeheader()
        for row in results:
            writer.writerow(
                {
                    "variant": row.name,
                    "kmers": ",".join(str(k) for k in row.kmers),
                    "accuracy": f"{row.accuracy:.6f}",
                    "f1": f"{row.f1:.6f}",
                    "mcc": f"{row.mcc:.6f}",
                    "best_thr": f"{row.best_thr:.4f}",
                    "best_acc": f"{row.best_acc:.6f}",
                    "config_path": str(row.config_path),
                    "output_dir": str(row.output_dir),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run multiscale Stage1 ensemble evaluation for each ablation variant."
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (relative to repo root or training/).",
    )
    parser.add_argument(
        "--variants",
        type=_parse_variants,
        default=None,
        help="Comma-separated list of variants to run (default: all).",
    )
    parser.add_argument(
        "--kmers",
        type=_parse_kmers_arg,
        default=None,
        help="Comma-separated list of k-mers to evaluate (overrides config).",
    )
    parser.add_argument(
        "--skip_train",
        action="store_true",
        help="Skip training and only run ensemble evaluation.",
    )
    parser.add_argument(
        "--eval_batch_size",
        type=int,
        default=None,
        help="Override batch size for evaluation.",
    )
    parser.add_argument(
        "--debug_first_batch",
        action="store_true",
        help="Check ensemble output matches mean of per-k outputs on the first batch.",
    )
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--weighted",
        action="store_true",
        help="Use weighted ensemble (default).",
    )
    mode_group.add_argument(
        "--unweighted",
        action="store_true",
        help="Use unweighted ensemble.",
    )
    args = parser.parse_args()

    training_dir = REPO_ROOT / "gene_whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    initial_cwd = Path.cwd()
    config_path = _resolve_config_path(args.config, training_dir=training_dir, initial_cwd=initial_cwd)
    with config_path.open("r", encoding="utf-8") as handle:
        base_cfg = yaml.safe_load(handle) or {}
    if not isinstance(base_cfg, dict):
        raise SystemExit(f"Expected config to be a mapping, got {type(base_cfg).__name__}")

    variants = [("baseline", {})] + STAGE1_VARIANTS
    variants_by_name = {name: overrides for name, overrides in variants}
    if args.variants is None:
        selected = variants
    else:
        requested: list[tuple[str, dict[str, Any]]] = []
        seen = set()
        for name in args.variants:
            if name not in variants_by_name:
                raise SystemExit(f"Unknown variant: {name}")
            if name not in seen:
                requested.append((name, variants_by_name[name]))
                seen.add(name)
        if "baseline" not in seen:
            print("Note: adding baseline to variants for comparison.")
            requested.insert(0, ("baseline", {}))
        selected = requested

    artifacts_root = (training_dir / "../artifacts").resolve()
    output_root = artifacts_root / "ablation_multiscale"
    configs_dir = output_root / "configs"
    output_root.mkdir(parents=True, exist_ok=True)
    configs_dir.mkdir(parents=True, exist_ok=True)

    results: list[VariantResult] = []
    for name, overrides in selected:
        variant_slug = _slugify(name)
        variant_dir = output_root / variant_slug
        variant_dir.mkdir(parents=True, exist_ok=True)

        cfg_run = copy.deepcopy(base_cfg)
        cfg_run.update(overrides)
        cfg_run["use_multi_scale"] = True
        if args.kmers is not None:
            cfg_run["multi_scale_kmers"] = args.kmers
        kmers = cfg_run.get("multi_scale_kmers")
        if not isinstance(kmers, list) or not kmers:
            raise SystemExit("multi_scale_kmers must be set in config or via --kmers")

        cfg_run = _materialize_config_paths(cfg_run, base_dir=config_path.parent)
        cfg_run["stage1_ckpt_by_k"] = _resolve_stage1_ckpt_map(
            cfg_run,
            base_dir=config_path.parent,
            artifacts_root=artifacts_root,
            kmers=kmers,
        )

        config_out = configs_dir / f"{variant_slug}.yaml"
        with config_out.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(cfg_run, handle, sort_keys=False, allow_unicode=False)

        if not args.skip_train:
            train_cmd = [
                sys.executable,
                "scripts/train_stage1_multiscale.py",
                "--config",
                str(config_out),
            ]
            train_code = _run_command(train_cmd, cwd=REPO_ROOT, label=f"Training ({name})")
            if train_code != 0:
                return train_code

        eval_cmd = [
            sys.executable,
            "scripts/eval_stage1_multiscale.py",
            "--config",
            str(config_out),
        ]
        if args.kmers is not None:
            eval_cmd.extend(["--kmers", ",".join(str(k) for k in args.kmers)])
        if args.eval_batch_size is not None:
            eval_cmd.extend(["--batch_size", str(args.eval_batch_size)])
        if args.debug_first_batch:
            eval_cmd.append("--debug_first_batch")
        if args.unweighted:
            eval_cmd.append("--unweighted")
        elif args.weighted:
            eval_cmd.append("--weighted")

        eval_code, metrics, output_lines = _run_eval_and_capture(eval_cmd, cwd=REPO_ROOT)
        if eval_code != 0:
            return eval_code

        eval_log_path = variant_dir / "eval_output.txt"
        with eval_log_path.open("w", encoding="utf-8") as handle:
            handle.writelines(output_lines)

        for k in kmers:
            report_path = artifacts_root / f"stage1_report_k{k}.json"
            if report_path.exists():
                shutil.copy2(report_path, variant_dir / report_path.name)
            checkpoint_path = artifacts_root / "checkpoints" / f"stage1_k{k}.pt"
            if checkpoint_path.exists():
                shutil.copy2(checkpoint_path, variant_dir / checkpoint_path.name)

        results.append(
            VariantResult(
                name=name,
                kmers=list(kmers),
                accuracy=float(metrics.get("accuracy", float("nan"))),
                f1=float(metrics.get("f1", float("nan"))),
                mcc=float(metrics.get("mcc", float("nan"))),
                best_thr=float(metrics.get("best_thr", float("nan"))),
                best_acc=float(metrics.get("best_acc", float("nan"))),
                config_path=config_out,
                output_dir=variant_dir,
            )
        )

        print(f"Completed {name} at {time.strftime('%H:%M:%S')}")

    summary_path = output_root / "ablation_multiscale_summary.csv"
    _write_summary_csv(summary_path, results)
    print(f"Wrote summary CSV: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
