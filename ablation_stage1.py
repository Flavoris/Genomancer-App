from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml


def _parse_seeds(value: str) -> list[int]:
    parts: list[int] = []
    for chunk in value.replace(",", " ").split():
        if not chunk:
            continue
        try:
            parts.append(int(chunk))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid seed: {chunk!r}") from exc
    if not parts:
        raise argparse.ArgumentTypeError("Expected at least one seed")
    return parts


def _mean_std(values: Iterable[float]) -> tuple[float, float]:
    arr = np.asarray(list(values), dtype=np.float64)
    if arr.size == 0:
        return float("nan"), float("nan")
    mean = float(arr.mean())
    if arr.size == 1:
        return mean, 0.0
    return mean, float(arr.std(ddof=1))


@dataclass(frozen=True)
class VariantSummary:
    name: str
    seeds: list[int]
    auc_mean: float
    auc_std: float
    acc_mean: float
    acc_std: float
    mcc_mean: float
    mcc_std: float
    delta_auc_mean: float
    delta_acc_mean: float
    delta_mcc_mean: float


def _resolve_config_path(config_arg: str, *, training_dir: Path, initial_cwd: Path) -> Path:
    candidate = Path(config_arg)
    if candidate.is_absolute() and candidate.exists():
        return candidate

    # Prefer paths relative to the invocation directory (usually repo root).
    by_cwd = (initial_cwd / candidate).resolve()
    if by_cwd.exists():
        return by_cwd

    # Then fallback to the training directory (matches train_stage1.py behavior).
    by_training = (training_dir / candidate).resolve()
    if by_training.exists():
        return by_training

    raise FileNotFoundError(
        f"Config not found: tried {by_cwd} and {by_training} (arg={config_arg!r})"
    )


def _apply_quick_defaults(cfg: dict[str, Any]) -> None:
    if cfg.get("ablation_max_train_samples_stage1") is None:
        cfg["ablation_max_train_samples_stage1"] = 2048
    if cfg.get("ablation_max_val_samples_stage1") is None:
        cfg["ablation_max_val_samples_stage1"] = 2048
    if cfg.get("ablation_max_steps_stage1") is None:
        cfg["ablation_max_steps_stage1"] = 200
    if cfg.get("epochs") is None:
        cfg["epochs"] = 5


def _run_variant(
    *,
    base_cfg: dict[str, Any],
    variant_name: str,
    variant_overrides: dict[str, Any],
    config_path: Path,
    seeds: list[int],
    run_stage1_training: Any,
) -> dict[int, dict[str, float]]:
    results_by_seed: dict[int, dict[str, float]] = {}

    for seed in seeds:
        overrides: dict[str, Any] = {
            "_config_path": str(config_path),
            "seed": seed,
        }
        overrides.update(variant_overrides)

        report = run_stage1_training(base_cfg, overrides=overrides)
        try:
            results_by_seed[seed] = {
                "best_val_auc": float(report["best_val_auc"]),
                "best_val_acc": float(report["best_val_acc"]),
                "best_val_mcc": float(report["best_val_mcc"]),
            }
        except Exception as exc:
            raise RuntimeError(
                f"Unexpected training report for {variant_name} seed={seed}: {report}"
            ) from exc

    return results_by_seed


def _summarize_variant(
    *,
    name: str,
    seeds: list[int],
    results_by_seed: dict[int, dict[str, float]],
) -> tuple[float, float, float, float, float, float]:
    auc_mean, auc_std = _mean_std(results_by_seed[s]["best_val_auc"] for s in seeds)
    acc_mean, acc_std = _mean_std(results_by_seed[s]["best_val_acc"] for s in seeds)
    mcc_mean, mcc_std = _mean_std(results_by_seed[s]["best_val_mcc"] for s in seeds)
    return auc_mean, auc_std, acc_mean, acc_std, mcc_mean, mcc_std


def _format_pm(mean: float, std: float) -> str:
    return f"{mean:.4f}±{std:.4f}"


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


def _materialize_config_for_export(cfg: dict[str, Any], *, base_dir: Path) -> dict[str, Any]:
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


def _read_negligible_candidates_from_csv(
    path: Path,
    *,
    negligible_auc: float,
    negligible_mcc: float,
) -> list[str]:
    candidates: list[str] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = str(row.get("name") or "").strip()
            if not name or name == "baseline":
                continue
            try:
                delta_auc = float(row.get("delta_auc") or "nan")
                delta_mcc = float(row.get("delta_mcc") or "nan")
            except ValueError:
                continue
            if delta_auc < negligible_auc and delta_mcc < negligible_mcc:
                candidates.append(name)
    return candidates


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 1 ablation runner")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument(
        "--seeds",
        type=_parse_seeds,
        default=_parse_seeds("1337"),
        help="Comma-separated list of seeds (e.g. 1337 or 1337,2027,9001)",
    )
    args = parser.parse_args()

    initial_cwd = Path.cwd()
    repo_root = Path(__file__).resolve().parent
    training_dir = repo_root / "Gene Whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    # Ensure training modules import cleanly (directory name contains a space).
    sys.path.insert(0, str(training_dir))
    from train_stage1 import load_config, run_stage1_training  # type: ignore

    config_path = _resolve_config_path(args.config, training_dir=training_dir, initial_cwd=initial_cwd)
    cfg: dict[str, Any] = load_config(config_path)
    _apply_quick_defaults(cfg)

    # Resolve relative dataset paths consistently with the config file location.
    os.chdir(config_path.parent)

    seeds: list[int] = args.seeds

    ablations: list[tuple[str, dict[str, Any]]] = [
        ("no_attention_pool", {"use_attention_pool": False}),
        ("no_tcn", {"use_tcn": False}),
        ("no_postcnn_transformer", {"post_cnn_transformer_layers": 0}),
        ("no_engineered_features", {"stage1_use_engineered_features": False}),
        ("no_tnc", {"stage1_feature_enable_tnc": False}),
        ("no_pseeiip", {"stage1_feature_enable_pseeiip": False}),
    ]

    baseline_results = _run_variant(
        base_cfg=cfg,
        variant_name="baseline",
        variant_overrides={},
        config_path=config_path,
        seeds=seeds,
        run_stage1_training=run_stage1_training,
    )
    base_auc_mean, base_auc_std, base_acc_mean, base_acc_std, base_mcc_mean, base_mcc_std = _summarize_variant(
        name="baseline",
        seeds=seeds,
        results_by_seed=baseline_results,
    )

    summaries: list[VariantSummary] = [
        VariantSummary(
            name="baseline",
            seeds=seeds,
            auc_mean=base_auc_mean,
            auc_std=base_auc_std,
            acc_mean=base_acc_mean,
            acc_std=base_acc_std,
            mcc_mean=base_mcc_mean,
            mcc_std=base_mcc_std,
            delta_auc_mean=0.0,
            delta_acc_mean=0.0,
            delta_mcc_mean=0.0,
        )
    ]

    for name, overrides in ablations:
        variant_results = _run_variant(
            base_cfg=cfg,
            variant_name=name,
            variant_overrides=overrides,
            config_path=config_path,
            seeds=seeds,
            run_stage1_training=run_stage1_training,
        )
        auc_mean, auc_std, acc_mean, acc_std, mcc_mean, mcc_std = _summarize_variant(
            name=name,
            seeds=seeds,
            results_by_seed=variant_results,
        )

        delta_auc_mean, _ = _mean_std(
            baseline_results[s]["best_val_auc"] - variant_results[s]["best_val_auc"] for s in seeds
        )
        delta_acc_mean, _ = _mean_std(
            baseline_results[s]["best_val_acc"] - variant_results[s]["best_val_acc"] for s in seeds
        )
        delta_mcc_mean, _ = _mean_std(
            baseline_results[s]["best_val_mcc"] - variant_results[s]["best_val_mcc"] for s in seeds
        )

        summaries.append(
            VariantSummary(
                name=name,
                seeds=seeds,
                auc_mean=auc_mean,
                auc_std=auc_std,
                acc_mean=acc_mean,
                acc_std=acc_std,
                mcc_mean=mcc_mean,
                mcc_std=mcc_std,
                delta_auc_mean=float(delta_auc_mean),
                delta_acc_mean=float(delta_acc_mean),
                delta_mcc_mean=float(delta_mcc_mean),
            )
        )

    baseline, *ablation_rows = summaries
    ablation_rows.sort(key=lambda r: r.delta_auc_mean, reverse=True)
    ordered = [baseline] + ablation_rows

    artifacts_root = (training_dir / "../artifacts").resolve()
    artifacts_root.mkdir(parents=True, exist_ok=True)
    out_path = artifacts_root / "ablation_stage1.csv"
    negligible_auc = 0.005
    negligible_mcc = 0.005

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "seeds",
                "n_seeds",
                "auc_mean",
                "auc_std",
                "acc_mean",
                "acc_std",
                "mcc_mean",
                "mcc_std",
                "delta_auc",
                "delta_acc",
                "delta_mcc",
            ],
        )
        writer.writeheader()
        for row in ordered:
            writer.writerow(
                {
                    "name": row.name,
                    "seeds": ",".join(str(s) for s in row.seeds),
                    "n_seeds": len(row.seeds),
                    "auc_mean": f"{row.auc_mean:.6f}",
                    "auc_std": f"{row.auc_std:.6f}",
                    "acc_mean": f"{row.acc_mean:.6f}",
                    "acc_std": f"{row.acc_std:.6f}",
                    "mcc_mean": f"{row.mcc_mean:.6f}",
                    "mcc_std": f"{row.mcc_std:.6f}",
                    "delta_auc": f"{row.delta_auc_mean:.6f}",
                    "delta_acc": f"{row.delta_acc_mean:.6f}",
                    "delta_mcc": f"{row.delta_mcc_mean:.6f}",
                }
            )

    ablation_by_name = dict(ablations)
    negligible_rows = [
        row
        for row in ablation_rows
        if row.delta_auc_mean < negligible_auc and row.delta_mcc_mean < negligible_mcc
    ]
    candidates_to_remove = [row.name for row in negligible_rows]

    stage2_csv_path = artifacts_root / "ablation_stage2.csv"
    stage2_candidates = (
        _read_negligible_candidates_from_csv(
            stage2_csv_path,
            negligible_auc=negligible_auc,
            negligible_mcc=negligible_mcc,
        )
        if stage2_csv_path.exists()
        else None
    )
    overall_candidates = (
        sorted(set(candidates_to_remove).intersection(stage2_candidates)) if stage2_candidates is not None else None
    )
    if overall_candidates is not None:
        stage2_shared_overrides: dict[str, dict[str, Any]] = {
            "no_engineered_features": {"stage2_use_engineered_features": False},
        }
        overall_overrides: dict[str, Any] = {}
        for name in overall_candidates:
            overall_overrides.update(ablation_by_name.get(name, {}))
            overall_overrides.update(stage2_shared_overrides.get(name, {}))
        overall_cfg = _materialize_config_for_export(cfg, base_dir=config_path.parent)
        overall_cfg.update(overall_overrides)
        overall_override_path = artifacts_root / "minimal_both_stages_override.yaml"
        with overall_override_path.open("w", encoding="utf-8") as handle:
            handle.write("# Generated by ablation_stage1.py\n")
            handle.write(f"# Base config: {config_path}\n")
            handle.write("# Source: intersection of negligible components in stage1+stage2\n")
            handle.write(
                "# Negligible thresholds: "
                f"delta_auc < {negligible_auc:.6f} AND delta_mcc < {negligible_mcc:.6f}\n"
            )
            handle.write("# Candidates: " + (", ".join(overall_candidates) if overall_candidates else "<none>") + "\n")
            yaml.safe_dump(overall_cfg, handle, sort_keys=False)
    else:
        overall_override_path = None

    print("\n" + "=" * 80)
    print("STAGE 1 ABLATION SUMMARY")
    print("=" * 80)
    print(f"Config: {config_path}")
    print(f"Seeds: {','.join(str(s) for s in seeds)}")
    print(
        "Baseline  "
        f"AUC {_format_pm(baseline.auc_mean, baseline.auc_std)} | "
        f"ACC {_format_pm(baseline.acc_mean, baseline.acc_std)} | "
        f"MCC {_format_pm(baseline.mcc_mean, baseline.mcc_std)}"
    )
    print("-" * 80)
    print("Ablations (sorted by ΔAUC drop vs baseline)")
    for row in ablation_rows:
        print(
            f"{row.name:<24} "
            f"ΔAUC {row.delta_auc_mean:+.4f} | "
            f"ΔACC {row.delta_acc_mean:+.4f} | "
            f"ΔMCC {row.delta_mcc_mean:+.4f} | "
            f"AUC {_format_pm(row.auc_mean, row.auc_std)} | "
            f"ACC {_format_pm(row.acc_mean, row.acc_std)} | "
            f"MCC {_format_pm(row.mcc_mean, row.mcc_std)}"
        )
    print("-" * 80)
    print("Distillation suggestion (keep only what matters)")
    print(
        f"Negligible thresholds: ΔAUC < {negligible_auc:.4f} AND ΔMCC < {negligible_mcc:.4f}"
    )
    if candidates_to_remove:
        print("Candidates to remove:")
        for row in negligible_rows:
            print(f"- {row.name} (ΔAUC {row.delta_auc_mean:+.4f}, ΔMCC {row.delta_mcc_mean:+.4f})")
    else:
        print("Candidates to remove: <none>")
    if stage2_candidates is None:
        print("Overall candidates (stage 1 + stage 2): <stage 2 CSV not found; run ablation_stage2.py as well>")
    else:
        print(
            "Overall candidates (stage 1 + stage 2): "
            + (", ".join(overall_candidates) if overall_candidates else "<none>")
        )
        if overall_override_path is not None:
            print(f"Wrote YAML override: {overall_override_path}")
    print("-" * 80)
    print(f"Wrote CSV: {out_path}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
