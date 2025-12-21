from __future__ import annotations

import argparse
import copy
import csv
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

VARIANTS: list[tuple[str, dict[str, Any]]] = [
    ("baseline", {}),
    ("no_alibi", {"use_alibi": False}),
    ("no_attention_pool", {"use_attention_pool": False}),
    ("no_tcn", {"use_tcn": False}),
    ("no_postcnn_transformer", {"post_cnn_transformer_layers": 0}),
    (
        "no_engineered_features",
        {"stage1_use_engineered_features": False, "stage2_use_engineered_features": False},
    ),
    ("no_tnc", {"stage1_feature_enable_tnc": False}),
    ("no_pstnp", {"stage1_feature_enable_pstnp": False}),
    ("no_pseeiip", {"stage1_feature_enable_pseeiip": False}),
    ("no_cksnap", {"stage1_feature_enable_cksnap": False}),
]

DEFAULT_SCORE_WEIGHTS: dict[str, float] = {
    "pretrain": 0.15,
    "stage1": 0.35,
    "stage2": 0.35,
    "transfer": 0.15,
    "cost": 0.10,
}

TRANSFER_EFFICIENCY_MODES = {"per_pretrain_second", "per_total_second", "per_deploy_param"}


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


def _parse_variants(value: str) -> list[str]:
    parts = [chunk.strip() for chunk in value.replace(",", " ").split()]
    parts = [chunk for chunk in parts if chunk]
    if not parts:
        raise argparse.ArgumentTypeError("Expected at least one variant name")
    return parts


def _mean_std(values: Iterable[float]) -> tuple[float, float]:
    arr = np.asarray(list(values), dtype=np.float64)
    if arr.size == 0:
        return float("nan"), float("nan")
    mean = float(arr.mean())
    if arr.size == 1:
        return mean, 0.0
    return mean, float(arr.std(ddof=1))


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


def _coerce_seed_list(value: Any) -> list[int] | None:
    if isinstance(value, list):
        seeds: list[int] = []
        for item in value:
            try:
                seeds.append(int(item))
            except (TypeError, ValueError):
                continue
        return seeds or None
    return None


def _resolve_score_weights(ablation_cfg: dict[str, Any]) -> dict[str, float]:
    weights = dict(DEFAULT_SCORE_WEIGHTS)
    raw = ablation_cfg.get("score_weights")
    if isinstance(raw, dict):
        for key in weights:
            if key in raw:
                try:
                    weights[key] = float(raw[key])
                except (TypeError, ValueError):
                    continue
    return weights


def _normalize_value(value: float, min_val: float, max_val: float, *, higher_better: bool) -> float:
    if not math.isfinite(value):
        return 0.0
    if max_val <= min_val:
        norm = 0.5
    else:
        norm = (value - min_val) / (max_val - min_val)
    if not higher_better:
        norm = 1.0 - norm
    return max(0.0, min(1.0, norm))


def _metric_ranges(rows: list[dict[str, Any]], *, seeds: list[int], key: str) -> dict[int, tuple[float, float]]:
    ranges: dict[int, tuple[float, float]] = {}
    for seed in seeds:
        values = [
            float(row[key])
            for row in rows
            if row.get("seed") == seed and math.isfinite(float(row.get(key, float("nan"))))
        ]
        if values:
            ranges[seed] = (min(values), max(values))
        else:
            ranges[seed] = (0.0, 0.0)
    return ranges


def _score_rows(
    rows: list[dict[str, Any]],
    *,
    seeds: list[int],
    weights: dict[str, float],
    ranges: dict[str, dict[int, tuple[float, float]]] | None = None,
) -> dict[str, dict[int, tuple[float, float]]]:
    if ranges is None:
        ranges = {
            "stage1_mcc": _metric_ranges(rows, seeds=seeds, key="stage1_mcc"),
            "stage2_mcc": _metric_ranges(rows, seeds=seeds, key="stage2_mcc"),
            "transfer_efficiency": _metric_ranges(rows, seeds=seeds, key="transfer_efficiency"),
            "mlm_loss": _metric_ranges(rows, seeds=seeds, key="mlm_loss"),
            "cost": _metric_ranges(rows, seeds=seeds, key="cost"),
        }

    for row in rows:
        seed = int(row["seed"])
        stage1_min, stage1_max = ranges["stage1_mcc"][seed]
        stage2_min, stage2_max = ranges["stage2_mcc"][seed]
        transfer_min, transfer_max = ranges["transfer_efficiency"][seed]
        loss_min, loss_max = ranges["mlm_loss"][seed]
        cost_min, cost_max = ranges["cost"][seed]

        norm_stage1 = _normalize_value(
            float(row["stage1_mcc"]), stage1_min, stage1_max, higher_better=True
        )
        norm_stage2 = _normalize_value(
            float(row["stage2_mcc"]), stage2_min, stage2_max, higher_better=True
        )
        norm_transfer = _normalize_value(
            float(row["transfer_efficiency"]), transfer_min, transfer_max, higher_better=True
        )
        norm_loss = _normalize_value(
            float(row["mlm_loss"]), loss_min, loss_max, higher_better=False
        )
        norm_cost = _normalize_value(
            float(row["cost"]), cost_min, cost_max, higher_better=False
        )

        score = (
            weights.get("pretrain", 0.0) * norm_loss
            + weights.get("stage1", 0.0) * norm_stage1
            + weights.get("stage2", 0.0) * norm_stage2
            + weights.get("transfer", 0.0) * norm_transfer
            + weights.get("cost", 0.0) * norm_cost
        )

        row["norm_mlm_loss"] = norm_loss
        row["norm_stage1_mcc"] = norm_stage1
        row["norm_stage2_mcc"] = norm_stage2
        row["norm_transfer_efficiency"] = norm_transfer
        row["norm_cost"] = norm_cost
        row["score"] = score

    return ranges


def _finite(values: Iterable[float]) -> list[float]:
    return [float(v) for v in values if math.isfinite(float(v))]


def _summarize_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_variant: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_variant.setdefault(str(row["variant"]), []).append(row)

    summary_rows: list[dict[str, Any]] = []
    for variant, variant_rows in by_variant.items():
        seeds = sorted({int(row["seed"]) for row in variant_rows})
        summary_rows.append(
            {
                "variant": variant,
                "seeds": ",".join(str(s) for s in seeds),
                "n_seeds": len(seeds),
                "mlm_loss_mean": _mean_std(_finite(row["mlm_loss"] for row in variant_rows))[0],
                "mlm_loss_std": _mean_std(_finite(row["mlm_loss"] for row in variant_rows))[1],
                "stage1_mcc_mean": _mean_std(_finite(row["stage1_mcc"] for row in variant_rows))[0],
                "stage1_mcc_std": _mean_std(_finite(row["stage1_mcc"] for row in variant_rows))[1],
                "stage1_auc_mean": _mean_std(_finite(row["stage1_auc"] for row in variant_rows))[0],
                "stage1_auc_std": _mean_std(_finite(row["stage1_auc"] for row in variant_rows))[1],
                "stage2_mcc_mean": _mean_std(_finite(row["stage2_mcc"] for row in variant_rows))[0],
                "stage2_mcc_std": _mean_std(_finite(row["stage2_mcc"] for row in variant_rows))[1],
                "stage2_auc_mean": _mean_std(_finite(row["stage2_auc"] for row in variant_rows))[0],
                "stage2_auc_std": _mean_std(_finite(row["stage2_auc"] for row in variant_rows))[1],
                "transfer_efficiency_mean": _mean_std(
                    _finite(row["transfer_efficiency"] for row in variant_rows)
                )[0],
                "transfer_efficiency_std": _mean_std(
                    _finite(row["transfer_efficiency"] for row in variant_rows)
                )[1],
                "deploy_params_mean": _mean_std(
                    _finite(row["deploy_params"] for row in variant_rows)
                )[0],
                "deploy_params_std": _mean_std(
                    _finite(row["deploy_params"] for row in variant_rows)
                )[1],
                "train_seconds_total_mean": _mean_std(
                    _finite(row["train_seconds_total"] for row in variant_rows)
                )[0],
                "train_seconds_total_std": _mean_std(
                    _finite(row["train_seconds_total"] for row in variant_rows)
                )[1],
                "cost_mean": _mean_std(_finite(row["cost"] for row in variant_rows))[0],
                "cost_std": _mean_std(_finite(row["cost"] for row in variant_rows))[1],
                "score_mean": _mean_std(_finite(row["score"] for row in variant_rows))[0],
                "score_std": _mean_std(_finite(row["score"] for row in variant_rows))[1],
            }
        )

    return summary_rows


def _format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6f}"
    return str(value)


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_value(row.get(key)) for key in fieldnames})


def _pareto_frontier(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates = []
    for row in summary_rows:
        try:
            s1 = float(row["stage1_mcc_mean"])
            s2 = float(row["stage2_mcc_mean"])
            neg_loss = -float(row["mlm_loss_mean"])
            neg_cost = -float(row["cost_mean"])
            if not all(math.isfinite(v) for v in (s1, s2, neg_loss, neg_cost)):
                continue
            candidates.append((row, s1, s2, neg_loss, neg_cost))
        except (TypeError, ValueError):
            continue

    frontier: list[dict[str, Any]] = []
    for row, s1, s2, neg_loss, neg_cost in candidates:
        dominated = False
        for other_row, o1, o2, o_neg_loss, o_neg_cost in candidates:
            if other_row is row:
                continue
            if (
                o1 >= s1
                and o2 >= s2
                and o_neg_loss >= neg_loss
                and o_neg_cost >= neg_cost
                and (o1 > s1 or o2 > s2 or o_neg_loss > neg_loss or o_neg_cost > neg_cost)
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(row)
    frontier.sort(key=lambda r: float(r["score_mean"]), reverse=True)
    return frontier


def _expect_float(report: dict[str, Any], key: str, label: str) -> float:
    if key not in report:
        raise RuntimeError(f"Missing {key} in {label} report: {report}")
    return float(report[key])


def _run_full_pipeline(
    *,
    base_cfg: dict[str, Any],
    variant_name: str,
    variant_overrides: dict[str, Any],
    seed: int,
    config_path: Path,
    compute_transfer_efficiency: bool,
    transfer_efficiency_mode: str,
    run_mlm_pretrain: Any,
    run_stage1_training: Any,
    run_stage2_training: Any,
) -> dict[str, Any]:
    cfg_run = copy.deepcopy(base_cfg)
    cfg_run["seed"] = seed
    ablation_cfg = cfg_run.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
        cfg_run["ablation"] = ablation_cfg
    ablation_cfg["enabled"] = True

    base_overrides: dict[str, Any] = {"_config_path": str(config_path), "seed": seed}
    pretrain_overrides = dict(base_overrides)
    pretrain_overrides.update(variant_overrides)

    pretrain_report = run_mlm_pretrain(cfg_run, overrides=pretrain_overrides)
    encoder_ckpt_path = str(pretrain_report.get("encoder_ckpt_path", ""))
    mlm_loss = _expect_float(pretrain_report, "mlm_loss", f"{variant_name} pretrain")
    pretrain_seconds = _expect_float(pretrain_report, "pretrain_seconds", f"{variant_name} pretrain")
    pretrain_params = int(pretrain_report.get("params", 0))

    stage1_transfer_mode = variant_overrides.get(
        "stage1_mlm_transfer_mode",
        cfg_run.get("stage1_mlm_transfer_mode", "embed_only"),
    )
    stage1_overrides = dict(base_overrides)
    stage1_overrides.update(variant_overrides)
    stage1_overrides.update(
        {
            "mlm_encoder_path": encoder_ckpt_path,
            "stage1_load_mlm_weights": True,
            "stage1_mlm_transfer_mode": stage1_transfer_mode,
        }
    )
    stage1_report = run_stage1_training(cfg_run, overrides=stage1_overrides)
    stage1_mcc = _expect_float(stage1_report, "stage1_mcc", f"{variant_name} stage1")
    stage1_auc = _expect_float(stage1_report, "stage1_auc", f"{variant_name} stage1")
    stage1_acc = _expect_float(stage1_report, "stage1_acc", f"{variant_name} stage1")
    stage1_seconds = float(stage1_report.get("stage1_seconds", 0.0))
    stage1_params = int(stage1_report.get("params", 0))

    stage2_transfer_mode = variant_overrides.get(
        "stage2_mlm_transfer_mode",
        cfg_run.get("stage2_mlm_transfer_mode", "embed_only"),
    )
    stage2_overrides = dict(base_overrides)
    stage2_overrides.update(variant_overrides)
    stage2_overrides.update(
        {
            "mlm_encoder_path": encoder_ckpt_path,
            "stage2_load_mlm_weights": True,
            "stage2_mlm_transfer_mode": stage2_transfer_mode,
        }
    )
    stage2_report = run_stage2_training(cfg_run, overrides=stage2_overrides)
    stage2_mcc = _expect_float(stage2_report, "stage2_mcc", f"{variant_name} stage2")
    stage2_auc = _expect_float(stage2_report, "stage2_auc", f"{variant_name} stage2")
    stage2_acc = _expect_float(stage2_report, "stage2_acc", f"{variant_name} stage2")
    stage2_seconds = float(stage2_report.get("stage2_seconds", 0.0))
    stage2_params = int(stage2_report.get("params", 0))

    deploy_params = max(stage1_params, stage2_params)
    train_seconds_total = pretrain_seconds + stage1_seconds + stage2_seconds

    stage1_mcc_scratch = float("nan")
    stage1_auc_scratch = float("nan")
    stage1_acc_scratch = float("nan")
    stage2_mcc_scratch = float("nan")
    stage2_auc_scratch = float("nan")
    stage2_acc_scratch = float("nan")
    transfer_gain_mcc = 0.0
    transfer_gain_auc = 0.0
    transfer_efficiency = 0.0

    if compute_transfer_efficiency:
        stage1_scratch_overrides = dict(base_overrides)
        stage1_scratch_overrides.update(variant_overrides)
        stage1_scratch_overrides["stage1_load_mlm_weights"] = False
        stage1_scratch_overrides["stage1_mlm_transfer_mode"] = "none"
        stage1_scratch_report = run_stage1_training(cfg_run, overrides=stage1_scratch_overrides)
        stage1_mcc_scratch = _expect_float(
            stage1_scratch_report, "stage1_mcc", f"{variant_name} stage1 scratch"
        )
        stage1_auc_scratch = _expect_float(
            stage1_scratch_report, "stage1_auc", f"{variant_name} stage1 scratch"
        )
        stage1_acc_scratch = _expect_float(
            stage1_scratch_report, "stage1_acc", f"{variant_name} stage1 scratch"
        )

        stage2_scratch_overrides = dict(base_overrides)
        stage2_scratch_overrides.update(variant_overrides)
        stage2_scratch_overrides["stage2_load_mlm_weights"] = False
        stage2_scratch_overrides["stage2_mlm_transfer_mode"] = "none"
        stage2_scratch_report = run_stage2_training(cfg_run, overrides=stage2_scratch_overrides)
        stage2_mcc_scratch = _expect_float(
            stage2_scratch_report, "stage2_mcc", f"{variant_name} stage2 scratch"
        )
        stage2_auc_scratch = _expect_float(
            stage2_scratch_report, "stage2_auc", f"{variant_name} stage2 scratch"
        )
        stage2_acc_scratch = _expect_float(
            stage2_scratch_report, "stage2_acc", f"{variant_name} stage2 scratch"
        )

        transfer_gain_mcc = 0.5 * (
            (stage1_mcc - stage1_mcc_scratch) + (stage2_mcc - stage2_mcc_scratch)
        )
        transfer_gain_auc = 0.5 * (
            (stage1_auc - stage1_auc_scratch) + (stage2_auc - stage2_auc_scratch)
        )
        if transfer_efficiency_mode == "per_deploy_param":
            denom = max(float(deploy_params), 1e-6)
        elif transfer_efficiency_mode in {"per_pretrain_second", "per_total_second"}:
            denom = max(train_seconds_total, 1e-9)
        else:
            raise ValueError(f"Unknown transfer_efficiency_mode: {transfer_efficiency_mode}")
        transfer_efficiency = transfer_gain_mcc / denom

    cost = float(deploy_params)

    return {
        "variant": variant_name,
        "seed": seed,
        "mlm_loss": mlm_loss,
        "pretrain_seconds": pretrain_seconds,
        "pretrain_params": pretrain_params,
        "encoder_ckpt_path": encoder_ckpt_path,
        "stage1_acc": stage1_acc,
        "stage1_auc": stage1_auc,
        "stage1_mcc": stage1_mcc,
        "stage1_seconds": stage1_seconds,
        "stage1_params": stage1_params,
        "stage2_acc": stage2_acc,
        "stage2_auc": stage2_auc,
        "stage2_mcc": stage2_mcc,
        "stage2_seconds": stage2_seconds,
        "stage2_params": stage2_params,
        "deploy_params": deploy_params,
        "train_seconds_total": train_seconds_total,
        "stage1_acc_scratch": stage1_acc_scratch,
        "stage1_auc_scratch": stage1_auc_scratch,
        "stage1_mcc_scratch": stage1_mcc_scratch,
        "stage2_acc_scratch": stage2_acc_scratch,
        "stage2_auc_scratch": stage2_auc_scratch,
        "stage2_mcc_scratch": stage2_mcc_scratch,
        "transfer_gain_mcc": transfer_gain_mcc,
        "transfer_gain_auc": transfer_gain_auc,
        "transfer_efficiency": transfer_efficiency,
        "cost": cost,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Full pipeline ablation runner")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--variants", type=_parse_variants, help="Comma-separated variants list")
    parser.add_argument("--seeds", type=_parse_seeds, help="Comma-separated list of seeds")
    args = parser.parse_args()

    initial_cwd = Path.cwd()
    repo_root = Path(__file__).resolve().parent
    training_dir = repo_root / "Gene Whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    sys.path.insert(0, str(training_dir))
    from pretrain_mlm import run_mlm_pretrain  # type: ignore
    from train_stage1 import load_config, run_stage1_training  # type: ignore
    from train_stage2 import run_stage2_training  # type: ignore

    config_path = _resolve_config_path(args.config, training_dir=training_dir, initial_cwd=initial_cwd)
    cfg: dict[str, Any] = load_config(config_path)

    ablation_cfg = cfg.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
        cfg["ablation"] = ablation_cfg
    base_variant_name = ablation_cfg.get("distill_base_variant", "baseline")
    if base_variant_name is None:
        base_variant_name = "baseline"
    base_variant_name = str(base_variant_name).strip() or "baseline"

    if args.seeds is None:
        seeds = _coerce_seed_list(ablation_cfg.get("seeds")) or [1337]
    else:
        seeds = args.seeds

    variants_by_name = dict(VARIANTS)
    if base_variant_name not in variants_by_name:
        raise SystemExit(
            "Unknown ablation.distill_base_variant: "
            f"{base_variant_name!r}. Choose from: {sorted(variants_by_name)}"
        )
    if args.variants is None:
        selected_variants = VARIANTS
    else:
        requested = []
        seen = set()
        for name in args.variants:
            if name == "diet_combo":
                raise SystemExit("diet_combo is generated automatically; do not pass it in --variants.")
            if name not in variants_by_name:
                raise SystemExit(f"Unknown variant: {name}")
            if name not in seen:
                requested.append(name)
                seen.add(name)
        if base_variant_name not in seen:
            print(f"Note: adding {base_variant_name} to variants for scoring and distillation.")
            requested.insert(0, base_variant_name)
        selected_variants = [(name, variants_by_name[name]) for name in requested]

    compute_transfer_efficiency = bool(ablation_cfg.get("compute_transfer_efficiency", False))
    transfer_efficiency_mode = str(ablation_cfg.get("transfer_efficiency_mode", "per_total_second"))
    if transfer_efficiency_mode not in TRANSFER_EFFICIENCY_MODES:
        raise SystemExit(
            "Unknown ablation.transfer_efficiency_mode: "
            f"{transfer_efficiency_mode!r}. Choose from: {sorted(TRANSFER_EFFICIENCY_MODES)}"
        )
    score_weights = _resolve_score_weights(ablation_cfg)
    negligible_score_drop = float(ablation_cfg.get("negligible_score_drop", 0.01))
    min_cost_improvement = float(ablation_cfg.get("min_cost_improvement", 0.03))

    os.chdir(config_path.parent)

    run_rows: list[dict[str, Any]] = []
    for name, overrides in selected_variants:
        for seed in seeds:
            run_rows.append(
                _run_full_pipeline(
                    base_cfg=cfg,
                    variant_name=name,
                    variant_overrides=overrides,
                    seed=seed,
                    config_path=config_path,
                    compute_transfer_efficiency=compute_transfer_efficiency,
                    transfer_efficiency_mode=transfer_efficiency_mode,
                    run_mlm_pretrain=run_mlm_pretrain,
                    run_stage1_training=run_stage1_training,
                    run_stage2_training=run_stage2_training,
                )
            )

    ranges = _score_rows(run_rows, seeds=seeds, weights=score_weights)
    summary_rows = _summarize_rows(run_rows)

    summary_by_name = {row["variant"]: row for row in summary_rows}
    base_summary = summary_by_name.get(base_variant_name)
    base_score = float("nan")
    base_cost = float("nan")
    if base_summary is not None:
        base_score = float(base_summary.get("score_mean", float("nan")))
        base_cost = float(base_summary.get("cost_mean", float("nan")))

    removable: list[dict[str, Any]] = []
    removable_variants: list[str] = []
    if math.isfinite(base_score) and math.isfinite(base_cost):
        for name, overrides in selected_variants:
            if name == base_variant_name:
                continue
            summary = summary_by_name.get(name)
            if summary is None:
                continue
            score_mean = float(summary["score_mean"])
            cost_mean = float(summary["cost_mean"])
            if not math.isfinite(score_mean) or not math.isfinite(cost_mean):
                continue
            score_drop = base_score - score_mean
            if score_mean >= base_score:
                score_drop = 0.0
            cost_improvement = 0.0
            if base_cost > 0:
                cost_improvement = (base_cost - cost_mean) / base_cost
            if score_drop <= negligible_score_drop and cost_improvement >= min_cost_improvement:
                removable.append(
                    {
                        "name": name,
                        "score_drop": score_drop,
                        "cost_improvement": cost_improvement,
                        "score_mean": score_mean,
                        "cost_mean": cost_mean,
                        "delta_score": score_mean - base_score,
                        "delta_cost": cost_mean - base_cost,
                    }
                )
                removable_variants.append(name)

    minimal_overrides: dict[str, Any] = {}
    for name in removable_variants:
        minimal_overrides.update(variants_by_name[name])

    artifacts_root = (training_dir / "../artifacts").resolve()
    artifacts_root.mkdir(parents=True, exist_ok=True)

    minimal_cfg = _materialize_config_for_export(minimal_overrides, base_dir=config_path.parent)
    minimal_override_path = artifacts_root / "minimal_pipeline_override.yaml"
    with minimal_override_path.open("w", encoding="utf-8") as handle:
        handle.write("# Generated by ablation_all.py\n")
        handle.write(f"# Base config: {config_path}\n")
        handle.write(
            "# Criteria: score_drop <= "
            f"{negligible_score_drop:.6f} AND cost_improvement >= {min_cost_improvement:.6f}\n"
        )
        handle.write(
            "# Removable: " + (", ".join(r["name"] for r in removable) if removable else "none") + "\n"
        )
        yaml.safe_dump(minimal_cfg, handle, sort_keys=False)

    diet_seed = seeds[0]
    diet_row = _run_full_pipeline(
        base_cfg=cfg,
        variant_name="diet_combo",
        variant_overrides=minimal_overrides,
        seed=diet_seed,
        config_path=config_path,
        compute_transfer_efficiency=compute_transfer_efficiency,
        transfer_efficiency_mode=transfer_efficiency_mode,
        run_mlm_pretrain=run_mlm_pretrain,
        run_stage1_training=run_stage1_training,
        run_stage2_training=run_stage2_training,
    )
    run_rows.append(diet_row)
    # Reuse variant ranges so diet_combo is comparable without reshuffling scores.
    _score_rows([diet_row], seeds=[diet_seed], weights=score_weights, ranges=ranges)

    summary_rows = _summarize_rows(run_rows)

    run_csv_path = artifacts_root / "ablation_all_runs.csv"
    summary_csv_path = artifacts_root / "ablation_all_summary.csv"

    run_fieldnames = [
        "variant",
        "seed",
        "mlm_loss",
        "pretrain_seconds",
        "pretrain_params",
        "encoder_ckpt_path",
        "stage1_acc",
        "stage1_auc",
        "stage1_mcc",
        "stage1_seconds",
        "stage1_params",
        "stage2_acc",
        "stage2_auc",
        "stage2_mcc",
        "stage2_seconds",
        "stage2_params",
        "deploy_params",
        "train_seconds_total",
        "stage1_acc_scratch",
        "stage1_auc_scratch",
        "stage1_mcc_scratch",
        "stage2_acc_scratch",
        "stage2_auc_scratch",
        "stage2_mcc_scratch",
        "transfer_gain_mcc",
        "transfer_gain_auc",
        "transfer_efficiency",
        "cost",
        "norm_mlm_loss",
        "norm_stage1_mcc",
        "norm_stage2_mcc",
        "norm_transfer_efficiency",
        "norm_cost",
        "score",
    ]
    summary_fieldnames = [
        "variant",
        "seeds",
        "n_seeds",
        "mlm_loss_mean",
        "mlm_loss_std",
        "stage1_mcc_mean",
        "stage1_mcc_std",
        "stage1_auc_mean",
        "stage1_auc_std",
        "stage2_mcc_mean",
        "stage2_mcc_std",
        "stage2_auc_mean",
        "stage2_auc_std",
        "transfer_efficiency_mean",
        "transfer_efficiency_std",
        "deploy_params_mean",
        "deploy_params_std",
        "train_seconds_total_mean",
        "train_seconds_total_std",
        "cost_mean",
        "cost_std",
        "score_mean",
        "score_std",
    ]

    _write_csv(run_csv_path, run_rows, run_fieldnames)
    _write_csv(summary_csv_path, summary_rows, summary_fieldnames)

    summary_rows_sorted = sorted(
        summary_rows, key=lambda r: float(r["score_mean"]) if math.isfinite(float(r["score_mean"])) else -1e9, reverse=True
    )
    top_rows = summary_rows_sorted[:5]
    pareto_rows = _pareto_frontier(summary_rows)

    print("\n" + "=" * 80)
    print("ABLATION PIPELINE SUMMARY")
    print("=" * 80)
    print(f"Config: {config_path}")
    print(f"Seeds: {','.join(str(s) for s in seeds)}")
    variant_names = [name for name, _ in selected_variants] + ["diet_combo"]
    print(f"Variants: {', '.join(variant_names)}")
    print("-" * 80)
    print("Top 5 by composite score (mean across seeds)")
    for row in top_rows:
        print(
            f"{row['variant']:<24} "
            f"score {float(row['score_mean']):.4f} | "
            f"s1_mcc {float(row['stage1_mcc_mean']):.4f} | "
            f"s2_mcc {float(row['stage2_mcc_mean']):.4f} | "
            f"mlm_loss {float(row['mlm_loss_mean']):.4f} | "
            f"cost {float(row['cost_mean']):.0f}"
        )
    if pareto_rows:
        print("-" * 80)
        print("Pareto frontier (mean across seeds)")
        for row in pareto_rows:
            print(
                f"{row['variant']:<24} "
                f"s1_mcc {float(row['stage1_mcc_mean']):.4f} | "
                f"s2_mcc {float(row['stage2_mcc_mean']):.4f} | "
                f"mlm_loss {float(row['mlm_loss_mean']):.4f} | "
                f"cost {float(row['cost_mean']):.0f}"
            )
    print("-" * 80)
    print("Distillation suggestion")
    if math.isfinite(base_score) and math.isfinite(base_cost):
        print(f"Base ({base_variant_name}) score {base_score:.4f} | cost {base_cost:.0f}")
    print(
        "Criteria: score_drop <= "
        f"{negligible_score_drop:.4f} AND cost_improvement >= {min_cost_improvement:.4f}"
    )
    if removable:
        for row in removable:
            print(
                f"- {row['name']:<24} "
                f"score {row['score_mean']:.4f} | "
                f"cost {row['cost_mean']:.0f} | "
                f"delta_score {row['delta_score']:+.4f} | "
                f"delta_cost {row['delta_cost']:+.0f}"
            )
    else:
        print("Removable components: none")
    print("Minimal overrides (diet_combo):")
    if minimal_overrides:
        formatted_overrides = yaml.safe_dump(minimal_overrides, sort_keys=False).rstrip()
        for line in formatted_overrides.splitlines():
            print(f"  {line}")
    else:
        print("  none")
    print(f"Minimal override: {minimal_override_path}")
    print(f"Diet combo seed: {diet_seed}")
    print("-" * 80)
    print(f"Wrote CSV: {run_csv_path}")
    print(f"Wrote CSV: {summary_csv_path}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
