#!/usr/bin/env python3
from __future__ import annotations

__test__ = False

import argparse
from pathlib import Path
from typing import Any

import yaml


def _load_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"config.yaml not found at {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"Expected config to be a mapping, got {type(cfg).__name__}")
    return cfg


def _resolve_config_path(config_arg: str, training_dir: Path) -> Path:
    config_path = Path(config_arg).expanduser()
    if not config_path.is_absolute():
        config_path = (Path.cwd() / config_path).resolve()
    if config_path.exists():
        return config_path
    fallback = (training_dir / config_arg).resolve()
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"config.yaml not found at {config_arg}")


def _resolve_path(base_dir: Path, value: str | Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _get_kmers(cfg: dict[str, Any]) -> list[int]:
    kmers = cfg.get("multi_scale_kmers")
    if not isinstance(kmers, list) or not kmers:
        raise ValueError(f"multi_scale_kmers must be a non-empty list, got {kmers!r}")
    if not all(isinstance(k, int) for k in kmers):
        raise ValueError(f"multi_scale_kmers must be list[int], got {kmers!r}")
    return kmers


def main() -> int:
    parser = argparse.ArgumentParser(description="Check multiscale artifacts exist")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (defaults to training/config.yaml).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    training_dir = repo_root / "gene_whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training dir not found: {training_dir}")

    config_path = _resolve_config_path(args.config, training_dir)
    cfg = _load_config(config_path)
    kmers = _get_kmers(cfg)

    stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k")
    if not isinstance(stage1_ckpt_by_k, dict):
        raise ValueError("stage1_ckpt_by_k missing or invalid in config")

    vocab_cache_dir = cfg.get("vocab_cache_dir")
    if not vocab_cache_dir:
        raise ValueError("vocab_cache_dir missing from config")

    base_dir = config_path.parent
    vocab_dir = _resolve_path(base_dir, vocab_cache_dir)

    missing: list[str] = []

    for k in kmers:
        ckpt_value = stage1_ckpt_by_k.get(k)
        if ckpt_value is None:
            ckpt_value = stage1_ckpt_by_k.get(str(k))
        if ckpt_value is None:
            missing.append(f"k={k} missing stage1_ckpt_by_k entry")
        else:
            ckpt_path = _resolve_path(base_dir, ckpt_value)
            if not ckpt_path.exists():
                missing.append(f"k={k} missing checkpoint: {ckpt_path}")

        vocab_path = vocab_dir / f"k{k}_vocab.json"
        if not vocab_path.exists():
            missing.append(f"k={k} missing vocab: {vocab_path}")

    if missing:
        print("FAIL: multiscale artifacts missing")
        for entry in missing:
            print(f"- {entry}")
        return 1

    print(f"PASS: verified {len(kmers)} k-mer artifact sets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
