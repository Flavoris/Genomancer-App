#!/usr/bin/env python3
__test__ = False

from pathlib import Path

import yaml


def _load_config() -> tuple[dict, Path]:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    cfg_path = repo_root / "gene_whisperer" / "training" / "config.yaml"
    if not cfg_path.exists():
        raise FileNotFoundError(f"config.yaml not found at {cfg_path}")
    with cfg_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    return cfg, cfg_path


def _resolve_cfg_path(cfg_path: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (cfg_path.parent / path).resolve()
    return path


def _ckpt_for_k(mapping: dict, k: int):
    if k in mapping:
        return mapping[k]
    return mapping.get(str(k))


def main() -> None:
    cfg, cfg_path = _load_config()

    use_multi_scale = cfg.get("use_multi_scale")
    assert use_multi_scale is True, f"use_multi_scale must be True, got {use_multi_scale!r}"

    multi_scale_kmers = cfg.get("multi_scale_kmers")
    assert isinstance(multi_scale_kmers, list) and multi_scale_kmers, (
        f"multi_scale_kmers must be a non-empty list, got {multi_scale_kmers!r}"
    )
    assert all(isinstance(k, int) for k in multi_scale_kmers), (
        f"multi_scale_kmers must be list[int], got {multi_scale_kmers!r}"
    )

    mlm_kmer = cfg.get("mlm_kmer")
    stage1_load_mlm_weights = bool(cfg.get("stage1_load_mlm_weights", False))
    print(f"mlm_kmer: {mlm_kmer}")
    print(f"stage1_load_mlm_weights: {stage1_load_mlm_weights}")
    if stage1_load_mlm_weights and any(k != mlm_kmer for k in multi_scale_kmers):
        print(
            "WARNING: multi_scale_kmers includes k != mlm_kmer while "
            "stage1_load_mlm_weights is True"
        )

    stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k")
    assert isinstance(stage1_ckpt_by_k, dict), (
        f"stage1_ckpt_by_k must be a mapping, got {stage1_ckpt_by_k!r}"
    )
    for k in multi_scale_kmers:
        ckpt = _ckpt_for_k(stage1_ckpt_by_k, k)
        assert ckpt is not None, f"stage1_ckpt_by_k missing entry for k={k}"
        assert isinstance(ckpt, str) and ckpt.endswith(".pt"), (
            f"stage1_ckpt_by_k[{k}] must be a .pt path string, got {ckpt!r}"
        )

    vocab_cache_dir = cfg.get("vocab_cache_dir")
    assert isinstance(vocab_cache_dir, str) and vocab_cache_dir, (
        f"vocab_cache_dir must be a non-empty string, got {vocab_cache_dir!r}"
    )
    vocab_path = _resolve_cfg_path(cfg_path, vocab_cache_dir)
    try:
        vocab_path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise AssertionError(f"vocab_cache_dir is not creatable: {vocab_path}") from exc
    assert vocab_path.exists() and vocab_path.is_dir(), (
        f"vocab_cache_dir does not exist or is not a directory: {vocab_path}"
    )


if __name__ == "__main__":
    main()
