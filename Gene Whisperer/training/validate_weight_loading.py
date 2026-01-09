"""Validate that correct MLM weights are loaded for each k-mer."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import torch
import yaml


def _expected_vocab_size(kmer: int) -> int:
    return (4 ** kmer) + 3


def _load_config(config_path: str) -> dict:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")
    with config_file.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def validate_weight_loading(
    config_path: str = "config.yaml",
    kmers: Optional[Iterable[int]] = None,
) -> bool:
    try:
        cfg = _load_config(config_path)
    except FileNotFoundError as exc:
        print(f"FAIL: {exc}")
        return False

    print("MLM Weight Loading Validation")
    print("=" * 60)

    try:
        from train_stage1 import get_mlm_checkpoint_for_kmer
    except Exception as exc:  # pragma: no cover - import errors should be obvious
        print(f"FAIL: get_mlm_checkpoint_for_kmer import error: {exc}")
        return False

    kmers_to_check = list(kmers) if kmers is not None else [3, 4, 5, 6]
    all_passed = True

    for kmer in kmers_to_check:
        ckpt_path = get_mlm_checkpoint_for_kmer(cfg, kmer)

        if ckpt_path is None:
            print(f"WARN k={kmer}: no checkpoint found")
            continue

        if not ckpt_path.exists():
            print(f"FAIL k={kmer}: checkpoint path not found: {ckpt_path}")
            all_passed = False
            continue

        try:
            checkpoint = torch.load(ckpt_path, map_location="cpu")
        except Exception as exc:
            print(f"FAIL k={kmer}: failed to load checkpoint: {exc}")
            all_passed = False
            continue

        # Only validate metadata when present in the checkpoint.
        ckpt_kmer = checkpoint.get("kmer") if isinstance(checkpoint, dict) else None
        if ckpt_kmer is not None and ckpt_kmer != kmer:
            print(f"FAIL k={kmer}: checkpoint is for k={ckpt_kmer}")
            all_passed = False
            continue

        ckpt_vocab = checkpoint.get("vocab_size") if isinstance(checkpoint, dict) else None
        expected_vocab = _expected_vocab_size(kmer)
        if ckpt_vocab is not None and ckpt_vocab != expected_vocab:
            print(
                f"FAIL k={kmer}: vocab mismatch (expected {expected_vocab}, got {ckpt_vocab})"
            )
            all_passed = False
            continue

        val_loss = checkpoint.get("best_val_loss", "N/A") if isinstance(checkpoint, dict) else "N/A"
        print(f"OK k={kmer}: {ckpt_path.name} (val_loss={val_loss})")

    print("=" * 60)
    if all_passed:
        print("WEIGHT LOADING VALIDATION PASSED")
    else:
        print("SOME VALIDATIONS FAILED")

    return all_passed


if __name__ == "__main__":
    validate_weight_loading()
