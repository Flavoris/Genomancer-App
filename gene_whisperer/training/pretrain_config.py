"""Configuration helpers for MLM pretraining."""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import List

import yaml


@dataclass
class PretrainConfig:
    fasta_paths: List[Path]
    tokenizer_path: Path
    vocab_size: int
    max_length: int
    window_size: int
    mask_prob: float
    max_bases_per_file: int | None
    tokenizer_max_bases: int | None
    tokenizer_max_sequences: int | None
    tokenizer_window_size: int
    tokenizer_min_freq: int
    tokenizer_max_token_length: int
    tokenizer_retrain_if_mismatch: bool
    sample_by_length: bool
    mask_ambiguous_tokens: bool
    min_masked_tokens: int
    min_maskable_tokens: int
    resample_attempts: int
    batch_size: int
    epochs: int
    min_epochs: int
    early_stopping_patience: int
    early_stopping_min_delta: float
    save_best_only: bool
    samples_per_epoch: int
    log_interval: int
    num_workers: int
    warmup_ratio: float
    min_lr_ratio: float
    max_grad_norm: float
    use_amp: bool
    lr: float
    weight_decay: float
    grad_accum_steps: int
    embedding_dim: int
    num_layers: int
    num_heads: int
    ff_dim: int
    dropout: float
    seed: int
    output_dir: Path


def _resolve_path(path_str: str, base_dir: Path) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _optional_int(value: object) -> int | None:
    if value is None:
        return None
    return int(value)


def load_pretrain_config(path: Path) -> PretrainConfig:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    config_dir = path.parent
    mlm = data["mlm"]
    model = data["model"]
    training = data["training"]

    return PretrainConfig(
        fasta_paths=[_resolve_path(str(p), config_dir) for p in mlm["fasta_paths"]],
        tokenizer_path=_resolve_path(str(mlm["tokenizer_path"]), config_dir),
        vocab_size=int(mlm["vocab_size"]),
        max_length=int(mlm["max_length"]),
        window_size=int(mlm["window_size"]),
        mask_prob=float(mlm.get("mask_prob", 0.15)),
        max_bases_per_file=_optional_int(mlm.get("max_bases_per_file")),
        tokenizer_max_bases=_optional_int(mlm.get("tokenizer_max_bases")),
        tokenizer_max_sequences=_optional_int(mlm.get("tokenizer_max_sequences")),
        tokenizer_window_size=max(64, int(mlm.get("tokenizer_window_size", 2048))),
        tokenizer_min_freq=max(2, int(mlm.get("tokenizer_min_freq", 8))),
        tokenizer_max_token_length=max(2, int(mlm.get("tokenizer_max_token_length", 16))),
        tokenizer_retrain_if_mismatch=bool(mlm.get("tokenizer_retrain_if_mismatch", True)),
        sample_by_length=bool(mlm.get("sample_by_length", True)),
        mask_ambiguous_tokens=bool(mlm.get("mask_ambiguous_tokens", False)),
        min_masked_tokens=max(1, int(mlm.get("min_masked_tokens", 1))),
        min_maskable_tokens=max(1, int(mlm.get("min_maskable_tokens", mlm.get("min_masked_tokens", 1)))),
        resample_attempts=max(1, int(mlm.get("resample_attempts", 8))),
        batch_size=int(training["batch_size"]),
        epochs=max(1, int(training["epochs"])),
        min_epochs=max(1, int(training.get("min_epochs", 1))),
        early_stopping_patience=max(1, int(training.get("early_stopping_patience", 10))),
        early_stopping_min_delta=max(0.0, float(training.get("early_stopping_min_delta", 0.0))),
        save_best_only=bool(training.get("save_best_only", True)),
        samples_per_epoch=max(1, int(training.get("samples_per_epoch", 4096))),
        log_interval=max(1, int(training.get("log_interval", 50))),
        num_workers=max(0, int(training.get("num_workers", 0))),
        warmup_ratio=min(max(float(training.get("warmup_ratio", 0.03)), 0.0), 0.5),
        min_lr_ratio=min(max(float(training.get("min_lr_ratio", 0.05)), 0.0), 1.0),
        max_grad_norm=max(0.0, float(training.get("max_grad_norm", 1.0))),
        use_amp=bool(training.get("use_amp", True)),
        lr=float(training["lr"]),
        weight_decay=float(training["weight_decay"]),
        grad_accum_steps=max(1, int(training.get("grad_accum_steps", 1))),
        embedding_dim=int(model["embedding_dim"]),
        num_layers=int(model["num_layers"]),
        num_heads=int(model["num_heads"]),
        ff_dim=int(model["ff_dim"]),
        dropout=float(model.get("dropout", 0.1)),
        seed=int(training.get("seed", 1337)),
        output_dir=_resolve_path(str(training["output_dir"]), config_dir),
    )


def sample_tokenizer_corpus(
    sequences: List[str],
    seed: int,
    max_bases: int | None,
    max_sequences: int | None,
    window_size: int,
) -> List[str]:
    if max_bases is None and max_sequences is None:
        return sequences

    if not sequences:
        return []

    rng = random.Random(seed)
    indices = list(range(len(sequences)))
    rng.shuffle(indices)

    sampled: List[str] = []
    consumed_bases = 0
    sampled_sequences = 0

    for idx in indices:
        if max_sequences is not None and sampled_sequences >= max_sequences:
            break
        if max_bases is not None and consumed_bases >= max_bases:
            break

        seq = sequences[idx]
        if not seq:
            continue

        span = min(window_size, len(seq))
        if span <= 0:
            continue

        if len(seq) == span:
            window = seq
        else:
            start = rng.randint(0, len(seq) - span)
            window = seq[start : start + span]

        if max_bases is not None:
            remaining = max_bases - consumed_bases
            if remaining <= 0:
                break
            if len(window) > remaining:
                window = window[:remaining]

        if not window:
            continue

        sampled.append(window)
        consumed_bases += len(window)
        sampled_sequences += 1

    return sampled
