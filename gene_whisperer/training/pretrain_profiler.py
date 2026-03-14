"""Preflight profiling helpers for MLM pretraining."""
from __future__ import annotations

import json
import random
from bisect import bisect_right
from pathlib import Path
from typing import Any, Mapping, Sequence

import torch

from gene_whisperer.datasets.fasta import load_fasta_sequences
from gene_whisperer.datasets.mlm_dataset import MLMDataset
from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.pretrain_config import PretrainConfig
from gene_whisperer.training.pretrain_metrics import summarize_series


def load_sequences_with_sources(
    paths: Sequence[Path],
    max_bases_per_file: int | None = None,
    verbose: bool = False,
) -> tuple[list[str], list[str], dict[str, dict[str, Any]]]:
    sequences: list[str] = []
    sequence_sources: list[str] = []
    source_summary: dict[str, dict[str, Any]] = {}

    for path in paths:
        source_name = path.name
        source_sequences = load_fasta_sequences(
            [path],
            max_bases_per_file=max_bases_per_file,
            verbose=verbose,
        )
        source_bases = sum(len(sequence) for sequence in source_sequences)
        source_summary[source_name] = {
            "configured_path": str(path),
            "sequence_count": len(source_sequences),
            "bases": source_bases,
        }
        sequences.extend(source_sequences)
        sequence_sources.extend([source_name] * len(source_sequences))

    return sequences, sequence_sources, source_summary


def _sample_tokenizer_corpus_with_source_coverage(
    sequences: Sequence[str],
    sequence_sources: Sequence[str],
    seed: int,
    max_bases: int | None,
    max_sequences: int | None,
    window_size: int,
) -> tuple[list[str], dict[str, Any]]:
    cleaned_sequences: list[str] = []
    cleaned_sources: list[str] = []
    for sequence, source in zip(sequences, sequence_sources):
        if not sequence:
            continue
        cleaned_sequences.append(sequence)
        cleaned_sources.append(source)

    if not cleaned_sequences:
        return [], {"sampled_sequences": 0, "sampled_bases": 0, "sources": {}}

    if max_bases is None and max_sequences is None:
        source_windows: dict[str, int] = {}
        source_bases: dict[str, int] = {}
        for sequence, source in zip(cleaned_sequences, cleaned_sources):
            source_windows[source] = source_windows.get(source, 0) + 1
            source_bases[source] = source_bases.get(source, 0) + len(sequence)
        return list(cleaned_sequences), _format_coverage_report(source_windows, source_bases)

    rng = random.Random(seed)
    cumulative_lengths: list[int] = []
    total_bases = 0
    for sequence in cleaned_sequences:
        total_bases += len(sequence)
        cumulative_lengths.append(total_bases)

    sampled_windows: list[str] = []
    source_window_counts: dict[str, int] = {}
    source_base_counts: dict[str, int] = {}
    consumed_bases = 0

    while True:
        if max_sequences is not None and len(sampled_windows) >= max_sequences:
            break
        if max_bases is not None and consumed_bases >= max_bases:
            break

        sequence_target = rng.randrange(total_bases)
        sequence_idx = bisect_right(cumulative_lengths, sequence_target)
        sequence_idx = min(sequence_idx, len(cleaned_sequences) - 1)
        sequence = cleaned_sequences[sequence_idx]
        source = cleaned_sources[sequence_idx]

        remaining_bases = None if max_bases is None else max_bases - consumed_bases
        span = min(window_size, len(sequence))
        if remaining_bases is not None:
            span = min(span, remaining_bases)
        if span <= 0:
            break

        if len(sequence) == span:
            window = sequence
        else:
            start = rng.randint(0, len(sequence) - span)
            window = sequence[start : start + span]

        sampled_windows.append(window)
        consumed_bases += len(window)
        source_window_counts[source] = source_window_counts.get(source, 0) + 1
        source_base_counts[source] = source_base_counts.get(source, 0) + len(window)

    return sampled_windows, _format_coverage_report(source_window_counts, source_base_counts)


def _format_coverage_report(
    source_window_counts: Mapping[str, int],
    source_base_counts: Mapping[str, int],
) -> dict[str, Any]:
    sampled_sequences = sum(source_window_counts.values())
    sampled_bases = sum(source_base_counts.values())
    sources: dict[str, Any] = {}
    for source_name in sorted(source_window_counts):
        source_sequences = source_window_counts[source_name]
        source_bases = source_base_counts.get(source_name, 0)
        sources[source_name] = {
            "sampled_windows": source_sequences,
            "sampled_bases": source_bases,
            "window_share": (
                source_sequences / sampled_sequences if sampled_sequences else 0.0
            ),
            "base_share": source_bases / sampled_bases if sampled_bases else 0.0,
        }
    return {
        "sampled_sequences": sampled_sequences,
        "sampled_bases": sampled_bases,
        "sources": sources,
    }


def build_preflight_profile(
    config: PretrainConfig,
    sequences: Sequence[str],
    tokenizer: BPETokenizer,
    *,
    sequence_sources: Sequence[str] | None = None,
    source_summary: Mapping[str, Any] | None = None,
    sample_count: int = 512,
) -> dict[str, Any]:
    effective_samples = max(1, min(sample_count, config.samples_per_epoch))
    dataset = MLMDataset(
        sequences=sequences,
        tokenizer=tokenizer,
        window_size=config.window_size,
        max_length=config.max_length,
        mask_prob=config.mask_prob,
        num_samples=effective_samples,
        seed=config.seed,
        sample_by_length=config.sample_by_length,
        mask_ambiguous_tokens=config.mask_ambiguous_tokens,
        min_masked_tokens=config.min_masked_tokens,
        min_maskable_tokens=config.min_maskable_tokens,
        min_tokenized_tokens=config.min_tokenized_tokens,
        resample_attempts=config.resample_attempts,
    )

    tokenized_lengths: list[int] = []
    maskable_counts: list[int] = []
    masked_counts: list[int] = []
    unk_counts: list[int] = []

    for index in range(effective_samples):
        sample = dataset[index]
        tokenized_count = int(sample["tokenized_count"].item())
        maskable_count = int(sample["maskable_count"].item())
        masked_count = int(sample["masked_count"].item())
        active_ids = sample["input_ids"][:tokenized_count]
        unk_count = int(torch.sum(active_ids == tokenizer.unk_token_id).item())

        tokenized_lengths.append(tokenized_count)
        maskable_counts.append(maskable_count)
        masked_counts.append(masked_count)
        unk_counts.append(unk_count)

    total_tokenized = int(sum(tokenized_lengths))
    total_maskable = int(sum(maskable_counts))
    total_masked = int(sum(masked_counts))
    total_unk = int(sum(unk_counts))

    tokenizer_windows: list[str] = []
    tokenizer_coverage = {
        "sampled_sequences": 0,
        "sampled_bases": 0,
        "sources": {},
    }
    if sequence_sources is not None and len(sequence_sources) == len(sequences):
        tokenizer_windows, tokenizer_coverage = _sample_tokenizer_corpus_with_source_coverage(
            sequences=sequences,
            sequence_sources=sequence_sources,
            seed=config.seed,
            max_bases=config.tokenizer_max_bases,
            max_sequences=config.tokenizer_max_sequences,
            window_size=config.tokenizer_window_size,
        )
    tokenizer_window_lengths = [len(window) for window in tokenizer_windows]

    return {
        "sequence_count": len(sequences),
        "total_bases": sum(len(sequence) for sequence in sequences),
        "source_summary": dict(source_summary or {}),
        "tokenizer": {
            "vocab_size": len(tokenizer.vocab),
            "metadata": dict(tokenizer.metadata),
            "corpus_coverage": tokenizer_coverage,
            "window_length_summary": summarize_series(
                tokenizer_window_lengths,
                "tokenizer_window_length",
            ),
        },
        "mlm_windows": {
            "sample_count": effective_samples,
            "tokenized_tokens_total": total_tokenized,
            "maskable_tokens_total": total_maskable,
            "masked_tokens_total": total_masked,
            "unk_tokens_total": total_unk,
            "maskable_token_rate": total_maskable / total_tokenized if total_tokenized else 0.0,
            "masked_token_rate": total_masked / total_tokenized if total_tokenized else 0.0,
            "unk_token_rate": total_unk / total_tokenized if total_tokenized else 0.0,
            "windows_below_min_tokenized": sum(
                count < config.min_tokenized_tokens for count in tokenized_lengths
            ),
            "windows_below_min_maskable": sum(
                count < config.min_maskable_tokens for count in maskable_counts
            ),
            "windows_with_zero_masked_targets": sum(count == 0 for count in masked_counts),
            "tokenized_length_summary": summarize_series(
                tokenized_lengths,
                "tokenized_length",
            ),
            "maskable_count_summary": summarize_series(
                maskable_counts,
                "maskable_count",
            ),
            "masked_count_summary": summarize_series(
                masked_counts,
                "masked_count",
            ),
        },
    }


def write_json_report(path: Path, report: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
