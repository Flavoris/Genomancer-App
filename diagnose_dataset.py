from __future__ import annotations

import argparse
import csv
import random
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable, Sequence

import yaml

REQUIRED_COLUMNS = ("sequence", "is_promoter")


@dataclass(frozen=True)
class DatasetStats:
    total_samples: int
    promoter_count: int
    non_promoter_count: int
    min_length: int
    max_length: int
    mean_length: float


def load_config(config_path: Path) -> dict:
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    return config


def resolve_dataset_path(config_path: Path, data_path: str) -> Path:
    candidate = Path(data_path)
    if candidate.is_absolute():
        return candidate
    return (config_path.parent / candidate).resolve()


def _parse_label(raw_label: str, row_index: int, path: Path) -> int:
    if raw_label is None:
        raise ValueError(f"Missing is_promoter value at row {row_index} in {path}")
    try:
        value = float(str(raw_label).strip())
    except ValueError as exc:
        raise ValueError(
            f"Invalid is_promoter value '{raw_label}' at row {row_index} in {path}"
        ) from exc
    if value not in (0.0, 1.0):
        raise ValueError(f"Unexpected is_promoter value '{raw_label}' at row {row_index} in {path}")
    return int(value)


def load_tsv(path: Path, delimiter: str) -> tuple[list[str], list[int]]:
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")
    sequences: list[str] = []
    labels: list[int] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
        if missing:
            raise ValueError(f"Missing columns in {path}: {missing}")
        for row_index, row in enumerate(reader, start=2):
            sequence = str(row.get("sequence", "")).strip()
            label = _parse_label(row.get("is_promoter"), row_index, path)
            sequences.append(sequence)
            labels.append(label)
    return sequences, labels


def compute_dataset_stats(sequences: Sequence[str], labels: Sequence[int]) -> DatasetStats:
    if len(sequences) != len(labels):
        raise ValueError("Sequences and labels must have the same length")

    promoter_count = sum(1 for label in labels if label == 1)
    non_promoter_count = sum(1 for label in labels if label == 0)
    unknown_count = len(labels) - promoter_count - non_promoter_count
    if unknown_count:
        raise ValueError(f"Unexpected label values encountered: {unknown_count}")

    lengths = [len(sequence) for sequence in sequences]
    if lengths:
        min_length = min(lengths)
        max_length = max(lengths)
        mean_length = sum(lengths) / len(lengths)
    else:
        min_length = 0
        max_length = 0
        mean_length = 0.0

    return DatasetStats(
        total_samples=len(sequences),
        promoter_count=promoter_count,
        non_promoter_count=non_promoter_count,
        min_length=min_length,
        max_length=max_length,
        mean_length=mean_length,
    )


def find_sequence_overlap(train_sequences: Iterable[str], val_sequences: Iterable[str]) -> set[str]:
    train_set = set(train_sequences)
    val_set = set(val_sequences)
    return train_set.intersection(val_set)


def similarity(a: str, b: str) -> float:
    return SequenceMatcher(None, a, b).ratio()


def sample_sequences(
    sequences: Sequence[str],
    sample_size: int,
    rng: random.Random,
) -> list[str]:
    if not sequences or sample_size <= 0:
        return []
    if sample_size >= len(sequences):
        return list(sequences)
    return rng.sample(list(sequences), sample_size)


def count_high_similarity_pairs(sequences: Sequence[str], threshold: float) -> int:
    high_similarity_pairs = 0
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            if similarity(sequences[i], sequences[j]) > threshold:
                high_similarity_pairs += 1
    return high_similarity_pairs


def format_ratio(train_total: int, val_total: int) -> str:
    combined = train_total + val_total
    if combined == 0:
        return "0.00 (train 0.0%, val 0.0%)"
    train_ratio = train_total / combined
    return f"{train_ratio:.2f} (train {train_ratio * 100:.1f}%, val {(1 - train_ratio) * 100:.1f}%)"


def print_dataset_report(
    train_stats: DatasetStats,
    val_stats: DatasetStats,
    overlap_count: int,
    high_similarity_pairs: int,
) -> None:
    print("Dataset diagnostics")
    print("===================")
    print(f"Train samples: {train_stats.total_samples}")
    print(f"Validation samples: {val_stats.total_samples}")
    print(f"Train/val split ratio: {format_ratio(train_stats.total_samples, val_stats.total_samples)}")
    print()
    print("Class balance")
    print("-------------")
    print(
        "Train promoters: {promoters} | non-promoters: {non_promoters}".format(
            promoters=train_stats.promoter_count,
            non_promoters=train_stats.non_promoter_count,
        )
    )
    print(
        "Validation promoters: {promoters} | non-promoters: {non_promoters}".format(
            promoters=val_stats.promoter_count,
            non_promoters=val_stats.non_promoter_count,
        )
    )
    print()
    print("Sequence length distribution (bp)")
    print("---------------------------------")
    print(
        "Train min: {min_len} | max: {max_len} | mean: {mean_len:.2f}".format(
            min_len=train_stats.min_length,
            max_len=train_stats.max_length,
            mean_len=train_stats.mean_length,
        )
    )
    print(
        "Validation min: {min_len} | max: {max_len} | mean: {mean_len:.2f}".format(
            min_len=val_stats.min_length,
            max_len=val_stats.max_length,
            mean_len=val_stats.mean_length,
        )
    )
    print()
    print(f"Sequences in both train and val: {overlap_count}")
    if overlap_count > 0:
        print("WARNING: Data leakage detected!")
    print(f"High similarity pairs in sample: {high_similarity_pairs}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnose stage1 dataset splits.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("gene_whisperer/configs/finetune.yaml"),
        help="Path to config.yaml containing stage1 dataset locations.",
    )
    parser.add_argument(
        "--train",
        type=Path,
        default=None,
        help="Override stage1 training TSV path.",
    )
    parser.add_argument(
        "--val",
        type=Path,
        default=None,
        help="Override stage1 validation TSV path.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=500,
        help="Number of sequences to sample for similarity checks.",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.85,
        help="Similarity threshold for counting near-duplicates.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1337,
        help="Random seed for sampling sequences.",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    config_path = args.config
    config = load_config(config_path)
    delimiter = config.get("delimiter", "\t")

    train_path = args.train
    if train_path is None:
        stage1_train = config.get("stage1_train")
        if not stage1_train:
            raise ValueError("stage1_train is not configured in the config file")
        train_path = resolve_dataset_path(config_path, stage1_train)

    val_path = args.val
    if val_path is None:
        stage1_val = config.get("stage1_val")
        if not stage1_val:
            raise ValueError("stage1_val is not configured in the config file")
        val_path = resolve_dataset_path(config_path, stage1_val)

    train_sequences, train_labels = load_tsv(train_path, delimiter)
    val_sequences, val_labels = load_tsv(val_path, delimiter)

    train_stats = compute_dataset_stats(train_sequences, train_labels)
    val_stats = compute_dataset_stats(val_sequences, val_labels)

    overlap = find_sequence_overlap(train_sequences, val_sequences)

    rng = random.Random(args.seed)
    sample = sample_sequences(train_sequences, args.sample_size, rng)
    high_similarity_pairs = count_high_similarity_pairs(sample, args.similarity_threshold)

    print_dataset_report(
        train_stats=train_stats,
        val_stats=val_stats,
        overlap_count=len(overlap),
        high_similarity_pairs=high_similarity_pairs,
    )


if __name__ == "__main__":
    main()
