from __future__ import annotations

import csv
import random
from pathlib import Path
import sys

import pytest

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from diagnose_dataset import (
    compute_dataset_stats,
    count_high_similarity_pairs,
    find_sequence_overlap,
    load_tsv,
    sample_sequences,
)


def _write_tsv(path: Path, rows: list[tuple[str, int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "is_promoter"])
        for sequence, label in rows:
            writer.writerow([sequence, label])


def test_load_tsv_and_stats(tmp_path: Path) -> None:
    dataset_path = tmp_path / "stage1_train.tsv"
    _write_tsv(
        dataset_path,
        [
            ("AAAA", 1),
            ("AAAAT", 0),
            ("GGGG", 1),
        ],
    )

    sequences, labels = load_tsv(dataset_path, delimiter="\t")
    stats = compute_dataset_stats(sequences, labels)

    assert stats.total_samples == 3
    assert stats.promoter_count == 2
    assert stats.non_promoter_count == 1
    assert stats.min_length == 4
    assert stats.max_length == 5
    assert stats.mean_length == pytest.approx(4.33, rel=1e-2)


def test_find_sequence_overlap() -> None:
    train_sequences = ["AAAA", "CCCC", "GGGG"]
    val_sequences = ["TTTT", "CCCC"]
    overlap = find_sequence_overlap(train_sequences, val_sequences)
    assert overlap == {"CCCC"}


def test_similarity_count_is_deterministic() -> None:
    sequences = ["AAAAAA", "AAAAAA", "CCCCCC"]
    sample = sample_sequences(sequences, sample_size=len(sequences), rng=random.Random(0))
    count = count_high_similarity_pairs(sample, threshold=0.85)
    assert count == 1
