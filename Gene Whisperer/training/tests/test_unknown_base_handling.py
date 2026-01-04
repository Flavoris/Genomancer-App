"""Tests for unknown-base handling in genome preprocessing."""
from __future__ import annotations

from genome_io import find_acgt_segments
from pretrain_mlm import extract_windows


def test_find_acgt_segments_breaks_on_unknowns() -> None:
    sequence = "AAAACNNNNGGTTNCC"
    segments = find_acgt_segments(sequence, min_length=2)
    assert segments == [(0, 5), (9, 13), (14, 16)]


def test_extract_windows_skip_strategy() -> None:
    sequence = "ACGTNNACGTAC"
    windows = extract_windows(sequence, window_size=4, stride=1, unknown_base_strategy="skip")
    assert windows == ["ACGT", "ACGT", "CGTA", "GTAC"]
    assert all("N" not in window for window in windows)
