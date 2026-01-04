"""FASTA ingestion and unknown-base sanitization utilities."""
from __future__ import annotations

import random
from pathlib import Path
from typing import List, Tuple

ALLOWED_BASES = {"A", "C", "G", "T"}
BASE_CHOICES = ("A", "C", "G", "T")


def _normalize_sequence_line(line: str) -> str:
    """Uppercase and map any non-ACGT character to N."""
    line = line.strip().upper()
    if not line:
        return ""
    return "".join(base if base in ALLOWED_BASES else "N" for base in line)


def read_fasta_records(path: Path, *, keep_header: bool = False) -> List[Tuple[str, str]]:
    """
    Read a FASTA file into ordered (record_id, sequence) tuples.

    - Treats lines starting with '>' as headers
    - Concatenates multi-line sequences per record
    - Uppercases and maps non-ACGT to 'N'
    - When keep_header=True, use the full header line as record_id
    """
    records: List[Tuple[str, str]] = []
    current_id: str | None = None
    current_parts: List[str] = []
    fallback_index = 0

    def flush() -> None:
        if current_id is not None:
            records.append((current_id, "".join(current_parts)))

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if header:
                    current_id = header if keep_header else header.split()[0]
                else:
                    current_id = f"seq_{fallback_index}"
                    fallback_index += 1
                current_parts = []
                continue
            if current_id is None:
                current_id = f"seq_{fallback_index}"
                fallback_index += 1
                current_parts = []
            current_parts.append(_normalize_sequence_line(line))

    flush()
    return records


def find_acgt_segments(sequence: str, min_length: int) -> List[Tuple[int, int]]:
    """Return (start, end) segments of A/C/G/T runs with length >= min_length."""
    if min_length <= 0:
        raise ValueError(f"min_length must be > 0, got {min_length}")
    segments: List[Tuple[int, int]] = []
    start: int | None = None

    for idx, base in enumerate(sequence):
        if base in ALLOWED_BASES:
            if start is None:
                start = idx
        else:
            if start is not None:
                if idx - start >= min_length:
                    segments.append((start, idx))
                start = None

    if start is not None and len(sequence) - start >= min_length:
        segments.append((start, len(sequence)))

    return segments


def sanitize_unknown_bases(
    sequence: str,
    strategy: str = "random",
    *,
    rng: random.Random | None = None,
) -> str:
    """
    Replace any non-ACGT base with a sanitized value.

    Currently supported strategies:
      - "random": replace each unknown base with random A/C/G/T
    """
    if rng is None:
        rng = random
    strategy = (strategy or "").strip().lower()
    if strategy != "random":
        raise ValueError(f"Unknown unknown-base strategy: {strategy!r}")
    sequence = sequence.upper()
    return "".join(base if base in ALLOWED_BASES else rng.choice(BASE_CHOICES) for base in sequence)
