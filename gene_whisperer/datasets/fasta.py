"""FASTA loading utilities."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List


def _normalize(seq: str) -> str:
    return "".join(ch for ch in seq.upper().strip() if ch in "ACGTN")


def _iter_fasta_file(path: Path) -> Iterator[str]:
    buffer: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buffer:
                    yield _normalize("".join(buffer))
                    buffer = []
                continue
            buffer.append(line)
    if buffer:
        yield _normalize("".join(buffer))


def iter_fasta_sequences(paths: Iterable[Path]) -> Iterator[str]:
    for path in paths:
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if child.is_file():
                    yield from _iter_fasta_file(child)
        else:
            yield from _iter_fasta_file(path)


def load_fasta_sequences(paths: Iterable[Path]) -> List[str]:
    return list(iter_fasta_sequences(paths))
