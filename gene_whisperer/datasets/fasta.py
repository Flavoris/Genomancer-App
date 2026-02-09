"""FASTA loading utilities."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List


def _normalize(seq: str) -> str:
    return "".join(ch for ch in seq.upper().strip() if ch in "ACGTN")


def _iter_fasta_file(
    path: Path,
    max_bases_per_file: int | None = None,
) -> Iterator[str]:
    buffer: List[str] = []
    consumed_bases = 0

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buffer:
                    yield _normalize("".join(buffer))
                    buffer = []
                continue

            if max_bases_per_file is not None:
                remaining = max_bases_per_file - consumed_bases
                if remaining <= 0:
                    break
                if len(line) > remaining:
                    line = line[:remaining]

            consumed_bases += len(line)
            buffer.append(line)

    if buffer:
        yield _normalize("".join(buffer))


def iter_fasta_sequences(
    paths: Iterable[Path],
    max_bases_per_file: int | None = None,
    verbose: bool = False,
) -> Iterator[str]:
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"FASTA path not found: {path}")

        if path.is_dir():
            files = [child for child in sorted(path.iterdir()) if child.is_file()]
        else:
            files = [path]

        for file_path in files:
            sequence_count = 0
            base_count = 0
            if verbose:
                print(f"Loading FASTA: {file_path}", flush=True)

            for sequence in _iter_fasta_file(
                file_path,
                max_bases_per_file=max_bases_per_file,
            ):
                sequence_count += 1
                base_count += len(sequence)
                yield sequence

            if verbose:
                print(
                    f"Loaded FASTA: {file_path} sequences={sequence_count} bases={base_count}",
                    flush=True,
                )


def load_fasta_sequences(
    paths: Iterable[Path],
    max_bases_per_file: int | None = None,
    verbose: bool = False,
) -> List[str]:
    return list(
        iter_fasta_sequences(
            paths,
            max_bases_per_file=max_bases_per_file,
            verbose=verbose,
        )
    )
