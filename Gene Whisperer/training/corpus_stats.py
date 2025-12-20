#!/usr/bin/env python3
"""
Compute 3-mer coverage statistics over a DNA corpus.

Samples N windows from either:
  1) A genome FASTA (sliding windows with stride), or
  2) A line-delimited corpus file (optionally TSV with a leading "sequence" header).

For the sampled windows, it:
  - counts occurrences for each of the 64 real 3-mers (AAA..TTT)
  - prints per-3-mer frequencies
  - reports min/max and flags missing 3-mers
"""

from __future__ import annotations

import argparse
import random
import sys
from itertools import product
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None

BASES: Tuple[str, ...] = ("A", "C", "G", "T")
BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}
TRI_KMER_VOCAB: List[str] = ["".join(kmer) for kmer in product(BASES, repeat=3)]

ALLOWED_BASES = set(BASES)
REV_COMP_MAP = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of an A/C/G/T-only DNA sequence."""
    return sequence.translate(REV_COMP_MAP)[::-1]


def read_fasta_sequence(path: Path) -> str:
    sequence_parts: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)


def _resolve_cfg_path(path_value: str | Path | None, *, base_dir: Path) -> Path:
    if path_value is None or str(path_value).strip() == "":
        raise ValueError("Expected a non-empty path value")
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def iter_fasta_windows(
    fasta_path: Path,
    window_size: int,
    stride: int,
    *,
    include_reverse_complements: bool,
) -> Iterator[str]:
    sequence = read_fasta_sequence(fasta_path)
    if len(sequence) < window_size:
        return
    max_start = len(sequence) - window_size + 1
    for start in range(0, max_start, stride):
        window = sequence[start : start + window_size]
        if not set(window) <= ALLOWED_BASES:
            continue
        yield window
        if include_reverse_complements:
            yield reverse_complement(window)


def iter_corpus_windows(
    corpus_path: Path,
    *,
    delimiter: str,
    include_reverse_complements: bool,
) -> Iterator[str]:
    with corpus_path.open("r", encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            first_field = line.split(delimiter, 1)[0].strip()
            if line_num == 1 and first_field.lower() == "sequence":
                continue
            window = first_field.upper()
            if not window or not set(window) <= ALLOWED_BASES:
                continue
            yield window
            if include_reverse_complements:
                yield reverse_complement(window)


def reservoir_sample(
    iterable: Iterable[str],
    sample_size: int,
    rng: random.Random,
) -> Tuple[List[str], int]:
    """
    Uniformly sample up to sample_size items from iterable without storing all items.

    Returns: (sampled_items, total_items_seen)
    """
    if sample_size <= 0:
        raise ValueError("sample_size must be > 0 for reservoir sampling")

    reservoir: List[str] = []
    total_seen = 0
    for item in iterable:
        total_seen += 1
        if len(reservoir) < sample_size:
            reservoir.append(item)
            continue
        j = rng.randrange(total_seen)
        if j < sample_size:
            reservoir[j] = item
    return reservoir, total_seen


def count_3mers(windows: Iterable[str]) -> Tuple[List[int], int, int]:
    """
    Count all 3-mers (AAA..TTT) across windows.

    Returns: (counts[64], total_tokens, windows_used)
    """
    counts = [0] * 64
    total_tokens = 0
    windows_used = 0
    for window in windows:
        if len(window) < 3:
            continue
        windows_used += 1
        for i in range(len(window) - 2):
            a = BASE_TO_INT[window[i]]
            b = BASE_TO_INT[window[i + 1]]
            c = BASE_TO_INT[window[i + 2]]
            idx = (a << 4) + (b << 2) + c
            counts[idx] += 1
            total_tokens += 1
    return counts, total_tokens, windows_used


def print_report(
    *,
    counts: List[int],
    total_tokens: int,
    windows_used: int,
    total_windows_seen: Optional[int],
    source_label: str,
) -> List[str]:
    if len(counts) != 64:
        raise ValueError(f"Expected 64 3-mer counts, got {len(counts)}")

    print("=== Corpus 3-mer Coverage ===")
    print(f"Source: {source_label}")
    if total_windows_seen is not None:
        print(f"Windows seen: {total_windows_seen}")
    print(f"Windows used: {windows_used}")
    print(f"3-mer tokens counted: {total_tokens}")
    print()
    print(f"{'3-mer':<5} {'count':>12} {'freq':>12}")
    for idx, kmer in enumerate(TRI_KMER_VOCAB):
        count = counts[idx]
        freq = (count / total_tokens) if total_tokens > 0 else 0.0
        print(f"{kmer:<5} {count:>12} {freq:>11.6%}")

    min_count = min(counts) if counts else 0
    max_count = max(counts) if counts else 0
    min_kmers = [TRI_KMER_VOCAB[i] for i, c in enumerate(counts) if c == min_count]
    max_kmers = [TRI_KMER_VOCAB[i] for i, c in enumerate(counts) if c == max_count]
    missing = [TRI_KMER_VOCAB[i] for i, c in enumerate(counts) if c == 0]

    sorted_idx = sorted(range(64), key=lambda i: (counts[i], TRI_KMER_VOCAB[i]))
    bottom = [TRI_KMER_VOCAB[i] for i in sorted_idx[:8]]
    top = [TRI_KMER_VOCAB[i] for i in reversed(sorted_idx[-8:])]

    print()
    print("=== Summary ===")
    print(f"Min count: {min_count} ({', '.join(min_kmers)})")
    print(f"Max count: {max_count} ({', '.join(max_kmers)})")
    print(f"Bottom-8:  {', '.join(bottom)}")
    print(f"Top-8:     {', '.join(top)}")
    if missing:
        print(f"MISSING:   {', '.join(missing)}")
    else:
        print("MISSING:   (none)")

    return missing


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify 3-mer coverage over a corpus.")
    parser.add_argument("--config", type=str, default=None, help="Optional YAML config (uses mlm_* keys).")
    parser.add_argument("--fasta", type=str, default=None, help="Genome FASTA path (overrides config).")
    parser.add_argument("--corpus", type=str, default=None, help="Line-delimited corpus/TSV path.")
    parser.add_argument("--delimiter", type=str, default="\t", help="Delimiter for --corpus files (default: tab).")
    parser.add_argument("--window-size", type=int, default=None, help="Window size for FASTA mode.")
    parser.add_argument("--stride", type=int, default=None, help="Stride for FASTA mode.")
    parser.add_argument("--num-windows", type=int, default=10_000, help="Sample size N (0 => use all).")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for sampling (default: config seed or 1337).")
    parser.add_argument(
        "--include-reverse-complements",
        action="store_true",
        help="Include reverse-complement windows (overrides config default).",
    )
    parser.add_argument(
        "--fail-on-missing",
        action="store_true",
        help="Exit with code 1 if any 3-mer never appears in the sample.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent

    cfg = {}
    if args.config is not None:
        if yaml is None:
            raise RuntimeError("PyYAML is required for --config")
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = (script_dir / config_path).resolve()
        with config_path.open("r", encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle) or {}

    seed = int(args.seed if args.seed is not None else cfg.get("seed", 1337))

    include_rc_default = bool(
        cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False))
    )
    include_rc = bool(args.include_reverse_complements or include_rc_default)

    num_windows = int(args.num_windows)

    if args.corpus and args.fasta:
        print("Error: Provide only one of --corpus or --fasta.", file=sys.stderr)
        return 2

    if args.corpus:
        corpus_path = Path(args.corpus)
        if not corpus_path.is_absolute():
            corpus_path = (script_dir / corpus_path).resolve()
        if not corpus_path.exists():
            print(f"Error: corpus not found: {corpus_path}", file=sys.stderr)
            return 2

        windows_iter = iter_corpus_windows(
            corpus_path,
            delimiter=args.delimiter,
            include_reverse_complements=include_rc,
        )
        source_label = str(corpus_path)
    else:
        fasta_value = args.fasta if args.fasta is not None else cfg.get("mlm_fasta_path")
        if fasta_value is None:
            print("Error: Provide --fasta or --corpus (or set mlm_fasta_path in --config).", file=sys.stderr)
            return 2
        fasta_path = _resolve_cfg_path(fasta_value, base_dir=script_dir)
        if not fasta_path.exists():
            print(f"Error: FASTA not found: {fasta_path}", file=sys.stderr)
            return 2

        window_size = int(args.window_size if args.window_size is not None else cfg.get("mlm_window_size", 81))
        stride = int(args.stride if args.stride is not None else cfg.get("mlm_stride", 20))

        windows_iter = iter_fasta_windows(
            fasta_path,
            window_size,
            stride,
            include_reverse_complements=include_rc,
        )
        source_label = f"{fasta_path} (window_size={window_size}, stride={stride}, rc={include_rc})"

    rng = random.Random(seed)

    if num_windows <= 0:
        counts, total_tokens, windows_used = count_3mers(windows_iter)
        missing = print_report(
            counts=counts,
            total_tokens=total_tokens,
            windows_used=windows_used,
            total_windows_seen=None,
            source_label=source_label,
        )
    else:
        sampled, total_seen = reservoir_sample(windows_iter, num_windows, rng)
        counts, total_tokens, windows_used = count_3mers(sampled)
        missing = print_report(
            counts=counts,
            total_tokens=total_tokens,
            windows_used=windows_used,
            total_windows_seen=total_seen,
            source_label=source_label,
        )

    if args.fail_on_missing and missing:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
