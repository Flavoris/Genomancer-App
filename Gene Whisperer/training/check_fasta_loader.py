#!/usr/bin/env python3
"""Quick FASTA loader sanity check for the bundled genomes."""
from __future__ import annotations

import random
from pathlib import Path

from genome_io import read_fasta_records, sanitize_unknown_bases

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"


def _summarize(path: Path, strategy: str) -> int:
    records = read_fasta_records(path)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    total_bp = sum(len(sanitize_unknown_bases(seq, strategy)) for _, seq in records)
    top_names = [record_id for record_id, _ in records[:3]]
    print(f"{path.name}: records={len(records)} top3={top_names} total_bp={total_bp}")
    return len(records)


def main() -> None:
    random.seed(1337)
    genomes = [
        DATA_DIR / "GCF_000005845.2_ASM584v2_genomic.fna",
        DATA_DIR / "B_subtilis_genome",
        DATA_DIR / "S_cerevisiae_genome",
    ]

    for path in genomes:
        if not path.exists():
            raise FileNotFoundError(f"Expected genome file not found: {path}")

    strategy = "random"
    _summarize(genomes[0], strategy)
    _summarize(genomes[1], strategy)
    yeast_records = _summarize(genomes[2], strategy)

    print(f"S_cerevisiae multi-record: {yeast_records > 1}")
    if yeast_records <= 1:
        raise AssertionError("Expected >1 FASTA record for S. cerevisiae")

if __name__ == "__main__":
    main()
