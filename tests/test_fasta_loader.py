from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.datasets.fasta import load_fasta_sequences


def test_load_fasta_sequences_with_base_cap(tmp_path: Path) -> None:
    fasta_path = tmp_path / "genome.fna"
    fasta_path.write_text(
        ">seq1\n" + "A" * 10 + "\n" + "C" * 10 + "\n",
        encoding="utf-8",
    )

    full = load_fasta_sequences([fasta_path])
    capped = load_fasta_sequences([fasta_path], max_bases_per_file=12)

    assert full == ["A" * 10 + "C" * 10]
    assert capped == ["A" * 10 + "C" * 2]


def test_load_fasta_sequences_missing_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fna"
    try:
        load_fasta_sequences([missing])
    except FileNotFoundError as exc:
        assert "FASTA path not found" in str(exc)
    else:
        raise AssertionError("Expected FileNotFoundError for missing FASTA path")
