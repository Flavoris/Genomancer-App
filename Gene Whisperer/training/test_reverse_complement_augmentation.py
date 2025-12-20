from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from build_pretrain_corpus import build_pretrain_corpus
from pretrain_mlm import prepare_dataset


def _write_fasta(path: Path, sequence: str) -> None:
    path.write_text(f">test\n{sequence}\n", encoding="utf-8")


def test_build_pretrain_corpus_include_reverse_complements_doubles(tmp_path: Path) -> None:
    fasta_path = tmp_path / "toy.fna"
    # Two windows: AAAA, CCCC (seq_len=4, stride=4)
    _write_fasta(fasta_path, "AAAACCCC")

    out_path = tmp_path / "corpus.txt"
    stats = build_pretrain_corpus(
        fasta_path=fasta_path,
        output_path=out_path,
        seq_len=4,
        num_samples=None,
        seed=123,
        sampling="sliding",
        stride=4,
        include_reverse_complements=True,
        shuffle=False,
    )

    lines = out_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4
    assert lines == ["AAAA", "CCCC", "TTTT", "GGGG"]
    assert stats.reverse_complement_samples == 2
    assert stats.reverse_complement_fraction == pytest.approx(0.5, abs=1e-9)


def test_prepare_dataset_include_reverse_complements_doubles(tmp_path: Path) -> None:
    fasta_path = tmp_path / "toy.fna"
    _write_fasta(fasta_path, "AAAACCCC")

    vocab_path = tmp_path / "vocab.json"
    cfg = {
        "mlm_fasta_path": str(fasta_path),
        "mlm_window_size": 4,
        "mlm_stride": 4,
        "mlm_kmer": 3,
        "mlm_vocab_path": str(vocab_path),
        "include_reverse_complements": True,
    }

    token_tensors, _ = prepare_dataset(cfg)
    assert len(token_tensors) == 4
