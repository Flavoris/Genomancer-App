"""Tests for BPE MLM padding behavior in pretraining utilities."""
from __future__ import annotations

from pathlib import Path

from bpe_tokenizer import DNABPETokenizer
from pretrain_mlm import prepare_dataset


def _write_fasta(path: Path, sequence: str) -> None:
    path.write_text(f">seq1\n{sequence}\n", encoding="utf-8")


def _build_vocab(path: Path, sequence: str) -> None:
    tokenizer = DNABPETokenizer(vocab_size=64)
    # Small corpus is fine for unit tests; merges will stop early if frequency is low.
    tokenizer.train([sequence], min_frequency=2, verbose=False)
    tokenizer.save(str(path))


def test_prepare_dataset_uses_mlm_max_token_len(tmp_path: Path) -> None:
    fasta_path = tmp_path / "test.fna"
    vocab_path = tmp_path / "bpe_vocab.json"

    sequence = "ATCG" * 25  # 100 bp
    _write_fasta(fasta_path, sequence)
    _build_vocab(vocab_path, sequence)

    cfg = {
        "mlm_fasta_path": str(fasta_path),
        "mlm_window_size": len(sequence),
        "mlm_stride": len(sequence),
        "mlm_unknown_base_strategy": "random",
        "bpe_vocab_path": str(vocab_path),
        "mlm_max_token_len": 16,
    }

    token_tensors, _vocab = prepare_dataset(cfg)
    assert token_tensors, "Expected at least one tokenized window"
    assert token_tensors[0].numel() == 16, "BPE windows should pad to mlm_max_token_len"
