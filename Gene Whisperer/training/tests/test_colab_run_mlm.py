"""Tests for the Colab MLM runner script."""
from __future__ import annotations

from pathlib import Path


def _repo_root() -> Path:
    """Resolve repository root from this test location."""
    return Path(__file__).resolve().parents[3]


def test_colab_run_mlm_supports_auto_multi_kmer() -> None:
    """Ensure the script supports auto multi-kmer selection."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_run_mlm.sh"
    content = script_path.read_text(encoding="utf-8")

    assert "--single_kmer" in content
    assert "--multi_kmer" in content
    assert "multi_kmer_pretrain_enabled" in content


def test_colab_run_mlm_builds_missing_vocabs() -> None:
    """Ensure the script can build missing vocabularies."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_run_mlm.sh"
    content = script_path.read_text(encoding="utf-8")

    assert "build_multi_kmer_vocabs.py" in content
    assert "skip_vocab_build" in content
