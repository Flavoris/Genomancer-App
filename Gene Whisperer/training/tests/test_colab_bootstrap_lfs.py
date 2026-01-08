"""Sanity checks for Colab bootstrap LFS handling."""
from __future__ import annotations

from pathlib import Path


def _repo_root() -> Path:
    """Resolve repository root from this test location."""
    return Path(__file__).resolve().parents[3]


def test_lfs_filters_disabled() -> None:
    """Ensure LFS filters are not enforced in .gitattributes."""
    gitattributes = (_repo_root() / ".gitattributes").read_text(encoding="utf-8")
    assert "filter=lfs" not in gitattributes


def test_bootstrap_supports_skip_lfs_pull() -> None:
    """Ensure the bootstrap script supports skipping LFS pulls."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_bootstrap.sh"
    content = script_path.read_text(encoding="utf-8")
    assert "SKIP_LFS_PULL" in content
    assert "git lfs pull" in content


def test_bootstrap_repairs_reference_genomes() -> None:
    """Ensure the bootstrap script can repair small reference genomes."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_bootstrap.sh"
    content = script_path.read_text(encoding="utf-8")
    assert "ensure_reference_genomes.py" in content
    assert "SKIP_REFERENCE_GENOME_FIX" in content
