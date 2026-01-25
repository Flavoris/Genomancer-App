"""Tests for the Colab MLM runner script."""
from __future__ import annotations

from pathlib import Path


def _repo_root() -> Path:
    """Resolve repository root from this test location."""
    return Path(__file__).resolve().parents[3]


def test_colab_run_mlm_script_exists() -> None:
    """Ensure the colab MLM runner script exists."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_run_mlm.sh"
    assert script_path.exists(), f"Script not found: {script_path}"


def test_colab_run_mlm_has_training_command() -> None:
    """Ensure the script contains a training command invocation."""
    script_path = _repo_root() / "Gene Whisperer" / "scripts" / "colab_run_mlm.sh"
    content = script_path.read_text(encoding="utf-8")

    # Script should reference the training directory and python
    assert "python" in content
    assert "pretrain" in content.lower() or "train" in content.lower()
