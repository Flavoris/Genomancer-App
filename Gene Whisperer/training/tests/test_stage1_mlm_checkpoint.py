"""Tests for k-mer aware MLM checkpoint selection."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from train_stage1 import get_mlm_checkpoint_for_kmer


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("stub", encoding="utf-8")
    return path


def test_get_mlm_checkpoint_prefers_k_specific(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    custom_path = _touch(script_dir / "mlm_k3_custom.pt")
    _touch(tmp_path / "artifacts" / "mlm_encoder_k3.pt")

    cfg = {"mlm_encoder_ckpt_by_k": {"3": "mlm_k3_custom.pt"}}
    ckpt_path = get_mlm_checkpoint_for_kmer(cfg, 3, script_dir=script_dir)

    assert ckpt_path == custom_path.resolve()


def test_get_mlm_checkpoint_reads_nested_mapping(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    nested_path = _touch(script_dir / "nested_k4.pt")
    cfg = {"checkpoints": {"mlm_encoder_ckpt_by_k": {4: "nested_k4.pt"}}}

    ckpt_path = get_mlm_checkpoint_for_kmer(cfg, 4, script_dir=script_dir)

    assert ckpt_path == nested_path.resolve()


def test_get_mlm_checkpoint_uses_default_pattern(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    default_path = _touch(tmp_path / "artifacts" / "mlm_encoder_k5.pt")

    ckpt_path = get_mlm_checkpoint_for_kmer({}, 5, script_dir=script_dir)

    assert ckpt_path == default_path.resolve()


def test_get_mlm_checkpoint_falls_back_to_generic(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    generic_path = _touch(script_dir / "generic.pt")
    cfg = {"mlm_encoder_checkpoint": "generic.pt"}

    ckpt_path = get_mlm_checkpoint_for_kmer(cfg, 6, script_dir=script_dir)

    assert ckpt_path == generic_path.resolve()


def test_get_mlm_checkpoint_returns_none_when_missing(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    ckpt_path = get_mlm_checkpoint_for_kmer({}, 3, script_dir=script_dir)

    assert ckpt_path is None
