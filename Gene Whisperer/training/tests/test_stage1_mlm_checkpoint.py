"""Tests for BPE MLM checkpoint selection."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from train_stage1 import get_mlm_checkpoint


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("stub", encoding="utf-8")
    return path


def test_get_mlm_checkpoint_uses_explicit_config(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    explicit_path = _touch(script_dir / "custom_mlm.pt")
    cfg = {"mlm_encoder_checkpoint": "custom_mlm.pt"}

    ckpt_path = get_mlm_checkpoint(cfg, script_dir=script_dir)

    assert ckpt_path == explicit_path.resolve()


def test_get_mlm_checkpoint_uses_mlm_encoder_path(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    path = _touch(script_dir / "encoder.pt")
    cfg = {"mlm_encoder_path": "encoder.pt"}

    ckpt_path = get_mlm_checkpoint(cfg, script_dir=script_dir)

    assert ckpt_path == path.resolve()


def test_get_mlm_checkpoint_uses_default_bpe_path(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    default_path = _touch(tmp_path / "artifacts" / "mlm_encoder_bpe.pt")

    ckpt_path = get_mlm_checkpoint({}, script_dir=script_dir)

    assert ckpt_path == default_path.resolve()


def test_get_mlm_checkpoint_returns_none_when_missing(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    ckpt_path = get_mlm_checkpoint({}, script_dir=script_dir)

    assert ckpt_path is None


def test_get_mlm_checkpoint_prefers_explicit_over_default(tmp_path: Path) -> None:
    script_dir = tmp_path / "training"
    script_dir.mkdir()

    explicit_path = _touch(script_dir / "my_encoder.pt")
    _touch(tmp_path / "artifacts" / "mlm_encoder_bpe.pt")

    cfg = {"mlm_encoder_checkpoint": "my_encoder.pt"}
    ckpt_path = get_mlm_checkpoint(cfg, script_dir=script_dir)

    assert ckpt_path == explicit_path.resolve()
