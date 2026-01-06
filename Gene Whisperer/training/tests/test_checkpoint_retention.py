"""Smoke tests for MLM checkpoint retention."""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from checkpoint_utils import prune_step_checkpoints


def _touch(path: Path) -> None:
    """Create a dummy checkpoint file."""
    path.write_text("checkpoint", encoding="utf-8")


def test_prune_keeps_latest_checkpoints() -> None:
    """Keep only the newest step checkpoints."""
    with tempfile.TemporaryDirectory() as tmpdir:
        checkpoints_dir = Path(tmpdir)
        for step in (10, 20, 30):
            _touch(checkpoints_dir / f"checkpoint_step_{step}.pt")

        removed = prune_step_checkpoints(checkpoints_dir, keep_last=2)

        remaining = {path.name for path in checkpoints_dir.glob("checkpoint_step_*.pt")}
        assert remaining == {"checkpoint_step_20.pt", "checkpoint_step_30.pt"}
        assert {path.name for path in removed} == {"checkpoint_step_10.pt"}

    print("test_prune_keeps_latest_checkpoints: PASSED")


def test_prune_preserves_keep_paths() -> None:
    """Keep explicit checkpoints even when older."""
    with tempfile.TemporaryDirectory() as tmpdir:
        checkpoints_dir = Path(tmpdir)
        keep_path = checkpoints_dir / "checkpoint_step_5.pt"
        _touch(keep_path)
        _touch(checkpoints_dir / "checkpoint_step_10.pt")
        _touch(checkpoints_dir / "checkpoint_step_15.pt")

        removed = prune_step_checkpoints(
            checkpoints_dir,
            keep_last=1,
            keep_paths=[keep_path],
        )

        remaining = {path.name for path in checkpoints_dir.glob("checkpoint_step_*.pt")}
        assert remaining == {"checkpoint_step_5.pt", "checkpoint_step_15.pt"}
        assert {path.name for path in removed} == {"checkpoint_step_10.pt"}

    print("test_prune_preserves_keep_paths: PASSED")
