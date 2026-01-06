"""Checkpoint pruning helpers for training workflows."""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple

_STEP_CHECKPOINT_RE = re.compile(r"checkpoint_step_(\d+)\.pt$")


def _parse_checkpoint_step(path: Path) -> Optional[int]:
    """Extract the step number from a checkpoint filename."""
    match = _STEP_CHECKPOINT_RE.match(path.name)
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def prune_step_checkpoints(
    checkpoints_dir: Path,
    keep_last: int,
    *,
    keep_paths: Optional[Iterable[Path]] = None,
    logger: Optional[logging.Logger] = None,
) -> List[Path]:
    """
    Remove older step checkpoints while keeping the most recent ones.

    Args:
        checkpoints_dir: Directory containing checkpoint_step_*.pt files.
        keep_last: Number of most recent step checkpoints to keep (min 1).
        keep_paths: Additional checkpoint paths to preserve.
        logger: Optional logger for warnings.

    Returns:
        List of paths that were removed.
    """
    if keep_last < 1:
        keep_last = 1

    step_checkpoints: List[Tuple[int, Path]] = []
    for path in checkpoints_dir.glob("checkpoint_step_*.pt"):
        step = _parse_checkpoint_step(path)
        if step is not None:
            step_checkpoints.append((step, path))

    if len(step_checkpoints) <= keep_last:
        return []

    step_checkpoints.sort(key=lambda item: item[0])
    keep_set: Set[Path] = {
        path.resolve() for _, path in step_checkpoints[-keep_last:]
    }
    if keep_paths:
        keep_set.update({Path(path).resolve() for path in keep_paths})

    removed: List[Path] = []
    for _, path in step_checkpoints[:-keep_last]:
        if path.resolve() in keep_set:
            continue
        try:
            path.unlink()
            removed.append(path)
        except FileNotFoundError:
            continue
        except OSError as exc:
            if logger is not None:
                logger.warning("Failed to remove old checkpoint %s: %s", path, exc)

    return removed
