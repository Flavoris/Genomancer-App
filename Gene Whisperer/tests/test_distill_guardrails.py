"""
Unit tests for distillation guardrails.

This module tests the `apply_distill_guardrails` function from train_stage1.py
to ensure distillation is only enabled in safe scenarios.
"""

import sys
from pathlib import Path

# Add training directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "training"))

from train_stage1 import apply_distill_guardrails


def test_disabled_returns_false():
    """When distillation is disabled, guardrail should return False."""
    enabled, reason = apply_distill_guardrails(
        distill_enabled=False,
        teacher_kmer=6,
        student_kmer=4,
        teacher_ckpt_path=Path("/tmp/teacher.pt"),
        student_ckpt_path=Path("/tmp/student.pt"),
    )
    assert enabled is False, f"Expected False when disabled, got {enabled}"
    assert "disabled" in reason.lower() or "not enabled" in reason.lower(), \
        f"Expected reason to mention 'disabled', got: {reason}"
    print("PASS: test_disabled_returns_false")


def test_same_kmer_disables():
    """When teacher k-mer equals student k-mer, distillation should be disabled."""
    enabled, reason = apply_distill_guardrails(
        distill_enabled=True,
        teacher_kmer=6,
        student_kmer=6,
        teacher_ckpt_path=Path("/tmp/teacher.pt"),
        student_ckpt_path=Path("/tmp/student.pt"),
    )
    assert enabled is False, f"Expected False for same k-mer, got {enabled}"
    assert "6" in reason, f"Expected reason to mention k-mer value, got: {reason}"
    print("PASS: test_same_kmer_disables")


def test_same_checkpoint_path_disables():
    """When teacher and student checkpoint paths are the same, distillation should be disabled."""
    same_path = Path("/tmp/same_checkpoint.pt")
    enabled, reason = apply_distill_guardrails(
        distill_enabled=True,
        teacher_kmer=6,
        student_kmer=4,
        teacher_ckpt_path=same_path,
        student_ckpt_path=same_path,
    )
    assert enabled is False, f"Expected False for same checkpoint path, got {enabled}"
    assert "path" in reason.lower() or "checkpoint" in reason.lower(), \
        f"Expected reason to mention path/checkpoint, got: {reason}"
    print("PASS: test_same_checkpoint_path_disables")


def test_different_kmer_and_paths_enables():
    """When k-mers differ and paths differ, distillation should be enabled."""
    enabled, reason = apply_distill_guardrails(
        distill_enabled=True,
        teacher_kmer=6,
        student_kmer=4,
        teacher_ckpt_path=Path("/tmp/teacher.pt"),
        student_ckpt_path=Path("/tmp/student.pt"),
    )
    assert enabled is True, f"Expected True for valid distillation setup, got {enabled}"
    assert reason == "", f"Expected empty reason, got: {reason}"
    print("PASS: test_different_kmer_and_paths_enables")


if __name__ == "__main__":
    print("Running distillation guardrail tests...\n")

    test_disabled_returns_false()
    test_same_kmer_disables()
    test_same_checkpoint_path_disables()
    test_different_kmer_and_paths_enables()

    print("\nOK")
