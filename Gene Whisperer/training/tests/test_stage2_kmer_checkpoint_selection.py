"""Smoke tests for Stage 2 k-mer-based checkpoint auto-selection.

Verifies that:
1. Stage 2 with k=4 auto-selects stage1_ckpt_by_k[4]
2. Stage 2 with k=6 auto-selects stage1_ckpt_by_k[6]
3. Explicit --stage1_ckpt overrides auto-selection
4. Missing checkpoint in stage1_ckpt_by_k falls back gracefully

Usage:
    python tests/test_stage2_kmer_checkpoint_selection.py

    # Or with pytest
    pytest tests/test_stage2_kmer_checkpoint_selection.py -v
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def _resolve_optional_path(value: Any, *, base_dir: Path) -> Optional[Path]:
    """Copy of the helper from train_stage2.py for testing."""
    if value is None:
        return None
    if isinstance(value, Path):
        candidate = value
    else:
        text = str(value).strip()
        if not text:
            return None
        candidate = Path(text)
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    return candidate


def resolve_stage1_ckpt_for_kmer(
    cfg: Dict[str, Any],
    stage2_kmer: int,
    explicit_stage1_ckpt: Optional[str],
    base_dir: Path,
) -> tuple[Optional[Path], bool]:
    """
    Resolve the Stage 1 checkpoint path based on k-mer.
    Returns (path, was_auto_selected).

    This mirrors the logic from train_stage2.py.
    """
    stage1_ckpt_path = _resolve_optional_path(explicit_stage1_ckpt, base_dir=base_dir)
    stage1_ckpt_auto_selected = False

    if stage1_ckpt_path is None:
        stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k")
        if isinstance(stage1_ckpt_by_k, dict):
            # Keys may be int or str depending on YAML parsing
            ckpt_for_k = stage1_ckpt_by_k.get(stage2_kmer) or stage1_ckpt_by_k.get(str(stage2_kmer))
            if ckpt_for_k:
                stage1_ckpt_path = _resolve_optional_path(ckpt_for_k, base_dir=base_dir)
                if stage1_ckpt_path and stage1_ckpt_path.exists():
                    stage1_ckpt_auto_selected = True
                elif stage1_ckpt_path:
                    # File doesn't exist, reset to None
                    stage1_ckpt_path = None

    return stage1_ckpt_path, stage1_ckpt_auto_selected


def test_k4_selects_k4_checkpoint():
    """Test that Stage 2 with k=4 auto-selects stage1_k4.pt."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ckpt_k4 = tmpdir / "stage1_k4.pt"
        ckpt_k6 = tmpdir / "stage1_k6.pt"

        # Create dummy checkpoint files
        ckpt_k4.write_bytes(b"dummy_k4")
        ckpt_k6.write_bytes(b"dummy_k6")

        cfg = {
            "stage1_ckpt_by_k": {
                4: str(ckpt_k4),
                6: str(ckpt_k6),
            }
        }

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=4, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is not None, "k=4 should resolve to a path"
        assert path == ckpt_k4, f"k=4 should select stage1_k4.pt, got {path}"
        assert auto_selected is True, "Should be marked as auto-selected"

    print("test_k4_selects_k4_checkpoint: PASSED")


def test_k6_selects_k6_checkpoint():
    """Test that Stage 2 with k=6 auto-selects stage1_k6.pt."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ckpt_k4 = tmpdir / "stage1_k4.pt"
        ckpt_k6 = tmpdir / "stage1_k6.pt"

        # Create dummy checkpoint files
        ckpt_k4.write_bytes(b"dummy_k4")
        ckpt_k6.write_bytes(b"dummy_k6")

        cfg = {
            "stage1_ckpt_by_k": {
                4: str(ckpt_k4),
                6: str(ckpt_k6),
            }
        }

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=6, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is not None, "k=6 should resolve to a path"
        assert path == ckpt_k6, f"k=6 should select stage1_k6.pt, got {path}"
        assert auto_selected is True, "Should be marked as auto-selected"

    print("test_k6_selects_k6_checkpoint: PASSED")


def test_explicit_ckpt_overrides_auto_selection():
    """Test that --stage1_ckpt CLI arg overrides auto-selection."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ckpt_k4 = tmpdir / "stage1_k4.pt"
        ckpt_k6 = tmpdir / "stage1_k6.pt"
        ckpt_custom = tmpdir / "custom_checkpoint.pt"

        # Create dummy checkpoint files
        ckpt_k4.write_bytes(b"dummy_k4")
        ckpt_k6.write_bytes(b"dummy_k6")
        ckpt_custom.write_bytes(b"dummy_custom")

        cfg = {
            "stage1_ckpt_by_k": {
                4: str(ckpt_k4),
                6: str(ckpt_k6),
            }
        }

        # Even though k=6, the explicit path should win
        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=6, explicit_stage1_ckpt=str(ckpt_custom), base_dir=tmpdir
        )

        assert path is not None, "Explicit path should resolve"
        assert path == ckpt_custom, f"Explicit path should override, got {path}"
        assert auto_selected is False, "Should NOT be marked as auto-selected"

    print("test_explicit_ckpt_overrides_auto_selection: PASSED")


def test_missing_kmer_in_mapping_returns_none():
    """Test that a missing k-mer in stage1_ckpt_by_k returns None gracefully."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ckpt_k4 = tmpdir / "stage1_k4.pt"

        # Only create k4 checkpoint
        ckpt_k4.write_bytes(b"dummy_k4")

        cfg = {
            "stage1_ckpt_by_k": {
                4: str(ckpt_k4),
                # k=6 is NOT in the mapping
            }
        }

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=6, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is None, f"k=6 (not in mapping) should return None, got {path}"
        assert auto_selected is False, "Should NOT be marked as auto-selected"

    print("test_missing_kmer_in_mapping_returns_none: PASSED")


def test_nonexistent_checkpoint_file_falls_back():
    """Test that a non-existent checkpoint file in stage1_ckpt_by_k falls back to None."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        cfg = {
            "stage1_ckpt_by_k": {
                4: str(tmpdir / "nonexistent_k4.pt"),  # File doesn't exist
                6: str(tmpdir / "nonexistent_k6.pt"),  # File doesn't exist
            }
        }

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=4, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is None, f"Non-existent file should return None, got {path}"
        assert auto_selected is False, "Should NOT be marked as auto-selected"

    print("test_nonexistent_checkpoint_file_falls_back: PASSED")


def test_string_keys_in_mapping():
    """Test that string keys (from YAML parsing) also work."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ckpt_k4 = tmpdir / "stage1_k4.pt"

        ckpt_k4.write_bytes(b"dummy_k4")

        # YAML sometimes parses integer keys as strings
        cfg = {
            "stage1_ckpt_by_k": {
                "4": str(ckpt_k4),  # String key
            }
        }

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=4, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is not None, "String key '4' should also match"
        assert path == ckpt_k4, f"String key should work, got {path}"
        assert auto_selected is True

    print("test_string_keys_in_mapping: PASSED")


def test_no_mapping_configured():
    """Test that no stage1_ckpt_by_k config returns None gracefully."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        cfg = {}  # No stage1_ckpt_by_k at all

        path, auto_selected = resolve_stage1_ckpt_for_kmer(
            cfg, stage2_kmer=6, explicit_stage1_ckpt=None, base_dir=tmpdir
        )

        assert path is None, "No mapping should return None"
        assert auto_selected is False

    print("test_no_mapping_configured: PASSED")


def run_all_tests():
    """Run all k-mer checkpoint selection smoke tests."""
    print("=" * 60)
    print("Running Stage 2 K-mer Checkpoint Selection Smoke Tests")
    print("=" * 60)

    tests = [
        ("k=4 selects k4 checkpoint", test_k4_selects_k4_checkpoint),
        ("k=6 selects k6 checkpoint", test_k6_selects_k6_checkpoint),
        ("explicit --stage1_ckpt overrides", test_explicit_ckpt_overrides_auto_selection),
        ("missing k-mer returns None", test_missing_kmer_in_mapping_returns_none),
        ("nonexistent file falls back", test_nonexistent_checkpoint_file_falls_back),
        ("string keys in mapping", test_string_keys_in_mapping),
        ("no mapping configured", test_no_mapping_configured),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            print(f"\n[TEST] {name}")
            test_func()
            print(f"[PASS] {name}")
            passed += 1
        except Exception as e:
            print(f"[FAIL] {name}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
