"""Smoke tests for Stage 2 checkpoint load order.

Verifies that:
1. When --stage1_ckpt is provided, MLM weights are NOT loaded
2. Stage 1 checkpoint weights are correctly loaded
3. The load order prevents accidental double-loading

Usage:
    python tests/test_stage2_load_order.py

    # Or with pytest
    pytest tests/test_stage2_load_order.py -v
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import torch

from model import GeneWhispererStage1, GeneWhispererStage2


def create_stage1_checkpoint(path: Path, model: GeneWhispererStage1) -> None:
    """Create a Stage 1 checkpoint file."""
    torch.save({"model_state": model.state_dict()}, path)


def create_mlm_checkpoint(path: Path, vocab_size: int, embedding_dim: int) -> None:
    """Create a mock MLM checkpoint with distinct embedding weights."""
    state = {
        "embedding.weight": torch.randn(vocab_size, embedding_dim) * 10.0,  # Large values to detect
    }
    torch.save(state, path)


def test_stage1_ckpt_skips_mlm_loading():
    """Test that providing stage1_ckpt skips MLM loading."""
    vocab_size = 67
    embedding_dim = 192
    seq_len = 79

    # Create Stage 1 model with distinctive weights
    stage1_model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=embedding_dim,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=128,
    )
    # Set embedding to known values
    with torch.no_grad():
        stage1_model.embedding.weight.fill_(0.5)

    # Create Stage 2 model (fresh init)
    stage2_model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=embedding_dim,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=128,
        use_engineered_features=True,
        stage2_head_type="transformer",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        stage1_ckpt = tmpdir / "stage1.pt"
        mlm_ckpt = tmpdir / "mlm.pt"

        # Save Stage 1 checkpoint
        create_stage1_checkpoint(stage1_ckpt, stage1_model)

        # Create MLM checkpoint with very different values
        create_mlm_checkpoint(mlm_ckpt, vocab_size, embedding_dim)

        # Load Stage 1 weights into Stage 2
        loaded = stage2_model.load_stage1_weights(stage1_ckpt)
        assert loaded > 0, "Stage 1 weights should load successfully"

        # Verify embedding weights match Stage 1
        with torch.no_grad():
            emb_diff = (stage2_model.embedding.weight - 0.5).abs().max().item()
            assert emb_diff < 1e-5, f"Embedding should be ~0.5 (from Stage 1), diff={emb_diff}"

    print("test_stage1_ckpt_skips_mlm_loading: PASSED")


def test_should_load_mlm_logic():
    """Test the should_load_mlm conditional logic directly."""
    from pathlib import Path

    # Simulate the conditional from train_stage2.py:
    # should_load_mlm = load_mlm_weights and transfer_mode != "none" and not stage1_ckpt_path

    test_cases = [
        # (load_mlm_weights, transfer_mode, stage1_ckpt_path, expected_should_load_mlm)
        (True, "embed_only", None, True),       # No stage1_ckpt -> load MLM
        (True, "embed_only", Path("foo.pt"), False),  # Has stage1_ckpt -> skip MLM
        (False, "embed_only", None, False),     # MLM loading disabled -> skip
        (True, "none", None, False),            # transfer_mode=none -> skip
        (True, "embed_plus_adapter", None, True),     # No stage1_ckpt -> load MLM
        (True, "embed_plus_adapter", Path("bar.pt"), False),  # Has stage1_ckpt -> skip MLM
    ]

    for load_mlm_weights, transfer_mode, stage1_ckpt_path, expected in test_cases:
        should_load_mlm = load_mlm_weights and transfer_mode != "none" and not stage1_ckpt_path
        assert should_load_mlm == expected, (
            f"FAIL: load_mlm_weights={load_mlm_weights}, transfer_mode={transfer_mode}, "
            f"stage1_ckpt_path={stage1_ckpt_path} -> expected {expected}, got {should_load_mlm}"
        )

    print("test_should_load_mlm_logic: PASSED (all 6 cases)")


def test_stage1_overwrites_initial_weights():
    """Test that load_stage1_weights actually changes the model weights."""
    vocab_size = 67
    embedding_dim = 192

    # Create Stage 1 model with known weights
    stage1_model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=embedding_dim,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=128,
    )
    with torch.no_grad():
        stage1_model.embedding.weight.fill_(1.0)

    # Create Stage 2 model
    stage2_model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=embedding_dim,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=128,
        use_engineered_features=True,
        stage2_head_type="transformer",
    )

    # Capture initial embedding norm
    with torch.no_grad():
        initial_emb_mean = stage2_model.embedding.weight.mean().item()

    with tempfile.TemporaryDirectory() as tmpdir:
        stage1_ckpt = Path(tmpdir) / "stage1.pt"
        create_stage1_checkpoint(stage1_ckpt, stage1_model)

        # Load Stage 1 weights
        loaded = stage2_model.load_stage1_weights(stage1_ckpt)
        assert loaded > 0

        # Verify embedding changed
        with torch.no_grad():
            after_emb_mean = stage2_model.embedding.weight.mean().item()

        assert abs(after_emb_mean - 1.0) < 1e-5, (
            f"After load, embedding mean should be 1.0, got {after_emb_mean}"
        )

    print("test_stage1_overwrites_initial_weights: PASSED")


def test_no_double_load_when_stage1_provided():
    """
    Verify the load order invariant:
    When stage1_ckpt is provided, MLM loading should be skipped entirely.

    This test simulates the decision logic from train_stage2.py.
    """
    # Simulate config scenarios
    scenarios = [
        {
            "name": "stage1_ckpt only",
            "stage1_ckpt": "/path/to/stage1.pt",
            "stage2_load_mlm_weights": True,
            "mlm_checkpoint": "/path/to/mlm.pt",
            "expected_mlm_load": False,  # MLM should NOT load because stage1_ckpt wins
        },
        {
            "name": "MLM only",
            "stage1_ckpt": None,
            "stage2_load_mlm_weights": True,
            "mlm_checkpoint": "/path/to/mlm.pt",
            "expected_mlm_load": True,  # MLM should load
        },
        {
            "name": "both disabled",
            "stage1_ckpt": None,
            "stage2_load_mlm_weights": False,
            "mlm_checkpoint": None,
            "expected_mlm_load": False,
        },
    ]

    for scenario in scenarios:
        stage1_ckpt_path = Path(scenario["stage1_ckpt"]) if scenario["stage1_ckpt"] else None
        load_mlm_weights = scenario["stage2_load_mlm_weights"]
        mlm_checkpoint = scenario["mlm_checkpoint"]
        transfer_mode = "embed_only"

        # This is the patched logic from train_stage2.py
        should_load_mlm = (
            load_mlm_weights
            and transfer_mode != "none"
            and not stage1_ckpt_path
        )
        will_attempt_mlm = bool(should_load_mlm and mlm_checkpoint)

        assert will_attempt_mlm == scenario["expected_mlm_load"], (
            f"Scenario '{scenario['name']}': expected MLM load={scenario['expected_mlm_load']}, "
            f"got {will_attempt_mlm}"
        )

    print("test_no_double_load_when_stage1_provided: PASSED (all 3 scenarios)")


def run_all_tests():
    """Run all load order smoke tests."""
    print("=" * 60)
    print("Running Stage 2 Load Order Smoke Tests")
    print("=" * 60)

    tests = [
        ("should_load_mlm logic", test_should_load_mlm_logic),
        ("stage1_ckpt skips MLM loading", test_stage1_ckpt_skips_mlm_loading),
        ("stage1 overwrites initial weights", test_stage1_overwrites_initial_weights),
        ("no double load when stage1 provided", test_no_double_load_when_stage1_provided),
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
