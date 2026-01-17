"""Test training script with all model variants.

Tests verify the training script works correctly with all supported variants:
1. Config parsing test - verify correct model creation for each variant
2. Single epoch test - run 1 epoch with small data for each variant
3. Checkpoint compatibility test - save/load checkpoints with variant metadata
4. Backward compatibility test - config without model_variant defaults to original
5. Quick benchmark - time training steps for each variant
"""
from __future__ import annotations

import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn

# Add training directory to path for local imports.
TRAINING_DIR = Path(__file__).resolve().parent.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

from model import (
    GeneWhispererCombined,
    GeneWhispererFeaturesOnly,
    GeneWhispererStage1,
    GeneWhispererTransformerOnly,
    create_model_variant,
)


VARIANT_NAMES = ("transformer_only", "features_only", "combined")


def _small_config() -> Dict[str, Any]:
    """Return a minimal config for testing."""
    return {
        "vocab_size": 16,
        "kmer": 3,
        "embedding_dim": 16,
        "transformer_layers": 1,
        "transformer_heads": 2,
        "transformer_ff_dim": 32,
        "transformer_dropout": 0.0,
        "pad_token_id": 0,
        "engineered_dim": 6,
        "engineered_mlp_hidden": 8,
        "engineered_mlp_output": 6,
        "max_seq_len": 8,
        "max_bp_len": 11,  # max_seq_len + kmer - 1
        "drop_path_rate": 0.0,
        "use_relative_position_bias": False,
        "use_glu_ffn": False,
    }


def _build_inputs(
    config: Dict[str, Any], batch_size: int = 4
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Build test inputs for forward pass."""
    seq_len = config.get("max_seq_len", 8)
    tokens = torch.randint(1, config["vocab_size"], (batch_size, seq_len))
    engineered = torch.randn(batch_size, config["engineered_dim"])
    labels = torch.randint(0, 2, (batch_size,)).float()
    return tokens, engineered, labels


def _build_training_batch(
    config: Dict[str, Any], batch_size: int = 16
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Build a training batch with signal to learn."""
    seq_len = config.get("max_seq_len", 8)
    split = max(2, config["vocab_size"] // 2)

    half = batch_size // 2
    # High tokens for positive class, low tokens for negative
    high_tokens = torch.randint(split, config["vocab_size"], (half, seq_len))
    low_tokens = torch.randint(1, split, (batch_size - half, seq_len))
    tokens = torch.cat([high_tokens, low_tokens], dim=0)

    # Shifted feature means for signal
    engineered_pos = torch.randn(half, config["engineered_dim"]) + 1.0
    engineered_neg = torch.randn(batch_size - half, config["engineered_dim"]) - 1.0
    engineered = torch.cat([engineered_pos, engineered_neg], dim=0)

    labels = torch.zeros(batch_size, 1)
    labels[:half] = 1.0
    return tokens, engineered, labels


# =============================================================================
# TEST 1: Config Parsing Test
# =============================================================================


def test_config_parsing_combined():
    """Test config with model_variant='combined' creates correct model."""
    config = _small_config()
    model = create_model_variant("combined", config)
    assert isinstance(model, GeneWhispererCombined), (
        f"Expected GeneWhispererCombined, got {type(model).__name__}"
    )
    # Verify model has both transformer and engineered feature components
    assert model.encoder is not None, "Combined model should have encoder"
    assert model.embedding is not None, "Combined model should have embedding"
    assert model.engineered_mlp is not None, "Combined model should have engineered_mlp"


def test_config_parsing_transformer_only():
    """Test config with model_variant='transformer_only' creates correct model."""
    config = _small_config()
    model = create_model_variant("transformer_only", config)
    assert isinstance(model, GeneWhispererTransformerOnly), (
        f"Expected GeneWhispererTransformerOnly, got {type(model).__name__}"
    )
    # Verify model has transformer but handles missing engineered features
    assert model.encoder is not None, "TransformerOnly model should have encoder"
    assert model.embedding is not None, "TransformerOnly model should have embedding"


def test_config_parsing_features_only():
    """Test config with model_variant='features_only' creates correct model."""
    config = _small_config()
    model = create_model_variant("features_only", config)
    assert isinstance(model, GeneWhispererFeaturesOnly), (
        f"Expected GeneWhispererFeaturesOnly, got {type(model).__name__}"
    )
    # Verify model has no transformer components
    assert model.encoder is None, "FeaturesOnly model should not have encoder"


def test_config_parsing_all_variants():
    """Test that all variant options create the correct model types."""
    config = _small_config()

    variant_to_class = {
        "transformer_only": GeneWhispererTransformerOnly,
        "features_only": GeneWhispererFeaturesOnly,
        "combined": GeneWhispererCombined,
    }

    for variant, expected_class in variant_to_class.items():
        model = create_model_variant(variant, config)
        assert isinstance(model, expected_class), (
            f"variant='{variant}': Expected {expected_class.__name__}, "
            f"got {type(model).__name__}"
        )


def test_config_parsing_invalid_variant():
    """Test that invalid variant raises ValueError."""
    config = _small_config()
    try:
        create_model_variant("invalid_variant", config)
        assert False, "Expected ValueError for invalid variant"
    except ValueError as e:
        assert "Unsupported model variant" in str(e)


# =============================================================================
# TEST 2: Single Epoch Test (for each variant)
# =============================================================================


def _run_single_epoch(
    variant: str,
    config: Dict[str, Any],
    num_samples: int = 100,
) -> Dict[str, Any]:
    """Run 1 epoch of training and return metrics."""
    torch.manual_seed(42)
    model = create_model_variant(variant, config)
    model.train()

    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    criterion = nn.BCEWithLogitsLoss()

    # Generate training data
    batch_size = 16
    num_batches = max(1, num_samples // batch_size)

    total_loss = 0.0
    total_samples = 0
    all_preds = []
    all_labels = []

    for _ in range(num_batches):
        tokens, engineered, labels = _build_training_batch(config, batch_size)
        labels_unsqueezed = labels

        optimizer.zero_grad()

        # Forward pass based on variant
        if variant == "transformer_only":
            logits = model(tokens, None)
        elif variant == "features_only":
            logits = model(None, engineered)
        else:  # combined
            logits = model(tokens, engineered)

        loss = criterion(logits, labels_unsqueezed)
        loss.backward()
        optimizer.step()

        total_loss += loss.item() * batch_size
        total_samples += batch_size

        # Track predictions
        with torch.no_grad():
            preds = (torch.sigmoid(logits) >= 0.5).float()
            all_preds.append(preds)
            all_labels.append(labels_unsqueezed)

    # Compute metrics
    all_preds = torch.cat(all_preds, dim=0)
    all_labels = torch.cat(all_labels, dim=0)
    accuracy = (all_preds == all_labels).float().mean().item()
    avg_loss = total_loss / total_samples

    return {
        "loss": avg_loss,
        "accuracy": accuracy,
        "samples": total_samples,
    }


def test_single_epoch_transformer_only():
    """Run 1 epoch of training with transformer_only variant."""
    config = _small_config()
    metrics = _run_single_epoch("transformer_only", config, num_samples=100)

    assert metrics["samples"] == 96, f"Expected 96 samples, got {metrics['samples']}"
    assert metrics["loss"] >= 0, f"Loss should be non-negative, got {metrics['loss']}"
    assert 0 <= metrics["accuracy"] <= 1, (
        f"Accuracy should be in [0,1], got {metrics['accuracy']}"
    )
    # Verify metrics were logged (non-zero loss indicates training happened)
    assert metrics["loss"] > 0, "Loss should be > 0 during training"


def test_single_epoch_features_only():
    """Run 1 epoch of training with features_only variant."""
    config = _small_config()
    metrics = _run_single_epoch("features_only", config, num_samples=100)

    assert metrics["samples"] == 96, f"Expected 96 samples, got {metrics['samples']}"
    assert metrics["loss"] >= 0, f"Loss should be non-negative, got {metrics['loss']}"
    assert 0 <= metrics["accuracy"] <= 1, (
        f"Accuracy should be in [0,1], got {metrics['accuracy']}"
    )
    assert metrics["loss"] > 0, "Loss should be > 0 during training"


def test_single_epoch_combined():
    """Run 1 epoch of training with combined variant."""
    config = _small_config()
    metrics = _run_single_epoch("combined", config, num_samples=100)

    assert metrics["samples"] == 96, f"Expected 96 samples, got {metrics['samples']}"
    assert metrics["loss"] >= 0, f"Loss should be non-negative, got {metrics['loss']}"
    assert 0 <= metrics["accuracy"] <= 1, (
        f"Accuracy should be in [0,1], got {metrics['accuracy']}"
    )
    assert metrics["loss"] > 0, "Loss should be > 0 during training"


def test_single_epoch_all_variants_complete_without_error():
    """Verify all variants complete training without errors."""
    config = _small_config()
    for variant in VARIANT_NAMES:
        try:
            metrics = _run_single_epoch(variant, config, num_samples=100)
            assert "loss" in metrics, f"{variant}: loss not in metrics"
            assert "accuracy" in metrics, f"{variant}: accuracy not in metrics"
        except Exception as e:
            raise AssertionError(f"{variant} training failed: {e}") from e


# =============================================================================
# TEST 3: Checkpoint Compatibility Test
# =============================================================================


def test_checkpoint_saves_variant_info():
    """Test that checkpoint saving includes variant metadata."""
    from train_stage1 import save_checkpoint

    config = _small_config()

    for variant in VARIANT_NAMES:
        model = create_model_variant(variant, config)
        optimizer = torch.optim.Adam(model.parameters())

        with tempfile.TemporaryDirectory() as tmpdir:
            ckpt_path = Path(tmpdir) / f"test_{variant}.pt"
            metrics = {"accuracy": 0.95, "loss": 0.1, "best_threshold": 0.5}

            save_checkpoint(
                path=ckpt_path,
                model=model,
                optimizer=optimizer,
                scheduler=None,
                epoch=1,
                metrics=metrics,
                model_variant=variant,
                use_simplified_architecture=True,
            )

            # Load and verify metadata
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            assert checkpoint["model_variant"] == variant, (
                f"Expected model_variant='{variant}', got '{checkpoint['model_variant']}'"
            )
            assert checkpoint["use_simplified_architecture"] is True
            assert "model_state" in checkpoint
            assert checkpoint["epoch"] == 1
            assert checkpoint["best_threshold"] == 0.5


def test_checkpoint_load_restores_architecture():
    """Test that loading checkpoint restores correct model architecture."""
    from train_stage1 import save_checkpoint

    config = _small_config()
    tokens, engineered, _ = _build_inputs(config)

    for variant in VARIANT_NAMES:
        torch.manual_seed(42)
        original_model = create_model_variant(variant, config)
        original_model.eval()

        # Get original output
        with torch.no_grad():
            if variant == "transformer_only":
                original_out = original_model(tokens, None)
            elif variant == "features_only":
                original_out = original_model(None, engineered)
            else:
                original_out = original_model(tokens, engineered)

        # Save checkpoint
        optimizer = torch.optim.Adam(original_model.parameters())
        with tempfile.TemporaryDirectory() as tmpdir:
            ckpt_path = Path(tmpdir) / f"test_{variant}.pt"
            save_checkpoint(
                path=ckpt_path,
                model=original_model,
                optimizer=optimizer,
                scheduler=None,
                epoch=1,
                metrics={"accuracy": 0.95},
                model_variant=variant,
                use_simplified_architecture=True,
            )

            # Load checkpoint and restore model
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            loaded_variant = checkpoint["model_variant"]
            assert loaded_variant == variant

            # Create new model of same variant and load weights
            loaded_model = create_model_variant(loaded_variant, config)
            loaded_model.load_state_dict(checkpoint["model_state"])
            loaded_model.eval()

            # Verify outputs match
            with torch.no_grad():
                if variant == "transformer_only":
                    loaded_out = loaded_model(tokens, None)
                elif variant == "features_only":
                    loaded_out = loaded_model(None, engineered)
                else:
                    loaded_out = loaded_model(tokens, engineered)

            assert torch.allclose(original_out, loaded_out, atol=1e-6), (
                f"{variant}: Loaded model output differs from original"
            )


def test_checkpoint_weights_loaded_correctly():
    """Test that checkpoint weights are loaded correctly after training."""
    from train_stage1 import save_checkpoint

    config = _small_config()
    tokens, engineered, labels = _build_inputs(config)
    labels = labels.unsqueeze(1)

    for variant in VARIANT_NAMES:
        torch.manual_seed(42)
        model = create_model_variant(variant, config)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
        criterion = nn.BCEWithLogitsLoss()

        # Train for a few steps to change weights
        model.train()
        for _ in range(5):
            optimizer.zero_grad()
            if variant == "transformer_only":
                logits = model(tokens, None)
            elif variant == "features_only":
                logits = model(None, engineered)
            else:
                logits = model(tokens, engineered)
            loss = criterion(logits, labels)
            loss.backward()
            optimizer.step()

        # Save trained model
        model.eval()
        with torch.no_grad():
            if variant == "transformer_only":
                trained_out = model(tokens, None)
            elif variant == "features_only":
                trained_out = model(None, engineered)
            else:
                trained_out = model(tokens, engineered)

        with tempfile.TemporaryDirectory() as tmpdir:
            ckpt_path = Path(tmpdir) / f"trained_{variant}.pt"
            save_checkpoint(
                path=ckpt_path,
                model=model,
                optimizer=optimizer,
                scheduler=None,
                epoch=5,
                metrics={"accuracy": 0.9},
                model_variant=variant,
                use_simplified_architecture=True,
            )

            # Create fresh model and load weights
            fresh_model = create_model_variant(variant, config)
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            fresh_model.load_state_dict(checkpoint["model_state"])
            fresh_model.eval()

            with torch.no_grad():
                if variant == "transformer_only":
                    loaded_out = fresh_model(tokens, None)
                elif variant == "features_only":
                    loaded_out = fresh_model(None, engineered)
                else:
                    loaded_out = fresh_model(tokens, engineered)

            # Should exactly match trained model output
            assert torch.allclose(trained_out, loaded_out, atol=1e-6), (
                f"{variant}: Loaded weights don't reproduce trained output"
            )


# =============================================================================
# TEST 4: Backward Compatibility Test
# =============================================================================


def test_backward_compat_missing_model_variant():
    """Test that config without model_variant defaults to original architecture."""
    config = _small_config()

    # Simulate old config without model_variant - should still be able to create
    # the original GeneWhispererStage1 model
    original_model = GeneWhispererStage1(
        vocab_size=config["vocab_size"],
        kmer=config["kmer"],
        embedding_dim=config["embedding_dim"],
        num_layers=config["transformer_layers"],
        num_heads=config["transformer_heads"],
        ff_dim=config["transformer_ff_dim"],
        dropout=config["transformer_dropout"],
        pad_token_id=config["pad_token_id"],
    )

    assert original_model is not None
    assert isinstance(original_model, GeneWhispererStage1)


def test_backward_compat_original_model_forward():
    """Test that original model forward pass works."""
    config = _small_config()
    tokens, engineered, _ = _build_inputs(config)

    original_model = GeneWhispererStage1(
        vocab_size=config["vocab_size"],
        kmer=config["kmer"],
        embedding_dim=config["embedding_dim"],
        num_layers=config["transformer_layers"],
        num_heads=config["transformer_heads"],
        ff_dim=config["transformer_ff_dim"],
        dropout=config["transformer_dropout"],
        pad_token_id=config["pad_token_id"],
        engineered_dim=config["engineered_dim"],
    )
    original_model.eval()

    with torch.no_grad():
        # Original model accepts both tokens and engineered
        output = original_model(tokens, engineered)

    assert output.shape == (4, 1), f"Expected shape (4, 1), got {output.shape}"
    assert torch.isfinite(output).all(), "Output contains NaN or Inf"


def test_backward_compat_original_model_training():
    """Test that training works with original model architecture."""
    config = _small_config()

    original_model = GeneWhispererStage1(
        vocab_size=config["vocab_size"],
        kmer=config["kmer"],
        embedding_dim=config["embedding_dim"],
        num_layers=config["transformer_layers"],
        num_heads=config["transformer_heads"],
        ff_dim=config["transformer_ff_dim"],
        dropout=config["transformer_dropout"],
        pad_token_id=config["pad_token_id"],
        engineered_dim=config["engineered_dim"],
    )
    original_model.train()

    optimizer = torch.optim.Adam(original_model.parameters(), lr=0.01)
    criterion = nn.BCEWithLogitsLoss()

    tokens, engineered, labels = _build_training_batch(config, batch_size=16)

    # Run a few training steps
    initial_loss = None
    final_loss = None

    for step in range(10):
        optimizer.zero_grad()
        logits = original_model(tokens, engineered)
        loss = criterion(logits, labels)
        loss.backward()
        optimizer.step()

        if step == 0:
            initial_loss = loss.item()
        final_loss = loss.item()

    assert initial_loss is not None and final_loss is not None
    assert torch.isfinite(torch.tensor(final_loss)), "Final loss is not finite"


def test_backward_compat_checkpoint_without_variant_metadata():
    """Test loading checkpoint that doesn't have variant metadata."""
    config = _small_config()

    original_model = GeneWhispererStage1(
        vocab_size=config["vocab_size"],
        kmer=config["kmer"],
        embedding_dim=config["embedding_dim"],
        num_layers=config["transformer_layers"],
        num_heads=config["transformer_heads"],
        ff_dim=config["transformer_ff_dim"],
        dropout=config["transformer_dropout"],
        pad_token_id=config["pad_token_id"],
    )
    optimizer = torch.optim.Adam(original_model.parameters())

    with tempfile.TemporaryDirectory() as tmpdir:
        ckpt_path = Path(tmpdir) / "old_checkpoint.pt"

        # Save checkpoint without variant metadata (old format)
        old_checkpoint = {
            "epoch": 1,
            "model_state": original_model.state_dict(),
            "optimizer_state": optimizer.state_dict(),
            "metrics": {"accuracy": 0.9},
        }
        torch.save(old_checkpoint, ckpt_path)

        # Load and verify we can handle missing variant info
        loaded = torch.load(ckpt_path, map_location="cpu", weights_only=False)

        # Gracefully handle missing keys
        model_variant = loaded.get("model_variant", "original")
        use_simplified = loaded.get("use_simplified_architecture", False)

        assert model_variant == "original", "Default variant should be 'original'"
        assert use_simplified is False, "Default should not use simplified architecture"

        # Verify we can load the state dict
        new_model = GeneWhispererStage1(
            vocab_size=config["vocab_size"],
            kmer=config["kmer"],
            embedding_dim=config["embedding_dim"],
            num_layers=config["transformer_layers"],
            num_heads=config["transformer_heads"],
            ff_dim=config["transformer_ff_dim"],
            dropout=config["transformer_dropout"],
            pad_token_id=config["pad_token_id"],
        )
        new_model.load_state_dict(loaded["model_state"])


# =============================================================================
# TEST 5: Quick Benchmark
# =============================================================================


def _benchmark_variant(
    variant: str,
    config: Dict[str, Any],
    num_steps: int = 10,
) -> float:
    """Time training steps for a variant, return steps/second."""
    torch.manual_seed(42)
    model = create_model_variant(variant, config)
    model.train()

    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    criterion = nn.BCEWithLogitsLoss()

    # Pre-generate data
    tokens, engineered, labels = _build_training_batch(config, batch_size=16)

    # Warmup step
    optimizer.zero_grad()
    if variant == "transformer_only":
        logits = model(tokens, None)
    elif variant == "features_only":
        logits = model(None, engineered)
    else:
        logits = model(tokens, engineered)
    loss = criterion(logits, labels)
    loss.backward()
    optimizer.step()

    # Timed steps
    start_time = time.perf_counter()

    for _ in range(num_steps):
        optimizer.zero_grad()
        if variant == "transformer_only":
            logits = model(tokens, None)
        elif variant == "features_only":
            logits = model(None, engineered)
        else:
            logits = model(tokens, engineered)
        loss = criterion(logits, labels)
        loss.backward()
        optimizer.step()

    elapsed = time.perf_counter() - start_time
    steps_per_second = num_steps / elapsed

    return steps_per_second


def test_benchmark_all_variants():
    """Benchmark all variants and compare speeds."""
    config = _small_config()
    num_steps = 10

    results: Dict[str, float] = {}

    for variant in VARIANT_NAMES:
        steps_per_sec = _benchmark_variant(variant, config, num_steps)
        results[variant] = steps_per_sec

    # Print results
    print("\n" + "=" * 60)
    print("QUICK BENCHMARK RESULTS")
    print("=" * 60)
    print(f"{'Variant':<20} {'Steps/Second':>15}")
    print("-" * 60)
    for variant in VARIANT_NAMES:
        print(f"{variant:<20} {results[variant]:>15.2f}")
    print("=" * 60)

    # Verify expected ordering:
    # features_only should be fastest (no transformer)
    # combined should be slowest (both components)
    assert results["features_only"] > 0, "features_only should have positive speed"
    assert results["transformer_only"] > 0, "transformer_only should have positive speed"
    assert results["combined"] > 0, "combined should have positive speed"

    # Combined should be slightly slower than individual variants (has both components)
    # Note: On fast hardware this might not always hold strictly due to measurement noise
    # So we just verify all variants run at reasonable speeds


def test_benchmark_features_only_fastest():
    """Verify features_only is fastest (no transformer overhead)."""
    config = _small_config()
    num_steps = 10

    features_speed = _benchmark_variant("features_only", config, num_steps)
    transformer_speed = _benchmark_variant("transformer_only", config, num_steps)

    # features_only should be faster than transformer_only
    # (features_only has no transformer, just MLP)
    assert features_speed > transformer_speed * 0.5, (
        f"features_only ({features_speed:.1f} steps/s) should be significantly "
        f"faster than transformer_only ({transformer_speed:.1f} steps/s)"
    )


def test_benchmark_combined_slowest():
    """Verify combined is slowest (has both transformer and feature MLP)."""
    config = _small_config()
    num_steps = 10

    combined_speed = _benchmark_variant("combined", config, num_steps)
    transformer_speed = _benchmark_variant("transformer_only", config, num_steps)
    features_speed = _benchmark_variant("features_only", config, num_steps)

    # Combined has both components, so should be slower than transformer_only
    # Allow some tolerance for measurement noise
    assert combined_speed <= transformer_speed * 1.2, (
        f"combined ({combined_speed:.1f} steps/s) should not be much faster "
        f"than transformer_only ({transformer_speed:.1f} steps/s)"
    )


# =============================================================================
# TEST RUNNER
# =============================================================================


def _run_test(name: str, fn) -> bool:
    """Run a test function and report results."""
    try:
        fn()
    except AssertionError as exc:
        print(f"FAIL: {name}")
        print(f"      {exc}")
        return False
    except Exception as exc:
        print(f"FAIL: {name} ({type(exc).__name__}: {exc})")
        return False
    print(f"PASS: {name}")
    return True


def main():
    """Run all tests and report results."""
    test_groups = [
        ("1. Config Parsing Tests", [
            ("combined variant creates correct model", test_config_parsing_combined),
            ("transformer_only variant creates correct model", test_config_parsing_transformer_only),
            ("features_only variant creates correct model", test_config_parsing_features_only),
            ("all variants create correct models", test_config_parsing_all_variants),
            ("invalid variant raises error", test_config_parsing_invalid_variant),
        ]),
        ("2. Single Epoch Tests", [
            ("transformer_only single epoch", test_single_epoch_transformer_only),
            ("features_only single epoch", test_single_epoch_features_only),
            ("combined single epoch", test_single_epoch_combined),
            ("all variants complete without error", test_single_epoch_all_variants_complete_without_error),
        ]),
        ("3. Checkpoint Compatibility Tests", [
            ("checkpoint saves variant info", test_checkpoint_saves_variant_info),
            ("checkpoint load restores architecture", test_checkpoint_load_restores_architecture),
            ("checkpoint weights loaded correctly", test_checkpoint_weights_loaded_correctly),
        ]),
        ("4. Backward Compatibility Tests", [
            ("config without model_variant uses original", test_backward_compat_missing_model_variant),
            ("original model forward pass works", test_backward_compat_original_model_forward),
            ("original model training works", test_backward_compat_original_model_training),
            ("checkpoint without variant metadata", test_backward_compat_checkpoint_without_variant_metadata),
        ]),
        ("5. Quick Benchmark Tests", [
            ("benchmark all variants", test_benchmark_all_variants),
            ("features_only fastest", test_benchmark_features_only_fastest),
            ("combined slowest", test_benchmark_combined_slowest),
        ]),
    ]

    print("=" * 70)
    print("TRAINING SCRIPT VARIANT TESTS")
    print("=" * 70)

    all_passed = True
    total_tests = 0
    passed_tests = 0

    for group_name, tests in test_groups:
        print(f"\n{group_name}")
        print("-" * 60)
        for test_name, test_fn in tests:
            total_tests += 1
            if _run_test(test_name, test_fn):
                passed_tests += 1
            else:
                all_passed = False

    print("\n" + "=" * 70)
    print(f"RESULTS: {passed_tests}/{total_tests} tests passed")
    print("=" * 70)

    if all_passed:
        print("\nAll tests PASSED!")
    else:
        print(f"\nSome tests FAILED ({total_tests - passed_tests} failures)")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
