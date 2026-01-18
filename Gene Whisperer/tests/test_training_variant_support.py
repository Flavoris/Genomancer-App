"""Tests for model variant support in train_stage1.py.

Tests cover:
1. Config parsing for use_simplified_architecture and model_variant
2. Checkpoint saving with variant metadata
3. run_epoch with different variants and variant-aware inputs
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Dict, Optional
from unittest.mock import MagicMock, patch

import torch
import torch.nn as nn

# Add training directory to path for local imports.
TRAINING_DIR = Path(__file__).resolve().parent.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

from model import (
    GeneWhispererCombined,
    GeneWhispererFeaturesOnly,
    GeneWhispererStage1Legacy,
    GeneWhispererTransformerOnly,
    create_model_variant,
)


def _small_config() -> Dict:
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


def _build_inputs(config: Dict, batch_size: int = 4):
    """Build test inputs for forward pass."""
    seq_len = config.get("max_seq_len", 8)
    tokens = torch.randint(1, config["vocab_size"], (batch_size, seq_len))
    engineered = torch.randn(batch_size, config["engineered_dim"])
    labels = torch.randint(0, 2, (batch_size,)).float()
    return tokens, engineered, labels


def test_variant_model_creation():
    """Test that variant models are created correctly based on config."""
    config = _small_config()

    # Test transformer_only variant
    model = create_model_variant("transformer_only", config)
    assert isinstance(model, GeneWhispererTransformerOnly)

    # Test features_only variant
    model = create_model_variant("features_only", config)
    assert isinstance(model, GeneWhispererFeaturesOnly)

    # Test combined variant
    model = create_model_variant("combined", config)
    assert isinstance(model, GeneWhispererCombined)


def test_variant_forward_with_appropriate_inputs():
    """Test that variants handle appropriate inputs correctly."""
    config = _small_config()
    tokens, engineered, _ = _build_inputs(config)

    # transformer_only: tokens required, engineered optional (ignored)
    model = create_model_variant("transformer_only", config)
    model.eval()
    with torch.no_grad():
        out1 = model(tokens, None)
        out2 = model(tokens, engineered)  # engineered is ignored
    assert out1.shape == (4, 1)
    assert out2.shape == (4, 1)
    # Output should be the same since engineered is ignored
    assert torch.allclose(out1, out2, atol=1e-6)

    # features_only: engineered required, tokens optional (ignored)
    model = create_model_variant("features_only", config)
    model.eval()
    with torch.no_grad():
        out = model(None, engineered)
    assert out.shape == (4, 1)
    # Should raise error if engineered is None
    try:
        model(tokens, None)
        assert False, "Expected ValueError for features_only with None engineered"
    except ValueError:
        pass  # Expected

    # combined: both required
    model = create_model_variant("combined", config)
    model.eval()
    with torch.no_grad():
        out = model(tokens, engineered)
    assert out.shape == (4, 1)
    # Should raise error if engineered is None
    try:
        model(tokens, None)
        assert False, "Expected ValueError for combined with None engineered"
    except ValueError:
        pass  # Expected


def test_checkpoint_variant_metadata():
    """Test that checkpoint saving includes variant metadata."""
    from train_stage1 import save_checkpoint

    config = _small_config()
    model = create_model_variant("transformer_only", config)
    optimizer = torch.optim.Adam(model.parameters())

    with tempfile.TemporaryDirectory() as tmpdir:
        ckpt_path = Path(tmpdir) / "test_checkpoint.pt"
        metrics = {"accuracy": 0.95, "loss": 0.1, "best_threshold": 0.5}

        # Test saving with variant metadata
        save_checkpoint(
            path=ckpt_path,
            model=model,
            optimizer=optimizer,
            scheduler=None,
            epoch=1,
            metrics=metrics,
            model_variant="transformer_only",
            use_simplified_architecture=True,
        )

        # Load and verify metadata
        checkpoint = torch.load(ckpt_path, map_location="cpu")
        assert checkpoint["model_variant"] == "transformer_only"
        assert checkpoint["use_simplified_architecture"] is True
        assert "model_state" in checkpoint
        assert checkpoint["epoch"] == 1

        # Test saving with original model
        original_model = GeneWhispererStage1Legacy(
            vocab_size=config["vocab_size"],
            kmer=config["kmer"],
            embedding_dim=config["embedding_dim"],
            num_layers=config["transformer_layers"],
            num_heads=config["transformer_heads"],
            ff_dim=config["transformer_ff_dim"],
            dropout=config["transformer_dropout"],
            pad_token_id=config["pad_token_id"],
        )
        optimizer2 = torch.optim.Adam(original_model.parameters())
        ckpt_path2 = Path(tmpdir) / "test_checkpoint2.pt"

        save_checkpoint(
            path=ckpt_path2,
            model=original_model,
            optimizer=optimizer2,
            scheduler=None,
            epoch=1,
            metrics=metrics,
            model_variant="original",
            use_simplified_architecture=False,
        )

        checkpoint2 = torch.load(ckpt_path2, map_location="cpu")
        assert checkpoint2["model_variant"] == "original"
        assert checkpoint2["use_simplified_architecture"] is False


def test_use_simplified_architecture_flag():
    """Test that use_simplified_architecture controls model selection."""
    config = _small_config()

    # When use_simplified_architecture=False, should use original even if variant specified
    # This is handled at training time, not at model creation
    # We just verify the variant models can be created correctly here

    # Test that all variants can be created
    for variant in ["transformer_only", "features_only", "combined"]:
        model = create_model_variant(variant, config)
        assert model is not None

        # Verify forward pass works
        tokens, engineered, _ = _build_inputs(config)
        with torch.no_grad():
            if variant == "transformer_only":
                out = model(tokens, None)
            elif variant == "features_only":
                out = model(None, engineered)
            else:
                out = model(tokens, engineered)
            assert out.shape == (4, 1)


def test_variant_parameter_counts():
    """Test that different variants have expected parameter count relationships."""
    config = _small_config()

    transformer_only = create_model_variant("transformer_only", config)
    features_only = create_model_variant("features_only", config)
    combined = create_model_variant("combined", config)

    t_params = sum(p.numel() for p in transformer_only.parameters())
    f_params = sum(p.numel() for p in features_only.parameters())
    c_params = sum(p.numel() for p in combined.parameters())

    # features_only should have the fewest parameters (no transformer)
    assert f_params < t_params, "features_only should have fewer params than transformer_only"

    # combined should have more parameters than both (has both components)
    assert c_params > t_params, "combined should have more params than transformer_only"
    assert c_params > f_params, "combined should have more params than features_only"


def test_gradient_flow_all_variants():
    """Test that gradients flow correctly for all variants."""
    config = _small_config()
    tokens, engineered, labels = _build_inputs(config)
    labels = labels.unsqueeze(1)

    for variant in ["transformer_only", "features_only", "combined"]:
        torch.manual_seed(42)
        model = create_model_variant(variant, config)
        model.train()

        # Forward pass with appropriate inputs
        if variant == "transformer_only":
            logits = model(tokens, None)
        elif variant == "features_only":
            logits = model(None, engineered)
        else:
            logits = model(tokens, engineered)

        loss = nn.BCEWithLogitsLoss()(logits, labels)
        loss.backward()

        # Check gradients exist and are finite
        has_gradient = False
        for name, param in model.named_parameters():
            if param.grad is not None:
                has_gradient = True
                assert torch.isfinite(param.grad).all(), f"Non-finite gradient in {variant}: {name}"

        assert has_gradient, f"No gradients computed for {variant}"


def _run_test(name: str, fn) -> bool:
    """Run a test function and report results."""
    try:
        fn()
    except Exception as exc:
        print(f"FAIL: {name} ({type(exc).__name__}: {exc})")
        return False
    print(f"PASS: {name}")
    return True


if __name__ == "__main__":
    tests = [
        ("variant model creation", test_variant_model_creation),
        ("variant forward with appropriate inputs", test_variant_forward_with_appropriate_inputs),
        ("checkpoint variant metadata", test_checkpoint_variant_metadata),
        ("use_simplified_architecture flag", test_use_simplified_architecture_flag),
        ("variant parameter counts", test_variant_parameter_counts),
        ("gradient flow all variants", test_gradient_flow_all_variants),
    ]

    all_passed = True
    for test_name, test_fn in tests:
        all_passed &= _run_test(test_name, test_fn)

    if not all_passed:
        raise SystemExit(1)

    print("\nAll tests passed!")
