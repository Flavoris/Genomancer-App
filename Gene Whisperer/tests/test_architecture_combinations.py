"""Test all 4 combinations of architectural features (RelativePositionBias + GLU FFN)."""

import sys
from pathlib import Path

import torch

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "training"))

from model import GeneWhispererStage1


def test_architecture_combinations():
    """Test all 4 combinations of architectural features."""
    configs = [
        {"use_relative_position_bias": False, "use_glu_ffn": False},
        {"use_relative_position_bias": True, "use_glu_ffn": False},
        {"use_relative_position_bias": False, "use_glu_ffn": True},
        {"use_relative_position_bias": True, "use_glu_ffn": True},
    ]

    for cfg in configs:
        print(f"\nTesting config: {cfg}")

        model = GeneWhispererStage1(
            vocab_size=4099,
            kmer=6,
            embedding_dim=256,
            num_layers=2,
            num_heads=8,
            ff_dim=1024,
            engineered_dim=288,
            post_cnn_transformer_layers=2,
            **cfg
        )

        # Test forward pass
        tokens = torch.randint(0, 4096, (2, 76))
        features = torch.randn(2, 288)
        output = model(tokens, features)

        assert output.shape == (2, 1), f"Expected (2, 1), got {output.shape}"

        # Count parameters
        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

        print(f"  Output shape: {output.shape}")
        print(f"  Total params: {total_params:,}")
        print(f"  Trainable params: {trainable_params:,}")

        # Check that relative position bias is created when enabled
        if cfg["use_relative_position_bias"]:
            assert model._full_encoder.relative_position_bias is not None, \
                "DNAEncoder should have relative_position_bias"
            assert model.post_cnn_transformer.relative_position_bias is not None, \
                "PostCNNTransformerAdapter should have relative_position_bias"
        else:
            assert model._full_encoder.relative_position_bias is None, \
                "DNAEncoder should NOT have relative_position_bias"
            assert model.post_cnn_transformer.relative_position_bias is None, \
                "PostCNNTransformerAdapter should NOT have relative_position_bias"

        # Check that GLU FFN is used when enabled
        if cfg["use_glu_ffn"]:
            # Check that the FFN has gate_proj (GLU signature)
            layer = model.post_cnn_transformer.layers[0]
            assert hasattr(layer.ffn, "gate_proj"), "FFN should be GLUFFN with gate_proj"
        else:
            # Check that the FFN is a Sequential (standard FFN)
            layer = model.post_cnn_transformer.layers[0]
            assert isinstance(layer.ffn, torch.nn.Sequential), \
                "FFN should be standard Sequential"

        print(f"  Config {cfg} passed!")

    print("\n" + "=" * 60)
    print("All 4 architecture combinations passed!")
    print("=" * 60)


def test_backward_pass():
    """Test that gradients flow through all architecture variants."""
    configs = [
        {"use_relative_position_bias": False, "use_glu_ffn": False},
        {"use_relative_position_bias": True, "use_glu_ffn": True},
    ]

    for cfg in configs:
        print(f"\nTesting backward pass for: {cfg}")

        model = GeneWhispererStage1(
            vocab_size=4099,
            kmer=6,
            embedding_dim=256,
            num_layers=2,
            num_heads=8,
            ff_dim=1024,
            engineered_dim=288,
            post_cnn_transformer_layers=2,
            **cfg
        )

        tokens = torch.randint(0, 4096, (2, 76))
        features = torch.randn(2, 288)

        output = model(tokens, features)
        loss = output.sum()
        loss.backward()

        # Check that gradients exist
        has_grad = False
        for name, param in model.named_parameters():
            if param.grad is not None and param.grad.abs().sum() > 0:
                has_grad = True
                break

        assert has_grad, f"No gradients computed for config {cfg}"
        print(f"  Backward pass OK for {cfg}")

    print("\nBackward pass tests passed!")


def test_load_pretrained_compatibility():
    """Test that load_pretrained handles missing keys gracefully."""
    import tempfile
    import os

    print("\nTesting load_pretrained backward compatibility...")

    # Create a model WITHOUT new features and save it
    old_model = GeneWhispererStage1(
        vocab_size=4099,
        kmer=6,
        embedding_dim=256,
        num_layers=2,
        num_heads=8,
        ff_dim=1024,
        engineered_dim=288,
        post_cnn_transformer_layers=2,
        use_relative_position_bias=False,
        use_glu_ffn=False,
    )

    # Save the old model
    with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
        torch.save(old_model.state_dict(), f.name)
        checkpoint_path = f.name

    try:
        # Create a NEW model WITH new features
        new_model = GeneWhispererStage1(
            vocab_size=4099,
            kmer=6,
            embedding_dim=256,
            num_layers=2,
            num_heads=8,
            ff_dim=1024,
            engineered_dim=288,
            post_cnn_transformer_layers=2,
            use_relative_position_bias=True,
            use_glu_ffn=True,
        )

        # Load old checkpoint into new model (should work with strict=False)
        new_model.load_pretrained(checkpoint_path, strict=False)

        # Verify forward pass still works
        tokens = torch.randint(0, 4096, (2, 76))
        features = torch.randn(2, 288)
        output = new_model(tokens, features)

        assert output.shape == (2, 1), f"Expected (2, 1), got {output.shape}"
        print("  load_pretrained with missing keys: OK")

        # Test that strict=True raises error for missing keys
        new_model2 = GeneWhispererStage1(
            vocab_size=4099,
            kmer=6,
            embedding_dim=256,
            num_layers=2,
            num_heads=8,
            ff_dim=1024,
            engineered_dim=288,
            post_cnn_transformer_layers=2,
            use_relative_position_bias=True,
            use_glu_ffn=True,
        )

        try:
            new_model2.load_pretrained(checkpoint_path, strict=True)
            assert False, "Should have raised RuntimeError for strict=True"
        except RuntimeError as e:
            assert "strict=True" in str(e)
            print("  load_pretrained with strict=True correctly raises error: OK")

    finally:
        os.unlink(checkpoint_path)

    print("\nload_pretrained compatibility tests passed!")


if __name__ == "__main__":
    test_architecture_combinations()
    test_backward_pass()
    test_load_pretrained_compatibility()
    print("\n" + "=" * 60)
    print("ALL TESTS PASSED!")
    print("=" * 60)
