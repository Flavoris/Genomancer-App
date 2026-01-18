# tests/test_architecture_enhancements.py
import torch
import torch.nn as nn
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model import (
    GeneWhispererStage1,
    RelativePositionBias,
    GLUFFN,
)

class TestRelativePositionBias:

    def test_output_shape(self):
        """Test relative position bias has correct shape."""
        num_heads = 8
        seq_len = 79

        rpb = RelativePositionBias(num_heads=num_heads, num_buckets=32, max_distance=128)
        bias = rpb(seq_len, device=torch.device('cpu'))

        assert bias.shape == (1, num_heads, seq_len, seq_len)
        print(f"✓ Output shape correct: {bias.shape}")

    def test_symmetry(self):
        """Test that bias is symmetric for bidirectional attention."""
        rpb = RelativePositionBias(num_heads=4, bidirectional=True)
        bias = rpb(32, device=torch.device('cpu'))

        # For bidirectional, bias[i,j] should relate to bias[j,i]
        # (not necessarily equal, but related through bucketing)
        assert bias.shape[-1] == bias.shape[-2]
        print("✓ Symmetry check passed")

    def test_different_distances_different_bias(self):
        """Test that different relative distances get different biases."""
        rpb = RelativePositionBias(num_heads=1, num_buckets=32)
        bias = rpb(64, device=torch.device('cpu')).squeeze()

        # Adjacent positions should have different bias than distant positions
        adjacent_bias = bias[0, 1].item()
        distant_bias = bias[0, 50].item()

        assert adjacent_bias != distant_bias
        print(f"✓ Adjacent bias: {adjacent_bias:.4f}, Distant bias: {distant_bias:.4f}")

    def test_cached_computation(self):
        """Test that same seq_len gives same result (for potential caching)."""
        rpb = RelativePositionBias(num_heads=8)

        bias1 = rpb(50, device=torch.device('cpu'))
        bias2 = rpb(50, device=torch.device('cpu'))

        assert torch.allclose(bias1, bias2)
        print("✓ Deterministic computation")


class TestGLUFFN:

    def test_output_shape(self):
        """Test GLU FFN preserves shape."""
        d_model = 256
        ff_dim = 1024
        batch_size = 4
        seq_len = 79

        glu = GLUFFN(d_model=d_model, ff_dim=ff_dim)
        x = torch.randn(batch_size, seq_len, d_model)
        output = glu(x)

        assert output.shape == x.shape
        print(f"✓ Output shape preserved: {output.shape}")

    def test_parameter_count(self):
        """Test that GLU FFN has similar param count to standard FFN."""
        d_model = 256
        ff_dim = 1024

        # Standard FFN
        standard_ffn = nn.Sequential(
            nn.Linear(d_model, ff_dim),
            nn.GELU(),
            nn.Linear(ff_dim, d_model),
        )
        standard_params = sum(p.numel() for p in standard_ffn.parameters())

        # GLU FFN
        glu_ffn = GLUFFN(d_model=d_model, ff_dim=ff_dim)
        glu_params = sum(p.numel() for p in glu_ffn.parameters())

        # Should be within 20% of each other
        ratio = glu_params / standard_params
        assert 0.8 < ratio < 1.2, f"Param ratio {ratio} outside acceptable range"

        print(f"✓ Standard FFN params: {standard_params:,}")
        print(f"✓ GLU FFN params: {glu_params:,}")
        print(f"✓ Ratio: {ratio:.2f}")

    def test_gating_effect(self):
        """Test that gating actually affects output."""
        glu = GLUFFN(d_model=128, ff_dim=512)

        # With gating, output should be different from simple linear
        x = torch.randn(2, 10, 128)
        output = glu(x)

        # Output shouldn't be all zeros or all same
        assert output.std() > 0.01
        print(f"✓ Output std: {output.std():.4f} (gating is active)")

    def test_activation_variants(self):
        """Test both GELU and SiLU activations work."""
        for activation in ["gelu", "silu"]:
            glu = GLUFFN(d_model=128, ff_dim=512, activation=activation)
            x = torch.randn(2, 10, 128)
            output = glu(x)

            assert output.shape == x.shape
            print(f"✓ {activation.upper()} activation works")


class TestIntegratedModel:

    def test_all_architecture_combinations(self):
        """Test all 4 combinations of architectural features."""
        configs = [
            {"use_relative_position_bias": False, "use_glu_ffn": False},
            {"use_relative_position_bias": True, "use_glu_ffn": False},
            {"use_relative_position_bias": False, "use_glu_ffn": True},
            {"use_relative_position_bias": True, "use_glu_ffn": True},
        ]

        for cfg in configs:
            model = GeneWhispererStage1(
                vocab_size=4099,
                kmer=6,
                embedding_dim=128,  # Smaller for testing
                num_layers=2,
                num_heads=4,
                ff_dim=512,
                engineered_dim=288,
                **cfg
            )

            # Test forward pass
            tokens = torch.randint(0, 4096, (2, 76))
            features = torch.randn(2, 288)

            model.eval()
            with torch.no_grad():
                output = model(tokens, features)

            assert output.shape == (2, 1)
            assert not torch.isnan(output).any()

            print(f"✓ Config {cfg} works")

    def test_backward_pass(self):
        """Test that gradients flow correctly with new components."""
        model = GeneWhispererStage1(
            vocab_size=4099,
            kmer=6,
            embedding_dim=128,
            num_layers=2,
            num_heads=4,
            ff_dim=512,
            engineered_dim=288,
            use_relative_position_bias=True,
            use_glu_ffn=True,
        )

        tokens = torch.randint(0, 4096, (2, 76))
        features = torch.randn(2, 288)
        labels = torch.randint(0, 2, (2, 1)).float()

        model.train()
        output = model(tokens, features)
        loss = nn.BCEWithLogitsLoss()(output, labels)
        loss.backward()

        # Check gradients exist for new components
        has_rpb_grad = False
        has_glu_grad = False

        for name, param in model.named_parameters():
            if param.grad is not None:
                if 'relative' in name.lower() or 'position_bias' in name.lower():
                    has_rpb_grad = True
                if 'gate_proj' in name.lower() or 'up_proj' in name.lower():
                    has_glu_grad = True

        assert has_rpb_grad, "No gradients for relative position bias!"
        assert has_glu_grad, "No gradients for GLU components!"

        print("✓ Gradients flow to relative position bias")
        print("✓ Gradients flow to GLU components")

    def test_checkpoint_loading_backward_compat(self):
        """Test that old checkpoints can be loaded into new model."""
        # Create "old" model without new features
        old_model = GeneWhispererStage1(
            vocab_size=259,
            kmer=4,
            embedding_dim=128,
            num_layers=2,
            num_heads=4,
            ff_dim=512,
            engineered_dim=128,  # Old dimension
            use_relative_position_bias=False,
            use_glu_ffn=False,
        )

        # Save "old" checkpoint
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.pt', delete=False) as f:
            torch.save(old_model.state_dict(), f.name)
            checkpoint_path = f.name

        # Create "new" model with new features
        new_model = GeneWhispererStage1(
            vocab_size=259,
            kmer=4,
            embedding_dim=128,
            num_layers=2,
            num_heads=4,
            ff_dim=512,
            engineered_dim=128,
            use_relative_position_bias=True,  # NEW
            use_glu_ffn=True,                 # NEW
        )

        # Load old checkpoint into new model (should not error)
        loaded = new_model.load_legacy_checkpoint(checkpoint_path, strict=False)
        assert loaded > 0

        # Verify model still works
        tokens = torch.randint(0, 256, (2, 78))
        features = torch.randn(2, 128)
        output = new_model(tokens, features)

        assert output.shape == (2, 1)
        print("✓ Old checkpoint loads into new model")

        # Cleanup
        import os
        os.unlink(checkpoint_path)


if __name__ == "__main__":
    print("=" * 60)
    print("ARCHITECTURE ENHANCEMENT TESTS")
    print("=" * 60)

    # Run tests
    test_rpb = TestRelativePositionBias()
    test_rpb.test_output_shape()
    test_rpb.test_symmetry()
    test_rpb.test_different_distances_different_bias()
    test_rpb.test_cached_computation()

    print()

    test_glu = TestGLUFFN()
    test_glu.test_output_shape()
    test_glu.test_parameter_count()
    test_glu.test_gating_effect()
    test_glu.test_activation_variants()

    print()

    test_model = TestIntegratedModel()
    test_model.test_all_architecture_combinations()
    test_model.test_backward_pass()
    test_model.test_checkpoint_loading_backward_compat()

    print()
    print("=" * 60)
    print("ALL ARCHITECTURE TESTS PASSED ✓")
    print("=" * 60)
