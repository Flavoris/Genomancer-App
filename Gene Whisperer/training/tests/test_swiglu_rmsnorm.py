"""Tests for SwiGLU and RMSNorm modules.

These modules are modern transformer components used in LLaMA, PaLM, and Mistral.
"""

import pytest
import torch

import sys
from pathlib import Path

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from model import RMSNorm, SwiGLU, PreNormTransformerLayer, DNAEncoder


class TestRMSNorm:
    """Tests for RMSNorm module."""

    def test_output_shape(self):
        """RMSNorm preserves input shape."""
        norm = RMSNorm(dim=64)
        x = torch.randn(2, 10, 64)
        out = norm(x)
        assert out.shape == x.shape

    def test_normalization(self):
        """RMSNorm normalizes by RMS, not mean."""
        norm = RMSNorm(dim=64)
        x = torch.randn(2, 10, 64) * 10  # Large values
        out = norm(x)

        # Compute RMS manually
        rms = torch.sqrt(x.pow(2).mean(dim=-1, keepdim=True) + norm.eps)
        expected = (x / rms) * norm.weight

        assert torch.allclose(out, expected, atol=1e-5)

    def test_weight_initialization(self):
        """RMSNorm weight is initialized to ones."""
        norm = RMSNorm(dim=128)
        assert torch.allclose(norm.weight, torch.ones(128))

    def test_gradient_flow(self):
        """Gradients flow through RMSNorm."""
        norm = RMSNorm(dim=32)
        x = torch.randn(4, 8, 32, requires_grad=True)
        out = norm(x)
        loss = out.sum()
        loss.backward()

        assert x.grad is not None
        assert norm.weight.grad is not None

    def test_different_dimensions(self):
        """RMSNorm works with different dimensions."""
        for dim in [32, 64, 128, 256, 384]:
            norm = RMSNorm(dim=dim)
            x = torch.randn(2, 10, dim)
            out = norm(x)
            assert out.shape == (2, 10, dim)

    def test_batch_independence(self):
        """Each batch element is normalized independently."""
        norm = RMSNorm(dim=64)
        x1 = torch.randn(1, 10, 64)
        x2 = torch.randn(1, 10, 64) * 100  # Very different scale

        out1 = norm(x1)
        out2 = norm(x2)

        # RMS of outputs should be similar (both normalized)
        rms1 = torch.sqrt(out1.pow(2).mean())
        rms2 = torch.sqrt(out2.pow(2).mean())

        # Both should be around 1.0 (normalized)
        assert abs(rms1.item() - rms2.item()) < 0.5


class TestSwiGLU:
    """Tests for SwiGLU module."""

    def test_output_shape(self):
        """SwiGLU preserves input/output dimensions."""
        swiglu = SwiGLU(in_features=64, out_features=64)
        x = torch.randn(2, 10, 64)
        out = swiglu(x)
        assert out.shape == x.shape

    def test_custom_hidden_dim(self):
        """SwiGLU respects custom hidden dimension."""
        swiglu = SwiGLU(in_features=64, hidden_features=128, out_features=64)
        # Check internal dimensions
        assert swiglu.w1.out_features == 128
        assert swiglu.w2.in_features == 128
        assert swiglu.w3.out_features == 128

    def test_default_hidden_dim(self):
        """SwiGLU default hidden is 8/3 * in_features, rounded to 64."""
        swiglu = SwiGLU(in_features=384)
        expected = ((int(384 * 8 / 3) + 63) // 64) * 64
        assert swiglu.w1.out_features == expected

    def test_no_bias_by_default(self):
        """SwiGLU has no bias by default."""
        swiglu = SwiGLU(in_features=64)
        assert swiglu.w1.bias is None
        assert swiglu.w2.bias is None
        assert swiglu.w3.bias is None

    def test_with_bias(self):
        """SwiGLU can use bias if specified."""
        swiglu = SwiGLU(in_features=64, bias=True)
        assert swiglu.w1.bias is not None
        assert swiglu.w2.bias is not None
        assert swiglu.w3.bias is not None

    def test_gradient_flow(self):
        """Gradients flow through SwiGLU."""
        swiglu = SwiGLU(in_features=32, out_features=32)
        x = torch.randn(4, 8, 32, requires_grad=True)
        out = swiglu(x)
        loss = out.sum()
        loss.backward()

        assert x.grad is not None
        assert swiglu.w1.weight.grad is not None
        assert swiglu.w2.weight.grad is not None
        assert swiglu.w3.weight.grad is not None

    def test_swiglu_formula(self):
        """SwiGLU follows the formula: (Swish(W1*x) * W3*x) * W2."""
        torch.manual_seed(42)
        swiglu = SwiGLU(in_features=32, hidden_features=64, out_features=32, dropout=0.0)
        swiglu.eval()

        x = torch.randn(2, 5, 32)
        out = swiglu(x)

        # Manual computation
        w1_out = torch.nn.functional.silu(swiglu.w1(x))
        w3_out = swiglu.w3(x)
        expected = swiglu.w2(w1_out * w3_out)

        assert torch.allclose(out, expected, atol=1e-5)

    def test_dropout(self):
        """SwiGLU applies dropout during training."""
        swiglu = SwiGLU(in_features=64, dropout=0.5)
        swiglu.train()

        x = torch.randn(4, 10, 64)
        out1 = swiglu(x)
        out2 = swiglu(x)

        # Outputs should differ due to dropout
        assert not torch.allclose(out1, out2)

    def test_no_dropout_in_eval(self):
        """SwiGLU has no dropout during eval."""
        swiglu = SwiGLU(in_features=64, dropout=0.5)
        swiglu.eval()

        x = torch.randn(4, 10, 64)
        out1 = swiglu(x)
        out2 = swiglu(x)

        # Outputs should be identical
        assert torch.allclose(out1, out2)


class TestPreNormTransformerLayerWithSwiGLU:
    """Tests for PreNormTransformerLayer with SwiGLU and RMSNorm."""

    def test_swiglu_ffn_type(self):
        """PreNormTransformerLayer uses SwiGLU when ffn_type='swiglu'."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,
            ffn_type="swiglu",
        )
        assert isinstance(layer.ffn, SwiGLU)

    def test_rmsnorm_type(self):
        """PreNormTransformerLayer uses RMSNorm when norm_type='rmsnorm'."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,
            norm_type="rmsnorm",
        )
        assert isinstance(layer.norm1, RMSNorm)
        assert isinstance(layer.norm2, RMSNorm)

    def test_layernorm_default(self):
        """PreNormTransformerLayer uses LayerNorm by default."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,
        )
        assert isinstance(layer.norm1, torch.nn.LayerNorm)
        assert isinstance(layer.norm2, torch.nn.LayerNorm)

    def test_combined_swiglu_rmsnorm(self):
        """SwiGLU and RMSNorm work together."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,
            ffn_type="swiglu",
            norm_type="rmsnorm",
        )

        x = torch.randn(2, 10, 64)
        out = layer(x)

        assert out.shape == x.shape
        assert isinstance(layer.ffn, SwiGLU)
        assert isinstance(layer.norm1, RMSNorm)

    def test_ffn_mult_parameter(self):
        """ffn_mult controls SwiGLU hidden dimension."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,  # This is ignored when using swiglu
            ffn_type="swiglu",
            ffn_mult=3.0,  # 3x expansion
        )
        # Hidden should be ceil(64 * 3.0 / 64) * 64 = 192
        expected = ((int(64 * 3.0) + 63) // 64) * 64
        assert layer.ffn.w1.out_features == expected

    def test_gradient_flow_combined(self):
        """Gradients flow through combined SwiGLU + RMSNorm."""
        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=256,
            ffn_type="swiglu",
            norm_type="rmsnorm",
        )

        x = torch.randn(2, 10, 64, requires_grad=True)
        out = layer(x)
        loss = out.sum()
        loss.backward()

        assert x.grad is not None


class TestDNAEncoderWithSwiGLU:
    """Tests for DNAEncoder with SwiGLU and RMSNorm."""

    def test_encoder_with_swiglu_rmsnorm(self):
        """DNAEncoder works with SwiGLU and RMSNorm."""
        encoder = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=64,
            num_layers=2,
            num_heads=4,
            ff_dim=256,
            ffn_type="swiglu",
            norm_type="rmsnorm",
        )

        tokens = torch.randint(0, 67, (2, 20))
        out = encoder(tokens)

        assert out.shape == (2, 20, 64)

    def test_encoder_layers_use_correct_types(self):
        """All encoder layers use SwiGLU and RMSNorm when specified."""
        encoder = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=64,
            num_layers=4,
            num_heads=4,
            ff_dim=256,
            ffn_type="swiglu",
            norm_type="rmsnorm",
        )

        for layer in encoder.layers:
            assert isinstance(layer.ffn, SwiGLU)
            assert isinstance(layer.norm1, RMSNorm)
            assert isinstance(layer.norm2, RMSNorm)

    def test_backward_compatibility(self):
        """DNAEncoder still works with legacy use_glu_ffn parameter."""
        encoder = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=64,
            num_layers=2,
            num_heads=4,
            ff_dim=256,
            use_glu_ffn=True,
            glu_activation="silu",
        )

        tokens = torch.randint(0, 67, (2, 20))
        out = encoder(tokens)

        assert out.shape == (2, 20, 64)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
