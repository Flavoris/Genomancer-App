"""Tests for GeneWhispererV6 modern architecture.

Tests the V6 model which integrates:
- Rotary Position Embeddings (RoPE)
- SwiGLU FFN
- RMSNorm
- Optional positional motif bias
- Optional Grouped Query Attention (GQA)
"""

import pytest
import torch

import sys
from pathlib import Path

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from model_v6 import (
    RMSNorm,
    SwiGLU,
    ModernTransformerBlock,
    GeneWhispererV6,
    create_v6_model,
)
from model import create_model_variant


class TestModernTransformerBlock:
    """Tests for ModernTransformerBlock."""

    def test_output_shape(self):
        """Block preserves input shape."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4)
        x = torch.randn(2, 10, 64)
        out = block(x)
        assert out.shape == x.shape

    def test_with_padding_mask(self):
        """Block handles padding mask correctly."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4)
        x = torch.randn(2, 10, 64)
        mask = torch.zeros(2, 10, dtype=torch.bool)
        mask[0, 8:] = True  # Last 2 positions are padding
        mask[1, 5:] = True  # Last 5 positions are padding

        out = block(x, mask=mask)
        assert out.shape == x.shape

    def test_rope_enabled(self):
        """Block works with RoPE enabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_rope=True)
        assert block.rotary is not None
        x = torch.randn(2, 10, 64)
        out = block(x)
        assert out.shape == x.shape

    def test_rope_disabled(self):
        """Block works with RoPE disabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_rope=False)
        assert block.rotary is None
        x = torch.randn(2, 10, 64)
        out = block(x)
        assert out.shape == x.shape

    def test_swiglu_enabled(self):
        """Block uses SwiGLU FFN when enabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_swiglu=True)
        assert isinstance(block.ffn, SwiGLU)

    def test_swiglu_disabled(self):
        """Block uses standard FFN when SwiGLU disabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_swiglu=False)
        assert not isinstance(block.ffn, SwiGLU)
        assert isinstance(block.ffn, torch.nn.Sequential)

    def test_rmsnorm_enabled(self):
        """Block uses RMSNorm when enabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_rmsnorm=True)
        assert isinstance(block.norm1, RMSNorm)
        assert isinstance(block.norm2, RMSNorm)

    def test_rmsnorm_disabled(self):
        """Block uses LayerNorm when RMSNorm disabled."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4, use_rmsnorm=False)
        assert isinstance(block.norm1, torch.nn.LayerNorm)
        assert isinstance(block.norm2, torch.nn.LayerNorm)

    def test_positional_bias_enabled(self):
        """Block creates motif bias when enabled."""
        block = ModernTransformerBlock(
            embed_dim=64, num_heads=4, use_positional_bias=True
        )
        assert block.motif_bias is not None

    def test_positional_bias_disabled(self):
        """Block has no motif bias when disabled."""
        block = ModernTransformerBlock(
            embed_dim=64, num_heads=4, use_positional_bias=False
        )
        assert block.motif_bias is None

    def test_gqa_same_heads(self):
        """GQA with same Q and KV heads (standard attention)."""
        block = ModernTransformerBlock(
            embed_dim=64, num_heads=4, num_kv_heads=4
        )
        assert block.kv_groups == 1
        x = torch.randn(2, 10, 64)
        out = block(x)
        assert out.shape == x.shape

    def test_gqa_fewer_kv_heads(self):
        """GQA with fewer KV heads (grouped query attention)."""
        block = ModernTransformerBlock(
            embed_dim=64, num_heads=8, num_kv_heads=2
        )
        assert block.kv_groups == 4
        x = torch.randn(2, 10, 64)
        out = block(x)
        assert out.shape == x.shape

    def test_gqa_invalid_configuration(self):
        """GQA raises error for invalid head configuration."""
        with pytest.raises(ValueError, match="divisible"):
            ModernTransformerBlock(embed_dim=64, num_heads=8, num_kv_heads=3)

    def test_gradient_flow(self):
        """Gradients flow through the block."""
        block = ModernTransformerBlock(embed_dim=64, num_heads=4)
        x = torch.randn(2, 10, 64, requires_grad=True)
        out = block(x)
        loss = out.sum()
        loss.backward()
        assert x.grad is not None


class TestGeneWhispererV6:
    """Tests for GeneWhispererV6 model."""

    def test_output_shape_without_features(self):
        """V6 produces correct output shape without engineered features."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            use_engineered_features=False,
        )
        tokens = torch.randint(0, 100, (2, 10))
        out = model(tokens)
        assert out.shape == (2, 1)

    def test_output_shape_with_features(self):
        """V6 produces correct output shape with engineered features."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            engineered_dim=32,
            use_engineered_features=True,
        )
        tokens = torch.randint(0, 100, (2, 10))
        features = torch.randn(2, 32)
        out = model(tokens, features)
        assert out.shape == (2, 1)

    def test_return_logits(self):
        """V6 returns both probs and logits when requested."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            use_engineered_features=False,
        )
        tokens = torch.randint(0, 100, (2, 10))
        probs, logits = model(tokens, return_logits=True)
        assert probs.shape == (2, 1)
        assert logits.shape == (2, 1)
        assert torch.all((probs >= 0) & (probs <= 1))

    def test_padding_handling(self):
        """V6 handles padding correctly."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            pad_token_id=99,
            use_engineered_features=False,
        )
        # Create tokens with padding
        tokens = torch.randint(0, 98, (2, 10))
        tokens[0, 8:] = 99  # Pad last 2 tokens
        tokens[1, 5:] = 99  # Pad last 5 tokens

        out = model(tokens)
        assert out.shape == (2, 1)
        assert not torch.isnan(out).any()

    def test_all_modern_features_enabled(self):
        """V6 works with all modern features enabled."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            use_rope=True,
            use_swiglu=True,
            use_rmsnorm=True,
            use_positional_bias=True,
            use_engineered_features=True,
            engineered_dim=32,
        )
        tokens = torch.randint(0, 100, (2, 10))
        features = torch.randn(2, 32)
        out = model(tokens, features)
        assert out.shape == (2, 1)

    def test_gqa_configuration(self):
        """V6 works with GQA (fewer KV heads)."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=8,
            num_kv_heads=2,  # GQA with 4 groups
            use_engineered_features=False,
        )
        tokens = torch.randint(0, 100, (2, 10))
        out = model(tokens)
        assert out.shape == (2, 1)

    def test_encoder_property(self):
        """V6 encoder property returns self for compatibility."""
        model = GeneWhispererV6(vocab_size=100, embed_dim=64, num_layers=2, num_heads=4)
        assert model.encoder is model

    def test_gradient_flow(self):
        """Gradients flow through the entire V6 model."""
        model = GeneWhispererV6(
            vocab_size=100,
            embed_dim=64,
            num_layers=2,
            num_heads=4,
            use_engineered_features=False,
        )
        tokens = torch.randint(0, 100, (2, 10))
        out = model(tokens)
        loss = out.sum()
        loss.backward()

        # Check gradients exist for key components
        assert model.embed.weight.grad is not None
        assert model.classifier[0].weight.grad is not None


class TestCreateV6Model:
    """Tests for create_v6_model factory function."""

    def test_basic_config(self):
        """Factory creates model from basic config."""
        config = {
            "vocab_size": 4096,
            "embedding_dim": 384,
            "transformer_layers": 12,
            "transformer_heads": 12,
            "max_token_len": 24,
        }
        model = create_v6_model(config)
        assert isinstance(model, GeneWhispererV6)
        assert model.embed_dim == 384

    def test_with_gqa_config(self):
        """Factory creates GQA model when configured."""
        config = {
            "vocab_size": 4096,
            "embedding_dim": 384,
            "transformer_layers": 12,
            "transformer_heads": 12,
            "use_gqa": True,
            "num_kv_heads": 4,
        }
        model = create_v6_model(config)
        assert isinstance(model, GeneWhispererV6)
        # Check first layer has GQA
        assert model.layers[0].num_kv_heads == 4

    def test_with_disabled_features(self):
        """Factory creates model with features disabled."""
        config = {
            "vocab_size": 100,
            "embedding_dim": 64,
            "transformer_layers": 2,
            "transformer_heads": 4,
            "use_rope": False,
            "use_swiglu": False,
            "use_rmsnorm": False,
            "use_positional_bias": False,
        }
        model = create_v6_model(config)
        assert isinstance(model, GeneWhispererV6)
        # Check first layer settings
        assert model.layers[0].rotary is None
        assert not isinstance(model.layers[0].ffn, SwiGLU)

    def test_with_motif_regions(self):
        """Factory creates model with custom motif regions."""
        config = {
            "vocab_size": 100,
            "embedding_dim": 64,
            "transformer_layers": 2,
            "transformer_heads": 4,
            "use_positional_bias": True,
            "motif_regions": {
                "tata_box": [45, 55],
                "tss": [75, 81],
            },
        }
        model = create_v6_model(config)
        # First layer should have motif bias
        assert model.layers[0].motif_bias is not None


class TestModelVariantFactory:
    """Tests for create_model_variant with v6 variant."""

    def test_create_v6_variant(self):
        """Factory creates V6 model with 'v6' variant."""
        config = {
            "vocab_size": 100,
            "embedding_dim": 64,
            "transformer_layers": 2,
            "transformer_heads": 4,
        }
        model = create_model_variant("v6", config)
        assert isinstance(model, GeneWhispererV6)

    def test_create_v6_case_insensitive(self):
        """Factory handles case-insensitive variant names."""
        config = {
            "vocab_size": 100,
            "embedding_dim": 64,
            "transformer_layers": 2,
            "transformer_heads": 4,
        }
        model1 = create_model_variant("V6", config)
        model2 = create_model_variant("v6", config)
        model3 = create_model_variant(" v6 ", config)  # With whitespace
        assert all(isinstance(m, GeneWhispererV6) for m in [model1, model2, model3])


class TestV6ParameterCount:
    """Tests for V6 model parameter efficiency."""

    def test_gqa_reduces_parameters(self):
        """GQA reduces parameter count compared to standard attention."""
        config_base = {
            "vocab_size": 4096,
            "embedding_dim": 384,
            "transformer_layers": 12,
            "transformer_heads": 12,
            "use_gqa": False,
        }
        config_gqa = {
            **config_base,
            "use_gqa": True,
            "num_kv_heads": 4,
        }

        model_base = create_v6_model(config_base)
        model_gqa = create_v6_model(config_gqa)

        params_base = sum(p.numel() for p in model_base.parameters())
        params_gqa = sum(p.numel() for p in model_gqa.parameters())

        # GQA should have fewer parameters
        assert params_gqa < params_base

    def test_reasonable_parameter_count(self):
        """V6 with default config has reasonable parameter count."""
        config = {
            "vocab_size": 4096,
            "embedding_dim": 384,
            "transformer_layers": 12,
            "transformer_heads": 12,
            "engineered_dim": 288,
            "stage1_use_engineered_features": True,
        }
        model = create_v6_model(config)
        params = sum(p.numel() for p in model.parameters())

        # Should be in reasonable range (15-30M parameters)
        assert 10_000_000 < params < 50_000_000
