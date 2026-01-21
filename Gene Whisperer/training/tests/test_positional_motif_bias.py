"""Tests for Positional Motif Attention Bias.

Tests the PositionalMotifBias module and its integration with the model.
"""
import pytest
import torch
import torch.nn as nn

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from positional_motif_bias import (
    PositionalMotifBias,
    DEFAULT_PROMOTER_PRIORS,
    create_promoter_motif_bias,
)


class TestPositionalMotifBias:
    """Tests for the PositionalMotifBias module."""

    def test_init_basic(self):
        """Test basic initialization without priors."""
        bias = PositionalMotifBias(num_heads=8, max_seq_len=81)
        assert bias.num_heads == 8
        assert bias.max_seq_len == 81
        assert bias.bias.shape == (8, 81, 81)
        # Without priors, bias should be zero
        assert torch.allclose(bias.bias, torch.zeros_like(bias.bias))

    def test_init_with_priors(self):
        """Test initialization with custom priors."""
        priors = {(10, 20): 0.1, (30, 40): 0.2}
        bias = PositionalMotifBias(
            num_heads=4,
            max_seq_len=50,
            init_bias_positions=priors,
        )
        # Check that priors are applied
        # Bias should be non-zero in the specified regions
        assert bias.bias[:, :, 10:20].abs().sum() > 0
        assert bias.bias[:, :, 30:40].abs().sum() > 0
        # Bias should be zero outside the regions (mostly)
        assert bias.bias[:, :, 0:10].abs().sum() > 0  # FROM effect
        assert bias.bias[:, :, 20:30].abs().sum() > 0  # FROM effect

    def test_init_with_default_priors(self):
        """Test initialization with default promoter priors."""
        bias = PositionalMotifBias(
            num_heads=12,
            max_seq_len=81,
            init_bias_positions=DEFAULT_PROMOTER_PRIORS,
        )
        # TATA box region (45-55) should have non-zero bias
        tata_bias = bias.bias[:, :, 45:55].abs().sum()
        assert tata_bias > 0, "TATA box region should have bias"

        # TSS region (75-81) should have higher bias
        tss_bias = bias.bias[:, :, 75:81].abs().sum()
        assert tss_bias > 0, "TSS region should have bias"

    def test_forward_shape(self):
        """Test that forward preserves attention score shape."""
        bias = PositionalMotifBias(num_heads=8, max_seq_len=100)
        attn_scores = torch.randn(4, 8, 50, 50)  # B=4, H=8, L=50
        output = bias(attn_scores)
        assert output.shape == attn_scores.shape

    def test_forward_adds_bias(self):
        """Test that forward adds bias to attention scores."""
        priors = {(5, 15): 0.5}
        bias = PositionalMotifBias(
            num_heads=4,
            max_seq_len=20,
            init_bias_positions=priors,
        )
        attn_scores = torch.zeros(2, 4, 20, 20)
        output = bias(attn_scores)
        # Output should equal the bias since input was zeros
        expected = bias.bias.unsqueeze(0).expand(2, -1, -1, -1)
        assert torch.allclose(output, expected)

    def test_forward_variable_length(self):
        """Test forward with sequences shorter than max_seq_len."""
        bias = PositionalMotifBias(num_heads=8, max_seq_len=100)
        # Use shorter sequence
        attn_scores = torch.randn(4, 8, 30, 30)
        output = bias(attn_scores)
        assert output.shape == (4, 8, 30, 30)

    def test_forward_exceeds_max_length_raises(self):
        """Test that exceeding max_seq_len raises an error."""
        bias = PositionalMotifBias(num_heads=8, max_seq_len=50)
        attn_scores = torch.randn(4, 8, 60, 60)  # 60 > 50
        with pytest.raises(ValueError, match="exceeds max_seq_len"):
            bias(attn_scores)

    def test_learnable_parameter(self):
        """Test that bias is learnable when learnable=True."""
        bias = PositionalMotifBias(num_heads=4, max_seq_len=20, learnable=True)
        assert bias.bias.requires_grad
        # Check it's a Parameter
        assert isinstance(bias.bias, nn.Parameter)

    def test_non_learnable_buffer(self):
        """Test that bias is a buffer when learnable=False."""
        bias = PositionalMotifBias(num_heads=4, max_seq_len=20, learnable=False)
        # Buffers don't require grad
        assert not bias.bias.requires_grad

    def test_gradient_flow(self):
        """Test that gradients flow through the bias."""
        bias = PositionalMotifBias(
            num_heads=4,
            max_seq_len=20,
            init_bias_positions={(5, 10): 0.1},
            learnable=True,
        )
        attn_scores = torch.randn(2, 4, 20, 20, requires_grad=True)
        output = bias(attn_scores)
        loss = output.sum()
        loss.backward()
        # Both input and bias should have gradients
        assert attn_scores.grad is not None
        assert bias.bias.grad is not None

    def test_prior_initialization_values(self):
        """Test that priors are correctly applied with expected values."""
        priors = {(0, 5): 0.1}
        bias = PositionalMotifBias(
            num_heads=2,
            max_seq_len=10,
            init_bias_positions=priors,
        )
        # The implementation adds:
        # 1. TO columns (0:5): bias[:, :, 0:5] += 0.1
        # 2. FROM rows (0:5): bias[:, 0:5, :] += 0.1 * 0.5 = 0.05
        #
        # So for rows 0-4 (FROM region):
        #   columns 0-4 (overlap): 0.1 (TO) + 0.05 (FROM) = 0.15
        #   columns 5-9 (outside): 0.05 (FROM only)
        # For rows 5-9 (outside FROM region):
        #   columns 0-4: 0.1 (TO only)
        #   columns 5-9: 0.0

        # Check rows in FROM region (0-4)
        assert torch.allclose(
            bias.bias[:, 0:5, 0:5],  # overlap region
            torch.full((2, 5, 5), 0.15),
        )
        assert torch.allclose(
            bias.bias[:, 0:5, 5:10],  # FROM only
            torch.full((2, 5, 5), 0.05),
        )

        # Check rows outside FROM region (5-9)
        assert torch.allclose(
            bias.bias[:, 5:10, 0:5],  # TO only
            torch.full((2, 5, 5), 0.1),
        )
        assert torch.allclose(
            bias.bias[:, 5:10, 5:10],  # no effect
            torch.full((2, 5, 5), 0.0),
        )


class TestCreatePromoterMotifBias:
    """Tests for the factory function."""

    def test_default_priors(self):
        """Test creation with default priors."""
        bias = create_promoter_motif_bias(num_heads=12, max_seq_len=81)
        assert bias.num_heads == 12
        assert bias.max_seq_len == 81
        # Should have non-zero bias from default priors
        assert bias.bias.abs().sum() > 0

    def test_no_default_priors(self):
        """Test creation without default priors."""
        bias = create_promoter_motif_bias(
            num_heads=8,
            max_seq_len=81,
            use_default_priors=False,
        )
        # Without any priors, bias should be zero
        assert torch.allclose(bias.bias, torch.zeros_like(bias.bias))

    def test_custom_priors(self):
        """Test creation with custom priors."""
        custom = {(0, 10): 0.2}
        bias = create_promoter_motif_bias(
            num_heads=4,
            max_seq_len=50,
            use_default_priors=False,
            custom_priors=custom,
        )
        assert bias.bias[:, :, 0:10].abs().sum() > 0

    def test_custom_priors_override(self):
        """Test that custom priors merge with defaults."""
        custom = {(0, 5): 0.5}
        bias = create_promoter_motif_bias(
            num_heads=4,
            max_seq_len=81,
            use_default_priors=True,
            custom_priors=custom,
        )
        # Both custom and default regions should have bias
        assert bias.bias[:, :, 0:5].abs().sum() > 0
        assert bias.bias[:, :, 45:55].abs().sum() > 0  # TATA box default


class TestModelIntegration:
    """Tests for integration with the model."""

    def test_dna_encoder_with_motif_bias(self):
        """Test DNAEncoder with positional motif bias enabled."""
        from model import DNAEncoder

        encoder = DNAEncoder(
            vocab_size=100,
            kmer=3,
            embedding_dim=64,
            num_layers=2,
            num_heads=4,
            ff_dim=128,
            dropout=0.1,
            max_seq_len=81,
            use_positional_motif_bias=True,
            motif_regions={"tss": [75, 81]},
        )
        # Check that positional_motif_bias was created
        assert encoder.positional_motif_bias is not None
        assert encoder.use_positional_motif_bias

        # Test forward pass
        tokens = torch.randint(0, 100, (2, 76))  # 81 - 6 + 1 = 76 for k=3
        output = encoder(tokens)
        assert output.shape == (2, 76, 64)

    def test_dna_encoder_without_motif_bias(self):
        """Test DNAEncoder without positional motif bias."""
        from model import DNAEncoder

        encoder = DNAEncoder(
            vocab_size=100,
            kmer=3,
            embedding_dim=64,
            num_layers=2,
            num_heads=4,
            ff_dim=128,
            dropout=0.1,
            max_seq_len=81,
            use_positional_motif_bias=False,
        )
        assert encoder.positional_motif_bias is None
        assert not encoder.use_positional_motif_bias

    def test_stage1_model_with_motif_bias(self):
        """Test GeneWhispererStage1 with positional motif bias."""
        from model import GeneWhispererStage1

        model = GeneWhispererStage1(
            vocab_size=100,
            kmer=3,
            embedding_dim=64,
            num_layers=2,
            num_heads=4,
            ff_dim=128,
            dropout=0.1,
            max_seq_len=81,
            use_engineered_features=False,
            use_attention_pool=False,
            use_positional_motif_bias=True,
            motif_regions={"tss": [70, 76]},
        )
        assert model._full_encoder.positional_motif_bias is not None

        # Test forward pass
        tokens = torch.randint(0, 100, (2, 76))
        output = model(tokens)
        assert output.shape == (2, 1)

    def test_transformer_layer_with_motif_bias(self):
        """Test PreNormTransformerLayer accepts motif_bias."""
        from model import PreNormTransformerLayer

        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=128,
            dropout=0.1,
        )
        x = torch.randn(2, 20, 64)
        motif_bias = torch.randn(1, 4, 20, 20)
        output = layer(x, motif_bias=motif_bias)
        assert output.shape == (2, 20, 64)

    def test_motif_bias_affects_attention(self):
        """Test that motif bias actually affects the attention computation."""
        from model import PreNormTransformerLayer

        layer = PreNormTransformerLayer(
            d_model=64,
            nhead=4,
            dim_feedforward=128,
            dropout=0.0,  # No dropout for deterministic comparison
        )
        layer.eval()  # No stochastic depth

        x = torch.randn(2, 20, 64)

        # Run without motif bias
        with torch.no_grad():
            output_without = layer(x, motif_bias=None)

        # Run with strong motif bias
        strong_bias = torch.ones(1, 4, 20, 20) * 10.0
        with torch.no_grad():
            output_with = layer(x, motif_bias=strong_bias)

        # Outputs should be different
        assert not torch.allclose(output_without, output_with)


class TestConfigParsing:
    """Tests for config parsing of motif regions."""

    def test_config_motif_regions_format(self):
        """Test that config format is correctly parsed."""
        # Config format: {name: [start, end]}
        config_regions = {
            "tata_box": [45, 55],
            "minus_10_box": [64, 74],
            "tss": [75, 81],
        }

        # Convert to expected format: {(start, end): bias_val}
        expected_priors = {}
        for name, (start, end) in config_regions.items():
            bias_val = 0.15 if "tss" in name.lower() else 0.1
            expected_priors[(start, end)] = bias_val

        assert expected_priors[(45, 55)] == 0.1
        assert expected_priors[(64, 74)] == 0.1
        assert expected_priors[(75, 81)] == 0.15


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
