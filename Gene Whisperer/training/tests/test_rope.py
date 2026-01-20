"""Tests for Rotary Position Embedding (RoPE) implementation."""
import pytest
import torch

from rope import RotaryEmbedding, RoPEAttention, apply_rotary_emb, rotate_half


class TestRotateHalf:
    """Tests for the rotate_half function."""

    def test_rotate_half_basic(self):
        """Test basic rotate_half operation."""
        x = torch.tensor([[1.0, 2.0, 3.0, 4.0]])
        result = rotate_half(x)
        expected = torch.tensor([[-3.0, -4.0, 1.0, 2.0]])
        assert torch.allclose(result, expected)

    def test_rotate_half_preserves_shape(self):
        """Test that rotate_half preserves tensor shape."""
        x = torch.randn(2, 4, 8, 64)
        result = rotate_half(x)
        assert result.shape == x.shape

    def test_rotate_half_double_application(self):
        """Test that applying rotate_half twice gives -x."""
        x = torch.randn(2, 4, 8, 64)
        result = rotate_half(rotate_half(x))
        assert torch.allclose(result, -x)


class TestRotaryEmbedding:
    """Tests for the RotaryEmbedding class."""

    def test_init(self):
        """Test RotaryEmbedding initialization."""
        rope = RotaryEmbedding(dim=64, max_seq_len=128)
        assert rope.dim == 64
        assert rope.max_seq_len == 128
        assert rope.inv_freq.shape == (32,)  # dim // 2

    def test_forward_shape(self):
        """Test that forward returns correct shapes."""
        rope = RotaryEmbedding(dim=64, max_seq_len=128)
        cos, sin = rope(seq_len=32)
        assert cos.shape == (1, 1, 32, 64)
        assert sin.shape == (1, 1, 32, 64)

    def test_forward_cache_extension(self):
        """Test that cache is extended for longer sequences."""
        rope = RotaryEmbedding(dim=64, max_seq_len=32)
        assert rope.max_seq_len == 32

        # Request longer sequence - cache should extend
        cos, sin = rope(seq_len=64)
        assert rope.max_seq_len == 64
        assert cos.shape == (1, 1, 64, 64)

    def test_cos_sin_normalized(self):
        """Test that cos and sin values are in valid range."""
        rope = RotaryEmbedding(dim=64, max_seq_len=128)
        cos, sin = rope(seq_len=32)
        assert cos.min() >= -1.0 and cos.max() <= 1.0
        assert sin.min() >= -1.0 and sin.max() <= 1.0


class TestApplyRotaryEmb:
    """Tests for the apply_rotary_emb function."""

    def test_apply_rotary_emb_shape(self):
        """Test that apply_rotary_emb preserves shapes."""
        batch_size, num_heads, seq_len, head_dim = 2, 8, 32, 64
        q = torch.randn(batch_size, num_heads, seq_len, head_dim)
        k = torch.randn(batch_size, num_heads, seq_len, head_dim)

        rope = RotaryEmbedding(dim=head_dim, max_seq_len=128)
        cos, sin = rope(seq_len)

        q_rot, k_rot = apply_rotary_emb(q, k, cos, sin)
        assert q_rot.shape == q.shape
        assert k_rot.shape == k.shape

    def test_apply_rotary_emb_different_positions(self):
        """Test that different positions get different rotations."""
        batch_size, num_heads, seq_len, head_dim = 1, 1, 4, 8
        q = torch.ones(batch_size, num_heads, seq_len, head_dim)
        k = torch.ones(batch_size, num_heads, seq_len, head_dim)

        rope = RotaryEmbedding(dim=head_dim, max_seq_len=8)
        cos, sin = rope(seq_len)

        q_rot, k_rot = apply_rotary_emb(q, k, cos, sin)

        # Different positions should have different rotations
        assert not torch.allclose(q_rot[0, 0, 0], q_rot[0, 0, 1])
        assert not torch.allclose(q_rot[0, 0, 1], q_rot[0, 0, 2])


class TestRoPEAttention:
    """Tests for the RoPEAttention module."""

    def test_init(self):
        """Test RoPEAttention initialization."""
        attn = RoPEAttention(embed_dim=256, num_heads=8)
        assert attn.embed_dim == 256
        assert attn.num_heads == 8
        assert attn.head_dim == 32

    def test_init_invalid_dims(self):
        """Test that invalid dimensions raise error."""
        with pytest.raises(ValueError):
            RoPEAttention(embed_dim=100, num_heads=8)  # 100 not divisible by 8

    def test_forward_shape(self):
        """Test that forward returns correct shape."""
        batch_size, seq_len, embed_dim = 2, 32, 256
        attn = RoPEAttention(embed_dim=embed_dim, num_heads=8)

        x = torch.randn(batch_size, seq_len, embed_dim)
        output = attn(x)
        assert output.shape == (batch_size, seq_len, embed_dim)

    def test_forward_with_mask(self):
        """Test forward with padding mask."""
        batch_size, seq_len, embed_dim = 2, 32, 256
        attn = RoPEAttention(embed_dim=embed_dim, num_heads=8)

        x = torch.randn(batch_size, seq_len, embed_dim)
        mask = torch.zeros(batch_size, seq_len, dtype=torch.bool)
        mask[:, -5:] = True  # Mask last 5 positions

        output = attn(x, mask=mask)
        assert output.shape == (batch_size, seq_len, embed_dim)

    def test_forward_deterministic(self):
        """Test that forward is deterministic in eval mode."""
        attn = RoPEAttention(embed_dim=256, num_heads=8)
        attn.eval()

        x = torch.randn(2, 32, 256)
        out1 = attn(x)
        out2 = attn(x)
        assert torch.allclose(out1, out2)


class TestRoPEIntegration:
    """Integration tests for RoPE with model components."""

    def test_rope_with_dna_encoder(self):
        """Test RoPE integration with DNAEncoder."""
        # Import here to avoid circular imports during test collection
        import sys
        sys.path.insert(0, str(__file__).rsplit("/tests", 1)[0])
        from model import DNAEncoder

        encoder = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=192,
            num_layers=2,
            num_heads=8,
            ff_dim=384,
            dropout=0.1,
            max_seq_len=64,
            use_rope=True,
            rope_base=10000.0,
        )

        batch_size, seq_len = 2, 32
        tokens = torch.randint(0, 64, (batch_size, seq_len))

        output = encoder(tokens)
        assert output.shape == (batch_size, seq_len, 192)

    def test_rope_disabled(self):
        """Test that model works with RoPE disabled."""
        import sys
        sys.path.insert(0, str(__file__).rsplit("/tests", 1)[0])
        from model import DNAEncoder

        encoder = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=192,
            num_layers=2,
            num_heads=8,
            ff_dim=384,
            dropout=0.1,
            max_seq_len=64,
            use_rope=False,
        )

        batch_size, seq_len = 2, 32
        tokens = torch.randint(0, 64, (batch_size, seq_len))

        output = encoder(tokens)
        assert output.shape == (batch_size, seq_len, 192)

    def test_rope_vs_no_rope_different_outputs(self):
        """Test that RoPE produces different outputs than no RoPE."""
        import sys
        sys.path.insert(0, str(__file__).rsplit("/tests", 1)[0])
        from model import DNAEncoder

        torch.manual_seed(42)
        encoder_rope = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=192,
            num_layers=2,
            num_heads=8,
            ff_dim=384,
            dropout=0.0,  # Disable dropout for determinism
            max_seq_len=64,
            use_rope=True,
        )

        torch.manual_seed(42)
        encoder_no_rope = DNAEncoder(
            vocab_size=67,
            kmer=3,
            embedding_dim=192,
            num_layers=2,
            num_heads=8,
            ff_dim=384,
            dropout=0.0,
            max_seq_len=64,
            use_rope=False,
        )

        encoder_rope.eval()
        encoder_no_rope.eval()

        tokens = torch.randint(0, 64, (1, 16))

        with torch.no_grad():
            out_rope = encoder_rope(tokens)
            out_no_rope = encoder_no_rope(tokens)

        # Outputs should be different (RoPE changes attention computation)
        assert not torch.allclose(out_rope, out_no_rope, atol=1e-3)


class TestRoPERelativePosition:
    """Tests verifying RoPE captures relative position information."""

    def test_relative_position_invariance(self):
        """Test that RoPE attention patterns are relative-position aware."""
        rope = RotaryEmbedding(dim=64, max_seq_len=32)

        # Get rotations for positions 0-7
        cos1, sin1 = rope(8)

        # Create Q at position 2 and K at position 5
        # The relative distance is 3
        q_pos2 = torch.randn(1, 1, 1, 64)  # Single query at position 2
        k_pos5 = torch.randn(1, 1, 1, 64)  # Single key at position 5

        # Apply rotations for these specific positions
        cos_q = cos1[:, :, 2:3, :]
        sin_q = sin1[:, :, 2:3, :]
        cos_k = cos1[:, :, 5:6, :]
        sin_k = sin1[:, :, 5:6, :]

        q_rot = (q_pos2 * cos_q) + (rotate_half(q_pos2) * sin_q)
        k_rot = (k_pos5 * cos_k) + (rotate_half(k_pos5) * sin_k)

        # Attention score between rotated Q and K
        score1 = torch.matmul(q_rot, k_rot.transpose(-2, -1))

        # Now shift both positions by same amount (e.g., pos 4 and 7)
        # Relative distance is still 3
        cos_q2 = cos1[:, :, 4:5, :]
        sin_q2 = sin1[:, :, 4:5, :]
        cos_k2 = cos1[:, :, 7:8, :]
        sin_k2 = sin1[:, :, 7:8, :]

        q_rot2 = (q_pos2 * cos_q2) + (rotate_half(q_pos2) * sin_q2)
        k_rot2 = (k_pos5 * cos_k2) + (rotate_half(k_pos5) * sin_k2)

        score2 = torch.matmul(q_rot2, k_rot2.transpose(-2, -1))

        # Scores should be equal because relative position is the same
        assert torch.allclose(score1, score2, atol=1e-5)
