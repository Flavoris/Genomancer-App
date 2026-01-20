"""Rotary Position Embedding (RoPE) implementation for Gene Whisperer.

Rotary Position Embedding from RoFormer/LLaMA. Instead of adding position
embeddings, RoPE rotates query and key vectors based on their positions.
This captures relative positions better than absolute embeddings, which is
important for DNA motif detection.

Reference: https://arxiv.org/abs/2104.09864 (RoFormer)
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class RotaryEmbedding(nn.Module):
    """Rotary Position Embedding from RoFormer/LLaMA.

    Instead of adding position embeddings, RoPE rotates query and key vectors
    based on their positions. This captures relative positions better than
    absolute embeddings, which is important for DNA motif detection.

    Args:
        dim: Dimension of the embedding (typically head_dim)
        max_seq_len: Maximum sequence length to cache
        base: Base for the frequency computation (default: 10000.0)
    """

    def __init__(self, dim: int, max_seq_len: int = 512, base: float = 10000.0):
        super().__init__()
        self.dim = dim
        self.max_seq_len = max_seq_len
        self.base = base

        # Compute inverse frequencies
        inv_freq = 1.0 / (base ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer("inv_freq", inv_freq, persistent=False)

        # Build cos/sin cache
        self._build_cache(max_seq_len)

    def _build_cache(self, seq_len: int) -> None:
        """Precompute cos and sin for all positions."""
        self.max_seq_len = seq_len
        t = torch.arange(seq_len, device=self.inv_freq.device, dtype=self.inv_freq.dtype)
        freqs = torch.outer(t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)

        self.register_buffer("cos_cached", emb.cos(), persistent=False)
        self.register_buffer("sin_cached", emb.sin(), persistent=False)

    def forward(self, seq_len: int) -> tuple[torch.Tensor, torch.Tensor]:
        """Get cos and sin for sequence length.

        Args:
            seq_len: Current sequence length

        Returns:
            Tuple of (cos, sin) tensors with shape (1, 1, L, D)
        """
        if seq_len > self.max_seq_len:
            self._build_cache(seq_len)

        return (
            self.cos_cached[:seq_len].unsqueeze(0).unsqueeze(0),  # (1, 1, L, D)
            self.sin_cached[:seq_len].unsqueeze(0).unsqueeze(0),
        )


def rotate_half(x: torch.Tensor) -> torch.Tensor:
    """Rotate half the hidden dims.

    Splits the input tensor in half along the last dimension and
    swaps the halves with negation on the first half.

    Args:
        x: Input tensor of shape (..., D)

    Returns:
        Rotated tensor of shape (..., D)
    """
    x1, x2 = x.chunk(2, dim=-1)
    return torch.cat((-x2, x1), dim=-1)


def apply_rotary_emb(
    q: torch.Tensor,
    k: torch.Tensor,
    cos: torch.Tensor,
    sin: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Apply rotary embedding to query and key tensors.

    Args:
        q: Query tensor of shape (B, H, L, D)
        k: Key tensor of shape (B, H, L, D)
        cos: Cosine tensor of shape (1, 1, L, D)
        sin: Sine tensor of shape (1, 1, L, D)

    Returns:
        Tuple of (rotated_q, rotated_k) with same shapes as inputs
    """
    q_embed = (q * cos) + (rotate_half(q) * sin)
    k_embed = (k * cos) + (rotate_half(k) * sin)
    return q_embed, k_embed


class RoPEAttention(nn.Module):
    """Multi-head attention with Rotary Position Embedding.

    This module implements multi-head self-attention using RoPE for
    position encoding. RoPE is applied to query and key projections
    before computing attention scores.

    Args:
        embed_dim: Total embedding dimension
        num_heads: Number of attention heads
        dropout: Dropout probability
        max_seq_len: Maximum sequence length for RoPE cache
        rope_base: Base for RoPE frequency computation
    """

    def __init__(
        self,
        embed_dim: int,
        num_heads: int,
        dropout: float = 0.1,
        max_seq_len: int = 512,
        rope_base: float = 10000.0,
    ):
        super().__init__()
        if embed_dim % num_heads != 0:
            raise ValueError(
                f"embed_dim ({embed_dim}) must be divisible by num_heads ({num_heads})"
            )

        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads

        self.q_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.k_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.v_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.out_proj = nn.Linear(embed_dim, embed_dim, bias=False)

        self.rotary = RotaryEmbedding(self.head_dim, max_seq_len, rope_base)
        self.dropout = nn.Dropout(dropout)
        self.scale = self.head_dim ** -0.5

    def forward(
        self,
        x: torch.Tensor,
        mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward pass with RoPE attention.

        Args:
            x: Input tensor of shape (B, L, D)
            mask: Optional padding mask of shape (B, L), True = pad

        Returns:
            Output tensor of shape (B, L, D)
        """
        B, L, D = x.shape

        # Project to Q, K, V
        q = self.q_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)

        # Apply rotary embeddings
        cos, sin = self.rotary(L)
        q, k = apply_rotary_emb(q, k, cos, sin)

        # Attention scores
        attn = torch.matmul(q, k.transpose(-2, -1)) * self.scale

        if mask is not None:
            attn = attn.masked_fill(mask.unsqueeze(1).unsqueeze(2), float("-inf"))

        attn = F.softmax(attn, dim=-1)
        attn = self.dropout(attn)

        out = torch.matmul(attn, v)
        out = out.transpose(1, 2).contiguous().view(B, L, D)

        return self.out_proj(out)
