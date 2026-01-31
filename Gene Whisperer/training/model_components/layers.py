"""Core transformer layers for Gene Whisperer models."""
from __future__ import annotations

import math
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from .norm import RMSNorm
from .ffn import SwiGLU, GLUFFN

# Import RoPE from the training directory
import sys
from pathlib import Path
_training_dir = Path(__file__).resolve().parent.parent
if str(_training_dir) not in sys.path:
    sys.path.insert(0, str(_training_dir))
from rope import RotaryEmbedding, apply_rotary_emb


class StochasticDepth(nn.Module):
    """Stochastic Depth (Layer Drop) for regularization.

    Randomly drops entire layers during training, which:
    1. Acts as strong regularization
    2. Reduces effective network depth, speeding up training
    3. Creates an implicit ensemble of different depth networks
    """

    def __init__(self, drop_prob: float = 0.0):
        super().__init__()
        self.drop_prob = drop_prob

    def forward(self, x: torch.Tensor, residual: torch.Tensor) -> torch.Tensor:
        if not self.training or self.drop_prob == 0.0:
            return x + residual

        keep_prob = 1.0 - self.drop_prob
        if torch.rand(1).item() < self.drop_prob:
            return x
        return x + residual / keep_prob


class RelativePositionBias(nn.Module):
    """Relative position bias for attention, similar to T5/ALiBi.

    Adds learnable bias based on relative distance between query and key positions.
    Uses bucketed distances to handle long sequences efficiently.

    Args:
        num_heads: Number of attention heads
        max_distance: Maximum relative distance to consider
        num_buckets: Number of distance buckets (for efficiency)
        bidirectional: Whether to use separate buckets for forward/backward
    """

    def __init__(
        self,
        num_heads: int = 8,
        max_distance: int = 128,
        num_buckets: int = 32,
        bidirectional: bool = True,
    ):
        super().__init__()
        self.num_heads = num_heads
        self.num_buckets = num_buckets
        self.max_distance = max_distance
        self.bidirectional = bidirectional

        self.relative_attention_bias = nn.Embedding(num_buckets, num_heads)
        nn.init.normal_(self.relative_attention_bias.weight, std=0.02)

    def _relative_position_bucket(
        self,
        relative_position: torch.Tensor,
    ) -> torch.Tensor:
        """Convert relative positions to bucket indices using T5-style bucketing."""
        relative_buckets = torch.zeros_like(relative_position)

        if self.bidirectional:
            num_buckets = self.num_buckets
            n = num_buckets // 2
            relative_buckets += (relative_position > 0).long() * n
            relative_position = torch.abs(relative_position)
        else:
            relative_position = torch.clamp(-relative_position, min=0)
            num_buckets = self.num_buckets
            n = num_buckets

        max_exact = n // 2
        is_small = relative_position < max_exact

        relative_position_if_large = max_exact + (
            torch.log(relative_position.float() / max_exact)
            / math.log(self.max_distance / max_exact)
            * (n - max_exact)
        ).long()

        relative_position_if_large = torch.clamp(
            relative_position_if_large, max=n - 1
        )

        relative_buckets += torch.where(
            is_small, relative_position, relative_position_if_large
        )

        return relative_buckets

    def forward(self, seq_len: int, device: torch.device) -> torch.Tensor:
        """Compute relative position bias matrix."""
        context_position = torch.arange(seq_len, device=device)[:, None]
        memory_position = torch.arange(seq_len, device=device)[None, :]
        relative_position = memory_position - context_position
        relative_position_bucket = self._relative_position_bucket(relative_position)
        values = self.relative_attention_bias(relative_position_bucket)
        values = values.permute(2, 0, 1).unsqueeze(0)
        return values


class PreNormTransformerLayer(nn.Module):
    """Pre-norm transformer layer with stochastic depth.

    Supports optional GLU FFN, RoPE, SwiGLU FFN, and RMSNorm.
    """

    def __init__(
        self,
        d_model: int,
        nhead: int,
        dim_feedforward: int,
        dropout: float = 0.1,
        drop_path: float = 0.0,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
        use_rope: bool = False,
        rope_base: float = 10000.0,
        max_seq_len: int = 512,
        ffn_type: str = None,
        norm_type: str = "layernorm",
        ffn_mult: float = None,
    ):
        super().__init__()
        self.d_model = d_model
        self.nhead = nhead
        self.head_dim = d_model // nhead
        self.use_rope = use_rope

        # Normalization layers
        if norm_type == "rmsnorm":
            self.norm1 = RMSNorm(d_model)
            self.norm2 = RMSNorm(d_model)
        else:
            self.norm1 = nn.LayerNorm(d_model)
            self.norm2 = nn.LayerNorm(d_model)

        # Initialize RoPE if enabled
        if use_rope:
            self.rotary = RotaryEmbedding(self.head_dim, max_seq_len, rope_base)
        else:
            self.rotary = None

        # Attention projections
        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)

        self.attn_dropout = nn.Dropout(dropout)
        self.resid_dropout = nn.Dropout(dropout)

        # Determine FFN type
        if ffn_type is None:
            ffn_type = "glu" if use_glu_ffn else "gelu"

        # FFN: standard GELU, GLU-style, or SwiGLU
        if ffn_type == "swiglu":
            swiglu_mult = ffn_mult if ffn_mult is not None else 2.67
            self.ffn = SwiGLU(
                in_features=d_model,
                hidden_features=int(d_model * swiglu_mult),
                out_features=d_model,
                bias=False,
                dropout=dropout,
            )
        elif ffn_type == "glu":
            self.ffn = GLUFFN(
                d_model=d_model,
                ff_dim=dim_feedforward,
                dropout=dropout,
                activation=glu_activation,
            )
        else:
            self.ffn = nn.Sequential(
                nn.Linear(d_model, dim_feedforward),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(dim_feedforward, d_model),
                nn.Dropout(dropout),
            )

        self.drop_path_attn = StochasticDepth(drop_path)
        self.drop_path_ffn = StochasticDepth(drop_path)

        self.scale = self.head_dim ** -0.5
        self._init_weights()

    def _init_weights(self):
        """Scaled initialization for better training dynamics."""
        for module in [self.q_proj, self.k_proj, self.v_proj, self.out_proj]:
            nn.init.xavier_uniform_(module.weight, gain=1.0 / math.sqrt(2))
            if module.bias is not None:
                nn.init.zeros_(module.bias)

    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
        position_bias: Optional[torch.Tensor] = None,
        motif_bias: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        B, L, D = x.shape

        normed = self.norm1(x)

        q = self.q_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)
        k = self.k_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)
        v = self.v_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)

        if self.use_rope and self.rotary is not None:
            cos, sin = self.rotary(L)
            q, k = apply_rotary_emb(q, k, cos, sin)

        attn_weights = (q @ k.transpose(-2, -1)) * self.scale

        if position_bias is not None and not self.use_rope:
            attn_weights = attn_weights + position_bias

        if motif_bias is not None:
            attn_weights = attn_weights + motif_bias

        if key_padding_mask is not None:
            attn_weights = attn_weights.masked_fill(
                key_padding_mask.unsqueeze(1).unsqueeze(2), float('-inf')
            )

        attn_weights = F.softmax(attn_weights, dim=-1)
        attn_weights = self.attn_dropout(attn_weights)

        attn_out = (attn_weights @ v).transpose(1, 2).reshape(B, L, D)
        attn_out = self.out_proj(attn_out)

        x = self.drop_path_attn(x, self.resid_dropout(attn_out))

        ffn_out = self.ffn(self.norm2(x))
        x = self.drop_path_ffn(x, ffn_out)

        return x
