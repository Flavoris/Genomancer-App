"""GeneWhispererV6: Modern LLM Architecture for DNA Classification.

This module implements the V6 architecture which integrates state-of-the-art
techniques from modern language models (LLaMA, Mistral, PaLM):
- Rotary Position Embeddings (RoPE) instead of absolute
- SwiGLU FFN instead of GELU
- RMSNorm instead of LayerNorm
- Optional positional motif bias for biological priors
- Optional Grouped Query Attention (GQA) for parameter reduction

Expected improvement over V5: +2-4% accuracy
"""
from __future__ import annotations

import math
from typing import Optional, Tuple, Union
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

from rope import RotaryEmbedding, apply_rotary_emb
from positional_motif_bias import PositionalMotifBias, DEFAULT_PROMOTER_PRIORS


class RMSNorm(nn.Module):
    """Root Mean Square Layer Normalization.

    Simpler and faster than LayerNorm - no mean centering, no bias term.
    Used in LLaMA, Mistral, and other modern LLMs.
    """

    def __init__(self, dim: int, eps: float = 1e-6):
        super().__init__()
        self.eps = eps
        self.weight = nn.Parameter(torch.ones(dim))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        rms = torch.sqrt(x.pow(2).mean(dim=-1, keepdim=True) + self.eps)
        return (x / rms) * self.weight


class SwiGLU(nn.Module):
    """SwiGLU Feed-Forward Network (LLaMA/PaLM style).

    SwiGLU(x) = (Swish(W1*x) * W3*x) * W2

    Better gradient flow than GELU. Uses smaller expansion ratio (2.67x vs 4x)
    due to gating mechanism.
    """

    def __init__(
        self,
        in_features: int,
        hidden_features: int = None,
        out_features: int = None,
        bias: bool = False,
        dropout: float = 0.0,
    ):
        super().__init__()
        out_features = out_features or in_features
        hidden_features = hidden_features or int(in_features * 8 / 3)
        # Round to multiple of 64 for efficiency
        hidden_features = ((hidden_features + 63) // 64) * 64

        self.w1 = nn.Linear(in_features, hidden_features, bias=bias)
        self.w2 = nn.Linear(hidden_features, out_features, bias=bias)
        self.w3 = nn.Linear(in_features, hidden_features, bias=bias)
        self.dropout = nn.Dropout(dropout) if dropout > 0 else nn.Identity()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.w2(F.silu(self.w1(x)) * self.w3(x)))


class ModernTransformerBlock(nn.Module):
    """Modern Transformer Block with RoPE, SwiGLU, RMSNorm, and optional GQA.

    This block integrates all modern LLM improvements:
    - Pre-norm architecture with RMSNorm
    - Rotary Position Embeddings applied to Q and K
    - SwiGLU feed-forward network
    - Optional Grouped Query Attention for parameter efficiency
    - Optional positional motif bias for biological priors

    Args:
        embed_dim: Model dimension
        num_heads: Number of query heads
        num_kv_heads: Number of key/value heads (for GQA). If None, uses num_heads.
        ff_mult: FFN expansion multiplier (2.67 for SwiGLU)
        dropout: Dropout probability
        max_seq_len: Maximum sequence length for RoPE
        use_rope: Whether to use Rotary Position Embedding
        use_swiglu: Whether to use SwiGLU FFN
        use_rmsnorm: Whether to use RMSNorm
        use_positional_bias: Whether to use positional motif bias
        positional_bias_priors: Dict of (start, end) -> bias_value for motif regions
    """

    def __init__(
        self,
        embed_dim: int = 384,
        num_heads: int = 12,
        num_kv_heads: int = None,
        ff_mult: float = 2.67,
        dropout: float = 0.1,
        max_seq_len: int = 128,
        use_rope: bool = True,
        use_swiglu: bool = True,
        use_rmsnorm: bool = True,
        use_positional_bias: bool = False,
        positional_bias_priors: dict = None,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.num_kv_heads = num_kv_heads or num_heads
        self.head_dim = embed_dim // num_heads
        self.use_rope = use_rope
        self.use_positional_bias = use_positional_bias

        # Validate GQA configuration
        if num_heads % self.num_kv_heads != 0:
            raise ValueError(
                f"num_heads ({num_heads}) must be divisible by "
                f"num_kv_heads ({self.num_kv_heads})"
            )
        self.kv_groups = num_heads // self.num_kv_heads

        # Normalization layers
        NormClass = RMSNorm if use_rmsnorm else nn.LayerNorm
        self.norm1 = NormClass(embed_dim)
        self.norm2 = NormClass(embed_dim)

        # Attention projections
        self.q_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.k_proj = nn.Linear(embed_dim, self.num_kv_heads * self.head_dim, bias=False)
        self.v_proj = nn.Linear(embed_dim, self.num_kv_heads * self.head_dim, bias=False)
        self.out_proj = nn.Linear(embed_dim, embed_dim, bias=False)

        # RoPE for position encoding
        if use_rope:
            self.rotary = RotaryEmbedding(self.head_dim, max_seq_len)
        else:
            self.rotary = None

        # Positional motif bias for biological priors
        if use_positional_bias:
            self.motif_bias = PositionalMotifBias(
                num_heads=num_heads,
                max_seq_len=max_seq_len,
                init_bias_positions=positional_bias_priors or DEFAULT_PROMOTER_PRIORS,
                learnable=True,
            )
        else:
            self.motif_bias = None

        # FFN: SwiGLU or standard GELU
        if use_swiglu:
            self.ffn = SwiGLU(
                in_features=embed_dim,
                hidden_features=int(embed_dim * ff_mult),
                out_features=embed_dim,
                bias=False,
                dropout=dropout,
            )
        else:
            ff_dim = int(embed_dim * 4)
            self.ffn = nn.Sequential(
                nn.Linear(embed_dim, ff_dim),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(ff_dim, embed_dim),
            )

        self.attn_dropout = nn.Dropout(dropout)
        self.resid_dropout = nn.Dropout(dropout)
        self.scale = self.head_dim ** -0.5

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Initialize projections with scaled initialization."""
        for module in [self.q_proj, self.k_proj, self.v_proj, self.out_proj]:
            nn.init.xavier_uniform_(module.weight, gain=1.0 / math.sqrt(2))

    def _repeat_kv(self, x: torch.Tensor, n_rep: int) -> torch.Tensor:
        """Repeat KV heads for Grouped Query Attention."""
        if n_rep == 1:
            return x
        B, num_kv_heads, L, head_dim = x.shape
        x = x[:, :, None, :, :].expand(B, num_kv_heads, n_rep, L, head_dim)
        return x.reshape(B, num_kv_heads * n_rep, L, head_dim)

    def forward(
        self,
        x: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Forward pass.

        Args:
            x: Input tensor (B, L, D)
            mask: Optional padding mask (B, L), True = pad positions

        Returns:
            Output tensor (B, L, D)
        """
        B, L, D = x.shape

        # Pre-norm for attention
        normed = self.norm1(x)

        # Project to Q, K, V
        q = self.q_proj(normed).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(normed).view(B, L, self.num_kv_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(normed).view(B, L, self.num_kv_heads, self.head_dim).transpose(1, 2)

        # Apply RoPE if enabled
        if self.use_rope and self.rotary is not None:
            cos, sin = self.rotary(L)
            q, k = apply_rotary_emb(q, k, cos, sin)

        # Repeat KV for GQA
        k = self._repeat_kv(k, self.kv_groups)
        v = self._repeat_kv(v, self.kv_groups)

        # Scaled dot-product attention
        attn_weights = (q @ k.transpose(-2, -1)) * self.scale

        # Add positional motif bias if enabled
        if self.use_positional_bias and self.motif_bias is not None:
            bias = self.motif_bias.bias[:, :L, :L].unsqueeze(0)
            attn_weights = attn_weights + bias

        # Apply padding mask
        if mask is not None:
            attn_weights = attn_weights.masked_fill(
                mask.unsqueeze(1).unsqueeze(2), float('-inf')
            )

        attn_weights = F.softmax(attn_weights, dim=-1)
        attn_weights = self.attn_dropout(attn_weights)

        # Apply attention to values
        attn_out = (attn_weights @ v).transpose(1, 2).reshape(B, L, D)
        attn_out = self.out_proj(attn_out)

        # Residual connection
        x = x + self.resid_dropout(attn_out)

        # Pre-norm FFN with residual
        x = x + self.ffn(self.norm2(x))

        return x


class GeneWhispererV6(nn.Module):
    """Gene Whisperer v6: Modern LLM Architecture for DNA Classification.

    Improvements over v5:
    - Rotary Position Embeddings (RoPE) instead of absolute
    - SwiGLU FFN instead of GELU
    - RMSNorm instead of LayerNorm
    - Optional positional motif bias
    - Optional Grouped Query Attention (GQA)

    Architecture:
        Input -> Embedding -> [RoPE+SwiGLU Transformer x N] -> Pool -> [+Features] -> Classifier
    """

    def __init__(
        self,
        vocab_size: int,
        embed_dim: int = 384,
        num_layers: int = 12,
        num_heads: int = 12,
        num_kv_heads: int = None,
        ff_mult: float = 2.67,
        dropout: float = 0.1,
        max_seq_len: int = 128,
        pad_token_id: int = None,
        # Engineered features
        engineered_dim: int = 288,
        use_engineered_features: bool = True,
        # Modern options
        use_rope: bool = True,
        use_swiglu: bool = True,
        use_rmsnorm: bool = True,
        use_positional_bias: bool = True,
        positional_bias_priors: dict = None,
    ):
        super().__init__()
        self.pad_token_id = pad_token_id
        self.use_engineered = use_engineered_features and engineered_dim > 0
        self.embed_dim = embed_dim

        # Token embedding (no position embedding with RoPE)
        self.embed = nn.Embedding(vocab_size, embed_dim, padding_idx=pad_token_id)
        self.embed_dropout = nn.Dropout(dropout)

        # Initialize embedding
        nn.init.normal_(self.embed.weight, std=0.02)
        if pad_token_id is not None:
            nn.init.zeros_(self.embed.weight[pad_token_id])

        # Build transformer layers
        self.layers = nn.ModuleList([
            ModernTransformerBlock(
                embed_dim=embed_dim,
                num_heads=num_heads,
                num_kv_heads=num_kv_heads,
                ff_mult=ff_mult,
                dropout=dropout,
                max_seq_len=max_seq_len,
                use_rope=use_rope,
                use_swiglu=use_swiglu,
                use_rmsnorm=use_rmsnorm,
                use_positional_bias=use_positional_bias and (i == 0),  # Only first layer
                positional_bias_priors=positional_bias_priors,
            )
            for i in range(num_layers)
        ])

        # Final normalization
        NormClass = RMSNorm if use_rmsnorm else nn.LayerNorm
        self.final_norm = NormClass(embed_dim)

        # Engineered features processing
        if self.use_engineered:
            self.eng_proj = nn.Sequential(
                nn.Linear(engineered_dim, embed_dim),
                NormClass(embed_dim),
                SwiGLU(embed_dim, int(embed_dim * 2), embed_dim) if use_swiglu else nn.Sequential(
                    nn.Linear(embed_dim, embed_dim * 4),
                    nn.GELU(),
                    nn.Linear(embed_dim * 4, embed_dim),
                ),
                nn.Dropout(dropout),
            )
            classifier_dim = embed_dim * 2
        else:
            classifier_dim = embed_dim

        # Classifier head
        self.classifier = nn.Sequential(
            nn.Linear(classifier_dim, embed_dim),
            NormClass(embed_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(embed_dim, 1),
        )

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: torch.Tensor = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            tokens: (B, L) token indices
            engineered_features: Optional (B, engineered_dim) features
            return_logits: If True, return (probs, logits) tuple

        Returns:
            (B, 1) logits, or (probs, logits) if return_logits=True
        """
        # Embed tokens
        x = self.embed(tokens)
        x = self.embed_dropout(x)

        # Create padding mask
        mask = None
        if self.pad_token_id is not None:
            mask = tokens.eq(self.pad_token_id)

        # Pass through transformer layers
        for layer in self.layers:
            x = layer(x, mask)

        x = self.final_norm(x)

        # Pool (mean, excluding padding)
        if mask is not None:
            x = x.masked_fill(mask.unsqueeze(-1), 0.0)
            lengths = (~mask).sum(dim=1, keepdim=True).float()
            x = x.sum(dim=1) / lengths.clamp(min=1)
        else:
            x = x.mean(dim=1)

        # Add engineered features
        if self.use_engineered and engineered_features is not None:
            eng = self.eng_proj(engineered_features)
            x = torch.cat([x, eng], dim=-1)

        logits = self.classifier(x)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits

    def load_pretrained_weights(
        self,
        checkpoint_path: str,
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ):
        """Load pretrained weights, handling architecture differences.

        Args:
            checkpoint_path: Path to pretrained checkpoint file
            strict: Whether to require exact weight match (default False)
            transfer_mode: One of ["embed_only", "embed_plus_adapter", "none"]
                          Controls which weights to load. For V6, both embed modes
                          load all compatible weights.
        """
        if transfer_mode == "none":
            print("Skipping weight loading (transfer_mode=none)")
            return

        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")

        state_dict = torch.load(checkpoint_path, map_location='cpu')
        if 'model_state_dict' in state_dict:
            state_dict = state_dict['model_state_dict']

        # Filter compatible weights
        model_dict = self.state_dict()
        compatible = {k: v for k, v in state_dict.items()
                     if k in model_dict and model_dict[k].shape == v.shape}

        model_dict.update(compatible)
        self.load_state_dict(model_dict, strict=False)

        print(f"Loaded {len(compatible)}/{len(state_dict)} weights from checkpoint")
        if transfer_mode == "embed_plus_adapter":
            print("V6 architecture: loaded all compatible encoder weights")

    @property
    def encoder(self):
        """Compatibility property - returns self since V6 is the encoder."""
        return self


def create_v6_model(config: dict) -> GeneWhispererV6:
    """Factory function to create GeneWhispererV6 from config dict.

    Args:
        config: Configuration dictionary with model parameters

    Returns:
        Configured GeneWhispererV6 instance
    """
    vocab_size = int(config.get("vocab_size", 4099))
    embed_dim = int(config.get("embedding_dim", 384))
    num_layers = int(config.get("transformer_layers", 12))
    num_heads = int(config.get("transformer_heads", 12))
    dropout = float(config.get("transformer_dropout", 0.1))
    pad_token_id = config.get("pad_token_id")
    if pad_token_id is not None:
        pad_token_id = int(pad_token_id)

    kmer = int(config.get("kmer", 6))
    max_bp_len = int(config.get("max_bp_len", 81))
    max_seq_len = int(config.get("max_seq_len", max_bp_len - kmer + 1))

    engineered_dim = int(config.get("engineered_dim", 288))
    use_engineered = bool(config.get("stage1_use_engineered_features", True))

    # V6 specific options
    use_rope = bool(config.get("use_rope", True))
    use_swiglu = bool(config.get("use_swiglu", True))
    use_rmsnorm = bool(config.get("use_rmsnorm", True))
    use_positional_bias = bool(config.get("use_positional_bias", True))
    use_gqa = bool(config.get("use_gqa", False))
    num_kv_heads = config.get("num_kv_heads")
    if num_kv_heads is not None and use_gqa:
        num_kv_heads = int(num_kv_heads)
    else:
        num_kv_heads = None

    ff_mult = float(config.get("ffn_mult", 2.67))

    # Build positional bias priors from motif_regions config
    positional_bias_priors = None
    motif_regions = config.get("motif_regions")
    if motif_regions and use_positional_bias:
        positional_bias_priors = {}
        for name, (start, end) in motif_regions.items():
            bias_val = 0.15 if "tss" in name.lower() else 0.1
            positional_bias_priors[(start, end)] = bias_val

    return GeneWhispererV6(
        vocab_size=vocab_size,
        embed_dim=embed_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        num_kv_heads=num_kv_heads,
        ff_mult=ff_mult,
        dropout=dropout,
        max_seq_len=max_seq_len,
        pad_token_id=pad_token_id,
        engineered_dim=engineered_dim,
        use_engineered_features=use_engineered,
        use_rope=use_rope,
        use_swiglu=use_swiglu,
        use_rmsnorm=use_rmsnorm,
        use_positional_bias=use_positional_bias,
        positional_bias_priors=positional_bias_priors,
    )
