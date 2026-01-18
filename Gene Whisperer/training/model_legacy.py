"""Legacy CNN/TCN and fusion components for Gene Whisperer.

This module preserves the pre-simplification building blocks so that older
architectures can still be instantiated or inspected.
"""
from __future__ import annotations

from typing import Optional, Tuple

import torch  # pyright: ignore[reportMissingImports]
import torch.nn as nn  # pyright: ignore[reportMissingImports]

from model import PreNormTransformerLayer, RelativePositionBias


class CausalConv1d(nn.Module):
    """
    Causal convolution that ensures output at time t only depends on inputs up to time t.

    Uses left-padding to maintain causality - critical for DNA sequence modeling
    where we want to preserve the natural order of nucleotides.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        dilation: int = 1,
        groups: int = 1,
    ):
        super().__init__()
        self.padding = (kernel_size - 1) * dilation
        self.conv = nn.Conv1d(
            in_channels,
            out_channels,
            kernel_size,
            padding=self.padding,
            dilation=dilation,
            groups=groups,
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, L) input tensor
        Returns:
            (B, C, L) output tensor with same length (causal)
        """
        out = self.conv(x)
        # Remove future padding to maintain causality
        if self.padding > 0:
            out = out[:, :, :-self.padding]
        return out


class TCNBlock(nn.Module):
    """
    Temporal Convolutional Block with dilated causal convolutions.

    Key features from iPro-TCN:
    - Dilated convolutions for exponentially growing receptive field
    - Residual connections for gradient flow
    - Weight normalization for training stability
    - Dropout for regularization
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int = 3,
        dilation: int = 1,
        dropout: float = 0.2,
    ):
        super().__init__()

        # First causal convolution
        self.conv1 = CausalConv1d(in_channels, out_channels, kernel_size, dilation)
        self.bn1 = nn.BatchNorm1d(out_channels)

        # Second causal convolution
        self.conv2 = CausalConv1d(out_channels, out_channels, kernel_size, dilation)
        self.bn2 = nn.BatchNorm1d(out_channels)

        self.relu = nn.GELU()
        self.dropout = nn.Dropout(dropout)

        # Residual connection with 1x1 conv if dimensions differ
        self.downsample = (
            nn.Conv1d(in_channels, out_channels, 1)
            if in_channels != out_channels else None
        )

        # Initialize weights
        self._init_weights()

    def _init_weights(self) -> None:
        """Initialize weights using He initialization."""
        for module in self.modules():
            if isinstance(module, nn.Conv1d):
                nn.init.kaiming_normal_(
                    module.weight, mode="fan_out", nonlinearity="relu"
                )
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, L) input tensor
        Returns:
            (B, C, L) output tensor
        """
        residual = x

        # First conv block
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.dropout(out)

        # Second conv block
        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)
        out = self.dropout(out)

        # Residual connection
        if self.downsample is not None:
            residual = self.downsample(residual)

        return self.relu(out + residual)


class TCNEncoder(nn.Module):
    """
    Multi-level TCN encoder with exponentially growing receptive field.

    Architecture inspired by iPro-TCN:
    - Multiple TCN blocks with increasing dilation rates (1, 2, 4, 8, ...)
    - Each level doubles the receptive field
    - Total receptive field = kernel_size * (2^num_levels - 1)

    For DNA sequences of 81bp, 4 levels with kernel_size=3 gives:
    receptive field = 3 * (2^4 - 1) = 45bp coverage per position
    """

    def __init__(
        self,
        in_channels: int,
        hidden_channels: int = 256,
        num_levels: int = 4,
        kernel_size: int = 3,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.num_levels = num_levels

        layers = []
        for i in range(num_levels):
            dilation = 2 ** i  # Exponential dilation: 1, 2, 4, 8, ...
            in_ch = in_channels if i == 0 else hidden_channels
            layers.append(
                TCNBlock(
                    in_channels=in_ch,
                    out_channels=hidden_channels,
                    kernel_size=kernel_size,
                    dilation=dilation,
                    dropout=dropout,
                )
            )

        self.network = nn.ModuleList(layers)
        self.out_channels = hidden_channels

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, L) input tensor
        Returns:
            (B, hidden_channels, L) output tensor
        """
        for layer in self.network:
            x = layer(x)
        return x


class MultiScaleCNN(nn.Module):
    """
    Multi-scale parallel CNN for capturing DNA motifs of different lengths.

    Promoter motifs vary in size:
    - TATA box: ~7bp
    - -10 element: ~6bp
    - -35 element: ~6bp
    - Initiator: ~7bp
    - Downstream elements: ~15-20bp

    This module uses parallel convolutions with different kernel sizes
    to capture all these patterns simultaneously.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels_per_scale: int = 64,
        kernel_sizes: Tuple[int, ...] = (3, 5, 7, 9, 15),
        dropout: float = 0.1,
    ):
        super().__init__()
        self.kernel_sizes = kernel_sizes
        self.num_scales = len(kernel_sizes)

        # Parallel convolution branches for different motif lengths
        self.convs = nn.ModuleList([
            nn.Conv1d(
                in_channels,
                out_channels_per_scale,
                kernel_size=ks,
                padding=ks // 2,  # Same padding to maintain length
            )
            for ks in kernel_sizes
        ])

        total_channels = out_channels_per_scale * self.num_scales
        self.bn = nn.BatchNorm1d(total_channels)
        self.activation = nn.GELU()
        self.dropout = nn.Dropout(dropout)

        # Channel reduction to manage complexity
        self.channel_reduce = nn.Conv1d(total_channels, total_channels // 2, 1)
        self.bn_reduce = nn.BatchNorm1d(total_channels // 2)

        self.out_channels = total_channels // 2

        self._init_weights()

    def _init_weights(self) -> None:
        """Initialize convolution weights."""
        for module in self.modules():
            if isinstance(module, nn.Conv1d):
                nn.init.kaiming_normal_(
                    module.weight, mode="fan_out", nonlinearity="relu"
                )
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, L) input tensor
        Returns:
            (B, out_channels, L) multi-scale features
        """
        # Apply parallel convolutions
        conv_outputs = [conv(x) for conv in self.convs]

        # Concatenate along channel dimension
        out = torch.cat(conv_outputs, dim=1)  # (B, num_scales * out_per_scale, L)

        out = self.bn(out)
        out = self.activation(out)
        out = self.dropout(out)

        # Reduce channels
        out = self.channel_reduce(out)
        out = self.bn_reduce(out)
        out = self.activation(out)

        return out


class MultiScaleTCN(nn.Module):
    """
    Combined Multi-Scale CNN + TCN architecture.

    This hybrid approach combines:
    1. Multi-Scale CNN: Captures motifs of different lengths in parallel
    2. TCN: Captures long-range dependencies with dilated convolutions

    Architecture:
    Input → Multi-Scale CNN → TCN → Output

    This is more powerful than either alone because:
    - Multi-Scale CNN extracts diverse local patterns
    - TCN builds hierarchical long-range representations
    """

    def __init__(
        self,
        in_channels: int,
        multiscale_out_per_branch: int = 64,
        multiscale_kernels: Tuple[int, ...] = (3, 5, 7, 9, 15),
        tcn_hidden: int = 256,
        tcn_levels: int = 4,
        tcn_kernel: int = 3,
        dropout: float = 0.2,
    ):
        super().__init__()

        # Multi-scale CNN for local pattern extraction
        self.multiscale = MultiScaleCNN(
            in_channels=in_channels,
            out_channels_per_scale=multiscale_out_per_branch,
            kernel_sizes=multiscale_kernels,
            dropout=dropout,
        )

        # TCN for long-range dependencies
        self.tcn = TCNEncoder(
            in_channels=self.multiscale.out_channels,
            hidden_channels=tcn_hidden,
            num_levels=tcn_levels,
            kernel_size=tcn_kernel,
            dropout=dropout,
        )

        self.out_channels = tcn_hidden

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, L) input tensor
        Returns:
            (B, tcn_hidden, L) hierarchical features
        """
        # Extract multi-scale local patterns
        x = self.multiscale(x)

        # Build long-range dependencies
        x = self.tcn(x)

        return x


class GatedAttentionFusion(nn.Module):
    """
    Attention-based fusion of sequence and engineered features with residual.

    Uses cross-attention where engineered features attend to sequence features,
    plus a gating mechanism to control the contribution of each modality.
    This is more expressive than simple concatenation.

    Key improvement: Residual connection from sequence features ensures
    the model doesn't lose important sequence information during fusion.
    """

    def __init__(
        self,
        seq_dim: int,
        eng_dim: int,
        hidden_dim: int = 256,
        num_heads: int = 4,
        dropout: float = 0.15,  # Reduced from 0.2
    ):
        super().__init__()
        self.seq_dim = seq_dim
        self.eng_dim = eng_dim
        self.hidden_dim = hidden_dim

        # Project both to same dimension for attention
        self.seq_proj = nn.Linear(seq_dim, hidden_dim)
        self.eng_proj = nn.Linear(eng_dim, hidden_dim)

        # Cross-attention: engineered features query sequence features
        self.cross_attn = nn.MultiheadAttention(
            embed_dim=hidden_dim,
            num_heads=num_heads,
            dropout=dropout,
            batch_first=True,
        )

        # Self-attention for sequence features to enhance them
        self.seq_self_attn = nn.MultiheadAttention(
            embed_dim=hidden_dim,
            num_heads=num_heads,
            dropout=dropout,
            batch_first=True,
        )

        # Gating mechanism to balance sequence vs engineered contributions
        self.gate_seq = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Sigmoid(),
        )
        self.gate_eng = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Sigmoid(),
        )

        # Final projection after fusion with residual pathway
        self.out_proj = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )

        # Residual projection for sequence features (preserves seq info)
        self.seq_residual_proj = nn.Linear(seq_dim, hidden_dim)

        self.output_dim = hidden_dim

    def forward(
        self,
        seq_features: torch.Tensor,
        eng_features: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            seq_features: (B, seq_dim) pooled sequence features
            eng_features: (B, eng_dim) processed engineered features
        Returns:
            (B, hidden_dim) fused representation
        """
        # Residual from original sequence features
        seq_residual = self.seq_residual_proj(seq_features)  # (B, hidden_dim)

        # Project to common dimension
        seq_h = self.seq_proj(seq_features)  # (B, hidden_dim)
        eng_h = self.eng_proj(eng_features)  # (B, hidden_dim)

        # Cross-attention: engineered features attend to sequence
        # Reshape for attention: (B, 1, hidden_dim) as single token
        seq_h_expanded = seq_h.unsqueeze(1)  # (B, 1, hidden_dim)
        eng_h_expanded = eng_h.unsqueeze(1)  # (B, 1, hidden_dim)

        # Engineered queries sequence (learn what sequence info is relevant for eng features)
        attn_out, _ = self.cross_attn(
            query=eng_h_expanded,
            key=seq_h_expanded,
            value=seq_h_expanded,
        )
        attn_out = attn_out.squeeze(1)  # (B, hidden_dim)

        # Concatenate for gating
        concat = torch.cat([seq_h, eng_h], dim=1)  # (B, hidden_dim * 2)

        # Compute gates
        g_seq = self.gate_seq(concat)  # (B, hidden_dim)
        g_eng = self.gate_eng(concat)  # (B, hidden_dim)

        # Apply gates with sequence residual
        gated_seq = g_seq * (seq_h + seq_residual)  # Add residual for seq preservation
        gated_eng = g_eng * (eng_h + attn_out)  # Add attention-enhanced eng features

        # Fuse
        fused = torch.cat([gated_seq, gated_eng], dim=1)  # (B, hidden_dim * 2)
        return self.out_proj(fused)


class PostCNNTransformerAdapter(nn.Module):
    """
    Adapter to apply pretrained transformer layers after CNN/TCN.

    This bridges the dimension mismatch between CNN/TCN output (tcn_hidden)
    and the pretrained transformer (embedding_dim), allowing us to reuse
    the pretrained transformer weights from MLM.

    Architecture:
    CNN/TCN output (tcn_hidden) → Projection → Pretrained Transformer → Output
    """

    def __init__(
        self,
        input_dim: int,
        transformer_dim: int,
        num_layers: int = 3,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        drop_path_rate: float = 0.1,  # Stochastic depth
        use_glu_ffn: bool = False,  # Use GLU-style FFN
        glu_activation: str = "gelu",  # Activation for GLU: "gelu" or "silu"
        use_relative_position_bias: bool = False,  # Relative position bias
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.transformer_dim = transformer_dim
        self.use_relative_position_bias = use_relative_position_bias

        # Projection to transformer dimension
        self.input_proj = nn.Linear(input_dim, transformer_dim)
        self.input_norm = nn.LayerNorm(transformer_dim)

        # Relative position bias (shared across all layers)
        if use_relative_position_bias:
            self.relative_position_bias = RelativePositionBias(
                num_heads=num_heads,
                max_distance=relative_position_max_distance,
                num_buckets=relative_position_num_buckets,
                bidirectional=True,
            )
        else:
            self.relative_position_bias = None

        # Transformer layers (can be initialized from pretrained DNAEncoder)
        if num_layers < 0:
            raise ValueError(f"num_layers must be >= 0, got {num_layers}")
        dpr = (
            [x.item() for x in torch.linspace(0, drop_path_rate, num_layers)]
            if num_layers > 0
            else []
        )
        self.layers = nn.ModuleList([
            PreNormTransformerLayer(
                d_model=transformer_dim,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
                drop_path=dpr[i],
                use_glu_ffn=use_glu_ffn,
                glu_activation=glu_activation,
            )
            for i in range(num_layers)
        ])

        self.final_norm = nn.LayerNorm(transformer_dim)

        # Project back to original dimension for pooling
        self.output_proj = nn.Linear(transformer_dim, input_dim)
        self.output_norm = nn.LayerNorm(input_dim)

    def load_pretrained_layers(self, encoder: "DNAEncoder", num_layers: int) -> int:
        """
        Load transformer layers from a pretrained DNAEncoder.

        Args:
            encoder: Pretrained DNAEncoder
            num_layers: Number of layers to copy (from the top of encoder)
        Returns:
            Number of layers successfully loaded
        """
        loaded = 0
        encoder_layers = list(encoder.layers)

        # Copy from the TOP of the pretrained encoder (most general layers)
        # to our adapter layers
        for i in range(min(num_layers, len(self.layers), len(encoder_layers))):
            src_idx = len(encoder_layers) - num_layers + i  # Take from top
            if src_idx >= 0:
                self.layers[i].load_state_dict(encoder_layers[src_idx].state_dict())
                loaded += 1

        return loaded

    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, input_dim) features from CNN/TCN
            key_padding_mask: Optional (B, L) padding mask
        Returns:
            (B, L, input_dim) contextualized features
        """
        if len(self.layers) == 0:
            return x

        B, L, _ = x.shape

        # Project to transformer dimension
        x = self.input_proj(x)
        x = self.input_norm(x)

        # Compute relative position bias once (shared across all layers)
        position_bias = None
        if self.relative_position_bias is not None:
            position_bias = self.relative_position_bias(L, x.device)

        # Apply transformer layers
        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask, position_bias=position_bias)

        x = self.final_norm(x)

        # Project back to input dimension
        x = self.output_proj(x)
        x = self.output_norm(x)

        return x


__all__ = [
    "CausalConv1d",
    "TCNBlock",
    "TCNEncoder",
    "MultiScaleCNN",
    "MultiScaleTCN",
    "GatedAttentionFusion",
    "PostCNNTransformerAdapter",
]
