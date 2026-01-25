"""Legacy CNN/TCN components for Gene Whisperer.

These classes are kept for backward compatibility with older checkpoints. The
simplified transformer-only Stage 1 architecture performs better in ablation
studies (84.1% accuracy vs 83.3%) while using 38% fewer parameters.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Tuple, Union

import torch  # pyright: ignore[reportMissingImports]
import torch.nn as nn  # pyright: ignore[reportMissingImports]

LOGGER = logging.getLogger(__name__)

_LEGACY_WARNING = (
    "GeneWhispererStage1Legacy is deprecated. Use GeneWhispererStage1 (simplified) instead. "
    "The simplified model achieves 84.1% accuracy vs 83.3% with 38% fewer parameters."
)


def _warn_legacy() -> None:
    warnings.warn(_LEGACY_WARNING, DeprecationWarning)


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
        _warn_legacy()
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
        _warn_legacy()

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
            if in_channels != out_channels
            else None
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
        _warn_legacy()
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
        _warn_legacy()
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
        _warn_legacy()

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
        _warn_legacy()
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
        _warn_legacy()
        from model import PreNormTransformerLayer, RelativePositionBias

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

        _, seq_len, _ = x.shape

        # Project to transformer dimension
        x = self.input_proj(x)
        x = self.input_norm(x)

        # Compute relative position bias once (shared across all layers)
        position_bias = None
        if self.relative_position_bias is not None:
            position_bias = self.relative_position_bias(seq_len, x.device)

        # Apply transformer layers
        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask, position_bias=position_bias)

        x = self.final_norm(x)

        # Project back to input dimension
        x = self.output_proj(x)
        x = self.output_norm(x)

        return x


class GeneWhispererStage1Legacy(nn.Module):
    """
    Legacy Stage 1 classifier with CNN/TCN backbone and attention fusion.

    Architecture: Embedding → CNN/TCN → Transformer → Pooling → Attention Fusion → Classifier

    Deprecated: use GeneWhispererStage1 (simplified transformer-only backbone).
    """

    def __init__(
        self,
        vocab_size: int = 67,
        kmer: int = 3,
        embedding_dim: int = 192,
        num_layers: int = 6,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        pad_token_id: Optional[int] = None,
        encoder: Optional["DNAEncoder"] = None,
        engineered_dim: int = 288,
        use_engineered_features: bool = True,
        use_attention_pool: bool = True,
        # TCN parameters
        use_tcn: bool = True,
        tcn_hidden: int = 256,
        tcn_levels: int = 4,
        tcn_kernel: int = 3,
        # Multi-scale CNN parameters
        multiscale_channels: int = 64,
        multiscale_kernels: Tuple[int, ...] = (3, 5, 7, 9, 15),
        # LSTM parameters (kept for compatibility but optional)
        lstm_hidden: int = 192,
        # New architecture parameters
        post_cnn_transformer_layers: int = 3,
        engineered_mlp_hidden: int = 512,
        engineered_mlp_output: int = 256,
        engineered_mlp_pre_norm: bool = True,
        engineered_mlp_input_dropout: float = 0.05,
        engineered_mlp_use_gated: bool = True,
        engineered_mlp_use_residual: bool = True,
        fusion_hidden: int = 384,
        # Stochastic depth (drop path) rate for encoder
        drop_path_rate: float = 0.1,
        # Maximum sequence length for positional embeddings
        max_seq_len: int = 256,
        # Relative position bias parameters
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        # GLU FFN parameters
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
    ):
        super().__init__()
        _warn_legacy()
        from model import AttentionPooling, DNAEncoder, EngineeredFeatureMLP

        self.use_tcn = use_tcn
        self.pad_token_id = pad_token_id
        self.max_seq_len = int(max_seq_len)
        self.pos_embedding = nn.Embedding(self.max_seq_len, embedding_dim)

        # Step 1: Embedding layer (from DNAEncoder, but we only use embedding)
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_token_id)
        self.embed_dropout = nn.Dropout(dropout)

        # Store the full encoder for MLM weight loading compatibility
        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )
        # Copy embedding weights from encoder
        self.embedding.weight = self._full_encoder.embedding.weight

        # Step 2: CNN/TCN for local pattern extraction (FIRST)
        if use_tcn:
            self.feature_extractor = MultiScaleTCN(
                in_channels=embedding_dim,
                multiscale_out_per_branch=multiscale_channels,
                multiscale_kernels=multiscale_kernels,
                tcn_hidden=tcn_hidden,
                tcn_levels=tcn_levels,
                tcn_kernel=tcn_kernel,
                dropout=dropout,
            )
            cnn_out_dim = tcn_hidden
        else:
            self.conv = nn.Conv1d(embedding_dim, 256, kernel_size=7, padding=3)
            self.conv_bn = nn.BatchNorm1d(256)
            self.conv_act = nn.GELU()
            cnn_out_dim = 256

        # Step 3: Transformer adapter for global contextualization (AFTER CNN)
        # Uses same dimension as pretrained encoder for weight transfer
        self.post_cnn_transformer = PostCNNTransformerAdapter(
            input_dim=cnn_out_dim,
            transformer_dim=embedding_dim,  # Match pretrained encoder dimension
            num_layers=post_cnn_transformer_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
        )

        # Step 4: Attention pooling
        self.use_attention_pool = use_attention_pool
        if use_attention_pool:
            self.pool = AttentionPooling(cnn_out_dim, num_heads=4)

        seq_out_dim = cnn_out_dim

        # Step 5: Engineered features MLP branch
        self.use_engineered_features = use_engineered_features and engineered_dim > 0
        self.engineered_dim = engineered_dim

        if self.use_engineered_features:
            self.engineered_mlp = EngineeredFeatureMLP(
                input_dim=engineered_dim,
                hidden_dim=engineered_mlp_hidden,
                output_dim=engineered_mlp_output,
                dropout=dropout,
                use_pre_norm=engineered_mlp_pre_norm,
                input_dropout=engineered_mlp_input_dropout,
                use_gated=engineered_mlp_use_gated,
                use_residual=engineered_mlp_use_residual,
            )
            eng_out_dim = self.engineered_mlp.output_dim

            # Step 6: Attention-based fusion
            self.fusion = GatedAttentionFusion(
                seq_dim=seq_out_dim,
                eng_dim=eng_out_dim,
                hidden_dim=fusion_hidden,
                num_heads=4,
                dropout=dropout,
            )
            classifier_in = self.fusion.output_dim
        else:
            classifier_in = seq_out_dim

        # Step 7: Classifier head with reduced dropout (outputs logits)
        # Original had 0.4/0.3 dropout which was too aggressive
        classifier_hidden = 512 if use_tcn else 384
        classifier_dropout = dropout * 1.5  # Scale with base dropout (0.15 -> 0.225)

        self.classifier = nn.Sequential(
            nn.Linear(classifier_in, classifier_hidden),
            nn.LayerNorm(classifier_hidden),
            nn.GELU(),
            nn.Dropout(classifier_dropout),  # Was 0.4, now ~0.22
            nn.Linear(classifier_hidden, classifier_hidden // 2),
            nn.LayerNorm(classifier_hidden // 2),
            nn.GELU(),
            nn.Dropout(classifier_dropout * 0.8),  # Was 0.3, now ~0.18
            nn.Linear(classifier_hidden // 2, 1),
        )

    @property
    def encoder(self) -> "DNAEncoder":
        """Compatibility property for MLM weight loading."""
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """
        Load pretrained MLM weights into the model.

        This properly transfers:
        1. Embedding weights (directly used in forward pass)
        2. Transformer layer weights (into post_cnn_transformer adapter) when requested

        Args:
            checkpoint_path: Path to pretrained encoder checkpoint
            strict: Whether to require exact match
            transfer_mode: One of ["embed_only", "embed_plus_adapter", "none"]
        """
        if transfer_mode == "none":
            LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
            return
        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")
        # Load into the full encoder first
        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

        # Transfer transformer layers to post_cnn_transformer
        if transfer_mode == "embed_plus_adapter":
            num_layers = len(self.post_cnn_transformer.layers)
            loaded = self.post_cnn_transformer.load_pretrained_layers(self._full_encoder, num_layers)
            LOGGER.info("Transferred %d transformer layers to post-CNN adapter", loaded)

    def load_pretrained(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
    ) -> None:
        """
        Load checkpoint with backward compatibility for architecture enhancements.

        New layers (relative position bias, GLU FFN components) will be randomly
        initialized if not present in the checkpoint. This allows loading older
        checkpoints into models with new architectural features enabled.

        Args:
            checkpoint_path: Path to checkpoint file
            strict: If False (default), allows missing keys for new components
        """
        from model import _load_checkpoint_file

        state_dict = _load_checkpoint_file(Path(checkpoint_path))

        # Handle nested state dict (e.g., from training checkpoints)
        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]

        # Get current model state
        model_state = self.state_dict()

        # Filter to only keys that exist in both and have matching shapes
        filtered_state = {}
        for key, value in state_dict.items():
            if key in model_state:
                if model_state[key].shape == value.shape:
                    filtered_state[key] = value
                else:
                    LOGGER.warning(
                        "Shape mismatch for %s: checkpoint %s vs model %s (skipping)",
                        key,
                        value.shape,
                        model_state[key].shape,
                    )

        # Log which keys are new (in model but not in checkpoint)
        new_keys = set(model_state.keys()) - set(filtered_state.keys())
        if new_keys:
            # Group by component for cleaner logging
            rel_pos_keys = [k for k in new_keys if "relative_position" in k]
            glu_keys = [k for k in new_keys if "gate_proj" in k or "up_proj" in k or "down_proj" in k]
            other_keys = [k for k in new_keys if k not in rel_pos_keys and k not in glu_keys]

            if rel_pos_keys:
                LOGGER.info(
                    "New relative position bias parameters (randomly initialized): %d keys",
                    len(rel_pos_keys),
                )
            if glu_keys:
                LOGGER.info("New GLU FFN parameters (randomly initialized): %d keys", len(glu_keys))
            if other_keys:
                LOGGER.info("Other new parameters (randomly initialized): %s", other_keys)

        # Log which keys are in checkpoint but not in model (removed components)
        missing_keys = set(state_dict.keys()) - set(model_state.keys())
        if missing_keys:
            LOGGER.info("Checkpoint keys not in model (ignored): %d keys", len(missing_keys))

        # Load the filtered state dict
        # When strict=True, fail if model has parameters not covered by checkpoint
        if strict and new_keys:
            raise RuntimeError(
                f"strict=True but {len(new_keys)} model parameters not in checkpoint: "
                f"{list(new_keys)[:5]}{'...' if len(new_keys) > 5 else ''}"
            )
        self.load_state_dict(filtered_state, strict=False)
        LOGGER.info(
            "Loaded %d/%d parameters from checkpoint (new: %d, ignored: %d)",
            len(filtered_state),
            len(state_dict),
            len(new_keys),
            len(missing_keys),
        )

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        """Generate padding mask from tokens."""
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def _sequence_backbone_from_embeds(
        self,
        embeds: torch.Tensor,
        padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Process pre-computed embeddings through CNN/TCN → Transformer → Pool.

        This method enables embedding-space mixup and distillation by allowing
        forward passes that start from embeddings rather than token indices.

        Args:
            embeds: (B, L, embedding_dim) pre-computed embeddings
            padding_mask: Optional (B, L) padding mask (True = pad)
        Returns:
            (B, D) pooled sequence features
        """
        # Add positional embeddings
        batch_size, seq_len, _ = embeds.shape
        if seq_len > self.pos_embedding.num_embeddings:
            raise ValueError(
                f"Sequence token length {seq_len} exceeds max_seq_len={self.pos_embedding.num_embeddings}. "
                f"Increase max_bp_len or pass a larger max_seq_len."
            )
        pos_ids = torch.arange(seq_len, device=embeds.device).unsqueeze(0).expand(batch_size, seq_len)
        x = embeds + self.pos_embedding(pos_ids)
        x = self.embed_dropout(x)

        # CNN/TCN (transpose for conv: B, L, D -> B, D, L)
        x = x.transpose(1, 2)

        if self.use_tcn:
            x = self.feature_extractor(x)  # (B, tcn_hidden, L)
        else:
            x = self.conv(x)
            x = self.conv_bn(x)
            x = self.conv_act(x)

        # Transpose back: (B, D, L) -> (B, L, D)
        x = x.transpose(1, 2)

        # Post-CNN Transformer for global context
        x = self.post_cnn_transformer(x, key_padding_mask=padding_mask)

        # Pooling
        if self.use_attention_pool:
            pooled = self.pool(x, padding_mask)
        else:
            pooled = torch.max(x, dim=1).values

        return pooled

    def _sequence_backbone(
        self,
        tokens: torch.Tensor,
    ) -> torch.Tensor:
        """
        Process sequence through embedding → CNN/TCN → Transformer → Pool.

        Args:
            tokens: (B, L) token indices
        Returns:
            (B, D) pooled sequence features
        """
        embeds = self.embedding(tokens)
        padding_mask = self._get_padding_mask(tokens)
        return self._sequence_backbone_from_embeds(embeds, padding_mask)

    def extract_features(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Extract sequence and fused features."""
        seq_features = self._sequence_backbone(tokens)

        if engineered_features is not None:
            assert engineered_features.size(-1) == self.engineered_dim, (
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )
        if self.use_engineered_features and engineered_features is not None:
            eng_features = self.engineered_mlp(engineered_features)
            fused = self.fusion(seq_features, eng_features)
            return seq_features, fused

        return seq_features, seq_features

    def extract_features_from_embeds(
        self,
        embeds: torch.Tensor,
        padding_mask: Optional[torch.Tensor] = None,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Extract sequence and fused features from pre-computed embeddings.

        This method enables embedding-space mixup and distillation by allowing
        feature extraction that starts from embeddings rather than token indices.

        Args:
            embeds: (B, L, embedding_dim) pre-computed embeddings
            padding_mask: Optional (B, L) padding mask (True = pad)
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        Returns:
            Tuple of (seq_features, fused_features) both (B, D)
        """
        seq_features = self._sequence_backbone_from_embeds(embeds, padding_mask)

        if engineered_features is not None:
            assert engineered_features.size(-1) == self.engineered_dim, (
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )
        if self.use_engineered_features and engineered_features is not None:
            eng_features = self.engineered_mlp(engineered_features)
            fused = self.fusion(seq_features, eng_features)
            return seq_features, fused

        return seq_features, seq_features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """
        Forward pass with legacy architecture order.

        Args:
            tokens: (B, L) k-mer token indices
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        return_logits: If True, also return sigmoid probabilities alongside logits
        Returns:
            (B, 1) promoter logits, or tuple of (probs, logits) if return_logits=True
        """
        _, fused = self.extract_features(tokens, engineered_features)

        logits = self.classifier(fused)
        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits

    def forward_from_embeds(
        self,
        embeds: torch.Tensor,
        padding_mask: Optional[torch.Tensor] = None,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """
        Forward pass from pre-computed embeddings (parity path for mixup/distillation).

        This method enables embedding-space mixup and distillation by allowing
        forward passes that start from embeddings rather than token indices.
        Behavior is identical to forward() except it takes embeddings directly.

        Args:
            embeds: (B, L, embedding_dim) pre-computed embeddings
            padding_mask: Optional (B, L) padding mask (True = pad)
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        return_logits: If True, also return sigmoid probabilities alongside logits
        Returns:
            (B, 1) promoter logits, or tuple of (probs, logits) if return_logits=True
        """
        _, fused = self.extract_features_from_embeds(embeds, padding_mask, engineered_features)

        logits = self.classifier(fused)
        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits


__all__ = [
    "CausalConv1d",
    "TCNBlock",
    "TCNEncoder",
    "MultiScaleCNN",
    "MultiScaleTCN",
    "GatedAttentionFusion",
    "PostCNNTransformerAdapter",
    "GeneWhispererStage1Legacy",
]
