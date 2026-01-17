"""Enhanced model definitions for Gene Whisperer with TCN and multi-scale CNN.

Key improvements from DNABERT, msBERT, Procables, and iPro-TCN papers:
1. Multi-scale k-mer encoding (k=3,4,5,6) for capturing patterns at different scales
2. Temporal Convolutional Network (TCN) with dilated causal convolutions
3. Multi-scale parallel CNN for capturing different motif lengths (3-15bp)
4. Attention pooling instead of max pooling for position-aware aggregation
5. CoreML-compatible design for iPhone deployment
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import List, Optional, Tuple, Union

import torch  # pyright: ignore[reportMissingImports]
import torch.nn as nn  # pyright: ignore[reportMissingImports]
import torch.nn.functional as F  # pyright: ignore[reportMissingImports]

LOGGER = logging.getLogger(__name__)


def _load_checkpoint_file(checkpoint_path: Path) -> dict:
    """
    Load a checkpoint file with robust error handling and format detection.

    Supports:
    - Standard PyTorch .pt/.pth files
    - Safetensors format (.safetensors)
    - Handles corrupted/incomplete files gracefully

    Raises:
        FileNotFoundError: If checkpoint file doesn't exist
        ValueError: If file format is invalid or corrupted
    """
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")

    suffix = checkpoint_path.suffix.lower()

    # Handle safetensors format
    if suffix == ".safetensors":
        try:
            from safetensors.torch import load_file
            return load_file(str(checkpoint_path))
        except ImportError:
            raise ImportError(
                f"Checkpoint is in safetensors format ({checkpoint_path}), "
                "but safetensors is not installed. Install with: pip install safetensors"
            )
        except Exception as e:
            raise ValueError(f"Failed to load safetensors checkpoint {checkpoint_path}: {e}")

    # Validate file appears to be a valid PyTorch checkpoint
    try:
        with open(checkpoint_path, "rb") as f:
            header = f.read(8)
    except Exception as e:
        raise ValueError(f"Cannot read checkpoint file {checkpoint_path}: {e}")

    # PyTorch checkpoints use pickle which starts with specific bytes
    # Common magic numbers: 0x80 (pickle protocol 2+), ZIP files start with PK
    if len(header) < 2:
        raise ValueError(
            f"Checkpoint file is too small to be valid: {checkpoint_path} ({len(header)} bytes)"
        )

    # Check for obvious non-checkpoint files
    if header[:4] == b"PK\x03\x04":
        # ZIP file - could be a safetensors in disguise or compressed checkpoint
        pass  # Let torch.load try
    elif header[0:1] not in (b"\x80", b"P"):
        # Not pickle protocol or ZIP, check if it might be text
        try:
            sample = header.decode("utf-8", errors="strict")
            if sample.isprintable() or sample.startswith("{"):
                raise ValueError(
                    f"Checkpoint file appears to be a text/JSON file, not a PyTorch checkpoint: "
                    f"{checkpoint_path}\nFirst bytes: {header[:50]!r}"
                )
        except UnicodeDecodeError:
            pass  # Binary file, let torch.load try

    # Attempt to load with torch
    try:
        return torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    except Exception as e:
        error_msg = str(e)
        if "invalid load key" in error_msg:
            raise ValueError(
                f"Checkpoint file is corrupted or not a valid PyTorch checkpoint: "
                f"{checkpoint_path}\n"
                f"This can happen if:\n"
                f"  1. The download was incomplete\n"
                f"  2. The file is in a different format (e.g., safetensors)\n"
                f"  3. The file was corrupted during transfer\n"
                f"Original error: {e}"
            )
        raise ValueError(f"Failed to load checkpoint {checkpoint_path}: {e}")


class AttentionPooling(nn.Module):
    """
    Learned attention-based pooling for sequence aggregation.
    
    Better than max/mean pooling for position-sensitive tasks like
    promoter detection, as it learns which positions are important.
    """
    
    def __init__(self, hidden_dim: int, num_heads: int = 4):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        
        # Query vector is learned, keys/values come from sequence
        self.query = nn.Parameter(torch.randn(1, num_heads, 1, self.head_dim) * 0.02)
        self.key_proj = nn.Linear(hidden_dim, hidden_dim)
        self.value_proj = nn.Linear(hidden_dim, hidden_dim)
        self.out_proj = nn.Linear(hidden_dim, hidden_dim)
        
        self.scale = self.head_dim ** -0.5
    
    def forward(
        self, 
        x: torch.Tensor, 
        mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, D) sequence features
            mask: Optional (B, L) padding mask (True = pad)
        Returns:
            (B, D) pooled representation
        """
        B, L, D = x.shape
        
        # Project to keys and values
        keys = self.key_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        values = self.value_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        
        # Expand query for batch
        query = self.query.expand(B, -1, -1, -1)  # (B, H, 1, head_dim)
        
        # Attention scores
        attn = (query @ keys.transpose(-2, -1)) * self.scale  # (B, H, 1, L)
        
        if mask is not None:
            mask = mask.unsqueeze(1).unsqueeze(2)  # (B, 1, 1, L)
            attn = attn.masked_fill(mask, float('-inf'))
        
        attn = F.softmax(attn, dim=-1)
        
        # Weighted sum
        out = (attn @ values).squeeze(2)  # (B, H, head_dim)
        out = out.reshape(B, D)
        
        return self.out_proj(out)


# =============================================================================
# TCN (Temporal Convolutional Network) Components
# Based on iPro-TCN paper for promoter prediction
# =============================================================================

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
    
    def _init_weights(self):
        """Initialize weights using He initialization."""
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
    
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


# =============================================================================
# Multi-Scale CNN Components
# Based on PROCABLES paper for capturing different motif lengths
# =============================================================================

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
    
    def _init_weights(self):
        """Initialize convolution weights."""
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
    
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


class DNAEncoder(nn.Module):
    """
    Enhanced transformer encoder for DNA k-mer tokens.
    
    Key improvements over original:
    - Pre-norm architecture for training stability
    - Stochastic depth for regularization (linearly increasing drop rate)
    - Scaled embedding initialization (critical for MLM pretraining)
    - Optional gradient checkpointing for memory efficiency
    
    IMPORTANT FOR MLM PRETRAINING:
    - Use scaled initialization for embeddings (std = 0.02)
    - Deeper models (8+ layers) learn better representations
    - Larger FF dimension (4x embedding_dim) helps capacity
    """
    
    def __init__(
        self,
        vocab_size: int,
        kmer: int,
        embedding_dim: int = 192,
        num_layers: int = 6,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        pad_token_id: Optional[int] = None,
        max_seq_len: int = 256,
        drop_path_rate: float = 0.1,  # Stochastic depth rate
        layer_scale_init: float = 1.0,  # Layer scale initialization
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,  # Use GLU-style FFN
        glu_activation: str = "gelu",  # Activation for GLU: "gelu" or "silu"
    ):
        super().__init__()
        self.embedding_dim = embedding_dim
        self.vocab_size = vocab_size
        self.k = kmer
        self.pad_token_id = pad_token_id
        self.num_heads = num_heads
        self.use_relative_position_bias = use_relative_position_bias

        # Embedding with scaled initialization (BERT-style)
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_token_id)
        self._init_embedding()

        self.embed_dropout = nn.Dropout(dropout)

        # Embedding layer norm for stability (like BERT)
        self.embed_norm = nn.LayerNorm(embedding_dim)

        # Positional encoding - ALWAYS use position embeddings for input
        # This is critical for MLM: when all tokens are [MASK], the model needs
        # positional information in the embeddings themselves, not just in attention.
        self.pos_embedding = nn.Embedding(max_seq_len, embedding_dim)
        nn.init.normal_(self.pos_embedding.weight, std=0.02)

        # Relative position bias (T5/ALiBi style) - shared across all layers
        if use_relative_position_bias:
            self.relative_position_bias = RelativePositionBias(
                num_heads=num_heads,
                max_distance=relative_position_max_distance,
                num_buckets=relative_position_num_buckets,
                bidirectional=True,
            )
        else:
            self.relative_position_bias = None

        # Pre-norm transformer layers with linearly increasing stochastic depth
        # This drops more aggressively in later layers (which are more task-specific)
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, num_layers)]
        self.layers = nn.ModuleList([
            PreNormTransformerLayer(
                d_model=embedding_dim,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
                drop_path=dpr[i],
                use_glu_ffn=use_glu_ffn,
                glu_activation=glu_activation,
            )
            for i in range(num_layers)
        ])

        self.final_norm = nn.LayerNorm(embedding_dim)
    
    def _init_embedding(self):
        """
        Initialize embeddings with truncated normal distribution.
        This is critical for good MLM pretraining - standard deviation of 0.02
        is the BERT default and works well for small vocabularies.
        """
        nn.init.trunc_normal_(self.embedding.weight, std=0.02, a=-0.04, b=0.04)
        if self.pad_token_id is not None:
            # Ensure padding embedding is zero
            nn.init.zeros_(self.embedding.weight[self.pad_token_id])
    
    def load_mlm_weights(self, checkpoint: Union[str, Path, dict], strict: bool = False) -> None:
        """Load pretrained MLM encoder weights with flexible matching."""
        if isinstance(checkpoint, (str, Path)):
            checkpoint_path = Path(checkpoint)
            state_dict = _load_checkpoint_file(checkpoint_path)
        else:
            state_dict = checkpoint

        # Handle nested checkpoint format (e.g., {"model_state_dict": ...})
        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]
        elif "state_dict" in state_dict:
            state_dict = state_dict["state_dict"]

        # Try to load what we can
        model_dict = self.state_dict()
        pretrained_dict = {}

        for k, v in state_dict.items():
            if k in model_dict and model_dict[k].shape == v.shape:
                pretrained_dict[k] = v

        model_dict.update(pretrained_dict)
        self.load_state_dict(model_dict)
        LOGGER.info("Loaded %d/%d weights from MLM checkpoint", len(pretrained_dict), len(state_dict))
    
    def freeze_bottom_layers(self, n: int) -> None:
        """Freeze the first n transformer layers."""
        for i, layer in enumerate(self.layers):
            if i < n:
                for param in layer.parameters():
                    param.requires_grad = False
    
    def unfreeze_all_layers(self) -> None:
        """Unfreeze all transformer layers."""
        for layer in self.layers:
            for param in layer.parameters():
                param.requires_grad = True
    
    def forward(
        self, 
        token_ids: torch.LongTensor,
        key_padding_mask: Optional[torch.Tensor] = None,
        return_all_hidden_states: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, List[torch.Tensor]]]:
        """
        Args:
            token_ids: (B, L) token indices
            key_padding_mask: Optional (B, L) boolean mask (True = padding)
            return_all_hidden_states: If True, return tuple of (output, hidden_states_list)
        Returns:
            (B, L, D) encoded representations, optionally with hidden states
        """
        B, L = token_ids.shape
        
        x = self.embedding(token_ids)
        
        # ALWAYS add positional embeddings to input
        # This is critical for MLM: when inputs are [MASK], positions must be encoded
        # in the embeddings, not just in attention.
        positions = torch.arange(L, device=token_ids.device)
        x = x + self.pos_embedding(positions).unsqueeze(0)
        
        # Apply embedding layer norm for stability
        x = self.embed_norm(x)
        x = self.embed_dropout(x)
        
        # Build padding mask if not provided
        if key_padding_mask is None and self.pad_token_id is not None:
            key_padding_mask = token_ids.eq(self.pad_token_id)

        # Compute relative position bias once (shared across all layers)
        position_bias = None
        if self.relative_position_bias is not None:
            position_bias = self.relative_position_bias(L, token_ids.device)

        # Process through layers
        hidden_states = [] if return_all_hidden_states else None
        for layer in self.layers:
            if return_all_hidden_states:
                hidden_states.append(x)
            x = layer(x, key_padding_mask=key_padding_mask, position_bias=position_bias)
        
        output = self.final_norm(x)
        
        if return_all_hidden_states:
            hidden_states.append(output)
            return output, hidden_states
        
        return output


class StochasticDepth(nn.Module):
    """
    Stochastic Depth (Layer Drop) for regularization.

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

        # During training, randomly skip this layer
        keep_prob = 1.0 - self.drop_prob
        if torch.rand(1).item() < self.drop_prob:
            return x  # Skip the residual (drop this layer)

        # Scale residual to compensate for dropped layers
        return x + residual / keep_prob


class GLUFFN(nn.Module):
    """
    Gated Linear Unit Feed-Forward Network.

    Uses gating mechanism: output = Linear(GELU(gate) * up)
    where gate and up are separate linear projections.

    This is used in modern transformers like PaLM, LLaMA, and others
    for improved expressiveness over standard FFN.

    Note: GLU effectively doubles the parameters in the up-projection,
    so we use ff_dim * 2/3 for each branch to keep param count similar
    to standard FFN.
    """

    def __init__(
        self,
        d_model: int,
        ff_dim: int,
        dropout: float = 0.1,
        activation: str = "gelu",  # or "silu" for SwiGLU
    ):
        super().__init__()
        # To maintain similar param count, each branch gets 2/3 of ff_dim
        hidden_dim = int(ff_dim * 2 / 3)
        # Round to multiple of 8 for efficiency
        hidden_dim = ((hidden_dim + 7) // 8) * 8

        self.gate_proj = nn.Linear(d_model, hidden_dim, bias=False)
        self.up_proj = nn.Linear(d_model, hidden_dim, bias=False)
        self.down_proj = nn.Linear(hidden_dim, d_model, bias=False)

        self.dropout = nn.Dropout(dropout)

        if activation == "gelu":
            self.activation = nn.GELU()
        elif activation == "silu":
            self.activation = nn.SiLU()  # SwiGLU uses SiLU
        else:
            self.activation = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (batch, seq_len, d_model)
        Returns:
            Output tensor of shape (batch, seq_len, d_model)
        """
        # Gate branch: activation(gate_proj(x))
        gate = self.activation(self.gate_proj(x))
        # Up branch: up_proj(x)
        up = self.up_proj(x)
        # Gated combination
        hidden = gate * up
        # Down projection
        output = self.down_proj(hidden)
        return self.dropout(output)


class RelativePositionBias(nn.Module):
    """
    Relative position bias for attention, similar to T5/ALiBi.

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

        # Learnable bias for each bucket and head
        self.relative_attention_bias = nn.Embedding(num_buckets, num_heads)

        # Initialize with small values for stable training
        nn.init.normal_(self.relative_attention_bias.weight, std=0.02)

    def _relative_position_bucket(
        self,
        relative_position: torch.Tensor,
    ) -> torch.Tensor:
        """
        Convert relative positions to bucket indices.
        Uses logarithmic bucketing for distant positions (T5-style).

        For bidirectional:
        - Half the buckets are for negative (backward) positions
        - Half are for positive (forward) positions
        - Within each half: first quarter is exact, rest is log-bucketed
        """
        relative_buckets = torch.zeros_like(relative_position)

        if self.bidirectional:
            num_buckets = self.num_buckets
            n = num_buckets // 2
            # Offset for positive positions
            relative_buckets += (relative_position > 0).long() * n
            relative_position = torch.abs(relative_position)
        else:
            # Clamp to non-negative for unidirectional
            relative_position = torch.clamp(-relative_position, min=0)
            num_buckets = self.num_buckets
            n = num_buckets

        # First half of buckets are for exact positions
        max_exact = n // 2
        is_small = relative_position < max_exact

        # Second half uses logarithmic bucketing
        # Maps [max_exact, max_distance] -> [max_exact, n-1]
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
        """
        Compute relative position bias matrix.

        Args:
            seq_len: Sequence length
            device: Device to create tensors on

        Returns:
            Tensor of shape (1, num_heads, seq_len, seq_len)
        """
        # Create position indices
        context_position = torch.arange(seq_len, device=device)[:, None]
        memory_position = torch.arange(seq_len, device=device)[None, :]

        # Relative positions: (seq_len, seq_len)
        # Positive = memory is ahead of context, negative = behind
        relative_position = memory_position - context_position

        # Convert to bucket indices
        relative_position_bucket = self._relative_position_bucket(relative_position)

        # Look up bias values: (seq_len, seq_len, num_heads)
        values = self.relative_attention_bias(relative_position_bucket)

        # Reshape to (1, num_heads, seq_len, seq_len)
        values = values.permute(2, 0, 1).unsqueeze(0)

        return values


class PreNormTransformerLayer(nn.Module):
    """
    Pre-norm transformer layer with stochastic depth.

    IMPORTANT FIX: The attention mask handling was causing interference with learning.
    Now properly handles padding masks without interfering with attention computation.

    Supports optional GLU FFN (Gated Linear Unit) for improved expressiveness.
    """

    def __init__(
        self,
        d_model: int,
        nhead: int,
        dim_feedforward: int,
        dropout: float = 0.1,
        drop_path: float = 0.0,  # Stochastic depth probability
        use_glu_ffn: bool = False,  # Use GLU-style FFN
        glu_activation: str = "gelu",  # Activation for GLU: "gelu" or "silu"
    ):
        super().__init__()
        self.d_model = d_model
        self.nhead = nhead
        self.head_dim = d_model // nhead

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

        # Use custom attention for tighter control over masking behavior
        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)

        self.attn_dropout = nn.Dropout(dropout)
        self.resid_dropout = nn.Dropout(dropout)

        # FFN: either standard or GLU-style
        if use_glu_ffn:
            self.ffn = GLUFFN(
                d_model=d_model,
                ff_dim=dim_feedforward,
                dropout=dropout,
                activation=glu_activation,
            )
        else:
            # Standard FFN
            self.ffn = nn.Sequential(
                nn.Linear(d_model, dim_feedforward),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(dim_feedforward, d_model),
                nn.Dropout(dropout),
            )

        # Stochastic depth for regularization
        self.drop_path_attn = StochasticDepth(drop_path)
        self.drop_path_ffn = StochasticDepth(drop_path)

        self.scale = self.head_dim ** -0.5

        # Initialize projections with scaled initialization
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
    ) -> torch.Tensor:
        B, L, D = x.shape

        # Pre-norm
        normed = self.norm1(x)

        # Project to Q, K, V
        q = self.q_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)
        k = self.k_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)
        v = self.v_proj(normed).view(B, L, self.nhead, self.head_dim).transpose(1, 2)

        # Scaled dot-product attention
        attn_weights = (q @ k.transpose(-2, -1)) * self.scale  # (B, H, L, L)

        # Add relative position bias BEFORE softmax (T5/ALiBi style)
        if position_bias is not None:
            attn_weights = attn_weights + position_bias

        # Apply padding mask
        if key_padding_mask is not None:
            # key_padding_mask: (B, L), True = pad
            # Expand to (B, 1, 1, L) and mask with -inf
            attn_weights = attn_weights.masked_fill(
                key_padding_mask.unsqueeze(1).unsqueeze(2), float('-inf')
            )

        attn_weights = F.softmax(attn_weights, dim=-1)
        attn_weights = self.attn_dropout(attn_weights)
        
        # Apply attention to values
        attn_out = (attn_weights @ v).transpose(1, 2).reshape(B, L, D)
        attn_out = self.out_proj(attn_out)
        
        # Stochastic depth on attention residual
        x = self.drop_path_attn(x, self.resid_dropout(attn_out))
        
        # Pre-norm FFN with residual + stochastic depth
        ffn_out = self.ffn(self.norm2(x))
        x = self.drop_path_ffn(x, ffn_out)
        
        return x


# =============================================================================
# Engineered Features Processing and Attention Fusion
# =============================================================================

class EngineeredFeatureMLP(nn.Module):
    """
    Separate MLP branch for processing engineered features.
    
    This gives engineered features (TNC, PseEIIP) their own
    representation space before fusion with sequence features, allowing
    the model to learn better feature interactions.

    Includes optional pre-norm, gated hidden block, and residual projection
    to stabilize training on larger feature spaces.
    """
    
    def __init__(
        self,
        input_dim: int = 288,
        hidden_dim: int = 512,
        output_dim: int = 256,
        dropout: float = 0.3,
        use_pre_norm: bool = True,
        input_dropout: float = 0.05,
        use_gated: bool = True,
        use_residual: bool = True,
    ):
        super().__init__()
        self.use_pre_norm = use_pre_norm
        self.use_gated = use_gated
        self.use_residual = use_residual

        self.input_norm = nn.LayerNorm(input_dim) if use_pre_norm else None
        self.pre_norm_scale = (
            nn.Parameter(torch.tensor(0.0)) if use_pre_norm else None
        )
        self.input_dropout = nn.Dropout(input_dropout) if input_dropout > 0 else None

        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
            nn.LayerNorm(output_dim),
        )
        if use_gated:
            self.gate_proj = nn.Linear(hidden_dim, hidden_dim)
            # Initialized to a no-op so old checkpoints still load cleanly.
            self.gate_scale = nn.Parameter(torch.tensor(0.0))
        else:
            self.gate_proj = None
            self.gate_scale = None

        if use_residual:
            self.residual_proj = (
                nn.Identity()
                if input_dim == output_dim
                else nn.Linear(input_dim, output_dim)
            )
            # Initialized to a no-op for backwards compatibility.
            self.residual_scale = nn.Parameter(torch.tensor(0.0))
        else:
            self.residual_proj = None
            self.residual_scale = None

        self.output_dim = output_dim
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, input_dim) engineered features
        Returns:
            (B, output_dim) processed features
        """
        residual_in = x
        if self.input_norm is not None:
            normed = self.input_norm(x)
            scale = torch.tanh(self.pre_norm_scale)
            x = x + scale * (normed - x)

        if self.input_dropout is not None:
            x = self.input_dropout(x)

        # Keep layer ordering explicit to insert gated block without renaming weights.
        x = self.mlp[0](x)
        x = self.mlp[1](x)
        x = self.mlp[2](x)
        x = self.mlp[3](x)

        if self.gate_proj is not None:
            gate = torch.sigmoid(self.gate_proj(x))
            scale = torch.tanh(self.gate_scale)
            x = x * (1.0 + scale * (gate - 0.5) * 2.0)

        x = self.mlp[4](x)
        x = self.mlp[5](x)
        x = self.mlp[6](x)
        x = self.mlp[7](x)
        x = self.mlp[8](x)
        x = self.mlp[9](x)

        if self.residual_proj is not None:
            residual = self.residual_proj(residual_in)
            scale = torch.tanh(self.residual_scale)
            x = x + scale * residual

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

        # Linearly increasing stochastic depth
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
    
    def load_pretrained_layers(self, encoder: 'DNAEncoder', num_layers: int) -> int:
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

        B, L, D = x.shape

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


class GeneWhispererStage1(nn.Module):
    """
    Enhanced Stage 1 classifier with improved architecture order and attention fusion.
    
    NEW Architecture: Embedding → CNN/TCN → Transformer → Pooling → Attention Fusion → Classifier
    
    Key improvements:
    1. CNN/TCN FIRST: Extract local patterns before global contextualization
    2. Transformer AFTER: Add global context to local features
    3. Separate MLP for engineered features: Learn dedicated representations
    4. Attention-based fusion: Gated cross-attention instead of concatenation
    
    This architecture processes information in a more logical order:
    - Local patterns (motifs) are extracted first by CNN/TCN
    - Global context is then added by transformer
    - Engineered features are processed separately and fused via attention
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
        encoder: Optional[DNAEncoder] = None,
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
            kmer=kmer,
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
    def encoder(self):
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
        
        # The embedding is already shared (self.embedding.weight = self._full_encoder.embedding.weight)
        
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
        state_dict = _load_checkpoint_file(Path(checkpoint_path))

        # Handle nested state dict (e.g., from training checkpoints)
        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]

        # Get current model state
        model_state = self.state_dict()

        # Filter to only keys that exist in both and have matching shapes
        filtered_state = {}
        for k, v in state_dict.items():
            if k in model_state:
                if model_state[k].shape == v.shape:
                    filtered_state[k] = v
                else:
                    LOGGER.warning(
                        "Shape mismatch for %s: checkpoint %s vs model %s (skipping)",
                        k, v.shape, model_state[k].shape
                    )

        # Log which keys are new (in model but not in checkpoint)
        new_keys = set(model_state.keys()) - set(filtered_state.keys())
        if new_keys:
            # Group by component for cleaner logging
            rel_pos_keys = [k for k in new_keys if "relative_position" in k]
            glu_keys = [k for k in new_keys if "gate_proj" in k or "up_proj" in k or "down_proj" in k]
            other_keys = [k for k in new_keys if k not in rel_pos_keys and k not in glu_keys]

            if rel_pos_keys:
                LOGGER.info("New relative position bias parameters (randomly initialized): %d keys", len(rel_pos_keys))
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
            len(filtered_state), len(state_dict), len(new_keys), len(missing_keys)
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
        B, L, D = embeds.shape
        if L > self.pos_embedding.num_embeddings:
            raise ValueError(
                f"Sequence token length {L} exceeds max_seq_len={self.pos_embedding.num_embeddings}. "
                f"Increase max_bp_len or pass a larger max_seq_len."
            )
        pos_ids = torch.arange(L, device=embeds.device).unsqueeze(0).expand(B, L)
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
        Forward pass with new architecture order.

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


# =============================================================================
# Simplified Stage 1: Transformer-Only Backbone
# =============================================================================

class GeneWhispererSimplified(nn.Module):
    """
    Simplified Stage 1 classifier that keeps the pretrained DNAEncoder intact and
    removes CNN/TCN and post-CNN adapters that underperformed in ablations.

    Architecture rationale:
    - Preserve the full transformer stack (defaults to 12 layers) to retain MLM-pretrained context.
    - Pool directly over transformer outputs (attention or masked mean).
    - Retain the engineered feature MLP for non-sequence signals.
    - Fuse via concatenation to avoid heavier gated attention.
    - Use a small two-layer classifier head to reduce parameter count.

    The constructor accepts the same configuration parameters as
    GeneWhispererStage1 for drop-in compatibility; CNN/TCN-related arguments are
    accepted but ignored in simplified mode. Set use_simplified=False to
    delegate to GeneWhispererStage1 for a full model comparison.
    """

    def __init__(
        self,
        vocab_size: int = 67,
        kmer: int = 3,
        embedding_dim: int = 192,
        num_layers: int = 12,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        pad_token_id: Optional[int] = None,
        encoder: Optional[DNAEncoder] = None,
        engineered_dim: int = 288,
        use_engineered_features: bool = True,
        use_attention_pool: bool = True,
        use_simplified: bool = True,
        # TCN parameters (accepted for compatibility)
        use_tcn: bool = True,
        tcn_hidden: int = 256,
        tcn_levels: int = 4,
        tcn_kernel: int = 3,
        # Multi-scale CNN parameters (accepted for compatibility)
        multiscale_channels: int = 64,
        multiscale_kernels: Tuple[int, ...] = (3, 5, 7, 9, 15),
        # LSTM parameters (accepted for compatibility)
        lstm_hidden: int = 192,
        # Post-CNN transformer parameters (accepted for compatibility)
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

        self.use_simplified = use_simplified
        self._fallback_model: Optional[GeneWhispererStage1] = None

        if not use_simplified:
            self._fallback_model = GeneWhispererStage1(
                vocab_size=vocab_size,
                kmer=kmer,
                embedding_dim=embedding_dim,
                num_layers=num_layers,
                num_heads=num_heads,
                ff_dim=ff_dim,
                dropout=dropout,
                pad_token_id=pad_token_id,
                encoder=encoder,
                engineered_dim=engineered_dim,
                use_engineered_features=use_engineered_features,
                use_attention_pool=use_attention_pool,
                use_tcn=use_tcn,
                tcn_hidden=tcn_hidden,
                tcn_levels=tcn_levels,
                tcn_kernel=tcn_kernel,
                multiscale_channels=multiscale_channels,
                multiscale_kernels=multiscale_kernels,
                lstm_hidden=lstm_hidden,
                post_cnn_transformer_layers=post_cnn_transformer_layers,
                engineered_mlp_hidden=engineered_mlp_hidden,
                engineered_mlp_output=engineered_mlp_output,
                engineered_mlp_pre_norm=engineered_mlp_pre_norm,
                engineered_mlp_input_dropout=engineered_mlp_input_dropout,
                engineered_mlp_use_gated=engineered_mlp_use_gated,
                engineered_mlp_use_residual=engineered_mlp_use_residual,
                fusion_hidden=fusion_hidden,
                drop_path_rate=drop_path_rate,
                max_seq_len=max_seq_len,
                use_relative_position_bias=use_relative_position_bias,
                relative_position_num_buckets=relative_position_num_buckets,
                relative_position_max_distance=relative_position_max_distance,
                use_glu_ffn=use_glu_ffn,
                glu_activation=glu_activation,
            )
            return

        self.pad_token_id = pad_token_id
        self.use_attention_pool = use_attention_pool
        self.engineered_dim = engineered_dim
        self.use_engineered_features = use_engineered_features and engineered_dim > 0

        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            kmer=kmer,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )
        self.embedding = self._full_encoder.embedding
        encoder_dim = self._full_encoder.embedding_dim

        if use_attention_pool:
            if encoder_dim % 4 != 0:
                raise ValueError(
                    f"embedding_dim ({encoder_dim}) must be divisible by 4 for attention pooling"
                )
            self.pool = AttentionPooling(encoder_dim, num_heads=4)
        else:
            self.pool = None

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
            classifier_in = encoder_dim + self.engineered_mlp.output_dim
        else:
            self.engineered_mlp = None
            classifier_in = encoder_dim

        classifier_hidden = encoder_dim
        self.classifier = nn.Sequential(
            nn.Linear(classifier_in, classifier_hidden),
            nn.LayerNorm(classifier_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(classifier_hidden, 1),
        )

    @property
    def encoder(self):
        """Compatibility property for MLM weight loading."""
        if self._fallback_model is not None:
            return self._fallback_model.encoder
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """
        Load pretrained MLM weights into the transformer encoder.

        In simplified mode, transformer-only weights are loaded directly into
        DNAEncoder. The transfer_mode flag is accepted for compatibility with
        GeneWhispererStage1; "embed_plus_adapter" behaves like "embed_only".
        """
        if self._fallback_model is not None:
            self._fallback_model.load_pretrained_weights(
                checkpoint_path, strict=strict, transfer_mode=transfer_mode
            )
            return

        if transfer_mode == "none":
            LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
            return
        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")

        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)
        if transfer_mode == "embed_plus_adapter":
            LOGGER.info("Simplified model ignores adapter transfer; loaded encoder weights only.")

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        """Generate padding mask from tokens."""
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def _encode_tokens(
        self,
        tokens: torch.Tensor,
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        """
        Encode tokens with the full transformer and return token features + mask.
        """
        seq_len = tokens.size(1)
        max_len = self._full_encoder.pos_embedding.num_embeddings
        if seq_len > max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds max_seq_len={max_len}. "
                "Increase max_seq_len or shorten inputs."
            )
        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        return token_features, padding_mask

    def _mean_pool(
        self,
        token_features: torch.Tensor,
        padding_mask: Optional[torch.Tensor],
    ) -> torch.Tensor:
        if padding_mask is None:
            return token_features.mean(dim=1)
        keep_mask = (~padding_mask).unsqueeze(-1).type_as(token_features)
        summed = (token_features * keep_mask).sum(dim=1)
        denom = keep_mask.sum(dim=1).clamp(min=1.0)
        return summed / denom

    def _pool_sequence(
        self,
        token_features: torch.Tensor,
        padding_mask: Optional[torch.Tensor],
    ) -> torch.Tensor:
        if self.use_attention_pool:
            return self.pool(token_features, padding_mask)
        return self._mean_pool(token_features, padding_mask)

    def extract_features(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Extract pooled sequence features and optional fused features.

        Returns:
            (seq_features, fused_features) both shaped (B, D).
        """
        if self._fallback_model is not None:
            return self._fallback_model.extract_features(tokens, engineered_features)

        token_features, padding_mask = self._encode_tokens(tokens)
        seq_features = self._pool_sequence(token_features, padding_mask)

        if engineered_features is not None:
            assert engineered_features.size(-1) == self.engineered_dim, (
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )

        if self.use_engineered_features and engineered_features is not None:
            eng_features = self.engineered_mlp(engineered_features)
            fused = torch.cat([seq_features, eng_features], dim=-1)
            return seq_features, fused

        return seq_features, seq_features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """
        Forward pass for simplified Stage 1 classification.

        Args:
            tokens: (B, L) k-mer token indices
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        return_logits: If True, also return sigmoid probabilities alongside logits
        Returns:
            (B, 1) promoter logits, or tuple of (probs, logits) if return_logits=True
        """
        if self._fallback_model is not None:
            return self._fallback_model(
                tokens, engineered_features=engineered_features, return_logits=return_logits
            )

        _, fused = self.extract_features(tokens, engineered_features)
        logits = self.classifier(fused)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits


# =============================================================================
# Stage 1 Variant Models (Transformer-only, Features-only, Combined)
# =============================================================================

def _mean_pool_sequence(
    token_features: torch.Tensor,
    padding_mask: Optional[torch.Tensor],
) -> torch.Tensor:
    """Mean-pool sequence features with optional padding mask."""
    if padding_mask is None:
        return token_features.mean(dim=1)
    keep_mask = (~padding_mask).unsqueeze(-1).type_as(token_features)
    summed = (token_features * keep_mask).sum(dim=1)
    denom = keep_mask.sum(dim=1).clamp(min=1.0)
    return summed / denom


def _encode_from_embeds(
    encoder: DNAEncoder,
    embeds: torch.Tensor,
    padding_mask: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """Run a DNAEncoder stack starting from pre-computed embeddings."""
    batch_size, seq_len, embed_dim = embeds.shape
    if embed_dim != encoder.embedding_dim:
        raise ValueError(
            f"Embedding dim mismatch: got {embed_dim}, expected {encoder.embedding_dim}"
        )
    if seq_len > encoder.pos_embedding.num_embeddings:
        raise ValueError(
            f"Sequence length {seq_len} exceeds max_seq_len={encoder.pos_embedding.num_embeddings}. "
            "Increase max_seq_len or shorten inputs."
        )

    positions = torch.arange(seq_len, device=embeds.device)
    x = embeds + encoder.pos_embedding(positions).unsqueeze(0)
    x = encoder.embed_norm(x)
    x = encoder.embed_dropout(x)

    position_bias = None
    if encoder.relative_position_bias is not None:
        position_bias = encoder.relative_position_bias(seq_len, embeds.device)

    for layer in encoder.layers:
        x = layer(x, key_padding_mask=padding_mask, position_bias=position_bias)

    return encoder.final_norm(x)


class GeneWhispererTransformerOnly(nn.Module):
    """
    Variant A: Transformer-only Stage 1 model with mean pooling and a 2-layer classifier.

    Parameter count (defaults; vocab_size=4099, embedding_dim=384, num_layers=12,
    num_heads=12, ff_dim=1536, max_seq_len=256): 23,115,649
    """

    def __init__(
        self,
        vocab_size: int = 4099,
        kmer: int = 6,
        embedding_dim: int = 384,
        num_layers: int = 12,
        num_heads: int = 12,
        ff_dim: int = 1536,
        dropout: float = 0.12,
        pad_token_id: Optional[int] = 4098,
        max_seq_len: int = 256,
        encoder: Optional[DNAEncoder] = None,
        classifier_hidden: Optional[int] = None,
        classifier_dropout: Optional[float] = None,
        drop_path_rate: float = 0.1,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
    ):
        super().__init__()
        self.pad_token_id = pad_token_id

        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            kmer=kmer,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )
        self.embedding = self._full_encoder.embedding

        hidden_dim = classifier_hidden or embedding_dim
        classifier_drop = dropout if classifier_dropout is None else classifier_dropout
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(classifier_drop),
            nn.Linear(hidden_dim, 1),
        )

    @property
    def encoder(self) -> DNAEncoder:
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """Load pretrained MLM weights into the transformer encoder."""
        if transfer_mode == "none":
            LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
            return
        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")
        if transfer_mode == "embed_plus_adapter":
            LOGGER.info("Transformer-only variant ignores adapter transfer; loading encoder weights only.")
        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is None:
            return None
        return tokens.eq(self.pad_token_id)

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        pooled = _mean_pool_sequence(token_features, padding_mask)
        logits = self.classifier(pooled)
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
        token_features = _encode_from_embeds(self._full_encoder, embeds, padding_mask)
        pooled = _mean_pool_sequence(token_features, padding_mask)
        logits = self.classifier(pooled)
        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits
        return logits


class GeneWhispererFeaturesOnly(nn.Module):
    """
    Variant B: Engineered-features-only classifier with a 3-layer MLP.

    Parameter count (defaults; engineered_dim=288, hidden_dim=512, second_hidden_dim=384): 345,345
    """

    def __init__(
        self,
        engineered_dim: int = 288,
        hidden_dim: int = 512,
        second_hidden_dim: int = 384,
        dropout: float = 0.12,
    ):
        super().__init__()
        if engineered_dim <= 0:
            raise ValueError(f"engineered_dim must be > 0, got {engineered_dim}")
        self.engineered_dim = engineered_dim

        self.classifier = nn.Sequential(
            nn.Linear(engineered_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, second_hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(second_hidden_dim, 1),
        )

    @property
    def encoder(self) -> Optional[DNAEncoder]:
        return None

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """No-op for features-only model (no transformer weights)."""
        LOGGER.info("Features-only variant has no MLM weights to load; skipping.")

    def _validate_engineered(self, engineered_features: Optional[torch.Tensor]) -> torch.Tensor:
        if engineered_features is None:
            raise ValueError("Engineered features are required for features-only variant.")
        if engineered_features.size(-1) != self.engineered_dim:
            raise ValueError(
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )
        return engineered_features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        engineered = self._validate_engineered(engineered_features)
        logits = self.classifier(engineered)
        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits
        return logits


class GeneWhispererCombined(nn.Module):
    """
    Variant C: Transformer + engineered features with mean pooling and concatenation.

    Parameter count (defaults; vocab_size=4099, embedding_dim=384, num_layers=12,
    engineered_dim=288, engineered_mlp_hidden=512, engineered_mlp_output=384, max_seq_len=256): 23,608,065
    """

    def __init__(
        self,
        vocab_size: int = 4099,
        kmer: int = 6,
        embedding_dim: int = 384,
        num_layers: int = 12,
        num_heads: int = 12,
        ff_dim: int = 1536,
        dropout: float = 0.12,
        pad_token_id: Optional[int] = 4098,
        max_seq_len: int = 256,
        encoder: Optional[DNAEncoder] = None,
        engineered_dim: int = 288,
        engineered_mlp_hidden: int = 512,
        engineered_mlp_output: int = 384,
        classifier_hidden: Optional[int] = None,
        classifier_dropout: Optional[float] = None,
        drop_path_rate: float = 0.1,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
    ):
        super().__init__()
        self.pad_token_id = pad_token_id
        self.engineered_dim = engineered_dim

        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            kmer=kmer,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )
        self.embedding = self._full_encoder.embedding

        self.engineered_mlp = nn.Sequential(
            nn.Linear(engineered_dim, engineered_mlp_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(engineered_mlp_hidden, engineered_mlp_output),
            nn.GELU(),
        )

        hidden_dim = classifier_hidden or embedding_dim
        classifier_drop = dropout if classifier_dropout is None else classifier_dropout
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim + engineered_mlp_output, hidden_dim),
            nn.GELU(),
            nn.Dropout(classifier_drop),
            nn.Linear(hidden_dim, 1),
        )

    @property
    def encoder(self) -> DNAEncoder:
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """Load pretrained MLM weights into the transformer encoder."""
        if transfer_mode == "none":
            LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
            return
        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")
        if transfer_mode == "embed_plus_adapter":
            LOGGER.info("Combined variant ignores adapter transfer; loading encoder weights only.")
        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is None:
            return None
        return tokens.eq(self.pad_token_id)

    def _validate_engineered(self, engineered_features: Optional[torch.Tensor]) -> torch.Tensor:
        if engineered_features is None:
            raise ValueError("Engineered features are required for combined variant.")
        if engineered_features.size(-1) != self.engineered_dim:
            raise ValueError(
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )
        return engineered_features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        pooled = _mean_pool_sequence(token_features, padding_mask)
        engineered = self._validate_engineered(engineered_features)
        engineered_proj = self.engineered_mlp(engineered)
        fused = torch.cat([pooled, engineered_proj], dim=-1)
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
        token_features = _encode_from_embeds(self._full_encoder, embeds, padding_mask)
        pooled = _mean_pool_sequence(token_features, padding_mask)
        engineered = self._validate_engineered(engineered_features)
        engineered_proj = self.engineered_mlp(engineered)
        fused = torch.cat([pooled, engineered_proj], dim=-1)
        logits = self.classifier(fused)
        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits
        return logits


def create_model_variant(variant: str, config: dict) -> nn.Module:
    """
    Factory for transformer-only, features-only, or combined model variants.
    """
    if not isinstance(config, dict):
        raise ValueError("config must be a dict")

    normalized = variant.strip().lower()
    supported = {"transformer_only", "features_only", "combined"}
    if normalized not in supported:
        raise ValueError(f"Unsupported model variant: {variant}")

    vocab_size = int(config.get("vocab_size", 67))
    kmer = int(config.get("kmer", 3))
    embedding_dim = int(config.get("embedding_dim", 192))
    num_layers = int(config.get("transformer_layers", 12))
    num_heads = int(config.get("transformer_heads", 8))
    ff_dim = int(config.get("transformer_ff_dim", embedding_dim * 4))
    dropout = float(config.get("transformer_dropout", 0.1))
    pad_token_id = config.get("pad_token_id")
    if pad_token_id is not None:
        pad_token_id = int(pad_token_id)
    max_seq_len = int(
        config.get(
            "max_seq_len",
            int(config.get("max_bp_len", 81)) - kmer + 1,
        )
    )
    engineered_dim = int(config.get("engineered_dim", 288))
    engineered_mlp_hidden = int(config.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(config.get("engineered_mlp_output", 128))

    drop_path_rate = float(config.get("drop_path_rate", 0.1))
    use_relative_position_bias = bool(config.get("use_relative_position_bias", False))
    relative_position_num_buckets = int(config.get("relative_position_num_buckets", 32))
    relative_position_max_distance = int(config.get("relative_position_max_distance", 128))
    use_glu_ffn = bool(config.get("use_glu_ffn", False))
    glu_activation = str(config.get("glu_activation", "gelu"))

    if normalized == "transformer_only":
        return GeneWhispererTransformerOnly(
            vocab_size=vocab_size,
            kmer=kmer,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )
    if normalized == "features_only":
        return GeneWhispererFeaturesOnly(
            engineered_dim=engineered_dim,
            hidden_dim=engineered_mlp_hidden,
            second_hidden_dim=engineered_mlp_output,
            dropout=dropout,
        )
    return GeneWhispererCombined(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=dropout,
        pad_token_id=pad_token_id,
        max_seq_len=max_seq_len,
        engineered_dim=engineered_dim,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        drop_path_rate=drop_path_rate,
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=relative_position_num_buckets,
        relative_position_max_distance=relative_position_max_distance,
        use_glu_ffn=use_glu_ffn,
        glu_activation=glu_activation,
    )


# =============================================================================
# Stage 2 Heads: Transformer, BiLSTM, and Combined
# =============================================================================

class StrengthTransformerHead(nn.Module):
    """
    Small transformer encoder head for Stage 2 strength classification.

    Applies a lightweight transformer on token-level features followed by
    attention pooling. This captures global dependencies specifically for
    the strong/weak classification task.
    """

    def __init__(
        self,
        input_dim: int,
        n_layers: int = 2,
        n_heads: int = 4,
        ff_dim: Optional[int] = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        ff_dim = ff_dim or (2 * input_dim)

        # Transformer layers
        self.layers = nn.ModuleList([
            PreNormTransformerLayer(
                d_model=input_dim,
                nhead=n_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
                drop_path=0.0,  # No stochastic depth in small head
            )
            for _ in range(n_layers)
        ])

        self.norm = nn.LayerNorm(input_dim)
        self.pool = AttentionPooling(input_dim, num_heads=n_heads)
        self.output_dim = input_dim

    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, D) token-level features
            key_padding_mask: Optional (B, L) padding mask
        Returns:
            (B, D) pooled representation
        """
        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask)

        x = self.norm(x)
        pooled = self.pool(x, key_padding_mask)
        return pooled


class StrengthBiLSTMHead(nn.Module):
    """
    Bidirectional LSTM head for Stage 2 strength classification.

    Captures sequential dependencies with BiLSTM followed by attention pooling.
    BiLSTM can capture different patterns than transformer attention.
    """

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 128,
        num_layers: int = 1,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim

        # BiLSTM (outputs 2 * hidden_dim due to bidirectional)
        self.lstm = nn.LSTM(
            input_size=input_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if num_layers > 1 else 0.0,
        )

        lstm_out_dim = 2 * hidden_dim  # Bidirectional
        self.norm = nn.LayerNorm(lstm_out_dim)
        self.pool = AttentionPooling(lstm_out_dim, num_heads=4)
        self.output_dim = lstm_out_dim

    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, D) token-level features
            key_padding_mask: Optional (B, L) padding mask (not used by LSTM, kept for API consistency)
        Returns:
            (B, 2*hidden_dim) pooled representation
        """
        # LSTM forward
        lstm_out, _ = self.lstm(x)  # (B, L, 2*hidden_dim)
        lstm_out = self.norm(lstm_out)

        # Attention pooling
        pooled = self.pool(lstm_out, key_padding_mask)
        return pooled


class CombinedStrengthHead(nn.Module):
    """
    Combined Transformer + BiLSTM head for Stage 2.

    Runs both heads and concatenates or averages their outputs.
    Order is configurable: 'transformer_first' or 'bilstm_first'.
    """

    def __init__(
        self,
        input_dim: int,
        transformer_layers: int = 2,
        transformer_heads: int = 4,
        transformer_ff_dim: Optional[int] = None,
        bilstm_hidden: int = 128,
        bilstm_layers: int = 1,
        dropout: float = 0.1,
        combine_mode: str = "concat",  # "concat" or "avg"
        order: str = "transformer_first",  # "transformer_first" or "bilstm_first"
    ):
        super().__init__()
        self.combine_mode = combine_mode
        self.order = order

        self.transformer_head = StrengthTransformerHead(
            input_dim=input_dim,
            n_layers=transformer_layers,
            n_heads=transformer_heads,
            ff_dim=transformer_ff_dim,
            dropout=dropout,
        )

        self.bilstm_head = StrengthBiLSTMHead(
            input_dim=input_dim,
            hidden_dim=bilstm_hidden,
            num_layers=bilstm_layers,
            dropout=dropout,
        )

        if combine_mode == "concat":
            self.output_dim = self.transformer_head.output_dim + self.bilstm_head.output_dim
        else:  # avg
            # For averaging, project to common dimension
            self.proj_transformer = nn.Linear(self.transformer_head.output_dim, input_dim)
            self.proj_bilstm = nn.Linear(self.bilstm_head.output_dim, input_dim)
            self.output_dim = input_dim

    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, D) token-level features
            key_padding_mask: Optional (B, L) padding mask
        Returns:
            (B, output_dim) pooled representation
        """
        trans_out = self.transformer_head(x, key_padding_mask)
        bilstm_out = self.bilstm_head(x, key_padding_mask)

        if self.combine_mode == "concat":
            if self.order == "transformer_first":
                return torch.cat([trans_out, bilstm_out], dim=-1)
            else:
                return torch.cat([bilstm_out, trans_out], dim=-1)
        else:  # avg
            trans_proj = self.proj_transformer(trans_out)
            bilstm_proj = self.proj_bilstm(bilstm_out)
            return (trans_proj + bilstm_proj) / 2.0


class GeneWhispererStage2(nn.Module):
    """
    Stage 2 model for strong vs weak promoter classification.

    Reuses the same backbone as GeneWhispererStage1:
    - Embedding layer
    - CNN/TCN for local pattern extraction
    - Post-CNN Transformer for global context

    Instead of immediately pooling, returns token-level features (B, L, D)
    and applies a configurable "strength head":
    - "transformer": StrengthTransformerHead (2 layers, 4 heads)
    - "bilstm": StrengthBiLSTMHead (BiLSTM hidden=128)
    - "both": CombinedStrengthHead (transformer + BiLSTM)

    Final classifier is a small MLP -> logits for binary strong/weak.
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
        encoder: Optional[DNAEncoder] = None,
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
        # LSTM parameters (kept for compatibility)
        lstm_hidden: int = 192,
        # Post-CNN transformer parameters
        post_cnn_transformer_layers: int = 3,
        # Stage 2 head configuration
        stage2_head_type: str = "transformer",  # "transformer", "bilstm", or "both"
        stage2_transformer_layers: int = 2,
        stage2_transformer_heads: int = 4,
        stage2_transformer_ff_dim: Optional[int] = None,
        stage2_bilstm_hidden: int = 128,
        stage2_bilstm_layers: int = 1,
        stage2_combine_mode: str = "concat",  # For "both": "concat" or "avg"
        stage2_combine_order: str = "transformer_first",
        # Engineered features MLP parameters
        engineered_mlp_hidden: int = 512,
        engineered_mlp_output: int = 256,
        engineered_mlp_pre_norm: bool = True,
        engineered_mlp_input_dropout: float = 0.05,
        engineered_mlp_use_gated: bool = True,
        engineered_mlp_use_residual: bool = True,
        fusion_hidden: int = 384,
        # Stochastic depth (drop path) rate for encoder
        drop_path_rate: float = 0.1,
        # Relative position bias parameters
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        # GLU FFN parameters
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
    ):
        super().__init__()

        self.use_tcn = use_tcn
        self.pad_token_id = pad_token_id
        self.stage2_head_type = stage2_head_type

        # Step 1: Embedding layer (from DNAEncoder, but we only use embedding)
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_token_id)
        self.embed_dropout = nn.Dropout(dropout)

        # Store the full encoder for MLM weight loading compatibility
        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            kmer=kmer,
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

        # Step 2: CNN/TCN for local pattern extraction
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

        # Step 3: Transformer adapter for global contextualization
        self.post_cnn_transformer = PostCNNTransformerAdapter(
            input_dim=cnn_out_dim,
            transformer_dim=embedding_dim,
            num_layers=post_cnn_transformer_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
        )

        # Token-level feature dimension after backbone
        self.token_feature_dim = cnn_out_dim

        # Step 4: Stage 2 strength head (configurable)
        if stage2_head_type == "transformer":
            self.strength_head = StrengthTransformerHead(
                input_dim=cnn_out_dim,
                n_layers=stage2_transformer_layers,
                n_heads=stage2_transformer_heads,
                ff_dim=stage2_transformer_ff_dim,
                dropout=dropout,
            )
        elif stage2_head_type == "bilstm":
            self.strength_head = StrengthBiLSTMHead(
                input_dim=cnn_out_dim,
                hidden_dim=stage2_bilstm_hidden,
                num_layers=stage2_bilstm_layers,
                dropout=dropout,
            )
        elif stage2_head_type == "both":
            self.strength_head = CombinedStrengthHead(
                input_dim=cnn_out_dim,
                transformer_layers=stage2_transformer_layers,
                transformer_heads=stage2_transformer_heads,
                transformer_ff_dim=stage2_transformer_ff_dim,
                bilstm_hidden=stage2_bilstm_hidden,
                bilstm_layers=stage2_bilstm_layers,
                dropout=dropout,
                combine_mode=stage2_combine_mode,
                order=stage2_combine_order,
            )
        else:
            raise ValueError(f"Unknown stage2_head_type: {stage2_head_type}. Must be 'transformer', 'bilstm', or 'both'")

        head_out_dim = self.strength_head.output_dim

        # Step 5: Engineered features processing (optional)
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

            # Fusion of sequence head output and engineered features
            self.fusion = GatedAttentionFusion(
                seq_dim=head_out_dim,
                eng_dim=self.engineered_mlp.output_dim,
                hidden_dim=fusion_hidden,
                num_heads=4,
                dropout=dropout,
            )
            classifier_in = self.fusion.output_dim
        else:
            classifier_in = head_out_dim

        # Step 6: Final classifier (small MLP -> logits)
        classifier_dropout = dropout * 1.2
        self.classifier = nn.Sequential(
            nn.Linear(classifier_in, classifier_in // 2),
            nn.LayerNorm(classifier_in // 2),
            nn.GELU(),
            nn.Dropout(classifier_dropout),
            nn.Linear(classifier_in // 2, 1),
        )

    @property
    def encoder(self):
        """Compatibility property for MLM weight loading."""
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """
        Load pretrained MLM weights into the model backbone.

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

        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)

        if transfer_mode == "embed_plus_adapter":
            num_layers = len(self.post_cnn_transformer.layers)
            loaded = self.post_cnn_transformer.load_pretrained_layers(self._full_encoder, num_layers)
            LOGGER.info("Transferred %d transformer layers to post-CNN adapter", loaded)

    def load_stage1_weights(
        self,
        checkpoint_path: Union[str, Path],
        device: Optional[torch.device] = None,
    ) -> int:
        """
        Load Stage 1 checkpoint weights into the backbone.

        Loads compatible backbone weights (embedding, CNN/TCN, post_cnn_transformer)
        and skips classifier/head weights that don't match.

        Args:
            checkpoint_path: Path to Stage 1 checkpoint
            device: Device to load weights to
        Returns:
            Number of tensors successfully loaded
        """
        checkpoint_path = Path(checkpoint_path)
        if not checkpoint_path.exists():
            LOGGER.warning("Stage 1 checkpoint not found: %s", checkpoint_path)
            return 0

        if device is None:
            device = next(self.parameters()).device

        checkpoint = _load_checkpoint_file(checkpoint_path)

        # Handle different checkpoint formats
        if isinstance(checkpoint, dict) and "model_state" in checkpoint:
            state = checkpoint["model_state"]
        elif isinstance(checkpoint, dict):
            state = checkpoint
        else:
            raise ValueError(f"Unsupported checkpoint format: {type(checkpoint)}")

        # Normalize state dict (remove "module." prefix if present)
        if any(k.startswith("module.") for k in state):
            state = {k.removeprefix("module."): v for k, v in state.items()}

        model_state = self.state_dict()

        compatible: dict = {}
        skipped: list = []

        for key, value in state.items():
            target = model_state.get(key)
            if target is None:
                skipped.append(key)
                continue
            if getattr(target, "shape", None) != getattr(value, "shape", None):
                skipped.append(key)
                continue
            compatible[key] = value

        missing = [k for k in model_state.keys() if k not in compatible]
        self.load_state_dict(compatible, strict=False)

        loaded = len(compatible)
        LOGGER.info(
            "Loaded %d compatible tensors from Stage 1 (skipped=%d, missing=%d)",
            loaded,
            len(skipped),
            len(missing),
        )

        if skipped:
            preview = ", ".join(skipped[:5])
            suffix = "…" if len(skipped) > 5 else ""
            LOGGER.info("Skipped keys (preview): %s%s", preview, suffix)

        return loaded

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        """Generate padding mask from tokens."""
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def _backbone_forward(
        self,
        tokens: torch.Tensor,
    ) -> torch.Tensor:
        """
        Process sequence through backbone: embedding → CNN/TCN → Transformer.

        Returns token-level features (B, L, D) instead of pooling.

        Args:
            tokens: (B, L) token indices
        Returns:
            (B, L, D) token-level features
        """
        # Step 1: Embedding
        x = self.embedding(tokens)  # (B, L, embedding_dim)
        x = self.embed_dropout(x)

        # Step 2: CNN/TCN (transpose for conv: B, L, D -> B, D, L)
        x = x.transpose(1, 2)

        if self.use_tcn:
            x = self.feature_extractor(x)  # (B, tcn_hidden, L)
        else:
            x = self.conv(x)
            x = self.conv_bn(x)
            x = self.conv_act(x)

        # Transpose back: (B, D, L) -> (B, L, D)
        x = x.transpose(1, 2)

        # Step 3: Post-CNN Transformer for global context
        padding_mask = self._get_padding_mask(tokens)
        x = self.post_cnn_transformer(x, key_padding_mask=padding_mask)

        return x  # (B, L, D) token-level features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward pass for Stage 2 strength classification.

        Args:
            tokens: (B, L) k-mer token indices
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        Returns:
            (B, 1) strength logits (1 = strong, 0 = weak)
        """
        # Get token-level features from backbone
        token_features = self._backbone_forward(tokens)  # (B, L, D)

        # Apply strength head
        padding_mask = self._get_padding_mask(tokens)
        pooled = self.strength_head(token_features, padding_mask)  # (B, head_out_dim)

        # Fuse with engineered features if available
        if self.use_engineered_features and engineered_features is not None:
            eng_processed = self.engineered_mlp(engineered_features)
            fused = self.fusion(pooled, eng_processed)
        else:
            fused = pooled

        # Final classification (logits)
        return self.classifier(fused)


class MultiScaleEnsemble(nn.Module):
    """
    Multi-scale k-mer ensemble following msBERT approach.

    Trains separate models for k=3,4,5,6 and averages predictions
    via soft voting (probability averaging).

    For CoreML/iPhone deployment, we can either:
    1. Export each model separately and average in Swift
    2. Use a smaller subset (k=3,4) for efficiency
    """

    def __init__(
        self,
        models: List[nn.Module],
        weights: Optional[List[float]] = None,
        learnable_weights: bool = False,
    ):
        super().__init__()
        if not models:
            raise ValueError("MultiScaleEnsemble requires at least one model")
        self.models = nn.ModuleList(models)
        self.learnable_weights = learnable_weights

        if weights is not None and len(weights) != len(models):
            raise ValueError("Weights length must match number of models")

        if learnable_weights:
            if weights is None:
                init_logits = torch.zeros(len(models), dtype=torch.float32)
            else:
                weights_tensor = torch.tensor(weights, dtype=torch.float32)
                weight_sum = weights_tensor.sum()
                if weight_sum.item() == 0.0:
                    raise ValueError("Weights must sum to a non-zero value")
                normed = weights_tensor / weight_sum
                eps = torch.finfo(normed.dtype).eps
                normed = normed.clamp(min=eps)
                init_logits = torch.log(normed)
            self.weights = nn.Parameter(init_logits)
        else:
            if weights is None:
                self.weights = None
            else:
                weights_tensor = torch.tensor(weights, dtype=torch.float32)
                weight_sum = weights_tensor.sum()
                if weight_sum.item() == 0.0:
                    raise ValueError("Weights must sum to a non-zero value")
                self.register_buffer("weights", weights_tensor / weight_sum)

    def forward(self, batch: dict, return_logits: bool = False) -> torch.Tensor:
        """
        Args:
            batch: Dict mapping k-mer size to (tokens, engineered_features)
            return_logits: If True, return logits instead of probabilities
        Returns:
            Ensemble probabilities by default; logits if return_logits=True
        """
        outputs = []
        for model in self.models:
            encoder = getattr(model, "encoder", None)
            kmer = getattr(encoder, "k", None)
            if kmer is None:
                raise ValueError("Model missing encoder k-mer metadata")
            if kmer not in batch:
                raise KeyError(f"No batch entry for k-mer {kmer}")
            tokens, engineered = batch[kmer]
            outputs.append(model(tokens, engineered))

        stacked = torch.stack(outputs, dim=0)
        probs = torch.sigmoid(stacked)

        if self.learnable_weights:
            weights = torch.softmax(self.weights, dim=0)
            weights = weights.to(device=probs.device, dtype=probs.dtype)
            probs = (weights.view(-1, 1, 1) * probs).sum(dim=0)
        elif self.weights is None:
            probs = probs.mean(dim=0)
        else:
            weights = self.weights.to(device=probs.device, dtype=probs.dtype)
            probs = (weights.view(-1, 1, 1) * probs).sum(dim=0)

        if return_logits:
            eps = torch.finfo(probs.dtype).eps
            probs = probs.clamp(min=eps, max=1 - eps)
            return torch.logit(probs)

        return probs


def _post_cnn_checksums(model: GeneWhispererStage1) -> dict:
    checksums = {}
    for name, param in model.post_cnn_transformer.named_parameters():
        checksums[name] = param.detach().float().mean().item()
    return checksums


def _diff_checksums(before: dict, after: dict, atol: float = 1e-12) -> list:
    changed = []
    for name, before_val in before.items():
        after_val = after.get(name)
        if after_val is None:
            changed.append(name)
            continue
        if abs(after_val - before_val) > atol:
            changed.append(name)
    return changed


def _count_parameters(model: nn.Module) -> int:
    return sum(param.numel() for param in model.parameters())


if __name__ == "__main__":
    import argparse

    try:
        import yaml
    except ImportError as exc:
        raise SystemExit("PyYAML is required for the debug utility") from exc

    parser = argparse.ArgumentParser(
        description="Verify MLM transfer modes for post-CNN transformer adapter weights."
    )
    parser.add_argument(
        "--mlm-checkpoint",
        help="Path to pretrained MLM encoder checkpoint",
    )
    parser.add_argument(
        "--config",
        default=str(Path(__file__).resolve().parent / "config.yaml"),
        help="Path to config.yaml for model settings",
    )
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    cfg = {}
    if config_path.exists():
        with config_path.open("r", encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle) or {}
    else:
        LOGGER.warning("Config not found at %s; using model defaults", config_path)

    max_bp_len = cfg.get("max_bp_len")
    if max_bp_len is not None:
        max_seq_len = int(max_bp_len) - int(cfg.get("kmer", 3)) + 1
    else:
        max_seq_len = int(cfg.get("max_seq_len", 256))

    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.15))

    model = GeneWhispererStage1(
        vocab_size=int(cfg.get("vocab_size", 67)),
        kmer=int(cfg.get("kmer", 3)),
        embedding_dim=int(cfg.get("embedding_dim", 192)),
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=cfg.get("pad_token_id"),
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 256)),
        tcn_levels=int(cfg.get("tcn_levels", 4)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        engineered_mlp_pre_norm=bool(cfg.get("engineered_mlp_pre_norm", True)),
        engineered_mlp_input_dropout=float(cfg.get("engineered_mlp_input_dropout", 0.05)),
        engineered_mlp_use_gated=bool(cfg.get("engineered_mlp_use_gated", True)),
        engineered_mlp_use_residual=bool(cfg.get("engineered_mlp_use_residual", True)),
        fusion_hidden=int(cfg.get("fusion_hidden", 256)),
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        max_seq_len=max_seq_len,
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", False)),
        relative_position_num_buckets=int(cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", False)),
        glu_activation=str(cfg.get("glu_activation", "gelu")),
    )

    simplified_model = GeneWhispererSimplified(
        vocab_size=int(cfg.get("vocab_size", 67)),
        kmer=int(cfg.get("kmer", 3)),
        embedding_dim=int(cfg.get("embedding_dim", 192)),
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=cfg.get("pad_token_id"),
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_simplified=True,
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        engineered_mlp_pre_norm=bool(cfg.get("engineered_mlp_pre_norm", True)),
        engineered_mlp_input_dropout=float(cfg.get("engineered_mlp_input_dropout", 0.05)),
        engineered_mlp_use_gated=bool(cfg.get("engineered_mlp_use_gated", True)),
        engineered_mlp_use_residual=bool(cfg.get("engineered_mlp_use_residual", True)),
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        max_seq_len=max_seq_len,
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", False)),
        relative_position_num_buckets=int(cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", False)),
        glu_activation=str(cfg.get("glu_activation", "gelu")),
    )

    stage1_params = _count_parameters(model)
    simplified_params = _count_parameters(simplified_model)
    reduction_pct = 0.0
    if stage1_params > 0:
        reduction_pct = (1.0 - (simplified_params / stage1_params)) * 100.0

    print(f"GeneWhispererStage1 params: {stage1_params:,}")
    print(f"GeneWhispererSimplified params: {simplified_params:,}")
    print(f"Reduction: {reduction_pct:.2f}%")

    if args.mlm_checkpoint is None:
        print("No MLM checkpoint provided; skipping transfer-mode verification.")
        raise SystemExit(0)

    checkpoint_path = Path(args.mlm_checkpoint).expanduser().resolve()
    if not checkpoint_path.exists():
        raise SystemExit(f"Checkpoint not found: {checkpoint_path}")

    before = _post_cnn_checksums(model)
    model.load_pretrained_weights(checkpoint_path, strict=False, transfer_mode="embed_only")
    after_embed_only = _post_cnn_checksums(model)
    embed_only_changes = _diff_checksums(before, after_embed_only)
    if embed_only_changes:
        raise SystemExit(
            "embed_only changed post_cnn_transformer weights: "
            + ", ".join(embed_only_changes[:5])
        )
    print("embed_only: post_cnn_transformer weights unchanged")

    model.load_pretrained_weights(checkpoint_path, strict=False, transfer_mode="embed_plus_adapter")
    after_adapter = _post_cnn_checksums(model)
    adapter_changes = _diff_checksums(after_embed_only, after_adapter)
    if not adapter_changes:
        raise SystemExit("embed_plus_adapter did not change post_cnn_transformer weights")
    print("embed_plus_adapter: post_cnn_transformer weights changed")
