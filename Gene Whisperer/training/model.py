"""Enhanced model definitions for Gene Whisperer with TCN and multi-scale CNN.

Key improvements from DNABERT, msBERT, Procables, and iPro-TCN papers:
1. Multi-scale k-mer encoding (k=3,4,5,6) for capturing patterns at different scales
2. Temporal Convolutional Network (TCN) with dilated causal convolutions
3. Multi-scale parallel CNN for capturing different motif lengths (3-15bp)
4. Attention pooling instead of max pooling for position-aware aggregation
5. ALiBi-style relative positional encoding for better generalization
6. CoreML-compatible design for iPhone deployment
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import List, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

LOGGER = logging.getLogger(__name__)


class ALiBiPositionalBias(nn.Module):
    """
    ALiBi (Attention with Linear Biases) positional encoding.
    
    More memory-efficient than learned positional embeddings and
    generalizes better to different sequence lengths. Works well
    for DNA sequences where relative position matters.
    """
    
    def __init__(self, num_heads: int, max_seq_len: int = 512):
        super().__init__()
        self.num_heads = num_heads
        # Slopes decrease geometrically for each head
        slopes = self._get_slopes(num_heads)
        self.register_buffer("slopes", torch.tensor(slopes).view(num_heads, 1, 1))
        
        # Pre-compute relative position matrix
        positions = torch.arange(max_seq_len)
        rel_pos = positions.unsqueeze(0) - positions.unsqueeze(1)  # (L, L)
        self.register_buffer("rel_pos_template", rel_pos.float())
    
    def _get_slopes(self, n: int) -> List[float]:
        """Generate head-specific slopes using geometric series."""
        def get_slopes_power_of_2(n: int) -> List[float]:
            start = 2 ** (-(2 ** -(math.log2(n) - 3)))
            ratio = start
            return [start * (ratio ** i) for i in range(n)]
        
        if math.log2(n).is_integer():
            return get_slopes_power_of_2(n)
        else:
            closest_power = 2 ** math.floor(math.log2(n))
            return (
                get_slopes_power_of_2(closest_power) +
                self._get_slopes(2 * closest_power)[0::2][:n - closest_power]
            )
    
    def forward(self, seq_len: int) -> torch.Tensor:
        """
        Returns ALiBi bias matrix of shape (num_heads, seq_len, seq_len).
        """
        rel_pos = self.rel_pos_template[:seq_len, :seq_len]  # (L, L)
        alibi = self.slopes * rel_pos.unsqueeze(0)  # (H, L, L)
        return alibi


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
    - ALiBi positional encoding for better length generalization
    - Pre-norm architecture for training stability
    - Optional gradient checkpointing for memory efficiency
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
        use_alibi: bool = True,
        pad_token_id: Optional[int] = None,
        max_seq_len: int = 256,
    ):
        super().__init__()
        self.embedding_dim = embedding_dim
        self.k = kmer
        self.pad_token_id = pad_token_id
        self.num_heads = num_heads
        self.use_alibi = use_alibi
        
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_token_id)
        self.embed_dropout = nn.Dropout(dropout)
        
        # ALiBi for relative position encoding
        if use_alibi:
            self.alibi = ALiBiPositionalBias(num_heads, max_seq_len)
        else:
            # Fallback to learned positional embedding
            self.pos_embedding = nn.Embedding(max_seq_len, embedding_dim)
            self.alibi = None
        
        # Pre-norm transformer layers
        self.layers = nn.ModuleList([
            PreNormTransformerLayer(
                d_model=embedding_dim,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
            )
            for _ in range(num_layers)
        ])
        
        self.final_norm = nn.LayerNorm(embedding_dim)
    
    def load_mlm_weights(self, checkpoint: Union[str, Path, dict], strict: bool = False) -> None:
        """Load pretrained MLM encoder weights with flexible matching."""
        if isinstance(checkpoint, (str, Path)):
            state_dict = torch.load(checkpoint, map_location="cpu")
        else:
            state_dict = checkpoint
        
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
    ) -> torch.Tensor:
        """
        Args:
            token_ids: (B, L) token indices
            key_padding_mask: Optional (B, L) boolean mask (True = padding)
        Returns:
            (B, L, D) encoded representations
        """
        B, L = token_ids.shape
        
        x = self.embedding(token_ids)
        
        if not self.use_alibi and hasattr(self, 'pos_embedding'):
            positions = torch.arange(L, device=token_ids.device)
            x = x + self.pos_embedding(positions).unsqueeze(0)
        
        x = self.embed_dropout(x)
        
        # Build padding mask if not provided
        if key_padding_mask is None and self.pad_token_id is not None:
            key_padding_mask = token_ids.eq(self.pad_token_id)
        
        # Get ALiBi bias
        attn_bias = None
        if self.use_alibi:
            attn_bias = self.alibi(L)  # (H, L, L)
        
        # Process through layers
        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask, attn_bias=attn_bias)
        
        return self.final_norm(x)


class PreNormTransformerLayer(nn.Module):
    """Pre-norm transformer layer with optional ALiBi bias."""
    
    def __init__(
        self,
        d_model: int,
        nhead: int,
        dim_feedforward: int,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        
        self.self_attn = nn.MultiheadAttention(
            d_model, nhead, dropout=dropout, batch_first=True
        )
        
        self.ffn = nn.Sequential(
            nn.Linear(d_model, dim_feedforward),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim_feedforward, d_model),
            nn.Dropout(dropout),
        )
        
        self.dropout = nn.Dropout(dropout)
    
    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
        attn_bias: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        # Pre-norm self-attention with residual
        normed = self.norm1(x)
        B, L, _ = x.shape
        
        # For ALiBi, we need to add the bias to attention logits
        # PyTorch's MHA doesn't directly support this, so we use attn_mask
        attn_mask = None
        if attn_bias is not None:
            # attn_bias is (H, L, L), we need (B*H, L, L) for MHA
            attn_mask = attn_bias.unsqueeze(0).expand(B, -1, -1, -1)
            attn_mask = attn_mask.reshape(B * attn_bias.size(0), attn_bias.size(1), attn_bias.size(2))
        
        # Convert key_padding_mask to float to match attn_mask type
        # This avoids the deprecation warning about mismatched mask types
        float_key_padding_mask = None
        if key_padding_mask is not None:
            float_key_padding_mask = key_padding_mask.float().masked_fill(key_padding_mask, float('-inf'))
        
        attn_out, _ = self.self_attn(
            normed, normed, normed,
            key_padding_mask=float_key_padding_mask,
            attn_mask=attn_mask,
            need_weights=False,
        )
        x = x + self.dropout(attn_out)
        
        # Pre-norm FFN with residual
        x = x + self.ffn(self.norm2(x))
        
        return x


# =============================================================================
# Engineered Features Processing and Attention Fusion
# =============================================================================

class EngineeredFeatureMLP(nn.Module):
    """
    Separate MLP branch for processing engineered features.
    
    This gives engineered features (TNC, PSTNP, PseEIIP, CKSNAP) their own
    representation space before fusion with sequence features, allowing
    the model to learn better feature interactions.
    """
    
    def __init__(
        self,
        input_dim: int = 208,
        hidden_dim: int = 256,
        output_dim: int = 128,
        dropout: float = 0.3,
    ):
        super().__init__()
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
        self.output_dim = output_dim
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, input_dim) engineered features
        Returns:
            (B, output_dim) processed features
        """
        return self.mlp(x)


class GatedAttentionFusion(nn.Module):
    """
    Attention-based fusion of sequence and engineered features.
    
    Uses cross-attention where engineered features attend to sequence features,
    plus a gating mechanism to control the contribution of each modality.
    This is more expressive than simple concatenation.
    """
    
    def __init__(
        self,
        seq_dim: int,
        eng_dim: int,
        hidden_dim: int = 256,
        num_heads: int = 4,
        dropout: float = 0.2,
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
        
        # Gating mechanism to balance sequence vs engineered contributions
        self.gate_seq = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.Sigmoid(),
        )
        self.gate_eng = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.Sigmoid(),
        )
        
        # Final projection after fusion
        self.out_proj = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        
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
        B = seq_features.size(0)
        
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
        
        # Apply gates
        gated_seq = g_seq * seq_h
        gated_eng = g_eng * (eng_h + attn_out)  # Add attention-enhanced eng features
        
        # Fuse
        fused = torch.cat([gated_seq, gated_eng], dim=1)  # (B, hidden_dim * 2)
        return self.out_proj(fused)


class LightweightTransformerEncoder(nn.Module):
    """
    Lightweight transformer encoder for post-CNN/TCN contextualization.
    
    This processes the CNN/TCN features to add global context via self-attention.
    Uses fewer layers than the full DNAEncoder since CNN/TCN already extracted
    local features.
    """
    
    def __init__(
        self,
        d_model: int,
        num_layers: int = 3,
        num_heads: int = 4,
        ff_dim: int = 256,
        dropout: float = 0.15,
        use_alibi: bool = True,
        max_seq_len: int = 256,
    ):
        super().__init__()
        self.d_model = d_model
        self.use_alibi = use_alibi
        
        if use_alibi:
            self.alibi = ALiBiPositionalBias(num_heads, max_seq_len)
        else:
            self.pos_embedding = nn.Embedding(max_seq_len, d_model)
            self.alibi = None
        
        self.layers = nn.ModuleList([
            PreNormTransformerLayer(
                d_model=d_model,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                dropout=dropout,
            )
            for _ in range(num_layers)
        ])
        
        self.final_norm = nn.LayerNorm(d_model)
    
    def forward(
        self,
        x: torch.Tensor,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, L, D) features from CNN/TCN
            key_padding_mask: Optional (B, L) padding mask
        Returns:
            (B, L, D) contextualized features
        """
        B, L, D = x.shape
        
        if not self.use_alibi and hasattr(self, 'pos_embedding'):
            positions = torch.arange(L, device=x.device)
            x = x + self.pos_embedding(positions).unsqueeze(0)
        
        attn_bias = None
        if self.use_alibi:
            attn_bias = self.alibi(L)
        
        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask, attn_bias=attn_bias)
        
        return self.final_norm(x)


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
        use_alibi: bool = True,
        pad_token_id: Optional[int] = None,
        encoder: Optional[DNAEncoder] = None,
        engineered_dim: int = 208,
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
        engineered_mlp_hidden: int = 256,
        engineered_mlp_output: int = 128,
        fusion_hidden: int = 256,
    ):
        super().__init__()
        
        self.use_tcn = use_tcn
        self.pad_token_id = pad_token_id
        
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
            use_alibi=use_alibi,
            pad_token_id=pad_token_id,
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
        
        # Step 3: Lightweight transformer for global contextualization (AFTER CNN)
        self.post_cnn_transformer = LightweightTransformerEncoder(
            d_model=cnn_out_dim,
            num_layers=post_cnn_transformer_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            use_alibi=use_alibi,
        )
        
        # Step 4: Attention pooling
        self.use_attention_pool = use_attention_pool
        if use_attention_pool:
            self.pool = AttentionPooling(cnn_out_dim, num_heads=4)
        
        seq_out_dim = cnn_out_dim
        
        # Step 5: Engineered features MLP branch
        self.use_engineered_features = use_engineered_features and engineered_dim > 0
        self.engineered_dim = engineered_dim if self.use_engineered_features else 0
        
        if self.use_engineered_features:
            self.engineered_mlp = EngineeredFeatureMLP(
                input_dim=engineered_dim,
                hidden_dim=engineered_mlp_hidden,
                output_dim=engineered_mlp_output,
                dropout=dropout,
            )
            eng_out_dim = engineered_mlp_output
            
            # Step 6: Attention-based fusion
            self.fusion = GatedAttentionFusion(
                seq_dim=seq_out_dim,
                eng_dim=eng_out_dim,
                hidden_dim=fusion_hidden,
                num_heads=4,
                dropout=dropout,
            )
            classifier_in = fusion_hidden
        else:
            classifier_in = seq_out_dim
        
        # Step 7: Classifier head
        hidden_dim = 512 if use_tcn else 384
        self.classifier = nn.Sequential(
            nn.Linear(classifier_in, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.4),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid(),
        )
    
    @property
    def encoder(self):
        """Compatibility property for MLM weight loading."""
        return self._full_encoder
    
    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        """Generate padding mask from tokens."""
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None
    
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
        
        # Step 4: Pooling
        if self.use_attention_pool:
            pooled = self.pool(x, padding_mask)
        else:
            pooled = torch.max(x, dim=1).values
        
        return pooled
    
    def extract_features(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Extract sequence and fused features."""
        seq_features = self._sequence_backbone(tokens)
        
        if self.use_engineered_features and engineered_features is not None:
            eng_features = self.engineered_mlp(engineered_features)
            fused = self.fusion(seq_features, eng_features)
            return seq_features, fused
        
        return seq_features, seq_features
    
    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward pass with new architecture order.
        
        Args:
            tokens: (B, L) k-mer token indices
            engineered_features: Optional (B, engineered_dim) hand-crafted features
        Returns:
            (B, 1) promoter probability
        """
        _, fused = self.extract_features(tokens, engineered_features)
        return self.classifier(fused)


class MultiScaleEnsemble(nn.Module):
    """
    Multi-scale k-mer ensemble following msBERT approach.
    
    Trains separate models for k=3,4,5,6 and averages predictions.
    This captures patterns at different scales - crucial for 
    improving beyond 80% accuracy.
    
    For CoreML/iPhone deployment, we can either:
    1. Export each model separately and average in Swift
    2. Use a smaller subset (k=3,4) for efficiency
    """
    
    def __init__(self, models: List[nn.Module]):
        super().__init__()
        if not models:
            raise ValueError("MultiScaleEnsemble requires at least one model")
        self.models = nn.ModuleList(models)
    
    def forward(self, batch: dict) -> torch.Tensor:
        """
        Args:
            batch: Dict mapping k-mer size to (tokens, engineered_features)
        Returns:
            Averaged predictions
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
        return stacked.mean(dim=0)


