"""DNA encoder and attention pooling modules."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

from .layers import PreNormTransformerLayer, RelativePositionBias

# Import from training directory
import sys
_training_dir = Path(__file__).resolve().parent.parent
if str(_training_dir) not in sys.path:
    sys.path.insert(0, str(_training_dir))
from positional_motif_bias import PositionalMotifBias, DEFAULT_PROMOTER_PRIORS

LOGGER = logging.getLogger(__name__)


def _load_checkpoint_file(checkpoint_path: Path) -> dict:
    """Load a checkpoint file with robust error handling and format detection."""
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")

    suffix = checkpoint_path.suffix.lower()

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

    try:
        with open(checkpoint_path, "rb") as f:
            header = f.read(8)
    except Exception as e:
        raise ValueError(f"Cannot read checkpoint file {checkpoint_path}: {e}")

    if len(header) < 2:
        raise ValueError(
            f"Checkpoint file is too small to be valid: {checkpoint_path} ({len(header)} bytes)"
        )

    if header[:4] == b"PK\x03\x04":
        pass
    elif header[0:1] not in (b"\x80", b"P"):
        try:
            sample = header.decode("utf-8", errors="strict")
            if sample.isprintable() or sample.startswith("{"):
                raise ValueError(
                    f"Checkpoint file appears to be a text/JSON file, not a PyTorch checkpoint: "
                    f"{checkpoint_path}\nFirst bytes: {header[:50]!r}"
                )
        except UnicodeDecodeError:
            pass

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
    """Learned attention-based pooling for sequence aggregation.

    Better than max/mean pooling for position-sensitive tasks like
    promoter detection, as it learns which positions are important.
    """

    def __init__(self, hidden_dim: int, num_heads: int = 4):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads

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

        keys = self.key_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)
        values = self.value_proj(x).view(B, L, self.num_heads, self.head_dim).transpose(1, 2)

        query = self.query.expand(B, -1, -1, -1)

        attn = (query @ keys.transpose(-2, -1)) * self.scale

        if mask is not None:
            mask = mask.unsqueeze(1).unsqueeze(2)
            attn = attn.masked_fill(mask, float('-inf'))

        attn = F.softmax(attn, dim=-1)

        out = (attn @ values).squeeze(2)
        out = out.reshape(B, D)

        return self.out_proj(out)


class DNAEncoder(nn.Module):
    """Enhanced transformer encoder for DNA k-mer tokens.

    Key improvements:
    - Pre-norm architecture for training stability
    - Stochastic depth for regularization (linearly increasing drop rate)
    - Scaled embedding initialization (critical for MLM pretraining)
    - Optional gradient checkpointing for memory efficiency
    """

    def __init__(
        self,
        vocab_size: int,
        embedding_dim: int = 192,
        num_layers: int = 6,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        pad_token_id: Optional[int] = None,
        max_seq_len: int = 256,
        drop_path_rate: float = 0.1,
        layer_scale_init: float = 1.0,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
        use_rope: bool = False,
        rope_base: float = 10000.0,
        ffn_type: str = None,
        norm_type: str = "layernorm",
        ffn_mult: float = None,
        use_positional_motif_bias: bool = False,
        motif_regions: Optional[dict] = None,
    ):
        super().__init__()
        self.embedding_dim = embedding_dim
        self.vocab_size = vocab_size
        self.pad_token_id = pad_token_id
        self.num_heads = num_heads
        self.use_relative_position_bias = use_relative_position_bias
        self.use_rope = use_rope
        self.use_positional_motif_bias = use_positional_motif_bias

        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_token_id)
        self._init_embedding()

        self.embed_dropout = nn.Dropout(dropout)
        self.embed_norm = nn.LayerNorm(embedding_dim)

        if not use_rope:
            self.pos_embedding = nn.Embedding(max_seq_len, embedding_dim)
            nn.init.normal_(self.pos_embedding.weight, std=0.02)
        else:
            self.pos_embedding = nn.Embedding(max_seq_len, embedding_dim)
            nn.init.normal_(self.pos_embedding.weight, std=0.02)

        if use_relative_position_bias:
            self.relative_position_bias = RelativePositionBias(
                num_heads=num_heads,
                max_distance=relative_position_max_distance,
                num_buckets=relative_position_num_buckets,
                bidirectional=True,
            )
        else:
            self.relative_position_bias = None

        if use_positional_motif_bias:
            motif_priors = None
            if motif_regions:
                motif_priors = {}
                for name, (start, end) in motif_regions.items():
                    bias_val = 0.15 if "tss" in name.lower() else 0.1
                    motif_priors[(start, end)] = bias_val
            else:
                motif_priors = DEFAULT_PROMOTER_PRIORS

            self.positional_motif_bias = PositionalMotifBias(
                num_heads=num_heads,
                max_seq_len=max_seq_len,
                init_bias_positions=motif_priors,
                learnable=True,
            )
        else:
            self.positional_motif_bias = None

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
                use_rope=use_rope,
                rope_base=rope_base,
                max_seq_len=max_seq_len,
                ffn_type=ffn_type,
                norm_type=norm_type,
                ffn_mult=ffn_mult,
            )
            for i in range(num_layers)
        ])

        self.final_norm = nn.LayerNorm(embedding_dim)

    def _init_embedding(self):
        """Initialize embeddings with truncated normal distribution."""
        nn.init.trunc_normal_(self.embedding.weight, std=0.02, a=-0.04, b=0.04)
        if self.pad_token_id is not None:
            nn.init.zeros_(self.embedding.weight[self.pad_token_id])

    def load_mlm_weights(self, checkpoint: Union[str, Path, dict], strict: bool = False) -> None:
        """Load pretrained MLM encoder weights with flexible matching."""
        if isinstance(checkpoint, (str, Path)):
            checkpoint_path = Path(checkpoint)
            state_dict = _load_checkpoint_file(checkpoint_path)
        else:
            state_dict = checkpoint

        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]
        elif "state_dict" in state_dict:
            state_dict = state_dict["state_dict"]

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

        if not self.use_rope:
            positions = torch.arange(L, device=token_ids.device)
            x = x + self.pos_embedding(positions).unsqueeze(0)

        x = self.embed_norm(x)
        x = self.embed_dropout(x)

        if key_padding_mask is None and self.pad_token_id is not None:
            key_padding_mask = token_ids.eq(self.pad_token_id)

        position_bias = None
        if self.relative_position_bias is not None:
            position_bias = self.relative_position_bias(L, token_ids.device)

        motif_bias = None
        if self.positional_motif_bias is not None:
            motif_bias = self.positional_motif_bias.bias[:, :L, :L].unsqueeze(0)

        hidden_states = [] if return_all_hidden_states else None
        for layer in self.layers:
            if return_all_hidden_states:
                hidden_states.append(x)
            x = layer(
                x,
                key_padding_mask=key_padding_mask,
                position_bias=position_bias,
                motif_bias=motif_bias,
            )

        output = self.final_norm(x)

        if return_all_hidden_states:
            hidden_states.append(output)
            return output, hidden_states

        return output
