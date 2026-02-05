"""Transformer encoder backbone used for DNA sequence modeling."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn


@dataclass
class TransformerConfig:
    vocab_size: int
    embedding_dim: int = 256
    num_layers: int = 6
    num_heads: int = 8
    ff_dim: int = 1024
    max_seq_len: int = 256
    dropout: float = 0.1
    pad_token_id: int = 0


class TransformerEncoder(nn.Module):
    """BERT-style transformer encoder with learned positional embeddings."""

    def __init__(self, config: TransformerConfig) -> None:
        super().__init__()
        self.config = config
        self.token_embed = nn.Embedding(
            config.vocab_size, config.embedding_dim, padding_idx=config.pad_token_id
        )
        self.pos_embed = nn.Embedding(config.max_seq_len, config.embedding_dim)
        self.embed_dropout = nn.Dropout(config.dropout)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=config.embedding_dim,
            nhead=config.num_heads,
            dim_feedforward=config.ff_dim,
            dropout=config.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=config.num_layers)
        self.final_norm = nn.LayerNorm(config.embedding_dim)

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        batch_size, seq_len = input_ids.shape
        positions = torch.arange(seq_len, device=input_ids.device).unsqueeze(0)
        positions = positions.expand(batch_size, seq_len)

        hidden = self.token_embed(input_ids) + self.pos_embed(positions)
        hidden = self.embed_dropout(hidden)

        if attention_mask is None:
            key_padding_mask = None
        else:
            key_padding_mask = attention_mask == 0

        hidden = self.encoder(hidden, src_key_padding_mask=key_padding_mask)
        return self.final_norm(hidden)

    def pool(self, hidden: torch.Tensor, attention_mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """Mean pooling over non-padding tokens."""
        if attention_mask is None:
            return hidden.mean(dim=1)
        mask = attention_mask.unsqueeze(-1).float()
        summed = (hidden * mask).sum(dim=1)
        denom = mask.sum(dim=1).clamp(min=1.0)
        return summed / denom
