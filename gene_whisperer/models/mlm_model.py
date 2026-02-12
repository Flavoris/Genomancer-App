"""Masked language model head for DNA BPE tokens."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn

from gene_whisperer.models.transformer import TransformerConfig, TransformerEncoder


@dataclass
class MLMConfig:
    transformer: TransformerConfig
    tie_weights: bool = True
    head_dropout: float = 0.1


class DNAMLMModel(nn.Module):
    """BERT-style MLM model for DNA sequences."""

    def __init__(self, config: MLMConfig) -> None:
        super().__init__()
        self.config = config
        self.encoder = TransformerEncoder(config.transformer)
        self.mlm_dense = nn.Linear(
            config.transformer.embedding_dim,
            config.transformer.embedding_dim,
        )
        self.mlm_activation = nn.GELU()
        self.mlm_norm = nn.LayerNorm(config.transformer.embedding_dim)
        self.mlm_dropout = nn.Dropout(config.head_dropout)
        self.lm_head = nn.Linear(
            config.transformer.embedding_dim, config.transformer.vocab_size, bias=False
        )
        self.lm_bias = nn.Parameter(torch.zeros(config.transformer.vocab_size))
        if config.tie_weights:
            self.lm_head.weight = self.encoder.token_embed.weight

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        hidden = self.encoder(input_ids, attention_mask=attention_mask)
        hidden = self.mlm_dense(hidden)
        hidden = self.mlm_activation(hidden)
        hidden = self.mlm_norm(hidden)
        hidden = self.mlm_dropout(hidden)
        return self.lm_head(hidden) + self.lm_bias
