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


class DNAMLMModel(nn.Module):
    """BERT-style MLM model for DNA sequences."""

    def __init__(self, config: MLMConfig) -> None:
        super().__init__()
        self.config = config
        self.encoder = TransformerEncoder(config.transformer)
        self.lm_head = nn.Linear(
            config.transformer.embedding_dim, config.transformer.vocab_size, bias=False
        )
        self.output_norm = nn.LayerNorm(config.transformer.embedding_dim)
        if config.tie_weights:
            self.lm_head.weight = self.encoder.token_embed.weight

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        hidden = self.encoder(input_ids, attention_mask=attention_mask)
        hidden = self.output_norm(hidden)
        return self.lm_head(hidden)
