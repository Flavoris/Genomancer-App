"""Promoter classifier with transformer + engineered feature fusion."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn

from gene_whisperer.models.transformer import TransformerConfig, TransformerEncoder


@dataclass
class PromoterConfig:
    transformer: TransformerConfig
    engineered_dim: int = 0
    engineered_hidden: int = 128
    fusion_hidden: int = 256
    dropout: float = 0.1
    use_engineered_features: bool = True


class PromoterModel(nn.Module):
    """Shared encoder with two heads: promoter vs non, strong vs weak."""

    def __init__(self, config: PromoterConfig) -> None:
        super().__init__()
        self.config = config
        self.encoder = TransformerEncoder(config.transformer)
        self.dropout = nn.Dropout(config.dropout)

        if config.use_engineered_features and config.engineered_dim > 0:
            self.engineered_proj = nn.Sequential(
                nn.LayerNorm(config.engineered_dim),
                nn.Linear(config.engineered_dim, config.engineered_hidden),
                nn.GELU(),
                nn.Dropout(config.dropout),
            )
            fusion_input = config.transformer.embedding_dim + config.engineered_hidden
        else:
            self.engineered_proj = None
            fusion_input = config.transformer.embedding_dim

        self.fusion = nn.Sequential(
            nn.Linear(fusion_input, config.fusion_hidden),
            nn.GELU(),
            nn.Dropout(config.dropout),
        )
        self.stage1_head = nn.Linear(config.fusion_hidden, 1)
        self.stage2_head = nn.Linear(config.fusion_hidden, 1)

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
        engineered_features: Optional[torch.Tensor] = None,
        task: str = "stage1",
    ) -> torch.Tensor:
        hidden = self.encoder(input_ids, attention_mask=attention_mask)
        pooled = self.encoder.pool(hidden, attention_mask=attention_mask)
        pooled = self.dropout(pooled)

        if self.engineered_proj is not None:
            if engineered_features is None:
                raise ValueError("engineered_features required when use_engineered_features is True")
            engineered_out = self.engineered_proj(engineered_features)
            pooled = torch.cat([pooled, engineered_out], dim=-1)

        fused = self.fusion(pooled)

        if task == "stage1":
            return self.stage1_head(fused)
        if task == "stage2":
            return self.stage2_head(fused)
        raise ValueError(f"Unknown task: {task}")
