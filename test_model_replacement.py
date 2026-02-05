"""Smoke tests for the new Gene Whisperer models."""
from __future__ import annotations

import sys
from pathlib import Path

import torch

ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.models.mlm_model import DNAMLMModel, MLMConfig
from gene_whisperer.models.promoter_model import PromoterConfig, PromoterModel
from gene_whisperer.models.transformer import TransformerConfig


def _count_parameters(model: torch.nn.Module) -> int:
    return sum(param.numel() for param in model.parameters())


def test_promoter_forward_shapes() -> None:
    transformer_cfg = TransformerConfig(
        vocab_size=256,
        embedding_dim=64,
        num_layers=2,
        num_heads=4,
        ff_dim=128,
        max_seq_len=32,
        dropout=0.1,
        pad_token_id=0,
    )
    model = PromoterModel(
        PromoterConfig(
            transformer=transformer_cfg,
            engineered_dim=16,
            engineered_hidden=8,
            fusion_hidden=32,
            dropout=0.1,
            use_engineered_features=True,
        )
    )
    input_ids = torch.randint(0, 255, (2, 32))
    attention_mask = torch.ones((2, 32), dtype=torch.long)
    engineered = torch.randn(2, 16)

    stage1_logits = model(
        input_ids,
        attention_mask=attention_mask,
        engineered_features=engineered,
        task="stage1",
    )
    stage2_logits = model(
        input_ids,
        attention_mask=attention_mask,
        engineered_features=engineered,
        task="stage2",
    )

    assert stage1_logits.shape == (2, 1)
    assert stage2_logits.shape == (2, 1)


def test_mlm_forward_shape() -> None:
    transformer_cfg = TransformerConfig(
        vocab_size=128,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        max_seq_len=16,
        dropout=0.1,
        pad_token_id=0,
    )
    model = DNAMLMModel(MLMConfig(transformer=transformer_cfg))
    input_ids = torch.randint(0, 127, (2, 16))
    attention_mask = torch.ones((2, 16), dtype=torch.long)
    logits = model(input_ids, attention_mask=attention_mask)
    assert logits.shape == (2, 16, 128)


def test_parameter_budget_mobile_friendly() -> None:
    transformer_cfg = TransformerConfig(
        vocab_size=4096,
        embedding_dim=256,
        num_layers=6,
        num_heads=8,
        ff_dim=1024,
        max_seq_len=128,
        dropout=0.1,
        pad_token_id=0,
    )
    model = PromoterModel(
        PromoterConfig(
            transformer=transformer_cfg,
            engineered_dim=200,
            engineered_hidden=128,
            fusion_hidden=256,
            dropout=0.1,
            use_engineered_features=True,
        )
    )
    total_params = _count_parameters(model)
    assert total_params < 12_000_000, (
        f"Model too large for mobile target: {total_params:,} params"
    )
