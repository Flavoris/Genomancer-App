from __future__ import annotations

import sys
from pathlib import Path

import torch

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.models.mlm_model import DNAMLMModel, MLMConfig
from gene_whisperer.models.transformer import TransformerConfig


def _build_transformer_config() -> TransformerConfig:
    return TransformerConfig(
        vocab_size=128,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        max_seq_len=16,
        dropout=0.1,
        pad_token_id=0,
    )


def test_mlm_model_ties_embeddings_and_emits_logits() -> None:
    model = DNAMLMModel(MLMConfig(transformer=_build_transformer_config(), tie_weights=True))
    input_ids = torch.randint(1, 127, (2, 16))
    attention_mask = torch.ones((2, 16), dtype=torch.long)

    logits = model(input_ids, attention_mask=attention_mask)

    assert model.lm_head.weight.data_ptr() == model.encoder.token_embed.weight.data_ptr()
    assert logits.shape == (2, 16, 128)
    assert torch.isfinite(logits).all()


def test_mlm_model_initializes_untied_head_bias_to_zero() -> None:
    model = DNAMLMModel(MLMConfig(transformer=_build_transformer_config(), tie_weights=False))

    assert model.lm_head.weight.data_ptr() != model.encoder.token_embed.weight.data_ptr()
    assert torch.count_nonzero(model.lm_bias) == 0
