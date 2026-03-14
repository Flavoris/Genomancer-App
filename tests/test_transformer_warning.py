from __future__ import annotations

import sys
import warnings
from pathlib import Path

import torch

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.models.transformer import TransformerConfig, TransformerEncoder


def test_transformer_encoder_avoids_nested_tensor_norm_first_warning() -> None:
    config = TransformerConfig(
        vocab_size=128,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        max_seq_len=16,
        dropout=0.1,
        pad_token_id=0,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = TransformerEncoder(config)
        assert hasattr(model, "embed_norm")
        assert torch.count_nonzero(model.token_embed.weight[config.pad_token_id]) == 0
        input_ids = torch.randint(0, 127, (2, 16))
        attention_mask = torch.ones((2, 16), dtype=torch.long)
        _ = model(input_ids, attention_mask=attention_mask)

    nested_tensor_warning = [
        warning
        for warning in caught
        if "enable_nested_tensor is True" in str(warning.message)
    ]
    assert not nested_tensor_warning
