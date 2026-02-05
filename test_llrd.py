"""LLRD sanity checks for the new Gene Whisperer model."""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.models.promoter_model import PromoterConfig, PromoterModel
from gene_whisperer.models.transformer import TransformerConfig
from gene_whisperer.training.train_utils import get_parameter_groups_with_llrd


def _get_group_lr(param_groups: list[dict], group_name: str) -> float:
    for group in param_groups:
        if group.get("name") == group_name:
            return float(group["lr"])
    raise AssertionError(f"Missing param group: {group_name}")


def test_llrd_param_groups() -> None:
    base_lr = 2e-5
    layer_decay = 0.9
    weight_decay = 0.01
    num_layers = 3

    transformer_cfg = TransformerConfig(
        vocab_size=128,
        embedding_dim=32,
        num_layers=num_layers,
        num_heads=4,
        ff_dim=64,
        max_seq_len=32,
        dropout=0.1,
        pad_token_id=0,
    )
    model = PromoterModel(
        PromoterConfig(
            transformer=transformer_cfg,
            engineered_dim=0,
            use_engineered_features=False,
            fusion_hidden=32,
            dropout=0.1,
        )
    )

    param_groups = get_parameter_groups_with_llrd(
        model=model,
        base_lr=base_lr,
        weight_decay=weight_decay,
        layer_decay=layer_decay,
    )

    embedding_lr = _get_group_lr(param_groups, "embedding")
    expected_embedding_lr = base_lr * (layer_decay ** num_layers)
    assert math.isclose(embedding_lr, expected_embedding_lr, rel_tol=1e-6), (
        f"embedding lr mismatch: got {embedding_lr}, expected {expected_embedding_lr}"
    )

    for layer_idx in range(num_layers):
        layer_lr = _get_group_lr(param_groups, f"layer_{layer_idx}")
        expected_lr = base_lr * (layer_decay ** (num_layers - layer_idx - 1))
        assert math.isclose(layer_lr, expected_lr, rel_tol=1e-6), (
            f"layer {layer_idx} lr mismatch: got {layer_lr}, expected {expected_lr}"
        )

    top_lr = _get_group_lr(param_groups, "top")
    assert math.isclose(top_lr, base_lr, rel_tol=1e-6), (
        f"top lr mismatch: got {top_lr}, expected {base_lr}"
    )


def test_llrd_decay_ratio_math() -> None:
    base_lr = 2e-5
    layer_decay = 0.9
    num_layers = 12

    embedding_lr = base_lr * (layer_decay ** num_layers)
    ratio = base_lr / embedding_lr
    expected_ratio = (1 / layer_decay) ** num_layers

    assert math.isclose(ratio, expected_ratio, rel_tol=1e-6), (
        f"ratio mismatch: got {ratio}, expected {expected_ratio}"
    )
