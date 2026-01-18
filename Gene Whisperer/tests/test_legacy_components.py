"""Tests for legacy model components kept for backward compatibility."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import torch

# Add training directory to path for local imports.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model_legacy import (  # noqa: E402
    CausalConv1d,
    GatedAttentionFusion,
    GeneWhispererStage1Legacy,
    MultiScaleCNN,
    MultiScaleTCN,
    PostCNNTransformerAdapter,
    TCNBlock,
    TCNEncoder,
)


def _assert_deprecated(factory) -> None:
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always", DeprecationWarning)
        factory()
    assert any(
        isinstance(item.message, DeprecationWarning) for item in recorded
    ), "Expected a DeprecationWarning from legacy component"


def test_legacy_components_emit_deprecation_warning() -> None:
    _assert_deprecated(lambda: CausalConv1d(in_channels=4, out_channels=8, kernel_size=3))
    _assert_deprecated(lambda: TCNBlock(in_channels=4, out_channels=8, kernel_size=3, dropout=0.0))
    _assert_deprecated(
        lambda: TCNEncoder(in_channels=4, hidden_channels=8, num_levels=1, kernel_size=3, dropout=0.0)
    )
    _assert_deprecated(
        lambda: MultiScaleCNN(in_channels=4, out_channels_per_scale=2, kernel_sizes=(3,), dropout=0.0)
    )
    _assert_deprecated(
        lambda: MultiScaleTCN(
            in_channels=4,
            multiscale_out_per_branch=2,
            multiscale_kernels=(3,),
            tcn_hidden=8,
            tcn_levels=1,
            tcn_kernel=3,
            dropout=0.0,
        )
    )
    _assert_deprecated(
        lambda: PostCNNTransformerAdapter(
            input_dim=8,
            transformer_dim=8,
            num_layers=0,
            num_heads=2,
            ff_dim=16,
            dropout=0.0,
        )
    )
    _assert_deprecated(
        lambda: GatedAttentionFusion(
            seq_dim=8,
            eng_dim=6,
            hidden_dim=8,
            num_heads=2,
            dropout=0.0,
        )
    )


def test_legacy_stage1_init_warns_and_runs_forward() -> None:
    _assert_deprecated(
        lambda: GeneWhispererStage1Legacy(
            vocab_size=16,
            kmer=3,
            embedding_dim=8,
            num_layers=1,
            num_heads=2,
            ff_dim=16,
            dropout=0.0,
            pad_token_id=0,
            engineered_dim=0,
            use_engineered_features=False,
            use_attention_pool=False,
            use_tcn=False,
            post_cnn_transformer_layers=0,
            max_seq_len=16,
            drop_path_rate=0.0,
        )
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        model = GeneWhispererStage1Legacy(
            vocab_size=16,
            kmer=3,
            embedding_dim=8,
            num_layers=1,
            num_heads=2,
            ff_dim=16,
            dropout=0.0,
            pad_token_id=0,
            engineered_dim=0,
            use_engineered_features=False,
            use_attention_pool=False,
            use_tcn=False,
            post_cnn_transformer_layers=0,
            max_seq_len=16,
            drop_path_rate=0.0,
        )
    tokens = torch.zeros((2, 8), dtype=torch.long)
    logits = model(tokens)
    assert logits.shape == (2, 1)
