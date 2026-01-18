"""Smoke tests for post-CNN transformer disabling.

Ensures `post_cnn_transformer_layers=0` behaves as an identity pass-through and
does not break the Stage 1 forward path.
"""
from __future__ import annotations

import sys
from pathlib import Path

import torch

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent))

from model import GeneWhispererStage1Legacy


def test_stage1_forward_with_zero_post_cnn_transformer_layers():
    model = GeneWhispererStage1Legacy(
        vocab_size=67,
        kmer=3,
        embedding_dim=32,
        num_layers=1,
        num_heads=2,
        ff_dim=64,
        dropout=0.0,
        engineered_dim=0,
        use_engineered_features=False,
        use_attention_pool=False,
        use_tcn=True,
        tcn_hidden=32,
        tcn_levels=1,
        tcn_kernel=3,
        multiscale_channels=8,
        multiscale_kernels=(3, 5),
        post_cnn_transformer_layers=0,
    )

    tokens = torch.randint(0, 67, (2, 64))
    out = model(tokens, None)
    assert out.shape == (2, 1)
