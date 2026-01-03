"""Coverage for Stage 1 engineered feature dimension wiring."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model import GeneWhispererStage1


def _base_stage1_kwargs() -> dict:
    return {
        "vocab_size": 50,
        "kmer": 3,
        "embedding_dim": 16,
        "num_layers": 1,
        "num_heads": 4,
        "ff_dim": 32,
        "dropout": 0.0,
        "pad_token_id": 0,
        "use_tcn": False,
        "use_attention_pool": False,
        "post_cnn_transformer_layers": 1,
    }


@pytest.mark.parametrize(
    "engineered_dim, mlp_hidden, mlp_output, fusion_hidden",
    [
        (128, 256, 128, 256),
        (288, 512, 256, 384),
    ],
)
def test_stage1_engineered_dims_forward(
    engineered_dim: int,
    mlp_hidden: int,
    mlp_output: int,
    fusion_hidden: int,
) -> None:
    model = GeneWhispererStage1(
        **_base_stage1_kwargs(),
        engineered_dim=engineered_dim,
        engineered_mlp_hidden=mlp_hidden,
        engineered_mlp_output=mlp_output,
        fusion_hidden=fusion_hidden,
        use_engineered_features=True,
    )
    model.eval()

    tokens = torch.randint(1, 50, (2, 12))
    engineered = torch.randn(2, engineered_dim)

    with torch.no_grad():
        output = model(tokens, engineered)

    assert output.shape == (2, 1)
    assert model.engineered_mlp.output_dim == mlp_output
    assert model.fusion.output_dim == fusion_hidden


def test_stage1_engineered_defaults_match_config() -> None:
    model = GeneWhispererStage1(
        **_base_stage1_kwargs(),
        engineered_dim=288,
        use_engineered_features=True,
    )

    first_linear = model.engineered_mlp.mlp[0]
    assert first_linear.in_features == 288
    assert first_linear.out_features == 512
    assert model.engineered_mlp.output_dim == 256
    assert model.fusion.output_dim == 384


def test_engineered_mlp_improvements_default_noops() -> None:
    model = GeneWhispererStage1(
        **_base_stage1_kwargs(),
        engineered_dim=128,
        use_engineered_features=True,
    )
    mlp = model.engineered_mlp

    assert mlp.input_norm is not None
    assert mlp.input_dropout is not None
    assert mlp.gate_proj is not None
    assert mlp.residual_proj is not None
    assert torch.allclose(mlp.pre_norm_scale, torch.zeros_like(mlp.pre_norm_scale))
    assert torch.allclose(mlp.gate_scale, torch.zeros_like(mlp.gate_scale))
    assert torch.allclose(mlp.residual_scale, torch.zeros_like(mlp.residual_scale))


def test_engineered_mlp_improvements_can_disable() -> None:
    model = GeneWhispererStage1(
        **_base_stage1_kwargs(),
        engineered_dim=128,
        use_engineered_features=True,
        engineered_mlp_pre_norm=False,
        engineered_mlp_input_dropout=0.0,
        engineered_mlp_use_gated=False,
        engineered_mlp_use_residual=False,
    )
    mlp = model.engineered_mlp

    assert mlp.input_norm is None
    assert mlp.input_dropout is None
    assert mlp.gate_proj is None
    assert mlp.residual_proj is None
