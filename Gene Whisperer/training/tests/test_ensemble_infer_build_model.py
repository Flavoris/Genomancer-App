"""Tests for ensemble_infer.build_model architecture alignment."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TRAINING_DIR = PROJECT_ROOT / "training"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(TRAINING_DIR) not in sys.path:
    sys.path.insert(0, str(TRAINING_DIR))

from bpe_tokenizer import DNABPETokenizer
from ensemble_infer import build_model
from model import RMSNorm, SwiGLU


def _make_vocab() -> DNABPETokenizer:
    vocab = DNABPETokenizer(vocab_size=20)
    vocab.train(["ACGTACGTAC", "GGGGCCCCAA", "ATATATAT"])
    return vocab


def test_build_model_applies_stage1_arch_params() -> None:
    vocab = _make_vocab()
    cfg = {
        "embedding_dim": 32,
        "transformer_layers": 2,
        "transformer_heads": 4,
        "transformer_ff_dim": 64,
        "transformer_dropout": 0.1,
        "engineered_dim": 0,
        "stage1_use_engineered_features": False,
        "max_token_len": 20,
        "drop_path_rate": 0.05,
        "stage1_drop_path_rate": 0.2,
        "use_relative_position_bias": False,
        "use_glu_ffn": False,
        "glu_activation": "gelu",
        "use_rope": True,
        "rope_base": 12345.0,
        "ffn_type": "swiglu",
        "norm_type": "rmsnorm",
        "ffn_mult": 2.0,
        "simplified_model": {
            "pooling_type": "attention",
            "classifier_hidden": 24,
            "classifier_dropout": 0.3,
        },
    }

    model = build_model(cfg, vocab, torch.device("cpu"))
    encoder = model._full_encoder
    layer = encoder.layers[-1]

    assert model.use_attention_pool is True
    assert model.pool is not None
    assert encoder.use_rope is True
    assert layer.rotary is not None
    assert layer.rotary.base == pytest.approx(12345.0)
    assert isinstance(layer.ffn, SwiGLU)
    assert isinstance(layer.norm1, RMSNorm)
    assert layer.drop_path_attn.drop_prob == pytest.approx(0.2)
    assert model.classifier[0].out_features == 24
    assert model.classifier[3].p == pytest.approx(0.3)


def test_build_model_uses_simplified_defaults() -> None:
    vocab = _make_vocab()
    cfg = {
        "embedding_dim": 16,
        "transformer_layers": 1,
        "transformer_heads": 4,
        "transformer_ff_dim": 32,
        "transformer_dropout": 0.05,
        "engineered_dim": 0,
        "stage1_use_engineered_features": False,
        "max_token_len": 12,
    }

    model = build_model(cfg, vocab, torch.device("cpu"))

    assert model.use_attention_pool is True
    assert model.pool is not None
    assert model.classifier[0].out_features == 256
    assert model.classifier[3].p == pytest.approx(0.15)
