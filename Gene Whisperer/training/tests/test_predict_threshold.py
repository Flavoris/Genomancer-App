"""Tests for loading inference thresholds from Stage 1 checkpoints."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TRAINING_DIR = PROJECT_ROOT / "training"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(TRAINING_DIR) not in sys.path:
    sys.path.insert(0, str(TRAINING_DIR))

from dataset import KmerVocabulary
from gene_whisperer.predict import _build_model, load_stage1_ensemble


def _write_config(tmp_path: Path) -> tuple[dict, Path]:
    cfg = {
        "embedding_dim": 12,
        "transformer_layers": 2,
        "transformer_heads": 3,
        "transformer_ff_dim": 24,
        "transformer_dropout": 0.1,
        "use_tcn": True,
        "tcn_hidden": 8,
        "tcn_levels": 1,
        "tcn_kernel": 3,
        "multiscale_channels": 4,
        "multiscale_kernels": [3],
        "lstm_hidden": 8,
        "post_cnn_transformer_layers": 1,
        "engineered_mlp_hidden": 8,
        "engineered_mlp_output": 4,
        "fusion_hidden": 8,
        "engineered_dim": 288,
        "stage1_use_engineered_features": True,
        "use_attention_pool": True,
    }
    config_path = tmp_path / "config.yaml"
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle)
    return cfg, config_path


def _create_checkpoint(
    tmp_path: Path,
    best_threshold: float,
) -> tuple[Path, Path, Path]:
    cfg, config_path = _write_config(tmp_path)
    vocab = KmerVocabulary.build_from_sequences(["ACGTACGT"], k=3)
    vocab_path = tmp_path / "k3_vocab.json"
    vocab.save(vocab_path)

    model = _build_model(cfg, vocab, device=torch.device("cpu"))
    checkpoint = {
        "model_state": model.state_dict(),
        "best_threshold": best_threshold,
    }
    checkpoint_path = tmp_path / "stage1_k3.pt"
    torch.save(checkpoint, checkpoint_path)
    return config_path, checkpoint_path, vocab_path


def test_load_stage1_ensemble_uses_checkpoint_threshold(tmp_path: Path) -> None:
    config_path, checkpoint_path, vocab_path = _create_checkpoint(tmp_path, best_threshold=0.37)

    predictor = load_stage1_ensemble(
        config_path=config_path,
        checkpoints_and_vocabs={
            3: {"checkpoint": str(checkpoint_path), "vocab": str(vocab_path)}
        },
        device=torch.device("cpu"),
    )

    assert predictor.threshold == pytest.approx(0.37)


def test_load_stage1_ensemble_allows_threshold_override(tmp_path: Path) -> None:
    config_path, checkpoint_path, vocab_path = _create_checkpoint(tmp_path, best_threshold=0.37)

    predictor = load_stage1_ensemble(
        config_path=config_path,
        checkpoints_and_vocabs={
            3: {"checkpoint": str(checkpoint_path), "vocab": str(vocab_path)}
        },
        device=torch.device("cpu"),
        threshold=0.5,
    )

    assert predictor.threshold == 0.5
