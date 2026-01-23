"""Tests for ensemble_infer.verify_checkpoint_compatibility."""
from __future__ import annotations

import sys
from pathlib import Path

import torch
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TRAINING_DIR = PROJECT_ROOT / "training"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(TRAINING_DIR) not in sys.path:
    sys.path.insert(0, str(TRAINING_DIR))

from dataset import KmerVocabulary
from ensemble_infer import build_model, verify_checkpoint_compatibility


def test_verify_checkpoint_compatibility(tmp_path: Path) -> None:
    vocab = KmerVocabulary.build_from_sequences(["ACGTACGT"], k=3)
    vocab_path = tmp_path / "k3_vocab.json"
    vocab.save(vocab_path)

    cfg = {
        "embedding_dim": 16,
        "transformer_layers": 1,
        "transformer_heads": 4,
        "transformer_ff_dim": 32,
        "transformer_dropout": 0.05,
        "engineered_dim": 0,
        "stage1_use_engineered_features": False,
        "max_bp_len": 12,
        "vocab_cache_dir": str(tmp_path),
        "simplified_model": {
            "pooling_type": "attention",
            "classifier_hidden": 16,
            "classifier_dropout": 0.1,
        },
    }
    config_path = tmp_path / "config.yaml"
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle)

    model = build_model(cfg, vocab, torch.device("cpu"))
    checkpoint_path = tmp_path / "stage1_k3.pt"
    torch.save({"model_state": model.state_dict()}, checkpoint_path)

    assert verify_checkpoint_compatibility(
        str(checkpoint_path),
        config_path=str(config_path),
    )
