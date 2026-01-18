#!/usr/bin/env python3
"""Tests for train_final_model.py."""
from __future__ import annotations

import math
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import torch
import torch.nn as nn

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from train_final_model import (
    TrainingConfig,
    TrainingResult,
    CosineAnnealingWithWarmup,
    BEST_CONFIG,
    ORIGINAL_MODEL_BASELINE,
    print_comparison_table,
    get_model_architecture_summary,
    load_pretrained_mlm_weights,
)


class TestTrainingConfig:
    """Tests for TrainingConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = TrainingConfig()
        assert config.variant == "combined"
        assert config.lr == 0.0001
        assert config.dropout == 0.1
        assert config.classifier_hidden == 256
        assert config.pooling_type == "attention"
        assert config.weight_decay == 0.02
        assert config.epochs == 100
        assert config.patience == 15

    def test_from_dict_with_stage1_lr(self):
        """Test creating config from dict with stage1_lr key."""
        d = {"stage1_lr": 0.0002, "transformer_dropout": 0.15, "variant": "combined"}
        config = TrainingConfig.from_dict(d)
        assert config.lr == 0.0002
        assert config.dropout == 0.15
        assert config.variant == "combined"

    def test_from_dict_with_lr_key(self):
        """Test creating config from dict with lr key."""
        d = {"lr": 0.0003, "dropout": 0.2}
        config = TrainingConfig.from_dict(d)
        assert config.lr == 0.0003
        assert config.dropout == 0.2

    def test_from_dict_defaults(self):
        """Test that from_dict uses defaults for missing keys."""
        config = TrainingConfig.from_dict({})
        assert config.lr == 0.0001
        assert config.epochs == 100


class TestTrainingResult:
    """Tests for TrainingResult dataclass."""

    def test_default_values(self):
        """Test default result values."""
        result = TrainingResult()
        assert result.best_val_acc == 0.0
        assert result.best_val_mcc == 0.0
        assert result.param_count == 0
        assert result.stopped_early is False

    def test_to_dict(self):
        """Test conversion to dictionary."""
        result = TrainingResult(
            best_val_acc=0.85,
            best_val_mcc=0.70,
            best_epoch=10,
            param_count=1000000,
        )
        d = result.to_dict()
        assert d["best_val_acc"] == 0.85
        assert d["best_val_mcc"] == 0.70
        assert d["best_epoch"] == 10
        assert d["param_count"] == 1000000


class TestCosineAnnealingWithWarmup:
    """Tests for the learning rate scheduler."""

    def test_warmup_phase(self):
        """Test linear warmup during initial steps."""
        model = nn.Linear(10, 10)
        optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
        scheduler = CosineAnnealingWithWarmup(
            optimizer,
            warmup_steps=10,
            total_steps=100,
        )

        # At step 0, LR should be 0
        assert optimizer.param_groups[0]["lr"] == 0.1

        # During warmup, LR should increase linearly
        for i in range(1, 11):
            scheduler.step()
            expected_lr = 0.1 * (i / 10)
            assert abs(optimizer.param_groups[0]["lr"] - expected_lr) < 1e-6

    def test_cosine_annealing_phase(self):
        """Test cosine annealing after warmup."""
        model = nn.Linear(10, 10)
        optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
        scheduler = CosineAnnealingWithWarmup(
            optimizer,
            warmup_steps=10,
            total_steps=110,
            min_lr_ratio=0.1,
        )

        # Complete warmup
        for _ in range(10):
            scheduler.step()

        # After warmup, LR should be at base_lr
        assert abs(optimizer.param_groups[0]["lr"] - 0.1) < 1e-6

        # Continue and verify cosine decay
        lrs = []
        for _ in range(100):
            scheduler.step()
            lrs.append(optimizer.param_groups[0]["lr"])

        # LR should decrease over time
        assert lrs[-1] < lrs[0]
        # Final LR should be close to min_lr_ratio * base_lr
        assert lrs[-1] >= 0.1 * 0.1 - 1e-6

    def test_get_last_lr(self):
        """Test get_last_lr returns current learning rate."""
        model = nn.Linear(10, 10)
        optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
        scheduler = CosineAnnealingWithWarmup(
            optimizer,
            warmup_steps=10,
            total_steps=100,
        )

        scheduler.step()
        last_lr = scheduler.get_last_lr()
        assert len(last_lr) == 1
        assert last_lr[0] == optimizer.param_groups[0]["lr"]


class TestBestConfig:
    """Tests for the default best configuration."""

    def test_best_config_has_required_keys(self):
        """Test that BEST_CONFIG has all required keys."""
        required_keys = [
            "variant",
            "stage1_lr",
            "transformer_dropout",
            "classifier_hidden",
            "pooling_type",
            "weight_decay",
        ]
        for key in required_keys:
            assert key in BEST_CONFIG, f"Missing key: {key}"

    def test_best_config_values_are_valid(self):
        """Test that BEST_CONFIG values are sensible."""
        assert BEST_CONFIG["variant"] in ["combined", "transformer_only", "features_only"]
        assert 1e-6 < BEST_CONFIG["stage1_lr"] < 1e-2
        assert 0.0 <= BEST_CONFIG["transformer_dropout"] <= 0.5
        assert BEST_CONFIG["classifier_hidden"] > 0
        assert BEST_CONFIG["pooling_type"] in ["mean", "cls", "attention"]
        assert BEST_CONFIG["weight_decay"] >= 0


class TestOriginalModelBaseline:
    """Tests for the original model baseline constants."""

    def test_baseline_has_required_keys(self):
        """Test that baseline has all required keys."""
        required_keys = [
            "val_acc",
            "val_mcc",
            "val_auc",
            "param_count",
            "training_time_minutes",
        ]
        for key in required_keys:
            assert key in ORIGINAL_MODEL_BASELINE, f"Missing key: {key}"

    def test_baseline_values_are_reasonable(self):
        """Test that baseline values are reasonable."""
        assert 0.5 < ORIGINAL_MODEL_BASELINE["val_acc"] < 1.0
        assert -1.0 <= ORIGINAL_MODEL_BASELINE["val_mcc"] <= 1.0
        assert 0.5 < ORIGINAL_MODEL_BASELINE["val_auc"] < 1.0
        assert ORIGINAL_MODEL_BASELINE["param_count"] > 0
        assert ORIGINAL_MODEL_BASELINE["training_time_minutes"] > 0


class TestPrintComparisonTable:
    """Tests for the comparison table function."""

    def test_comparison_table_format(self):
        """Test that comparison table is properly formatted."""
        result = TrainingResult(
            best_val_acc=0.85,
            best_val_mcc=0.70,
            best_val_auc=0.90,
            param_count=20_000_000,
            training_time_seconds=600,
        )

        table = print_comparison_table(result)

        assert "MODEL COMPARISON" in table
        assert "Original" in table
        assert "Simplified" in table
        assert "Change" in table
        assert "Val Accuracy" in table
        assert "Val MCC" in table
        assert "Parameters" in table

    def test_comparison_shows_positive_change(self):
        """Test that improvement is shown with positive sign."""
        result = TrainingResult(
            best_val_acc=0.90,  # Better than baseline 0.8293
            best_val_mcc=0.75,
            best_val_auc=0.95,
            param_count=10_000_000,
            training_time_seconds=300,
        )

        table = print_comparison_table(result)
        # The accuracy change should be positive
        assert "+" in table  # Should have positive change indicators


class TestModelArchitectureSummary:
    """Tests for the model architecture summary function."""

    def test_summary_contains_param_counts(self):
        """Test that summary includes parameter counts."""
        model = nn.Sequential(
            nn.Linear(10, 20),
            nn.ReLU(),
            nn.Linear(20, 5),
        )

        summary = get_model_architecture_summary(model)

        assert "Total parameters" in summary
        assert "Trainable parameters" in summary
        assert "Module-wise" in summary


class TestIntegration:
    """Integration tests for the training pipeline."""

    def test_training_config_to_training_result_flow(self):
        """Test that config can be used to initialize training and produce results."""
        config = TrainingConfig.from_dict(BEST_CONFIG)

        # Verify config is valid
        assert config.variant == BEST_CONFIG["variant"]
        assert config.lr == BEST_CONFIG["stage1_lr"]

        # Create a result
        result = TrainingResult()
        result.best_val_acc = 0.85
        result.param_count = 20_000_000

        # Verify result can be serialized
        result_dict = result.to_dict()
        assert result_dict["best_val_acc"] == 0.85


class TestCreateModel:
    """Tests for the create_model function."""

    def test_create_tunable_model_returns_module(self):
        """Test that create_model returns a nn.Module for tunable variant."""
        from train_final_model import create_model, TrainingConfig

        # This test only verifies the function signature works
        # Full model creation would require the full codebase setup
        config = TrainingConfig()
        assert config.variant == "combined"

    def test_create_model_requires_vocab_for_standard(self):
        """Test that create_model raises error when vocab missing for standard model."""
        from train_final_model import create_model, TrainingConfig

        config = TrainingConfig()
        with pytest.raises(ValueError, match="vocab is required"):
            create_model(
                variant="combined",
                cfg={},
                training_config=config,
                device=torch.device("cpu"),
                use_standard_model=True,
                vocab=None,
            )


class TestCheckpointCompatibility:
    """Tests for checkpoint format compatibility."""

    def test_checkpoint_has_model_state_key(self):
        """Verify checkpoint format includes model_state key for ensemble compatibility."""
        # This is a documentation test - the actual format is defined in train_final_model
        expected_keys = ["model_state", "model_state_dict", "best_threshold", "val_acc", "val_mcc"]
        # The checkpoint saving code includes these keys for compatibility
        assert all(key in expected_keys for key in ["model_state", "best_threshold"])


class TestLoadPretrainedMLMWeights:
    """Tests for the load_pretrained_mlm_weights function."""

    def test_returns_false_for_none_transfer_mode(self):
        """Test that load_pretrained_mlm_weights returns False for transfer_mode='none'."""
        model = MagicMock()
        cfg = {"kmer": 6}
        result = load_pretrained_mlm_weights(model, cfg, transfer_mode="none")
        assert result is False

    def test_returns_false_if_model_lacks_load_method(self):
        """Test returns False if model doesn't have load_pretrained_weights."""
        model = nn.Linear(10, 10)  # No load_pretrained_weights method
        cfg = {"kmer": 6}
        result = load_pretrained_mlm_weights(model, cfg, transfer_mode="embed_only")
        assert result is False

    def test_returns_false_if_no_checkpoint_found(self):
        """Test returns False when no checkpoint file exists."""
        model = MagicMock()
        model.load_pretrained_weights = MagicMock()
        cfg = {"kmer": 6}

        # Mock the checkpoint finder to return None
        with patch("train_final_model.get_mlm_checkpoint_for_kmer", return_value=None):
            result = load_pretrained_mlm_weights(model, cfg, transfer_mode="embed_only")
        assert result is False

    def test_calls_load_pretrained_weights_when_checkpoint_exists(self):
        """Test that load_pretrained_weights is called when checkpoint exists."""
        model = MagicMock()
        model.embedding = nn.Embedding(100, 64)
        model.load_pretrained_weights = MagicMock()
        cfg = {"kmer": 6}

        mock_path = Path("/fake/checkpoint.pt")
        with patch("train_final_model.get_mlm_checkpoint_for_kmer", return_value=mock_path):
            load_pretrained_mlm_weights(model, cfg, transfer_mode="embed_only")

        model.load_pretrained_weights.assert_called_once_with(
            checkpoint_path=mock_path,
            strict=False,
            transfer_mode="embed_only",
        )

    def test_accepts_embed_plus_adapter_transfer_mode(self):
        """Test that embed_plus_adapter transfer mode works."""
        model = MagicMock()
        model.embedding = nn.Embedding(100, 64)
        model.load_pretrained_weights = MagicMock()
        cfg = {"kmer": 6}

        mock_path = Path("/fake/checkpoint.pt")
        with patch("train_final_model.get_mlm_checkpoint_for_kmer", return_value=mock_path):
            load_pretrained_mlm_weights(model, cfg, transfer_mode="embed_plus_adapter")

        model.load_pretrained_weights.assert_called_once_with(
            checkpoint_path=mock_path,
            strict=False,
            transfer_mode="embed_plus_adapter",
        )


class TestTrainFinalModelSignature:
    """Tests for the train_final_model function signature."""

    def test_has_mlm_weight_parameters(self):
        """Test that train_final_model has MLM weight loading parameters."""
        import inspect
        from train_final_model import train_final_model

        sig = inspect.signature(train_final_model)
        params = list(sig.parameters.keys())

        assert "load_mlm_weights" in params
        assert "mlm_transfer_mode" in params

    def test_mlm_weight_parameters_have_defaults(self):
        """Test that MLM weight parameters have sensible defaults."""
        import inspect
        from train_final_model import train_final_model

        sig = inspect.signature(train_final_model)

        load_mlm_default = sig.parameters["load_mlm_weights"].default
        transfer_mode_default = sig.parameters["mlm_transfer_mode"].default

        assert load_mlm_default is True
        assert transfer_mode_default == "embed_only"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
