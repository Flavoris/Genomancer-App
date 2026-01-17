"""
Tests for the variant comparison experiment runner.
"""
import csv
import json
import tempfile
from dataclasses import asdict
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

import pytest
import torch
import torch.nn as nn

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "training"))

from run_variant_comparison import (
    VariantResult,
    count_parameters,
    create_model_for_variant,
    save_results_csv,
    save_results_json,
    print_summary_table,
    print_statistical_comparison,
    VARIANTS,
)


class TestVariantResult:
    """Tests for VariantResult dataclass."""

    def test_default_values(self):
        """Test default values for VariantResult."""
        result = VariantResult(variant="test")
        assert result.variant == "test"
        assert result.best_val_acc == 0.0
        assert result.best_val_mcc == 0.0
        assert result.best_val_auc == 0.0
        assert result.best_epoch == 0
        assert result.param_count == 0
        assert result.training_time_seconds == 0.0
        assert result.training_curves == {}

    def test_with_values(self):
        """Test VariantResult with custom values."""
        result = VariantResult(
            variant="combined",
            best_val_acc=0.95,
            best_val_mcc=0.85,
            best_val_auc=0.92,
            best_epoch=42,
            param_count=1_000_000,
            training_time_seconds=3600.0,
        )
        assert result.variant == "combined"
        assert result.best_val_acc == 0.95
        assert result.best_epoch == 42

    def test_asdict(self):
        """Test that VariantResult can be converted to dict."""
        result = VariantResult(variant="test", best_val_acc=0.9)
        d = asdict(result)
        assert d["variant"] == "test"
        assert d["best_val_acc"] == 0.9


class TestCountParameters:
    """Tests for count_parameters function."""

    def test_simple_linear(self):
        """Test counting parameters in a simple linear layer."""
        model = nn.Linear(10, 5)
        count = count_parameters(model)
        # 10*5 weights + 5 bias = 55
        assert count == 55

    def test_frozen_parameters(self):
        """Test that frozen parameters are not counted."""
        model = nn.Linear(10, 5)
        for param in model.parameters():
            param.requires_grad = False
        count = count_parameters(model)
        assert count == 0

    def test_mixed_frozen(self):
        """Test model with mix of frozen and trainable params."""
        model = nn.Sequential(
            nn.Linear(10, 5),  # 55 params
            nn.Linear(5, 3),   # 18 params
        )
        # Freeze first layer
        for param in model[0].parameters():
            param.requires_grad = False
        count = count_parameters(model)
        assert count == 18


class TestVariants:
    """Tests for variant constants."""

    def test_all_variants_present(self):
        """Test that all expected variants are defined."""
        assert "transformer_only" in VARIANTS
        assert "features_only" in VARIANTS
        assert "combined" in VARIANTS
        assert "original" in VARIANTS
        assert len(VARIANTS) == 4


class TestSaveResults:
    """Tests for results saving functions."""

    def test_save_results_csv(self, tmp_path):
        """Test saving results to CSV."""
        results = [
            VariantResult(
                variant="transformer_only",
                best_val_acc=0.90,
                best_val_mcc=0.80,
                best_val_auc=0.88,
                param_count=1000,
                training_time_seconds=100.0,
                best_epoch=10,
                best_threshold=0.45,
            ),
            VariantResult(
                variant="combined",
                best_val_acc=0.92,
                best_val_mcc=0.82,
                best_val_auc=0.90,
                param_count=2000,
                training_time_seconds=200.0,
                best_epoch=15,
                best_threshold=0.55,
            ),
        ]
        csv_path = tmp_path / "results.csv"
        save_results_csv(results, csv_path)

        assert csv_path.exists()
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 2
            assert rows[0]["Variant"] == "transformer_only"
            assert rows[0]["Val Acc"] == "0.9000"
            assert rows[1]["Variant"] == "combined"

    def test_save_results_json(self, tmp_path):
        """Test saving results to JSON."""
        results = [
            VariantResult(
                variant="features_only",
                best_val_acc=0.85,
                training_curves={"val_acc": [0.80, 0.82, 0.85]},
            ),
        ]
        json_path = tmp_path / "results.json"
        save_results_json(results, json_path)

        assert json_path.exists()
        with open(json_path) as f:
            data = json.load(f)
            assert "timestamp" in data
            assert "variants" in data
            assert len(data["variants"]) == 1
            assert data["variants"][0]["variant"] == "features_only"
            assert data["variants"][0]["training_curves"]["val_acc"] == [0.80, 0.82, 0.85]


class TestPrintFunctions:
    """Tests for print/display functions."""

    def test_print_summary_table(self, capsys):
        """Test that summary table prints without error."""
        results = [
            VariantResult(variant="test", best_val_acc=0.9, param_count=1000),
        ]
        print_summary_table(results)
        captured = capsys.readouterr()
        assert "VARIANT COMPARISON" in captured.out
        assert "test" in captured.out
        assert "0.9000" in captured.out

    def test_print_statistical_comparison(self, capsys):
        """Test that statistical comparison prints without error."""
        results = [
            VariantResult(variant="transformer_only", best_val_acc=0.88),
            VariantResult(variant="features_only", best_val_acc=0.75),
            VariantResult(variant="combined", best_val_acc=0.92, param_count=2000),
            VariantResult(variant="original", best_val_acc=0.90, param_count=10000),
        ]
        print_statistical_comparison(results)
        captured = capsys.readouterr()
        assert "STATISTICAL ANALYSIS" in captured.out
        # Combined > transformer by > 1%, so engineered features help
        assert "Engineered features help" in captured.out or "combined" in captured.out


class TestCreateModelForVariant:
    """Tests for model creation."""

    @pytest.fixture
    def small_config(self):
        """Minimal config for testing."""
        return {
            "vocab_size": 67,
            "kmer": 3,
            "embedding_dim": 32,
            "transformer_layers": 1,
            "transformer_heads": 2,
            "transformer_ff_dim": 64,
            "transformer_dropout": 0.1,
            "max_bp_len": 20,
            "engineered_dim": 16,
            "engineered_mlp_hidden": 32,
            "engineered_mlp_output": 16,
            "pad_token_id": 66,
            "drop_path_rate": 0.0,
            "use_relative_position_bias": False,
            "use_glu_ffn": False,
            # Original model specific
            "use_tcn": False,
            "tcn_hidden": 32,
            "tcn_levels": 2,
            "tcn_kernel": 3,
            "multiscale_channels": 16,
            "multiscale_kernels": [3, 5],
            "use_lstm": False,
            "lstm_hidden": 32,
            "post_cnn_transformer_layers": 1,
            "fusion_hidden": 32,
            "use_attention_pool": True,
        }

    def test_create_transformer_only(self, small_config):
        """Test creating transformer_only variant."""
        device = torch.device("cpu")
        model = create_model_for_variant("transformer_only", small_config, device)
        assert model is not None
        params = count_parameters(model)
        assert params > 0

    def test_create_features_only(self, small_config):
        """Test creating features_only variant."""
        device = torch.device("cpu")
        model = create_model_for_variant("features_only", small_config, device)
        assert model is not None
        params = count_parameters(model)
        assert params > 0
        # Features only should have far fewer params
        transformer_model = create_model_for_variant("transformer_only", small_config, device)
        assert params < count_parameters(transformer_model)

    def test_create_combined(self, small_config):
        """Test creating combined variant."""
        device = torch.device("cpu")
        model = create_model_for_variant("combined", small_config, device)
        assert model is not None
        params = count_parameters(model)
        # Combined should have more params than transformer_only
        transformer_model = create_model_for_variant("transformer_only", small_config, device)
        assert params > count_parameters(transformer_model)

    def test_create_original(self, small_config):
        """Test creating original variant."""
        device = torch.device("cpu")
        model = create_model_for_variant("original", small_config, device)
        assert model is not None
        params = count_parameters(model)
        assert params > 0


class TestStatisticalAnalysis:
    """Tests for statistical comparison logic."""

    def test_engineered_features_help(self, capsys):
        """Test detection when engineered features help."""
        results = [
            VariantResult(variant="transformer_only", best_val_acc=0.85),
            VariantResult(variant="combined", best_val_acc=0.92),
        ]
        print_statistical_comparison(results)
        captured = capsys.readouterr()
        assert "Engineered features help" in captured.out

    def test_transformer_helps(self, capsys):
        """Test detection when transformer helps."""
        results = [
            VariantResult(variant="features_only", best_val_acc=0.75),
            VariantResult(variant="combined", best_val_acc=0.90),
        ]
        print_statistical_comparison(results)
        captured = capsys.readouterr()
        assert "Transformer helps" in captured.out

    def test_synergy_detection(self, capsys):
        """Test synergy detection when combined > max of individual."""
        results = [
            VariantResult(variant="transformer_only", best_val_acc=0.85),
            VariantResult(variant="features_only", best_val_acc=0.80),
            VariantResult(variant="combined", best_val_acc=0.92),
        ]
        print_statistical_comparison(results)
        captured = capsys.readouterr()
        # Synergy: 0.92 > 0.85 + 0.01 threshold
        assert "Synergy detected" in captured.out

    def test_best_variant_identified(self, capsys):
        """Test that best variant is identified."""
        results = [
            VariantResult(variant="transformer_only", best_val_acc=0.85),
            VariantResult(variant="combined", best_val_acc=0.92),
            VariantResult(variant="original", best_val_acc=0.90),
        ]
        print_statistical_comparison(results)
        captured = capsys.readouterr()
        assert "Best variant: combined" in captured.out


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
