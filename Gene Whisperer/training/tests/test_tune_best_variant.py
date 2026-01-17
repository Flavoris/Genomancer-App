#!/usr/bin/env python3
"""Tests for tune_best_variant.py hyperparameter tuning script."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from tune_best_variant import (
    TuningResult,
    HYPERPARAM_SPACE,
    generate_search_configs,
    get_best_variant,
    load_comparison_results,
    print_best_config_yaml,
    save_results_csv,
    save_results_json,
)


class TestTuningResult:
    """Tests for TuningResult dataclass."""

    def test_default_values(self):
        """Test TuningResult has correct default values."""
        result = TuningResult(
            config_id=1,
            lr=1e-5,
            dropout=0.1,
            classifier_hidden=256,
            pooling_type="mean",
            weight_decay=0.01,
        )
        assert result.config_id == 1
        assert result.lr == 1e-5
        assert result.dropout == 0.1
        assert result.classifier_hidden == 256
        assert result.pooling_type == "mean"
        assert result.weight_decay == 0.01
        assert result.best_val_acc == 0.0
        assert result.best_val_mcc == 0.0
        assert result.best_epoch == 0
        assert result.training_time_seconds == 0.0

    def test_to_dict(self):
        """Test TuningResult converts to dict correctly."""
        result = TuningResult(
            config_id=5,
            lr=2e-5,
            dropout=0.15,
            classifier_hidden=512,
            pooling_type="attention",
            weight_decay=0.02,
            best_val_acc=0.85,
            best_val_mcc=0.70,
            best_epoch=15,
            training_time_seconds=120.5,
        )
        d = result.to_dict()
        assert d["config_id"] == 5
        assert d["lr"] == 2e-5
        assert d["dropout"] == 0.15
        assert d["classifier_hidden"] == 512
        assert d["pooling_type"] == "attention"
        assert d["weight_decay"] == 0.02
        assert d["best_val_acc"] == 0.85
        assert d["best_val_mcc"] == 0.70
        assert d["best_epoch"] == 15
        assert d["training_time_seconds"] == 120.5


class TestHyperparameterSpace:
    """Tests for hyperparameter search space."""

    def test_space_contains_required_params(self):
        """Test search space has all required parameters."""
        assert "lr" in HYPERPARAM_SPACE
        assert "dropout" in HYPERPARAM_SPACE
        assert "classifier_hidden" in HYPERPARAM_SPACE
        assert "pooling_type" in HYPERPARAM_SPACE
        assert "weight_decay" in HYPERPARAM_SPACE

    def test_lr_values(self):
        """Test learning rate options."""
        expected = [1e-5, 2e-5, 5e-5, 1e-4]
        assert HYPERPARAM_SPACE["lr"] == expected

    def test_dropout_values(self):
        """Test dropout options."""
        expected = [0.1, 0.15, 0.2]
        assert HYPERPARAM_SPACE["dropout"] == expected

    def test_classifier_hidden_values(self):
        """Test classifier hidden dim options."""
        expected = [128, 256, 512]
        assert HYPERPARAM_SPACE["classifier_hidden"] == expected

    def test_pooling_type_values(self):
        """Test pooling type options."""
        expected = ["mean", "cls", "attention"]
        assert HYPERPARAM_SPACE["pooling_type"] == expected

    def test_weight_decay_values(self):
        """Test weight decay options."""
        expected = [0.01, 0.02, 0.05]
        assert HYPERPARAM_SPACE["weight_decay"] == expected


class TestGenerateSearchConfigs:
    """Tests for configuration generation."""

    def test_random_search_count(self):
        """Test random search generates correct number of configs."""
        configs = generate_search_configs("random", num_trials=10, seed=42)
        assert len(configs) == 10

    def test_random_search_has_all_keys(self):
        """Test random search configs have all required keys."""
        configs = generate_search_configs("random", num_trials=5, seed=42)
        required_keys = {"lr", "dropout", "classifier_hidden", "pooling_type", "weight_decay"}
        for cfg in configs:
            assert set(cfg.keys()) == required_keys

    def test_random_search_values_in_space(self):
        """Test random search values come from search space."""
        configs = generate_search_configs("random", num_trials=20, seed=42)
        for cfg in configs:
            assert cfg["lr"] in HYPERPARAM_SPACE["lr"]
            assert cfg["dropout"] in HYPERPARAM_SPACE["dropout"]
            assert cfg["classifier_hidden"] in HYPERPARAM_SPACE["classifier_hidden"]
            assert cfg["pooling_type"] in HYPERPARAM_SPACE["pooling_type"]
            assert cfg["weight_decay"] in HYPERPARAM_SPACE["weight_decay"]

    def test_random_search_reproducible(self):
        """Test random search is reproducible with same seed."""
        configs1 = generate_search_configs("random", num_trials=5, seed=123)
        configs2 = generate_search_configs("random", num_trials=5, seed=123)
        assert configs1 == configs2

    def test_random_search_different_seeds(self):
        """Test random search differs with different seeds."""
        configs1 = generate_search_configs("random", num_trials=20, seed=1)
        configs2 = generate_search_configs("random", num_trials=20, seed=2)
        # Very unlikely to be identical
        assert configs1 != configs2

    def test_grid_search_count(self):
        """Test grid search generates all combinations."""
        configs = generate_search_configs("grid", num_trials=0, seed=42)
        expected_count = (
            len(HYPERPARAM_SPACE["lr"]) *
            len(HYPERPARAM_SPACE["dropout"]) *
            len(HYPERPARAM_SPACE["classifier_hidden"]) *
            len(HYPERPARAM_SPACE["pooling_type"]) *
            len(HYPERPARAM_SPACE["weight_decay"])
        )
        # 4 * 3 * 3 * 3 * 3 = 324
        assert expected_count == 324
        assert len(configs) == expected_count

    def test_grid_search_unique_configs(self):
        """Test grid search produces unique configurations."""
        configs = generate_search_configs("grid", num_trials=0, seed=42)
        # Convert to tuples for hashing
        config_tuples = [tuple(sorted(c.items())) for c in configs]
        assert len(config_tuples) == len(set(config_tuples))


class TestGetBestVariant:
    """Tests for best variant selection."""

    def test_gets_best_by_accuracy(self):
        """Test selects variant with highest accuracy."""
        results = {
            "variants": [
                {"variant": "transformer_only", "best_val_acc": 0.80},
                {"variant": "combined", "best_val_acc": 0.85},
                {"variant": "features_only", "best_val_acc": 0.78},
            ]
        }
        best = get_best_variant(results)
        assert best == "combined"

    def test_empty_variants_raises(self):
        """Test raises error for empty variants list."""
        with pytest.raises(ValueError, match="No variant results found"):
            get_best_variant({"variants": []})

    def test_missing_variants_raises(self):
        """Test raises error for missing variants key."""
        with pytest.raises(ValueError, match="No variant results found"):
            get_best_variant({})


class TestLoadComparisonResults:
    """Tests for loading comparison results."""

    def test_loads_json_file(self):
        """Test loads JSON file correctly."""
        test_data = {
            "timestamp": "2024-01-01",
            "variants": [
                {"variant": "combined", "best_val_acc": 0.83}
            ]
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(test_data, f)
            f.flush()

            results = load_comparison_results(Path(f.name))
            assert results["timestamp"] == "2024-01-01"
            assert len(results["variants"]) == 1


class TestSaveResults:
    """Tests for saving results."""

    def test_save_results_csv(self):
        """Test saves results to CSV correctly."""
        results = [
            TuningResult(1, 1e-5, 0.1, 256, "mean", 0.01, 0.80, 0.65, 10, 100.0),
            TuningResult(2, 2e-5, 0.15, 512, "attention", 0.02, 0.82, 0.68, 12, 120.0),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "results.csv"
            save_results_csv(results, output_path)

            assert output_path.exists()
            with open(output_path) as f:
                content = f.read()
                assert "config_id" in content
                assert "lr" in content
                assert "best_val_acc" in content
                assert "0.8" in content or "0.80" in content

    def test_save_results_json(self):
        """Test saves results to JSON correctly."""
        results = [
            TuningResult(1, 1e-5, 0.1, 256, "mean", 0.01, 0.80, 0.65, 10, 100.0),
        ]
        best_result = results[0]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "results.json"
            save_results_json(results, best_result, "combined", output_path)

            assert output_path.exists()
            with open(output_path) as f:
                data = json.load(f)
                assert data["variant"] == "combined"
                assert "best_config" in data
                assert "all_results" in data
                assert data["best_config"]["lr"] == 1e-5


class TestPrintBestConfigYaml:
    """Tests for YAML config generation."""

    def test_generates_yaml_snippet(self):
        """Test generates valid YAML-like snippet."""
        result = TuningResult(
            config_id=1,
            lr=2e-5,
            dropout=0.15,
            classifier_hidden=256,
            pooling_type="attention",
            weight_decay=0.02,
            best_val_acc=0.8307,
            best_val_mcc=0.6644,
            best_epoch=16,
        )

        yaml_str = print_best_config_yaml(result, "combined")

        assert "stage1_lr: 2e-05" in yaml_str or "stage1_lr: 0.00002" in yaml_str
        assert "transformer_dropout: 0.15" in yaml_str
        assert "classifier_hidden: 256" in yaml_str
        assert 'pooling_type: "attention"' in yaml_str
        assert "weight_decay: 0.02" in yaml_str
        assert "Val Accuracy: 0.8307" in yaml_str
        assert "combined" in yaml_str

    def test_includes_comments(self):
        """Test YAML snippet includes helpful comments."""
        result = TuningResult(1, 1e-5, 0.1, 128, "mean", 0.01, 0.75, 0.55, 5, 50.0)
        yaml_str = print_best_config_yaml(result, "transformer_only")

        assert "Best Hyperparameters" in yaml_str
        assert "Variant:" in yaml_str
        assert "Val Accuracy:" in yaml_str
        assert "Val MCC:" in yaml_str


class TestCreateTunableModel:
    """Tests for tunable model creation (requires torch)."""

    @pytest.fixture
    def mock_cfg(self):
        """Provide a mock configuration dict."""
        return {
            "kmer": 6,
            "max_bp_len": 81,
            "vocab_size": 4099,
            "embedding_dim": 384,
            "transformer_layers": 2,  # Reduced for testing
            "transformer_heads": 8,
            "transformer_ff_dim": 768,
            "pad_token_id": 4098,
            "drop_path_rate": 0.1,
            "use_relative_position_bias": False,
            "use_glu_ffn": False,
            "engineered_dim": 288,
            "engineered_mlp_hidden": 256,
            "engineered_mlp_output": 128,
        }

    def test_creates_features_only_model(self, mock_cfg):
        """Test creates features-only model."""
        try:
            import torch
            from tune_best_variant import create_tunable_model

            hp_config = {
                "lr": 1e-5,
                "dropout": 0.1,
                "classifier_hidden": 256,
                "pooling_type": "mean",
                "weight_decay": 0.01,
            }

            device = torch.device("cpu")
            model = create_tunable_model("features_only", mock_cfg, hp_config, device)

            # Test forward pass
            batch_size = 2
            eng_features = torch.randn(batch_size, 288)
            dummy_tokens = torch.zeros(batch_size, 10, dtype=torch.long)

            output = model(dummy_tokens, engineered_features=eng_features)
            assert output.shape == (batch_size, 1)
        except ImportError:
            pytest.skip("PyTorch not available")

    def test_creates_transformer_only_model(self, mock_cfg):
        """Test creates transformer-only model with different pooling types."""
        try:
            import torch
            from tune_best_variant import create_tunable_model

            device = torch.device("cpu")

            for pooling_type in ["mean", "cls", "attention"]:
                hp_config = {
                    "lr": 1e-5,
                    "dropout": 0.1,
                    "classifier_hidden": 256,
                    "pooling_type": pooling_type,
                    "weight_decay": 0.01,
                }

                model = create_tunable_model("transformer_only", mock_cfg, hp_config, device)

                batch_size = 2
                seq_len = 76  # max_bp_len - kmer + 1
                tokens = torch.randint(0, 4099, (batch_size, seq_len))

                output = model(tokens)
                assert output.shape == (batch_size, 1), f"Failed for pooling_type={pooling_type}"
        except ImportError:
            pytest.skip("PyTorch not available")

    def test_creates_combined_model(self, mock_cfg):
        """Test creates combined model."""
        try:
            import torch
            from tune_best_variant import create_tunable_model

            hp_config = {
                "lr": 2e-5,
                "dropout": 0.15,
                "classifier_hidden": 512,
                "pooling_type": "attention",
                "weight_decay": 0.02,
            }

            device = torch.device("cpu")
            model = create_tunable_model("combined", mock_cfg, hp_config, device)

            batch_size = 2
            seq_len = 76
            tokens = torch.randint(0, 4099, (batch_size, seq_len))
            eng_features = torch.randn(batch_size, 288)

            output = model(tokens, engineered_features=eng_features)
            assert output.shape == (batch_size, 1)
        except ImportError:
            pytest.skip("PyTorch not available")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
