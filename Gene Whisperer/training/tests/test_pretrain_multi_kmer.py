"""Tests for pretrain_multi_kmer.py module."""

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pretrain_multi_kmer import (
    get_kmer_config,
    load_config,
    validate_vocabs_exist,
)


class TestLoadConfig:
    """Tests for load_config function."""

    def test_load_valid_config(self, tmp_path: Path):
        """Test loading a valid YAML config."""
        config_content = """
seed: 1337
mlm_kmer: 6
mlm_epochs: 100
pretrain_kmers: [3, 4, 5, 6]
kmer_pretrain_overrides:
  3:
    mlm_epochs: 50
    mlm_lr: 0.0003
"""
        config_file = tmp_path / "config.yaml"
        config_file.write_text(config_content)

        cfg = load_config(str(config_file))

        assert cfg["seed"] == 1337
        assert cfg["mlm_kmer"] == 6
        assert cfg["pretrain_kmers"] == [3, 4, 5, 6]
        assert cfg["kmer_pretrain_overrides"][3]["mlm_epochs"] == 50

    def test_load_missing_config(self):
        """Test error when config file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            load_config("/nonexistent/config.yaml")


class TestGetKmerConfig:
    """Tests for get_kmer_config function."""

    def test_vocab_size_calculation(self):
        """Test that vocab size is calculated correctly for each k-mer."""
        cfg = {"vocab_cache_dir": "../artifacts/vocabs"}

        # k=3: 4^3 + 3 = 64 + 3 = 67
        overrides_k3 = get_kmer_config(cfg, 3)
        assert overrides_k3["vocab_size"] == 67

        # k=4: 4^4 + 3 = 256 + 3 = 259
        overrides_k4 = get_kmer_config(cfg, 4)
        assert overrides_k4["vocab_size"] == 259

        # k=5: 4^5 + 3 = 1024 + 3 = 1027
        overrides_k5 = get_kmer_config(cfg, 5)
        assert overrides_k5["vocab_size"] == 1027

        # k=6: 4^6 + 3 = 4096 + 3 = 4099
        overrides_k6 = get_kmer_config(cfg, 6)
        assert overrides_k6["vocab_size"] == 4099

    def test_mlm_kmer_set(self):
        """Test that mlm_kmer is set correctly."""
        cfg = {"vocab_cache_dir": "../artifacts/vocabs"}

        for k in [3, 4, 5, 6]:
            overrides = get_kmer_config(cfg, k)
            assert overrides["mlm_kmer"] == k

    def test_vocab_path_generation(self):
        """Test that vocab paths are generated correctly."""
        cfg = {"vocab_cache_dir": "../artifacts/vocabs"}

        overrides_k3 = get_kmer_config(cfg, 3)
        assert "mlm_k3_vocab.json" in overrides_k3["mlm_vocab_path"]

        overrides_k6 = get_kmer_config(cfg, 6)
        assert "mlm_k6_vocab.json" in overrides_k6["mlm_vocab_path"]

    def test_encoder_path_generation(self):
        """Test that encoder paths are generated correctly."""
        cfg = {"vocab_cache_dir": "../artifacts/vocabs"}

        overrides = get_kmer_config(cfg, 4)
        assert "mlm_encoder_k4.pt" in overrides["mlm_encoder_path"]

    def test_kmer_overrides_applied(self):
        """Test that k-mer specific overrides are applied."""
        cfg = {
            "vocab_cache_dir": "../artifacts/vocabs",
            "mlm_epochs": 200,
            "mlm_lr": 0.0002,
            "kmer_pretrain_overrides": {
                3: {
                    "mlm_epochs": 100,
                    "mlm_lr": 0.0003,
                },
                6: {
                    "mlm_epochs": 200,
                    "mlm_lr": 0.0002,
                },
            },
        }

        # k=3 should have custom epochs and lr
        overrides_k3 = get_kmer_config(cfg, 3)
        assert overrides_k3["mlm_epochs"] == 100
        assert overrides_k3["mlm_lr"] == 0.0003

        # k=6 should have default epochs from override
        overrides_k6 = get_kmer_config(cfg, 6)
        assert overrides_k6["mlm_epochs"] == 200

        # k=4 (no override) should not have mlm_epochs in overrides
        overrides_k4 = get_kmer_config(cfg, 4)
        assert "mlm_epochs" not in overrides_k4

    def test_checkpoint_paths_from_config(self):
        """Test that checkpoint paths can be loaded from config."""
        cfg = {
            "vocab_cache_dir": "../artifacts/vocabs",
            "checkpoints": {
                "mlm_encoder_ckpt_by_k": {
                    3: "../artifacts/custom_k3.pt",
                    6: "../artifacts/custom_k6.pt",
                },
            },
        }

        overrides_k3 = get_kmer_config(cfg, 3)
        assert overrides_k3["mlm_encoder_path"] == "../artifacts/custom_k3.pt"

        # k=4 not in config, should use default
        overrides_k4 = get_kmer_config(cfg, 4)
        assert "mlm_encoder_k4.pt" in overrides_k4["mlm_encoder_path"]


class TestValidateVocabsExist:
    """Tests for validate_vocabs_exist function."""

    def test_all_vocabs_exist(self, tmp_path: Path):
        """Test when all vocabularies exist."""
        vocab_dir = tmp_path / "vocabs"
        vocab_dir.mkdir()

        # Create vocab files
        for k in [3, 4, 5, 6]:
            vocab_file = vocab_dir / f"mlm_k{k}_vocab.json"
            vocab_file.write_text(json.dumps({"k": k, "itos": []}))

        cfg = {"vocab_cache_dir": str(vocab_dir)}
        missing = validate_vocabs_exist(cfg, [3, 4, 5, 6])
        assert missing == []

    def test_some_vocabs_missing(self, tmp_path: Path):
        """Test when some vocabularies are missing."""
        vocab_dir = tmp_path / "vocabs"
        vocab_dir.mkdir()

        # Only create k=3 and k=4
        for k in [3, 4]:
            vocab_file = vocab_dir / f"mlm_k{k}_vocab.json"
            vocab_file.write_text(json.dumps({"k": k, "itos": []}))

        cfg = {"vocab_cache_dir": str(vocab_dir)}
        missing = validate_vocabs_exist(cfg, [3, 4, 5, 6])
        assert missing == [5, 6]

    def test_all_vocabs_missing(self, tmp_path: Path):
        """Test when all vocabularies are missing."""
        vocab_dir = tmp_path / "vocabs"
        vocab_dir.mkdir()

        cfg = {"vocab_cache_dir": str(vocab_dir)}
        missing = validate_vocabs_exist(cfg, [3, 4, 5, 6])
        assert missing == [3, 4, 5, 6]


class TestPretrainIntegration:
    """Integration tests for multi-kmer pretraining."""

    @pytest.fixture
    def mock_config(self, tmp_path: Path):
        """Create a mock config file."""
        config_content = """
seed: 1337
mlm_kmer: 6
mlm_epochs: 2
mlm_batch_size: 4
mlm_lr: 0.001
vocab_cache_dir: {vocab_dir}
pretrain_kmers: [3, 4]
kmer_pretrain_overrides:
  3:
    mlm_epochs: 1
    mlm_lr: 0.002
  4:
    mlm_epochs: 1
    mlm_lr: 0.0015
checkpoints:
  mlm_encoder_ckpt_by_k:
    3: {artifacts}/mlm_encoder_k3.pt
    4: {artifacts}/mlm_encoder_k4.pt
""".format(
            vocab_dir=str(tmp_path / "vocabs"),
            artifacts=str(tmp_path / "artifacts"),
        )

        config_file = tmp_path / "config.yaml"
        config_file.write_text(config_content)

        # Create directories
        (tmp_path / "vocabs").mkdir()
        (tmp_path / "artifacts").mkdir()

        return str(config_file)

    def test_config_loading_and_override(self, mock_config: str):
        """Test that config loads correctly with overrides."""
        cfg = load_config(mock_config)

        # Base config
        assert cfg["mlm_kmer"] == 6
        assert cfg["pretrain_kmers"] == [3, 4]

        # Get k=3 config with overrides
        overrides_k3 = get_kmer_config(cfg, 3)
        assert overrides_k3["mlm_kmer"] == 3
        assert overrides_k3["vocab_size"] == 67
        assert overrides_k3["mlm_epochs"] == 1
        assert overrides_k3["mlm_lr"] == 0.002

        # Get k=4 config with overrides
        overrides_k4 = get_kmer_config(cfg, 4)
        assert overrides_k4["mlm_kmer"] == 4
        assert overrides_k4["vocab_size"] == 259
        assert overrides_k4["mlm_epochs"] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
