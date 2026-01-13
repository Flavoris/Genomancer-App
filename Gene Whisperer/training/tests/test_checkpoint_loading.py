"""Tests for checkpoint loading with robust error handling."""

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).parent.parent))

from model import _load_checkpoint_file


class TestLoadCheckpointFile:
    """Tests for the _load_checkpoint_file helper function."""

    def test_load_valid_pytorch_checkpoint(self, tmp_path: Path) -> None:
        """Valid PyTorch checkpoint should load successfully."""
        checkpoint_path = tmp_path / "valid.pt"
        state_dict = {"layer.weight": torch.randn(10, 10), "layer.bias": torch.randn(10)}
        torch.save(state_dict, checkpoint_path)

        loaded = _load_checkpoint_file(checkpoint_path)

        assert "layer.weight" in loaded
        assert "layer.bias" in loaded
        assert loaded["layer.weight"].shape == (10, 10)

    def test_load_nested_checkpoint(self, tmp_path: Path) -> None:
        """Checkpoint with nested model_state_dict should load."""
        checkpoint_path = tmp_path / "nested.pt"
        state_dict = {
            "model_state_dict": {"encoder.weight": torch.randn(5, 5)},
            "optimizer_state_dict": {},
            "epoch": 10,
        }
        torch.save(state_dict, checkpoint_path)

        loaded = _load_checkpoint_file(checkpoint_path)

        assert "model_state_dict" in loaded
        assert loaded["epoch"] == 10

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        """Missing file should raise FileNotFoundError."""
        missing_path = tmp_path / "nonexistent.pt"

        with pytest.raises(FileNotFoundError, match="Checkpoint file not found"):
            _load_checkpoint_file(missing_path)

    def test_corrupted_file_raises_value_error(self, tmp_path: Path) -> None:
        """Corrupted checkpoint should raise ValueError with helpful message."""
        corrupt_path = tmp_path / "corrupt.pt"
        corrupt_path.write_bytes(b"invalid pickle data that starts with v")

        with pytest.raises(ValueError, match="corrupted|not a valid"):
            _load_checkpoint_file(corrupt_path)

    def test_text_file_raises_value_error(self, tmp_path: Path) -> None:
        """Text file should be detected and raise ValueError."""
        text_path = tmp_path / "text.pt"
        text_path.write_text('{"this": "is json"}', encoding="utf-8")

        with pytest.raises(ValueError, match="text|JSON"):
            _load_checkpoint_file(text_path)

    def test_empty_file_raises_value_error(self, tmp_path: Path) -> None:
        """Empty file should raise ValueError."""
        empty_path = tmp_path / "empty.pt"
        empty_path.write_bytes(b"")

        with pytest.raises(ValueError, match="too small"):
            _load_checkpoint_file(empty_path)

    def test_very_small_file_raises_value_error(self, tmp_path: Path) -> None:
        """File too small to be valid should raise ValueError."""
        tiny_path = tmp_path / "tiny.pt"
        tiny_path.write_bytes(b"x")

        with pytest.raises(ValueError, match="too small"):
            _load_checkpoint_file(tiny_path)

    def test_safetensors_without_library_raises_import_error(
        self, tmp_path: Path
    ) -> None:
        """Safetensors file without library should raise ImportError or ValueError."""
        safetensors_path = tmp_path / "model.safetensors"
        safetensors_path.write_bytes(b"safetensors format data")

        # Will raise ImportError if safetensors not installed, or ValueError if parsing fails
        with pytest.raises((ImportError, ValueError)):
            _load_checkpoint_file(safetensors_path)


class TestLoadCheckpointFileIntegration:
    """Integration tests using actual model checkpoints."""

    def test_save_and_load_roundtrip(self, tmp_path: Path) -> None:
        """Save and load should preserve tensor values."""
        checkpoint_path = tmp_path / "roundtrip.pt"
        original = {
            "weight": torch.tensor([[1.0, 2.0], [3.0, 4.0]]),
            "config": {"hidden_dim": 256},
        }
        torch.save(original, checkpoint_path)

        loaded = _load_checkpoint_file(checkpoint_path)

        assert torch.allclose(loaded["weight"], original["weight"])
        assert loaded["config"]["hidden_dim"] == 256
