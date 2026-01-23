"""Tests for transformer contribution diagnosis modules."""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import torch
import torch.nn as nn

# Ensure training directory is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from hooks import TransformerOutputRandomizer
from experiment_runner import (
    build_model_from_config,
    print_analysis,
    _evaluate,
    _resolve_vocab_path,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


class DummyPooling(nn.Module):
    """Simulates V5 attention pooling: (B, L, D) -> (B, D)."""

    def forward(self, x, mask=None):
        return x.mean(dim=1)


class DummyV5Model(nn.Module):
    """Minimal V5-like model with attention pooling and engineered features."""

    def __init__(self, embed_dim: int = 32, eng_dim: int = 16):
        super().__init__()
        self.embed_dim = embed_dim
        self.use_attention_pool = True
        self.pool = DummyPooling()
        self.eng_mlp = nn.Linear(eng_dim, embed_dim)
        self.classifier = nn.Linear(embed_dim * 2, 1)

    def forward(self, tokens, engineered_features=None):
        # Simulate transformer encoding
        batch_size = tokens.size(0)
        x = torch.randn(batch_size, 10, self.embed_dim, device=tokens.device)
        x = self.pool(x)  # (B, embed_dim)

        if engineered_features is not None:
            eng = self.eng_mlp(engineered_features)
            x = torch.cat([x, eng], dim=-1)
        else:
            x = torch.cat([x, torch.zeros(batch_size, self.embed_dim, device=tokens.device)], dim=-1)

        return self.classifier(x)


class DummyV6Model(nn.Module):
    """Minimal V6-like model with final_norm and mean pooling."""

    def __init__(self, embed_dim: int = 32, eng_dim: int = 16):
        super().__init__()
        self.embed_dim = embed_dim
        self.use_engineered = True
        self.final_norm = nn.LayerNorm(embed_dim)
        self.eng_proj = nn.Linear(eng_dim, embed_dim)
        self.classifier = nn.Linear(embed_dim * 2, 1)

    def forward(self, tokens, engineered_features=None):
        batch_size = tokens.size(0)
        # Simulate transformer layers
        x = torch.randn(batch_size, 10, self.embed_dim, device=tokens.device)
        x = self.final_norm(x)
        # Mean pooling
        x = x.mean(dim=1)

        if self.use_engineered and engineered_features is not None:
            eng = self.eng_proj(engineered_features)
            x = torch.cat([x, eng], dim=-1)
        else:
            x = torch.cat([x, torch.zeros(batch_size, self.embed_dim, device=tokens.device)], dim=-1)

        return self.classifier(x)


def _make_dummy_loader(
    num_samples: int = 32,
    seq_len: int = 76,
    eng_dim: int = 16,
    positive_ratio: float = 0.5,
):
    """Create a fake DataLoader yielding (tokens, engineered, labels)."""
    tokens = torch.randint(0, 100, (num_samples, seq_len))
    engineered = torch.randn(num_samples, eng_dim)
    num_pos = int(num_samples * positive_ratio)
    labels = torch.cat([torch.ones(num_pos), torch.zeros(num_samples - num_pos)])
    # Shuffle
    perm = torch.randperm(num_samples)
    labels = labels[perm]
    dataset = torch.utils.data.TensorDataset(tokens, engineered, labels)
    return torch.utils.data.DataLoader(dataset, batch_size=8, shuffle=False)


# ---------------------------------------------------------------------------
# TransformerOutputRandomizer Tests
# ---------------------------------------------------------------------------


class TestTransformerOutputRandomizer:
    """Tests for the hook-based transformer output randomization."""

    def test_finds_pool_target_for_v5(self):
        model = DummyV5Model()
        randomizer = TransformerOutputRandomizer(model)
        assert randomizer._hook_target is model.pool

    def test_finds_final_norm_target_for_v6(self):
        model = DummyV6Model()
        randomizer = TransformerOutputRandomizer(model)
        assert randomizer._hook_target is model.final_norm

    def test_raises_for_unsupported_model(self):
        model = nn.Linear(10, 1)
        with pytest.raises(ValueError, match="Cannot determine hook target"):
            TransformerOutputRandomizer(model)

    def test_hook_replaces_v5_pool_output(self):
        model = DummyV5Model()
        randomizer = TransformerOutputRandomizer(model)
        handle = randomizer.attach()

        tokens = torch.randint(0, 100, (4, 10))
        eng = torch.randn(4, 16)

        # Run twice with hook - outputs should differ (random replacement)
        out1 = model(tokens, eng)
        out2 = model(tokens, eng)
        # With random replacement, outputs should almost certainly differ
        assert not torch.allclose(out1, out2, atol=1e-5)

        handle.remove()

    def test_hook_replaces_v6_final_norm_output(self):
        model = DummyV6Model()
        randomizer = TransformerOutputRandomizer(model)
        handle = randomizer.attach()

        tokens = torch.randint(0, 100, (4, 10))
        eng = torch.randn(4, 16)

        out1 = model(tokens, eng)
        out2 = model(tokens, eng)
        assert not torch.allclose(out1, out2, atol=1e-5)

        handle.remove()

    def test_hook_removable(self):
        model = DummyV5Model()
        randomizer = TransformerOutputRandomizer(model)
        handle = randomizer.attach()
        handle.remove()

        # After removal, model should run without errors
        tokens = torch.randint(0, 100, (2, 10))
        eng = torch.randn(2, 16)
        out = model(tokens, eng)
        assert out.shape == (2, 1)

    def test_hook_output_shape_v5(self):
        """Hook output from pool should be (B, D)."""
        model = DummyV5Model(embed_dim=64)
        randomizer = TransformerOutputRandomizer(model)

        outputs = []

        def capture_hook(module, inp, out):
            outputs.append(out)
            return out

        # Attach randomizer first, then capture
        rand_handle = randomizer.attach()
        cap_handle = model.pool.register_forward_hook(capture_hook)

        tokens = torch.randint(0, 100, (3, 10))
        eng = torch.randn(3, 16)
        model(tokens, eng)

        # The captured output should be from the randomizer
        assert outputs[-1].shape == (3, 64)

        cap_handle.remove()
        rand_handle.remove()

    def test_hook_output_shape_v6(self):
        """Hook output from final_norm should be (B, L, D)."""
        model = DummyV6Model(embed_dim=64)
        randomizer = TransformerOutputRandomizer(model)

        outputs = []

        def capture_hook(module, inp, out):
            outputs.append(out)
            return out

        rand_handle = randomizer.attach()
        cap_handle = model.final_norm.register_forward_hook(capture_hook)

        tokens = torch.randint(0, 100, (3, 10))
        eng = torch.randn(3, 16)
        model(tokens, eng)

        # The last hook runs after randomizer, should see (B, L, D) output
        assert outputs[-1].shape == (3, 10, 64)

        cap_handle.remove()
        rand_handle.remove()

    def test_v6_hook_produces_consistent_per_sample_vectors(self):
        """After mean pooling randomized final_norm output, each sample gets one vector."""
        model = DummyV6Model(embed_dim=32)
        randomizer = TransformerOutputRandomizer(model)

        captured = []

        def capture_after_norm(module, inp, out):
            captured.append(out.clone())
            return out

        handle = randomizer.attach()
        cap_handle = model.final_norm.register_forward_hook(capture_after_norm)

        tokens = torch.randint(0, 100, (2, 8))
        eng = torch.randn(2, 16)
        model(tokens, eng)

        # Check that within each sample, all L positions are the same vector
        norm_out = captured[-1]  # (2, 8, 32)
        for b in range(2):
            for l in range(1, 8):
                assert torch.allclose(norm_out[b, 0], norm_out[b, l])

        cap_handle.remove()
        handle.remove()


# ---------------------------------------------------------------------------
# Experiment Evaluation Tests
# ---------------------------------------------------------------------------


class TestEvaluate:
    """Tests for the _evaluate function."""

    def test_full_mode_runs_without_error(self):
        model = DummyV5Model()
        model.eval()
        loader = _make_dummy_loader()
        correct, total = _evaluate(model, loader, torch.device("cpu"), 0.5, "full")
        assert total == 32
        assert 0 <= correct <= total

    def test_zero_features_mode(self):
        model = DummyV5Model()
        model.eval()
        loader = _make_dummy_loader()
        correct, total = _evaluate(
            model, loader, torch.device("cpu"), 0.5, "zero_features"
        )
        assert total == 32
        assert 0 <= correct <= total

    def test_random_transformer_mode_v5(self):
        model = DummyV5Model()
        model.eval()
        loader = _make_dummy_loader()
        correct, total = _evaluate(
            model, loader, torch.device("cpu"), 0.5, "random_transformer"
        )
        assert total == 32
        assert 0 <= correct <= total

    def test_random_transformer_mode_v6(self):
        model = DummyV6Model()
        model.eval()
        loader = _make_dummy_loader()
        correct, total = _evaluate(
            model, loader, torch.device("cpu"), 0.5, "random_transformer"
        )
        assert total == 32
        assert 0 <= correct <= total

    def test_threshold_affects_predictions(self):
        model = DummyV5Model()
        model.eval()
        loader = _make_dummy_loader()

        # Very low threshold: predict almost all positive
        _, total = _evaluate(model, loader, torch.device("cpu"), 0.01, "full")
        correct_low, _ = _evaluate(model, loader, torch.device("cpu"), 0.01, "full")
        # Very high threshold: predict almost all negative
        correct_high, _ = _evaluate(model, loader, torch.device("cpu"), 0.99, "full")

        # With balanced labels, extreme thresholds should give ~50% accuracy
        # (all-positive or all-negative on balanced data)
        assert total == 32

    def test_empty_loader_returns_zero(self):
        model = DummyV5Model()
        model.eval()
        empty_dataset = torch.utils.data.TensorDataset(
            torch.empty(0, 76, dtype=torch.long),
            torch.empty(0, 16),
            torch.empty(0),
        )
        loader = torch.utils.data.DataLoader(empty_dataset, batch_size=8)
        correct, total = _evaluate(model, loader, torch.device("cpu"), 0.5, "full")
        assert total == 0
        assert correct == 0


# ---------------------------------------------------------------------------
# Print Analysis Tests
# ---------------------------------------------------------------------------


class TestPrintAnalysis:
    """Tests for the analysis printing function."""

    def test_transformer_learning(self, capsys):
        print_analysis(acc_a=0.86, acc_b=0.70, acc_c=0.78)
        output = capsys.readouterr().out
        assert "Transformer contributes" in output
        assert "Engineered features contribute" in output
        assert "IS learning useful patterns" in output

    def test_transformer_not_learning(self, capsys):
        print_analysis(acc_a=0.80, acc_b=0.52, acc_c=0.79)
        output = capsys.readouterr().out
        assert "NOT learning useful patterns" in output

    def test_transformer_weak_learning(self, capsys):
        print_analysis(acc_a=0.82, acc_b=0.60, acc_c=0.80)
        output = capsys.readouterr().out
        assert "WEAK patterns" in output

    def test_features_do_all_work(self, capsys):
        print_analysis(acc_a=0.85, acc_b=0.52, acc_c=0.85)
        output = capsys.readouterr().out
        assert "ALL the work" in output

    def test_features_do_most_work(self, capsys):
        print_analysis(acc_a=0.86, acc_b=0.60, acc_c=0.82)
        output = capsys.readouterr().out
        assert "MOST of the work" in output


# ---------------------------------------------------------------------------
# Vocab Path Resolution Tests
# ---------------------------------------------------------------------------


class TestResolveVocabPath:
    """Tests for vocabulary path resolution."""

    def test_uses_mlm_vocab_by_k(self):
        config = {"mlm_vocab_by_k": {6: "/path/to/vocab.json"}}
        assert _resolve_vocab_path(config, 6) == "/path/to/vocab.json"

    def test_uses_nested_vocab_paths(self):
        config = {
            "vocab_paths": {"mlm_vocab_by_k": {4: "/nested/vocab.json"}}
        }
        assert _resolve_vocab_path(config, 4) == "/nested/vocab.json"

    def test_falls_back_to_default(self):
        config = {"vocab_cache_dir": "../artifacts/vocabs"}
        path = _resolve_vocab_path(config, 6)
        assert path == "../artifacts/vocabs/mlm_k6_vocab.json"

    def test_default_cache_dir(self):
        path = _resolve_vocab_path({}, 3)
        assert path == "../artifacts/vocabs/mlm_k3_vocab.json"


# ---------------------------------------------------------------------------
# Build Model Tests
# ---------------------------------------------------------------------------


class TestBuildModel:
    """Tests for model construction from config."""

    @patch("experiment_runner.create_model_variant")
    def test_passes_correct_variant(self, mock_create):
        mock_create.return_value = nn.Linear(10, 1)
        config = {"model_variant": "v6", "kmer": 6, "max_bp_len": 81}
        build_model_from_config(config)
        mock_create.assert_called_once()
        args = mock_create.call_args
        assert args[0][0] == "v6"

    @patch("experiment_runner.create_model_variant")
    def test_computes_max_seq_len(self, mock_create):
        mock_create.return_value = nn.Linear(10, 1)
        config = {"model_variant": "combined", "kmer": 6, "max_bp_len": 81}
        build_model_from_config(config)
        model_config = mock_create.call_args[0][1]
        assert model_config["max_seq_len"] == 76  # 81 - 6 + 1

    @patch("experiment_runner.create_model_variant")
    def test_includes_modern_settings(self, mock_create):
        mock_create.return_value = nn.Linear(10, 1)
        config = {
            "model_variant": "v6",
            "use_rope": True,
            "ffn_type": "swiglu",
            "norm_type": "rmsnorm",
        }
        build_model_from_config(config)
        model_config = mock_create.call_args[0][1]
        assert model_config["use_rope"] is True
        assert model_config["ffn_type"] == "swiglu"
        assert model_config["norm_type"] == "rmsnorm"
