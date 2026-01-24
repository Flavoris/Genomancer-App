"""Tests for MLM quality diagnosis modules."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import torch

# Ensure training directory is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dataset import KmerVocabulary
from mlm_quality_metrics import (
    EncoderConfig,
    MLMMetrics,
    build_mlm_model,
    compute_mlm_predictions,
    compute_silhouette,
    extract_embeddings,
    interpret_silhouette,
    load_encoder_config,
    mask_tokens_for_eval,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def k3_vocab() -> KmerVocabulary:
    """Create a small k=3 vocabulary for testing."""
    import itertools
    bases = ["A", "C", "G", "T"]
    tokens = ["".join(kmer) for kmer in itertools.product(bases, repeat=3)]
    return KmerVocabulary(k=3, tokens=tokens)


@pytest.fixture
def k4_vocab() -> KmerVocabulary:
    """Create a k=4 vocabulary for testing."""
    import itertools
    bases = ["A", "C", "G", "T"]
    tokens = ["".join(kmer) for kmer in itertools.product(bases, repeat=4)]
    return KmerVocabulary(k=4, tokens=tokens)


@pytest.fixture
def small_encoder_cfg() -> EncoderConfig:
    """Small encoder config for fast testing."""
    return EncoderConfig(
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        dropout=0.0,
        max_seq_len=128,
    )


@pytest.fixture
def sample_sequences() -> list[str]:
    """Sample DNA sequences for testing."""
    return [
        "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGT",
        "TTTAAAGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC",
    ]


# ---------------------------------------------------------------------------
# EncoderConfig Tests
# ---------------------------------------------------------------------------


class TestEncoderConfig:
    """Tests for encoder configuration loading."""

    def test_default_config(self):
        cfg = EncoderConfig()
        assert cfg.embedding_dim == 384
        assert cfg.num_layers == 12
        assert cfg.num_heads == 12
        assert cfg.ff_dim == 1536
        assert cfg.dropout == 0.1

    def test_load_from_metadata(self, tmp_path):
        metadata = {
            "embedding_dim": 256,
            "num_layers": 8,
            "heads": 8,
            "config_snapshot": {
                "transformer_ff_dim": 1024,
                "transformer_dropout": 0.1,
            },
        }
        meta_path = tmp_path / "meta.json"
        meta_path.write_text(json.dumps(metadata))

        cfg = load_encoder_config(meta_path)
        assert cfg.embedding_dim == 256
        assert cfg.num_layers == 8
        assert cfg.num_heads == 8
        assert cfg.ff_dim == 1024

    def test_load_missing_file_uses_defaults(self, tmp_path):
        cfg = load_encoder_config(tmp_path / "nonexistent.json")
        assert cfg.embedding_dim == 384
        assert cfg.num_layers == 12

    def test_load_partial_metadata(self, tmp_path):
        metadata = {"embedding_dim": 128}
        meta_path = tmp_path / "meta.json"
        meta_path.write_text(json.dumps(metadata))

        cfg = load_encoder_config(meta_path)
        assert cfg.embedding_dim == 128
        assert cfg.num_layers == 12  # Falls back to default


# ---------------------------------------------------------------------------
# Build Model Tests
# ---------------------------------------------------------------------------


class TestBuildModel:
    """Tests for model building and checkpoint loading."""

    def test_build_model_no_checkpoint(self, k3_vocab, small_encoder_cfg, tmp_path):
        missing_path = tmp_path / "missing.pt"
        encoder, loaded = build_mlm_model(k3_vocab, small_encoder_cfg, missing_path)
        assert not loaded
        assert encoder.embedding_dim == 32
        assert encoder.vocab_size == len(k3_vocab.itos)

    def test_build_model_with_checkpoint(self, k3_vocab, small_encoder_cfg, tmp_path):
        # Create a real encoder and save its state dict
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )
        ckpt_path = tmp_path / "encoder.pt"
        torch.save(encoder.state_dict(), ckpt_path)

        loaded_encoder, success = build_mlm_model(k3_vocab, small_encoder_cfg, ckpt_path)
        assert success
        assert loaded_encoder.embedding_dim == 32

    def test_build_model_nested_checkpoint(self, k3_vocab, small_encoder_cfg, tmp_path):
        """Handle checkpoints wrapped in model_state_dict key."""
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )
        ckpt_path = tmp_path / "encoder.pt"
        torch.save({"model_state_dict": encoder.state_dict()}, ckpt_path)

        loaded_encoder, success = build_mlm_model(k3_vocab, small_encoder_cfg, ckpt_path)
        assert success

    def test_build_model_shape_mismatch(self, k3_vocab, tmp_path):
        """Incompatible checkpoint returns loaded=False."""
        # Save a checkpoint with different dimensions
        fake_state = {"embedding.weight": torch.randn(10, 999)}
        ckpt_path = tmp_path / "bad.pt"
        torch.save(fake_state, ckpt_path)

        cfg = EncoderConfig(embedding_dim=32, num_layers=2, num_heads=4, ff_dim=64)
        encoder, loaded = build_mlm_model(k3_vocab, cfg, ckpt_path)
        assert not loaded


# ---------------------------------------------------------------------------
# Masking Tests
# ---------------------------------------------------------------------------


class TestMasking:
    """Tests for evaluation masking."""

    def test_mask_tokens_basic(self, k3_vocab):
        token_ids = torch.randint(0, 60, (4, 20))
        masked, labels = mask_tokens_for_eval(token_ids, k3_vocab, mask_prob=0.15)

        # Labels should be -100 for unmasked positions
        assert (labels == -100).sum() > 0
        # Some positions should have real labels
        assert (labels != -100).sum() > 0
        # Masked positions should have [MASK] token
        mask_positions = labels != -100
        assert (masked[mask_positions] == k3_vocab.mask_id).all()

    def test_mask_tokens_never_masks_pad(self, k3_vocab):
        # Create input with padding at the end
        token_ids = torch.randint(0, 60, (2, 20))
        token_ids[:, 15:] = k3_vocab.pad_id

        masked, labels = mask_tokens_for_eval(token_ids, k3_vocab, mask_prob=0.3)

        # Pad positions should never be masked (labels=-100)
        assert (labels[:, 15:] == -100).all()

    def test_mask_tokens_mask_prob(self, k3_vocab):
        token_ids = torch.randint(0, 60, (10, 50))
        masked, labels = mask_tokens_for_eval(token_ids, k3_vocab, mask_prob=0.15)

        # Approximately 15% should be masked (with some variance)
        total_tokens = 10 * 50
        masked_count = (labels != -100).sum().item()
        mask_ratio = masked_count / total_tokens
        assert 0.05 < mask_ratio < 0.30  # Allow wide range for small samples

    def test_mask_tokens_empty_batch(self, k3_vocab):
        # All-pad input
        token_ids = torch.full((2, 10), k3_vocab.pad_id)
        masked, labels = mask_tokens_for_eval(token_ids, k3_vocab, mask_prob=0.15)
        assert (labels == -100).all()


# ---------------------------------------------------------------------------
# MLM Prediction Tests
# ---------------------------------------------------------------------------


class TestMLMPredictions:
    """Tests for MLM prediction metrics computation."""

    def test_compute_predictions_returns_metrics(self, k3_vocab, small_encoder_cfg):
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        token_ids = torch.randint(0, 60, (8, 20))
        device = torch.device("cpu")

        metrics = compute_mlm_predictions(
            encoder, k3_vocab, token_ids, device,
            batch_size=4, mask_prob=0.15,
        )

        assert isinstance(metrics, MLMMetrics)
        assert metrics.k == 3
        assert 0.0 <= metrics.top1_accuracy <= 1.0
        assert 0.0 <= metrics.top5_accuracy <= 1.0
        assert metrics.top5_accuracy >= metrics.top1_accuracy
        assert metrics.avg_entropy >= 0.0
        assert metrics.num_masked_tokens > 0

    def test_predictions_entropy_positive(self, k3_vocab):
        """Entropy should be positive for non-trivial predictions."""
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        token_ids = torch.randint(0, 60, (16, 30))
        device = torch.device("cpu")

        metrics = compute_mlm_predictions(
            encoder, k3_vocab, token_ids, device,
            batch_size=8, mask_prob=0.20,
        )
        assert metrics.avg_entropy > 0.0

    def test_top5_geq_top1(self, k3_vocab):
        """Top-5 accuracy should always be >= top-1 accuracy."""
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        token_ids = torch.randint(0, 60, (20, 25))
        device = torch.device("cpu")

        metrics = compute_mlm_predictions(
            encoder, k3_vocab, token_ids, device,
            batch_size=10, mask_prob=0.15,
        )
        assert metrics.top5_accuracy >= metrics.top1_accuracy


# ---------------------------------------------------------------------------
# Embedding Tests
# ---------------------------------------------------------------------------


class TestEmbeddings:
    """Tests for embedding extraction."""

    def test_extract_embeddings_shape(self, k3_vocab):
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        token_ids = torch.randint(0, 60, (10, 20))
        device = torch.device("cpu")

        embeddings = extract_embeddings(encoder, token_ids, device, batch_size=4)
        assert embeddings.shape == (10, 32)

    def test_extract_embeddings_deterministic(self, k3_vocab):
        """Same input should produce same embeddings."""
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        token_ids = torch.randint(0, 60, (5, 20))
        device = torch.device("cpu")

        emb1 = extract_embeddings(encoder, token_ids, device, batch_size=5)
        emb2 = extract_embeddings(encoder, token_ids, device, batch_size=5)
        np.testing.assert_allclose(emb1, emb2, atol=1e-5)

    def test_extract_embeddings_handles_padding(self, k3_vocab):
        """Padding tokens should not contribute to the mean."""
        from model import DNAEncoder
        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        # Create input with varying padding
        token_ids = torch.randint(0, 60, (4, 20))
        token_ids[0, 10:] = k3_vocab.pad_id
        token_ids[1, 5:] = k3_vocab.pad_id

        device = torch.device("cpu")
        embeddings = extract_embeddings(encoder, token_ids, device, batch_size=4)

        # All embeddings should be finite
        assert np.all(np.isfinite(embeddings))
        assert embeddings.shape == (4, 32)


# ---------------------------------------------------------------------------
# Silhouette Score Tests
# ---------------------------------------------------------------------------


class TestSilhouette:
    """Tests for silhouette score computation."""

    def test_silhouette_separated_clusters(self):
        """Well-separated clusters should have high silhouette."""
        embeddings = np.vstack([
            np.random.randn(50, 10) + 5,  # Cluster 1
            np.random.randn(50, 10) - 5,  # Cluster 2
        ])
        labels = np.array([0] * 50 + [1] * 50)
        score = compute_silhouette(embeddings, labels)
        assert score > 0.5

    def test_silhouette_mixed_clusters(self):
        """Overlapping clusters should have low silhouette."""
        embeddings = np.random.randn(100, 10)
        labels = np.array([0] * 50 + [1] * 50)
        score = compute_silhouette(embeddings, labels)
        assert score < 0.5

    def test_silhouette_single_class(self):
        """Single class should return NaN."""
        embeddings = np.random.randn(50, 10)
        labels = np.zeros(50)
        score = compute_silhouette(embeddings, labels)
        assert np.isnan(score)

    def test_silhouette_subsampling(self):
        """Large inputs should be subsampled."""
        embeddings = np.vstack([
            np.random.randn(2000, 10) + 3,
            np.random.randn(2000, 10) - 3,
        ])
        labels = np.array([0] * 2000 + [1] * 2000)
        score = compute_silhouette(embeddings, labels, max_samples=500)
        assert not np.isnan(score)


# ---------------------------------------------------------------------------
# Interpretation Tests
# ---------------------------------------------------------------------------


class TestInterpretation:
    """Tests for score interpretation helpers."""

    def test_interpret_silhouette_ranges(self):
        assert interpret_silhouette(0.6) == "strong"
        assert interpret_silhouette(0.4) == "good"
        assert interpret_silhouette(0.2) == "moderate"
        assert interpret_silhouette(0.1) == "weak"
        assert interpret_silhouette(-0.1) == "poor"
        assert interpret_silhouette(float("nan")) == "N/A"


# ---------------------------------------------------------------------------
# Integration Tests
# ---------------------------------------------------------------------------


class TestIntegration:
    """End-to-end integration tests."""

    def test_full_pipeline_k3(self, k3_vocab, sample_sequences, tmp_path):
        """Full pipeline: tokenize → mask → predict → embed → silhouette."""
        from model import DNAEncoder

        encoder = DNAEncoder(
            vocab_size=len(k3_vocab.itos),
            kmer=3,
            embedding_dim=32,
            num_layers=2,
            num_heads=4,
            ff_dim=64,
            dropout=0.0,
            pad_token_id=k3_vocab.pad_id,
            max_seq_len=128,
            drop_path_rate=0.0,
        )

        # Tokenize
        max_bp_len = 81
        tokens = torch.stack([
            k3_vocab.tokenize(seq, max_bp_len) for seq in sample_sequences
        ])
        assert tokens.shape == (4, 79)  # 81 - 3 + 1 = 79

        device = torch.device("cpu")

        # MLM predictions
        metrics = compute_mlm_predictions(
            encoder, k3_vocab, tokens, device,
            batch_size=2, mask_prob=0.15,
        )
        assert metrics.num_samples == 4
        assert metrics.num_masked_tokens > 0

        # Embeddings
        embeddings = extract_embeddings(encoder, tokens, device, batch_size=2)
        assert embeddings.shape == (4, 32)

        # Silhouette (with dummy labels)
        labels = np.array([1, 1, 0, 0])
        score = compute_silhouette(embeddings, labels)
        assert isinstance(score, float)

    def test_different_kmer_values(self, sample_sequences):
        """Ensure pipeline works for multiple k values."""
        import itertools
        from model import DNAEncoder

        for k in [3, 4, 5]:
            tokens_list = ["".join(kmer) for kmer in itertools.product("ACGT", repeat=k)]
            vocab = KmerVocabulary(k=k, tokens=tokens_list)

            encoder = DNAEncoder(
                vocab_size=len(vocab.itos),
                kmer=k,
                embedding_dim=32,
                num_layers=2,
                num_heads=4,
                ff_dim=64,
                dropout=0.0,
                pad_token_id=vocab.pad_id,
                max_seq_len=128,
                drop_path_rate=0.0,
            )

            max_bp_len = 81
            tokens = torch.stack([
                vocab.tokenize(seq, max_bp_len) for seq in sample_sequences[:2]
            ])
            expected_len = max_bp_len - k + 1
            assert tokens.shape == (2, expected_len)

            metrics = compute_mlm_predictions(
                encoder, vocab, tokens, torch.device("cpu"),
                batch_size=2, mask_prob=0.15,
            )
            assert metrics.k == k
            assert metrics.num_masked_tokens > 0


# ---------------------------------------------------------------------------
# CLI Tests
# ---------------------------------------------------------------------------


class TestCLI:
    """Tests for the CLI entry point."""

    def test_load_validation_data(self, tmp_path):
        """Test loading validation TSV data."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from diagnose_mlm_quality import load_validation_data

        val_file = tmp_path / "val.tsv"
        val_file.write_text(
            "sequence\tis_promoter\n"
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\t1\n"
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\t0\n"
        )

        sequences, labels = load_validation_data(val_file)
        assert len(sequences) == 2
        assert labels[0] == 1
        assert labels[1] == 0

    def test_load_validation_data_max_samples(self, tmp_path):
        """Test max_samples parameter."""
        from diagnose_mlm_quality import load_validation_data

        lines = ["sequence\tis_promoter\n"]
        for i in range(50):
            seq = "A" * 81
            lines.append(f"{seq}\t{i % 2}\n")

        val_file = tmp_path / "val.tsv"
        val_file.write_text("".join(lines))

        sequences, labels = load_validation_data(val_file, max_samples=10)
        assert len(sequences) == 10
        assert len(labels) == 10

    def test_load_validation_missing_file(self, tmp_path):
        """Missing file raises FileNotFoundError."""
        from diagnose_mlm_quality import load_validation_data
        with pytest.raises(FileNotFoundError):
            load_validation_data(tmp_path / "missing.tsv")

    def test_load_validation_bad_columns(self, tmp_path):
        """Wrong column names raises ValueError."""
        from diagnose_mlm_quality import load_validation_data
        val_file = tmp_path / "bad.tsv"
        val_file.write_text("col1\tcol2\nfoo\tbar\n")
        with pytest.raises(ValueError, match="Expected columns"):
            load_validation_data(val_file)

    def test_tokenize_sequences(self, k3_vocab, sample_sequences):
        """Test sequence tokenization helper."""
        from diagnose_mlm_quality import tokenize_sequences
        tokens = tokenize_sequences(sample_sequences, k3_vocab, max_bp_len=81)
        assert tokens.shape == (4, 79)
        assert tokens.dtype == torch.long
