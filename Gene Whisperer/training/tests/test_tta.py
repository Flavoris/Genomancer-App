"""Tests for Test-Time Augmentation (TTA) functionality."""
import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bpe_tokenizer import DNABPETokenizer
from dataset import reverse_complement
from tta import (
    get_reverse_complement_tokens,
    predict_with_tta,
    predict_batch_with_tta,
    TTAWrapper,
)


class MockModel(torch.nn.Module):
    """Mock model for testing TTA."""

    def __init__(self, output_value: float = 0.7):
        super().__init__()
        self.output_value = output_value
        self._call_count = 0

    def forward(self, tokens, engineered_features=None):
        self._call_count += 1
        batch_size = tokens.shape[0]
        # Return logits (pre-sigmoid)
        logit = torch.log(torch.tensor(self.output_value / (1 - self.output_value)))
        return logit.expand(batch_size, 1)

    def eval(self):
        return self


@pytest.fixture
def vocab():
    """Create a BPE tokenizer for testing."""
    tok = DNABPETokenizer(vocab_size=30)
    tok.train(["ACGTACGTAC", "GGGGCCCCAA", "ATATATAT", "TTTTAAAACCCC"])
    return tok


class TestReverseComplement:
    """Test reverse complement functionality."""

    def test_reverse_complement_basic(self):
        """Test basic reverse complement."""
        assert reverse_complement("ACGT") == "ACGT"  # Palindrome
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GGGG") == "CCCC"
        assert reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_longer_sequence(self):
        """Test reverse complement of longer sequence."""
        seq = "ATGCATGC"
        rc = reverse_complement(seq)
        assert rc == "GCATGCAT"
        # Double RC should give original
        assert reverse_complement(rc) == seq

    def test_reverse_complement_invariant(self):
        """Test that double reverse complement returns original."""
        sequences = ["ACGTACGT", "AAAATTTT", "GCGCGCGC", "ATATATATATAT"]
        for seq in sequences:
            assert reverse_complement(reverse_complement(seq)) == seq


class TestGetReverseComplementTokens:
    """Test tokenized reverse complement."""

    def test_get_rc_tokens_shape(self, vocab):
        """Test that RC tokens have same shape as input."""
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)

        rc_tokens = get_reverse_complement_tokens(tokens, vocab, max_token_len)

        assert rc_tokens.shape == tokens.shape

    def test_get_rc_tokens_batch(self, vocab):
        """Test RC tokens with batch input."""
        max_token_len = 12
        sequences = ["ACGTACGTAC", "AAAATTTTGG"]
        tokens = torch.stack([vocab.tokenize_and_pad(s, max_token_len) for s in sequences])

        rc_tokens = get_reverse_complement_tokens(tokens, vocab, max_token_len)

        assert rc_tokens.shape == tokens.shape

    def test_double_rc_tokens_gives_original(self, vocab):
        """Test that double RC of tokens gives back original sequence."""
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)

        # Apply RC twice
        rc_tokens = get_reverse_complement_tokens(tokens, vocab, max_token_len)
        double_rc_tokens = get_reverse_complement_tokens(rc_tokens, vocab, max_token_len)

        # Decode and compare
        original_decoded = vocab.decode(tokens.squeeze().tolist())
        double_rc_decoded = vocab.decode(double_rc_tokens.squeeze().tolist())

        # Should be the same (or very close - may have padding differences)
        assert original_decoded[:len(seq)] == double_rc_decoded[:len(seq)]


class TestPredictWithTTA:
    """Test TTA prediction function."""

    def test_predict_with_tta_returns_probabilities(self, vocab):
        """Test that TTA returns valid probabilities."""
        model = MockModel(output_value=0.7)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        probs = predict_with_tta(
            model=model,
            tokens=tokens,
            engineered_features=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        assert probs.shape == (1, 1)
        assert 0 <= probs.item() <= 1

    def test_predict_with_tta_calls_model_twice(self, vocab):
        """Test that TTA calls model twice (forward + RC)."""
        model = MockModel(output_value=0.7)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        _ = predict_with_tta(
            model=model,
            tokens=tokens,
            engineered_features=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        assert model._call_count == 2

    def test_predict_with_tta_mean_aggregation(self, vocab):
        """Test mean aggregation of TTA predictions."""
        # Model returns same value for both orientations
        model = MockModel(output_value=0.7)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        probs = predict_with_tta(
            model=model,
            tokens=tokens,
            engineered_features=features,
            vocab=vocab,
            max_token_len=max_token_len,
            aggregation="mean",
        )

        # With same output for both, average should be 0.7
        assert abs(probs.item() - 0.7) < 0.01

    def test_predict_with_tta_geometric_mean_aggregation(self, vocab):
        """Test geometric mean aggregation."""
        model = MockModel(output_value=0.64)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        probs = predict_with_tta(
            model=model,
            tokens=tokens,
            engineered_features=features,
            vocab=vocab,
            max_token_len=max_token_len,
            aggregation="geometric_mean",
        )

        # sqrt(0.64 * 0.64) = 0.64
        assert abs(probs.item() - 0.64) < 0.01


class TestPredictBatchWithTTA:
    """Test batch TTA prediction."""

    def test_batch_tta_returns_three_tensors(self, vocab):
        """Test that batch TTA returns forward, RC, and averaged probs."""
        model = MockModel(output_value=0.7)
        max_token_len = 12
        sequences = ["ACGTACGTAC", "AAAATTTTGG"]
        tokens = torch.stack([vocab.tokenize_and_pad(s, max_token_len) for s in sequences])
        features = torch.randn(2, 288)

        prob_fwd, prob_rc, prob_avg = predict_batch_with_tta(
            model=model,
            tokens_batch=tokens,
            features_batch=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        assert prob_fwd.shape == (2, 1)
        assert prob_rc.shape == (2, 1)
        assert prob_avg.shape == (2, 1)

    def test_batch_tta_average_is_correct(self, vocab):
        """Test that averaged probability is correct."""
        model = MockModel(output_value=0.8)
        max_token_len = 12
        sequences = ["ACGTACGTAC"]
        tokens = torch.stack([vocab.tokenize_and_pad(s, max_token_len) for s in sequences])
        features = torch.randn(1, 288)

        prob_fwd, prob_rc, prob_avg = predict_batch_with_tta(
            model=model,
            tokens_batch=tokens,
            features_batch=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        # Average of 0.8 and 0.8 should be 0.8
        assert abs(prob_avg.item() - 0.8) < 0.01


class TestTTAWrapper:
    """Test TTA wrapper class."""

    def test_wrapper_enabled(self, vocab):
        """Test wrapper with TTA enabled."""
        model = MockModel(output_value=0.7)
        wrapper = TTAWrapper(
            model=model,
            vocab=vocab,
            max_token_len=12,
            enabled=True,
        )

        tokens = vocab.tokenize_and_pad("ACGTACGTAC", 12).unsqueeze(0)
        features = torch.randn(1, 288)

        probs = wrapper(tokens, features)

        assert probs.shape == (1, 1)
        assert model._call_count == 2  # Called for forward and RC

    def test_wrapper_disabled(self, vocab):
        """Test wrapper with TTA disabled."""
        model = MockModel(output_value=0.7)
        wrapper = TTAWrapper(
            model=model,
            vocab=vocab,
            max_token_len=12,
            enabled=False,
        )

        tokens = vocab.tokenize_and_pad("ACGTACGTAC", 12).unsqueeze(0)
        features = torch.randn(1, 288)

        probs = wrapper(tokens, features)

        assert probs.shape == (1, 1)
        assert model._call_count == 1  # Only forward pass

    def test_wrapper_predict_with_details(self, vocab):
        """Test wrapper returns detailed results."""
        model = MockModel(output_value=0.7)
        wrapper = TTAWrapper(
            model=model,
            vocab=vocab,
            max_token_len=12,
            enabled=True,
        )

        tokens = vocab.tokenize_and_pad("ACGTACGTAC", 12).unsqueeze(0)
        features = torch.randn(1, 288)

        result = wrapper.predict_with_details(tokens, features)

        assert "forward" in result
        assert "reverse_complement" in result
        assert "averaged" in result
        assert result["tta_enabled"] is True

    def test_wrapper_details_disabled(self, vocab):
        """Test wrapper details when TTA disabled."""
        model = MockModel(output_value=0.7)
        wrapper = TTAWrapper(
            model=model,
            vocab=vocab,
            max_token_len=12,
            enabled=False,
        )

        tokens = vocab.tokenize_and_pad("ACGTACGTAC", 12).unsqueeze(0)
        features = torch.randn(1, 288)

        result = wrapper.predict_with_details(tokens, features)

        assert result["tta_enabled"] is False
        # forward, rc, and averaged should all be the same
        assert torch.equal(result["forward"], result["averaged"])


class TestTTACorrectness:
    """Test TTA correctness properties."""

    def test_promoter_same_classification_both_strands(self, vocab):
        """A promoter should be classified the same on forward and RC strands.

        This is a key property: if TTA is working correctly, and the model
        has learned strand-invariant features, the forward and RC predictions
        should be similar.
        """
        # Use a model that always outputs high probability
        model = MockModel(output_value=0.9)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        prob_fwd, prob_rc, prob_avg = predict_batch_with_tta(
            model=model,
            tokens_batch=tokens,
            features_batch=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        # Both should be high probability (promoter)
        assert prob_fwd.item() > 0.5
        assert prob_rc.item() > 0.5
        assert prob_avg.item() > 0.5

    def test_tta_improves_confidence_for_consistent_predictions(self, vocab):
        """TTA should maintain confidence when both strands agree."""
        model = MockModel(output_value=0.85)
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)
        features = torch.randn(1, 288)

        prob_fwd, prob_rc, prob_avg = predict_batch_with_tta(
            model=model,
            tokens_batch=tokens,
            features_batch=features,
            vocab=vocab,
            max_token_len=max_token_len,
        )

        # Average should still be high when both agree
        assert prob_avg.item() >= 0.8
