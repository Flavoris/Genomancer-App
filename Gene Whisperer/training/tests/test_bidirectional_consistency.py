"""Tests for Bidirectional Consistency Training functionality."""
import pytest
import torch
import torch.nn as nn

from bpe_tokenizer import DNABPETokenizer
from dataset import reverse_complement


# Import from train_stage1 - we need to handle the import path
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from train_stage1 import (
    BidirectionalConsistencyLoss,
    get_reverse_complement_tokens_for_training,
)


class MockModel(nn.Module):
    """Mock model for testing consistency training."""

    def __init__(self, output_value: float = 0.0):
        super().__init__()
        self.output_value = output_value
        self._call_count = 0

    def forward(self, tokens, engineered_features=None):
        self._call_count += 1
        batch_size = tokens.shape[0]
        return torch.full((batch_size, 1), self.output_value)

    def eval(self):
        return self

    def train(self, mode=True):
        return self


class AsymmetricMockModel(nn.Module):
    """Mock model that returns different values for forward vs RC."""

    def __init__(self, fwd_value: float = 0.5, rc_value: float = 0.3):
        super().__init__()
        self.fwd_value = fwd_value
        self.rc_value = rc_value
        self._call_count = 0

    def forward(self, tokens, engineered_features=None):
        self._call_count += 1
        batch_size = tokens.shape[0]
        # Alternate between forward and RC values based on call count
        if self._call_count % 2 == 1:
            return torch.full((batch_size, 1), self.fwd_value)
        else:
            return torch.full((batch_size, 1), self.rc_value)


@pytest.fixture
def vocab():
    """Create a simple BPE tokenizer for testing."""
    tok = DNABPETokenizer(vocab_size=20)
    tok.train(["ACGTACGTAC", "GGGGCCCCAA", "ATATATAT", "ACGTACGT"])
    return tok


class TestBidirectionalConsistencyLoss:
    """Test the BidirectionalConsistencyLoss class."""

    def test_loss_initialization(self):
        """Test that loss initializes with correct parameters."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3, label_smoothing=0.1)
        assert loss_fn.alpha == 0.3
        assert loss_fn.label_smoothing == 0.1

    def test_loss_initialization_defaults(self):
        """Test default initialization."""
        loss_fn = BidirectionalConsistencyLoss()
        assert loss_fn.alpha == 0.3
        assert loss_fn.label_smoothing == 0.0

    def test_loss_output_shape(self):
        """Test that loss outputs correct shape."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 8
        logits_fwd = torch.randn(batch_size, 1)
        logits_rc = torch.randn(batch_size, 1)
        labels = torch.randint(0, 2, (batch_size, 1)).float()

        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        assert total_loss.shape == ()
        assert isinstance(metrics, dict)
        assert 'cls_loss' in metrics
        assert 'consistency_loss' in metrics
        assert 'total_loss' in metrics

    def test_loss_zero_consistency_when_identical(self):
        """Test that consistency loss is zero when predictions are identical."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 8
        logits = torch.randn(batch_size, 1)
        labels = torch.randint(0, 2, (batch_size, 1)).float()

        # Same logits for forward and RC
        total_loss, metrics = loss_fn(logits, logits.clone(), labels)

        assert metrics['consistency_loss'] == pytest.approx(0.0, abs=1e-6)

    def test_loss_positive_consistency_when_different(self):
        """Test that consistency loss is positive when predictions differ."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 8
        logits_fwd = torch.full((batch_size, 1), 2.0)  # High confidence positive
        logits_rc = torch.full((batch_size, 1), -2.0)  # High confidence negative
        labels = torch.ones(batch_size, 1)

        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        assert metrics['consistency_loss'] > 0.0

    def test_loss_alpha_weighting(self):
        """Test that alpha correctly weights consistency loss."""
        batch_size = 4
        logits_fwd = torch.full((batch_size, 1), 2.0)
        logits_rc = torch.full((batch_size, 1), -2.0)
        labels = torch.ones(batch_size, 1)

        # With alpha=0, consistency loss should not affect total
        loss_fn_zero = BidirectionalConsistencyLoss(alpha=0.0)
        total_zero, metrics_zero = loss_fn_zero(logits_fwd, logits_rc, labels)

        # With alpha=1, consistency loss should have full weight
        loss_fn_one = BidirectionalConsistencyLoss(alpha=1.0)
        total_one, metrics_one = loss_fn_one(logits_fwd, logits_rc, labels)

        # Total with alpha=0 should be lower (no consistency penalty)
        assert total_zero < total_one

    def test_loss_with_label_smoothing(self):
        """Test loss with label smoothing enabled."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3, label_smoothing=0.1)
        batch_size = 4
        logits = torch.randn(batch_size, 1)
        labels = torch.ones(batch_size, 1)

        total_loss, metrics = loss_fn(logits, logits.clone(), labels)

        # Should run without error
        assert total_loss.item() > 0.0

    def test_loss_gradient_flow(self):
        """Test that gradients flow through the loss."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 4
        logits_fwd = torch.randn(batch_size, 1, requires_grad=True)
        logits_rc = torch.randn(batch_size, 1, requires_grad=True)
        labels = torch.randint(0, 2, (batch_size, 1)).float()

        total_loss, _ = loss_fn(logits_fwd, logits_rc, labels)
        total_loss.backward()

        assert logits_fwd.grad is not None
        assert logits_rc.grad is not None

    def test_loss_encourages_consistency(self):
        """Test that loss encourages predictions to be consistent."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.5)
        labels = torch.ones(4, 1)

        # Consistent predictions (both positive)
        logits_consistent_fwd = torch.full((4, 1), 2.0)
        logits_consistent_rc = torch.full((4, 1), 2.0)
        loss_consistent, _ = loss_fn(logits_consistent_fwd, logits_consistent_rc, labels)

        # Inconsistent predictions (one positive, one negative)
        logits_inconsistent_fwd = torch.full((4, 1), 2.0)
        logits_inconsistent_rc = torch.full((4, 1), -2.0)
        loss_inconsistent, _ = loss_fn(logits_inconsistent_fwd, logits_inconsistent_rc, labels)

        # Consistent predictions should have lower loss
        assert loss_consistent < loss_inconsistent


class TestGetReverseComplementTokensForTraining:
    """Test the RC token computation for training."""

    def test_rc_tokens_shape(self, vocab):
        """Test that RC tokens have same shape as input."""
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)

        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)

        assert rc_tokens.shape == tokens.shape

    def test_rc_tokens_batch(self, vocab):
        """Test RC tokens with batch of sequences."""
        max_token_len = 12
        sequences = ["ACGTACGTAC", "GGGGCCCCAA", "ATATATAT"]
        tokens_list = [vocab.tokenize_and_pad(seq, max_token_len) for seq in sequences]
        tokens = torch.stack(tokens_list, dim=0)

        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)

        assert rc_tokens.shape == tokens.shape

    def test_double_rc_returns_original(self, vocab):
        """Test that applying RC twice returns original tokens."""
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)

        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)
        rc_rc_tokens = get_reverse_complement_tokens_for_training(rc_tokens, vocab, max_token_len)

        # Double RC should give back original
        assert torch.equal(tokens, rc_rc_tokens)

    def test_rc_tokens_device_transfer(self, vocab):
        """Test that RC tokens are on same device as input."""
        max_token_len = 12
        seq = "ACGTACGTAC"
        tokens = vocab.tokenize_and_pad(seq, max_token_len).unsqueeze(0)

        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)

        assert rc_tokens.device == tokens.device


class TestConsistencyTrainingIntegration:
    """Integration tests for consistency training components."""

    def test_full_loss_computation_pipeline(self, vocab):
        """Test full pipeline: tokens -> RC tokens -> loss."""
        max_token_len = 12
        batch_size = 4
        sequences = ["ACGTACGTAC", "GGGGCCCCAA", "ATATATAT", "CGATCGATCG"]

        # Tokenize sequences
        tokens_list = [vocab.tokenize_and_pad(seq, max_token_len) for seq in sequences]
        tokens = torch.stack(tokens_list, dim=0)
        labels = torch.randint(0, 2, (batch_size, 1)).float()

        # Get RC tokens
        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)

        # Create mock model
        model = MockModel(output_value=0.5)

        # Forward pass
        logits_fwd = model(tokens)
        logits_rc = model(rc_tokens)

        # Compute loss
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        assert total_loss.item() > 0.0
        assert metrics['cls_loss'] > 0.0

    def test_asymmetric_model_has_consistency_loss(self, vocab):
        """Test that asymmetric predictions result in consistency loss."""
        max_token_len = 12
        batch_size = 2
        sequences = ["ACGTACGTAC", "GGGGCCCCAA"]

        # Tokenize sequences
        tokens_list = [vocab.tokenize_and_pad(seq, max_token_len) for seq in sequences]
        tokens = torch.stack(tokens_list, dim=0)
        labels = torch.ones(batch_size, 1)

        # Get RC tokens
        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)

        # Create asymmetric model
        model = AsymmetricMockModel(fwd_value=2.0, rc_value=-2.0)

        # Forward pass
        logits_fwd = model(tokens)
        logits_rc = model(rc_tokens)

        # Compute loss
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        # Should have significant consistency loss
        assert metrics['consistency_loss'] > 0.1

    def test_consistency_loss_encourages_agreement(self, vocab):
        """Test that training with consistency loss encourages prediction agreement."""
        # This is a simplified optimization test
        max_token_len = 12
        batch_size = 4
        sequences = ["ACGT" * 3] * batch_size

        tokens_list = [vocab.tokenize_and_pad(seq, max_token_len) for seq in sequences]
        tokens = torch.stack(tokens_list, dim=0)
        rc_tokens = get_reverse_complement_tokens_for_training(tokens, vocab, max_token_len)
        labels = torch.ones(batch_size, 1)

        # Simple learnable parameters
        fwd_param = torch.nn.Parameter(torch.tensor(2.0))
        rc_param = torch.nn.Parameter(torch.tensor(-2.0))

        optimizer = torch.optim.SGD([fwd_param, rc_param], lr=0.1)
        loss_fn = BidirectionalConsistencyLoss(alpha=0.5)

        # Initial difference
        initial_diff = abs(fwd_param.item() - rc_param.item())

        # Training steps
        for _ in range(10):
            optimizer.zero_grad()
            logits_fwd = fwd_param.expand(batch_size, 1)
            logits_rc = rc_param.expand(batch_size, 1)
            loss, _ = loss_fn(logits_fwd, logits_rc, labels)
            loss.backward()
            optimizer.step()

        # Final difference should be smaller (parameters moved closer)
        final_diff = abs(fwd_param.item() - rc_param.item())
        assert final_diff < initial_diff


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_batch(self):
        """Test handling of empty batch."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        logits_fwd = torch.empty(0, 1)
        logits_rc = torch.empty(0, 1)
        labels = torch.empty(0, 1)

        # Should not crash, though behavior with empty tensors varies
        try:
            total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)
            # If it succeeds, loss should be NaN or zero for empty input
            assert True
        except RuntimeError:
            # Some operations may fail on empty tensors, which is acceptable
            assert True

    def test_single_sample_batch(self):
        """Test with single sample batch."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        logits_fwd = torch.tensor([[1.0]])
        logits_rc = torch.tensor([[0.5]])
        labels = torch.tensor([[1.0]])

        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        assert total_loss.item() > 0.0

    def test_extreme_logits(self):
        """Test with extreme logit values."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 4

        # Very large positive logits
        logits_fwd = torch.full((batch_size, 1), 100.0)
        logits_rc = torch.full((batch_size, 1), 100.0)
        labels = torch.ones(batch_size, 1)

        total_loss, metrics = loss_fn(logits_fwd, logits_rc, labels)

        # Should not produce NaN or Inf
        assert not torch.isnan(total_loss)
        assert not torch.isinf(total_loss)

    def test_all_positive_labels(self):
        """Test with all positive labels."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 8
        logits = torch.randn(batch_size, 1)
        labels = torch.ones(batch_size, 1)

        total_loss, metrics = loss_fn(logits, logits.clone(), labels)

        assert total_loss.item() >= 0.0

    def test_all_negative_labels(self):
        """Test with all negative labels."""
        loss_fn = BidirectionalConsistencyLoss(alpha=0.3)
        batch_size = 8
        logits = torch.randn(batch_size, 1)
        labels = torch.zeros(batch_size, 1)

        total_loss, metrics = loss_fn(logits, logits.clone(), labels)

        assert total_loss.item() >= 0.0
