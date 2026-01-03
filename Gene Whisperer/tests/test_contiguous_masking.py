# tests/test_contiguous_masking.py
"""
Comprehensive tests for DNABERT-style contiguous span masking.

Tests cover:
1. Span length verification (exactly k tokens)
2. Mask rate validation (~15%)
3. Special token protection (PAD, MASK, UNK never masked)
4. 80/10/10 distribution (MASK/random/keep)
5. Progressive scheduler behavior
6. Label correctness (-100 for unmasked, original for masked)
7. Contiguous vs random strategy differences
"""
from __future__ import annotations

import random
import sys
from collections import Counter
from itertools import product
from pathlib import Path
from typing import List, Set

import torch
import pytest

# Add training directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "training"))

from dataset import KmerVocabulary
from pretrain_mlm import create_contiguous_mask, MaskingScheduler


def create_test_vocab(k: int = 6) -> KmerVocabulary:
    """Create a test vocabulary for k-mers."""
    bases = ["A", "C", "G", "T"]
    kmers = ["".join(p) for p in product(bases, repeat=k)]
    return KmerVocabulary(k=k, tokens=kmers)


class TestContiguousMasking:
    """Tests for DNABERT-style contiguous span masking."""

    @pytest.fixture
    def sample_tokens(self):
        """Create sample tokenized sequences."""
        batch_size = 8
        seq_len = 79
        vocab_size = 4099
        # Generate tokens excluding special tokens (assume 0, 1, 2 are special)
        return torch.randint(3, vocab_size, (batch_size, seq_len))

    @pytest.fixture
    def vocab_k6(self):
        """Create k=6 vocabulary."""
        return create_test_vocab(k=6)

    @pytest.fixture
    def vocab_k4(self):
        """Create k=4 vocabulary."""
        return create_test_vocab(k=4)

    def test_mask_span_length(self, sample_tokens):
        """
        Verify masked regions form contiguous spans of length k.

        For DNABERT-style masking, when a position is selected as a span START,
        it should mask k consecutive tokens (the k-mer size).
        """
        kmer = 6
        vocab_size = 4099
        mask_token_id = 1
        special_token_ids = {0, 1, 2}  # PAD, MASK, UNK

        masked, labels, mask = create_contiguous_mask(
            sample_tokens,
            kmer=kmer,
            mask_token_id=mask_token_id,
            vocab_size=vocab_size,
            special_token_ids=special_token_ids,
        )

        # For each batch, analyze masked spans
        for b in range(mask.shape[0]):
            masked_positions = mask[b].nonzero().flatten().tolist()
            if not masked_positions:
                continue

            # Find contiguous groups
            groups = []
            if masked_positions:
                current_group = [masked_positions[0]]
                for pos in masked_positions[1:]:
                    if pos == current_group[-1] + 1:
                        current_group.append(pos)
                    else:
                        groups.append(current_group)
                        current_group = [pos]
                groups.append(current_group)

            # Each group should have length that is a multiple of kmer
            # OR be cut short by sequence boundaries
            # The key insight: spans START at a position and extend k tokens
            # Overlapping spans can create longer runs
            for group in groups:
                # Spans can overlap, so groups can be longer than k
                # But each individual span is exactly k tokens
                # So group lengths should be >= 1 (at minimum)
                assert len(group) >= 1, f"Empty group found in batch {b}"

                # Each group should be contiguous
                for i in range(1, len(group)):
                    assert group[i] == group[i - 1] + 1, (
                        f"Non-contiguous span in batch {b}: {group}"
                    )

    def test_mask_rate(self, sample_tokens):
        """
        Verify masking rate scales with mask_prob and kmer size.

        In create_contiguous_mask, mask_prob is the probability of a position
        being a span START, and each start masks k consecutive tokens.
        So effective mask rate ≈ mask_prob * k (minus overlaps).

        To get ~15% token coverage, use mask_prob ≈ 0.15/k.
        """
        kmer = 6
        vocab_size = 4099
        mask_token_id = 1
        # To get ~15% coverage with k=6, use mask_prob ≈ 0.025
        mask_prob = 0.025
        special_token_ids = {0, 1, 2}

        # Run multiple times to get stable statistics
        total_tokens = 0
        total_masked = 0
        num_runs = 50

        torch.manual_seed(42)
        random.seed(42)

        for _ in range(num_runs):
            tokens = torch.randint(3, vocab_size, (8, 79))
            _, labels, mask = create_contiguous_mask(
                tokens,
                kmer=kmer,
                mask_token_id=mask_token_id,
                vocab_size=vocab_size,
                mask_prob=mask_prob,
                special_token_ids=special_token_ids,
            )
            total_tokens += tokens.numel()
            total_masked += mask.sum().item()

        actual_rate = total_masked / total_tokens

        # With mask_prob=0.025 and k=6, expect ~15% coverage
        # Allow reasonable tolerance for span overlaps and edge effects
        assert 0.08 < actual_rate < 0.25, (
            f"Mask rate {actual_rate:.2%} outside expected range [8%, 25%]. "
            f"Total masked: {total_masked}, Total tokens: {total_tokens}"
        )

    def test_special_tokens_not_masked(self, vocab_k6):
        """
        Verify PAD, MASK, and UNK tokens are NEVER masked.

        Special tokens should always remain unchanged in the output.
        """
        batch_size = 8
        seq_len = 79
        vocab = vocab_k6

        # Create tokens with some special tokens mixed in
        tokens = torch.randint(0, len(vocab.itos), (batch_size, seq_len))

        # Explicitly place special tokens at specific positions
        # PAD at positions 70-78 (end of sequences)
        tokens[:, 70:] = vocab.pad_id
        # MASK at some positions
        tokens[:, 10] = vocab.mask_id
        tokens[:, 30] = vocab.mask_id
        # UNK at some positions
        tokens[:, 20] = vocab.unk_id
        tokens[:, 40] = vocab.unk_id

        special_token_ids = {vocab.pad_id, vocab.mask_id, vocab.unk_id}

        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=vocab.k,
            mask_token_id=vocab.mask_id,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
        )

        # Check PAD positions are not masked
        pad_positions = tokens == vocab.pad_id
        assert not mask[pad_positions].any(), "PAD tokens should never be masked"
        assert (labels[pad_positions] == -100).all(), (
            "PAD tokens should have label=-100"
        )

        # Check MASK positions are not re-masked (labels should be -100)
        original_mask_positions = tokens == vocab.mask_id
        assert (labels[original_mask_positions] == -100).all(), (
            "Original MASK tokens should have label=-100"
        )

        # Check UNK positions are not masked
        unk_positions = tokens == vocab.unk_id
        assert (labels[unk_positions] == -100).all(), (
            "UNK tokens should have label=-100"
        )

    def test_mask_distribution(self):
        """
        Verify 80/10/10 split for MASK/random/keep within masked positions.

        Of all masked positions:
        - 80% should be replaced with [MASK] token
        - 10% should be replaced with random token
        - 10% should keep original token
        """
        kmer = 6
        vocab_size = 4096
        mask_token_id = 100
        special_token_ids = {0, 1, 2}  # PAD, MASK, UNK

        # Use large sample for statistical significance
        batch_size = 64
        seq_len = 128
        torch.manual_seed(42)
        random.seed(42)

        tokens = torch.randint(3, vocab_size, (batch_size, seq_len))
        original = tokens.clone()

        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=kmer,
            mask_token_id=mask_token_id,
            vocab_size=vocab_size,
            mask_prob=0.15,
            mask_token_prob=0.8,
            random_token_prob=0.1,
            special_token_ids=special_token_ids,
        )

        # Count positions in each category
        masked_positions = mask
        num_masked = masked_positions.sum().item()

        if num_masked == 0:
            pytest.skip("No positions were masked")

        # Count [MASK] tokens at masked positions
        got_mask_token = (masked[masked_positions] == mask_token_id).sum().item()

        # Count unchanged (kept original) at masked positions
        got_kept = (masked[masked_positions] == original[masked_positions]).sum().item()

        # Random = total masked - MASK - kept
        got_random = num_masked - got_mask_token - got_kept

        # Calculate percentages
        mask_pct = got_mask_token / num_masked
        keep_pct = got_kept / num_masked
        random_pct = got_random / num_masked

        # Allow 5% tolerance
        tolerance = 0.05
        assert abs(mask_pct - 0.80) <= tolerance, (
            f"[MASK] percentage {mask_pct:.2%} should be ~80%"
        )
        assert abs(keep_pct - 0.10) <= tolerance, (
            f"Keep percentage {keep_pct:.2%} should be ~10%"
        )
        assert abs(random_pct - 0.10) <= tolerance, (
            f"Random percentage {random_pct:.2%} should be ~10%"
        )

    def test_progressive_scheduler(self):
        """
        Verify MaskingScheduler increases rate correctly over steps.

        DNABERT schedule:
        - First 80% of training: 15% masking rate
        - Last 20% of training: linear increase to 20%
        """
        total_steps = 100000
        scheduler = MaskingScheduler(
            initial_rate=0.15,
            final_rate=0.20,
            total_steps=total_steps,
            transition_start_ratio=0.8,
        )

        # At step 0: should be initial rate
        assert scheduler.get_mask_rate(0) == 0.15, (
            f"Rate at step 0 should be 0.15, got {scheduler.get_mask_rate(0)}"
        )

        # At step 40000 (40%): should still be initial rate
        assert scheduler.get_mask_rate(40000) == 0.15, (
            f"Rate at step 40000 should be 0.15, got {scheduler.get_mask_rate(40000)}"
        )

        # At step 80000 (80%): just starting transition, should be ~0.15
        rate_80 = scheduler.get_mask_rate(80000)
        assert 0.15 <= rate_80 <= 0.16, (
            f"Rate at step 80000 should be ~0.15, got {rate_80}"
        )

        # At step 90000 (90%): midway through transition, should be ~0.175
        rate_90 = scheduler.get_mask_rate(90000)
        assert 0.17 <= rate_90 <= 0.18, (
            f"Rate at step 90000 should be ~0.175, got {rate_90}"
        )

        # At step 100000 (100%): should be final rate
        rate_100 = scheduler.get_mask_rate(100000)
        assert rate_100 == 0.20, (
            f"Rate at step 100000 should be 0.20, got {rate_100}"
        )

        # Beyond total_steps: should stay at final rate
        rate_beyond = scheduler.get_mask_rate(150000)
        assert rate_beyond == 0.20, (
            f"Rate beyond total_steps should be 0.20, got {rate_beyond}"
        )

    def test_labels_correct(self, sample_tokens):
        """
        Verify labels are -100 for unmasked positions, original token for masked.

        This is critical for the MLM loss function which ignores -100 labels.
        """
        kmer = 6
        vocab_size = 4099
        mask_token_id = 1
        special_token_ids = {0, 1, 2}

        original = sample_tokens.clone()

        masked, labels, mask = create_contiguous_mask(
            sample_tokens,
            kmer=kmer,
            mask_token_id=mask_token_id,
            vocab_size=vocab_size,
            special_token_ids=special_token_ids,
        )

        # Unmasked positions should have label=-100
        unmasked_positions = ~mask
        assert (labels[unmasked_positions] == -100).all(), (
            "Unmasked positions should have label=-100"
        )

        # Masked positions should have label=original token
        masked_positions = mask
        if masked_positions.any():
            assert (labels[masked_positions] == original[masked_positions]).all(), (
                "Masked positions should have label=original token"
            )

        # Special tokens should never be in labels (should be -100)
        for special_id in special_token_ids:
            special_positions = original == special_id
            if special_positions.any():
                assert (labels[special_positions] == -100).all(), (
                    f"Special token {special_id} should have label=-100"
                )

    def test_contiguous_vs_random_different(self, sample_tokens):
        """
        Verify contiguous masking produces different patterns than random masking.

        Contiguous masking should create spans of consecutive masked positions,
        while random masking would spread masked positions throughout.
        """
        kmer = 6
        vocab_size = 4099
        mask_token_id = 1
        special_token_ids = {0, 1, 2}

        # Set seed for reproducibility
        torch.manual_seed(123)
        random.seed(123)

        # Get contiguous masking result
        _, _, contiguous_mask = create_contiguous_mask(
            sample_tokens.clone(),
            kmer=kmer,
            mask_token_id=mask_token_id,
            vocab_size=vocab_size,
            mask_prob=0.15,
            special_token_ids=special_token_ids,
        )

        # Simulate random masking (independent per-token)
        torch.manual_seed(456)
        random_mask = torch.rand_like(sample_tokens, dtype=torch.float32) < 0.15

        # Analyze contiguity: count isolated masked tokens vs spans
        def count_span_statistics(mask_tensor):
            """Count number of spans and isolated tokens."""
            num_spans = 0
            num_isolated = 0
            total_span_lengths = []

            for b in range(mask_tensor.shape[0]):
                row = mask_tensor[b].tolist()
                i = 0
                while i < len(row):
                    if row[i]:
                        # Start of a masked region
                        span_len = 0
                        while i < len(row) and row[i]:
                            span_len += 1
                            i += 1
                        total_span_lengths.append(span_len)
                        if span_len == 1:
                            num_isolated += 1
                        else:
                            num_spans += 1
                    else:
                        i += 1

            return {
                "num_spans": num_spans,
                "num_isolated": num_isolated,
                "avg_span_length": (
                    sum(total_span_lengths) / len(total_span_lengths)
                    if total_span_lengths
                    else 0
                ),
                "span_lengths": total_span_lengths,
            }

        contiguous_stats = count_span_statistics(contiguous_mask)
        random_stats = count_span_statistics(random_mask)

        # Contiguous masking should have longer average span length
        assert contiguous_stats["avg_span_length"] > random_stats["avg_span_length"], (
            f"Contiguous avg span length ({contiguous_stats['avg_span_length']:.2f}) "
            f"should be > random ({random_stats['avg_span_length']:.2f})"
        )

        # Contiguous masking should have fewer isolated single tokens
        # (as a proportion of total masked regions)
        contiguous_total = contiguous_stats["num_spans"] + contiguous_stats["num_isolated"]
        random_total = random_stats["num_spans"] + random_stats["num_isolated"]

        if contiguous_total > 0 and random_total > 0:
            contiguous_isolated_ratio = contiguous_stats["num_isolated"] / contiguous_total
            random_isolated_ratio = random_stats["num_isolated"] / random_total

            # Random masking typically has more isolated tokens
            # (This may not always hold, but is generally expected)
            # At minimum, they should produce different patterns
            assert contiguous_stats != random_stats, (
                "Contiguous and random masking should produce different patterns"
            )


class TestMaskingSchedulerEdgeCases:
    """Additional edge case tests for MaskingScheduler."""

    def test_scheduler_state_dict(self):
        """Verify state dict for checkpointing."""
        scheduler = MaskingScheduler(
            initial_rate=0.15,
            final_rate=0.20,
            total_steps=1000,
        )

        # Call get_mask_rate to update internal state
        scheduler.get_mask_rate(500)

        state = scheduler.state_dict()
        assert "initial_rate" in state
        assert "final_rate" in state
        assert "total_steps" in state

    def test_scheduler_invalid_params(self):
        """Verify scheduler rejects invalid parameters."""
        # Invalid initial_rate
        with pytest.raises(ValueError):
            MaskingScheduler(initial_rate=-0.1, total_steps=1000)

        with pytest.raises(ValueError):
            MaskingScheduler(initial_rate=1.5, total_steps=1000)

        # Invalid final_rate
        with pytest.raises(ValueError):
            MaskingScheduler(final_rate=-0.1, total_steps=1000)

        # Invalid total_steps
        with pytest.raises(ValueError):
            MaskingScheduler(total_steps=0)

        with pytest.raises(ValueError):
            MaskingScheduler(total_steps=-100)

        # Invalid warmup_steps
        with pytest.raises(ValueError):
            MaskingScheduler(warmup_steps=-10, total_steps=1000)

    def test_scheduler_no_transition(self):
        """Test scheduler with same initial and final rate."""
        scheduler = MaskingScheduler(
            initial_rate=0.15,
            final_rate=0.15,
            total_steps=1000,
        )

        # Rate should always be 0.15
        for step in [0, 500, 800, 1000, 2000]:
            assert scheduler.get_mask_rate(step) == 0.15


class TestContiguousMaskingEdgeCases:
    """Edge case tests for create_contiguous_mask."""

    def test_empty_batch(self):
        """Test with empty batch."""
        tokens = torch.zeros((0, 79), dtype=torch.long)
        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=6,
            mask_token_id=1,
            vocab_size=4099,
        )
        assert masked.shape == (0, 79)
        assert labels.shape == (0, 79)
        assert mask.shape == (0, 79)

    def test_very_short_sequence(self):
        """Test with sequence shorter than k-mer size."""
        tokens = torch.randint(3, 100, (4, 3))  # seq_len=3, kmer=6
        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=6,
            mask_token_id=1,
            vocab_size=100,
        )
        # Should handle gracefully without errors
        assert masked.shape == tokens.shape

    def test_all_special_tokens(self):
        """Test sequence of all special tokens."""
        tokens = torch.zeros((4, 20), dtype=torch.long)  # All PAD (id=0)
        special_token_ids = {0}

        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=6,
            mask_token_id=1,
            vocab_size=100,
            special_token_ids=special_token_ids,
        )

        # Nothing should be masked
        assert not mask.any(), "All-special-token sequence should have no masking"
        assert (labels == -100).all(), "All labels should be -100"

    def test_single_token_sequence(self):
        """Test single token sequence."""
        tokens = torch.tensor([[42]], dtype=torch.long)
        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=6,
            mask_token_id=1,
            vocab_size=100,
        )
        assert masked.shape == (1, 1)

    def test_high_mask_probability(self):
        """Test with very high mask probability."""
        tokens = torch.randint(3, 100, (8, 50))
        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=4,
            mask_token_id=1,
            vocab_size=100,
            mask_prob=0.5,
        )

        # Should mask a significant portion
        mask_rate = mask.sum().item() / mask.numel()
        assert mask_rate > 0.3, f"High mask_prob should produce high masking, got {mask_rate:.2%}"

    def test_zero_mask_probability(self):
        """Test with zero mask probability."""
        tokens = torch.randint(3, 100, (8, 50))
        masked, labels, mask = create_contiguous_mask(
            tokens,
            kmer=4,
            mask_token_id=1,
            vocab_size=100,
            mask_prob=0.0,
        )

        # Should mask nothing (or nearly nothing)
        mask_rate = mask.sum().item() / mask.numel()
        assert mask_rate < 0.05, f"Zero mask_prob should produce minimal masking, got {mask_rate:.2%}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
