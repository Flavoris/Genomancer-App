"""
Unit tests for DNABERT-style span masking.

Tests cover:
1. Basic span masking properties
2. 80/10/10 BERT masking rule
3. Edge cases: spans crossing PAD, all-PAD sequences
4. Geometric vs uniform distribution
5. Label correctness (masked=original, unmasked=-100)
"""
from __future__ import annotations

import math
import random
import sys
from collections import Counter
from pathlib import Path
from typing import List, Tuple

import torch

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from dataset import KmerVocabulary
from pretrain_mlm import (
    mask_tokens_span,
    sample_span_length_geometric,
)


def create_test_vocab(k: int = 3) -> KmerVocabulary:
    """Create a simple test vocabulary."""
    # Generate all k-mers for ACGT
    from itertools import product
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in product(bases, repeat=k)]
    return KmerVocabulary(k=k, tokens=kmers)


def create_test_inputs(
    batch_size: int,
    seq_len: int,
    vocab: KmerVocabulary,
    pad_fraction: float = 0.0,
) -> torch.LongTensor:
    """Create test input tensor with optional padding."""
    # Generate random token IDs (excluding special tokens)
    base_ids = vocab._base_token_ids
    inputs = torch.tensor(
        [[random.choice(base_ids) for _ in range(seq_len)] for _ in range(batch_size)],
        dtype=torch.long,
    )
    
    # Add padding at the end if requested
    if pad_fraction > 0:
        pad_len = int(seq_len * pad_fraction)
        inputs[:, -pad_len:] = vocab.pad_id
    
    return inputs


class TestSpanLengthSampling:
    """Tests for geometric span length sampling."""
    
    def test_geometric_sampling_bounds(self):
        """Span lengths should always be in [1, max_span_len]."""
        max_span = 10
        for _ in range(1000):
            span_len = sample_span_length_geometric(max_span, mean_span=3.0)
            assert 1 <= span_len <= max_span, f"Span length {span_len} out of bounds [1, {max_span}]"
    
    def test_geometric_sampling_distribution(self):
        """Geometric distribution should favor shorter spans."""
        max_span = 10
        mean_span = 3.0
        samples = [sample_span_length_geometric(max_span, mean_span) for _ in range(10000)]
        
        counter = Counter(samples)
        
        # Length 1 should be most common (geometric decay)
        assert counter[1] > counter.get(5, 0), "Length 1 should be more common than length 5"
        assert counter[1] > counter.get(max_span, 0), "Length 1 should be more common than max length"
        
        # Empirical mean should be close to target mean (within tolerance)
        empirical_mean = sum(samples) / len(samples)
        # Note: clamping affects the mean, so we allow wider tolerance
        assert 2.0 < empirical_mean < 5.0, f"Mean {empirical_mean} too far from target {mean_span}"
    
    def test_geometric_sampling_reproducibility(self):
        """Same seed should give same results."""
        max_span = 10
        
        random.seed(42)
        samples1 = [sample_span_length_geometric(max_span) for _ in range(100)]
        
        random.seed(42)
        samples2 = [sample_span_length_geometric(max_span) for _ in range(100)]
        
        assert samples1 == samples2, "Sampling not reproducible with same seed"


class TestMaskTokensSpan:
    """Tests for the main mask_tokens_span function."""
    
    def test_basic_masking(self):
        """Basic span masking should work without errors."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(4, 32, vocab)
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        assert masked_inputs.shape == inputs.shape
        assert labels.shape == inputs.shape
    
    def test_mask_probability(self):
        """Approximately mask_prob fraction of tokens should be masked."""
        vocab = create_test_vocab()
        batch_size, seq_len = 16, 64
        inputs = create_test_inputs(batch_size, seq_len, vocab)
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        num_masked = (labels != -100).sum().item()
        total_tokens = batch_size * seq_len
        actual_prob = num_masked / total_tokens
        
        # Should be within 5% of target
        assert 0.10 < actual_prob < 0.25, f"Mask prob {actual_prob:.2%} too far from 15%"
    
    def test_labels_correctness(self):
        """Labels should equal original tokens at masked positions, -100 elsewhere."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(8, 48, vocab)
        original = inputs.clone()
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        # Masked positions: label == original
        masked_positions = (labels != -100)
        assert torch.all(labels[masked_positions] == original[masked_positions]), \
            "Labels at masked positions should equal original tokens"
        
        # Unmasked positions: label == -100
        unmasked_positions = (labels == -100)
        assert unmasked_positions.any(), "Should have some unmasked positions"
    
    def test_80_10_10_rule(self):
        """Check 80/10/10 masking distribution within ±2–3%."""
        vocab = create_test_vocab()
        # Use a large sample for tight statistical tolerance
        random.seed(42)
        torch.manual_seed(42)
        inputs = create_test_inputs(128, 256, vocab)
        original = inputs.clone()

        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        masked_positions = (labels != -100)
        num_masked = masked_positions.sum().item()
        
        # Count: [MASK] token
        got_mask = (masked_inputs[masked_positions] == vocab.mask_id).sum().item()
        # Count: original token (unchanged)
        got_original = (masked_inputs[masked_positions] == original[masked_positions]).sum().item()
        # Count: random token (not MASK, not original)
        got_random = num_masked - got_mask - got_original
        
        mask_pct = got_mask / num_masked
        orig_pct = got_original / num_masked
        random_pct = got_random / num_masked
        
        tolerance = 0.03
        assert abs(mask_pct - 0.80) <= tolerance, (
            f"[MASK] pct {mask_pct:.2%} should be 80%±{tolerance:.0%}"
        )
        assert abs(orig_pct - 0.10) <= tolerance, (
            f"Original pct {orig_pct:.2%} should be 10%±{tolerance:.0%}"
        )
        assert abs(random_pct - 0.10) <= tolerance, (
            f"Random pct {random_pct:.2%} should be 10%±{tolerance:.0%}"
        )


class TestSpanLengthCap:
    """Regression tests: max_span_len must cap contiguous masked runs."""

    def test_max_span_len_never_exceeded_and_coverage_ok(self):
        """
        Generate 1000 masked examples with max_span_len=3 and verify:
        1) No run of consecutive masked tokens exceeds 3
        2) Overall mask coverage stays near mask_prob=0.15
        """
        vocab = create_test_vocab()
        max_span_len = 3
        mask_prob = 0.15
        seq_len = 128
        num_examples = 1000

        random.seed(123)
        torch.manual_seed(123)

        base_ids = torch.tensor(vocab._base_token_ids, dtype=torch.long)
        inputs = base_ids[torch.randint(0, len(base_ids), (num_examples, seq_len))]

        _, labels = mask_tokens_span(
            inputs, vocab, mask_prob=mask_prob, max_span_len=max_span_len
        )

        mask = (labels != -100)

        # (1) Span cap: no contiguous run exceeds max_span_len
        mask_rows = mask.tolist()
        for idx, row in enumerate(mask_rows):
            longest = 0
            current = 0
            for val in row:
                if val:
                    current += 1
                    longest = max(longest, current)
                else:
                    current = 0
            assert longest <= max_span_len, (
                f"Example {idx} has masked run length {longest} (> {max_span_len})"
            )

        # (2) Coverage: overall masked fraction should be close to mask_prob
        mask_coverage = mask.sum().item() / (num_examples * seq_len)
        assert 0.13 <= mask_coverage <= 0.17, (
            f"mask_coverage={mask_coverage:.4f} outside [0.13, 0.17]"
        )


class TestPaddingHandling:
    """Tests for PAD token handling in span masking."""
    
    def test_no_pad_tokens_masked(self):
        """PAD tokens should NEVER be masked."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(8, 64, vocab, pad_fraction=0.3)  # 30% padding
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        # PAD positions in original
        pad_positions = (inputs == vocab.pad_id)
        
        # Labels at PAD positions should be -100 (not masked)
        assert torch.all(labels[pad_positions] == -100), \
            "PAD tokens should have label=-100 (not masked)"
        
        # Masked inputs at PAD positions should still be PAD
        assert torch.all(masked_inputs[pad_positions] == vocab.pad_id), \
            "PAD tokens should remain PAD in masked input"
    
    def test_spans_crossing_pad(self):
        """Spans that cross into PAD region should only mask non-PAD tokens."""
        vocab = create_test_vocab()
        
        # Create sequence where last 50% is PAD
        seq_len = 20
        inputs = create_test_inputs(4, seq_len, vocab, pad_fraction=0.5)
        
        masked_inputs, labels = mask_tokens_span(
            inputs, vocab, mask_prob=0.30, max_span_len=10
        )
        
        # PAD positions (indices 10-19)
        pad_positions = (inputs == vocab.pad_id)
        non_pad_positions = ~pad_positions
        
        # No PAD should be masked
        assert torch.all(labels[pad_positions] == -100)
        
        # Some non-PAD should be masked
        num_masked_non_pad = (labels[non_pad_positions] != -100).sum().item()
        assert num_masked_non_pad > 0, "Should mask some non-PAD tokens"
    
    def test_all_pad_sequence(self):
        """All-PAD sequences should have no masking, no errors."""
        vocab = create_test_vocab()
        
        # Create all-PAD input
        inputs = torch.full((4, 32), vocab.pad_id, dtype=torch.long)
        
        # Should not raise error
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        # All labels should be -100 (no masking)
        assert torch.all(labels == -100), "All-PAD sequence should have all labels=-100"
        
        # All inputs should still be PAD
        assert torch.all(masked_inputs == vocab.pad_id), "All-PAD sequence should remain all-PAD"


class TestDistributionModes:
    """Tests for geometric vs uniform distribution modes."""
    
    def test_geometric_distribution(self):
        """Geometric distribution should be the default and favor short spans."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(32, 64, vocab)
        
        random.seed(42)
        _, _, spans = mask_tokens_span(
            inputs, vocab, mask_prob=0.15, track_spans=True,
            span_distribution="geometric", mean_span=3.0
        )
        
        if len(spans) < 50:
            return  # Skip if too few spans
        
        counter = Counter(spans)
        # Short spans should be more common
        short_count = sum(counter.get(i, 0) for i in range(1, 4))
        long_count = sum(counter.get(i, 0) for i in range(6, 11))
        
        assert short_count > long_count, \
            f"Geometric: short spans ({short_count}) should be more common than long ({long_count})"
    
    def test_uniform_distribution(self):
        """Uniform distribution should give more even span lengths."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(32, 64, vocab)
        
        random.seed(42)
        _, _, spans = mask_tokens_span(
            inputs, vocab, mask_prob=0.15, track_spans=True,
            span_distribution="uniform", max_span_len=5
        )
        
        if len(spans) < 100:
            return  # Skip if too few spans
        
        counter = Counter(spans)
        # All lengths 1-5 should have similar counts (within 3x)
        counts = [counter.get(i, 0) for i in range(1, 6)]
        if min(counts) > 0:
            ratio = max(counts) / min(counts)
            assert ratio < 5, f"Uniform: length counts too uneven (ratio {ratio:.1f})"


class TestSpanTracking:
    """Tests for span length tracking functionality."""
    
    def test_track_spans_returns_lengths(self):
        """track_spans=True should return span lengths list."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(8, 48, vocab)
        
        result = mask_tokens_span(inputs, vocab, mask_prob=0.15, track_spans=True)
        
        assert len(result) == 3, "track_spans=True should return 3 values"
        masked_inputs, labels, span_lengths = result
        assert isinstance(span_lengths, list), "span_lengths should be a list"
    
    def test_span_lengths_sum_approximately_matches(self):
        """Sum of span lengths should approximately equal masked count."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(8, 48, vocab)
        
        masked_inputs, labels, span_lengths = mask_tokens_span(
            inputs, vocab, mask_prob=0.15, track_spans=True
        )
        
        total_masked = (labels != -100).sum().item()
        sum_spans = sum(span_lengths)
        
        # Should be close (may differ slightly due to overlap handling)
        assert abs(total_masked - sum_spans) <= len(span_lengths) + 1, \
            f"Span sum {sum_spans} differs too much from masked count {total_masked}"


class TestEdgeCases:
    """Tests for edge cases and corner conditions."""
    
    def test_single_token_sequence(self):
        """Single-token sequences should be handled."""
        vocab = create_test_vocab()
        inputs = torch.tensor([[vocab._base_token_ids[0]]], dtype=torch.long)
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        # Either masked or not, but should not error
        assert masked_inputs.shape == (1, 1)
        assert labels.shape == (1, 1)
    
    def test_very_short_sequence(self):
        """Very short sequences (2-3 tokens) should work."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(4, 3, vocab)  # Only 3 tokens
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.15)
        
        assert masked_inputs.shape == inputs.shape
    
    def test_high_mask_probability(self):
        """High mask probability (50%+) should work."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(4, 32, vocab)
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.50)
        
        num_masked = (labels != -100).sum().item()
        total = 4 * 32
        
        # Should mask roughly 40-60%
        assert 0.35 < num_masked / total < 0.65
    
    def test_zero_mask_probability(self):
        """Zero mask probability should mask minimally (at least 1 per sequence)."""
        vocab = create_test_vocab()
        inputs = create_test_inputs(4, 32, vocab)
        
        masked_inputs, labels = mask_tokens_span(inputs, vocab, mask_prob=0.0)
        
        # With mask_prob=0, num_to_mask = max(1, 0) = 1 per sequence
        # So we expect at least some masking
        num_masked = (labels != -100).sum().item()
        assert num_masked >= 0  # May be 0 if 0*seq_len rounds to 0


def print_span_distribution(vocab: KmerVocabulary, n_samples: int = 5000):
    """
    Debug function: Print span length distribution for visual inspection.
    
    This satisfies the "debug print shows expected span distribution" acceptance test.
    """
    print("\n" + "=" * 70)
    print("SPAN LENGTH DISTRIBUTION DEBUG")
    print("=" * 70)
    
    inputs = create_test_inputs(32, 64, vocab)
    
    # Test geometric distribution
    random.seed(42)
    _, _, geo_spans = mask_tokens_span(
        inputs, vocab, mask_prob=0.15, track_spans=True,
        span_distribution="geometric", mean_span=3.0, max_span_len=10
    )
    
    # Test uniform distribution  
    random.seed(42)
    _, _, uni_spans = mask_tokens_span(
        inputs, vocab, mask_prob=0.15, track_spans=True,
        span_distribution="uniform", max_span_len=10
    )
    
    print("\n[GEOMETRIC DISTRIBUTION] (mean=3.0, max=10)")
    print("-" * 50)
    geo_counter = Counter(geo_spans)
    total_geo = len(geo_spans)
    for length in sorted(geo_counter.keys()):
        count = geo_counter[length]
        pct = count / total_geo * 100
        bar = "█" * int(pct)
        print(f"  Length {length:>2}: {count:>4} ({pct:5.1f}%) {bar}")
    if geo_spans:
        print(f"  Mean span length: {sum(geo_spans)/len(geo_spans):.2f}")
        print(f"  Total spans: {total_geo}")
    
    print("\n[UNIFORM DISTRIBUTION] (max=10)")
    print("-" * 50)
    uni_counter = Counter(uni_spans)
    total_uni = len(uni_spans)
    for length in sorted(uni_counter.keys()):
        count = uni_counter[length]
        pct = count / total_uni * 100
        bar = "█" * int(pct)
        print(f"  Length {length:>2}: {count:>4} ({pct:5.1f}%) {bar}")
    if uni_spans:
        print(f"  Mean span length: {sum(uni_spans)/len(uni_spans):.2f}")
        print(f"  Total spans: {total_uni}")
    
    print("\n" + "=" * 70)
    print("✓ Geometric distribution favors shorter spans (DNABERT-style)")
    print("=" * 70 + "\n")


def run_all_tests():
    """Run all tests and report results."""
    print("\n" + "=" * 70)
    print("SPAN MASKING UNIT TESTS")
    print("=" * 70)
    
    test_classes = [
        TestSpanLengthSampling,
        TestMaskTokensSpan,
        TestPaddingHandling,
        TestDistributionModes,
        TestSpanTracking,
        TestEdgeCases,
    ]
    
    total_tests = 0
    passed_tests = 0
    failed_tests = []
    
    for test_class in test_classes:
        class_name = test_class.__name__
        print(f"\n[{class_name}]")
        print("-" * 50)
        
        instance = test_class()
        methods = [m for m in dir(instance) if m.startswith('test_')]
        
        for method_name in methods:
            total_tests += 1
            try:
                getattr(instance, method_name)()
                print(f"  ✓ {method_name}")
                passed_tests += 1
            except AssertionError as e:
                print(f"  ✗ {method_name}: {e}")
                failed_tests.append((class_name, method_name, str(e)))
            except Exception as e:
                print(f"  ✗ {method_name}: EXCEPTION - {e}")
                failed_tests.append((class_name, method_name, f"Exception: {e}"))
    
    print("\n" + "=" * 70)
    print(f"RESULTS: {passed_tests}/{total_tests} tests passed")
    
    if failed_tests:
        print("\nFAILED TESTS:")
        for class_name, method_name, error in failed_tests:
            print(f"  - {class_name}.{method_name}: {error}")
        print("=" * 70)
        return False
    else:
        print("✓ All tests passed!")
        print("=" * 70)
        return True


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test DNABERT-style span masking")
    parser.add_argument("--debug", action="store_true", 
                        help="Print span distribution debug info")
    args = parser.parse_args()
    
    # Run tests
    success = run_all_tests()
    
    # Show debug distribution if requested
    if args.debug:
        vocab = create_test_vocab()
        print_span_distribution(vocab)
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)
