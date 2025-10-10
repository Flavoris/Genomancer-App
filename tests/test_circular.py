#!/usr/bin/env python3
"""
Unit tests for circular DNA digestion functionality.
Tests the compute_fragments function with circular topology.
"""

import pytest
import sys
import os

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fragment_calculator import compute_fragments, validate_fragment_total


class TestCircularNoChunks:
    """Test circular mode with no cuts."""
    
    def test_circular_no_cuts_returns_whole_circle(self):
        """With no cuts, circular DNA should return one intact circle."""
        seq_len = 10
        cuts = []
        frags = compute_fragments(cuts, seq_len, circular=True)
        
        assert len(frags) == 1
        assert frags[0]["length"] == 10
        assert frags[0]["wraps"] is True
        assert frags[0]["start"] == 0
        assert frags[0]["end"] == 0
        assert validate_fragment_total(frags, seq_len)
    
    def test_circular_no_cuts_large_plasmid(self):
        """Test with a larger plasmid size."""
        seq_len = 5000
        cuts = []
        frags = compute_fragments(cuts, seq_len, circular=True)
        
        assert len(frags) == 1
        assert frags[0]["length"] == 5000
        assert frags[0]["wraps"] is True
        assert validate_fragment_total(frags, seq_len)


class TestCircularSingleCut:
    """Test circular mode with a single cut."""
    
    def test_circular_one_cut_default_keeps_one_fragment(self):
        """By default, one cut in circular mode keeps one fragment (intact circle)."""
        seq_len = 10
        cuts = [3]
        frags = compute_fragments(cuts, seq_len, circular=True, circular_single_cut_linearizes=False)
        
        assert len(frags) == 1
        assert frags[0]["length"] == 10
        assert frags[0]["wraps"] is True
        assert validate_fragment_total(frags, seq_len)
    
    def test_circular_one_cut_linearizes_when_flag_true(self):
        """With flag set, one cut in circular mode yields two fragments."""
        seq_len = 10
        cuts = [3]
        frags = compute_fragments(cuts, seq_len, circular=True, circular_single_cut_linearizes=True)
        
        assert len(frags) == 2
        lengths = sorted(f["length"] for f in frags)
        assert lengths == [3, 7]
        assert all(f["wraps"] is False for f in frags)
        assert validate_fragment_total(frags, seq_len)
    
    def test_circular_one_cut_linearizes_different_positions(self):
        """Test single cut linearization at various positions."""
        seq_len = 20
        
        for cut_pos in [0, 5, 10, 15, 19]:
            frags = compute_fragments([cut_pos], seq_len, circular=True, circular_single_cut_linearizes=True)
            assert len(frags) == 2
            total_length = sum(f["length"] for f in frags)
            assert total_length == seq_len


class TestCircularTwoCuts:
    """Test circular mode with two cuts."""
    
    def test_circular_two_cuts_wrap_fragment_present(self):
        """Two cuts in circular mode should produce two fragments, one wrapping."""
        L = 12
        cuts = [3, 9]
        frags = compute_fragments(cuts, L, circular=True)
        
        assert len(frags) == 2
        lengths = [f["length"] for f in frags]
        assert sorted(lengths) == [6, 6]
        
        # One fragment should wrap
        wrap_count = sum(1 for f in frags if f["wraps"])
        assert wrap_count == 1
        
        # Find the wrapping fragment
        wrap_frag = next(f for f in frags if f["wraps"])
        assert wrap_frag["start"] == 9
        assert wrap_frag["end"] == 3
        assert wrap_frag["end"] < wrap_frag["start"]  # Characteristic of wrap
        
        assert validate_fragment_total(frags, L)
    
    def test_circular_two_cuts_asymmetric(self):
        """Test two cuts with unequal fragments."""
        L = 100
        cuts = [10, 80]
        frags = compute_fragments(cuts, L, circular=True)
        
        assert len(frags) == 2
        lengths = sorted(f["length"] for f in frags)
        assert lengths == [30, 70]
        assert validate_fragment_total(frags, L)
        
        # Verify the wrap fragment
        wrap_frag = next(f for f in frags if f["wraps"])
        assert wrap_frag["length"] == 30  # (100 - 80) + 10


class TestCircularMultipleCuts:
    """Test circular mode with three or more cuts."""
    
    def test_circular_three_cuts(self):
        """Three cuts in circular mode should produce three fragments."""
        L = 30
        cuts = [5, 15, 25]
        frags = compute_fragments(cuts, L, circular=True)
        
        assert len(frags) == 3
        lengths = sorted(f["length"] for f in frags)
        assert lengths == [10, 10, 10]
        
        # Exactly one should wrap
        wrap_count = sum(1 for f in frags if f["wraps"])
        assert wrap_count == 1
        
        assert validate_fragment_total(frags, L)
    
    def test_circular_many_cuts(self):
        """Test with many cuts."""
        L = 100
        cuts = [10, 20, 30, 50, 70, 90]
        frags = compute_fragments(cuts, L, circular=True)
        
        assert len(frags) == 6
        expected_lengths = [10, 10, 20, 20, 20, 20]
        lengths = sorted(f["length"] for f in frags)
        assert lengths == expected_lengths
        
        # Exactly one should wrap
        wrap_count = sum(1 for f in frags if f["wraps"])
        assert wrap_count == 1
        
        assert validate_fragment_total(frags, L)


class TestCircularDuplicateCuts:
    """Test that duplicate cut positions are properly deduplicated."""
    
    def test_multi_enzyme_same_pos_dedup(self):
        """Multiple enzymes cutting at same position should be deduplicated."""
        L = 10
        cuts = [2, 2, 7, 7, 7]  # Duplicates from different enzymes
        frags = compute_fragments(cuts, L, circular=True)
        
        # Should only count 2 unique positions
        assert len(frags) == 2
        lengths = sorted(f["length"] for f in frags)
        assert lengths == [5, 5]
        assert validate_fragment_total(frags, L)
    
    def test_duplicate_cuts_with_metadata(self):
        """Test deduplication preserves all enzyme metadata."""
        L = 20
        cuts = [5, 5, 15]  # Two enzymes cut at position 5
        metadata = {
            5: [
                {'enzyme': 'EcoRI', 'site': 'GAATTC', 'cut_index': 1, 'overhang_type': "5' overhang"},
                {'enzyme': 'BamHI', 'site': 'GGATCC', 'cut_index': 1, 'overhang_type': "5' overhang"}
            ],
            15: [
                {'enzyme': 'HindIII', 'site': 'AAGCTT', 'cut_index': 1, 'overhang_type': "5' overhang"}
            ]
        }
        
        frags = compute_fragments(cuts, L, circular=True, cut_metadata=metadata)
        
        assert len(frags) == 2
        
        # Check that position 5 has both enzymes in boundaries
        for frag in frags:
            if frag['boundaries']['left_cut'] and frag['boundaries']['left_cut']['pos'] == 5:
                enzymes = [e['enzyme'] for e in frag['boundaries']['left_cut']['enzymes']]
                assert 'EcoRI' in enzymes
                assert 'BamHI' in enzymes


class TestLinearMode:
    """Test that linear mode still works correctly (no regressions)."""
    
    def test_linear_no_cuts(self):
        """Linear DNA with no cuts returns one fragment."""
        seq_len = 10
        cuts = []
        frags = compute_fragments(cuts, seq_len, circular=False)
        
        assert len(frags) == 1
        assert frags[0]["length"] == 10
        assert frags[0]["wraps"] is False
        assert frags[0]["start"] == 0
        assert frags[0]["end"] == 10
        assert validate_fragment_total(frags, seq_len)
    
    def test_linear_one_cut(self):
        """Linear DNA with one cut yields two fragments."""
        seq_len = 10
        cuts = [3]
        frags = compute_fragments(cuts, seq_len, circular=False)
        
        assert len(frags) == 2
        lengths = [f["length"] for f in frags]
        assert lengths == [3, 7]
        assert all(f["wraps"] is False for f in frags)
        assert validate_fragment_total(frags, seq_len)
    
    def test_linear_two_cuts(self):
        """Linear DNA with two cuts yields three fragments."""
        seq_len = 12
        cuts = [3, 9]
        frags = compute_fragments(cuts, seq_len, circular=False)
        
        assert len(frags) == 3
        lengths = [f["length"] for f in frags]
        assert lengths == [3, 6, 3]
        assert all(f["wraps"] is False for f in frags)
        assert validate_fragment_total(frags, seq_len)
    
    def test_linear_multiple_cuts(self):
        """Linear DNA with multiple cuts."""
        seq_len = 100
        cuts = [10, 30, 60, 90]
        frags = compute_fragments(cuts, seq_len, circular=False)
        
        assert len(frags) == 5
        expected = [10, 20, 30, 30, 10]
        lengths = [f["length"] for f in frags]
        assert lengths == expected
        assert validate_fragment_total(frags, seq_len)


class TestBoundaryAnnotations:
    """Test that boundary annotations are correct."""
    
    def test_linear_boundaries(self):
        """Test boundary annotations in linear mode."""
        seq_len = 20
        cuts = [5, 15]
        metadata = {
            5: [{'enzyme': 'EcoRI', 'site': 'GAATTC', 'cut_index': 1, 'overhang_type': "5' overhang"}],
            15: [{'enzyme': 'BamHI', 'site': 'GGATCC', 'cut_index': 1, 'overhang_type': "5' overhang"}]
        }
        
        frags = compute_fragments(cuts, seq_len, circular=False, cut_metadata=metadata)
        
        assert len(frags) == 3
        
        # First fragment: START -> 5
        assert frags[0]['boundaries']['left_cut'] is None
        assert frags[0]['boundaries']['right_cut']['pos'] == 5
        assert frags[0]['boundaries']['right_cut']['enzymes'][0]['enzyme'] == 'EcoRI'
        
        # Middle fragment: 5 -> 15
        assert frags[1]['boundaries']['left_cut']['pos'] == 5
        assert frags[1]['boundaries']['right_cut']['pos'] == 15
        
        # Last fragment: 15 -> END
        assert frags[2]['boundaries']['left_cut']['pos'] == 15
        assert frags[2]['boundaries']['right_cut'] is None
    
    def test_circular_boundaries_wrap(self):
        """Test boundary annotations for wrap-around fragment in circular mode."""
        seq_len = 20
        cuts = [5, 15]
        metadata = {
            5: [{'enzyme': 'EcoRI', 'site': 'GAATTC', 'cut_index': 1, 'overhang_type': "5' overhang"}],
            15: [{'enzyme': 'BamHI', 'site': 'GGATCC', 'cut_index': 1, 'overhang_type': "5' overhang"}]
        }
        
        frags = compute_fragments(cuts, seq_len, circular=True, cut_metadata=metadata)
        
        assert len(frags) == 2
        
        # Find wrap fragment
        wrap_frag = next(f for f in frags if f["wraps"])
        
        # Should go from 15 to 5 (wrapping around)
        assert wrap_frag['start'] == 15
        assert wrap_frag['end'] == 5
        assert wrap_frag['boundaries']['left_cut']['pos'] == 15
        assert wrap_frag['boundaries']['right_cut']['pos'] == 5


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_zero_length_raises_error(self):
        """Sequence length of 0 should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence length must be positive"):
            compute_fragments([], 0, circular=False)
    
    def test_negative_length_raises_error(self):
        """Negative sequence length should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence length must be positive"):
            compute_fragments([], -5, circular=True)
    
    def test_cut_at_position_zero(self):
        """Test cut at position 0."""
        seq_len = 10
        cuts = [0, 5]
        
        # Linear mode
        frags = compute_fragments(cuts, seq_len, circular=False)
        assert len(frags) == 3
        assert frags[0]["length"] == 0  # Empty fragment at start
        
        # Circular mode
        frags = compute_fragments(cuts, seq_len, circular=True)
        assert len(frags) == 2
        assert validate_fragment_total(frags, seq_len)
    
    def test_cut_at_last_position(self):
        """Test cut at the last valid position."""
        seq_len = 10
        cuts = [9]
        
        frags = compute_fragments(cuts, seq_len, circular=False)
        assert len(frags) == 2
        assert frags[0]["length"] == 9
        assert frags[1]["length"] == 1
    
    def test_modulo_normalization(self):
        """Test that cut positions are normalized with modulo."""
        seq_len = 10
        cuts = [3, 13, 23]  # All equivalent to position 3 mod 10
        
        frags = compute_fragments(cuts, seq_len, circular=True, circular_single_cut_linearizes=True)
        
        # Should treat all as position 3
        assert len(frags) == 2


class TestFragmentOrdering:
    """Test that fragments are ordered consistently."""
    
    def test_linear_fragments_ordered_by_position(self):
        """Linear fragments should be ordered by start position."""
        seq_len = 30
        cuts = [20, 5, 15]  # Unordered input
        
        frags = compute_fragments(cuts, seq_len, circular=False)
        
        # Should be in order: [0-5], [5-15], [15-20], [20-30]
        starts = [f['start'] for f in frags]
        assert starts == sorted(starts)
    
    def test_circular_fragments_ordered_with_wrap_last(self):
        """Circular fragments should have wrap fragment last."""
        seq_len = 30
        cuts = [25, 10, 20]  # Unordered
        
        frags = compute_fragments(cuts, seq_len, circular=True)
        
        # Wrap fragment should be last
        assert frags[-1]["wraps"] is True
        
        # Non-wrap fragments should be ordered
        non_wrap = [f for f in frags if not f["wraps"]]
        starts = [f['start'] for f in non_wrap]
        assert starts == sorted(starts)


class TestValidateFragmentTotal:
    """Test the validation helper function."""
    
    def test_validate_correct_sum(self):
        """Valid fragments should pass validation."""
        frags = [
            {'length': 10},
            {'length': 20},
            {'length': 30}
        ]
        assert validate_fragment_total(frags, 60) is True
    
    def test_validate_incorrect_sum(self):
        """Invalid fragments should fail validation."""
        frags = [
            {'length': 10},
            {'length': 20}
        ]
        assert validate_fragment_total(frags, 100) is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

