#!/usr/bin/env python3
"""
Test suite for overhang correctness across all outputs.
Tests that the centralized compute_end_metadata function produces correct results.
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fragment_calculator import compute_end_metadata


class TestOverhangCorrectness:
    """Test that overhang metadata is computed correctly for common enzymes."""
    
    def test_ecori_5prime_overhang(self):
        """
        EcoRI cuts G^AATTC to produce 5' overhang AATT (4 bp).
        Recognition site: GAATTC
        Cut index: 1 (cuts after G)
        """
        dna = "AAAGAATTCGGG"
        cut_pos = 4  # Position after the G in GAATTC
        
        # Left end (fragment ending at this cut)
        result_left = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=True,
            circular=False
        )
        
        assert result_left['overhang_len'] == 4, "EcoRI should produce 4 bp overhang"
        assert result_left['sticky_seq'] == "AATT", f"EcoRI left end should be AATT, got {result_left['sticky_seq']}"
        assert result_left['end_bases'] == "AATT", "end_bases should match sticky_seq"
        assert result_left['polarity'] == "left"
        
        # Right end (fragment starting at this cut)
        result_right = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=False,
            circular=False
        )
        
        assert result_right['overhang_len'] == 4
        assert result_right['sticky_seq'] == "AATT", f"EcoRI right end should be AATT, got {result_right['sticky_seq']}"
        assert result_right['polarity'] == "right"
    
    def test_hindiii_5prime_overhang(self):
        """
        HindIII cuts A^AGCTT to produce 5' overhang AGCT (4 bp).
        Recognition site: AAGCTT
        Cut index: 1
        """
        dna = "TTTAAGCTTAAA"
        cut_pos = 4  # Position after first A in AAGCTT
        
        result = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="AAGCTT",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=True,
            circular=False
        )
        
        assert result['overhang_len'] == 4, "HindIII should produce 4 bp overhang"
        assert result['sticky_seq'] == "AGCT", f"HindIII should produce AGCT overhang, got {result['sticky_seq']}"
    
    def test_psti_3prime_overhang(self):
        """
        PstI cuts CTGCA^G to produce 3' overhang TGCA (4 bp).
        Recognition site: CTGCAG
        Cut index: 5 (cuts after CTGCA)
        """
        dna = "AAACTGCAGTTT"
        cut_pos = 8  # Position after CTGCA
        
        result = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="CTGCAG",
            cut_index=5,
            overhang_type="3' overhang",
            is_left_end=False,
            circular=False
        )
        
        assert result['overhang_len'] == 4, "PstI should produce 4 bp overhang"
        assert result['sticky_seq'] == "TGCA", f"PstI should produce TGCA overhang, got {result['sticky_seq']}"
    
    def test_smai_blunt(self):
        """
        SmaI cuts CCC^GGG to produce blunt end (0 bp overhang).
        Recognition site: CCCGGG
        Cut index: 3
        """
        dna = "AAACCCGGGAAA"
        cut_pos = 6  # Position after CCC
        
        result = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="CCCGGG",
            cut_index=3,
            overhang_type="Blunt",
            is_left_end=True,
            circular=False
        )
        
        assert result['overhang_len'] == 0, "SmaI should produce blunt end (0 bp)"
        assert result['sticky_seq'] == "", "Blunt end should have empty sticky_seq"
        assert result['end_bases'] == "", "Blunt end should have empty end_bases"
    
    def test_bamhi_5prime_overhang(self):
        """
        BamHI cuts G^GATCC to produce 5' overhang GATC (4 bp).
        Recognition site: GGATCC
        Cut index: 1
        """
        dna = "AAAGGATCCGGG"
        cut_pos = 4  # Position after first G in GGATCC
        
        result = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos,
            recognition_site="GGATCC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=True,
            circular=False
        )
        
        assert result['overhang_len'] == 4, "BamHI should produce 4 bp overhang"
        assert result['sticky_seq'] == "GATC", f"BamHI should produce GATC overhang, got {result['sticky_seq']}"
    
    def test_circular_wrapping(self):
        """Test that overhang extraction works correctly with circular DNA wrapping."""
        # Circular DNA where the recognition site wraps around
        dna = "AATTCGGGGG"  # Imagine this wraps: ...GGGGG+G^AATTC...
        
        # If we have a cut at position 1 (after wrapping G from end)
        result = compute_end_metadata(
            dna=dna,
            cut_pos=1,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=False,
            circular=True
        )
        
        assert result['overhang_len'] == 4, "EcoRI overhang should be 4 bp even in circular"
        # The sticky sequence might wrap, but length should be correct
        assert len(result['sticky_seq']) == 4, "Sticky sequence should be 4 bases"


class TestOverhangConsistency:
    """Test that overhang values are consistent across different use cases."""
    
    def test_left_and_right_same_enzyme(self):
        """
        When the same enzyme creates both ends of a fragment,
        both ends should have the same overhang length.
        """
        dna = "AAAGAATTCGGGGGAATTCAAA"
        cut_pos_1 = 4  # First EcoRI site at position 3, cuts after G at position 4
        cut_pos_2 = 14  # Second EcoRI site at position 13, cuts after G at position 14
        
        left_end = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos_1,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=False,
            circular=False
        )
        
        right_end = compute_end_metadata(
            dna=dna,
            cut_pos=cut_pos_2,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=True,
            circular=False
        )
        
        assert left_end['overhang_len'] == right_end['overhang_len'], \
            "Same enzyme should produce same overhang length on both ends"
        assert left_end['sticky_seq'] == right_end['sticky_seq'], \
            "EcoRI should produce compatible sticky ends"
    
    def test_polarity_assignment(self):
        """Test that polarity is correctly assigned based on is_left_end."""
        dna = "AAAGAATTCGGG"
        
        left = compute_end_metadata(
            dna=dna,
            cut_pos=4,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=True,
            circular=False
        )
        
        right = compute_end_metadata(
            dna=dna,
            cut_pos=4,
            recognition_site="GAATTC",
            cut_index=1,
            overhang_type="5' overhang",
            is_left_end=False,
            circular=False
        )
        
        assert left['polarity'] == "left"
        assert right['polarity'] == "right"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

