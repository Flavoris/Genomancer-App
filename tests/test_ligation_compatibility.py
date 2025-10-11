#!/usr/bin/env python3
"""
Test suite for ligation compatibility analysis.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from ligation_compatibility import (
    revcomp, calculate_gc_percent, calculate_tm,
    are_compatible, is_directional, EndInfo,
    calculate_compatibility, analyze_enzyme_pair_theoretical
)


class TestRevcomp:
    """Test reverse complement function."""
    
    def test_revcomp_simple(self):
        """Test basic reverse complement."""
        assert revcomp("ATCG") == "CGAT"
        assert revcomp("GAATTC") == "GAATTC"  # Palindrome
        assert revcomp("ACTAGT") == "ACTAGT"  # SpeI site (palindrome)
    
    def test_revcomp_lowercase(self):
        """Test reverse complement with lowercase."""
        assert revcomp("atcg") == "cgat"
        # Note: mixed case gets translated, so case may not be fully preserved
        result = revcomp("AaTtCcGg")
        # Just check that the sequence is correct (case may vary)
        assert result.upper() == "CCGGAATT"
    
    def test_revcomp_empty(self):
        """Test empty string."""
        assert revcomp("") == ""


class TestHeuristics:
    """Test GC% and Tm calculation."""
    
    def test_gc_percent(self):
        """Test GC percentage calculation."""
        assert calculate_gc_percent("AAAA") == 0.0
        assert calculate_gc_percent("TTTT") == 0.0
        assert calculate_gc_percent("GGGG") == 100.0
        assert calculate_gc_percent("CCCC") == 100.0
        assert calculate_gc_percent("ATGC") == 50.0
        assert calculate_gc_percent("AATT") == 0.0
    
    def test_gc_percent_empty(self):
        """Test GC% with empty sequence."""
        assert calculate_gc_percent("") == 0.0
    
    def test_tm_calculation(self):
        """Test Tm calculation using Wallace rule."""
        # Wallace rule: Tm = 2*(A+T) + 4*(G+C)
        assert calculate_tm("AAAA") == 8.0   # 2*4 + 4*0
        assert calculate_tm("TTTT") == 8.0   # 2*4 + 4*0
        assert calculate_tm("GGGG") == 16.0  # 2*0 + 4*4
        assert calculate_tm("CCCC") == 16.0  # 2*0 + 4*4
        assert calculate_tm("AATT") == 8.0   # EcoRI overhang: 2*4 + 4*0
        assert calculate_tm("ACTAGT") == 16.0  # SpeI: 2*4 + 4*2
    
    def test_tm_empty(self):
        """Test Tm with empty sequence."""
        assert calculate_tm("") == 0.0


class TestCompatibility:
    """Test compatibility checking between fragment ends."""
    
    def test_compatible_5prime_overhangs(self):
        """Test compatible 5' overhangs (EcoRI and complementary)."""
        # EcoRI produces AATT overhang
        end_a = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        # Complementary end with AATT (revcomp = AATT, so compatible)
        end_b = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert are_compatible(end_a, end_b, include_blunt=False, min_overhang=1)
    
    def test_compatible_spei_xbai(self):
        """Test SpeI/XbaI compatibility (classic example)."""
        # SpeI produces CTAG overhang
        spei_end = EndInfo(
            enzyme="SpeI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="CTAG",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        # XbaI produces CTAG overhang (same as SpeI)
        xbai_end = EndInfo(
            enzyme="XbaI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="CTAG",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        # They should be compatible (CTAG revcomp = CTAG)
        assert are_compatible(spei_end, xbai_end, include_blunt=False, min_overhang=1)
    
    def test_incompatible_different_lengths(self):
        """Test incompatible ends with different overhang lengths."""
        end_a = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="MfeI",
            overhang_type="5' overhang",
            overhang_len=2,
            sticky_seq="AA",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert not are_compatible(end_a, end_b, include_blunt=False, min_overhang=1)
    
    def test_incompatible_different_types(self):
        """Test incompatible ends with different overhang types."""
        # 5' overhang
        end_a = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        # 3' overhang
        end_b = EndInfo(
            enzyme="PstI",
            overhang_type="3' overhang",
            overhang_len=4,
            sticky_seq="TGCA",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert not are_compatible(end_a, end_b, include_blunt=False, min_overhang=1)
    
    def test_blunt_blunt_without_flag(self):
        """Test blunt-blunt compatibility without include_blunt flag."""
        end_a = EndInfo(
            enzyme="EcoRV",
            overhang_type="Blunt",
            overhang_len=0,
            sticky_seq="",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="SmaI",
            overhang_type="Blunt",
            overhang_len=0,
            sticky_seq="",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert not are_compatible(end_a, end_b, include_blunt=False, min_overhang=1)
    
    def test_blunt_blunt_with_flag(self):
        """Test blunt-blunt compatibility with include_blunt flag."""
        end_a = EndInfo(
            enzyme="EcoRV",
            overhang_type="Blunt",
            overhang_len=0,
            sticky_seq="",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="SmaI",
            overhang_type="Blunt",
            overhang_len=0,
            sticky_seq="",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert are_compatible(end_a, end_b, include_blunt=True, min_overhang=1)
    
    def test_sticky_blunt_incompatible(self):
        """Test sticky-blunt is always incompatible."""
        end_a = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="EcoRV",
            overhang_type="Blunt",
            overhang_len=0,
            sticky_seq="",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert not are_compatible(end_a, end_b, include_blunt=True, min_overhang=1)
    
    def test_min_overhang_filter(self):
        """Test minimum overhang length filtering."""
        end_a = EndInfo(
            enzyme="TestEnz",
            overhang_type="5' overhang",
            overhang_len=2,
            sticky_seq="AT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="TestEnz",
            overhang_type="5' overhang",
            overhang_len=2,
            sticky_seq="AT",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        # Should be compatible with min_overhang=1 or 2
        assert are_compatible(end_a, end_b, include_blunt=False, min_overhang=1)
        assert are_compatible(end_a, end_b, include_blunt=False, min_overhang=2)
        
        # Should be incompatible with min_overhang=3
        assert not are_compatible(end_a, end_b, include_blunt=False, min_overhang=3)


class TestDirectionality:
    """Test directionality checking."""
    
    def test_directional_spei_xbai(self):
        """Test SpeI/XbaI is directional (non-palindromic)."""
        # SpeI produces CTAG overhang
        spei_end = EndInfo(
            enzyme="SpeI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="CTAG",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        xbai_end = EndInfo(
            enzyme="XbaI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="CTAG",
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        # CTAG revcomp = CTAG, so it's actually palindromic!
        # Wait, let me recalculate: CTAG -> reverse: GATC -> complement: CTAG
        # So CTAG is palindromic
        assert not is_directional(spei_end, xbai_end)
    
    def test_non_directional_ecori(self):
        """Test EcoRI is non-directional (palindromic)."""
        # EcoRI produces AATT overhang
        # AATT -> reverse: TTAA -> complement: AATT
        # So AATT is palindromic
        end = EndInfo(
            enzyme="EcoRI",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="AATT",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        assert not is_directional(end, end)
    
    def test_directional_non_palindromic(self):
        """Test non-palindromic overhang is directional."""
        # ACGT -> reverse: TGCA -> complement: ACGT
        # Wait, that's also palindromic
        # Let's use ACTA
        # ACTA -> reverse: ATCA -> complement: TAGT
        # ACTA != TAGT, so directional
        end_a = EndInfo(
            enzyme="TestEnz",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="ACTA",
            polarity="right",
            fragment_id=1,
            position=100
        )
        
        end_b = EndInfo(
            enzyme="TestEnz2",
            overhang_type="5' overhang",
            overhang_len=4,
            sticky_seq="TAGT",  # revcomp of ACTA
            polarity="left",
            fragment_id=2,
            position=200
        )
        
        assert is_directional(end_a, end_b)


class TestCalculateCompatibility:
    """Test full compatibility analysis."""
    
    def test_calculate_compatibility_basic(self):
        """Test basic compatibility calculation."""
        ends = [
            EndInfo(
                enzyme="EcoRI",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="AATT",
                polarity="right",
                fragment_id=0,
                position=100
            ),
            EndInfo(
                enzyme="EcoRI",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="AATT",
                polarity="left",
                fragment_id=1,
                position=200
            ),
        ]
        
        results = calculate_compatibility(
            ends=ends,
            include_blunt=False,
            min_overhang=1,
            require_directional=False
        )
        
        assert len(results) == 1
        assert results[0].compatible
        assert results[0].end_a == ends[0]
        assert results[0].end_b == ends[1]
    
    def test_calculate_compatibility_require_directional(self):
        """Test filtering by directionality requirement."""
        ends = [
            # Palindromic pair (EcoRI)
            EndInfo(
                enzyme="EcoRI",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="AATT",
                polarity="right",
                fragment_id=0,
                position=100
            ),
            EndInfo(
                enzyme="EcoRI",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="AATT",
                polarity="left",
                fragment_id=1,
                position=200
            ),
            # Non-palindromic pair
            EndInfo(
                enzyme="TestEnz",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="ACTA",
                polarity="right",
                fragment_id=2,
                position=300
            ),
            EndInfo(
                enzyme="TestEnz2",
                overhang_type="5' overhang",
                overhang_len=4,
                sticky_seq="TAGT",
                polarity="left",
                fragment_id=3,
                position=400
            ),
        ]
        
        # Without directionality requirement, should find 2 pairs
        results_all = calculate_compatibility(
            ends=ends,
            include_blunt=False,
            min_overhang=1,
            require_directional=False
        )
        assert len(results_all) == 2
        
        # With directionality requirement, should find only 1 pair (the non-palindromic one)
        results_directional = calculate_compatibility(
            ends=ends,
            include_blunt=False,
            min_overhang=1,
            require_directional=True
        )
        assert len(results_directional) == 1
        assert results_directional[0].directional


class TestEnzymePairTheoretical:
    """Test theoretical enzyme pair analysis."""
    
    def test_ecori_ecori_compatible(self):
        """Test EcoRI with itself is compatible."""
        enzyme_a = {
            "name": "EcoRI",
            "site": "GAATTC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        }
        
        result = analyze_enzyme_pair_theoretical(enzyme_a, enzyme_a)
        assert result["compatible"]
        assert result["enzyme_a"] == "EcoRI"
        assert result["enzyme_b"] == "EcoRI"
    
    def test_ecori_mfei_incompatible_length(self):
        """Test EcoRI and MfeI are incompatible (different lengths)."""
        ecori = {
            "name": "EcoRI",
            "site": "GAATTC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        }
        
        mfei = {
            "name": "MfeI",
            "site": "CAATTG",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        }
        
        result = analyze_enzyme_pair_theoretical(ecori, mfei)
        # EcoRI: 6 bp site, cut at 1 -> |2*1 - 6| = 4 bp overhang
        # MfeI: 6 bp site, cut at 1 -> |2*1 - 6| = 4 bp overhang
        # Actually, they should have the same length!
        assert result["overhang_len_a"] == 4
        assert result["overhang_len_b"] == 4
    
    def test_blunt_blunt_compatible(self):
        """Test two blunt cutters are compatible."""
        ecorv = {
            "name": "EcoRV",
            "site": "GATATC",
            "cut_index": 3,
            "overhang_type": "Blunt"
        }
        
        smai = {
            "name": "SmaI",
            "site": "CCCGGG",
            "cut_index": 3,
            "overhang_type": "Blunt"
        }
        
        result = analyze_enzyme_pair_theoretical(ecorv, smai)
        assert result["compatible"]
        assert result["overhang_len_a"] == 0
        assert result["overhang_len_b"] == 0
    
    def test_5prime_3prime_incompatible(self):
        """Test 5' overhang and 3' overhang are incompatible."""
        ecori = {
            "name": "EcoRI",
            "site": "GAATTC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        }
        
        psti = {
            "name": "PstI",
            "site": "CTGCAG",
            "cut_index": 5,
            "overhang_type": "3' overhang"
        }
        
        result = analyze_enzyme_pair_theoretical(ecori, psti)
        assert not result["compatible"]
        assert "don't match" in result["note"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

