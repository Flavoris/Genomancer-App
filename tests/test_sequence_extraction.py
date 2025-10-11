#!/usr/bin/env python3
"""
Test suite for fragment sequence extraction features.
Tests the new functionality for extracting DNA sequences from fragments.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fragment_calculator import (
    compute_fragments_with_sequences,
    Fragment,
    EndInfo,
    slice_circular,
    elide_sequence,
    calculate_overhang_length
)


def test_linear_blunt_single_cut():
    """Test linear DNA with a single blunt-end cut (EcoRV)."""
    print("=" * 70)
    print("TEST: Linear DNA with single blunt cut (EcoRV)")
    print("=" * 70)
    
    # EcoRV cuts GATATC at position 3 (blunt)
    dna_sequence = "AAAAAGATATCGGGG"  # 15 bp
    cut_positions = [9]  # After GATATC
    
    cut_metadata = {
        9: [{
            'enzyme': 'EcoRV',
            'site': 'GATATC',
            'cut_index': 3,
            'overhang_type': 'Blunt'
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=False,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    # Fragment 1: AAAAAGATA
    frag1 = fragments[0]
    assert frag1.length == 9, f"Fragment 1 length: expected 9, got {frag1.length}"
    assert frag1.sequence == "AAAAAGATA", f"Fragment 1 seq: {frag1.sequence}"
    assert frag1.enzymes_at_ends[0] is None, "Fragment 1 left end should be None"
    assert frag1.enzymes_at_ends[1] is not None, "Fragment 1 right end should have enzyme"
    assert frag1.enzymes_at_ends[1].overhang_type == "Blunt", "Should be blunt"
    
    # Fragment 2: TCGGGG
    frag2 = fragments[1]
    assert frag2.length == 6, f"Fragment 2 length: expected 6, got {frag2.length}"
    assert frag2.sequence == "TCGGGG", f"Fragment 2 seq: {frag2.sequence}"
    assert frag2.enzymes_at_ends[0] is not None, "Fragment 2 left end should have enzyme"
    assert frag2.enzymes_at_ends[1] is None, "Fragment 2 right end should be None"
    
    # Verify total length
    total = sum(f.length for f in fragments)
    assert total == len(dna_sequence), f"Total length mismatch: {total} vs {len(dna_sequence)}"
    
    # Verify sequences concatenate correctly
    full_seq = frag1.sequence + frag2.sequence
    assert full_seq == dna_sequence, "Sequences don't concatenate correctly"
    
    print(f"✓ Fragment 1: {frag1.length} bp - {frag1.sequence}")
    print(f"✓ Fragment 2: {frag2.length} bp - {frag2.sequence}")
    print("✓ Test passed!")
    print()


def test_linear_5prime_overhang():
    """Test linear DNA with 5' overhang enzyme (EcoRI)."""
    print("=" * 70)
    print("TEST: Linear DNA with 5' overhang cut (EcoRI)")
    print("=" * 70)
    
    # EcoRI cuts G^AATTC (5' overhang of 4 bp)
    dna_sequence = "AAAAGAATTCGGGG"  # 14 bp
    cut_positions = [5]  # After first G in GAATTC
    
    cut_metadata = {
        5: [{
            'enzyme': 'EcoRI',
            'site': 'GAATTC',
            'cut_index': 1,
            'overhang_type': "5' overhang"
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=False,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    # Check overhang calculation
    frag1_right = fragments[0].enzymes_at_ends[1]
    assert frag1_right.overhang_len == 4, f"Expected 4 bp overhang, got {frag1_right.overhang_len}"
    
    print(f"✓ Fragment 1: {fragments[0].length} bp - {fragments[0].sequence}")
    print(f"  Right end: {frag1_right.enzyme}, overhang={frag1_right.overhang_len} bp")
    print(f"✓ Fragment 2: {fragments[1].length} bp - {fragments[1].sequence}")
    print("✓ Test passed!")
    print()


def test_linear_3prime_overhang():
    """Test linear DNA with 3' overhang enzyme (PstI)."""
    print("=" * 70)
    print("TEST: Linear DNA with 3' overhang cut (PstI)")
    print("=" * 70)
    
    # PstI cuts CTGCA^G (3' overhang of 4 bp)
    dna_sequence = "AAAACTGCAGGGGG"  # 14 bp
    cut_positions = [9]  # After CTGCA
    
    cut_metadata = {
        9: [{
            'enzyme': 'PstI',
            'site': 'CTGCAG',
            'cut_index': 5,
            'overhang_type': "3' overhang"
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=False,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    # Check overhang calculation
    frag1_right = fragments[0].enzymes_at_ends[1]
    assert frag1_right.overhang_len == 4, f"Expected 4 bp overhang, got {frag1_right.overhang_len}"
    assert frag1_right.overhang_type == "3' overhang", "Should be 3' overhang"
    
    print(f"✓ Fragment 1: {fragments[0].length} bp - {fragments[0].sequence}")
    print(f"  Right end: {frag1_right.enzyme}, overhang={frag1_right.overhang_len} bp (3' overhang)")
    print(f"✓ Fragment 2: {fragments[1].length} bp - {fragments[1].sequence}")
    print("✓ Test passed!")
    print()


def test_circular_no_cuts():
    """Test circular DNA with no cuts (intact plasmid)."""
    print("=" * 70)
    print("TEST: Circular DNA with no cuts")
    print("=" * 70)
    
    dna_sequence = "ATCGATCGATCGATCG"  # 16 bp plasmid
    cut_positions = []
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=True,
        cut_metadata={}
    )
    
    assert len(fragments) == 1, f"Expected 1 fragment, got {len(fragments)}"
    
    frag = fragments[0]
    assert frag.length == len(dna_sequence), "Length should equal full sequence"
    assert frag.sequence == dna_sequence, "Sequence should be complete"
    assert frag.wraps == True, "Should be marked as wrapping"
    assert frag.enzymes_at_ends[0] is None, "No left end enzyme"
    assert frag.enzymes_at_ends[1] is None, "No right end enzyme"
    
    print(f"✓ Intact plasmid: {frag.length} bp, wraps={frag.wraps}")
    print("✓ Test passed!")
    print()


def test_circular_single_cut_linearizes():
    """Test circular DNA with single cut in linearize mode."""
    print("=" * 70)
    print("TEST: Circular DNA with single cut (linearizes)")
    print("=" * 70)
    
    dna_sequence = "ATCGATCGATCGATCG"  # 16 bp
    cut_positions = [8]  # Cut in middle
    
    cut_metadata = {
        8: [{
            'enzyme': 'TestEnz',
            'site': 'ATCG',
            'cut_index': 2,
            'overhang_type': 'Blunt'
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=True,
        circular_single_cut_linearizes=True,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments in linearize mode, got {len(fragments)}"
    
    # Verify sequences concatenate to full plasmid
    combined = fragments[0].sequence + fragments[1].sequence
    assert combined == dna_sequence, "Linearized fragments should reconstruct plasmid"
    
    print(f"✓ Fragment 1: {fragments[0].length} bp - {fragments[0].sequence}")
    print(f"✓ Fragment 2: {fragments[1].length} bp - {fragments[1].sequence}")
    print(f"✓ Combined: {combined}")
    print("✓ Test passed!")
    print()


def test_circular_two_cuts_with_wrap():
    """Test circular DNA with two cuts, creating wrap-around fragment."""
    print("=" * 70)
    print("TEST: Circular DNA with two cuts (wrap-around)")
    print("=" * 70)
    
    dna_sequence = "AAAAGGGGCCCCTTTT"  # 16 bp
    cut_positions = [4, 12]  # Two cuts
    
    cut_metadata = {
        4: [{'enzyme': 'Enz1', 'site': 'AAAA', 'cut_index': 2, 'overhang_type': 'Blunt'}],
        12: [{'enzyme': 'Enz2', 'site': 'CCCC', 'cut_index': 2, 'overhang_type': 'Blunt'}]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=True,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    # One fragment should wrap
    wrap_frag = [f for f in fragments if f.wraps][0]
    assert wrap_frag.wraps == True, "Should have one wrap-around fragment"
    
    # Verify wrap-around sequence
    expected_wrap = dna_sequence[12:] + dna_sequence[0:4]  # "TTTT" + "AAAA"
    assert wrap_frag.sequence == expected_wrap, f"Wrap sequence: {wrap_frag.sequence} vs {expected_wrap}"
    
    # Verify all fragments sum to total
    total = sum(f.length for f in fragments)
    assert total == len(dna_sequence), "Fragments should sum to total length"
    
    print(f"✓ Fragment 1: {fragments[0].length} bp - {fragments[0].sequence} (wraps={fragments[0].wraps})")
    print(f"✓ Fragment 2: {fragments[1].length} bp - {fragments[1].sequence} (wraps={fragments[1].wraps})")
    print("✓ Test passed!")
    print()


def test_helper_functions():
    """Test helper functions like slice_circular and elide_sequence."""
    print("=" * 70)
    print("TEST: Helper functions")
    print("=" * 70)
    
    # Test slice_circular
    seq = "ABCDEFGH"
    
    # Normal slice
    assert slice_circular(seq, 0, 4) == "ABCD", "Normal slice failed"
    
    # Wrap-around slice
    assert slice_circular(seq, 6, 2) == "GH" + "AB", "Wrap-around slice failed"
    
    # Full circle
    assert slice_circular(seq, 0, 0) == seq, "Full circle slice failed"
    
    print("✓ slice_circular works correctly")
    
    # Test elide_sequence
    long_seq = "ATCGATCGATCGATCGATCGATCG"
    
    # No elision
    assert elide_sequence(long_seq, 0) == long_seq, "No elision failed"
    
    # With elision
    elided = elide_sequence(long_seq, 5)
    assert elided.startswith("ATCGA"), "Elision start failed"
    assert elided.endswith("ATCG"), "Elision end failed"
    assert "..." in elided, "Elision should contain ..."
    
    print(f"✓ elide_sequence works: '{long_seq}' -> '{elided}'")
    
    # Test calculate_overhang_length
    # EcoRI: GAATTC, cut at 1 -> |6 - 2*1| = 4
    assert calculate_overhang_length("GAATTC", 1) == 4, "EcoRI overhang failed"
    
    # Blunt: cut in middle
    assert calculate_overhang_length("GATATC", 3) == 0, "Blunt overhang failed"
    
    # PstI: CTGCAG, cut at 5 -> |6 - 2*5| = 4
    assert calculate_overhang_length("CTGCAG", 5) == 4, "PstI overhang failed"
    
    print("✓ calculate_overhang_length works correctly")
    print("✓ Test passed!")
    print()


def test_iupac_degenerate():
    """Test IUPAC degenerate base recognition."""
    print("=" * 70)
    print("TEST: IUPAC degenerate bases")
    print("=" * 70)
    
    # AccB2I recognizes RGCGCY (R=[AG], Y=[CT])
    # Should match AGCGCC or GGCGCT etc.
    dna_sequence = "AAAAAGCGCCGGGG"  # Contains AGCGCC
    cut_positions = [9]
    
    cut_metadata = {
        9: [{
            'enzyme': 'AccB2I',
            'site': 'RGCGCY',
            'cut_index': 5,
            'overhang_type': "3' overhang"
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=False,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    print(f"✓ IUPAC site RGCGCY matched AGCGCC")
    print(f"✓ Fragment 1: {fragments[0].length} bp")
    print(f"✓ Fragment 2: {fragments[1].length} bp")
    print("✓ Test passed!")
    print()


def test_type_iis_enzyme():
    """Test Type IIS enzyme with offset cutting (BsaI)."""
    print("=" * 70)
    print("TEST: Type IIS enzyme with offset cut (BsaI-like)")
    print("=" * 70)
    
    # BsaI recognizes GGTCTC and cuts 7 bp downstream (outside recognition site)
    # Sequence: AAAAGGTCTCNNNNGGG
    #                      ^--- cut here (position 14)
    dna_sequence = "AAAAGGTCTCAAAAGGG"  # 17 bp
    cut_positions = [14]  # 7 bp after start of GGTCTC
    
    cut_metadata = {
        14: [{
            'enzyme': 'BsaI',
            'site': 'GGTCTCNNNN',  # Recognition + offset
            'cut_index': 7,
            'overhang_type': "5' overhang"
        }]
    }
    
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_positions,
        circular=False,
        cut_metadata=cut_metadata
    )
    
    assert len(fragments) == 2, f"Expected 2 fragments, got {len(fragments)}"
    
    # Fragment 1 should be: AAAAGGTCTCAAAA (14 bp)
    assert fragments[0].length == 14, f"Fragment 1 length: {fragments[0].length}"
    assert fragments[0].sequence == "AAAAGGTCTCAAAA", f"Fragment 1: {fragments[0].sequence}"
    
    # Fragment 2 should be: GGG (3 bp)
    assert fragments[1].length == 3, f"Fragment 2 length: {fragments[1].length}"
    assert fragments[1].sequence == "GGG", f"Fragment 2: {fragments[1].sequence}"
    
    print(f"✓ Type IIS enzyme cut outside recognition site")
    print(f"✓ Fragment 1: {fragments[0].length} bp - {fragments[0].sequence}")
    print(f"✓ Fragment 2: {fragments[1].length} bp - {fragments[1].sequence}")
    print("✓ Test passed!")
    print()


def run_all_tests():
    """Run all test functions."""
    print("\n" + "=" * 70)
    print("RUNNING FRAGMENT SEQUENCE EXTRACTION TEST SUITE")
    print("=" * 70 + "\n")
    
    tests = [
        test_linear_blunt_single_cut,
        test_linear_5prime_overhang,
        test_linear_3prime_overhang,
        test_circular_no_cuts,
        test_circular_single_cut_linearizes,
        test_circular_two_cuts_with_wrap,
        test_helper_functions,
        test_iupac_degenerate,
        test_type_iis_enzyme
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"✗ FAILED: {test_func.__name__}")
            print(f"  Error: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ ERROR in {test_func.__name__}: {e}")
            failed += 1
    
    print("\n" + "=" * 70)
    print(f"TEST SUMMARY: {passed} passed, {failed} failed")
    print("=" * 70 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

