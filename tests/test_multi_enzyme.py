#!/usr/bin/env python3
"""
Test cases for the multi-enzyme restriction enzyme simulator.
"""

import sys
import os

# Add parent directory to path to import sim module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim import (
    find_cut_positions_linear,
    merge_cut_positions,
    fragments_linear,
    find_closest_enzyme_names
)


def test_single_enzyme_one_site():
    """Test: One enzyme / one site → two fragments."""
    print("Test 1: Single enzyme with one cut site")
    
    # Create a simple enzyme database
    db = {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1
        }
    }
    
    # Test sequence with one EcoRI site
    seq = "ATCGAATTCGGGATCCAAA"
    
    # Find cut positions
    cuts = find_cut_positions_linear(seq, "EcoRI", db)
    print(f"Cut positions: {cuts}")
    assert cuts == [4], f"Expected [4], got {cuts}"
    
    # Calculate fragments
    fragment_lengths = fragments_linear(len(seq), cuts)
    print(f"Fragment lengths: {fragment_lengths}")
    assert fragment_lengths == [4, 15], f"Expected [4, 15], got {fragment_lengths}"
    
    print("✓ Test 1 passed\n")


def test_two_enzymes_distinct_sites():
    """Test: Two enzymes with distinct sites → correct combined cuts and fragments."""
    print("Test 2: Two enzymes with distinct cut sites")
    
    # Create enzyme database
    db = {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1
        },
        "BamHI": {
            "sequence": "GGATCC",
            "cut_index": 1
        }
    }
    
    # Test sequence with both EcoRI and BamHI sites
    seq = "ATCGAATTCGGGATCCAAA"
    
    # Find cut positions for each enzyme
    eco_cuts = find_cut_positions_linear(seq, "EcoRI", db)
    bam_cuts = find_cut_positions_linear(seq, "BamHI", db)
    
    print(f"EcoRI cuts: {eco_cuts}")
    print(f"BamHI cuts: {bam_cuts}")
    
    assert eco_cuts == [4], f"Expected EcoRI cuts [4], got {eco_cuts}"
    assert bam_cuts == [11], f"Expected BamHI cuts [11], got {bam_cuts}"
    
    # Merge cut positions
    cuts_by_enzyme = {"EcoRI": eco_cuts, "BamHI": bam_cuts}
    combined_cuts = merge_cut_positions(cuts_by_enzyme, len(seq))
    print(f"Combined cuts: {combined_cuts}")
    assert combined_cuts == [4, 11], f"Expected [4, 11], got {combined_cuts}"
    
    # Calculate combined fragments
    fragment_lengths = fragments_linear(len(seq), combined_cuts)
    print(f"Fragment lengths: {fragment_lengths}")
    assert fragment_lengths == [4, 7, 8], f"Expected [4, 7, 8], got {fragment_lengths}"
    
    print("✓ Test 2 passed\n")


def test_two_enzymes_same_coordinate():
    """Test: Two enzymes cutting at the same coordinate → dedup works (count once)."""
    print("Test 3: Two enzymes cutting at the same coordinate")
    
    # Create enzyme database with enzymes that cut at the same position
    db = {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1
        },
        "Custom": {
            "sequence": "GAATTC",  # Same recognition sequence
            "cut_index": 1         # Same cut position
        }
    }
    
    # Test sequence with one site
    seq = "ATCGAATTCGGGATCCAAA"
    
    # Find cut positions for each enzyme
    eco_cuts = find_cut_positions_linear(seq, "EcoRI", db)
    custom_cuts = find_cut_positions_linear(seq, "Custom", db)
    
    print(f"EcoRI cuts: {eco_cuts}")
    print(f"Custom cuts: {custom_cuts}")
    
    assert eco_cuts == [4], f"Expected EcoRI cuts [4], got {eco_cuts}"
    assert custom_cuts == [4], f"Expected Custom cuts [4], got {custom_cuts}"
    
    # Merge cut positions (should deduplicate)
    cuts_by_enzyme = {"EcoRI": eco_cuts, "Custom": custom_cuts}
    combined_cuts = merge_cut_positions(cuts_by_enzyme, len(seq))
    print(f"Combined cuts: {combined_cuts}")
    assert combined_cuts == [4], f"Expected [4], got {combined_cuts}"
    assert len(combined_cuts) == 1, "Should have only one unique cut position"
    
    # Calculate combined fragments
    fragment_lengths = fragments_linear(len(seq), combined_cuts)
    print(f"Fragment lengths: {fragment_lengths}")
    assert fragment_lengths == [4, 15], f"Expected [4, 15], got {fragment_lengths}"
    
    print("✓ Test 3 passed\n")


def test_no_cuts():
    """Test: No cuts → single fragment equal to sequence length."""
    print("Test 4: No enzyme cuts")
    
    # Create enzyme database
    db = {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1
        }
    }
    
    # Test sequence with no cut sites
    seq = "ATCGATCGATCGATCG"
    
    # Find cut positions
    cuts = find_cut_positions_linear(seq, "EcoRI", db)
    print(f"Cut positions: {cuts}")
    assert cuts == [], f"Expected [], got {cuts}"
    
    # Calculate fragments
    fragment_lengths = fragments_linear(len(seq), cuts)
    print(f"Fragment lengths: {fragment_lengths}")
    assert fragment_lengths == [len(seq)], f"Expected [{len(seq)}], got {fragment_lengths}"
    
    print("✓ Test 4 passed\n")


def test_closest_enzyme_names():
    """Test the enzyme name similarity matching."""
    print("Test 5: Enzyme name similarity matching")
    
    available_names = ["EcoRI", "BamHI", "HindIII", "PstI", "NotI"]
    
    # Test exact match (case-insensitive)
    closest = find_closest_enzyme_names("ecori", available_names)
    print(f"Closest to 'ecori': {closest}")
    assert closest == ["EcoRI"], f"Expected ['EcoRI'], got {closest}"
    
    # Test partial match
    closest = find_closest_enzyme_names("Eco", available_names)
    print(f"Closest to 'Eco': {closest}")
    assert "EcoRI" in closest, "Should find EcoRI for 'Eco'"
    
    # Test no match
    closest = find_closest_enzyme_names("XYZ", available_names)
    print(f"Closest to 'XYZ': {closest}")
    assert len(closest) == 0, "Should find no matches for 'XYZ'"
    
    print("✓ Test 5 passed\n")


def run_all_tests():
    """Run all test cases."""
    print("Running multi-enzyme restriction enzyme simulator tests...")
    print("=" * 60)
    
    try:
        test_single_enzyme_one_site()
        test_two_enzymes_distinct_sites()
        test_two_enzymes_same_coordinate()
        test_no_cuts()
        test_closest_enzyme_names()
        
        print("=" * 60)
        print("🎉 All tests passed!")
        
    except AssertionError as e:
        print(f"❌ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
