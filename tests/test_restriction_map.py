#!/usr/bin/env python3
"""
Test suite for restriction map functionality.
Tests the build_restriction_map() function with various scenarios.
"""

import sys
import os
import subprocess

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fragment_calculator import build_restriction_map


def test_scale_and_markers_linear_basic():
    """Test basic linear map with scale and markers."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"},
        {"pos": 500, "enzyme": "HindIII", "overhang_type": "5' overhang", "site": "AAGCTT"}
    ]
    s = build_restriction_map(L, cuts, map_width=40, map_ticks=5, circular=False)
    
    assert "Length: 1000 bp" in s
    assert "Mode: linear" in s
    assert "All (" in s
    assert "^" in s or "|" in s
    print("✓ test_scale_and_markers_linear_basic passed")


def test_circular_wrap_note_present():
    """Test that circular mode with 2+ cuts shows wrap note."""
    L = 1200
    cuts = [
        {"pos": 100, "enzyme": "A", "overhang_type": "Blunt", "site": "NNNN"},
        {"pos": 900, "enzyme": "B", "overhang_type": "5' overhang", "site": "NNNN"}
    ]
    s = build_restriction_map(L, cuts, map_width=60, circular=True)
    
    assert "wrap spans" in s
    assert "900→100" in s
    print("✓ test_circular_wrap_note_present passed")


def test_group_by_position_lists_contributors():
    """Test position grouping mode lists all contributing enzymes."""
    L = 1000
    cuts = [
        {"pos": 250, "enzyme": "A", "overhang_type": "Blunt", "site": "X"},
        {"pos": 250, "enzyme": "B", "overhang_type": "5' overhang", "site": "Y"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, group_by="position")
    
    assert "pos 250:" in s
    assert "A" in s and "B" in s
    print("✓ test_group_by_position_lists_contributors passed")


def test_min_hits_filters_sparse_enzymes():
    """Test that min_hits filters out enzymes with too few cuts."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "A", "overhang_type": "Blunt", "site": "X"},
        {"pos": 200, "enzyme": "A", "overhang_type": "Blunt", "site": "X"},
        {"pos": 300, "enzyme": "B", "overhang_type": "Blunt", "site": "Y"}
    ]
    
    # With min_hits=1, both should appear
    s1 = build_restriction_map(L, cuts, map_width=50, map_min_hits=1)
    assert "A (2)" in s1
    assert "B (1)" in s1
    
    # With min_hits=2, only A should appear
    s2 = build_restriction_map(L, cuts, map_width=50, map_min_hits=2)
    assert "A (2)" in s2
    assert "B (1)" not in s2
    
    print("✓ test_min_hits_filters_sparse_enzymes passed")


def test_no_cuts_displays_message():
    """Test that zero cuts shows appropriate message."""
    L = 1000
    cuts = []
    s = build_restriction_map(L, cuts, map_width=50)
    
    assert "No cut sites for selected enzymes" in s
    assert "Length: 1000 bp" in s
    print("✓ test_no_cuts_displays_message passed")


def test_single_cut_linear():
    """Test single cut on linear DNA."""
    L = 1000
    cuts = [{"pos": 500, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"}]
    s = build_restriction_map(L, cuts, map_width=50, circular=False)
    
    assert "All (1)" in s
    assert "EcoRI (1)" in s
    assert "Mode: linear" in s
    print("✓ test_single_cut_linear passed")


def test_single_cut_circular():
    """Test single cut on circular DNA shows appropriate note."""
    L = 1000
    cuts = [{"pos": 500, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"}]
    s = build_restriction_map(L, cuts, map_width=50, circular=True)
    
    assert "All (1)" in s
    assert "Single cut on circular molecule" in s
    assert "Mode: circular" in s
    print("✓ test_single_cut_circular passed")


def test_shared_positions_marked_with_asterisk():
    """Test that shared cut positions are marked with * on All track."""
    L = 1000
    cuts = [
        {"pos": 500, "enzyme": "A", "overhang_type": "Blunt", "site": "X"},
        {"pos": 500, "enzyme": "B", "overhang_type": "5' overhang", "site": "Y"},
        {"pos": 700, "enzyme": "C", "overhang_type": "Blunt", "site": "Z"}
    ]
    s = build_restriction_map(L, cuts, map_width=60, circular=False)
    
    # Should have an asterisk for the shared position
    assert "*" in s
    print("✓ test_shared_positions_marked_with_asterisk passed")


def test_show_overhangs_flag():
    """Test that show_overhangs adds overhang labels."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"},
        {"pos": 500, "enzyme": "PstI", "overhang_type": "3' overhang", "site": "CTGCAG"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, show_overhangs=True)
    
    assert "5' overhang" in s
    assert "3' overhang" in s
    print("✓ test_show_overhangs_flag passed")


def test_show_sites_flag():
    """Test that show_sites adds recognition sequences."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, show_sites=True)
    
    assert "GAATTC" in s
    print("✓ test_show_sites_flag passed")


def test_site_truncation():
    """Test that very long sites are truncated."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "LongSite", "overhang_type": "Blunt", 
         "site": "AAAAAAAAAAAAAAAAAAAAAAAAAA"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, show_sites=True)
    
    # Should truncate with ...
    assert "..." in s
    print("✓ test_site_truncation passed")


def test_circular_origin_parameter():
    """Test that circular_origin parameter is reflected in output."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "A", "overhang_type": "Blunt", "site": "X"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, circular=True, circular_origin=250)
    
    assert "origin=250" in s
    print("✓ test_circular_origin_parameter passed")


def test_enzyme_sorting():
    """Test that enzymes are sorted by decreasing hit count, then alphabetically."""
    L = 2000
    cuts = [
        {"pos": 100, "enzyme": "Alpha", "overhang_type": "Blunt", "site": "X"},
        {"pos": 200, "enzyme": "Beta", "overhang_type": "Blunt", "site": "Y"},
        {"pos": 300, "enzyme": "Beta", "overhang_type": "Blunt", "site": "Y"},
        {"pos": 400, "enzyme": "Beta", "overhang_type": "Blunt", "site": "Y"},
        {"pos": 500, "enzyme": "Gamma", "overhang_type": "Blunt", "site": "Z"},
        {"pos": 600, "enzyme": "Gamma", "overhang_type": "Blunt", "site": "Z"}
    ]
    s = build_restriction_map(L, cuts, map_width=60)
    
    # Beta (3 hits) should come before Gamma (2 hits) and Alpha (1 hit)
    lines = s.split('\n')
    enzyme_lines = [l for l in lines if '(' in l and ')' in l and ':' in l]
    
    # Find positions of enzyme names in output
    beta_idx = None
    gamma_idx = None
    alpha_idx = None
    
    for i, line in enumerate(enzyme_lines):
        if "Beta" in line:
            beta_idx = i
        elif "Gamma" in line:
            gamma_idx = i
        elif "Alpha" in line:
            alpha_idx = i
    
    # Beta should come before Gamma, and Gamma before Alpha
    if beta_idx is not None and gamma_idx is not None:
        assert beta_idx < gamma_idx
    if gamma_idx is not None and alpha_idx is not None:
        assert gamma_idx < alpha_idx
    
    print("✓ test_enzyme_sorting passed")


def test_map_width_parameter():
    """Test that map_width controls the width of the map."""
    L = 1000
    cuts = [{"pos": 500, "enzyme": "EcoRI", "overhang_type": "5' overhang", "site": "GAATTC"}]
    
    s1 = build_restriction_map(L, cuts, map_width=40)
    s2 = build_restriction_map(L, cuts, map_width=80)
    
    # Extract the ruler line (should be all dashes and pipes)
    lines1 = s1.split('\n')
    lines2 = s2.split('\n')
    
    # Find ruler line (contains mostly dashes)
    ruler1 = [l for l in lines1 if l.count('-') > 10][0]
    ruler2 = [l for l in lines2 if l.count('-') > 10][0]
    
    # Width should be approximately as specified
    assert len(ruler1) <= 45  # Some tolerance
    assert len(ruler2) >= 75  # Some tolerance
    
    print("✓ test_map_width_parameter passed")


def test_degenerate_bases_in_site():
    """Test that degenerate bases in sites are displayed correctly."""
    L = 1000
    cuts = [
        {"pos": 100, "enzyme": "BsrI", "overhang_type": "5' overhang", "site": "ACTGGN"}
    ]
    s = build_restriction_map(L, cuts, map_width=50, show_sites=True)
    
    assert "ACTGGN" in s
    print("✓ test_degenerate_bases_in_site passed")


def test_cli_integration():
    """Test CLI integration with --print-map flag."""
    # Create a simple test sequence
    test_seq = "ATCGAATTCGGGATCCAAA"
    
    # Run with --print-map
    result = subprocess.run(
        ["python", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "BamHI", "--print-map"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.dirname(__file__))
    )
    
    assert result.returncode == 0
    assert "RESTRICTION MAP" in result.stdout
    assert "All (" in result.stdout
    print("✓ test_cli_integration passed")


def test_cli_circular_integration():
    """Test CLI integration with circular mode and --print-map."""
    test_seq = "ATCGAATTCGGGATCCAAA"
    
    result = subprocess.run(
        ["python", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "BamHI", 
         "--circular", "--print-map", "--map-width", "60"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.dirname(__file__))
    )
    
    assert result.returncode == 0
    assert "Mode: circular" in result.stdout
    assert "wrap spans" in result.stdout or "Single cut" in result.stdout
    print("✓ test_cli_circular_integration passed")


def test_print_map_only_flag():
    """Test that --print-map-only skips fragment table."""
    test_seq = "ATCGAATTCGGGATCCAAA"
    
    result = subprocess.run(
        ["python", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "--print-map-only"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.dirname(__file__))
    )
    
    assert result.returncode == 0
    assert "RESTRICTION MAP" in result.stdout
    # Should not have fragment details table
    assert "Fragment Details:" not in result.stdout
    print("✓ test_print_map_only_flag passed")


def run_all_tests():
    """Run all test functions."""
    test_functions = [
        test_scale_and_markers_linear_basic,
        test_circular_wrap_note_present,
        test_group_by_position_lists_contributors,
        test_min_hits_filters_sparse_enzymes,
        test_no_cuts_displays_message,
        test_single_cut_linear,
        test_single_cut_circular,
        test_shared_positions_marked_with_asterisk,
        test_show_overhangs_flag,
        test_show_sites_flag,
        test_site_truncation,
        test_circular_origin_parameter,
        test_enzyme_sorting,
        test_map_width_parameter,
        test_degenerate_bases_in_site,
        test_cli_integration,
        test_cli_circular_integration,
        test_print_map_only_flag
    ]
    
    print("\n" + "="*70)
    print("Running Restriction Map Tests")
    print("="*70 + "\n")
    
    passed = 0
    failed = 0
    
    for test_func in test_functions:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"✗ {test_func.__name__} FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {test_func.__name__} ERROR: {e}")
            failed += 1
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*70 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

