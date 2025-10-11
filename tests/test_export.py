#!/usr/bin/env python3
"""
Test suite for GenBank and CSV export functionality.
"""

import pytest
import os
import sys
import tempfile
import csv

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from exporters import (
    export_genbank, export_csv, gb_date_today, sanitize_locus_name,
    gb_loc_linear, gb_loc_wrap, wrap_origin
)


# ============================================================================
# HELPER FUNCTION TESTS
# ============================================================================

def test_sanitize_locus_name():
    """Test LOCUS name sanitization."""
    assert sanitize_locus_name("test.gbk") == "test"
    assert sanitize_locus_name("/path/to/my-file.gb") == "my_file"
    assert sanitize_locus_name("special!@#$chars.gbk") == "special____chars"
    
    # Test truncation to 16 chars
    long_name = "a" * 30 + ".gbk"
    assert len(sanitize_locus_name(long_name)) == 16


def test_gb_date_today():
    """Test GenBank date format."""
    date_str = gb_date_today()
    # Should be DD-MMM-YYYY format
    assert len(date_str) == 11  # e.g., "11-OCT-2025"
    assert date_str[2] == "-"
    assert date_str[6] == "-"
    # Month should be uppercase
    assert date_str[3:6].isupper()


def test_gb_loc_linear():
    """Test linear GenBank location strings."""
    # Simple linear location
    assert gb_loc_linear(0, 100, 1000) == "1..100"
    assert gb_loc_linear(50, 150, 1000) == "51..150"
    assert gb_loc_linear(999, 1000, 1000) == "1000..1000"


def test_gb_loc_wrap():
    """Test circular wrap GenBank location strings."""
    # Wrap around origin
    assert gb_loc_wrap(900, 100, 1000) == "join(901..1000,1..100)"
    assert gb_loc_wrap(50, 25, 100) == "join(51..100,1..25)"


def test_wrap_origin():
    """Test ORIGIN sequence formatting."""
    seq = "ATCG" * 20  # 80 bases
    formatted = wrap_origin(seq)
    
    lines = formatted.split('\n')
    # Should have 2 lines (60 nt/line)
    assert len(lines) == 2
    
    # Check first line format
    first_line = lines[0]
    assert first_line.startswith("        1")  # 9-char index field
    assert "atcgatcgat" in first_line  # Lowercase, 10-nt blocks


# ============================================================================
# GENBANK EXPORT TESTS
# ============================================================================

def test_genbank_export_linear_two_cuts():
    """Test GenBank export for linear DNA with two cuts."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gbk', delete=False) as f:
        temp_path = f.name
    
    try:
        # Linear sequence with two EcoRI sites
        sequence = "ATCGGAATTCGGGAATTCTA"  # 20 bp, cuts at pos 5 and 14
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 4
            },
            {
                'pos': 14,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 4
            }
        ]
        
        fragments = [
            {
                'index': 1,
                'start': 0,
                'end': 5,
                'length': 5,
                'wraps': False,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    }
                }
            },
            {
                'index': 2,
                'start': 5,
                'end': 14,
                'length': 9,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    },
                    'right_cut': {
                        'pos': 14,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    }
                }
            },
            {
                'index': 3,
                'start': 14,
                'end': 20,
                'length': 6,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 14,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    },
                    'right_cut': None
                }
            }
        ]
        
        export_genbank(
            sequence=sequence,
            cuts=cuts,
            fragments=fragments,
            path=temp_path,
            topology='linear',
            definition='Test linear digest',
            organism='Test organism'
        )
        
        # Read and validate the file
        with open(temp_path, 'r') as f:
            content = f.read()
        
        # Check key sections
        assert "LOCUS" in content
        assert "linear" in content
        assert str(len(sequence)) in content
        assert "DEFINITION  Test linear digest" in content
        assert "SOURCE      Test organism" in content
        assert "FEATURES" in content
        assert "source" in content
        assert "misc_feature" in content
        assert "EcoRI" in content
        assert "ORIGIN" in content
        assert "//" in content
        
        # Validate sequence length in LOCUS
        lines = content.split('\n')
        locus_line = [l for l in lines if l.startswith("LOCUS")][0]
        assert "20 bp" in locus_line
        
        # Count features (should have 1 source + 2 sites + 3 fragments = 6)
        feature_count = content.count("misc_feature")
        assert feature_count == 5  # 2 cuts + 3 fragments
        
        # Validate ORIGIN section has correct sequence (remove spaces for comparison)
        content_no_spaces = content.replace(' ', '').replace('\n', '')
        seq_no_spaces = sequence.replace(' ', '').replace('\n', '')
        assert seq_no_spaces.lower() in content_no_spaces.lower()
    
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def test_genbank_export_circular_wrap():
    """Test GenBank export for circular DNA with wrap fragment."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gbk', delete=False) as f:
        temp_path = f.name
    
    try:
        # Circular sequence with one cut near end
        sequence = "ATCGGAATTCTA"  # 12 bp, cut at pos 5
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 4
            }
        ]
        
        # In circular mode with one cut, we get one fragment that wraps
        fragments = [
            {
                'index': 1,
                'start': 5,
                'end': 5,  # Wraps around
                'length': 12,
                'wraps': True,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    },
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    }
                }
            }
        ]
        
        export_genbank(
            sequence=sequence,
            cuts=cuts,
            fragments=fragments,
            path=temp_path,
            topology='circular',
            definition='Test circular digest',
            organism='synthetic DNA'
        )
        
        # Read and validate the file
        with open(temp_path, 'r') as f:
            content = f.read()
        
        # Check circular-specific features
        assert "circular" in content
        assert "LOCUS" in content
        assert "FEATURES" in content
        
        # Check for source feature with circular note
        assert "source" in content
        assert '/note="circular"' in content
        
        # Check for join() in wrapped fragment
        assert "join(" in content
    
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def test_genbank_overhang_in_notes():
    """Test that overhang information appears in feature notes."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gbk', delete=False) as f:
        temp_path = f.name
    
    try:
        sequence = "ATCGGAATTCTA"
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 4
            }
        ]
        
        fragments = [
            {
                'index': 1,
                'start': 0,
                'end': 5,
                'length': 5,
                'wraps': False,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    }
                }
            },
            {
                'index': 2,
                'start': 5,
                'end': 12,
                'length': 7,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang"}]
                    },
                    'right_cut': None
                }
            }
        ]
        
        export_genbank(
            sequence=sequence,
            cuts=cuts,
            fragments=fragments,
            path=temp_path,
            topology='linear',
            definition='Test overhang notes',
            organism='synthetic DNA'
        )
        
        # Read and validate the file
        with open(temp_path, 'r') as f:
            content = f.read()
        
        # Check for overhang info in notes
        assert "overhang=5' overhang" in content
        assert "k=4" in content
        assert "site=GAATTC" in content
        assert "cut_index=1" in content
    
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


# ============================================================================
# CSV EXPORT TESTS
# ============================================================================

def test_csv_export_schema():
    """Test CSV export creates correct schema with all required columns."""
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test")
        
        sequence = "ATCGGAATTCTA"  # 12 bp
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 4
            }
        ]
        
        fragments = [
            {
                'index': 1,
                'start': 0,
                'end': 5,
                'length': 5,
                'wraps': False,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 4}]
                    }
                }
            },
            {
                'index': 2,
                'start': 5,
                'end': 12,
                'length': 7,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 4}]
                    },
                    'right_cut': None
                }
            }
        ]
        
        export_csv(
            prefix=prefix,
            cuts=cuts,
            fragments=fragments,
            topology='linear',
            dna_sequence=sequence
        )
        
        # Check fragments CSV
        frag_path = f"{prefix}_fragments.csv"
        assert os.path.exists(frag_path)
        
        with open(frag_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            # Check columns
            expected_cols = [
                'fragment_id', 'start_idx', 'end_idx', 'mode', 'length',
                'left_enzyme', 'left_overhang_type', 'left_overhang_len', 'left_end_bases',
                'right_enzyme', 'right_overhang_type', 'right_overhang_len', 'right_end_bases',
                'sequence'
            ]
            assert all(col in reader.fieldnames for col in expected_cols)
            
            # Check number of rows
            assert len(rows) == 2
            
            # Validate data
            assert rows[0]['fragment_id'] == '1'
            assert rows[0]['length'] == '5'
            assert rows[0]['sequence'] == 'ATCGG'
            assert rows[1]['fragment_id'] == '2'
            assert rows[1]['length'] == '7'
            assert rows[1]['sequence'] == 'AATTCTA'
        
        # Check cuts CSV
        cuts_path = f"{prefix}_cuts.csv"
        assert os.path.exists(cuts_path)
        
        with open(cuts_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            # Check columns
            expected_cols = [
                'cut_id', 'pos', 'enzyme', 'recognition_site',
                'cut_index', 'overhang_type', 'overhang_len'
            ]
            assert all(col in reader.fieldnames for col in expected_cols)
            
            # Check data
            assert len(rows) == 1
            assert rows[0]['enzyme'] == 'EcoRI'
            assert rows[0]['pos'] == '5'
            assert rows[0]['overhang_type'] == "5' overhang"


def test_csv_fragment_lengths_match():
    """Test that CSV fragment lengths match actual sequence lengths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test")
        
        # Use a longer sequence with multiple cuts
        sequence = "ATCGGAATTCGGGCTGCAGTTTAAGCTTAA"  # 30 bp
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 5
            },
            {
                'pos': 17,
                'enzyme': 'PstI',
                'recognition_site': 'CTGCAG',
                'cut_index': 5,
                'overhang_type': "3' overhang",
                'overhang_len': 5
            }
        ]
        
        fragments = [
            {
                'index': 1,
                'start': 0,
                'end': 5,
                'length': 5,
                'wraps': False,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 5}]
                    }
                }
            },
            {
                'index': 2,
                'start': 5,
                'end': 17,
                'length': 12,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 5}]
                    },
                    'right_cut': {
                        'pos': 17,
                        'enzymes': [{'enzyme': 'PstI', 'overhang_type': "3' overhang", 'overhang_len': 5}]
                    }
                }
            },
            {
                'index': 3,
                'start': 17,
                'end': 30,
                'length': 13,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': 17,
                        'enzymes': [{'enzyme': 'PstI', 'overhang_type': "3' overhang", 'overhang_len': 5}]
                    },
                    'right_cut': None
                }
            }
        ]
        
        export_csv(
            prefix=prefix,
            cuts=cuts,
            fragments=fragments,
            topology='linear',
            dna_sequence=sequence
        )
        
        # Read fragments CSV and verify lengths
        frag_path = f"{prefix}_fragments.csv"
        with open(frag_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Length should equal len(sequence)
                assert int(row['length']) == len(row['sequence'])


def test_csv_circular_wrap_fragment():
    """Test CSV export correctly handles circular wrap fragments."""
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test")
        
        sequence = "ATCGGAATTCTA"  # 12 bp, cut at 5
        
        cuts = [
            {
                'pos': 5,
                'enzyme': 'EcoRI',
                'recognition_site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang",
                'overhang_len': 5
            }
        ]
        
        # Single wrap fragment
        fragments = [
            {
                'index': 1,
                'start': 5,
                'end': 5,
                'length': 12,
                'wraps': True,
                'boundaries': {
                    'left_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 5}]
                    },
                    'right_cut': {
                        'pos': 5,
                        'enzymes': [{'enzyme': 'EcoRI', 'overhang_type': "5' overhang", 'overhang_len': 5}]
                    }
                }
            }
        ]
        
        export_csv(
            prefix=prefix,
            cuts=cuts,
            fragments=fragments,
            topology='circular',
            dna_sequence=sequence
        )
        
        # Read and validate
        frag_path = f"{prefix}_fragments.csv"
        with open(frag_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            assert len(rows) == 1
            row = rows[0]
            
            # Check mode is circular
            assert row['mode'] == 'circular'
            
            # Sequence should wrap: from pos 5 to end + start to pos 5
            expected_seq = sequence[5:] + sequence[:5]
            assert row['sequence'] == expected_seq
            assert int(row['length']) == len(expected_seq)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

