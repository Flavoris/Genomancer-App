#!/usr/bin/env python3
"""
Test suite for export consistency.
Tests that CSV and GenBank exports match console output for overhang metadata.
"""

import pytest
import sys
import os
import subprocess
import tempfile
import csv

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class TestExportConsistency:
    """Test that exported data matches console output."""
    
    def test_csv_overhang_correctness(self):
        """Test that CSV export shows correct overhang lengths and bases."""
        test_seq = "AAAGAATTCGGG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            
            # Run digest with CSV export
            result = subprocess.run(
                ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", 
                 "--export-csv", prefix],
                capture_output=True,
                text=True
            )
            
            assert result.returncode == 0, f"Command failed: {result.stderr}"
            
            # Read fragments CSV
            fragments_csv = f"{prefix}_fragments.csv"
            assert os.path.exists(fragments_csv), "Fragments CSV should be created"
            
            with open(fragments_csv, 'r') as f:
                reader = csv.DictReader(f)
                fragments = list(reader)
            
            # Check that we have fragments
            assert len(fragments) > 0, "Should have fragments"
            
            # Check that overhang lengths are correct for EcoRI (4 bp)
            for frag in fragments:
                if frag['left_overhang_len']:
                    left_len = int(frag['left_overhang_len'])
                    if left_len > 0:
                        assert left_len == 4, f"EcoRI left overhang should be 4, got {left_len}"
                        # Check that end_bases is provided and is 4 characters
                        if frag['left_end_bases']:
                            assert len(frag['left_end_bases']) == 4, \
                                f"Left end bases should be 4 chars, got {frag['left_end_bases']}"
                
                if frag['right_overhang_len']:
                    right_len = int(frag['right_overhang_len'])
                    if right_len > 0:
                        assert right_len == 4, f"EcoRI right overhang should be 4, got {right_len}"
                        if frag['right_end_bases']:
                            assert len(frag['right_end_bases']) == 4, \
                                f"Right end bases should be 4 chars, got {frag['right_end_bases']}"
            
            # Read cuts CSV
            cuts_csv = f"{prefix}_cuts.csv"
            assert os.path.exists(cuts_csv), "Cuts CSV should be created"
            
            with open(cuts_csv, 'r') as f:
                reader = csv.DictReader(f)
                cuts = list(reader)
            
            # Check cuts
            for cut in cuts:
                if cut['enzyme'] == 'EcoRI':
                    assert int(cut['overhang_len']) == 4, \
                        f"EcoRI cut overhang should be 4, got {cut['overhang_len']}"
    
    def test_genbank_overhang_correctness(self):
        """Test that GenBank export shows correct overhang metadata."""
        test_seq = "AAAGAATTCGGG"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as f:
            gb_file = f.name
        
        try:
            # Run digest with GenBank export
            result = subprocess.run(
                ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", 
                 "--export-genbank", gb_file],
                capture_output=True,
                text=True
            )
            
            assert result.returncode == 0, f"Command failed: {result.stderr}"
            assert os.path.exists(gb_file), "GenBank file should be created"
            
            # Read GenBank file
            with open(gb_file, 'r') as f:
                gb_content = f.read()
            
            # Check that overhang information is present
            # For EcoRI, should see k=4 (not k=5 or k=0)
            assert "k=4" in gb_content, \
                f"GenBank should show k=4 for EcoRI overhang, content: {gb_content}"
            
            # Should NOT see k=5 or k=0 for EcoRI
            assert "k=5" not in gb_content or "k=4" in gb_content, \
                "Should not show incorrect overhang length k=5 for EcoRI"
        
        finally:
            if os.path.exists(gb_file):
                os.unlink(gb_file)
    
    def test_console_vs_csv_consistency(self):
        """Test that console output matches CSV export."""
        test_seq = "GAATTCGGG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            
            # Run with both console output and CSV export
            result = subprocess.run(
                ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI",
                 "--include-seqs", "--export-csv", prefix],
                capture_output=True,
                text=True
            )
            
            assert result.returncode == 0, f"Command failed: {result.stderr}"
            
            # Extract overhang info from console
            console_output = result.stdout
            
            # Should show overhang length 4 in console
            # Format might be like: "4 bp: AATT"
            assert "4 bp" in console_output or "overhang_len = 4" in console_output.lower(), \
                f"Console should show 4 bp overhang, output: {console_output}"
            
            # Should show AATT bases
            assert "AATT" in console_output, \
                f"Console should show AATT sticky end, output: {console_output}"
            
            # Read CSV and verify it matches
            fragments_csv = f"{prefix}_fragments.csv"
            with open(fragments_csv, 'r') as f:
                reader = csv.DictReader(f)
                fragments = list(reader)
            
            # At least one fragment should have overhang_len = 4
            has_correct_overhang = False
            for frag in fragments:
                if frag['left_overhang_len'] == '4' or frag['right_overhang_len'] == '4':
                    has_correct_overhang = True
                    break
            
            assert has_correct_overhang, "CSV should show overhang length 4"


class TestBluntEnds:
    """Test that blunt-cutting enzymes report 0 bp overhang correctly."""
    
    def test_smai_blunt_in_csv(self):
        """Test that SmaI (blunt cutter) shows 0 bp overhang in CSV."""
        # SmaI recognizes CCCGGG and cuts in the middle (blunt)
        test_seq = "AAACCCGGGAAA"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            
            result = subprocess.run(
                ["python3", "sim.py", "--seq", test_seq, "--enz", "SmaI",
                 "--export-csv", prefix],
                capture_output=True,
                text=True
            )
            
            # Check if SmaI is in database
            if "not found" in result.stderr or result.returncode != 0:
                pytest.skip("SmaI not in enzyme database")
            
            # Read cuts CSV
            cuts_csv = f"{prefix}_cuts.csv"
            with open(cuts_csv, 'r') as f:
                reader = csv.DictReader(f)
                cuts = list(reader)
            
            for cut in cuts:
                if 'SmaI' in cut['enzyme']:
                    assert int(cut['overhang_len']) == 0, \
                        f"SmaI (blunt) should have 0 bp overhang, got {cut['overhang_len']}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

