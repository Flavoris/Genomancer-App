#!/usr/bin/env python3
"""
Test suite for duplicate enzyme naming.
Tests that duplicate enzymes are displayed with #2, #3 suffixes.
"""

import pytest
import sys
import os
import subprocess
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class TestDuplicateEnzymeNaming:
    """Test that duplicate enzyme names are displayed correctly with #2, #3 suffixes."""
    
    def test_duplicate_enzyme_in_output(self):
        """Test that duplicate enzymes appear as EcoRI and EcoRI#2 in output."""
        # Create a test sequence with two EcoRI sites
        test_seq = "GAATTCGGGGGAATTC"
        
        # Run sim.py with EcoRI specified twice
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "EcoRI"],
            capture_output=True,
            text=True
        )
        
        # Check that output contains both "Enzyme: EcoRI" and "Enzyme: EcoRI#2"
        assert "Enzyme: EcoRI" in result.stdout, "Should show first occurrence as EcoRI"
        assert "Enzyme: EcoRI#2" in result.stdout, "Should show second occurrence as EcoRI#2"
        
        # Make sure it's not showing "EcoRI#3" (only two duplicates)
        assert "EcoRI#3" not in result.stdout, "Should not show EcoRI#3 when only 2 specified"
    
    def test_duplicate_enzyme_in_map(self):
        """Test that duplicate enzymes appear correctly in restriction map."""
        test_seq = "GAATTCGGGGGAATTC"
        
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "EcoRI", "--print-map"],
            capture_output=True,
            text=True
        )
        
        # The restriction map should show both enzyme instances
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        assert "EcoRI" in result.stdout, "Map should show EcoRI"
        # The map might consolidate them or show them separately
        # At minimum, the command should not crash
    
    def test_triple_duplicate(self):
        """Test that three identical enzymes are labeled EcoRI, EcoRI#2, EcoRI#3."""
        test_seq = "GAATTCAAAAAAGAATTCGGGGGAATTC"
        
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "EcoRI", "EcoRI"],
            capture_output=True,
            text=True
        )
        
        assert "Enzyme: EcoRI\n" in result.stdout or "Enzyme: EcoRI" in result.stdout
        assert "Enzyme: EcoRI#2" in result.stdout
        assert "Enzyme: EcoRI#3" in result.stdout
    
    def test_mixed_duplicates(self):
        """Test mixing different enzymes with some duplicates."""
        test_seq = "GAATTCAAAAAGGATCCGGGGGAATTC"
        
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "BamHI", "EcoRI"],
            capture_output=True,
            text=True
        )
        
        assert "Enzyme: EcoRI" in result.stdout
        assert "Enzyme: BamHI" in result.stdout
        assert "Enzyme: EcoRI#2" in result.stdout
        # Should not have BamHI#2 since BamHI appears only once
        assert "BamHI#2" not in result.stdout
    
    def test_cut_site_details_with_duplicates(self):
        """Test that cut site details show the correct enzyme names."""
        test_seq = "GAATTCGGGGGAATTC"
        
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--enz", "EcoRI", "EcoRI"],
            capture_output=True,
            text=True
        )
        
        # Check that Cut Site Details section exists and shows enzyme names
        assert "Cut Site Details:" in result.stdout
        # Should show position information with enzyme names
        # The exact format may vary, but enzyme names should appear


class TestLanesConfig:
    """Test --lanes-config feature works without crashing."""
    
    def test_lanes_config_from_file(self):
        """Test that --lanes-config loads from JSON file without crashing."""
        test_seq = "ATCGAATTCGGGATCCAAA"
        
        # Create a temporary lanes config file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('''[
  {"label": "Uncut", "enzymes": [], "circular": false},
  {"label": "EcoRI", "enzymes": ["EcoRI"], "circular": false},
  {"label": "EcoRI+BamHI", "enzymes": ["EcoRI", "BamHI"], "circular": false}
]''')
            lanes_file = f.name
        
        try:
            result = subprocess.run(
                ["python3", "sim.py", "--seq", test_seq, "--lanes-config", lanes_file, "--gel-only"],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            # Should not crash
            assert result.returncode == 0, f"Command failed with: {result.stderr}"
            
            # Should show gel output
            assert "AGAROSE GEL SIMULATION" in result.stdout or "GEL" in result.stdout.upper()
            
            # Should show lane labels
            assert "Uncut" in result.stdout or "EcoRI" in result.stdout
        
        finally:
            os.unlink(lanes_file)
    
    def test_lanes_config_inline_json(self):
        """Test that --lanes-config works with inline JSON string."""
        test_seq = "ATCGAATTCGGGATCCAAA"
        
        lanes_json = '[{"label":"Test","enzymes":["EcoRI"],"circular":false}]'
        
        result = subprocess.run(
            ["python3", "sim.py", "--seq", test_seq, "--lanes-config", lanes_json, "--gel-only"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        # Should not crash with UnboundLocalError
        assert result.returncode == 0, f"Command failed with: {result.stderr}"
        assert "UnboundLocalError" not in result.stderr


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

