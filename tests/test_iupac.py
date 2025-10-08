#!/usr/bin/env python3
"""
Test suite for IUPAC degenerate base support in restriction enzyme simulator.
"""

import sys
import os
import unittest
import tempfile
import json

# Add parent directory to path to import sim module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim import iupac_to_regex, find_cut_sites, load_enzyme_database


class TestIUPACFunctionality(unittest.TestCase):
    """Test IUPAC degenerate base functionality."""
    
    def test_iupac_to_regex_basic(self):
        """Test basic IUPAC to regex conversion."""
        self.assertEqual(iupac_to_regex("A"), "A")
        self.assertEqual(iupac_to_regex("C"), "C")
        self.assertEqual(iupac_to_regex("G"), "G")
        self.assertEqual(iupac_to_regex("T"), "T")
    
    def test_iupac_to_regex_degenerate(self):
        """Test degenerate base conversion."""
        self.assertEqual(iupac_to_regex("R"), "[AG]")
        self.assertEqual(iupac_to_regex("Y"), "[CT]")
        self.assertEqual(iupac_to_regex("W"), "[AT]")
        self.assertEqual(iupac_to_regex("S"), "[GC]")
        self.assertEqual(iupac_to_regex("M"), "[AC]")
        self.assertEqual(iupac_to_regex("K"), "[GT]")
        self.assertEqual(iupac_to_regex("N"), "[ACGT]")
    
    def test_iupac_to_regex_mixed(self):
        """Test mixed standard and degenerate bases."""
        self.assertEqual(iupac_to_regex("AGR"), "AG[AG]")
        self.assertEqual(iupac_to_regex("ATW"), "AT[AT]")
        self.assertEqual(iupac_to_regex("AN"), "A[ACGT]")
        self.assertEqual(iupac_to_regex("GCNS"), "GC[ACGT][GC]")
    
    def test_iupac_to_regex_case_insensitive(self):
        """Test case insensitivity."""
        self.assertEqual(iupac_to_regex("r"), "[AG]")
        self.assertEqual(iupac_to_regex("AgR"), "AG[AG]")
        self.assertEqual(iupac_to_regex("atw"), "AT[AT]")
    
    def test_iupac_to_regex_invalid_chars(self):
        """Test validation of invalid characters."""
        with self.assertRaises(ValueError) as context:
            iupac_to_regex("AXG")
        self.assertIn("Invalid character 'X'", str(context.exception))
        
        with self.assertRaises(ValueError) as context:
            iupac_to_regex("ABG")
        self.assertIn("Invalid character 'B'", str(context.exception))


class TestIUPACCutSiteFinding(unittest.TestCase):
    """Test cut site finding with IUPAC patterns."""
    
    def test_iupac_simple_R(self):
        """Test enzyme site 'AGR' with cut_index=1 should match 'AGA' and 'AGG'."""
        # Test sequence with both AGA and AGG
        sequence = "ATAGATAGGCC"
        cut_positions = find_cut_sites(sequence, "AGR", 1)
        # Should find cuts at positions: 3 (AGA) and 7 (AGG)
        self.assertEqual(cut_positions, [3, 7])
    
    def test_iupac_W(self):
        """Test site 'ATW' matches 'ATA' and 'ATT'."""
        sequence = "ATATATTCC"
        cut_positions = find_cut_sites(sequence, "ATW", 2)
        # Should find cuts at positions: 2, 4, 6 (overlapping ATW matches)
        self.assertEqual(cut_positions, [2, 4, 6])
    
    def test_iupac_N(self):
        """Test site 'AN' matches 'AA', 'AC', 'AG', 'AT'."""
        sequence = "AATACGAT"
        cut_positions = find_cut_sites(sequence, "AN", 1)
        # Should find cuts at positions: 1 (AA), 2 (AC), 4 (AG), 7 (AT)
        self.assertEqual(cut_positions, [1, 2, 4, 7])
    
    def test_iupac_overlap(self):
        """Test overlapping matches."""
        sequence = "AAAAA"
        cut_positions = find_cut_sites(sequence, "AA", 1)
        # Should find overlapping matches at positions: 1,2,3,4
        self.assertEqual(cut_positions, [1, 2, 3, 4])
    
    def test_iupac_no_matches(self):
        """Test when no matches are found."""
        sequence = "CCCCCCC"
        cut_positions = find_cut_sites(sequence, "AGR", 1)
        self.assertEqual(cut_positions, [])
    
    def test_iupac_case_insensitive_sequence(self):
        """Test case insensitive matching on DNA sequence."""
        sequence = "atagataggcc"
        cut_positions = find_cut_sites(sequence, "AGR", 1)
        self.assertEqual(cut_positions, [3, 7])


class TestEnzymeDatabaseWithIUPAC(unittest.TestCase):
    """Test enzyme database loading with IUPAC sites."""
    
    def test_load_valid_iupac_enzyme(self):
        """Test loading enzyme database with valid IUPAC sites."""
        test_enzymes = [
            {"name": "TestR", "site": "AGR", "cut_index": 1},
            {"name": "TestW", "site": "ATW", "cut_index": 2},
            {"name": "TestN", "site": "AN", "cut_index": 1},
            {"name": "TestStandard", "site": "GATC", "cut_index": 2}
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_enzymes, f)
            temp_file = f.name
        
        try:
            # Mock the file loading by temporarily replacing the file path
            import sim
            original_open = open
            def mock_open(filename, mode='r'):
                if filename == "enzymes.json":
                    return original_open(temp_file, mode)
                return original_open(filename, mode)
            
            sim.open = mock_open
            
            enzymes = load_enzyme_database()
            
            self.assertIn("TestR", enzymes)
            self.assertIn("TestW", enzymes)
            self.assertIn("TestN", enzymes)
            self.assertIn("TestStandard", enzymes)
            
            self.assertEqual(enzymes["TestR"]["sequence"], "AGR")
            self.assertEqual(enzymes["TestW"]["sequence"], "ATW")
            
        finally:
            os.unlink(temp_file)
            sim.open = original_open
    
    def test_load_invalid_iupac_enzyme(self):
        """Test loading enzyme database with invalid characters."""
        test_enzymes = [
            {"name": "ValidEnzyme", "site": "AGR", "cut_index": 1},
            {"name": "InvalidEnzyme", "site": "AXG", "cut_index": 1},
            {"name": "AnotherValid", "site": "GATC", "cut_index": 2}
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_enzymes, f)
            temp_file = f.name
        
        try:
            # Mock the file loading
            import sim
            original_open = open
            def mock_open(filename, mode='r'):
                if filename == "enzymes.json":
                    return original_open(temp_file, mode)
                return original_open(filename, mode)
            
            sim.open = mock_open
            
            # Capture stdout to check warning messages
            import io
            import contextlib
            
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                enzymes = load_enzyme_database()
            
            output = f.getvalue()
            
            # Should load valid enzymes but skip invalid one
            self.assertIn("ValidEnzyme", enzymes)
            self.assertIn("AnotherValid", enzymes)
            self.assertNotIn("InvalidEnzyme", enzymes)
            
            # Should have warning message
            self.assertIn("Invalid character 'X'", output)
            
        finally:
            os.unlink(temp_file)
            sim.open = original_open
    
    def test_load_invalid_cut_index(self):
        """Test loading enzyme database with invalid cut_index."""
        test_enzymes = [
            {"name": "ValidEnzyme", "site": "AGR", "cut_index": 1},
            {"name": "InvalidCutIndex", "site": "AGR", "cut_index": 5},  # cut_index > len(site)
            {"name": "NegativeCutIndex", "site": "AGR", "cut_index": -1}  # negative cut_index
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_enzymes, f)
            temp_file = f.name
        
        try:
            # Mock the file loading
            import sim
            original_open = open
            def mock_open(filename, mode='r'):
                if filename == "enzymes.json":
                    return original_open(temp_file, mode)
                return original_open(filename, mode)
            
            sim.open = mock_open
            
            # Capture stdout to check warning messages
            import io
            import contextlib
            
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                enzymes = load_enzyme_database()
            
            output = f.getvalue()
            
            # Should only load valid enzyme
            self.assertIn("ValidEnzyme", enzymes)
            self.assertNotIn("InvalidCutIndex", enzymes)
            self.assertNotIn("NegativeCutIndex", enzymes)
            
            # Should have warning messages
            self.assertIn("cut_index 5 must be between 0 and 3", output)
            self.assertIn("cut_index -1 must be between 0 and 3", output)
            
        finally:
            os.unlink(temp_file)
            sim.open = original_open


class TestMultiEnzymeWithIUPAC(unittest.TestCase):
    """Test multi-enzyme behavior with IUPAC sites."""
    
    def test_multi_enzyme_merge_with_iupac(self):
        """Test that two enzymes producing shared cut coordinates dedupe correctly."""
        sequence = "ATAGATAGGCC"
        
        # Enzyme 1: AGR with cut_index=1 (matches AGA at pos 2, AGG at pos 6)
        cuts1 = find_cut_sites(sequence, "AGR", 1)
        self.assertEqual(cuts1, [3, 7])
        
        # Enzyme 2: ATA with cut_index=2 (matches ATA at pos 0 and 4)
        cuts2 = find_cut_sites(sequence, "ATA", 2)
        self.assertEqual(cuts2, [2, 6])
        
        # Combined cuts should dedupe position 6
        all_cuts = sorted(set(cuts1 + cuts2))
        self.assertEqual(all_cuts, [2, 3, 6, 7])


if __name__ == "__main__":
    unittest.main()
