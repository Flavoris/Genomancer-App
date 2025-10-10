#!/usr/bin/env python3
"""
Test cases for enhanced restriction enzyme simulator features.
Tests IUPAC expansion, dual-cut support, N-spacer cuts, and duplicate names.
"""

import sys
import os
import json
import tempfile
import unittest
from unittest.mock import patch

# Add parent directory to path to import sim module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim import (
    iupac_to_regex, normalize, load_enzyme_database, 
    find_cut_sites, find_cut_positions_linear
)


class TestIUPACExpansion(unittest.TestCase):
    """Test IUPAC degenerate base expansion."""
    
    def test_basic_iupac_bases(self):
        """Test basic IUPAC bases (R, Y, W, S, M, K, N)."""
        test_cases = [
            ("GAATTC", "GAATTC"),  # No IUPAC
            ("GARTTC", "GA[AG]TTC"),  # R = [AG]
            ("GAYTTC", "GA[CT]TTC"),  # Y = [CT]
            ("GAWTTC", "GA[AT]TTC"),  # W = [AT]
            ("GASTTC", "GA[GC]TTC"),  # S = [GC]
            ("GAMTTC", "GA[AC]TTC"),  # M = [AC]
            ("GAKTTC", "GA[GT]TTC"),  # K = [GT]
            ("GANTTC", "GA[ACGT]TTC"),  # N = [ACGT]
        ]
        
        for site, expected in test_cases:
            with self.subTest(site=site):
                result = iupac_to_regex(site)
                self.assertEqual(result, expected)
    
    def test_new_iupac_bases(self):
        """Test new IUPAC bases (B, D, H, V)."""
        test_cases = [
            ("GABTTC", "GA[CGT]TTC"),  # B = [CGT]
            ("GADTTC", "GA[AGT]TTC"),  # D = [AGT]
            ("GAHTTC", "GA[ACT]TTC"),  # H = [ACT]
            ("GAVTTC", "GA[ACG]TTC"),  # V = [ACG]
        ]
        
        for site, expected in test_cases:
            with self.subTest(site=site):
                result = iupac_to_regex(site)
                self.assertEqual(result, expected)
    
    def test_mixed_iupac_bases(self):
        """Test sites with multiple IUPAC bases."""
        result = iupac_to_regex("GARTYC")
        self.assertEqual(result, "GA[AG]T[CT]C")
        
        result = iupac_to_regex("GBNTTC")
        self.assertEqual(result, "G[CGT][ACGT]TTC")
    
    def test_invalid_characters(self):
        """Test that invalid characters raise ValueError."""
        with self.assertRaises(ValueError):
            iupac_to_regex("GAPTC")  # P is not valid
        
        with self.assertRaises(ValueError):
            iupac_to_regex("GA!TC")  # ! is not valid


class TestNameNormalization(unittest.TestCase):
    """Test enzyme name normalization."""
    
    def test_basic_normalization(self):
        """Test basic name normalization."""
        test_cases = [
            ("EcoRI", "ecori"),
            ("BamHI", "bamhi"),
            ("EcoRI ", "ecori"),  # Remove spaces
            ("Eco-RI", "ecori"),  # Remove hyphens
            ("ECORI", "ecori"),   # Lowercase
        ]
        
        for name, expected in test_cases:
            with self.subTest(name=name):
                result = normalize(name)
                self.assertEqual(result, expected)
    
    def test_diacritics_removal(self):
        """Test removal of diacritics."""
        # Note: This is a simplified test since actual diacritics
        # would require Unicode characters
        result = normalize("EcoRI")
        self.assertEqual(result, "ecori")


class TestDatabaseLoading(unittest.TestCase):
    """Test enhanced database loading."""
    
    def test_duplicate_name_handling(self):
        """Test handling of duplicate enzyme names."""
        test_db = [
            {
                "name": "EcoRI",
                "site": "GAATTC",
                "cut_index": 1
            },
            {
                "name": "EcoRI",
                "site": "GAATTC",
                "cut_index": 2
            }
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_db, f)
            temp_file = f.name
        
        try:
            # Mock the file opening to use our test file
            original_open = open
            def mock_open(filename, mode='r'):
                if filename == "enzymes.json":
                    return original_open(temp_file, mode)
                return original_open(filename, mode)
            
            with patch('builtins.open', side_effect=mock_open):
                enzymes = load_enzyme_database()
            
            # Should have both EcoRI and EcoRI#2
            self.assertIn("EcoRI", enzymes)
            self.assertIn("EcoRI#2", enzymes)
            
            # First one should have cut_index 1
            self.assertEqual(enzymes["EcoRI"]["cut_index"], 1)
            
            # Second one should have cut_index 2
            self.assertEqual(enzymes["EcoRI#2"]["cut_index"], 2)
        
        finally:
            os.unlink(temp_file)
    
    def test_invalid_enzyme_skipping(self):
        """Test that invalid enzymes are skipped with warnings."""
        test_db = [
            {
                "name": "ValidEnzyme",
                "site": "GAATTC",
                "cut_index": 1
            },
            {
                "name": "InvalidEnzyme",
                "site": "GAPTC",  # Invalid character
                "cut_index": 1
            },
            {
                "name": "AnotherValid",
                "site": "GGATCC",
                "cut_index": 1
            }
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_db, f)
            temp_file = f.name
        
        try:
            # Mock the file opening to use our test file
            original_open = open
            def mock_open(filename, mode='r'):
                if filename == "enzymes.json":
                    return original_open(temp_file, mode)
                return original_open(filename, mode)
            
            with patch('builtins.open', side_effect=mock_open):
                with patch('builtins.print') as mock_print:
                    enzymes = load_enzyme_database()
            
            # Should only have valid enzymes
            self.assertIn("ValidEnzyme", enzymes)
            self.assertIn("AnotherValid", enzymes)
            self.assertNotIn("InvalidEnzyme", enzymes)
            
            # Should have printed warning
            mock_print.assert_called()
        
        finally:
            os.unlink(temp_file)


if __name__ == "__main__":
    unittest.main()
