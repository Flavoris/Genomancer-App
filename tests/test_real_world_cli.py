#!/usr/bin/env python3
"""
Command-Line Interface Real-World Testing
Tests the actual CLI commands as documented in README
"""

import sys
import os
import subprocess
import tempfile
import shutil
import json
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestCLI:
    """Test CLI commands with real-world scenarios"""
    
    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()
        self.script_path = Path(__file__).parent.parent / "sim.py"
        
    def teardown_method(self):
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def run_sim(self, args):
        """Helper to run sim.py with arguments"""
        cmd = ["python", str(self.script_path)] + args
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode, result.stdout, result.stderr
    
    def test_basic_linear_digest(self):
        """Test basic linear DNA digestion"""
        # Direct sequence input
        returncode, stdout, stderr = self.run_sim([
            "--seq", "ATCGAATTCGGGATCCAAA",
            "--enz", "EcoRI"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "DIGESTION RESULTS" in stdout, "Should show digestion results"
        assert "EcoRI" in stdout, "Should mention EcoRI"
        assert "fragments" in stdout.lower(), "Should mention fragments"
        
        print("✓ Basic linear digest works")
    
    def test_multiple_enzymes(self):
        """Test multiple enzyme digestion"""
        returncode, stdout, stderr = self.run_sim([
            "--seq", "ATCGAATTCGGGATCCAAAGCTT",
            "--enz", "EcoRI", "BamHI", "HindIII"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "EcoRI" in stdout, "Should show EcoRI"
        assert "BamHI" in stdout, "Should show BamHI"
        assert "HindIII" in stdout, "Should show HindIII"
        
        print("✓ Multiple enzyme digest works")
    
    def test_circular_dna(self):
        """Test circular DNA mode"""
        # Create test file
        test_file = os.path.join(self.temp_dir, "plasmid.fasta")
        with open(test_file, 'w') as f:
            f.write(">test_plasmid\n")
            f.write("ATCGATCGATCGAATTCGGGGGGGATCCATATATATAT" * 20 + "\n")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_file,
            "--enz", "EcoRI", "BamHI",
            "--circular"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "circular" in stdout.lower(), "Should indicate circular topology"
        assert "wraps" in stdout.lower() or "Wraps" in stdout, "Should mention wrapping"
        
        print("✓ Circular DNA mode works")
    
    def test_sequence_extraction(self):
        """Test fragment sequence extraction"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--include-seqs"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "FRAGMENT SEQUENCES" in stdout, "Should show fragment sequences section"
        assert "seq:" in stdout, "Should show sequence data"
        
        print("✓ Sequence extraction works")
    
    def test_fasta_export(self):
        """Test FASTA file export"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        fasta_out = os.path.join(self.temp_dir, "fragments.fasta")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--fasta-out", fasta_out
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert os.path.exists(fasta_out), "FASTA file should be created"
        
        with open(fasta_out, 'r') as f:
            content = f.read()
            assert ">" in content, "Should have FASTA headers"
            assert "frag_" in content, "Should have fragment IDs"
        
        print("✓ FASTA export works")
    
    def test_restriction_map(self):
        """Test restriction map generation"""
        test_seq = "A" * 50 + "GAATTC" + "G" * 50 + "GGATCC" + "T" * 50
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--print-map"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "RESTRICTION MAP" in stdout, "Should show restriction map"
        assert "^" in stdout or "|" in stdout, "Should show cut markers"
        
        print("✓ Restriction map works")
    
    def test_gel_simulation(self):
        """Test gel simulation"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--simulate-gel"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "GEL SIMULATION" in stdout or "AGAROSE GEL" in stdout, "Should show gel section"
        assert "bp" in stdout, "Should show fragment sizes"
        
        print("✓ Gel simulation works")
    
    def test_compatibility_analysis(self):
        """Test ligation compatibility analysis"""
        test_seq = "AAAAGAATTCGGGGGATCCTTTTCTGCAGCCCCC"
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI", "PstI",
            "--compatibility"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "COMPATIBILITY" in stdout, "Should show compatibility analysis"
        
        print("✓ Compatibility analysis works")
    
    def test_genbank_export(self):
        """Test GenBank export"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        gb_out = os.path.join(self.temp_dir, "digest.gb")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--export-genbank", gb_out
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert os.path.exists(gb_out), "GenBank file should be created"
        
        with open(gb_out, 'r') as f:
            content = f.read()
            assert "LOCUS" in content, "Should have LOCUS line"
            assert "ORIGIN" in content, "Should have ORIGIN section"
        
        print("✓ GenBank export works")
    
    def test_csv_export(self):
        """Test CSV export"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        csv_prefix = os.path.join(self.temp_dir, "digest")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--export-csv", csv_prefix
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        
        frag_csv = csv_prefix + "_fragments.csv"
        cuts_csv = csv_prefix + "_cuts.csv"
        
        assert os.path.exists(frag_csv), "Fragments CSV should be created"
        assert os.path.exists(cuts_csv), "Cuts CSV should be created"
        
        with open(frag_csv, 'r') as f:
            content = f.read()
            assert "fragment_id" in content, "Should have headers"
        
        print("✓ CSV export works")
    
    def test_svg_graphics(self):
        """Test SVG graphics generation"""
        test_seq = "A" * 100 + "GAATTC" + "G" * 100
        svg_out = os.path.join(self.temp_dir, "map.svg")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI",
            "--out-svg", svg_out
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert os.path.exists(svg_out), "SVG file should be created"
        
        with open(svg_out, 'r') as f:
            content = f.read()
            assert "<svg" in content, "Should be valid SVG"
            assert "EcoRI" in content, "Should contain enzyme label"
        
        print("✓ SVG graphics works")
    
    def test_theoretical_compatibility(self):
        """Test theoretical enzyme compatibility"""
        returncode, stdout, stderr = self.run_sim([
            "--theoretical-enzymes", "EcoRI,BamHI,HindIII"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "THEORETICAL" in stdout or "theoretical" in stdout, "Should show theoretical analysis"
        
        print("✓ Theoretical compatibility works")
    
    def test_case_insensitive_enzymes(self):
        """Test case-insensitive enzyme names"""
        returncode, stdout, stderr = self.run_sim([
            "--seq", "ATCGAATTCGGGATCCAAA",
            "--enz", "ecori", "bamhi"  # lowercase
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "EcoRI" in stdout, "Should normalize enzyme names"
        
        print("✓ Case-insensitive enzymes work")
    
    def test_combined_output_formats(self):
        """Test combining multiple output formats"""
        test_seq = "AAAGAATTCGGGGGATCCTTTT"
        gb_out = os.path.join(self.temp_dir, "combined.gb")
        csv_prefix = os.path.join(self.temp_dir, "combined")
        svg_out = os.path.join(self.temp_dir, "combined.svg")
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--enz", "EcoRI", "BamHI",
            "--export-genbank", gb_out,
            "--export-csv", csv_prefix,
            "--out-svg", svg_out,
            "--print-map"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert os.path.exists(gb_out), "GenBank file should exist"
        assert os.path.exists(csv_prefix + "_fragments.csv"), "Fragments CSV should exist"
        assert os.path.exists(svg_out), "SVG file should exist"
        assert "RESTRICTION MAP" in stdout, "Should show restriction map"
        
        print("✓ Combined output formats work")
    
    def test_lanes_config_gel(self):
        """Test multi-lane gel configuration"""
        test_seq = "ATCGAATTCGGGATCCAAAGCTT"
        lanes_json = [
            {"label": "Control", "enzymes": []},
            {"label": "EcoRI", "enzymes": ["EcoRI"]},
            {"label": "Double", "enzymes": ["EcoRI", "BamHI"]}
        ]
        
        lanes_file = os.path.join(self.temp_dir, "lanes.json")
        with open(lanes_file, 'w') as f:
            json.dump(lanes_json, f)
        
        returncode, stdout, stderr = self.run_sim([
            "--seq", test_seq,
            "--lanes-config", lanes_file,
            "--gel-only"
        ])
        
        assert returncode == 0, f"Command failed: {stderr}"
        assert "Control" in stdout or "EcoRI" in stdout or "Double" in stdout, \
            "Should show lane labels"
        
        print("✓ Multi-lane gel configuration works")
    
    def test_error_handling_invalid_enzyme(self):
        """Test error handling for invalid enzyme"""
        returncode, stdout, stderr = self.run_sim([
            "--seq", "ATCGATCGATCG",
            "--enz", "InvalidEnzymeXYZ123"
        ])
        
        assert returncode != 0, "Should fail with invalid enzyme"
        output = stdout + stderr
        assert "not found" in output.lower() or "error" in output.lower(), \
            "Should show error message"
        
        print("✓ Error handling for invalid enzyme works")
    
    def test_no_cuts_found(self):
        """Test handling when no cuts are found"""
        returncode, stdout, stderr = self.run_sim([
            "--seq", "ATCGATCGATCGATCGATCG",
            "--enz", "NotI"  # 8bp cutter unlikely in short sequence
        ])
        
        assert returncode == 0, "Should succeed even with no cuts"
        assert "No cut sites" in stdout or "0" in stdout, "Should indicate no cuts"
        
        print("✓ No cuts found handling works")


def run_cli_tests():
    """Run all CLI tests"""
    test_class = TestCLI()
    
    test_methods = [method for method in dir(test_class) 
                   if method.startswith('test_') and callable(getattr(test_class, method))]
    
    print("=" * 80)
    print("COMMAND-LINE INTERFACE TESTING")
    print("=" * 80)
    print(f"Running {len(test_methods)} CLI test scenarios...\n")
    
    passed = 0
    failed = 0
    errors = []
    
    for test_name in test_methods:
        print(f"\n{'─' * 80}")
        print(f"Test: {test_name}")
        print('─' * 80)
        
        test_class.setup_method()
        
        try:
            test_method = getattr(test_class, test_name)
            test_method()
            passed += 1
            print(f"✓ {test_name} PASSED")
        except AssertionError as e:
            failed += 1
            error_msg = f"{test_name}: {str(e)}"
            errors.append(error_msg)
            print(f"✗ {test_name} FAILED: {e}")
        except Exception as e:
            failed += 1
            error_msg = f"{test_name}: {str(e)}"
            errors.append(error_msg)
            print(f"✗ {test_name} ERROR: {e}")
            import traceback
            traceback.print_exc()
        finally:
            test_class.teardown_method()
    
    print("\n" + "=" * 80)
    print("CLI TEST SUMMARY")
    print("=" * 80)
    print(f"Total tests: {len(test_methods)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Success rate: {passed / len(test_methods) * 100:.1f}%")
    
    if errors:
        print("\nFailed tests:")
        for error in errors:
            print(f"  ✗ {error}")
    
    return passed, failed


if __name__ == "__main__":
    passed, failed = run_cli_tests()
    sys.exit(0 if failed == 0 else 1)

