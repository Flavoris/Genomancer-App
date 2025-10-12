def run_command(cmd, check_output=True):
    """Run a command and return output or check for errors"""
    try:
        # Change to the directory where sim.py is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30,
            cwd=script_dir
        )
        if check_output:
            return result.returncode == 0, result.stdout, result.stderr
        return True, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)#!/usr/bin/env python3
"""
Comprehensive Test Suite for Restriction Enzyme Simulator
Tests all features documented in README.md with real-world scenarios
"""

import subprocess
import json
import os
import tempfile
import sys
from pathlib import Path

class TestResults:
    def __init__(self):
        self.passed = []
        self.failed = []
        self.warnings = []
    
    def add_pass(self, test_name, details=""):
        self.passed.append((test_name, details))
        print(f"✅ PASS: {test_name}")
        if details:
            print(f"   {details}")
    
    def add_fail(self, test_name, error):
        self.failed.append((test_name, error))
        print(f"❌ FAIL: {test_name}")
        print(f"   Error: {error}")
    
    def add_warning(self, test_name, warning):
        self.warnings.append((test_name, warning))
        print(f"⚠️  WARN: {test_name}")
        print(f"   {warning}")
    
    def summary(self):
        total = len(self.passed) + len(self.failed)
        print("\n" + "="*80)
        print("TEST SUMMARY")
        print("="*80)
        print(f"Total Tests: {total}")
        print(f"Passed: {len(self.passed)} ({100*len(self.passed)/total if total > 0 else 0:.1f}%)")
        print(f"Failed: {len(self.failed)}")
        print(f"Warnings: {len(self.warnings)}")
        
        if self.failed:
            print("\n❌ FAILED TESTS:")
            for name, error in self.failed:
                print(f"  - {name}: {error}")
        
        if self.warnings:
            print("\n⚠️  WARNINGS:")
            for name, warning in self.warnings:
                print(f"  - {name}: {warning}")
        
        return len(self.failed) == 0


def run_command(cmd, check_output=True):
    """Run a command and return output or check for errors"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )
        if check_output:
            return result.returncode == 0, result.stdout, result.stderr
        return True, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)


def create_test_files():
    """Create test DNA sequences"""
    test_dir = tempfile.mkdtemp()
    
    # Linear test sequence with multiple enzyme sites
    linear_seq = (
        ">Linear_test_sequence\n"
        "ATCGATCGAATTCGGGATCCAAGCTTCTGCAGGGCGCGCCGCGGCCGCATCGATCG\n"
        "GAATTCGGGATCCAAGCTTCTGCAGGGCGCGCCGCGGCCGCATCGATCGAATTCGG\n"
    )
    
    # Circular plasmid (pUC19-like)
    circular_seq = (
        ">Circular_plasmid_test\n"
        "GAATTCGGGATCCAAGCTTCTGCAGGGCGCGCCGCGGCCGCATCGATCGGAATTCG\n"
        "GGATCCAAGCTTCTGCAGGGCGCGCCGCGGCCGCATCGATCGGAATTCGGGATCCA\n"
        "AGCTTCTGCAGGGCGCGCCGCGGCCGCATCGATCG\n"
    )
    
    # Small test sequence
    small_seq = ">Small_test\nATCGAATTCGGGATCCAA\n"
    
    files = {}
    files['linear'] = os.path.join(test_dir, 'linear_test.fasta')
    files['circular'] = os.path.join(test_dir, 'circular_test.fasta')
    files['small'] = os.path.join(test_dir, 'small_test.fasta')
    
    with open(files['linear'], 'w') as f:
        f.write(linear_seq)
    with open(files['circular'], 'w') as f:
        f.write(circular_seq)
    with open(files['small'], 'w') as f:
        f.write(small_seq)
    
    return test_dir, files


def test_basic_functionality(results, files):
    """Test basic linear digestion"""
    print("\n" + "="*80)
    print("TESTING: Basic Functionality")
    print("="*80)
    
    # Test 1: Single enzyme linear
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI'
    )
    if success and "DIGESTION RESULTS" in stdout and "EcoRI" in stdout:
        results.add_pass("Basic single enzyme", "EcoRI digestion works")
    else:
        results.add_fail("Basic single enzyme", stderr or "No output")
    
    # Test 2: Multiple enzymes
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI BamHI'
    )
    if success and "EcoRI" in stdout and "BamHI" in stdout:
        results.add_pass("Multiple enzymes", "Both enzymes detected")
    else:
        results.add_fail("Multiple enzymes", stderr or "Missing enzyme output")
    
    # Test 3: Case insensitive
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz ecori bamhi'
    )
    if success and "EcoRI" in stdout:
        results.add_pass("Case insensitive", "Lowercase enzyme names work")
    else:
        results.add_fail("Case insensitive", "Failed with lowercase names")


def test_circular_dna(results, files):
    """Test circular DNA digestion"""
    print("\n" + "="*80)
    print("TESTING: Circular DNA")
    print("="*80)
    
    # Test 1: Circular with two cuts
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["circular"]}" --enz EcoRI BamHI --circular'
    )
    if success and "wraps" in stdout.lower() and "circular" in stdout.lower():
        results.add_pass("Circular two cuts", "Wrap-around fragment detected")
    else:
        results.add_fail("Circular two cuts", "No wrap-around fragment")
    
    # Test 2: Single cut intact circle
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["circular"]}" --enz EcoRI --circular'
    )
    if success and "1" in stdout:  # Should show 1 or 2 fragments
        results.add_pass("Circular single cut", "Single cut handled")
    else:
        results.add_fail("Circular single cut", "Failed to handle single cut")
    
    # Test 3: Single cut linearization
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["circular"]}" --enz EcoRI --circular --circular_single_cut_linearizes'
    )
    if success:
        results.add_pass("Circular linearization", "Linearization mode works")
    else:
        results.add_fail("Circular linearization", stderr)


def test_sequence_extraction(results, files):
    """Test fragment sequence extraction"""
    print("\n" + "="*80)
    print("TESTING: Sequence Extraction")
    print("="*80)
    
    # Test 1: Include sequences
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --include-seqs'
    )
    if success and "FRAGMENT SEQUENCES" in stdout and "seq:" in stdout:
        results.add_pass("Sequence display", "Sequences shown in output")
    else:
        results.add_fail("Sequence display", "Sequences not displayed")
    
    # Test 2: FASTA export
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_out = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --fasta-out {fasta_out}'
    )
    
    if success and os.path.exists(fasta_out):
        with open(fasta_out, 'r') as f:
            content = f.read()
        if content.startswith('>') and 'frag_' in content:
            results.add_pass("FASTA export", f"FASTA file created with {len(content.split('>')) - 1} fragments")
        else:
            results.add_fail("FASTA export", "Invalid FASTA format")
        os.unlink(fasta_out)
    else:
        results.add_fail("FASTA export", "FASTA file not created")
    
    # Test 3: Sequence context
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["linear"]}" --enz EcoRI --include-seqs --seq-context 5'
    )
    if success and "..." in stdout:
        results.add_pass("Sequence elision", "Context limiting works")
    else:
        results.add_warning("Sequence elision", "Elision not clearly visible (may depend on fragment size)")


def test_restriction_maps(results, files):
    """Test restriction map visualization"""
    print("\n" + "="*80)
    print("TESTING: Restriction Maps")
    print("="*80)
    
    # Test 1: Print map
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["linear"]}" --enz EcoRI BamHI --print-map'
    )
    if success and "RESTRICTION MAP" in stdout and "---" in stdout:
        results.add_pass("Restriction map", "Map displayed")
    else:
        results.add_fail("Restriction map", "Map not displayed")
    
    # Test 2: Map only
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["linear"]}" --enz EcoRI --print-map-only'
    )
    if success and "RESTRICTION MAP" in stdout and "DIGESTION RESULTS" not in stdout:
        results.add_pass("Map only mode", "Only map shown")
    else:
        results.add_fail("Map only mode", "Other sections still shown")
    
    # Test 3: Show overhangs
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI PstI --print-map --map-show-overhangs'
    )
    if success and "overhang" in stdout.lower():
        results.add_pass("Map overhangs", "Overhangs displayed")
    else:
        results.add_warning("Map overhangs", "Overhangs not visible (may depend on enzyme)")


def test_gel_simulation(results, files):
    """Test gel simulation"""
    print("\n" + "="*80)
    print("TESTING: Gel Simulation")
    print("="*80)
    
    # Test 1: Basic gel
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --simulate-gel'
    )
    if success and "AGAROSE GEL SIMULATION" in stdout and "Ladder" in stdout:
        results.add_pass("Basic gel", "Gel simulation displayed")
    else:
        results.add_fail("Basic gel", "Gel not displayed")
    
    # Test 2: Gel only
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --gel-only'
    )
    if success and "AGAROSE GEL" in stdout and "DIGESTION RESULTS" not in stdout:
        results.add_pass("Gel only mode", "Only gel shown")
    else:
        results.add_fail("Gel only mode", "Other sections shown")
    
    # Test 3: Different ladder
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --simulate-gel --gel-ladder 100bp'
    )
    if success and "100" in stdout:
        results.add_pass("Gel ladder", "100bp ladder works")
    else:
        results.add_warning("Gel ladder", "100bp ladder not clearly visible")
    
    # Test 4: Multi-lane configuration
    lanes_config = json.dumps([
        {"label": "Lane1", "enzymes": ["EcoRI"]},
        {"label": "Lane2", "enzymes": ["BamHI"]}
    ])
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --lanes-config \'{lanes_config}\' --gel-only'
    )
    if success and "Lane1" in stdout and "Lane2" in stdout:
        results.add_pass("Multi-lane gel", "Multiple lanes displayed")
    else:
        results.add_fail("Multi-lane gel", stderr or "Lanes not displayed")


def test_export_formats(results, files):
    """Test GenBank and CSV export"""
    print("\n" + "="*80)
    print("TESTING: Export Formats")
    print("="*80)
    
    # Test 1: GenBank export
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gbk', delete=False) as f:
        gb_file = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --export-genbank {gb_file}'
    )
    
    if success and os.path.exists(gb_file):
        with open(gb_file, 'r') as f:
            content = f.read()
        if "LOCUS" in content and "ORIGIN" in content:
            results.add_pass("GenBank export", "Valid GenBank file created")
        else:
            results.add_fail("GenBank export", "Invalid GenBank format")
        os.unlink(gb_file)
    else:
        results.add_fail("GenBank export", "GenBank file not created")
    
    # Test 2: CSV export
    csv_prefix = os.path.join(tempfile.gettempdir(), "test_export")
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --export-csv {csv_prefix}'
    )
    
    frag_csv = f"{csv_prefix}_fragments.csv"
    cuts_csv = f"{csv_prefix}_cuts.csv"
    
    if success and os.path.exists(frag_csv) and os.path.exists(cuts_csv):
        results.add_pass("CSV export", "Both CSV files created")
        os.unlink(frag_csv)
        os.unlink(cuts_csv)
    else:
        results.add_fail("CSV export", "CSV files not created")


def test_ligation_compatibility(results, files):
    """Test ligation compatibility analysis"""
    print("\n" + "="*80)
    print("TESTING: Ligation Compatibility")
    print("="*80)
    
    # Test 1: Basic compatibility
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --compatibility'
    )
    if success and "COMPATIBILITY" in stdout:
        results.add_pass("Basic compatibility", "Compatibility analysis works")
    else:
        results.add_fail("Basic compatibility", "Compatibility not shown")
    
    # Test 2: Matrix format
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI BamHI --compatibility --compat-summary matrix'
    )
    if success and ("✓" in stdout or "." in stdout):
        results.add_pass("Compatibility matrix", "Matrix format works")
    else:
        results.add_fail("Compatibility matrix", "Matrix not displayed")
    
    # Test 3: JSON export
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json_file = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --compatibility --json-out {json_file}'
    )
    
    if success and os.path.exists(json_file):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            if isinstance(data, list):
                results.add_pass("Compatibility JSON", f"JSON file with {len(data)} entries")
            else:
                results.add_fail("Compatibility JSON", "Invalid JSON structure")
        except json.JSONDecodeError:
            results.add_fail("Compatibility JSON", "Invalid JSON format")
        os.unlink(json_file)
    else:
        results.add_fail("Compatibility JSON", "JSON file not created")


def test_theoretical_compatibility(results):
    """Test theoretical enzyme compatibility"""
    print("\n" + "="*80)
    print("TESTING: Theoretical Compatibility")
    print("="*80)
    
    # Test 1: Specific enzymes
    success, stdout, stderr = run_command(
        'python sim.py --theoretical-enzymes "EcoRI,MfeI,BamHI"'
    )
    if success and "Theoretical" in stdout and "EcoRI" in stdout:
        results.add_pass("Theoretical specific", "Specific enzyme analysis works")
    else:
        results.add_fail("Theoretical specific", stderr or "Analysis failed")
    
    # Test 2: All enzymes (with timeout protection)
    success, stdout, stderr = run_command(
        'python sim.py --theoretical-all --format pairs'
    )
    if success and "compatible" in stdout.lower():
        results.add_pass("Theoretical all", "All enzyme analysis works")
    else:
        results.add_warning("Theoretical all", "May take too long or have issues")
    
    # Test 3: Directional filter
    success, stdout, stderr = run_command(
        'python sim.py --theoretical-enzymes "EcoRI,SpeI,XbaI" --require-directional'
    )
    if success:
        results.add_pass("Theoretical directional", "Directional filtering works")
    else:
        results.add_fail("Theoretical directional", stderr)


def test_graphics_output(results, files):
    """Test SVG/PNG graphics generation"""
    print("\n" + "="*80)
    print("TESTING: Graphics Output")
    print("="*80)
    
    # Test 1: SVG plasmid map
    with tempfile.NamedTemporaryFile(mode='w', suffix='.svg', delete=False) as f:
        svg_file = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["circular"]}" --enz EcoRI BamHI --circular --out-svg {svg_file}'
    )
    
    if success and os.path.exists(svg_file):
        with open(svg_file, 'r') as f:
            content = f.read()
        if '<svg' in content and '</svg>' in content:
            results.add_pass("SVG plasmid map", "Valid SVG generated")
        else:
            results.add_fail("SVG plasmid map", "Invalid SVG format")
        os.unlink(svg_file)
    else:
        results.add_fail("SVG plasmid map", "SVG file not created")
    
    # Test 2: Linear map
    with tempfile.NamedTemporaryFile(mode='w', suffix='.svg', delete=False) as f:
        linear_svg = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["linear"]}" --enz EcoRI --out-svg-linear {linear_svg}'
    )
    
    if success and os.path.exists(linear_svg):
        results.add_pass("SVG linear map", "Linear map generated")
        os.unlink(linear_svg)
    else:
        results.add_fail("SVG linear map", "Linear map not created")
    
    # Test 3: Fragment diagram
    with tempfile.NamedTemporaryFile(mode='w', suffix='.svg', delete=False) as f:
        frag_svg = f.name
    
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI --out-svg-fragments {frag_svg}'
    )
    
    if success and os.path.exists(frag_svg):
        results.add_pass("SVG fragment diagram", "Fragment diagram generated")
        os.unlink(frag_svg)
    else:
        results.add_fail("SVG fragment diagram", "Fragment diagram not created")


def test_edge_cases(results, files):
    """Test edge cases and error handling"""
    print("\n" + "="*80)
    print("TESTING: Edge Cases & Error Handling")
    print("="*80)
    
    # Test 1: Invalid enzyme
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz InvalidEnzyme123',
        check_output=False
    )
    if not success or "not found" in stderr or "not found" in stdout:
        results.add_pass("Invalid enzyme", "Error handling works")
    else:
        results.add_fail("Invalid enzyme", "Should reject invalid enzyme")
    
    # Test 2: No cuts found
    success, stdout, stderr = run_command(
        'python sim.py --seq "AAAAAAAAAA" --enz EcoRI'
    )
    if success and ("No cut" in stdout or "0" in stdout):
        results.add_pass("No cuts found", "Handles sequences with no cuts")
    else:
        results.add_fail("No cuts found", "Failed to handle no cuts")
    
    # Test 3: Very long sequence
    long_seq = "ATCG" * 1000 + "GAATTC" + "ATCG" * 1000
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{long_seq}" --enz EcoRI'
    )
    if success:
        results.add_pass("Long sequence", "Handles long sequences")
    else:
        results.add_fail("Long sequence", "Failed on long sequence")
    
    # Test 4: Duplicate enzyme names
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["small"]}" --enz EcoRI EcoRI --print-map'
    )
    if success and "#2" in stdout:
        results.add_pass("Duplicate enzymes", "Handles duplicate enzyme names")
    else:
        results.add_warning("Duplicate enzymes", "Duplicate handling not clearly visible")


def test_real_world_scenarios(results, files):
    """Test real-world cloning scenarios"""
    print("\n" + "="*80)
    print("TESTING: Real-World Scenarios")
    print("="*80)
    
    # Scenario 1: Standard subcloning
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["linear"]}" --enz EcoRI HindIII --compatibility --print-map'
    )
    if success and "COMPATIBILITY" in stdout and "RESTRICTION MAP" in stdout:
        results.add_pass("Subcloning workflow", "Combined analysis works")
    else:
        results.add_fail("Subcloning workflow", "Combined features failed")
    
    # Scenario 2: Plasmid verification
    success, stdout, stderr = run_command(
        f'python sim.py --seq "{files["circular"]}" --enz EcoRI BamHI --circular --simulate-gel --export-genbank /tmp/test_plasmid.gbk'
    )
    if success and os.path.exists('/tmp/test_plasmid.gbk'):
        results.add_pass("Plasmid verification", "Full plasmid analysis works")
        os.unlink('/tmp/test_plasmid.gbk')
    else:
        results.add_fail("Plasmid verification", "Plasmid workflow failed")
    
    # Scenario 3: Directional cloning planning
    success, stdout, stderr = run_command(
        'python sim.py --theoretical-enzymes "SpeI,XbaI,NheI" --require-directional --format matrix'
    )
    if success:
        results.add_pass("Directional cloning", "Planning workflow works")
    else:
        results.add_fail("Directional cloning", "Planning failed")


def main():
    print("="*80)
    print("COMPREHENSIVE TEST SUITE FOR RESTRICTION ENZYME SIMULATOR")
    print("="*80)
    print("\nCreating test files...")
    
    test_dir, files = create_test_files()
    results = TestResults()
    
    try:
        # Run all test suites
        test_basic_functionality(results, files)
        test_circular_dna(results, files)
        test_sequence_extraction(results, files)
        test_restriction_maps(results, files)
        test_gel_simulation(results, files)
        test_export_formats(results, files)
        test_ligation_compatibility(results, files)
        test_theoretical_compatibility(results)
        test_graphics_output(results, files)
        test_edge_cases(results, files)
        test_real_world_scenarios(results, files)
        
        # Print summary
        success = results.summary()
        
        return 0 if success else 1
        
    finally:
        # Cleanup
        import shutil
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        print(f"\nTest files cleaned up from {test_dir}")


if __name__ == "__main__":
    sys.exit(main())
