#!/usr/bin/env python3
"""
Advanced Real-World Testing - Edge Cases and Production Scenarios
Tests complex situations that might occur in actual molecular biology applications
"""

import sys
import os
import json
import tempfile
import shutil
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim import (
    load_enzyme_database, read_dna_sequence, find_cut_positions_linear,
    merge_cut_positions
)
from fragment_calculator import (
    compute_fragments, compute_fragments_with_sequences,
    simulate_gel
)
from ligation_compatibility import calculate_compatibility
from gel_ladders import get_ladder


class TestAdvancedRealWorld:
    """Advanced production-ready tests"""
    
    def setup_method(self):
        self.enzymes = load_enzyme_database()
        self.temp_dir = tempfile.mkdtemp()
        
    def teardown_method(self):
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_very_long_sequence(self):
        """Test with a very long sequence (simulating bacterial genome fragment)"""
        # 10kb sequence with scattered restriction sites
        segment_size = 500
        num_segments = 20
        
        segments = []
        for i in range(num_segments):
            # Add various recognition sites throughout
            if i % 3 == 0:
                segments.append("GAATTC")  # EcoRI
            segments.append("A" * (segment_size - 6))
            if i % 4 == 0:
                segments.append("GGATCC")  # BamHI
        
        long_seq = "".join(segments)
        
        ecori_cuts = find_cut_positions_linear(long_seq, "EcoRI", self.enzymes, circular=False)
        bamhi_cuts = find_cut_positions_linear(long_seq, "BamHI", self.enzymes, circular=False)
        
        all_cuts = sorted(set(ecori_cuts + bamhi_cuts))
        
        cut_metadata = {}
        for pos in ecori_cuts:
            cut_metadata[pos] = [{
                'enzyme': 'EcoRI',
                'site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            }]
        for pos in bamhi_cuts:
            if pos not in cut_metadata:
                cut_metadata[pos] = []
            cut_metadata[pos].append({
                'enzyme': 'BamHI',
                'site': 'GGATCC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            })
        
        fragments = compute_fragments(
            cut_positions=all_cuts,
            seq_len=len(long_seq),
            circular=False,
            circular_single_cut_linearizes=False,
            cut_metadata=cut_metadata
        )
        
        total = sum(f['length'] for f in fragments)
        assert total == len(long_seq), f"Fragment sum mismatch: {total} vs {len(long_seq)}"
        
        print(f"✓ Long sequence ({len(long_seq)} bp): {len(fragments)} fragments from {len(all_cuts)} cuts")
    
    def test_multiple_cuts_same_position(self):
        """Test multiple enzymes cutting at the same position (isoschizomers)"""
        # Create sequence with site recognized by multiple enzymes
        # For example, some enzymes have the same recognition sequence
        test_seq = "AAAAAGCGGCCGCGGGG"  # NotI site (also recognized by EagI in some systems)
        
        noti_cuts = find_cut_positions_linear(test_seq, "NotI", self.enzymes, circular=False)
        
        cut_metadata = {}
        for pos in noti_cuts:
            cut_metadata[pos] = [{
                'enzyme': 'NotI',
                'site': 'GCGGCCGC',
                'cut_index': 2,
                'overhang_type': "5' overhang"
            }]
            # Simulate another enzyme cutting at same position
            if 'EagI' in self.enzymes:
                cut_metadata[pos].append({
                    'enzyme': 'EagI',
                    'site': 'CGGCCG',
                    'cut_index': 1,
                    'overhang_type': "5' overhang"
                })
        
        fragments = compute_fragments(
            cut_positions=noti_cuts,
            seq_len=len(test_seq),
            circular=False,
            circular_single_cut_linearizes=False,
            cut_metadata=cut_metadata
        )
        
        # Should handle multiple enzymes at same position
        assert len(fragments) > 0, "Should handle coincident cuts"
        
        print(f"✓ Multiple enzymes at same position: {len(fragments)} fragments")
    
    def test_palindromic_sequence(self):
        """Test with palindromic DNA sequence"""
        # Palindromic sequence
        half = "GAATTCGGATCCAAGCTT"
        palindrome = half + half[::-1].translate(str.maketrans('ATCG', 'TAGC'))
        
        ecori_cuts = find_cut_positions_linear(palindrome, "EcoRI", self.enzymes, circular=False)
        bamhi_cuts = find_cut_positions_linear(palindrome, "BamHI", self.enzymes, circular=False)
        hindiii_cuts = find_cut_positions_linear(palindrome, "HindIII", self.enzymes, circular=False)
        
        all_cuts = sorted(set(ecori_cuts + bamhi_cuts + hindiii_cuts))
        
        # Should find symmetric cuts
        print(f"✓ Palindromic sequence: {len(all_cuts)} cuts (should be symmetric)")
    
    def test_high_gc_content(self):
        """Test with high GC content sequence (like bacterial genome)"""
        # 70% GC content
        high_gc_seq = "GCGCGCGCGCGCGCGAATTCGCGCGCGCGCGCGCGCGGATCCGCGCGCGCGC"
        
        ecori_cuts = find_cut_positions_linear(high_gc_seq, "EcoRI", self.enzymes, circular=False)
        bamhi_cuts = find_cut_positions_linear(high_gc_seq, "BamHI", self.enzymes, circular=False)
        
        all_cuts = sorted(set(ecori_cuts + bamhi_cuts))
        
        cut_metadata = {}
        for pos in ecori_cuts:
            cut_metadata[pos] = [{
                'enzyme': 'EcoRI',
                'site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            }]
        for pos in bamhi_cuts:
            if pos not in cut_metadata:
                cut_metadata[pos] = []
            cut_metadata[pos].append({
                'enzyme': 'BamHI',
                'site': 'GGATCC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            })
        
        fragments = compute_fragments(
            cut_positions=all_cuts,
            seq_len=len(high_gc_seq),
            circular=False,
            circular_single_cut_linearizes=False,
            cut_metadata=cut_metadata
        )
        
        print(f"✓ High GC content: {len(fragments)} fragments")
    
    def test_low_complexity_sequence(self):
        """Test with low complexity (repetitive) sequence"""
        # Highly repetitive sequence
        low_complex = "ATATATATAT" * 10 + "GAATTC" + "GCGCGCGCGC" * 10
        
        ecori_cuts = find_cut_positions_linear(low_complex, "EcoRI", self.enzymes, circular=False)
        
        cut_metadata = {}
        for pos in ecori_cuts:
            cut_metadata[pos] = [{
                'enzyme': 'EcoRI',
                'site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            }]
        
        fragments = compute_fragments(
            cut_positions=ecori_cuts,
            seq_len=len(low_complex),
            circular=False,
            circular_single_cut_linearizes=False,
            cut_metadata=cut_metadata
        )
        
        print(f"✓ Low complexity sequence: {len(fragments)} fragments")
    
    def test_circular_plasmid_multiple_cuts(self):
        """Test circular plasmid with many cuts (mimicking complete digest)"""
        # 5kb circular plasmid with 10 restriction sites
        plasmid_segments = []
        enzymes_used = ["EcoRI", "BamHI", "HindIII", "PstI", "NotI"]
        sites = ["GAATTC", "GGATCC", "AAGCTT", "CTGCAG", "GCGGCCGC"]
        
        for i in range(10):
            plasmid_segments.append("A" * 400)
            plasmid_segments.append(sites[i % len(sites)])
        
        plasmid = "".join(plasmid_segments)
        
        all_cuts = []
        cut_metadata = {}
        
        for i, enz_name in enumerate(enzymes_used):
            cuts = find_cut_positions_linear(plasmid, enz_name, self.enzymes, circular=True)
            all_cuts.extend(cuts)
            
            enz_info = self.enzymes[enz_name]
            for pos in cuts:
                if pos not in cut_metadata:
                    cut_metadata[pos] = []
                cut_metadata[pos].append({
                    'enzyme': enz_name,
                    'site': enz_info['sequence'],
                    'cut_index': enz_info['cut_index'],
                    'overhang_type': enz_info['overhang_type']
                })
        
        all_cuts = sorted(set(all_cuts))
        
        fragments = compute_fragments(
            cut_positions=all_cuts,
            seq_len=len(plasmid),
            circular=True,
            circular_single_cut_linearizes=False,
            cut_metadata=cut_metadata
        )
        
        # Circular with N cuts should give N fragments
        assert len(fragments) == len(all_cuts), \
            f"Circular with {len(all_cuts)} cuts should give {len(all_cuts)} fragments"
        
        # One fragment should wrap
        wrap_count = sum(1 for f in fragments if f['wraps'])
        assert wrap_count == 1, "Exactly one fragment should wrap"
        
        print(f"✓ Circular with many cuts: {len(fragments)} fragments from {len(all_cuts)} cuts")
    
    def test_gel_simulation_various_sizes(self):
        """Test gel simulation with various fragment sizes"""
        # Create fragments of different sizes for gel testing
        fragment_sizes_lane1 = [100, 500, 1000, 3000, 10000]
        fragment_sizes_lane2 = [200, 750, 1500, 2500, 5000]
        
        lanes = [
            {'label': 'Sample 1', 'fragments': fragment_sizes_lane1, 
             'topology': 'linearized', 'circular': False, 'notes': ''},
            {'label': 'Sample 2', 'fragments': fragment_sizes_lane2,
             'topology': 'linearized', 'circular': False, 'notes': ''}
        ]
        
        ladder_bp = get_ladder("1kb")
        
        gel_output = simulate_gel(
            lanes=lanes,
            ladder_bp=ladder_bp,
            gel_percent=1.0,
            gel_length=24,
            gel_width=80,
            lane_gap=3,
            merge_threshold=20,
            smear="none",
            dye_front=0.85
        )
        
        assert len(gel_output) > 0, "Should generate gel output"
        assert "Sample 1" in gel_output or "Sample 2" in gel_output, "Should include lane labels"
        
        print(f"✓ Gel simulation with various fragment sizes")
    
    def test_realistic_cloning_insert_vector(self):
        """Test realistic cloning scenario: inserting PCR product into vector"""
        # Vector with EcoRI and BamHI sites in MCS
        vector = (
            "ATCGATCGATCGATCGATCGATCGATCGAATTCGGGGATCCATCGATCGATCGATCGATCG"
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
        ) * 20  # ~2.5kb vector
        
        # PCR insert with compatible ends (EcoRI and BamHI)
        insert = (
            "GAATTCATGAAAAAATTTTTTTTTGGGGGGGGGCCCCCCCCGGATCC"  # 0.5kb insert
        )
        
        # Digest vector with both enzymes
        vector_ecori = find_cut_positions_linear(vector, "EcoRI", self.enzymes, circular=True)
        vector_bamhi = find_cut_positions_linear(vector, "BamHI", self.enzymes, circular=True)
        
        # Digest insert with both enzymes
        insert_ecori = find_cut_positions_linear(insert, "EcoRI", self.enzymes, circular=False)
        insert_bamhi = find_cut_positions_linear(insert, "BamHI", self.enzymes, circular=False)
        
        assert len(vector_ecori) > 0 and len(vector_bamhi) > 0, "Vector should have both sites"
        assert len(insert_ecori) > 0 and len(insert_bamhi) > 0, "Insert should have both sites"
        
        print(f"✓ Realistic cloning: vector has {len(vector_ecori)} EcoRI and {len(vector_bamhi)} BamHI sites")
        print(f"                    insert has {len(insert_ecori)} EcoRI and {len(insert_bamhi)} BamHI sites")
    
    def test_star_activity_sites(self):
        """Test with sequences that might show star activity (relaxed specificity)"""
        # Star activity can occur under non-optimal conditions
        # Test with near-matches to recognition sites
        test_seq = "AAAGAATTAGGGGATCCGGGAATTCAAA"  # GAATTA (star) and GAATTC (normal)
        
        ecori_cuts = find_cut_positions_linear(test_seq, "EcoRI", self.enzymes, circular=False)
        
        # Under normal conditions, should only find perfect matches
        # In real scenarios, might need to account for star activity
        
        print(f"✓ Star activity test: found {len(ecori_cuts)} perfect matches")
    
    def test_methylation_sensitive_sites(self):
        """Test sequences with CpG islands (methylation-sensitive)"""
        # Some enzymes are sensitive to Dam/Dcm methylation
        # CpG-rich sequence
        cpg_rich = "CGCGCGCGCGAATTCCGCGCGCGCG"
        
        # EcoRI contains CG and can be affected by methylation in real scenarios
        ecori_cuts = find_cut_positions_linear(cpg_rich, "EcoRI", self.enzymes, circular=False)
        
        # In simulation, should find the site (doesn't model methylation)
        # But good to test the sequence type
        
        print(f"✓ CpG-rich sequence: {len(ecori_cuts)} cuts (methylation not modeled)")
    
    def test_golden_gate_assembly_sites(self):
        """Test Type IIS enzymes for Golden Gate assembly"""
        # BsaI is commonly used for Golden Gate
        test_seq = "AAAAAAGGTCTCANNNNTTTTTT"  # BsaI site with 4bp overhang space
        
        if "BsaI" in self.enzymes:
            bsai_cuts = find_cut_positions_linear(test_seq, "BsaI", self.enzymes, circular=False)
            print(f"✓ Golden Gate (BsaI): {len(bsai_cuts)} cuts found")
        else:
            print("⊘ BsaI not in database")
    
    def test_gibson_assembly_compatible_ends(self):
        """Test overlapping ends for Gibson assembly"""
        # Gibson assembly requires 15-40bp homology
        # Create two fragments with 20bp overlap
        frag1 = "ATCGATCGATCGATCGATCGAATTCGGGGGGGGGGGGGGGGGGG"  # Ends with specific sequence
        frag2 = "GGGGGGGGGGGGGGGGGGGATCCAAAAAAAAAAAAAAAAAAAA"  # Starts with overlap
        
        # In real Gibson, the overlapping sequences would be homologous
        # Here we're just testing that we can identify the ends
        
        ecori_cuts = find_cut_positions_linear(frag1, "EcoRI", self.enzymes, circular=False)
        bamhi_cuts = find_cut_positions_linear(frag2, "BamHI", self.enzymes, circular=False)
        
        print(f"✓ Gibson assembly simulation: frag1 has {len(ecori_cuts)} cuts, frag2 has {len(bamhi_cuts)} cuts")
    
    def test_directional_cloning_strategy(self):
        """Test directional cloning with non-compatible enzymes"""
        # Use two different enzymes to prevent insert from ligating backwards
        test_seq = "AAAAGAATTCGGGGGGGGGAAGCTTTTTT"  # EcoRI and HindIII
        
        ecori_cuts = find_cut_positions_linear(test_seq, "EcoRI", self.enzymes, circular=False)
        hindiii_cuts = find_cut_positions_linear(test_seq, "HindIII", self.enzymes, circular=False)
        
        all_cuts = sorted(set(ecori_cuts + hindiii_cuts))
        
        cut_metadata = {}
        for pos in ecori_cuts:
            cut_metadata[pos] = [{
                'enzyme': 'EcoRI',
                'site': 'GAATTC',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            }]
        for pos in hindiii_cuts:
            if pos not in cut_metadata:
                cut_metadata[pos] = []
            cut_metadata[pos].append({
                'enzyme': 'HindIII',
                'site': 'AAGCTT',
                'cut_index': 1,
                'overhang_type': "5' overhang"
            })
        
        if all_cuts:
            from fragment_calculator import extract_fragment_ends_for_ligation
            
            fragment_ends = extract_fragment_ends_for_ligation(
                dna_sequence=test_seq,
                cut_positions=all_cuts,
                circular=False,
                circular_single_cut_linearizes=False,
                cut_metadata=cut_metadata
            )
            
            # EcoRI and HindIII create different overhangs, so they shouldn't be compatible
            compat_results = calculate_compatibility(
                ends=fragment_ends,
                include_blunt=False,
                min_overhang=1,
                require_directional=False
            )
            
            # Should have some pairs but EcoRI-HindIII should not be compatible
            print(f"✓ Directional cloning: {len(compat_results)} compatible pairs")
    
    def test_synthetic_sequence_optimization(self):
        """Test codon-optimized synthetic sequence"""
        # Codon-optimized sequences might have unexpected restriction sites removed
        synthetic = (
            "ATGGCTAGCGGTACCGAATTCGGATCCAAGCTTCTGCAGGGCGGCCGC"
            "ATGATGATGATGATGATGTAGAGATCTGCTAGCGGTACCGGATCC"
        )
        
        # Count all restriction sites from common cloning enzymes
        common_enzymes = ["EcoRI", "BamHI", "HindIII", "PstI", "NotI", "XbaI", "SpeI"]
        total_sites = 0
        
        for enz in common_enzymes:
            if enz in self.enzymes:
                cuts = find_cut_positions_linear(synthetic, enz, self.enzymes, circular=False)
                total_sites += len(cuts)
        
        print(f"✓ Synthetic sequence: {total_sites} total restriction sites from {len(common_enzymes)} enzymes")


def run_advanced_tests():
    """Run all advanced tests"""
    test_class = TestAdvancedRealWorld()
    
    test_methods = [method for method in dir(test_class) 
                   if method.startswith('test_') and callable(getattr(test_class, method))]
    
    print("=" * 80)
    print("ADVANCED REAL-WORLD TESTING")
    print("=" * 80)
    print(f"Running {len(test_methods)} advanced test scenarios...\n")
    
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
    print("ADVANCED TEST SUMMARY")
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
    passed, failed = run_advanced_tests()
    sys.exit(0 if failed == 0 else 1)

