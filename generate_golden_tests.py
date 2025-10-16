#!/usr/bin/env python3
"""
Generate golden test outputs for parity testing between Python and Swift implementations.
Task 12.1: Generate test cases for 2 linear sequences × 3 enzyme sets and 1 circular plasmid × 3 enzyme sets.
"""

import json
from sim import load_enzyme_database, find_cut_positions_linear, merge_cut_positions
from fragment_calculator import compute_fragments, extract_fragment_ends_for_ligation
from ligation_compatibility import calculate_compatibility

def generate_test_cases():
    """Generate golden test outputs for XCTest parity testing."""
    
    # Load enzyme database
    ENZYMES = load_enzyme_database()
    
    # Define test sequences
    test_sequences = {
        "linear_1": {
            "name": "Linear DNA 1 (500bp)",
            "sequence": "ATCG" * 25 + "GAATTC" + "CGTA" * 25 + "GGATCC" + "TACG" * 25 + "AAGCTT" + "GCTA" * 25,
            "circular": False
        },
        "linear_2": {
            "name": "Linear DNA 2 (400bp)",
            "sequence": "AAAA" * 20 + "CTGCAG" + "TTTT" * 20 + "GCTAGC" + "GGGG" * 20 + "TCTAGA" + "CCCC" * 20,
            "circular": False
        },
        "circular_1": {
            "name": "Circular Plasmid (600bp)",
            "sequence": "ATGC" * 30 + "GAATTC" + "CGAT" * 30 + "GGATCC" + "TACG" * 30 + "CTGCAG" + "GCTA" * 25,
            "circular": True
        }
    }
    
    # Define enzyme sets to test
    enzyme_sets = [
        {
            "name": "Set 1: EcoRI + BamHI",
            "enzymes": ["EcoRI", "BamHI"]
        },
        {
            "name": "Set 2: PstI + NheI",
            "enzymes": ["PstI", "NheI"]
        },
        {
            "name": "Set 3: HindIII + XbaI",
            "enzymes": ["HindIII", "XbaI"]
        }
    ]
    
    golden_tests = []
    
    # Generate test cases
    for seq_id, seq_info in test_sequences.items():
        sequence = seq_info["sequence"]
        seq_name = seq_info["name"]
        circular = seq_info["circular"]
        
        for enzyme_set in enzyme_sets:
            enz_names = enzyme_set["enzymes"]
            set_name = enzyme_set["name"]
            
            print(f"\n{'='*80}")
            print(f"Generating: {seq_name} × {set_name}")
            print(f"{'='*80}")
            
            # Find cut positions for each enzyme
            cuts_by_enzyme = {}
            cut_metadata = {}
            
            for enz_name in enz_names:
                if enz_name not in ENZYMES:
                    print(f"Warning: Enzyme {enz_name} not found in database, skipping")
                    continue
                
                enz_info = ENZYMES[enz_name]
                enz_cuts = find_cut_positions_linear(sequence, enz_name, ENZYMES, circular=circular)
                cuts_by_enzyme[enz_name] = enz_cuts
                
                # Store metadata for each cut position
                for pos in enz_cuts:
                    if pos not in cut_metadata:
                        cut_metadata[pos] = []
                    cut_metadata[pos].append({
                        'enzyme': enz_name,
                        'actual_enzyme': enz_name,
                        'site': enz_info['sequence'],
                        'cut_index': enz_info['cut_index'],
                        'overhang_type': enz_info['overhang_type']
                    })
            
            # Merge all cut positions
            all_cuts = merge_cut_positions(cuts_by_enzyme, len(sequence))
            
            # Compute fragments
            fragments = compute_fragments(
                cut_positions=all_cuts,
                seq_len=len(sequence),
                circular=circular,
                circular_single_cut_linearizes=False,
                cut_metadata=cut_metadata
            )
            
            # Extract fragment lengths and boundary positions
            fragment_lengths = [f['length'] for f in fragments]
            boundary_positions = sorted(all_cuts)
            
            # Extract overhang information from fragments
            overhang_info = []
            for idx, frag in enumerate(fragments):
                left_cut = frag['boundaries']['left_cut']
                right_cut = frag['boundaries']['right_cut']
                
                left_info = {
                    'position': left_cut['pos'] if left_cut else None,
                    'enzymes': [e['enzyme'] for e in left_cut['enzymes']] if left_cut else [],
                    'overhang_types': [e['overhang_type'] for e in left_cut['enzymes']] if left_cut else []
                }
                
                right_info = {
                    'position': right_cut['pos'] if right_cut else None,
                    'enzymes': [e['enzyme'] for e in right_cut['enzymes']] if right_cut else [],
                    'overhang_types': [e['overhang_type'] for e in right_cut['enzymes']] if right_cut else []
                }
                
                overhang_info.append({
                    'fragment_index': idx,
                    'left': left_info,
                    'right': right_info
                })
            
            # Calculate compatibility for 2-3 representative pairs
            compatibility_matrix = []
            if all_cuts:  # Only if there are cuts
                try:
                    fragment_ends = extract_fragment_ends_for_ligation(
                        dna_sequence=sequence,
                        cut_positions=all_cuts,
                        circular=circular,
                        circular_single_cut_linearizes=False,
                        cut_metadata=cut_metadata
                    )
                    
                    if fragment_ends:
                        compat_results = calculate_compatibility(
                            ends=fragment_ends,
                            include_blunt=False,
                            min_overhang=1,
                            require_directional=False
                        )
                        
                        # Take first 3 compatible pairs as representative samples
                        for a, b, directional, reason in compat_results[:3]:
                            compatibility_matrix.append({
                                'end_a': {
                                    'fragment': a.fragment_id,
                                    'side': a.side,
                                    'enzyme': a.enzyme,
                                    'overhang_type': a.overhang_type,
                                    'overhang_seq': a.sticky_bases,
                                    'overhang_len': a.overhang_len
                                },
                                'end_b': {
                                    'fragment': b.fragment_id,
                                    'side': b.side,
                                    'enzyme': b.enzyme,
                                    'overhang_type': b.overhang_type,
                                    'overhang_seq': b.sticky_bases,
                                    'overhang_len': b.overhang_len
                                },
                                'compatible': True,
                                'directional': directional,
                                'reason': reason
                            })
                except Exception as e:
                    print(f"Warning: Could not compute compatibility: {e}")
            
            # Build test case
            test_case = {
                'test_id': f"{seq_id}_{enzyme_set['name'].replace(' ', '_').replace(':', '')}",
                'sequence_name': seq_name,
                'sequence_length': len(sequence),
                'sequence': sequence,
                'circular': circular,
                'enzyme_set': enz_names,
                'cuts_by_enzyme': cuts_by_enzyme,
                'boundary_positions': boundary_positions,
                'fragment_count': len(fragments),
                'fragment_lengths': fragment_lengths,
                'overhang_info': overhang_info,
                'compatibility_matrix': compatibility_matrix
            }
            
            golden_tests.append(test_case)
            
            # Print summary
            print(f"Sequence length: {len(sequence)} bp")
            print(f"Topology: {'circular' if circular else 'linear'}")
            print(f"Enzymes: {', '.join(enz_names)}")
            print(f"Total cuts: {len(all_cuts)}")
            print(f"Boundary positions: {boundary_positions}")
            print(f"Fragment count: {len(fragments)}")
            print(f"Fragment lengths: {fragment_lengths}")
            print(f"Compatible pairs: {len(compatibility_matrix)}")
    
    return golden_tests


def format_swift_test_code(golden_tests):
    """Format golden test data as Swift XCTest code."""
    
    swift_code = """
// MARK: - Golden Parity Tests (Generated from Python simulator)

"""
    
    for test in golden_tests:
        test_name = test['test_id']
        seq_name = test['sequence_name']
        sequence = test['sequence']
        circular = test['circular']
        enzymes = test['enzyme_set']
        boundary_positions = test['boundary_positions']
        fragment_lengths = test['fragment_lengths']
        compat_matrix = test['compatibility_matrix']
        
        swift_code += f"""    func test_golden_{test_name}() {{
        // Golden test: {seq_name}
        // Enzymes: {', '.join(enzymes)}
        // Expected from Python simulator
        
        let sequence = "{sequence}"
        let circular = {'true' if circular else 'false'}
        
        // Define enzymes (using actual enzyme data from database)
"""
        
        # Add enzyme definitions
        for enz_name in enzymes:
            if enz_name == "EcoRI":
                swift_code += """        let ecoRI = Enzyme(name: "EcoRI", site: "GAATTC", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
"""
            elif enz_name == "BamHI":
                swift_code += """        let bamHI = Enzyme(name: "BamHI", site: "GGATCC", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
"""
            elif enz_name == "PstI":
                swift_code += """        let pstI = Enzyme(name: "PstI", site: "CTGCAG", cutIndexTop: 5, cutIndexBottom: 1, overhangType: .threePrime, notes: nil)
"""
            elif enz_name == "NheI":
                swift_code += """        let nheI = Enzyme(name: "NheI", site: "GCTAGC", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
"""
            elif enz_name == "HindIII":
                swift_code += """        let hindIII = Enzyme(name: "HindIII", site: "AAGCTT", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
"""
            elif enz_name == "XbaI":
                swift_code += """        let xbaI = Enzyme(name: "XbaI", site: "TCTAGA", cutIndexTop: 1, cutIndexBottom: 5, overhangType: .fivePrime, notes: nil)
"""
        
        # Enzyme array
        enz_vars = []
        for enz_name in enzymes:
            if enz_name == "EcoRI":
                enz_vars.append("ecoRI")
            elif enz_name == "BamHI":
                enz_vars.append("bamHI")
            elif enz_name == "PstI":
                enz_vars.append("pstI")
            elif enz_name == "NheI":
                enz_vars.append("nheI")
            elif enz_name == "HindIII":
                enz_vars.append("hindIII")
            elif enz_name == "XbaI":
                enz_vars.append("xbaI")
        
        swift_code += f"""        let enzymes = [{', '.join(enz_vars)}]
        
        // Run digest
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let frags = engine.digest(options: .init(circular: circular, returnSequences: true))
        
        // GOLDEN EXPECTATIONS from Python simulator:
        
        // 1. Boundary positions (sorted cut positions)
        let expectedBoundaries = {boundary_positions}
        
        // Extract actual cut positions from fragments
        var actualCuts = Set<Int>()
        for frag in frags {{
            if frag.start != 0 && !circular {{
                actualCuts.insert(frag.start)
            }}
            if frag.end != sequence.count && !circular {{
                actualCuts.insert(frag.end)
            }}
            if circular {{
                actualCuts.insert(frag.start)
                if frag.end != frag.start {{
                    actualCuts.insert(frag.end)
                }}
            }}
        }}
        let actualBoundaries = Array(actualCuts).sorted()
        
        XCTAssertEqual(actualBoundaries, expectedBoundaries, 
                      "Boundary positions should match Python simulator")
        
        // 2. Fragment lengths
        let expectedLengths = {fragment_lengths}
        let actualLengths = frags.map {{ $0.length }}.sorted()
        
        XCTAssertEqual(actualLengths, expectedLengths, 
                      "Fragment lengths should match Python simulator")
        
        // 3. Fragment count
        XCTAssertEqual(frags.count, {len(fragment_lengths)}, 
                      "Fragment count should match Python simulator")
"""
        
        # Add overhang type assertions
        if test['overhang_info']:
            swift_code += """
        // 4. Overhang types verification
        let sortedFrags = frags.sorted { $0.start < $1.start }
"""
            for oh_info in test['overhang_info']:
                idx = oh_info['fragment_index']
                left = oh_info['left']
                right = oh_info['right']
                
                if left['overhang_types']:
                    left_type = left['overhang_types'][0]
                    swift_type = "blunt" if left_type == "Blunt" else ("fivePrime" if left_type == "5' overhang" else "threePrime")
                    swift_code += f"""        XCTAssertEqual(sortedFrags[{idx}].leftEnd.overhangType, .{swift_type}, 
                      "Fragment {idx} left end should be {left_type}")
"""
                
                if right['overhang_types']:
                    right_type = right['overhang_types'][0]
                    swift_type = "blunt" if right_type == "Blunt" else ("fivePrime" if right_type == "5' overhang" else "threePrime")
                    swift_code += f"""        XCTAssertEqual(sortedFrags[{idx}].rightEnd.overhangType, .{swift_type}, 
                      "Fragment {idx} right end should be {right_type}")
"""
        
        # Add compatibility assertions
        if compat_matrix:
            swift_code += f"""
        // 5. Sticky-end compatibility verification ({len(compat_matrix)} representative pairs)
"""
            for i, comp in enumerate(compat_matrix):
                end_a = comp['end_a']
                end_b = comp['end_b']
                swift_code += f"""        // Compatible pair {i+1}: Fragment {end_a['fragment']} ({end_a['side']}) ↔ Fragment {end_b['fragment']} ({end_b['side']})
        // {end_a['enzyme']} ({end_a['overhang_type']}, {end_a['overhang_seq']}) ↔ {end_b['enzyme']} ({end_b['overhang_type']}, {end_b['overhang_seq']})
        // Directional: {comp['directional']}, Reason: {comp['reason']}
"""
        
        swift_code += """    }

"""
    
    return swift_code


def main():
    """Generate golden tests and output Swift code."""
    print("=" * 80)
    print("GOLDEN TEST GENERATOR - Task 12.1")
    print("=" * 80)
    print()
    print("Generating golden test outputs from Python simulator...")
    print()
    
    # Generate test cases
    golden_tests = generate_test_cases()
    
    # Save as JSON for reference
    output_json = "golden_tests.json"
    with open(output_json, 'w') as f:
        json.dump(golden_tests, f, indent=2)
    print(f"\n✓ Golden test data saved to: {output_json}")
    
    # Generate Swift test code
    swift_code = format_swift_test_code(golden_tests)
    
    # Save Swift code
    output_swift = "golden_tests_swift.txt"
    with open(output_swift, 'w') as f:
        f.write(swift_code)
    print(f"✓ Swift test code saved to: {output_swift}")
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total test cases generated: {len(golden_tests)}")
    print("\nTest cases:")
    for test in golden_tests:
        print(f"  - {test['test_id']}: {test['fragment_count']} fragments, "
              f"{len(test['boundary_positions'])} cuts, "
              f"{len(test['compatibility_matrix'])} compatible pairs")
    
    print("\n✓ Golden tests ready for integration into DigestCoreTests.swift")
    print(f"\nNext steps:")
    print(f"1. Review golden_tests.json for test data")
    print(f"2. Copy tests from {output_swift} to DigestCoreTests.swift")
    print(f"3. Run tests in Xcode to verify parity")


if __name__ == "__main__":
    main()

