#!/usr/bin/env python3
"""
Fragment Sequence Extraction Module

This module provides reference implementation for extracting DNA sequences
from restriction enzyme digest fragments. It demonstrates how to:
1. Extract the actual sequence for each fragment
2. Handle both linear and circular DNA topologies
3. Manage sequence wrapping for circular plasmids
4. Provide overhang and end information

This serves as a reference for the iOS/Swift implementation.
"""

from typing import List, Dict, Optional, Tuple
from fragment_calculator import (
    compute_fragments_with_sequences,
    Fragment,
    EndInfo,
    slice_circular
)


def extract_fragment_sequence(
    dna_sequence: str,
    start_position: int,
    end_position: int,
    circular: bool = False
) -> str:
    """
    Extract a DNA sequence fragment between two positions.
    
    This is the core function for getting the actual DNA sequence
    of a fragment after restriction digest.
    
    Args:
        dna_sequence: The complete DNA sequence
        start_position: Start position (0-based, inclusive)
        end_position: End position (0-based, exclusive)
        circular: Whether the DNA is circular (handles wrap-around)
        
    Returns:
        The DNA sequence fragment in 5'→3' orientation
        
    Examples:
        >>> dna = "ATGCGAATTCGCTAGC"
        >>> extract_fragment_sequence(dna, 0, 5, False)
        'ATGCG'
        
        >>> # Circular DNA with wrap-around
        >>> dna = "GAATTCGGGCCC"  # 12 bp
        >>> extract_fragment_sequence(dna, 10, 4, True)  # wraps around
        'CCGAAT'
    """
    seq_len = len(dna_sequence)
    
    if seq_len == 0:
        return ""
    
    if not circular:
        # Linear mode: simple slice
        return dna_sequence[start_position:end_position]
    else:
        # Circular mode: handle wrap-around
        return slice_circular(dna_sequence, start_position, end_position)


def get_fragment_with_sequence(
    dna_sequence: str,
    enzymes: List[Dict],
    enzyme_db: Dict[str, Dict],
    circular: bool = False
) -> List[Dict]:
    """
    Perform restriction digest and return fragments with full sequence information.
    
    This function demonstrates the complete workflow for getting fragment sequences:
    1. Find all cut sites for the given enzymes
    2. Compute fragments with positions
    3. Extract actual DNA sequences for each fragment
    4. Include overhang and end information
    
    Args:
        dna_sequence: The complete DNA sequence to digest
        enzymes: List of enzyme names to use
        enzyme_db: Database of enzyme information (from enzymes.json)
        circular: Whether DNA is circular
        
    Returns:
        List of dictionaries, each containing:
            - index: Fragment number (0-based)
            - start: Start position
            - end: End position  
            - length: Fragment length in bp
            - sequence: Actual DNA sequence (5'→3')
            - left_end: Information about left end (overhang, enzyme, etc.)
            - right_end: Information about right end
            - wraps: Whether fragment wraps around (circular only)
            
    Example:
        >>> # Example with EcoRI cutting at position 6 in a linear sequence
        >>> dna = "ATGCGAATTCGCTAGC"  # EcoRI site at position 5 (G^AATTC)
        >>> enzymes = ["EcoRI"]
        >>> enzyme_db = {"EcoRI": {"site": "GAATTC", "cut_index": 1, ...}}
        >>> fragments = get_fragment_with_sequence(dna, enzymes, enzyme_db, False)
        >>> # Returns two fragments:
        >>> # Fragment 0: "ATGCG" (0-6)
        >>> # Fragment 1: "AATTCGCTAGC" (6-16)
    """
    from sim import find_cut_sites
    
    # Step 1: Find all cut sites
    cut_sites = []
    cut_metadata = {}
    
    for enzyme_name in enzymes:
        if enzyme_name not in enzyme_db:
            continue
            
        enzyme_info = enzyme_db[enzyme_name]
        site = enzyme_info.get('site', '')
        cut_index = enzyme_info.get('cut_index', 0)
        overhang_type = enzyme_info.get('overhang_type', 'Blunt')
        
        # Find all occurrences of this enzyme's recognition site
        positions = find_cut_sites(dna_sequence, site, cut_index, circular)
        
        for pos in positions:
            cut_sites.append(pos)
            if pos not in cut_metadata:
                cut_metadata[pos] = []
            cut_metadata[pos].append({
                'enzyme': enzyme_name,
                'site': site,
                'cut_index': cut_index,
                'overhang_type': overhang_type
            })
    
    # Step 2: Compute fragments with sequences
    fragments = compute_fragments_with_sequences(
        dna_sequence=dna_sequence,
        cut_positions=cut_sites,
        circular=circular,
        circular_single_cut_linearizes=True,
        cut_metadata=cut_metadata
    )
    
    # Step 3: Convert to dictionary format for easier use
    result = []
    for idx, frag in enumerate(fragments):
        left_end_info = None
        right_end_info = None
        
        if frag.enzymes_at_ends[0]:
            left_end = frag.enzymes_at_ends[0]
            left_end_info = {
                'enzyme': left_end.enzyme,
                'overhang_type': left_end.overhang_type,
                'overhang_len': left_end.overhang_len,
                'end_bases': left_end.end_bases,
                'recognition_site': left_end.recognition_site
            }
        
        if frag.enzymes_at_ends[1]:
            right_end = frag.enzymes_at_ends[1]
            right_end_info = {
                'enzyme': right_end.enzyme,
                'overhang_type': right_end.overhang_type,
                'overhang_len': right_end.overhang_len,
                'end_bases': right_end.end_bases,
                'recognition_site': right_end.recognition_site
            }
        
        result.append({
            'index': idx,
            'start': frag.start_idx,
            'end': frag.end_idx,
            'length': frag.length,
            'sequence': frag.sequence,
            'left_end': left_end_info,
            'right_end': right_end_info,
            'wraps': frag.wraps
        })
    
    return result


def format_sequence_display(
    sequence: str,
    line_length: int = 60,
    group_size: int = 10,
    show_positions: bool = True
) -> str:
    """
    Format a DNA sequence for readable display.
    
    Args:
        sequence: DNA sequence to format
        line_length: Number of bases per line
        group_size: Number of bases per group (space-separated)
        show_positions: Whether to show position markers
        
    Returns:
        Formatted sequence string
        
    Example:
        >>> seq = "ATGCGAATTCGCTAGCTAGCTAGCTAGCTAGCTAGC"
        >>> print(format_sequence_display(seq, line_length=30, group_size=10))
           1 ATGCGAATTC GCTAGCTAG CTAGCTAG
          31 CTAGCTAGC
    """
    lines = []
    pos = 0
    
    while pos < len(sequence):
        # Get line chunk
        chunk = sequence[pos:pos + line_length]
        
        # Split into groups
        groups = []
        for i in range(0, len(chunk), group_size):
            groups.append(chunk[i:i + group_size])
        
        # Format line with position
        if show_positions:
            line = f"{pos + 1:>4} {' '.join(groups)}"
        else:
            line = ' '.join(groups)
        
        lines.append(line)
        pos += line_length
    
    return '\n'.join(lines)


def print_fragment_details(fragment: Dict, show_full_sequence: bool = False):
    """
    Print detailed information about a fragment including its sequence.
    
    Args:
        fragment: Fragment dictionary from get_fragment_with_sequence()
        show_full_sequence: If True, show formatted full sequence; 
                           if False, show truncated sequence
    """
    print(f"\n{'='*70}")
    print(f"Fragment {fragment['index'] + 1}")
    print(f"{'='*70}")
    print(f"Position: [{fragment['start']} .. {fragment['end']})")
    print(f"Length: {fragment['length']} bp")
    
    if fragment.get('wraps', False):
        print(f"⚠ Wraps around origin (circular DNA)")
    
    # Left end
    if fragment['left_end']:
        left = fragment['left_end']
        print(f"\nLeft End (5'):")
        print(f"  Enzyme: {left['enzyme']}")
        print(f"  Recognition site: {left['recognition_site']}")
        print(f"  Overhang: {left['overhang_type']}")
        if left['end_bases']:
            print(f"  Sticky sequence: {left['end_bases']}")
    else:
        print(f"\nLeft End (5'): Natural terminus (blunt)")
    
    # Right end
    if fragment['right_end']:
        right = fragment['right_end']
        print(f"\nRight End (3'):")
        print(f"  Enzyme: {right['enzyme']}")
        print(f"  Recognition site: {right['recognition_site']}")
        print(f"  Overhang: {right['overhang_type']}")
        if right['end_bases']:
            print(f"  Sticky sequence: {right['end_bases']}")
    else:
        print(f"\nRight End (3'): Natural terminus (blunt)")
    
    # Sequence
    seq = fragment['sequence']
    print(f"\nSequence (5'→3'):")
    if show_full_sequence or len(seq) <= 120:
        # Show full formatted sequence
        print(format_sequence_display(seq))
    else:
        # Show truncated sequence
        print(f"  {seq[:40]} ... {seq[-40:]}")
        print(f"  (showing first and last 40 bp, use show_full_sequence=True for complete sequence)")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

def main():
    """
    Demonstrate fragment sequence extraction with example data.
    """
    import json
    import os
    
    print("Fragment Sequence Extraction - Reference Implementation")
    print("="*70)
    
    # Load enzyme database
    enzyme_db_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'enzymes.json')
    try:
        with open(enzyme_db_path, 'r') as f:
            enzyme_list = json.load(f)
            # Convert to dict for easier lookup
            enzyme_db = {e['name']: e for e in enzyme_list}
    except FileNotFoundError:
        print(f"Error: Could not find enzymes.json at {enzyme_db_path}")
        return
    
    # Example 1: Linear DNA with EcoRI
    print("\n" + "="*70)
    print("EXAMPLE 1: Linear DNA with EcoRI")
    print("="*70)
    
    dna_linear = "ATGCGAATTCGCTAGCTAGCTAGCTGAATTCGGG"
    enzymes = ["EcoRI"]
    
    print(f"\nDNA Sequence ({len(dna_linear)} bp):")
    print(format_sequence_display(dna_linear, line_length=40))
    print(f"\nEnzymes: {', '.join(enzymes)}")
    
    fragments = get_fragment_with_sequence(dna_linear, enzymes, enzyme_db, circular=False)
    print(f"\nResult: {len(fragments)} fragments")
    
    for frag in fragments:
        print_fragment_details(frag, show_full_sequence=True)
    
    # Example 2: Circular DNA with multiple enzymes
    print("\n\n" + "="*70)
    print("EXAMPLE 2: Circular Plasmid with EcoRI and BamHI")
    print("="*70)
    
    dna_circular = "ATGCGAATTCGCTAGCGGATCCTAGCTAGCTAGCTGAATTCGGG"
    enzymes = ["EcoRI", "BamHI"]
    
    print(f"\nPlasmid Sequence ({len(dna_circular)} bp):")
    print(format_sequence_display(dna_circular, line_length=40))
    print(f"\nEnzymes: {', '.join(enzymes)}")
    
    fragments = get_fragment_with_sequence(dna_circular, enzymes, enzyme_db, circular=True)
    print(f"\nResult: {len(fragments)} fragments")
    
    for frag in fragments:
        print_fragment_details(frag, show_full_sequence=True)
    
    # Example 3: Using actual test sequence
    test_seq_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_sequence.txt')
    if os.path.exists(test_seq_path):
        print("\n\n" + "="*70)
        print("EXAMPLE 3: Test Sequence from data/test_sequence.txt")
        print("="*70)
        
        with open(test_seq_path, 'r') as f:
            test_seq = f.read().strip().replace('\n', '').replace(' ', '').upper()
        
        enzymes = ["HaeIII", "HinfI"]
        
        print(f"\nSequence length: {len(test_seq)} bp")
        print(f"Enzymes: {', '.join(enzymes)}")
        
        fragments = get_fragment_with_sequence(test_seq, enzymes, enzyme_db, circular=False)
        print(f"\nResult: {len(fragments)} fragments")
        
        # Show summary
        print(f"\nFragment Summary:")
        for frag in fragments:
            seq = frag['sequence']
            if len(seq) > 40:
                seq_display = f"{seq[:20]}...{seq[-20:]}"
            else:
                seq_display = seq
            print(f"  Fragment {frag['index'] + 1}: {frag['length']} bp - {seq_display}")


if __name__ == "__main__":
    main()

