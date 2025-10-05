#!/usr/bin/env python3
"""
Restriction Enzyme Simulator - Phase 1
A simple tool to simulate restriction enzyme cutting on linear DNA sequences.
"""

import argparse
import sys
from typing import List, Tuple, Dict


def load_enzyme_database() -> Dict[str, Dict[str, any]]:
    """
    Load the hardcoded enzyme database with recognition sequences and cut indices.
    
    Returns:
        Dictionary mapping enzyme names to their properties
    """
    # Hardcoded enzyme database with recognition sequences and cut positions
    # Cut index indicates where the enzyme cuts (0-based position after the cut)
    # For example: EcoRI recognizes GAATTC and cuts between G^AATTC, so cut_index=1
    enzymes = {
        'EcoRI': {
            'sequence': 'GAATTC',
            'cut_index': 1  # Cuts between G and A
        },
        'BamHI': {
            'sequence': 'GGATCC',
            'cut_index': 1  # Cuts between G and G
        },
        'HindIII': {
            'sequence': 'AAGCTT',
            'cut_index': 1  # Cuts between A and A
        },
        'PstI': {
            'sequence': 'CTGCAG',
            'cut_index': 5  # Cuts between A and G
        },
        'NotI': {
            'sequence': 'GCGGCCGC',
            'cut_index': 2  # Cuts between C and G
        }
    }
    return enzymes


def read_dna_sequence(seq_input: str) -> str:
    """
    Read DNA sequence from either a string or a FASTA file.
    
    Args:
        seq_input: Either a DNA sequence string or a file path
        
    Returns:
        DNA sequence as uppercase string
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the sequence contains invalid characters
    """
    # Check if input looks like a file path (contains .fasta, .txt, or .fa)
    if any(ext in seq_input.lower() for ext in ['.fasta', '.txt', '.fa']):
        try:
            with open(seq_input, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"File '{seq_input}' not found")
        
        # Parse FASTA format - skip header lines (starting with >)
        sequence_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('>'):
                sequence_lines.append(line)
        
        content = ''.join(sequence_lines)
    else:
        # Treat as direct DNA sequence
        content = seq_input
    
    # Remove any whitespace and convert to uppercase
    sequence = ''.join(content.split()).upper()
    
    # Validate that sequence contains only valid DNA bases
    valid_bases = set('ATCG')
    if not all(base in valid_bases for base in sequence):
        invalid_chars = set(sequence) - valid_bases
        raise ValueError(f"Invalid DNA characters found: {invalid_chars}")
    
    return sequence


def find_cut_sites(dna_sequence: str, enzyme_sequence: str, cut_index: int) -> List[int]:
    """
    Find all cut sites for a given enzyme in the DNA sequence.
    
    Args:
        dna_sequence: The DNA sequence to search
        enzyme_sequence: The recognition sequence of the enzyme
        cut_index: 0-based index within the recognition site where the cut occurs
        
    Returns:
        List of positions where cuts occur (0-based indices)
    """
    cut_sites = []
    enzyme_len = len(enzyme_sequence)
    
    # Search for the enzyme recognition sequence (case-insensitive)
    for i in range(len(dna_sequence) - enzyme_len + 1):
        if dna_sequence[i:i + enzyme_len] == enzyme_sequence:
            # Add the cut position relative to the start of the match
            cut_sites.append(i + cut_index)
    
    return cut_sites


def calculate_fragment_lengths(dna_sequence: str, cut_positions: List[int]) -> List[int]:
    """
    Calculate fragment lengths for a linear DNA molecule after cutting.
    
    Args:
        dna_sequence: The original DNA sequence
        cut_positions: List of cut positions (0-based indices)
        
    Returns:
        List of fragment lengths
    """
    if not cut_positions:
        # No cuts found, return the full sequence as one fragment
        return [len(dna_sequence)]
    
    # Add start and end positions for easier calculation
    positions = [0] + sorted(cut_positions) + [len(dna_sequence)]
    
    # Calculate fragment lengths
    fragment_lengths = []
    for i in range(len(positions) - 1):
        fragment_length = positions[i + 1] - positions[i]
        fragment_lengths.append(fragment_length)
    
    return fragment_lengths


def main():
    """Main function to run the restriction enzyme simulator."""
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description='Restriction Enzyme Simulator - Phase 1',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
  python sim.py --seq my_dna.fasta --enz BamHI
        """
    )
    
    parser.add_argument('--seq', required=True,
                       help='DNA sequence (string) or path to FASTA/text file')
    parser.add_argument('--enz', required=True,
                       help='Enzyme name (EcoRI, BamHI, HindIII, PstI, NotI)')
    
    args = parser.parse_args()
    
    try:
        # Load enzyme database
        enzymes = load_enzyme_database()
        
        # Check if the requested enzyme exists
        if args.enz not in enzymes:
            print(f"Error: Enzyme '{args.enz}' not found in database.")
            print(f"Available enzymes: {', '.join(enzymes.keys())}")
            sys.exit(1)
        
        # Read DNA sequence
        print(f"Reading DNA sequence from: {args.seq}")
        dna_sequence = read_dna_sequence(args.seq)
        print(f"DNA sequence length: {len(dna_sequence)} bp")
        print(f"DNA sequence: {dna_sequence}")
        print()
        
        # Get enzyme information
        enzyme_info = enzymes[args.enz]
        recognition_seq = enzyme_info['sequence']
        cut_index = enzyme_info['cut_index']
        
        print(f"Enzyme: {args.enz}")
        print(f"Recognition sequence: {recognition_seq}")
        print(f"Cut position: after position {cut_index} in recognition sequence")
        print()
        
        # Find cut sites
        cut_sites = find_cut_sites(dna_sequence, recognition_seq, cut_index)
        
        if not cut_sites:
            print("No cut sites found!")
            print(f"Result: 1 fragment of {len(dna_sequence)} bp")
        else:
            print(f"Found {len(cut_sites)} cut site(s):")
            for i, site in enumerate(cut_sites, 1):
                print(f"  Cut site {i}: position {site}")
            print()
            
            # Calculate fragment lengths
            fragment_lengths = calculate_fragment_lengths(dna_sequence, cut_sites)
            
            print("Fragment analysis:")
            print(f"  Number of fragments: {len(fragment_lengths)}")
            for i, length in enumerate(fragment_lengths, 1):
                print(f"  Fragment {i}: {length} bp")
            
            print(f"\nTotal length: {sum(fragment_lengths)} bp")
    
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)


if __name__ == "__main__":
    main()
