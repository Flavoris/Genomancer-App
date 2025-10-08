#!/usr/bin/env python3
"""
Restriction Enzyme Simulator - Phase 2
A simple tool to simulate restriction enzyme cutting on linear DNA sequences.
Now supports multiple enzymes for combined digest analysis.
"""

import argparse
import json
import sys
from typing import List, Dict, Tuple


def load_enzyme_database() -> Dict[str, Dict[str, any]]:
    """
    Load enzyme database from enzymes.json if it exists, otherwise use built-in database.

    Returns:
        Dictionary mapping enzyme names to their properties (sequence, cut_index)
    """
    try:
        # Try to load from enzymes.json file
        with open("enzymes.json", "r") as file:
            enzyme_list = json.load(file)
        
        # Convert list format to dict format for compatibility
        enzymes = {}
        for enzyme in enzyme_list:
            # Validate required fields
            if not all(key in enzyme for key in ["name", "site", "cut_index"]):
                print(f"Warning: Skipping invalid enzyme entry: {enzyme}")
                continue
            
            # Convert JSON format to internal format
            enzymes[enzyme["name"]] = {
                "sequence": enzyme["site"],
                "cut_index": enzyme["cut_index"]
            }
        
        return enzymes
        
    except FileNotFoundError:
        # If enzymes.json doesn't exist, use built-in small database
        print("enzymes.json not found, using built-in enzyme database")
        
    except json.JSONDecodeError as e:
        print(f"Error parsing enzymes.json: {e}")
        print("Using built-in enzyme database instead")
    
    # Built-in enzyme database with recognition sequences and cut positions
    # Cut index indicates where the enzyme cuts (0-based position after the cut)
    # For example: EcoRI recognizes GAATTC and cuts between G^AATTC, so cut_index=1
    enzymes = {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1,  # Cuts between G and A
        },
        "BamHI": {
            "sequence": "GGATCC",
            "cut_index": 1,  # Cuts between G and G
        },
        "HindIII": {
            "sequence": "AAGCTT",
            "cut_index": 1,  # Cuts between A and A
        },
        "PstI": {
            "sequence": "CTGCAG",
            "cut_index": 5,  # Cuts between A and G
        },
        "NotI": {
            "sequence": "GCGGCCGC",
            "cut_index": 2,  # Cuts between C and G
        },
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
    if any(ext in seq_input.lower() for ext in [".fasta", ".txt", ".fa"]):
        try:
            with open(seq_input, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"File '{seq_input}' not found")

        # Parse FASTA format - skip header lines (starting with >)
        sequence_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence_lines.append(line)

        content = "".join(sequence_lines)
    else:
        # Treat as direct DNA sequence
        content = seq_input

    # Remove any whitespace and convert to uppercase
    sequence = "".join(content.split()).upper()

    # Validate that sequence contains only valid DNA bases
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in sequence):
        invalid_chars = set(sequence) - valid_bases
        raise ValueError(f"Invalid DNA characters found: {invalid_chars}")

    return sequence


def find_cut_sites(
    dna_sequence: str, enzyme_sequence: str, cut_index: int
) -> List[int]:
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
        if dna_sequence[i : i + enzyme_len] == enzyme_sequence:
            # Add the cut position relative to the start of the match
            cut_sites.append(i + cut_index)

    return cut_sites


def calculate_fragments(
    dna_sequence: str, cut_positions: List[int]
) -> List[Tuple[str, int]]:
    """
    Calculate fragments for a linear DNA molecule after cutting.

    Args:
        dna_sequence: The original DNA sequence
        cut_positions: List of cut positions (0-based indices)

    Returns:
        List of tuples containing (fragment_sequence, fragment_length)
    """
    if not cut_positions:
        # No cuts found, return the full sequence as one fragment
        return [(dna_sequence, len(dna_sequence))]

    # Add start and end positions for easier calculation
    positions = [0] + sorted(cut_positions) + [len(dna_sequence)]

    # Calculate fragments
    fragments = []
    for i in range(len(positions) - 1):
        start_pos = positions[i]
        end_pos = positions[i + 1]
        fragment_sequence = dna_sequence[start_pos:end_pos]
        fragment_length = end_pos - start_pos
        fragments.append((fragment_sequence, fragment_length))

    return fragments


def find_cut_positions_linear(seq: str, enzyme_name: str, db: Dict[str, Dict[str, any]]) -> List[int]:
    """
    Find all cut positions for a single enzyme in a linear DNA sequence.
    
    Args:
        seq: DNA sequence to search
        enzyme_name: Name of the enzyme
        db: Enzyme database
        
    Returns:
        List of cut positions (0-based indices)
    """
    if enzyme_name not in db:
        return []
    
    enzyme_info = db[enzyme_name]
    enzyme_sequence = enzyme_info["sequence"]
    cut_index = enzyme_info["cut_index"]
    
    return find_cut_sites(seq, enzyme_sequence, cut_index)


def merge_cut_positions(cuts_by_enzyme: Dict[str, List[int]], seq_len: int) -> List[int]:
    """
    Merge cut positions from multiple enzymes into a sorted, unique list.
    
    Args:
        cuts_by_enzyme: Dictionary mapping enzyme names to their cut positions
        seq_len: Length of the DNA sequence
        
    Returns:
        Sorted list of unique cut positions
    """
    all_cuts = set()
    
    for enzyme_cuts in cuts_by_enzyme.values():
        all_cuts.update(enzyme_cuts)
    
    # Filter out any cuts beyond sequence length and sort
    valid_cuts = [cut for cut in all_cuts if 0 <= cut <= seq_len]
    return sorted(valid_cuts)


def fragments_linear(seq_len: int, cuts: List[int]) -> List[int]:
    """
    Calculate fragment lengths for linear DNA after cutting.
    
    Args:
        seq_len: Length of the DNA sequence
        cuts: List of cut positions (0-based indices)
        
    Returns:
        List of fragment lengths
    """
    if not cuts:
        return [seq_len]
    
    # Add start and end positions
    positions = [0] + sorted(cuts) + [seq_len]
    
    # Calculate fragment lengths
    fragment_lengths = []
    for i in range(len(positions) - 1):
        fragment_length = positions[i + 1] - positions[i]
        fragment_lengths.append(fragment_length)
    
    return fragment_lengths


def find_closest_enzyme_names(requested_name: str, available_names: List[str], max_distance: int = 2) -> List[str]:
    """
    Find enzyme names that are similar to the requested name.
    
    Args:
        requested_name: The name that was requested
        available_names: List of available enzyme names
        max_distance: Maximum edit distance for similarity
        
    Returns:
        List of similar enzyme names
    """
    def edit_distance(s1: str, s2: str) -> int:
        """Calculate edit distance between two strings."""
        m, n = len(s1), len(s2)
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        
        for i in range(m + 1):
            dp[i][0] = i
        for j in range(n + 1):
            dp[0][j] = j
            
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if s1[i-1] == s2[j-1]:
                    dp[i][j] = dp[i-1][j-1]
                else:
                    dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
        
        return dp[m][n]
    
    similar_names = []
    requested_lower = requested_name.lower()
    
    for name in available_names:
        name_lower = name.lower()
        distance = edit_distance(requested_lower, name_lower)
        
        # Check exact match (case-insensitive)
        if name_lower == requested_lower:
            return [name]
        
        # Check if it's a substring match
        if requested_lower in name_lower or name_lower in requested_lower:
            similar_names.append(name)
        # Check edit distance
        elif distance <= max_distance:
            similar_names.append(name)
    
    return similar_names[:5]  # Return top 5 matches


def main():
    """Main function to run the restriction enzyme simulator."""
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Restriction Enzyme Simulator - Phase 2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
  python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI BamHI
  python sim.py --seq my_dna.fasta --enz BamHI HindIII
        """,
    )

    parser.add_argument(
        "--seq", required=True, help="DNA sequence (string) or path to FASTA/text file"
    )
    parser.add_argument(
        "--enz", required=True, nargs='+', help="One or more enzyme names (see available enzymes in output)"
    )

    args = parser.parse_args()

    try:
        # Handle empty enzyme list
        if not args.enz:
            print("Error: No enzymes specified.")
            print("Usage: python sim.py --seq <sequence> --enz <enzyme1> [enzyme2] ...")
            sys.exit(2)

        # Load enzyme database
        ENZYMES = load_enzyme_database()
        
        # Create case-insensitive lookup dictionary
        enzyme_lookup = {name.lower(): name for name in ENZYMES.keys()}
        available_names = list(ENZYMES.keys())

        # Validate and normalize enzyme names
        validated_enzymes = []
        for enzyme_name in args.enz:
            enzyme_lower = enzyme_name.lower()
            if enzyme_lower in enzyme_lookup:
                validated_enzymes.append(enzyme_lookup[enzyme_lower])
            else:
                # Find closest matches
                closest = find_closest_enzyme_names(enzyme_name, available_names)
                print(f"Error: Enzyme '{enzyme_name}' not found in database.")
                if closest:
                    print(f"Did you mean one of these? {', '.join(closest)}")
                else:
                    print("No similar enzyme names found.")
                print(f"Available enzymes: {', '.join(sorted(available_names))}")
                print(f"Total enzymes available: {len(available_names)}")
                sys.exit(2)

        # Read DNA sequence
        print(f"Reading DNA sequence from: {args.seq}")
        dna_sequence = read_dna_sequence(args.seq)
        
        if not dna_sequence:
            print("Error: Empty sequence after filtering.")
            sys.exit(1)
            
        print(f"DNA sequence length: {len(dna_sequence)} bp")
        print(f"DNA sequence: {dna_sequence}")
        print()

        # Process each enzyme individually
        cuts_by_enzyme = {}
        
        for enzyme_name in validated_enzymes:
            enzyme_info = ENZYMES[enzyme_name]
            recognition_seq = enzyme_info["sequence"]
            cut_index = enzyme_info["cut_index"]
            
            print(f"Enzyme: {enzyme_name}")
            print(f"Recognition sequence: {recognition_seq}")
            print(f"Cut position: after position {cut_index} in recognition sequence")
            
            # Find cut positions for this enzyme
            cut_positions = find_cut_positions_linear(dna_sequence, enzyme_name, ENZYMES)
            cuts_by_enzyme[enzyme_name] = cut_positions
            
            if cut_positions:
                print(f"Found cut positions: {cut_positions}")
            else:
                print("No cut sites found.")
            print()

        # Calculate combined digest
        print("COMBINED DIGEST SUMMARY")
        print("=" * 40)
        
        # Merge all cut positions
        all_cuts = merge_cut_positions(cuts_by_enzyme, len(dna_sequence))
        
        print(f"Total cuts (unique across enzymes): {len(all_cuts)}")
        if all_cuts:
            print(f"Cut positions: {all_cuts}")
        
        # Calculate fragment lengths for combined digest
        fragment_lengths = fragments_linear(len(dna_sequence), all_cuts)
        print(f"Fragment lengths: {fragment_lengths}")
        
        # Verify total length
        total_length = sum(fragment_lengths)
        print(f"Total length verification: {total_length} bp (expected: {len(dna_sequence)} bp)")
        
        if total_length != len(dna_sequence):
            print("Warning: Fragment lengths don't sum to original sequence length!")

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)


if __name__ == "__main__":
    main()
