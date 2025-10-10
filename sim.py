#!/usr/bin/env python3
"""
Restriction Enzyme Simulator - Phase 3
A tool to simulate restriction enzyme cutting on linear and circular DNA sequences.
Supports multiple enzymes for combined digest analysis and circular topology.
"""

import argparse
import json
import re
import sys
import unicodedata
from typing import List, Dict, Tuple
from fragment_calculator import compute_fragments, validate_fragment_total

# IUPAC degenerate base mapping
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[GC]", 
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "D": "[AGT]", 
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
}


def iupac_to_regex(site: str) -> str:
    """
    Convert IUPAC degenerate base notation to regex character classes.
    
    Args:
        site: Recognition site string that may contain IUPAC letters
        
    Returns:
        Regex pattern (no anchors) where IUPAC letters are expanded to character classes
        
    Raises:
        ValueError: If site contains invalid characters
    """
    allowed_chars = set("ACGT") | set("RYWSMKNBDHV")
    site_upper = site.upper()
    
    for char in site_upper:
        if char not in allowed_chars:
            raise ValueError(f"Invalid character '{char}' in recognition site '{site}'. "
                           f"Allowed characters: A,C,G,T,R,Y,W,S,M,K,B,D,H,V,N")
    
    return "".join(IUPAC[ch] for ch in site_upper)


def normalize(name: str) -> str:
    """
    Normalize enzyme name by removing diacritics, converting to lowercase, 
    and removing whitespace/hyphens.
    
    Args:
        name: Original enzyme name
        
    Returns:
        Normalized name for lookup purposes
    """
    # Remove diacritics (e.g., HF® -> HFR)
    normalized = unicodedata.normalize('NFKD', name)
    normalized = ''.join(c for c in normalized if not unicodedata.combining(c))
    
    # Convert to lowercase and remove whitespace/hyphens
    normalized = normalized.lower().replace(' ', '').replace('-', '')
    
    return normalized


def load_enzyme_database() -> Dict[str, Dict[str, any]]:
    """
    Load enzyme database from enzymes.json if it exists, otherwise use built-in database.
    Supports the new format with overhang_type.

    Returns:
        Dictionary mapping enzyme names to their properties (sequence, cut_index, overhang_type)
    """
    try:
        # Try to load from enzymes.json file
        with open("enzymes.json", "r") as file:
            enzyme_list = json.load(file)
        
        # Convert list format to dict format for compatibility
        enzymes = {}
        normalized_names = {}  # Track normalized names for duplicate detection
        
        for enzyme in enzyme_list:
            # Validate required fields
            if not all(key in enzyme for key in ["name", "site"]):
                print(f"Warning: Skipping invalid enzyme entry: {enzyme}")
                continue
            
            # Check for cut specification
            if "cut_index" not in enzyme:
                print(f"Warning: Skipping enzyme '{enzyme['name']}': missing cut_index")
                continue
            
            # Normalize and validate site
            site = enzyme["site"].upper()
            try:
                iupac_to_regex(site)
            except ValueError as e:
                print(f"Warning: Skipping enzyme '{enzyme['name']}': {e}")
                continue
            
            # Validate cut specification
            site_len = len(site)
            cut_index = enzyme["cut_index"]
            if not (0 <= cut_index <= site_len):
                print(f"Warning: Skipping enzyme '{enzyme['name']}': "
                      f"cut_index {cut_index} must be between 0 and {site_len}")
                continue
            
            # Get overhang_type, default to "Unknown" if missing (legacy support)
            overhang_type = enzyme.get("overhang_type", "Unknown")
            if overhang_type not in ["5' overhang", "3' overhang", "Blunt", "Unknown"]:
                print(f"Warning: Invalid overhang_type '{overhang_type}' for enzyme '{enzyme['name']}', setting to 'Unknown'")
                overhang_type = "Unknown"
            
            # Handle duplicate names by adding suffixes
            original_name = enzyme["name"]
            normalized_name = normalize(original_name)
            
            if normalized_name in normalized_names:
                # This is a duplicate - find the next available suffix
                counter = 2
                while f"{original_name}#{counter}" in enzymes:
                    counter += 1
                final_name = f"{original_name}#{counter}"
                normalized_names[normalized_name].append(final_name)
            else:
                final_name = original_name
                normalized_names[normalized_name] = [final_name]
            
            # Store enzyme with new schema
            enzyme_data = {
                "sequence": site,
                "cut_index": cut_index,
                "overhang_type": overhang_type
            }
            
            enzymes[final_name] = enzyme_data
        
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
            "overhang_type": "5' overhang"
        },
        "BamHI": {
            "sequence": "GGATCC",
            "cut_index": 1,  # Cuts between G and G
            "overhang_type": "5' overhang"
        },
        "HindIII": {
            "sequence": "AAGCTT",
            "cut_index": 1,  # Cuts between A and A
            "overhang_type": "5' overhang"
        },
        "PstI": {
            "sequence": "CTGCAG",
            "cut_index": 5,  # Cuts between A and G
            "overhang_type": "3' overhang"
        },
        "NotI": {
            "sequence": "GCGGCCGC",
            "cut_index": 2,  # Cuts between C and G
            "overhang_type": "5' overhang"
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

    # Validate that sequence contains only valid DNA bases or IUPAC characters
    valid_bases = set("ATCG") | set("RYWSMKNBDHV")
    if not all(base in valid_bases for base in sequence):
        invalid_chars = set(sequence) - valid_bases
        raise ValueError(f"Invalid DNA characters found: {invalid_chars}")

    return sequence


def find_cut_sites(
    dna_sequence: str, enzyme_sequence: str, cut_index: int
) -> List[int]:
    """
    Find all cut sites for a given enzyme in the DNA sequence using IUPAC pattern matching.

    Args:
        dna_sequence: The DNA sequence to search
        enzyme_sequence: The recognition sequence of the enzyme (may contain IUPAC letters)
        cut_index: 0-based index within the recognition site where the enzyme cuts

    Returns:
        List of break positions (0-based indices)
    """
    break_positions = []
    
    # Convert IUPAC site to regex pattern and compile with lookahead for overlapping matches
    regex_pattern = f"(?={iupac_to_regex(enzyme_sequence)})"
    pattern = re.compile(regex_pattern, flags=re.IGNORECASE)
    
    # Find all overlapping matches
    for match in pattern.finditer(dna_sequence):
        # Calculate break position using cut_index
        break_pos = match.start() + cut_index
        break_positions.append(break_pos)

    return break_positions


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
        List of break positions (0-based indices)
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
        description="Restriction Enzyme Simulator - Phase 3 (Linear and Circular DNA)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Linear digestion (default)
  python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
  python sim.py --seq my_dna.fasta --enz EcoRI BamHI
  
  # Circular digestion
  python sim.py --seq plasmid.fasta --enz EcoRI BamHI --circular
  python sim.py --seq plasmid.fasta --enz EcoRI --circular --circular_single_cut_linearizes
        """,
    )

    parser.add_argument(
        "--seq", required=True, help="DNA sequence (string) or path to FASTA/text file"
    )
    parser.add_argument(
        "--enz", required=True, nargs='+', help="One or more enzyme names (see available enzymes in output)"
    )
    parser.add_argument(
        "--circular", action="store_true", default=False,
        help="Treat DNA as circular (default: linear)"
    )
    parser.add_argument(
        "--circular_single_cut_linearizes", action="store_true", default=False,
        help="In circular mode, one cut yields two fragments instead of one intact circle (default: False)"
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
        
        # Create normalized lookup dictionary for case-insensitive matching
        normalized_lookup = {}
        for name in ENZYMES.keys():
            norm_name = normalize(name)
            if norm_name not in normalized_lookup:
                normalized_lookup[norm_name] = []
            normalized_lookup[norm_name].append(name)
        
        available_names = list(ENZYMES.keys())

        # Validate and normalize enzyme names
        validated_enzymes = []
        for enzyme_name in args.enz:
            norm_name = normalize(enzyme_name)
            if norm_name in normalized_lookup:
                variants = normalized_lookup[norm_name]
                if len(variants) == 1:
                    validated_enzymes.append(variants[0])
                else:
                    # Ambiguous base name - show variants and exit
                    print(f"Error: Multiple enzyme variants found for '{enzyme_name}':")
                    for variant in variants:
                        print(f"  - {variant}")
                    print("Please specify the exact enzyme name with suffix if needed.")
                    sys.exit(2)
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
        
        # Display topology mode
        topology = "circular" if args.circular else "linear"
        print(f"Topology: {topology}")
        if args.circular and args.circular_single_cut_linearizes:
            print("Single-cut behavior: linearizes (yields 2 fragments)")
        elif args.circular:
            print("Single-cut behavior: intact circle (yields 1 fragment)")
            
        print(f"DNA sequence length: {len(dna_sequence)} bp")
        print(f"DNA sequence: {dna_sequence}")
        print()

        # Process each enzyme individually and collect cut metadata
        cuts_by_enzyme = {}
        cut_metadata = {}  # Maps cut position -> list of enzyme metadata
        
        for enzyme_name in validated_enzymes:
            enzyme_info = ENZYMES[enzyme_name]
            recognition_seq = enzyme_info["sequence"]
            cut_index = enzyme_info["cut_index"]
            overhang_type = enzyme_info["overhang_type"]
            
            print(f"Enzyme: {enzyme_name}")
            print(f"Site:   {recognition_seq}")
            print(f"Cut @:  index {cut_index}")
            print(f"Overhang: {overhang_type}")
            
            # Find cut positions for this enzyme
            break_positions = find_cut_positions_linear(dna_sequence, enzyme_name, ENZYMES)
            cuts_by_enzyme[enzyme_name] = break_positions
            
            # Store metadata for each cut position
            for pos in break_positions:
                if pos not in cut_metadata:
                    cut_metadata[pos] = []
                cut_metadata[pos].append({
                    'enzyme': enzyme_name,
                    'site': recognition_seq,
                    'cut_index': cut_index,
                    'overhang_type': overhang_type
                })
            
            if break_positions:
                print(f"Matches at positions: {', '.join(map(str, break_positions))}")
            else:
                print("No cut sites found.")
            print()

        # Merge all cut positions
        all_cuts = merge_cut_positions(cuts_by_enzyme, len(dna_sequence))
        
        # Compute fragments using new calculator
        fragments = compute_fragments(
            cut_positions=all_cuts,
            seq_len=len(dna_sequence),
            circular=args.circular,
            circular_single_cut_linearizes=args.circular_single_cut_linearizes,
            cut_metadata=cut_metadata
        )
        
        # Display detailed fragment information
        print("=" * 80)
        print("DIGESTION RESULTS")
        print("=" * 80)
        print(f"Mode: {topology}")
        print(f"Sequence length: {len(dna_sequence)} bp")
        print(f"Total cuts: {len(all_cuts)}")
        if all_cuts:
            print(f"Cut positions: {', '.join(map(str, sorted(all_cuts)))}")
        print(f"Fragments generated: {len(fragments)}")
        print()
        
        # Display fragments table
        print("Fragment Details:")
        print("-" * 80)
        header = f"{'Index':<8} {'Start':<8} {'End':<8} {'Length':<10} {'Wraps':<8} {'Boundaries'}"
        print(header)
        print("-" * 80)
        
        for frag in fragments:
            idx = frag['index']
            start = frag['start']
            end = frag['end']
            length = frag['length']
            wraps = 'Yes' if frag['wraps'] else 'No'
            
            # Build boundary info string
            left_info = "START"
            right_info = "END"
            
            if frag['boundaries']['left_cut'] is not None:
                left_pos = frag['boundaries']['left_cut']['pos']
                left_enzymes = [e['enzyme'] for e in frag['boundaries']['left_cut']['enzymes']]
                left_info = f"{left_pos}({','.join(left_enzymes) if left_enzymes else '?'})"
            
            if frag['boundaries']['right_cut'] is not None:
                right_pos = frag['boundaries']['right_cut']['pos']
                right_enzymes = [e['enzyme'] for e in frag['boundaries']['right_cut']['enzymes']]
                right_info = f"{right_pos}({','.join(right_enzymes) if right_enzymes else '?'})"
            
            boundary_str = f"{left_info} -> {right_info}"
            
            print(f"{idx:<8} {start:<8} {end:<8} {length:<10} {wraps:<8} {boundary_str}")
        
        print("-" * 80)
        
        # Verify total length
        if not validate_fragment_total(fragments, len(dna_sequence)):
            total = sum(f['length'] for f in fragments)
            print(f"WARNING: Fragment lengths don't sum to sequence length! ({total} vs {len(dna_sequence)})")
        else:
            print(f"✓ Fragment lengths sum correctly to {len(dna_sequence)} bp")
        
        # Display cut site details
        if all_cuts:
            print()
            print("Cut Site Details:")
            print("-" * 80)
            for pos in sorted(all_cuts):
                enzymes_at_pos = cut_metadata.get(pos, [])
                for enz_meta in enzymes_at_pos:
                    print(f"  Position {pos}: {enz_meta['enzyme']} "
                          f"(site: {enz_meta['site']}, overhang: {enz_meta['overhang_type']})")
        
        print()

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)


if __name__ == "__main__":
    main()
