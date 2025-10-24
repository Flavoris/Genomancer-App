#!/usr/bin/env python3
"""
Restriction Enzyme Simulator - Phase 4
A tool to simulate restriction enzyme cutting on linear and circular DNA sequences.
Supports multiple enzymes for combined digest analysis and circular topology.
"""

import argparse
import json
import re
import sys
import unicodedata
from typing import List, Dict, Tuple
from fragment_calculator import (
    compute_fragments, validate_fragment_total, build_restriction_map, simulate_gel,
    compute_fragments_with_sequences, elide_sequence, extract_fragment_ends_for_ligation,
    compute_end_metadata
)
from gel_ladders import get_ladder
from graphics import render_plasmid_map, render_linear_map, render_fragment_diagram, svg_to_png
from ligation_compatibility import (
    calculate_compatibility, format_pairs_output, format_matrix_output,
    format_detailed_output, export_to_json,
    theoretical_end_from_enzyme, calculate_theoretical_compatibility,
    format_theoretical_pairs, format_theoretical_matrix, format_theoretical_detailed
)
from exporters import export_genbank, export_csv

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
    dna_sequence: str, enzyme_sequence: str, cut_index: int, circular: bool = False
) -> List[int]:
    """
    Find all cut sites for a given enzyme in the DNA sequence using IUPAC pattern matching.

    Args:
        dna_sequence: The DNA sequence to search
        enzyme_sequence: The recognition sequence of the enzyme (may contain IUPAC letters)
        cut_index: 0-based index within the recognition site where the enzyme cuts
        circular: If True, wrap cut positions using modulo for circular DNA (handles Type IIS edge cases)

    Returns:
        List of break positions (0-based indices)
    """
    break_positions = []
    seq_len = len(dna_sequence)
    
    # Convert IUPAC site to regex pattern and compile with lookahead for overlapping matches
    regex_pattern = f"(?={iupac_to_regex(enzyme_sequence)})"
    pattern = re.compile(regex_pattern, flags=re.IGNORECASE)
    
    # Find all overlapping matches
    for match in pattern.finditer(dna_sequence):
        # Calculate break position using cut_index
        break_pos = match.start() + cut_index
        
        # For circular DNA, validate and wrap cut positions
        # This handles Type IIS enzymes that cut outside the recognition site
        if circular:
            # Wrap around using modulo
            break_pos = break_pos % seq_len
        else:
            # For linear DNA, ensure cut position is within valid range
            if break_pos < 0 or break_pos > seq_len:
                continue  # Skip invalid cut positions
        
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


def find_cut_positions_linear(seq: str, enzyme_name: str, db: Dict[str, Dict[str, any]], circular: bool = False) -> List[int]:
    """
    Find all cut positions for a single enzyme in a DNA sequence.
    
    Args:
        seq: DNA sequence to search
        enzyme_name: Name of the enzyme
        db: Enzyme database
        circular: If True, treat DNA as circular (wraps cut positions with modulo)
        
    Returns:
        List of break positions (0-based indices)
    """
    if enzyme_name not in db:
        return []
    
    enzyme_info = db[enzyme_name]
    enzyme_sequence = enzyme_info["sequence"]
    cut_index = enzyme_info["cut_index"]
    
    return find_cut_sites(seq, enzyme_sequence, cut_index, circular=circular)


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
        "--seq", required=False, help="DNA sequence (string) or path to FASTA/text file"
    )
    parser.add_argument(
        "--enz", required=False, nargs='+', help="One or more enzyme names (see available enzymes in output)"
    )
    parser.add_argument(
        "--circular", action="store_true", default=False,
        help="Treat DNA as circular (default: linear)"
    )
    parser.add_argument(
        "--circular_single_cut_linearizes", action="store_true", default=False,
        help="In circular mode, one cut yields two fragments instead of one intact circle (default: False)"
    )
    
    # Restriction map options
    parser.add_argument(
        "--print-map", action="store_true", default=False,
        help="Print restriction map after digestion results"
    )
    parser.add_argument(
        "--print-map-only", action="store_true", default=False,
        help="Print only the restriction map (skip fragment table)"
    )
    parser.add_argument(
        "--map-width", type=int, default=80,
        help="Width of the restriction map in characters (default: 80)"
    )
    parser.add_argument(
        "--map-ticks", type=int, default=10,
        help="Number of tick marks on the map scale (default: 10)"
    )
    parser.add_argument(
        "--map-min-hits", type=int, default=1,
        help="Minimum number of cuts to show an enzyme (default: 1)"
    )
    parser.add_argument(
        "--map-group-by", choices=["enzyme", "position"], default="enzyme",
        help="Group cuts by enzyme or by position (default: enzyme)"
    )
    parser.add_argument(
        "--map-show-overhangs", action="store_true", default=False,
        help="Show overhang type labels in the map"
    )
    parser.add_argument(
        "--map-show-sites", action="store_true", default=False,
        help="Show recognition sequences in the map"
    )
    parser.add_argument(
        "--map-circular-origin", type=int, default=0,
        help="Origin position for circular DNA map (default: 0)"
    )
    
    # Gel simulation options
    parser.add_argument(
        "--simulate-gel", action="store_true", default=False,
        help="Print ASCII gel simulation after digestion results"
    )
    parser.add_argument(
        "--gel-only", action="store_true", default=False,
        help="Print only the gel simulation (skip digestion results and map)"
    )
    parser.add_argument(
        "--gel-percent", type=float, default=1.0,
        help="Agarose gel percentage (0.7-3.0, default: 1.0)"
    )
    parser.add_argument(
        "--gel-length", type=int, default=24,
        help="Gel height in rows (default: 24)"
    )
    parser.add_argument(
        "--gel-width", type=int, default=80,
        help="Gel width in characters (default: 80)"
    )
    parser.add_argument(
        "--gel-lane-gap", type=int, default=3,
        help="Spacing between lanes (default: 3)"
    )
    parser.add_argument(
        "--gel-ladder", type=str, default="1kb",
        help="Ladder preset: 100bp, 1kb, broad (default: 1kb)"
    )
    parser.add_argument(
        "--gel-topology", choices=["auto", "linearized", "native"], default="auto",
        help="Topology rendering mode (default: auto)"
    )
    parser.add_argument(
        "--gel-merge-threshold", type=int, default=20,
        help="Merge bands closer than this size (bp, default: 20)"
    )
    parser.add_argument(
        "--gel-smear", choices=["none", "light", "heavy"], default="none",
        help="Gel smear effect (default: none)"
    )
    parser.add_argument(
        "--gel-dye-front", type=float, default=0.85,
        help="Dye front position, 0-1 fraction down gel (default: 0.85)"
    )
    parser.add_argument(
        "--lanes-config", type=str, default=None,
        help="JSON file or string defining multiple lanes for gel"
    )
    
    # Graphics output options
    parser.add_argument(
        "--out-svg", type=str, default=None,
        help="Output SVG file path for plasmid map"
    )
    parser.add_argument(
        "--out-svg-linear", type=str, default=None,
        help="Output SVG file path for linear restriction map"
    )
    parser.add_argument(
        "--out-svg-fragments", type=str, default=None,
        help="Output SVG file path for fragment diagram"
    )
    parser.add_argument(
        "--png", action="store_true", default=False,
        help="Also export PNG alongside each SVG (requires cairosvg)"
    )
    parser.add_argument(
        "--theme", choices=["light", "dark"], default="light",
        help="Color theme for graphics (default: light)"
    )
    parser.add_argument(
        "--title", type=str, default=None,
        help="Custom title for graphics output"
    )
    parser.add_argument(
        "--show-sites", action="store_true", default=False,
        help="Include recognition sequences in graphics labels"
    )
    parser.add_argument(
        "--hide-overhangs", action="store_true", default=False,
        help="Hide overhang type badges in graphics"
    )
    parser.add_argument(
        "--origin", type=int, default=0,
        help="Origin position for circular plasmid map (default: 0)"
    )
    parser.add_argument(
        "--svg-width", type=int, default=None,
        help="Override width for linear map and fragment diagram (default: 900)"
    )
    parser.add_argument(
        "--svg-height", type=int, default=None,
        help="Override height for linear map and fragment diagram"
    )
    
    # Sequence output options
    parser.add_argument(
        "--include-seqs", action="store_true", default=False,
        help="Include DNA sequences in text output"
    )
    parser.add_argument(
        "--seq-context", type=int, default=0,
        help="Limit sequence display to N bases at each end (0 = show full sequence)"
    )
    parser.add_argument(
        "--fasta-out", type=str, default=None,
        help="Output fragment sequences to a FASTA file"
    )
    
    # Ligation compatibility options
    parser.add_argument(
        "--compatibility", action="store_true", default=False,
        help="Compute sticky-end compatibility for fragment ends"
    )
    parser.add_argument(
        "--compat-summary", choices=["pairs", "matrix", "detailed"], default="pairs",
        help="Format for compatibility output (default: pairs)"
    )
    parser.add_argument(
        "--require-directional", action="store_true", default=False,
        help="Filter to directional pairs only (non-palindromic overhangs)"
    )
    parser.add_argument(
        "--include-blunt", action="store_true", default=False,
        help="Include blunt-blunt as compatible (default: False)"
    )
    parser.add_argument(
        "--min-overhang", type=int, default=1,
        help="Minimum overhang length for sticky classification (default: 1)"
    )
    parser.add_argument(
        "--json-out", type=str, default=None,
        help="Output compatibility results to JSON file"
    )
    
    # Theoretical compatibility options (no digest required)
    parser.add_argument(
        "--theoretical-enzymes", type=str, default=None,
        help="Comma-separated enzyme names for theoretical compatibility analysis (no sequence required)"
    )
    parser.add_argument(
        "--theoretical-all", action="store_true", default=False,
        help="Compute theoretical compatibility matrix for all enzymes in database"
    )
    parser.add_argument(
        "--format", choices=["pairs", "matrix", "detailed"], default="pairs",
        help="Format for theoretical compatibility output (default: pairs)"
    )
    
    # Export options
    parser.add_argument(
        "--export-genbank", type=str, default=None,
        help="Export digest to GenBank file at specified path"
    )
    parser.add_argument(
        "--export-csv", type=str, default=None,
        help="Export digest to CSV files with specified prefix (creates prefix_fragments.csv and prefix_cuts.csv)"
    )
    parser.add_argument(
        "--gb-definition", type=str, default=None,
        help="GenBank DEFINITION line text (default: 'Restriction digest export')"
    )
    parser.add_argument(
        "--source", type=str, default="synthetic DNA",
        help="GenBank SOURCE organism name (default: 'synthetic DNA')"
    )
    parser.add_argument(
        "--topology", type=str, choices=["linear", "circular"], default=None,
        help="Override topology for export (defaults to current run mode)"
    )
    
    # Cloning planner options
    parser.add_argument(
        "--plan-cloning", type=str, default=None,
        help="Run multi-step cloning planner with spec file (JSON or YAML)"
    )
    parser.add_argument(
        "--max-steps", type=int, default=3,
        help="Maximum number of cloning steps to search (default: 3)"
    )
    parser.add_argument(
        "--prefer-typeIIS", action="store_true", default=False,
        help="Prefer Type IIS (Golden Gate) one-pot assemblies"
    )
    parser.add_argument(
        "--avoid-enzymes", type=str, default=None,
        help="Comma-separated list of enzymes to avoid in planning"
    )
    parser.add_argument(
        "--allow-enzymes", type=str, default=None,
        help="Comma-separated list of allowed enzymes (restrict search space)"
    )
    parser.add_argument(
        "--frame-check", action="store_true", default=False,
        help="Enforce reading-frame continuity across junctions"
    )
    parser.add_argument(
        "--export-plan", type=str, default=None,
        help="Export plan to directory (creates GenBank/CSV for each step)"
    )
    parser.add_argument(
        "--print-gels", action="store_true", default=False,
        help="Include simulated gels for each step in plan output"
    )

    args = parser.parse_args()

    try:
        # Load enzyme database first (needed for all modes)
        ENZYMES = load_enzyme_database()
        
        # ====================================================================
        # CLONING PLANNER MODE
        # ====================================================================
        if args.plan_cloning:
            print("=" * 80)
            print("MULTI-STEP CLONING PLANNER")
            print("=" * 80)
            print()
            
            # Import planner modules
            try:
                from planner import plan_from_spec
                from planner_utils import (
                    load_json_or_yaml, validate_spec, 
                    format_plan_detailed, format_plan_json
                )
            except ImportError as e:
                print(f"Error: Could not import planner modules: {e}")
                sys.exit(1)
            
            # Load specification
            print(f"Loading cloning specification: {args.plan_cloning}")
            try:
                spec = load_json_or_yaml(args.plan_cloning)
            except (FileNotFoundError, ValueError) as e:
                print(f"Error loading spec: {e}")
                sys.exit(1)
            
            # Validate specification
            valid, error_msg = validate_spec(spec)
            if not valid:
                print(f"Error: Invalid specification - {error_msg}")
                sys.exit(1)
            
            print("✓ Specification loaded and validated")
            print()
            
            # Parse options
            planner_options = {}
            
            if args.avoid_enzymes:
                planner_options['avoid_enzymes'] = [e.strip() for e in args.avoid_enzymes.split(',')]
                print(f"Avoiding enzymes: {', '.join(planner_options['avoid_enzymes'])}")
            
            if args.allow_enzymes:
                planner_options['allow_enzymes'] = [e.strip() for e in args.allow_enzymes.split(',')]
                print(f"Allowed enzymes: {', '.join(planner_options['allow_enzymes'])}")
            
            planner_options['prefer_typeIIS'] = args.prefer_typeIIS
            if args.prefer_typeIIS:
                print("Preference: Type IIS (Golden Gate) assemblies")
            
            planner_options['frame_check'] = args.frame_check
            if args.frame_check:
                print("Frame checking: Enabled")
            
            planner_options['avoid_internal_cuts'] = spec.get('constraints', {}).get('avoid_internal_cuts', True)
            planner_options['min_overhang'] = spec.get('constraints', {}).get('min_overhang', 4)
            planner_options['beam_width'] = 10  # Can be made configurable
            
            print(f"Max steps: {args.max_steps}")
            print()
            
            # Run planner
            print("Searching for optimal cloning plan...")
            print()
            
            plan = plan_from_spec(
                spec=spec,
                enzyme_db=ENZYMES,
                max_steps=args.max_steps,
                options=planner_options
            )
            
            # Display results
            if not plan.feasible:
                print("⚠ No feasible plan found")
                print(f"Reason: {plan.reason}")
                print()
                print("Suggestions:")
                print("  - Increase --max-steps to allow more complex strategies")
                print("  - Relax constraints in the specification")
                print("  - Check that enzymes cut at appropriate sites")
                sys.exit(1)
            
            print("✓ Plan found!")
            print()
            
            # Format and display plan
            from planner import format_plan_summary
            
            summary = format_plan_summary(plan, show_gels=args.print_gels)
            print(summary)
            
            # Show detailed protocol if requested
            if args.export_plan or args.print_gels:
                detailed = format_plan_detailed(plan, export_dir=args.export_plan)
                print()
                print(detailed)
            
            # Export plan if requested
            if args.export_plan:
                print()
                print(f"Exporting plan to: {args.export_plan}")
                
                from planner import export_plan_to_files
                export_plan_to_files(plan, args.export_plan, ENZYMES)
                
                # Also export JSON representation
                import os
                json_path = os.path.join(args.export_plan, "plan.json")
                with open(json_path, 'w') as f:
                    f.write(format_plan_json(plan))
                print(f"✓ Plan JSON saved to: {json_path}")
            
            # Exit after planning
            sys.exit(0)
        
        # ====================================================================
        # THEORETICAL COMPATIBILITY MODE (no sequence required)
        # ====================================================================
        if args.theoretical_enzymes or args.theoretical_all:
            print("=" * 80)
            print("THEORETICAL ENZYME COMPATIBILITY ANALYSIS (NO DIGEST)")
            print("=" * 80)
            print()
            
            # Determine which enzymes to analyze
            if args.theoretical_all:
                enzyme_names = list(ENZYMES.keys())
                print(f"Analyzing all {len(enzyme_names)} enzymes in database...")
                if len(enzyme_names) > 100:
                    print(f"Warning: Large output expected ({len(enzyme_names)} enzymes)")
                print()
            else:
                # Parse comma-separated enzyme names
                enzyme_names = [name.strip() for name in args.theoretical_enzymes.split(',')]
                print(f"Analyzing {len(enzyme_names)} enzymes: {', '.join(enzyme_names)}")
                print()
            
            # Build TheoreticalEnd objects for each enzyme
            theoretical_ends = []
            skipped = []
            
            for name in enzyme_names:
                try:
                    end = theoretical_end_from_enzyme(name, ENZYMES)
                    theoretical_ends.append(end)
                except (KeyError, ValueError) as e:
                    skipped.append((name, str(e)))
            
            # Report skipped enzymes
            if skipped:
                print("Skipped enzymes:")
                for name, reason in skipped:
                    print(f"  - {name}: {reason}")
                print()
            
            # Filter by min-overhang if specified
            if args.min_overhang > 1:
                before_count = len(theoretical_ends)
                theoretical_ends = [e for e in theoretical_ends if e.k >= args.min_overhang or e.k == 0]
                filtered_count = before_count - len(theoretical_ends)
                if filtered_count > 0:
                    print(f"Filtered out {filtered_count} enzymes with overhang < {args.min_overhang} bp")
                    print()
            
            # Calculate compatibility
            results = calculate_theoretical_compatibility(
                ends=theoretical_ends,
                include_blunt=args.include_blunt,
                min_overhang=args.min_overhang,
                require_directional=args.require_directional
            )
            
            # Format and display output
            if args.format == "pairs":
                output = format_theoretical_pairs(results)
                print(output)
            elif args.format == "matrix":
                output = format_theoretical_matrix(results, theoretical_ends)
                print(output)
            elif args.format == "detailed":
                output = format_theoretical_detailed(results)
                print(output)
            
            # Save JSON output if requested
            if args.json_out:
                json_data = []
                for a, b, directional, reason in results:
                    entry = {
                        "enzyme_a": a.enzyme,
                        "enzyme_b": b.enzyme,
                        "overhang_type": a.overhang_type,
                        "k": a.k,
                        "template_a": a.sticky_template,
                        "template_b": b.sticky_template,
                        "compatible": True,
                        "directional": directional,
                        "reason": reason,
                        "palindromic_a": a.is_palindromic,
                        "palindromic_b": b.is_palindromic
                    }
                    json_data.append(entry)
                
                with open(args.json_out, 'w') as f:
                    json.dump(json_data, f, indent=2)
                print(f"\n✓ Results saved to: {args.json_out}")
            
            # Exit after theoretical analysis
            sys.exit(0)
        
        # ====================================================================
        # REGULAR DIGEST MODE (requires sequence)
        # ====================================================================
        
        # Validate required arguments for digest mode
        if not args.seq:
            print("Error: --seq is required for digest mode.")
            print("For theoretical compatibility without a sequence, use --theoretical-enzymes or --theoretical-all")
            sys.exit(2)
        
        # Allow --lanes-config without --enz if lanes define their own enzymes
        if not args.enz and not args.lanes_config:
            print("Error: No enzymes specified.")
            print("Usage: python sim.py --seq <sequence> --enz <enzyme1> [enzyme2] ...")
            print("Or use --lanes-config to define enzymes per lane")
            sys.exit(2)

        # Create normalized lookup dictionary for case-insensitive matching
        normalized_lookup = {}
        for name in ENZYMES.keys():
            norm_name = normalize(name)
            if norm_name not in normalized_lookup:
                normalized_lookup[norm_name] = []
            normalized_lookup[norm_name].append(name)
        
        available_names = list(ENZYMES.keys())

        # Validate and normalize enzyme names, handling duplicates
        validated_enzymes = []
        validated_display_names = []
        enzyme_count = {}  # Track how many times each enzyme appears
        
        # Only validate enzymes if --enz was provided
        if args.enz:
            for enzyme_name in args.enz:
                norm_name = normalize(enzyme_name)
                if norm_name in normalized_lookup:
                    variants = normalized_lookup[norm_name]
                    if len(variants) == 1:
                        actual_enzyme = variants[0]
                        validated_enzymes.append(actual_enzyme)
                        
                        # Track duplicate count (for display names)
                        if actual_enzyme not in enzyme_count:
                            enzyme_count[actual_enzyme] = 1
                        else:
                            enzyme_count[actual_enzyme] += 1
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
            
            # Build a list of display names in order
            validated_display_names = []
            current_counts = {}
            for enz in validated_enzymes:
                if enz not in current_counts:
                    current_counts[enz] = 1
                    validated_display_names.append(enz)
                else:
                    current_counts[enz] += 1
                    validated_display_names.append(f"{enz}#{current_counts[enz]}")

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
        
        # Only process enzymes if --enz was provided (not using lanes-config only)
        if validated_enzymes:
            for enzyme_name, display_name in zip(validated_enzymes, validated_display_names):
                enzyme_info = ENZYMES[enzyme_name]
                recognition_seq = enzyme_info["sequence"]
                cut_index = enzyme_info["cut_index"]
                overhang_type = enzyme_info["overhang_type"]
                
                print(f"Enzyme: {display_name}")
                print(f"Site:   {recognition_seq}")
                print(f"Cut @:  index {cut_index}")
                print(f"Overhang: {overhang_type}")
                
                # Find cut positions for this enzyme
                break_positions = find_cut_positions_linear(dna_sequence, enzyme_name, ENZYMES, circular=args.circular)
                cuts_by_enzyme[display_name] = break_positions
                
                # Store metadata for each cut position with display name
                for pos in break_positions:
                    if pos not in cut_metadata:
                        cut_metadata[pos] = []
                    cut_metadata[pos].append({
                        'enzyme': display_name,  # Use display name for output
                        'actual_enzyme': enzyme_name,  # Keep actual enzyme for lookups
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
        
        # Compute fragments with sequences if requested
        fragments_with_seqs = None
        if args.include_seqs or args.fasta_out:
            fragments_with_seqs = compute_fragments_with_sequences(
                dna_sequence=dna_sequence,
                cut_positions=all_cuts,
                circular=args.circular,
                circular_single_cut_linearizes=args.circular_single_cut_linearizes,
                cut_metadata=cut_metadata
            )
        
        # Build cut_events list for restriction map
        cut_events = []
        for pos in sorted(all_cuts):
            enzymes_at_pos = cut_metadata.get(pos, [])
            for enz_meta in enzymes_at_pos:
                cut_events.append({
                    'pos': pos,
                    'enzyme': enz_meta['enzyme'],
                    'site': enz_meta['site'],
                    'overhang_type': enz_meta['overhang_type']
                })
        
        # Display detailed fragment information (unless --print-map-only or --gel-only is set)
        if not args.print_map_only and not args.gel_only:
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
            
            # Display sequences if requested
            if args.include_seqs and fragments_with_seqs:
                print()
                print("=" * 80)
                print("FRAGMENT SEQUENCES")
                print("=" * 80)
                
                for idx, frag in enumerate(fragments_with_seqs):
                    print(f"\n# Fragment {idx + 1} (length: {frag.length} bp)")
                    print(f"start={frag.start_idx}  end={frag.end_idx}  mode={'circular' if args.circular else 'linear'}")
                    
                    # Display end information
                    left_end, right_end = frag.enzymes_at_ends
                    
                    if left_end:
                        left_str = f"{left_end.enzyme} ({left_end.overhang_type}"
                        if left_end.overhang_len > 0:
                            left_str += f", {left_end.overhang_len} bp: {left_end.end_bases}"
                        left_str += ")"
                    else:
                        left_str = "START"
                    
                    if right_end:
                        right_str = f"{right_end.enzyme} ({right_end.overhang_type}"
                        if right_end.overhang_len > 0:
                            right_str += f", {right_end.overhang_len} bp: {right_end.end_bases}"
                        right_str += ")"
                    else:
                        right_str = "END"
                    
                    print(f"ends: left={left_str}, right={right_str}")
                    
                    # Display sequence (with optional elision)
                    display_seq = elide_sequence(frag.sequence, args.seq_context)
                    print(f"seq: {display_seq}")
                
                print()
            
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
        
        # Display restriction map if requested (skip with --gel-only)
        if (args.print_map or args.print_map_only) and not args.gel_only:
            if args.print_map_only:
                print("=" * 80)
                print("RESTRICTION MAP")
                print("=" * 80)
            else:
                print()
                print("=" * 80)
                print("RESTRICTION MAP")
                print("=" * 80)
            
            restriction_map = build_restriction_map(
                L=len(dna_sequence),
                cut_events=cut_events,
                circular=args.circular,
                map_width=args.map_width,
                map_ticks=args.map_ticks,
                map_min_hits=args.map_min_hits,
                group_by=args.map_group_by,
                show_overhangs=args.map_show_overhangs,
                show_sites=args.map_show_sites,
                circular_origin=args.map_circular_origin
            )
            
            print(restriction_map)
            print()
        
        # Display gel simulation if requested
        if args.simulate_gel or args.gel_only:
            # Load ladder
            try:
                ladder_bp = get_ladder(args.gel_ladder)
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            
            # Prepare lanes for gel simulation
            gel_lanes = []
            
            if args.lanes_config:
                # Load multi-lane configuration
                try:
                    # Try to parse as JSON string first
                    if args.lanes_config.strip().startswith('[') or args.lanes_config.strip().startswith('{'):
                        lanes_cfg_data = json.loads(args.lanes_config)
                    else:
                        # Try to load from file
                        with open(args.lanes_config, 'r') as f:
                            lanes_cfg_data = json.load(f)
                    
                    # Ensure it's a list
                    if isinstance(lanes_cfg_data, dict):
                        lanes_cfg_data = [lanes_cfg_data]
                    
                    # Validate schema
                    if not isinstance(lanes_cfg_data, list):
                        raise ValueError("lanes-config must be a JSON array of lane objects")
                    
                    lanes_data = lanes_cfg_data
                    
                    # Process each lane configuration
                    for lane_idx, lane_config in enumerate(lanes_data):
                        # Validate lane object
                        if not isinstance(lane_config, dict):
                            print(f"Warning: Lane {lane_idx} is not a valid object, skipping")
                            continue
                        
                        # Validate required fields
                        if 'enzymes' not in lane_config:
                            raise ValueError(f"Lane {lane_idx} is missing required field 'enzymes'")
                        if not isinstance(lane_config['enzymes'], list):
                            raise ValueError(f"Lane {lane_idx}: 'enzymes' must be an array of enzyme names")
                        
                        # Get lane parameters with fallback to global settings
                        lane_enzymes = lane_config.get('enzymes', args.enz if args.enz else [])
                        lane_label = lane_config.get('label', '+'.join(lane_enzymes) if lane_enzymes else f'Lane{lane_idx+1}')
                        lane_circular = lane_config.get('circular', args.circular)
                        lane_notes = lane_config.get('notes', '')
                        
                        # Compute fragments for this lane
                        lane_cuts_by_enzyme = {}
                        lane_cut_metadata = {}
                        
                        for enz_name in lane_enzymes:
                            if enz_name not in ENZYMES:
                                print(f"Warning: Enzyme '{enz_name}' not found, skipping in lane '{lane_label}'")
                                continue
                            
                            enz_info = ENZYMES[enz_name]
                            enz_cuts = find_cut_positions_linear(dna_sequence, enz_name, ENZYMES, circular=lane_circular)
                            lane_cuts_by_enzyme[enz_name] = enz_cuts
                            
                            for pos in enz_cuts:
                                if pos not in lane_cut_metadata:
                                    lane_cut_metadata[pos] = []
                                lane_cut_metadata[pos].append({
                                    'enzyme': enz_name,
                                    'site': enz_info['sequence'],
                                    'cut_index': enz_info['cut_index'],
                                    'overhang_type': enz_info['overhang_type']
                                })
                        
                        lane_all_cuts = merge_cut_positions(lane_cuts_by_enzyme, len(dna_sequence))
                        
                        # Compute fragments for this lane
                        lane_fragments = compute_fragments(
                            cut_positions=lane_all_cuts,
                            seq_len=len(dna_sequence),
                            circular=lane_circular,
                            circular_single_cut_linearizes=args.circular_single_cut_linearizes,
                            cut_metadata=lane_cut_metadata
                        )
                        
                        # Extract fragment sizes
                        fragment_sizes = [f['length'] for f in lane_fragments]
                        
                        # Determine topology for gel rendering
                        gel_topology = args.gel_topology
                        
                        # Handle circular special cases for gel rendering
                        if lane_circular:
                            num_cuts = len(lane_all_cuts)
                            if gel_topology == "auto":
                                if num_cuts == 0:
                                    # No cuts: render SC/OC forms
                                    gel_topology = "native"
                                elif num_cuts == 1:
                                    # One cut: linearized
                                    gel_topology = "linearized"
                        
                        gel_lanes.append({
                            'label': lane_label,
                            'fragments': fragment_sizes,
                            'topology': gel_topology,
                            'circular': lane_circular,
                            'notes': lane_notes
                        })
                
                except json.JSONDecodeError as e:
                    print("Error: Invalid JSON format in lanes-config")
                    print(f"  Details: {e}")
                    print()
                    print("Expected format: [{'label': 'Lane1', 'enzymes': ['EcoRI']}, {'label': 'Lane2', 'enzymes': ['HindIII']}]")
                    print()
                    print("Required fields:")
                    print("  - label: string (lane name)")
                    print("  - enzymes: array of enzyme names")
                    print()
                    print("Optional fields:")
                    print("  - circular: boolean (default: false)")
                    print("  - notes: string (additional information)")
                    sys.exit(1)
                except FileNotFoundError:
                    print(f"Error: Could not find lanes-config file: {args.lanes_config}")
                    sys.exit(1)
                except ValueError as e:
                    print(f"Error: {e}")
                    print()
                    print("Expected format: [{'label': 'Lane1', 'enzymes': ['EcoRI']}]")
                    sys.exit(1)
            
            else:
                # Use current digest as single lane
                fragment_sizes = [f['length'] for f in fragments]
                lane_label = '+'.join(validated_display_names)
                
                # Determine topology for gel rendering
                gel_topology = args.gel_topology
                
                # Handle circular special cases
                if args.circular:
                    num_cuts = len(all_cuts)
                    if gel_topology == "auto":
                        if num_cuts == 0:
                            gel_topology = "native"
                        elif num_cuts == 1:
                            gel_topology = "linearized"
                
                gel_lanes.append({
                    'label': lane_label,
                    'fragments': fragment_sizes,
                    'topology': gel_topology,
                    'circular': args.circular,
                    'notes': ''
                })
            
            # Generate gel simulation
            gel_output = simulate_gel(
                lanes=gel_lanes,
                ladder_bp=ladder_bp,
                gel_percent=args.gel_percent,
                gel_length=args.gel_length,
                gel_width=args.gel_width,
                lane_gap=args.gel_lane_gap,
                merge_threshold=args.gel_merge_threshold,
                smear=args.gel_smear,
                dye_front=args.gel_dye_front
            )
            
            # Print gel
            if args.gel_only:
                # Clear previous output and only show gel
                pass  # Already skipped above with --gel-only
            else:
                print()
            
            print(gel_output)
        
        # Generate graphics outputs if requested
        if args.out_svg or args.out_svg_linear or args.out_svg_fragments:
            # Build cut list with metadata for graphics
            cuts_for_graphics = []
            for pos in sorted(all_cuts):
                enzymes_at_pos = cut_metadata.get(pos, [])
                for enz_meta in enzymes_at_pos:
                    cuts_for_graphics.append({
                        'pos': pos,
                        'enzyme': enz_meta['enzyme'],
                        'site': enz_meta['site'],
                        'overhang_type': enz_meta['overhang_type']
                    })
            
            # Determine title
            graphics_title = args.title if args.title else (
                "Plasmid" if args.circular else "DNA"
            )
            
            # Generate plasmid map SVG
            if args.out_svg:
                try:
                    svg_content = render_plasmid_map(
                        L=len(dna_sequence),
                        cuts=cuts_for_graphics,
                        title=graphics_title,
                        origin=args.origin,
                        show_sites=args.show_sites,
                        show_overhangs=not args.hide_overhangs,
                        theme=args.theme
                    )
                    
                    with open(args.out_svg, 'w') as f:
                        f.write(svg_content)
                    print(f"\n✓ Plasmid map saved to: {args.out_svg}")
                    
                    # Generate PNG if requested
                    if args.png:
                        png_path = args.out_svg.rsplit('.', 1)[0] + '.png'
                        try:
                            svg_to_png(svg_content, png_path)
                            print(f"✓ PNG saved to: {png_path}")
                        except ImportError as e:
                            print(f"Warning: Could not generate PNG - {e}")
                        except Exception:
                            print("Warning: Could not generate PNG - Cairo library not found.")
                            print("  The SVG was saved successfully. To enable PNG export, install Cairo:")
                            print("  macOS: brew install cairo")
                            print("  Ubuntu: sudo apt-get install libcairo2-dev")
                
                except Exception as e:
                    print(f"Error generating plasmid map SVG: {e}")
            
            # Generate linear map SVG
            if args.out_svg_linear:
                try:
                    linear_title = args.title if args.title else "Restriction Map"
                    svg_width = args.svg_width if args.svg_width else 900
                    svg_height = args.svg_height if args.svg_height else 180
                    
                    svg_content = render_linear_map(
                        L=len(dna_sequence),
                        cuts=cuts_for_graphics,
                        title=linear_title,
                        width=svg_width,
                        height=svg_height,
                        show_sites=args.show_sites,
                        show_overhangs=not args.hide_overhangs,
                        theme=args.theme
                    )
                    
                    with open(args.out_svg_linear, 'w') as f:
                        f.write(svg_content)
                    print(f"✓ Linear map saved to: {args.out_svg_linear}")
                    
                    # Generate PNG if requested
                    if args.png:
                        png_path = args.out_svg_linear.rsplit('.', 1)[0] + '.png'
                        try:
                            svg_to_png(svg_content, png_path)
                            print(f"✓ PNG saved to: {png_path}")
                        except ImportError as e:
                            print(f"Warning: Could not generate PNG - {e}")
                        except Exception:
                            print("Warning: Could not generate PNG - Cairo library not found.")
                            print("  The SVG was saved successfully. To enable PNG export, install Cairo:")
                            print("  macOS: brew install cairo")
                
                except Exception as e:
                    print(f"Error generating linear map SVG: {e}")
            
            # Generate fragment diagram SVG
            if args.out_svg_fragments:
                try:
                    frag_title = args.title if args.title else "Fragments"
                    svg_width = args.svg_width if args.svg_width else 900
                    svg_height = args.svg_height if args.svg_height else 140
                    
                    svg_content = render_fragment_diagram(
                        fragments=fragments,
                        L=len(dna_sequence),
                        title=frag_title,
                        width=svg_width,
                        height=svg_height,
                        theme=args.theme,
                        annotate_sizes=True
                    )
                    
                    with open(args.out_svg_fragments, 'w') as f:
                        f.write(svg_content)
                    print(f"✓ Fragment diagram saved to: {args.out_svg_fragments}")
                    
                    # Generate PNG if requested
                    if args.png:
                        png_path = args.out_svg_fragments.rsplit('.', 1)[0] + '.png'
                        try:
                            svg_to_png(svg_content, png_path)
                            print(f"✓ PNG saved to: {png_path}")
                        except ImportError as e:
                            print(f"Warning: Could not generate PNG - {e}")
                        except Exception:
                            print("Warning: Could not generate PNG - Cairo library not found.")
                            print("  The SVG was saved successfully. To enable PNG export, install Cairo:")
                            print("  macOS: brew install cairo")
                
                except Exception as e:
                    print(f"Error generating fragment diagram SVG: {e}")
        
        # Generate FASTA output if requested
        if args.fasta_out and fragments_with_seqs:
            try:
                with open(args.fasta_out, 'w') as fasta_file:
                    for idx, frag in enumerate(fragments_with_seqs):
                        # Build FASTA header with fragment information
                        frag_id = f"frag_{idx+1:03d}"
                        
                        # Left end info
                        left_end, right_end = frag.enzymes_at_ends
                        if left_end and left_end.overhang_len > 0:
                            left_info = f"{left_end.enzyme}:{left_end.overhang_type[0]}p:{left_end.overhang_len}:{left_end.end_bases}"
                        elif left_end:
                            left_info = f"{left_end.enzyme}:blunt:0"
                        else:
                            left_info = "START"
                        
                        # Right end info
                        if right_end and right_end.overhang_len > 0:
                            right_info = f"{right_end.enzyme}:{right_end.overhang_type[0]}p:{right_end.overhang_len}:{right_end.end_bases}"
                        elif right_end:
                            right_info = f"{right_end.enzyme}:blunt:0"
                        else:
                            right_info = "END"
                        
                        # Write header
                        header = f">{frag_id}|len={frag.length}|start={frag.start_idx}|end={frag.end_idx}|left={left_info}|right={right_info}"
                        if frag.wraps:
                            header += "|wraps=True"
                        fasta_file.write(header + "\n")
                        
                        # Write sequence (wrapped at 80 characters)
                        seq = frag.sequence
                        for i in range(0, len(seq), 80):
                            fasta_file.write(seq[i:i+80] + "\n")
                
                print(f"\n✓ Fragment sequences saved to: {args.fasta_out}")
                print(f"  Total fragments: {len(fragments_with_seqs)}")
            
            except Exception as e:
                print(f"\nError writing FASTA file: {e}")
        
        # Perform ligation compatibility analysis if requested
        if args.compatibility:
            print()
            print("=" * 80)
            print("LIGATION COMPATIBILITY ANALYSIS")
            print("=" * 80)
            print()
            
            # Extract fragment ends for compatibility analysis
            try:
                fragment_ends = extract_fragment_ends_for_ligation(
                    dna_sequence=dna_sequence,
                    cut_positions=all_cuts,
                    circular=args.circular,
                    circular_single_cut_linearizes=args.circular_single_cut_linearizes,
                    cut_metadata=cut_metadata
                )
                
                if not fragment_ends:
                    print("No fragment ends to analyze (no cuts or intact circle).")
                    print()
                else:
                    # Calculate compatibility
                    compat_results = calculate_compatibility(
                        ends=fragment_ends,
                        include_blunt=args.include_blunt,
                        min_overhang=args.min_overhang,
                        require_directional=args.require_directional
                    )
                    
                    # Display results based on format
                    if args.compat_summary == "pairs":
                        output = format_pairs_output(compat_results)
                        print(output)
                    elif args.compat_summary == "matrix":
                        output = format_matrix_output(compat_results, fragment_ends)
                        print(output)
                    elif args.compat_summary == "detailed":
                        output = format_detailed_output(compat_results)
                        print(output)
                    
                    # Export to JSON if requested
                    if args.json_out:
                        try:
                            export_to_json(compat_results, args.json_out)
                            print(f"✓ Compatibility results saved to: {args.json_out}")
                            print(f"  Total compatible pairs: {len(compat_results)}")
                            print()
                        except Exception as e:
                            print(f"Error writing JSON file: {e}")
            
            except Exception as e:
                print(f"Error during compatibility analysis: {e}")
                import traceback
                traceback.print_exc()
        
        # Export to GenBank if requested
        if args.export_genbank:
            try:
                # Build cuts list for export
                export_cuts = []
                for pos in sorted(all_cuts):
                    enzymes_at_pos = cut_metadata.get(pos, [])
                    for enz_meta in enzymes_at_pos:
                        # Use centralized function to compute overhang metadata
                        end_meta = compute_end_metadata(
                            dna=dna_sequence,
                            cut_pos=pos,
                            recognition_site=enz_meta['site'],
                            cut_index=enz_meta['cut_index'],
                            overhang_type=enz_meta['overhang_type'],
                            is_left_end=True,  # Doesn't matter for just getting length
                            circular=args.circular
                        )
                        
                        export_cuts.append({
                            'pos': pos,
                            'enzyme': enz_meta['enzyme'],
                            'recognition_site': enz_meta['site'],
                            'cut_index': enz_meta['cut_index'],
                            'overhang_type': enz_meta['overhang_type'],
                            'overhang_len': end_meta['overhang_len']
                        })
                
                # Determine topology for export
                export_topology = args.topology if args.topology else ("circular" if args.circular else "linear")
                
                # Determine definition
                gb_definition = args.gb_definition if args.gb_definition else "Restriction digest export"
                
                export_genbank(
                    sequence=dna_sequence,
                    cuts=export_cuts,
                    fragments=fragments,
                    path=args.export_genbank,
                    topology=export_topology,
                    definition=gb_definition,
                    organism=args.source
                )
            
            except Exception as e:
                print(f"Error exporting GenBank file: {e}")
                import traceback
                traceback.print_exc()
        
        # Export to CSV if requested
        if args.export_csv:
            try:
                # Build cuts list for export (using centralized function)
                export_cuts = []
                for pos in sorted(all_cuts):
                    enzymes_at_pos = cut_metadata.get(pos, [])
                    for enz_meta in enzymes_at_pos:
                        # Use centralized function to compute overhang metadata
                        end_meta = compute_end_metadata(
                            dna=dna_sequence,
                            cut_pos=pos,
                            recognition_site=enz_meta['site'],
                            cut_index=enz_meta['cut_index'],
                            overhang_type=enz_meta['overhang_type'],
                            is_left_end=True,  # Doesn't matter for just getting length
                            circular=args.circular
                        )
                        
                        export_cuts.append({
                            'pos': pos,
                            'enzyme': enz_meta['enzyme'],
                            'recognition_site': enz_meta['site'],
                            'cut_index': enz_meta['cut_index'],
                            'overhang_type': enz_meta['overhang_type'],
                            'overhang_len': end_meta['overhang_len']
                        })
                
                # Determine topology for export
                export_topology = args.topology if args.topology else ("circular" if args.circular else "linear")
                
                export_csv(
                    prefix=args.export_csv,
                    cuts=export_cuts,
                    fragments=fragments,
                    topology=export_topology,
                    dna_sequence=dna_sequence
                )
            
            except Exception as e:
                print(f"Error exporting CSV files: {e}")
                import traceback
                traceback.print_exc()

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)


if __name__ == "__main__":
    main()
