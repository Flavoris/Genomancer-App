#!/usr/bin/env python3
"""
Fragment Calculator Module
Handles fragment computation for both linear and circular DNA digestion.
Also includes ASCII agarose gel simulation.
"""

import math
import random
from typing import List, Dict, Tuple, Optional, NamedTuple


# ============================================================================
# DATA STRUCTURES FOR FRAGMENT SEQUENCES
# ============================================================================

class EndInfo(NamedTuple):
    """Information about a fragment end created by enzyme cutting."""
    enzyme: str
    recognition_site: str
    cut_index: int
    overhang_type: str   # "5' overhang" | "3' overhang" | "Blunt"
    end_type: str        # Same as overhang_type for consistency
    overhang_len: int
    end_bases: str       # 5'->3' bases present at that end


class Fragment(NamedTuple):
    """Complete fragment information including sequence."""
    start_idx: int       # 0-based, inclusive
    end_idx: int         # 0-based, exclusive (or wraps for circular)
    length: int
    sequence: str        # 5'->3' orientation
    enzymes_at_ends: Tuple[Optional[EndInfo], Optional[EndInfo]]  # (left, right)
    wraps: bool          # True if fragment wraps around origin in circular mode


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# IUPAC degenerate base mapping (imported from sim.py)
IUPAC_MAP = {
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
        Regex pattern where IUPAC letters are expanded to character classes
    """
    site_upper = site.upper()
    return "".join(IUPAC_MAP.get(ch, ch) for ch in site_upper)


def slice_circular(seq: str, start: int, end: int) -> str:
    """
    Extract a sequence slice that may wrap around in circular DNA.
    
    Args:
        seq: Full DNA sequence
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        
    Returns:
        Sliced sequence, wrapping around if start >= end
    """
    n = len(seq)
    if n == 0:
        return ""
    
    start = start % n
    end = end % n
    
    if start < end:
        return seq[start:end]
    elif start > end:
        # Wrap around
        return seq[start:] + seq[:end]
    else:
        # start == end means full circle or empty
        return seq if start == 0 else ""


def elide_sequence(seq: str, context: int) -> str:
    """
    Elide a sequence to show only context bases at each end.
    
    Args:
        seq: DNA sequence to elide
        context: Number of bases to show at each end (0 = show all)
        
    Returns:
        Elided sequence string
    """
    if context <= 0 or len(seq) <= 2 * context:
        return seq
    return seq[:context] + "..." + seq[-context:]


def calculate_overhang_length(recognition_site: str, cut_index: int) -> int:
    """
    Calculate the length of the overhang based on recognition site and cut position.
    
    For Type II enzymes, the overhang length is determined by the symmetry of the cut.
    For palindromic sites, overhang = |2*cut_index - site_length|
    
    Args:
        recognition_site: The enzyme recognition sequence
        cut_index: Position where enzyme cuts within the site
        
    Returns:
        Length of overhang in bases
    """
    site_len = len(recognition_site)
    
    # For standard Type II enzymes with palindromic sites
    # The overhang is the difference between the top and bottom strand cuts
    overhang = abs(2 * cut_index - site_len)
    
    return overhang


def compute_end_metadata(
    dna: str,
    cut_pos: int,
    recognition_site: str,
    cut_index: int,
    overhang_type: str,
    is_left_end: bool,
    circular: bool = False
) -> Dict[str, any]:
    """
    SINGLE SOURCE OF TRUTH for overhang metadata computation.
    
    This function computes all end-related metadata for a fragment end created by
    an enzyme cut. All output functions (console, CSV, GenBank, JSON, graphics)
    MUST use this function to ensure consistency.
    
    Args:
        dna: Full DNA sequence
        cut_pos: Position where the enzyme cuts (0-based)
        recognition_site: Recognition sequence (e.g., "GAATTC" for EcoRI)
        cut_index: Cut position within the site (0-based)
        overhang_type: "5' overhang" | "3' overhang" | "Blunt"
        is_left_end: True if this is the left (5') end of a fragment
        circular: True if the sequence is circular
        
    Returns:
        Dictionary with keys:
            - overhang_len: int (length of overhang in bp)
            - sticky_seq: str (5'→3' sequence of protruding nucleotides)
            - polarity: str ("left" | "right")
            - end_bases: str (same as sticky_seq, for display)
            - overhang_type: str (passed through)
    
    Examples:
        >>> # EcoRI cuts G^AATTC to produce 5' AATT overhang (4 bp)
        >>> # Both left and right ends show the same canonical sticky sequence: AATT
        >>> compute_end_metadata("AAAGAATTCGGG", 4, "GAATTC", 1, "5' overhang", is_left_end=True)
        {'overhang_len': 4, 'sticky_seq': 'AATT', 'polarity': 'left', 'end_bases': 'AATT', 'overhang_type': "5' overhang"}
        
        >>> compute_end_metadata("AAAGAATTCGGG", 4, "GAATTC", 1, "5' overhang", is_left_end=False)
        {'overhang_len': 4, 'sticky_seq': 'AATT', 'polarity': 'right', 'end_bases': 'AATT', 'overhang_type': "5' overhang"}
        
        >>> # PstI cuts CTGCA^G to produce 3' TGCA overhang (4 bp)
        >>> # Both ends show the canonical sticky sequence: TGCA
        >>> compute_end_metadata("AAACTGCAGTTT", 8, "CTGCAG", 5, "3' overhang", is_left_end=True)
        {'overhang_len': 4, 'sticky_seq': 'TGCA', 'polarity': 'left', 'end_bases': 'TGCA', 'overhang_type': "3' overhang"}
        
        >>> # SmaI cuts CCC^GGG to produce blunt end (0 bp)
        >>> compute_end_metadata("AAACCCGGGAAA", 6, "CCCGGG", 3, "Blunt", is_left_end=True)
        {'overhang_len': 0, 'sticky_seq': '', 'polarity': 'left', 'end_bases': '', 'overhang_type': 'Blunt'}
    """
    # Calculate overhang length
    overhang_len = calculate_overhang_length(recognition_site, cut_index)
    
    # Determine polarity based on fragment end
    polarity = "left" if is_left_end else "right"
    
    # Handle blunt ends
    if overhang_type == "Blunt" or overhang_len == 0:
        return {
            'overhang_len': 0,
            'sticky_seq': '',
            'polarity': polarity,
            'end_bases': '',
            'overhang_type': overhang_type
        }
    
    # Extract the CANONICAL sticky sequence (the protruding nucleotides)
    # This is the same for both left and right ends - we always report
    # the actual sticky overhang sequence that protrudes, not the recessed part.
    #
    # For 5' overhang (e.g., EcoRI: G^AATTC):
    #   Top strand:    5'---G AATT---3'
    #   Bottom strand: 3'---CTTAA G---5'
    #   Canonical sticky sequence: AATT (the 4 bases after the cut)
    #   Both left and right fragment ends show: AATT
    #
    # For 3' overhang (e.g., PstI: CTGCA^G):
    #   Top strand:    5'---CTGCA G---3'
    #   Bottom strand: 3'---G ACGTC---5'
    #   Canonical sticky sequence: TGCA (the 4 bases before the cut)
    #   Both left and right fragment ends show: TGCA
    
    seq_len = len(dna)
    
    if overhang_type == "5' overhang":
        # For 5' overhang: use dna[cut_pos:cut_pos + overhang_len]
        # This gets the protruding nucleotides (AATT for EcoRI)
        start = cut_pos
        end = cut_pos + overhang_len
    else:  # 3' overhang
        # For 3' overhang: use dna[cut_pos - overhang_len:cut_pos]
        # This gets the protruding nucleotides (TGCA for PstI)
        start = cut_pos - overhang_len
        end = cut_pos
    
    # Extract the sequence (handle circular wrapping if needed)
    if circular:
        sticky_seq = slice_circular(dna, start, end)
    else:
        # Clamp to sequence boundaries
        start = max(0, min(start, seq_len))
        end = max(0, min(end, seq_len))
        sticky_seq = dna[start:end]
    
    return {
        'overhang_len': overhang_len,
        'sticky_seq': sticky_seq,
        'polarity': polarity,
        'end_bases': sticky_seq,  # Same as sticky_seq for display
        'overhang_type': overhang_type
    }


def extract_end_bases(
    seq: str,
    cut_pos: int,
    enzyme_meta: Dict,
    is_left_end: bool,
    circular: bool = False
) -> str:
    """
    Extract the bases at a fragment end based on the overhang type.
    
    IMPORTANT: Returns the CANONICAL sticky sequence that is the same for both
    left and right fragment ends created by the same cut. This ensures consistency
    with compute_end_metadata().
    
    Args:
        seq: Full DNA sequence
        cut_pos: Position of the cut
        enzyme_meta: Dictionary with enzyme information
        is_left_end: True if this is the left (5') end of the fragment
        circular: True if sequence is circular
        
    Returns:
        String of bases at the end (5'->3' orientation) - canonical sticky sequence
    """
    overhang_type = enzyme_meta.get('overhang_type', 'Blunt')
    recognition_site = enzyme_meta.get('site', '')
    cut_index = enzyme_meta.get('cut_index', 0)
    
    # Calculate overhang length
    overhang_len = calculate_overhang_length(recognition_site, cut_index)
    
    if overhang_type == "Blunt" or overhang_len == 0:
        return ""
    
    # Extract the CANONICAL sticky sequence (same for both left and right ends)
    # This matches the behavior of compute_end_metadata()
    seq_len = len(seq)
    
    if overhang_type == "5' overhang":
        # For 5' overhang: canonical sequence is after the cut
        # Example: EcoRI G^AATTC produces "AATT" for both ends
        start = cut_pos
        end = cut_pos + overhang_len
    else:  # 3' overhang
        # For 3' overhang: canonical sequence is before the cut
        # Example: PstI CTGCA^G produces "TGCA" for both ends
        start = cut_pos - overhang_len
        end = cut_pos
    
    # Handle circular wrapping
    if circular:
        return slice_circular(seq, start, end)
    else:
        # Clamp to sequence boundaries
        start = max(0, min(start, seq_len))
        end = max(0, min(end, seq_len))
        return seq[start:end]


def extract_sticky_seq(
    seq: str,
    cut_pos: int,
    enzyme_meta: Dict,
    is_left_end: bool,
    circular: bool = False
) -> str:
    """
    Extract the sticky sequence for ligation compatibility analysis.
    
    IMPORTANT: Returns the CANONICAL sticky sequence (the protruding nucleotides)
    that is the same for both left and right fragment ends created by the same cut.
    This ensures consistency with compute_end_metadata().
    
    The sticky sequence represents the single-stranded overhang in 5'→3' orientation.
    For ligation compatibility, two ends are compatible if their sticky sequences
    are reverse complements of each other.
    
    Examples:
      - EcoRI (G^AATTC): produces "AATT" for both left and right ends (5' overhang)
      - PstI (CTGCA^G): produces "TGCA" for both left and right ends (3' overhang)
      - SmaI (CCC^GGG): produces "" for both ends (blunt)
    
    Args:
        seq: Full DNA sequence
        cut_pos: Position of the cut
        enzyme_meta: Dictionary with enzyme information
        is_left_end: True if this is the left (5') end of the fragment
        circular: True if sequence is circular
        
    Returns:
        String representing the canonical sticky overhang sequence (5'→3')
    """
    overhang_type = enzyme_meta.get('overhang_type', 'Blunt')
    recognition_site = enzyme_meta.get('site', '')
    cut_index = enzyme_meta.get('cut_index', 0)
    
    # Calculate overhang length
    overhang_len = calculate_overhang_length(recognition_site, cut_index)
    
    if overhang_type == "Blunt" or overhang_len == 0:
        return ""
    
    seq_len = len(seq)
    
    # Extract the CANONICAL sticky sequence (same for both left and right ends)
    # This matches the behavior of compute_end_metadata()
    
    if overhang_type == "5' overhang":
        # For 5' overhang: canonical sequence is after the cut
        # Example: EcoRI G^AATTC produces "AATT" for both ends
        #   Top strand:    5'---G AATT---3'
        #   Bottom strand: 3'---CTTAA G---5'
        #   Canonical sticky sequence: AATT
        start = cut_pos
        end = cut_pos + overhang_len
    
    else:  # 3' overhang
        # For 3' overhang: canonical sequence is before the cut
        # Example: PstI CTGCA^G produces "TGCA" for both ends
        #   Top strand:    5'---CTGCA G---3'
        #   Bottom strand: 3'---G ACGTC---5'
        #   Canonical sticky sequence: TGCA
        start = cut_pos - overhang_len
        end = cut_pos
    
    # Handle circular wrapping
    if circular:
        return slice_circular(seq, start, end)
    else:
        # Clamp to sequence boundaries
        start = max(0, min(start, seq_len))
        end = max(0, min(end, seq_len))
        return seq[start:end]


def build_end_info(
    enzyme_meta: Dict,
    seq: str,
    cut_pos: int,
    is_left_end: bool,
    circular: bool = False
) -> EndInfo:
    """
    Build EndInfo object for a fragment end using the centralized compute_end_metadata function.
    
    IMPORTANT: This function uses compute_end_metadata() as the single source of truth
    for all overhang calculations. All other functions (extract_end_bases, extract_sticky_seq)
    have been synchronized to use the same logic.
    
    Args:
        enzyme_meta: Dictionary with enzyme metadata
        seq: Full DNA sequence
        cut_pos: Position of the cut
        is_left_end: True if this is the left end
        circular: True if sequence is circular
        
    Returns:
        EndInfo object
    """
    enzyme = enzyme_meta['enzyme']
    recognition_site = enzyme_meta['site']
    cut_index = enzyme_meta['cut_index']
    overhang_type = enzyme_meta['overhang_type']
    
    # Use the centralized function for consistency
    metadata = compute_end_metadata(
        dna=seq,
        cut_pos=cut_pos,
        recognition_site=recognition_site,
        cut_index=cut_index,
        overhang_type=overhang_type,
        is_left_end=is_left_end,
        circular=circular
    )
    
    return EndInfo(
        enzyme=enzyme,
        recognition_site=recognition_site,
        cut_index=cut_index,
        overhang_type=overhang_type,
        end_type=overhang_type,
        overhang_len=metadata['overhang_len'],
        end_bases=metadata['end_bases']
    )


# ============================================================================
# FRAGMENT COMPUTATION (LEGACY)
# ============================================================================

def compute_fragments(
    cut_positions: List[int],
    seq_len: int,
    circular: bool = False,
    circular_single_cut_linearizes: bool = False,
    cut_metadata: Dict[int, List[Dict]] = None
) -> List[Dict]:
    """
    Compute DNA fragments after restriction enzyme cutting.
    
    Args:
        cut_positions: List of cut position indices [0, seq_len)
        seq_len: Length of the DNA sequence
        circular: If True, treat DNA as circular; if False, treat as linear
        circular_single_cut_linearizes: If True and circular, one cut yields two fragments
        cut_metadata: Dictionary mapping cut position to list of enzyme metadata dicts
                     Each dict contains: {enzyme, site, cut_index, overhang_type}
        
    Returns:
        List of fragment dictionaries with keys:
            - index: Fragment number (0-based)
            - length: Fragment length in bp
            - start: Start position in sequence
            - end: End position in sequence (exclusive)
            - wraps: Boolean indicating if fragment wraps around origin (circular only)
            - boundaries: Dict with left_cut and right_cut information
    
    Raises:
        ValueError: If seq_len is 0 or negative
    """
    if seq_len <= 0:
        raise ValueError(f"Sequence length must be positive, got {seq_len}")
    
    # Default empty metadata
    if cut_metadata is None:
        cut_metadata = {}
    
    # Normalize and deduplicate cut positions
    ps = sorted(set(p % seq_len for p in cut_positions))
    n = len(ps)
    
    # Build position-to-enzyme mapping with metadata
    pos_to_enzymes = {}
    for pos in ps:
        # Get all enzymes that cut at this position (handling modulo)
        enzymes_at_pos = []
        for orig_pos in cut_positions:
            normalized = orig_pos % seq_len
            if normalized == pos and orig_pos in cut_metadata:
                enzymes_at_pos.extend(cut_metadata[orig_pos])
        
        # Deduplicate by enzyme name
        seen = set()
        unique_enzymes = []
        for enz_meta in enzymes_at_pos:
            if enz_meta['enzyme'] not in seen:
                seen.add(enz_meta['enzyme'])
                unique_enzymes.append(enz_meta)
        
        pos_to_enzymes[pos] = unique_enzymes
    
    if not circular:
        # Linear mode
        if n == 0:
            # No cuts - return entire sequence
            return [{
                'index': 0,
                'length': seq_len,
                'start': 0,
                'end': seq_len,
                'wraps': False,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': None
                }
            }]
        
        fragments = []
        
        # First fragment: from sequence start to first cut
        fragments.append({
            'index': 0,
            'length': ps[0],
            'start': 0,
            'end': ps[0],
            'wraps': False,
            'boundaries': {
                'left_cut': None,  # Sequence start
                'right_cut': {
                    'pos': ps[0],
                    'enzymes': pos_to_enzymes.get(ps[0], [])
                }
            }
        })
        
        # Middle fragments: between consecutive cuts
        for i in range(n - 1):
            start = ps[i]
            end = ps[i + 1]
            fragments.append({
                'index': i + 1,
                'length': end - start,
                'start': start,
                'end': end,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': start,
                        'enzymes': pos_to_enzymes.get(start, [])
                    },
                    'right_cut': {
                        'pos': end,
                        'enzymes': pos_to_enzymes.get(end, [])
                    }
                }
            })
        
        # Last fragment: from last cut to sequence end
        fragments.append({
            'index': n,
            'length': seq_len - ps[-1],
            'start': ps[-1],
            'end': seq_len,
            'wraps': False,
            'boundaries': {
                'left_cut': {
                    'pos': ps[-1],
                    'enzymes': pos_to_enzymes.get(ps[-1], [])
                },
                'right_cut': None  # Sequence end
            }
        })
        
        return fragments
    
    else:
        # Circular mode
        if n == 0:
            # No cuts - intact circle
            return [{
                'index': 0,
                'length': seq_len,
                'start': 0,
                'end': 0,
                'wraps': True,
                'boundaries': {
                    'left_cut': None,
                    'right_cut': None
                }
            }]
        
        if n == 1:
            # One cut
            if circular_single_cut_linearizes:
                # Two fragments: linearized plasmid split at the cut
                p0 = ps[0]
                return [
                    {
                        'index': 0,
                        'length': seq_len - p0,
                        'start': p0,
                        'end': seq_len,
                        'wraps': False,
                        'boundaries': {
                            'left_cut': {
                                'pos': p0,
                                'enzymes': pos_to_enzymes.get(p0, [])
                            },
                            'right_cut': {
                                'pos': p0,
                                'enzymes': pos_to_enzymes.get(p0, [])
                            }
                        }
                    },
                    {
                        'index': 1,
                        'length': p0,
                        'start': 0,
                        'end': p0,
                        'wraps': False,
                        'boundaries': {
                            'left_cut': {
                                'pos': p0,
                                'enzymes': pos_to_enzymes.get(p0, [])
                            },
                            'right_cut': {
                                'pos': p0,
                                'enzymes': pos_to_enzymes.get(p0, [])
                            }
                        }
                    }
                ]
            else:
                # One fragment - intact circle (single cut doesn't fragment)
                return [{
                    'index': 0,
                    'length': seq_len,
                    'start': 0,
                    'end': 0,
                    'wraps': True,
                    'boundaries': {
                        'left_cut': None,
                        'right_cut': None
                    }
                }]
        
        # Multiple cuts in circular mode
        fragments = []
        
        # Fragments between consecutive cuts
        for i in range(n - 1):
            start = ps[i]
            end = ps[i + 1]
            fragments.append({
                'index': i,
                'length': end - start,
                'start': start,
                'end': end,
                'wraps': False,
                'boundaries': {
                    'left_cut': {
                        'pos': start,
                        'enzymes': pos_to_enzymes.get(start, [])
                    },
                    'right_cut': {
                        'pos': end,
                        'enzymes': pos_to_enzymes.get(end, [])
                    }
                }
            })
        
        # Wrap-around fragment: from last cut, around origin, to first cut
        start = ps[-1]
        end = ps[0]
        length = (seq_len - start) + end
        fragments.append({
            'index': n - 1,
            'length': length,
            'start': start,
            'end': end,
            'wraps': True,
            'boundaries': {
                'left_cut': {
                    'pos': start,
                    'enzymes': pos_to_enzymes.get(start, [])
                },
                'right_cut': {
                    'pos': end,
                    'enzymes': pos_to_enzymes.get(end, [])
                }
            }
        })
        
        return fragments


def validate_fragment_total(fragments: List[Dict], expected_length: int) -> bool:
    """
    Validate that fragment lengths sum to the expected total sequence length.
    
    Args:
        fragments: List of fragment dictionaries from compute_fragments
        expected_length: Expected total length
        
    Returns:
        True if lengths sum correctly, False otherwise
    """
    total = sum(f['length'] for f in fragments)
    return total == expected_length


# ============================================================================
# FRAGMENT COMPUTATION WITH SEQUENCES
# ============================================================================

def compute_fragments_with_sequences(
    dna_sequence: str,
    cut_positions: List[int],
    circular: bool = False,
    circular_single_cut_linearizes: bool = False,
    cut_metadata: Dict[int, List[Dict]] = None
) -> List[Fragment]:
    """
    Compute DNA fragments with full sequence information after restriction enzyme cutting.
    
    Args:
        dna_sequence: The complete DNA sequence
        cut_positions: List of cut position indices [0, seq_len)
        circular: If True, treat DNA as circular; if False, treat as linear
        circular_single_cut_linearizes: If True and circular, one cut yields two fragments
        cut_metadata: Dictionary mapping cut position to list of enzyme metadata dicts
                     Each dict contains: {enzyme, site, cut_index, overhang_type}
        
    Returns:
        List of Fragment namedtuples with sequence and end information
    """
    seq_len = len(dna_sequence)
    
    if seq_len == 0:
        return []
    
    # Default empty metadata
    if cut_metadata is None:
        cut_metadata = {}
    
    # Normalize and deduplicate cut positions
    ps = sorted(set(p % seq_len for p in cut_positions))
    n = len(ps)
    
    # Build position-to-enzyme mapping with metadata
    pos_to_enzymes = {}
    for pos in ps:
        # Get all enzymes that cut at this position (handling modulo)
        enzymes_at_pos = []
        for orig_pos in cut_positions:
            normalized = orig_pos % seq_len
            if normalized == pos and orig_pos in cut_metadata:
                enzymes_at_pos.extend(cut_metadata[orig_pos])
        
        # Deduplicate by enzyme name
        seen = set()
        unique_enzymes = []
        for enz_meta in enzymes_at_pos:
            if enz_meta['enzyme'] not in seen:
                seen.add(enz_meta['enzyme'])
                unique_enzymes.append(enz_meta)
        
        pos_to_enzymes[pos] = unique_enzymes
    
    fragments = []
    
    if not circular:
        # ========== LINEAR MODE ==========
        if n == 0:
            # No cuts - return entire sequence
            frag = Fragment(
                start_idx=0,
                end_idx=seq_len,
                length=seq_len,
                sequence=dna_sequence,
                enzymes_at_ends=(None, None),
                wraps=False
            )
            return [frag]
        
        # First fragment: from sequence start to first cut
        frag_seq = dna_sequence[0:ps[0]]
        right_enzymes = pos_to_enzymes.get(ps[0], [])
        right_end = build_end_info(right_enzymes[0], dna_sequence, ps[0], False, circular) if right_enzymes else None
        
        fragments.append(Fragment(
            start_idx=0,
            end_idx=ps[0],
            length=ps[0],
            sequence=frag_seq,
            enzymes_at_ends=(None, right_end),
            wraps=False
        ))
        
        # Middle fragments: between consecutive cuts
        for i in range(n - 1):
            start = ps[i]
            end = ps[i + 1]
            frag_seq = dna_sequence[start:end]
            
            left_enzymes = pos_to_enzymes.get(start, [])
            right_enzymes = pos_to_enzymes.get(end, [])
            
            left_end = build_end_info(left_enzymes[0], dna_sequence, start, True, circular) if left_enzymes else None
            right_end = build_end_info(right_enzymes[0], dna_sequence, end, False, circular) if right_enzymes else None
            
            fragments.append(Fragment(
                start_idx=start,
                end_idx=end,
                length=end - start,
                sequence=frag_seq,
                enzymes_at_ends=(left_end, right_end),
                wraps=False
            ))
        
        # Last fragment: from last cut to sequence end
        start = ps[-1]
        frag_seq = dna_sequence[start:seq_len]
        left_enzymes = pos_to_enzymes.get(start, [])
        left_end = build_end_info(left_enzymes[0], dna_sequence, start, True, circular) if left_enzymes else None
        
        fragments.append(Fragment(
            start_idx=start,
            end_idx=seq_len,
            length=seq_len - start,
            sequence=frag_seq,
            enzymes_at_ends=(left_end, None),
            wraps=False
        ))
        
    else:
        # ========== CIRCULAR MODE ==========
        if n == 0:
            # No cuts - intact circle
            frag = Fragment(
                start_idx=0,
                end_idx=0,
                length=seq_len,
                sequence=dna_sequence,
                enzymes_at_ends=(None, None),
                wraps=True
            )
            return [frag]
        
        if n == 1:
            # One cut
            if circular_single_cut_linearizes:
                # Two fragments: linearized plasmid split at the cut
                p0 = ps[0]
                enzymes_at_p0 = pos_to_enzymes.get(p0, [])
                end_info = build_end_info(enzymes_at_p0[0], dna_sequence, p0, True, circular) if enzymes_at_p0 else None
                
                # Fragment 1: from cut to end of sequence
                frag1_seq = dna_sequence[p0:seq_len]
                fragments.append(Fragment(
                    start_idx=p0,
                    end_idx=seq_len,
                    length=seq_len - p0,
                    sequence=frag1_seq,
                    enzymes_at_ends=(end_info, end_info),
                    wraps=False
                ))
                
                # Fragment 2: from start to cut
                frag2_seq = dna_sequence[0:p0]
                fragments.append(Fragment(
                    start_idx=0,
                    end_idx=p0,
                    length=p0,
                    sequence=frag2_seq,
                    enzymes_at_ends=(end_info, end_info),
                    wraps=False
                ))
            else:
                # One fragment - intact circle (single cut doesn't fragment)
                frag = Fragment(
                    start_idx=0,
                    end_idx=0,
                    length=seq_len,
                    sequence=dna_sequence,
                    enzymes_at_ends=(None, None),
                    wraps=True
                )
                return [frag]
        else:
            # Multiple cuts in circular mode
            
            # Fragments between consecutive cuts
            for i in range(n - 1):
                start = ps[i]
                end = ps[i + 1]
                frag_seq = dna_sequence[start:end]
                
                left_enzymes = pos_to_enzymes.get(start, [])
                right_enzymes = pos_to_enzymes.get(end, [])
                
                left_end = build_end_info(left_enzymes[0], dna_sequence, start, True, circular) if left_enzymes else None
                right_end = build_end_info(right_enzymes[0], dna_sequence, end, False, circular) if right_enzymes else None
                
                fragments.append(Fragment(
                    start_idx=start,
                    end_idx=end,
                    length=end - start,
                    sequence=frag_seq,
                    enzymes_at_ends=(left_end, right_end),
                    wraps=False
                ))
            
            # Wrap-around fragment: from last cut, around origin, to first cut
            start = ps[-1]
            end = ps[0]
            length = (seq_len - start) + end
            frag_seq = slice_circular(dna_sequence, start, end)
            
            left_enzymes = pos_to_enzymes.get(start, [])
            right_enzymes = pos_to_enzymes.get(end, [])
            
            left_end = build_end_info(left_enzymes[0], dna_sequence, start, True, circular) if left_enzymes else None
            right_end = build_end_info(right_enzymes[0], dna_sequence, end, False, circular) if right_enzymes else None
            
            fragments.append(Fragment(
                start_idx=start,
                end_idx=end,
                length=length,
                sequence=frag_seq,
                enzymes_at_ends=(left_end, right_end),
                wraps=True
            ))
    
    return fragments


# ============================================================================
# LIGATION COMPATIBILITY SUPPORT
# ============================================================================

def extract_fragment_ends_for_ligation(
    dna_sequence: str,
    cut_positions: List[int],
    circular: bool = False,
    circular_single_cut_linearizes: bool = False,
    cut_metadata: Dict[int, List[Dict]] = None
) -> List:
    """
    Extract all fragment ends for ligation compatibility analysis.
    
    Returns a list of EndInfo objects compatible with ligation_compatibility module.
    
    Args:
        dna_sequence: The complete DNA sequence
        cut_positions: List of cut position indices
        circular: If True, treat DNA as circular
        circular_single_cut_linearizes: If True and circular, one cut yields two fragments
        cut_metadata: Dictionary mapping cut position to enzyme metadata
        
    Returns:
        List of EndInfo objects with sticky_seq for ligation analysis
    """
    # Import the ligation compatibility EndInfo here to avoid circular imports
    from ligation_compatibility import EndInfo as LigationEndInfo
    
    seq_len = len(dna_sequence)
    
    if seq_len == 0:
        return []
    
    # Default empty metadata
    if cut_metadata is None:
        cut_metadata = {}
    
    # Normalize and deduplicate cut positions
    ps = sorted(set(p % seq_len for p in cut_positions))
    n = len(ps)
    
    # Build position-to-enzyme mapping with metadata
    pos_to_enzymes = {}
    for pos in ps:
        # Get all enzymes that cut at this position
        enzymes_at_pos = []
        for orig_pos in cut_positions:
            normalized = orig_pos % seq_len
            if normalized == pos and orig_pos in cut_metadata:
                enzymes_at_pos.extend(cut_metadata[orig_pos])
        
        # Deduplicate by enzyme name
        seen = set()
        unique_enzymes = []
        for enz_meta in enzymes_at_pos:
            if enz_meta['enzyme'] not in seen:
                seen.add(enz_meta['enzyme'])
                unique_enzymes.append(enz_meta)
        
        pos_to_enzymes[pos] = unique_enzymes
    
    ends = []
    fragment_id = 0
    
    if not circular:
        # ========== LINEAR MODE ==========
        if n == 0:
            # No cuts - no sticky ends to analyze
            return []
        
        # Process each cut position
        for i, cut_pos in enumerate(ps):
            enzymes = pos_to_enzymes.get(cut_pos, [])
            
            for enzyme_meta in enzymes:
                # Left fragment end (right side of previous fragment)
                if i > 0 or cut_pos > 0:
                    sticky_seq = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, False, circular)
                    end = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq,
                        polarity="left",
                        fragment_id=i + 1,
                        position=cut_pos
                    )
                    ends.append(end)
                
                # Right fragment end (left side of next fragment)
                if i < n - 1 or cut_pos < seq_len:
                    sticky_seq = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, True, circular)
                    end = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq,
                        polarity="right",
                        fragment_id=i,
                        position=cut_pos
                    )
                    ends.append(end)
    
    else:
        # ========== CIRCULAR MODE ==========
        if n == 0:
            # No cuts - intact circle, no ends
            return []
        
        if n == 1:
            # One cut
            if circular_single_cut_linearizes:
                # Two fragments with same enzyme at both ends
                cut_pos = ps[0]
                enzymes = pos_to_enzymes.get(cut_pos, [])
                
                for enzyme_meta in enzymes:
                    # Both ends of the linearized plasmid
                    sticky_seq_left = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, True, circular)
                    sticky_seq_right = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, False, circular)
                    
                    end_left = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq_left,
                        polarity="left",
                        fragment_id=0,
                        position=cut_pos
                    )
                    
                    end_right = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq_right,
                        polarity="right",
                        fragment_id=0,
                        position=cut_pos
                    )
                    
                    ends.extend([end_left, end_right])
            else:
                # One fragment - intact circle, no ends to analyze
                return []
        else:
            # Multiple cuts - process all boundaries
            for i, cut_pos in enumerate(ps):
                enzymes = pos_to_enzymes.get(cut_pos, [])
                
                for enzyme_meta in enzymes:
                    # Left end of fragment starting at this cut
                    sticky_seq_left = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, True, circular)
                    end_left = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq_left,
                        polarity="left",
                        fragment_id=i,
                        position=cut_pos
                    )
                    ends.append(end_left)
                    
                    # Right end of fragment ending at this cut
                    sticky_seq_right = extract_sticky_seq(dna_sequence, cut_pos, enzyme_meta, False, circular)
                    end_right = LigationEndInfo(
                        enzyme=enzyme_meta['enzyme'],
                        overhang_type=enzyme_meta['overhang_type'],
                        overhang_len=calculate_overhang_length(enzyme_meta['site'], enzyme_meta['cut_index']),
                        sticky_seq=sticky_seq_right,
                        polarity="right",
                        fragment_id=(i - 1) % n,  # Wrap around for circular
                        position=cut_pos
                    )
                    ends.append(end_right)
    
    return ends


def _scale_position_to_column(pos: int, seq_len: int, map_width: int) -> int:
    """
    Scale a position in the sequence to a column on the map.
    
    Args:
        pos: Position in sequence [0, seq_len)
        seq_len: Total sequence length
        map_width: Width of the map in characters
        
    Returns:
        Column index [0, map_width)
    """
    if seq_len <= 0:
        return 0
    return int(round(pos / seq_len * (map_width - 1)))


def _make_scale_lines(seq_len: int, map_width: int, map_ticks: int, circular: bool, origin: int) -> Tuple[str, str]:
    """
    Create the scale lines (labels and ruler) for the restriction map.
    
    Args:
        seq_len: Length of the DNA sequence
        map_width: Width of the map in characters
        map_ticks: Number of tick marks to show
        circular: Whether the sequence is circular
        origin: Origin position for circular DNA
        
    Returns:
        Tuple of (label_line, ruler_line)
    """
    # Calculate tick positions
    tick_cols = []
    tick_labels = []
    
    for i in range(map_ticks):
        # Calculate the position along the sequence for this tick
        if circular:
            # For circular, ticks represent positions relative to origin
            pos = int(round(i * seq_len / (map_ticks - 1))) if map_ticks > 1 else 0
            actual_pos = (origin + pos) % seq_len
            label = str(actual_pos) if i < map_ticks - 1 else str(seq_len)
        else:
            # For linear, evenly spaced ticks
            pos = int(round(i * seq_len / (map_ticks - 1))) if map_ticks > 1 else 0
            label = str(pos)
        
        col = int(round(i * (map_width - 1) / (map_ticks - 1))) if map_ticks > 1 else 0
        tick_cols.append(col)
        tick_labels.append(label)
    
    # Build label line
    label_chars = [' '] * map_width
    for col, label in zip(tick_cols, tick_labels):
        # Place label starting at column
        for j, ch in enumerate(label):
            idx = col + j
            if idx < map_width:
                label_chars[idx] = ch
    
    # Build ruler line
    ruler_chars = ['-'] * map_width
    for col in tick_cols:
        if col < map_width:
            ruler_chars[col] = '|'
    
    return ''.join(label_chars), ''.join(ruler_chars)


def _group_cuts_by_enzyme(cut_events: List[Dict]) -> Dict[str, List[int]]:
    """
    Group cut events by enzyme name.
    
    Args:
        cut_events: List of cut event dictionaries with keys: pos, enzyme, overhang_type, site
        
    Returns:
        Dictionary mapping enzyme name to list of positions
    """
    enzyme_cuts = {}
    for event in cut_events:
        enzyme = event['enzyme']
        pos = event['pos']
        if enzyme not in enzyme_cuts:
            enzyme_cuts[enzyme] = []
        enzyme_cuts[enzyme].append(pos)
    return enzyme_cuts


def _group_cuts_by_position(cut_events: List[Dict]) -> Dict[int, List[Dict]]:
    """
    Group cut events by position.
    
    Args:
        cut_events: List of cut event dictionaries
        
    Returns:
        Dictionary mapping position to list of enzyme info dicts
    """
    pos_cuts = {}
    for event in cut_events:
        pos = event['pos']
        if pos not in pos_cuts:
            pos_cuts[pos] = []
        pos_cuts[pos].append({
            'enzyme': event['enzyme'],
            'overhang_type': event['overhang_type'],
            'site': event['site']
        })
    return pos_cuts


def _render_track_line(
    positions: List[int],
    seq_len: int,
    map_width: int,
    glyph: str = '^',
    shared_positions: set = None
) -> str:
    """
    Render a single track line with cut markers.
    
    Args:
        positions: List of cut positions
        seq_len: Sequence length
        map_width: Width of the map
        glyph: Character to use for markers
        shared_positions: Set of positions that are shared with other enzymes
        
    Returns:
        String representing the track line
    """
    if shared_positions is None:
        shared_positions = set()
    
    track = [' '] * map_width
    
    for pos in positions:
        col = _scale_position_to_column(pos, seq_len, map_width)
        if col < map_width:
            # Use * for shared positions, otherwise use the glyph
            if pos in shared_positions:
                track[col] = '*'
            else:
                track[col] = glyph
    
    return ''.join(track)


def build_restriction_map(
    L: int,
    cut_events: List[Dict],
    *,
    circular: bool = False,
    map_width: int = 80,
    map_ticks: int = 10,
    map_min_hits: int = 1,
    group_by: str = "enzyme",
    show_overhangs: bool = False,
    show_sites: bool = False,
    circular_origin: int = 0
) -> str:
    """
    Build a text-mode restriction map showing cut sites along the sequence.
    
    Args:
        L: Length of the DNA sequence
        cut_events: List of cut event dictionaries with keys:
                   - pos: int (cut position)
                   - enzyme: str (enzyme name)
                   - overhang_type: str (from enzymes.json)
                   - site: str (recognition sequence)
        circular: If True, render as circular topology
        map_width: Width of the map in characters (default 80)
        map_ticks: Number of tick marks on the scale (default 10)
        map_min_hits: Minimum number of cuts to show an enzyme (default 1)
        group_by: "enzyme" or "position" (default "enzyme")
        show_overhangs: Include overhang type labels (default False)
        show_sites: Include recognition sequence (default False)
        circular_origin: Origin position for circular mode (default 0)
        
    Returns:
        Multi-line ASCII map string
    """
    if L <= 0:
        return "Error: Sequence length must be positive"
    
    # Deduplicate cut events by position
    unique_positions = sorted(set(event['pos'] for event in cut_events))
    n_cuts = len(unique_positions)
    
    # Build position-to-enzymes mapping
    pos_to_events = _group_cuts_by_position(cut_events)
    
    # Find shared positions (where multiple enzymes cut)
    shared_positions = set(pos for pos, events in pos_to_events.items() if len(events) > 1)
    
    lines = []
    
    # Header line
    mode = "circular" if circular else "linear"
    origin_suffix = f" (origin={circular_origin})" if circular and circular_origin != 0 else ""
    lines.append(f"Length: {L} bp   Mode: {mode}{origin_suffix}")
    
    # Scale lines
    label_line, ruler_line = _make_scale_lines(L, map_width, map_ticks, circular, circular_origin)
    lines.append(label_line)
    lines.append(ruler_line)
    
    # Handle no cuts case
    if n_cuts == 0:
        lines.append("")
        lines.append("No cut sites for selected enzymes.")
        return '\n'.join(lines)
    
    # Combined "All cuts" track
    all_track = _render_track_line(unique_positions, L, map_width, glyph='|', shared_positions=shared_positions)
    lines.append(f"All ({n_cuts}): {all_track}")
    
    # Render tracks based on grouping mode
    if group_by == "enzyme":
        # Group by enzyme
        enzyme_cuts = _group_cuts_by_enzyme(cut_events)
        
        # Filter by minimum hits and sort
        filtered_enzymes = [(enz, positions) for enz, positions in enzyme_cuts.items() 
                           if len(set(positions)) >= map_min_hits]
        
        # Sort by decreasing hit count, then alphabetically
        filtered_enzymes.sort(key=lambda x: (-len(set(x[1])), x[0]))
        
        for enzyme, positions in filtered_enzymes:
            unique_pos = sorted(set(positions))
            n_hits = len(unique_pos)
            
            # Render track
            track_line = _render_track_line(unique_pos, L, map_width, glyph='^', shared_positions=shared_positions)
            
            # Build enzyme label
            enzyme_label = f"{enzyme} ({n_hits})"
            
            # Add optional annotations
            annotations = []
            
            if show_overhangs:
                # Collect unique overhang types for this enzyme
                overhang_types = set()
                for event in cut_events:
                    if event['enzyme'] == enzyme:
                        overhang_types.add(event['overhang_type'])
                if overhang_types:
                    overhang_str = ', '.join(sorted(overhang_types))
                    annotations.append(f"[{overhang_str}]")
            
            if show_sites:
                # Get site (use first occurrence)
                site = None
                for event in cut_events:
                    if event['enzyme'] == enzyme:
                        site = event['site']
                        break
                if site:
                    # Truncate if too long
                    site_str = site if len(site) <= 16 else site[:13] + "..."
                    annotations.append(f"({site_str})")
            
            # Combine label and track
            full_label = enzyme_label
            if annotations:
                full_label += " " + " ".join(annotations)
            
            lines.append(f"{full_label}: {track_line}")
    
    elif group_by == "position":
        # Group by position
        for pos in unique_positions:
            events = pos_to_events[pos]
            
            # Build enzyme list
            enzyme_names = [e['enzyme'] for e in events]
            enzyme_str = ', '.join(enzyme_names)
            
            # Build overhang list if requested
            annotations = []
            if show_overhangs:
                overhang_types = set(e['overhang_type'] for e in events)
                overhang_str = ', '.join(sorted(overhang_types))
                annotations.append(f"[{overhang_str}]")
            
            # Render single marker at this position
            track_line = _render_track_line([pos], L, map_width, glyph='|')
            
            label = f"pos {pos}: {enzyme_str}"
            if annotations:
                label += " " + " ".join(annotations)
            
            lines.append(f"{label}: {track_line}")
    
    # Add wrap-around note for circular with 2+ cuts
    if circular and n_cuts >= 2:
        first_cut = unique_positions[0]
        last_cut = unique_positions[-1]
        wrap_len = (L - last_cut) + first_cut
        lines.append("")
        lines.append(f"wrap spans: last_cut→first_cut = {last_cut}→{first_cut} (len = {wrap_len})")
    
    # Add single cut note for circular
    if circular and n_cuts == 1:
        lines.append("")
        lines.append("Single cut on circular molecule (linearization depends on analysis mode)")
    
    return '\n'.join(lines)


# ============================================================================
# AGAROSE GEL SIMULATION
# ============================================================================

def gel_coefficients(percent: float) -> Tuple[float, float]:
    """
    Get migration coefficients (a, b) for the log-linear model based on agarose %.
    
    Migration formula: y_norm = clamp01(a + b * log10(bp))
    Higher % resolves smaller bands better, slows large ones.
    
    Small fragments migrate further (higher y_norm, larger row number).
    Large fragments migrate less (lower y_norm, smaller row number).
    
    Typical range: log10(100) = 2, log10(10000) = 4
    We want small fragments near y=0.9 and large near y=0.1
    
    Args:
        percent: Agarose concentration (0.7 to 3.0)
        
    Returns:
        Tuple of (a, b) coefficients
    """
    # Lower % gels are better for larger fragments
    if percent <= 0.8:
        return (1.8, -0.42)  # Good for large fragments
    elif percent <= 1.0:
        return (1.7, -0.40)  # General purpose
    elif percent <= 1.2:
        return (1.65, -0.38)
    elif percent <= 1.5:
        return (1.6, -0.37)
    elif percent <= 2.0:
        return (1.55, -0.36)  # Good for small fragments
    else:
        return (1.5, -0.35)  # Best for very small fragments


def calculate_migration_row(
    bp: int,
    gel_percent: float,
    gel_length: int,
    dye_front: float
) -> int:
    """
    Calculate the row index for a DNA fragment based on size.
    
    Args:
        bp: Fragment size in base pairs
        gel_percent: Agarose concentration
        gel_length: Total gel height in rows
        dye_front: Dye front position (0-1 fraction down the gel)
        
    Returns:
        Row index where the band appears
    """
    if bp <= 0:
        return 0
    
    # Get migration coefficients
    a, b = gel_coefficients(gel_percent)
    
    # Calculate normalized position (0 = top, 1 = bottom)
    y_norm = a + b * math.log10(bp)
    y_norm = max(0.0, min(1.0, y_norm))  # Clamp to [0, 1]
    
    # Convert to row index (wells at top, bands migrate down)
    max_row = int((gel_length - 1) * dye_front)
    row = int(y_norm * max_row)
    
    # Ensure band is below wells (row 0-1 are wells)
    row = max(2, row)
    
    return row


def merge_bands(
    fragments: List[int],
    merge_threshold: int,
    gel_percent: float,
    gel_length: int,
    dye_front: float
) -> Dict[int, List[int]]:
    """
    Group fragments into bands based on merge threshold.
    
    Fragments that land on the same row or have size difference <= threshold
    are merged into one band.
    
    Args:
        fragments: List of fragment sizes in bp
        merge_threshold: Size difference threshold for merging (bp)
        gel_percent: Agarose concentration
        gel_length: Gel height in rows
        dye_front: Dye front position
        
    Returns:
        Dictionary mapping row index to list of fragment sizes in that band
    """
    if not fragments:
        return {}
    
    # Sort fragments by size
    sorted_fragments = sorted(fragments)
    
    # Calculate row for each fragment
    fragment_rows = [(bp, calculate_migration_row(bp, gel_percent, gel_length, dye_front)) 
                     for bp in sorted_fragments]
    
    # Group by row and merge threshold
    bands = {}  # row -> list of bp values
    
    for bp, row in fragment_rows:
        placed = False
        
        # Check if we can merge with existing band at this row
        if row in bands:
            # Check if size difference is within threshold
            for existing_bp in bands[row]:
                if abs(bp - existing_bp) <= merge_threshold:
                    bands[row].append(bp)
                    placed = True
                    break
        
        if not placed:
            # Try adjacent rows
            for adj_row in [row - 1, row, row + 1]:
                if adj_row in bands and adj_row >= 2:
                    for existing_bp in bands[adj_row]:
                        if abs(bp - existing_bp) <= merge_threshold:
                            bands[adj_row].append(bp)
                            placed = True
                            break
                if placed:
                    break
        
        if not placed:
            # Create new band
            bands[row] = [bp]
    
    return bands


def get_band_glyph(intensity: int) -> str:
    """
    Get the character glyph for a band based on intensity (number of merged fragments).
    
    Args:
        intensity: Number of fragments merged into this band
        
    Returns:
        Character to display
    """
    if intensity >= 4:
        return '█'
    elif intensity >= 3:
        return '▮'
    elif intensity >= 2:
        return '•'
    else:
        return '·'


def add_smear(
    canvas: List[List[str]],
    lane_col: int,
    bands_in_lane: Dict[int, List[int]],
    smear: str,
    gel_length: int,
    seed: Optional[int] = None
) -> None:
    """
    Add smearing/noise to a lane to simulate gel artifacts.
    
    Args:
        canvas: 2D array of characters [row][col]
        lane_col: Column index for this lane
        bands_in_lane: Dictionary of bands (row -> bp list)
        smear: Smear intensity ("none", "light", "heavy")
        gel_length: Height of gel
        seed: Random seed for deterministic smearing
    """
    if smear == "none" or not bands_in_lane:
        return
    
    if seed is not None:
        random.seed(seed)
    
    # Find lowest band (highest row number)
    max_band_row = max(bands_in_lane.keys())
    
    # Add tailing below bands
    if smear == "light":
        # Light smear: few dots below bands
        for row in range(max_band_row + 1, min(max_band_row + 5, gel_length)):
            if random.random() < 0.15:
                if canvas[row][lane_col] == ' ':
                    canvas[row][lane_col] = '.'
    
    elif smear == "heavy":
        # Heavy smear: more dots and commas
        for row in range(max_band_row + 1, min(max_band_row + 8, gel_length)):
            rand_val = random.random()
            if rand_val < 0.25:
                if canvas[row][lane_col] == ' ':
                    canvas[row][lane_col] = '.' if rand_val < 0.20 else ','


def render_circular_topology_bands(
    L: int,
    topology: str,
    gel_percent: float,
    gel_length: int,
    dye_front: float
) -> Dict[int, List[tuple]]:
    """
    Render plasmid topology forms (supercoiled, nicked, linear) for circular DNA with 0 cuts.
    
    Args:
        L: Plasmid size in bp
        topology: "native" or "auto"
        gel_percent: Agarose concentration
        gel_length: Gel height in rows
        dye_front: Dye front position
        
    Returns:
        Dictionary mapping row to list of (label, size) tuples
    """
    # Calculate where linear form would migrate
    linear_row = calculate_migration_row(L, gel_percent, gel_length, dye_front)
    
    # Delta for topology forms (about 12% of gel height)
    delta = int(0.12 * gel_length)
    delta = max(1, delta)
    
    # Supercoiled migrates faster (higher row number = further down)
    sc_row = linear_row + delta
    # Nicked/open circular migrates slower (lower row number = stays higher)
    oc_row = linear_row - delta
    
    # Clamp to valid range (below wells, above max)
    sc_row = max(2, min(gel_length - 1, sc_row))
    oc_row = max(2, min(gel_length - 1, oc_row))
    
    bands = {}
    bands[sc_row] = [('SC', L)]
    bands[oc_row] = [('OC', L)]
    
    return bands


def simulate_gel(
    lanes: List[dict],
    ladder_bp: List[int],
    *,
    gel_percent: float = 1.0,
    gel_length: int = 24,
    gel_width: int = 80,
    lane_gap: int = 3,
    merge_threshold: int = 20,
    smear: str = "none",
    dye_front: float = 0.85
) -> str:
    """
    Simulate an ASCII agarose gel electrophoresis.
    
    Args:
        lanes: List of lane dictionaries, each with:
               - label: str (lane name)
               - fragments: List[int] (fragment sizes in bp)
               - topology: str ("auto", "linearized", "native")
               - circular: bool (whether this is circular DNA)
        ladder_bp: List of ladder fragment sizes
        gel_percent: Agarose concentration (0.7-3.0)
        gel_length: Gel height in character rows
        gel_width: Total gel width in characters
        lane_gap: Spacing between lanes
        merge_threshold: Size difference for merging bands (bp)
        smear: Smear effect ("none", "light", "heavy")
        dye_front: Dye front position (0-1, fraction down gel)
        
    Returns:
        Multi-line ASCII gel string
    """
    # Validate inputs
    if gel_percent < 0.7 or gel_percent > 3.0:
        gel_percent = max(0.7, min(3.0, gel_percent))
    
    if gel_length < 10:
        gel_length = 10
    
    if gel_width < 40:
        gel_width = 40
    
    # Calculate lane widths
    total_lanes = 1 + len(lanes)  # ladder + sample lanes
    lane_width = (gel_width - (total_lanes - 1) * lane_gap) // total_lanes
    lane_width = max(1, lane_width)
    
    # Initialize canvas
    canvas = [[' ' for _ in range(gel_width)] for _ in range(gel_length)]
    
    # Track lane info for legend
    lane_info = []
    
    # Process ladder first
    ladder_label = "Ladder"
    ladder_bands = merge_bands(ladder_bp, merge_threshold, gel_percent, gel_length, dye_front)
    ladder_col = lane_width // 2
    lane_info.append({
        'label': ladder_label,
        'col': ladder_col,
        'bands': ladder_bands,
        'fragments': sorted(ladder_bp)
    })
    
    # Process sample lanes
    for lane_idx, lane in enumerate(lanes):
        lane_label = lane.get('label', f'Lane{lane_idx+1}')
        fragments = lane.get('fragments', [])
        topology = lane.get('topology', 'auto')
        circular = lane.get('circular', False)
        
        # Calculate column position
        lane_col = (lane_idx + 1) * (lane_width + lane_gap) + lane_width // 2
        
        # Handle circular topology special cases
        if circular and len(fragments) == 1 and topology in ['auto', 'native']:
            # Check if this is a 0-cut or 1-cut scenario
            # For 0 cuts: render SC/OC forms
            # For 1 cut (auto): single linearized band
            L = fragments[0]
            
            # If we have exactly the plasmid length, it could be 0 or 1 cut
            # The caller should tell us via the fragments list
            # Assume if fragments=[L], it's either 1-cut linearized or intact
            # For native topology with 0 cuts, render SC/OC
            if topology == 'native':
                # Render supercoiled and open circular
                bands = render_circular_topology_bands(L, topology, gel_percent, gel_length, dye_front)
            else:
                # topology == 'auto': single cut linearized
                bands = merge_bands(fragments, merge_threshold, gel_percent, gel_length, dye_front)
        else:
            # Normal fragment migration
            bands = merge_bands(fragments, merge_threshold, gel_percent, gel_length, dye_front)
        
        lane_info.append({
            'label': lane_label,
            'col': lane_col,
            'bands': bands,
            'fragments': fragments,
            'topology': topology,
            'circular': circular
        })
    
    # Draw wells at top (rows 0-1)
    for info in lane_info:
        col = info['col']
        if col > 0 and col < gel_width - 1:
            canvas[0][col-1] = '┏'
            canvas[0][col] = '━'
            canvas[0][col+1] = '┓'
    
    # Draw bands for all lanes
    for info in lane_info:
        col = info['col']
        bands = info['bands']
        
        for row, bp_list in bands.items():
            if row < gel_length and col < gel_width:
                intensity = len(bp_list)
                glyph = get_band_glyph(intensity)
                canvas[row][col] = glyph
        
        # Add smear if requested
        add_smear(canvas, col, bands, smear, gel_length, seed=42)
    
    # Draw dye front line
    dye_row = int((gel_length - 1) * dye_front)
    if 2 <= dye_row < gel_length:
        for c in range(gel_width):
            if canvas[dye_row][c] == ' ':
                canvas[dye_row][c] = '~'
    
    # Convert canvas to string
    gel_lines = [''.join(row) for row in canvas]
    
    # Build legend on the right side
    legend_lines = []
    legend_lines.append(f"  Agarose: {gel_percent}%")
    legend_lines.append(f"  Dye front: {int(dye_front*100)}%")
    legend_lines.append("")
    
    for idx, info in enumerate(lane_info):
        label = info['label']
        fragments = info.get('fragments', [])
        bands = info['bands']
        
        # Count merged bands
        band_labels = []
        for row in sorted(bands.keys()):
            bp_list = bands[row]
            if len(bp_list) > 1:
                # Merged band
                avg_size = int(sum(bp_list) / len(bp_list))
                band_labels.append(f"{avg_size} ({len(bp_list)}×)")
            elif bp_list and isinstance(bp_list[0], tuple):
                # Topology label (SC, OC)
                band_labels.append(f"{bp_list[0][1]} ({bp_list[0][0]})")
            else:
                # Single band
                band_labels.append(f"{bp_list[0]}")
        
        band_str = ", ".join(band_labels) if band_labels else "no cuts"
        legend_lines.append(f"  {label}: {band_str} bp")
    
    # Combine gel and legend
    result_lines = []
    result_lines.append("=" * gel_width)
    result_lines.append("AGAROSE GEL SIMULATION")
    result_lines.append("=" * gel_width)
    result_lines.append("")
    
    # Add gel canvas
    result_lines.extend(gel_lines)
    result_lines.append("")
    
    # Add legend
    result_lines.append("Legend:")
    result_lines.extend(legend_lines)
    result_lines.append("")
    
    return '\n'.join(result_lines)

