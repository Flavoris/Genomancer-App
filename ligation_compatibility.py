#!/usr/bin/env python3
"""
Ligation Compatibility Analysis Module
Analyzes sticky-end compatibility between fragment ends produced by restriction enzymes.
"""

from typing import List, Tuple, Dict, Optional, NamedTuple
import json


# ============================================================================
# DATA STRUCTURES
# ============================================================================

class EndInfo(NamedTuple):
    """Information about a fragment end with ligation compatibility data."""
    enzyme: str
    overhang_type: str      # "5' overhang" | "3' overhang" | "Blunt"
    overhang_len: int       # 0 for blunt
    sticky_seq: str         # 5'→3' sequence of the single-stranded overhang
    polarity: str           # "left" or "right" end of fragment
    fragment_id: int        # Fragment index
    position: int           # Position in the sequence where this end is


class CompatibilityResult(NamedTuple):
    """Result of compatibility check between two ends."""
    end_a: EndInfo
    end_b: EndInfo
    compatible: bool
    directional: bool
    gc_percent_a: float
    gc_percent_b: float
    tm_a: float
    tm_b: float
    note: str


# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def revcomp(s: str) -> str:
    """
    Calculate reverse complement of a DNA sequence.
    
    Args:
        s: DNA sequence (uppercase)
        
    Returns:
        Reverse complement sequence
    """
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]


def calculate_gc_percent(seq: str) -> float:
    """
    Calculate GC percentage of a sequence.
    
    Args:
        seq: DNA sequence
        
    Returns:
        GC percentage (0-100)
    """
    if not seq:
        return 0.0
    
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq)) * 100.0


def calculate_tm(seq: str) -> float:
    """
    Calculate rough melting temperature using Wallace rule.
    
    Tm ≈ 2*(A+T) + 4*(G+C)
    
    This is only accurate for short oligonucleotides (<14 nt).
    For longer sequences, more complex models would be needed.
    
    Args:
        seq: DNA sequence
        
    Returns:
        Estimated Tm in °C
    """
    if not seq:
        return 0.0
    
    seq_upper = seq.upper()
    at_count = seq_upper.count('A') + seq_upper.count('T')
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    
    return 2.0 * at_count + 4.0 * gc_count


def are_compatible(
    a: EndInfo,
    b: EndInfo,
    include_blunt: bool = False,
    min_overhang: int = 1
) -> bool:
    """
    Check if two fragment ends are compatible for ligation.
    
    Rules:
    1. Sticky↔sticky only; blunt↔blunt optional if include_blunt
    2. Sticky↔blunt is always incompatible
    3. Length must match
    4. Sequences must be complementary: sticky_seq_A == revcomp(sticky_seq_B)
    5. Overhang types must be the same (both 5', both 3', or both blunt)
    
    Args:
        a: First end
        b: Second end
        include_blunt: If True, blunt-blunt is compatible
        min_overhang: Minimum overhang length for sticky classification
        
    Returns:
        True if compatible, False otherwise
    """
    # Check if both are blunt
    if a.overhang_len == 0 and b.overhang_len == 0:
        return include_blunt
    
    # Check if one is blunt and one is sticky (always incompatible)
    if (a.overhang_len == 0) != (b.overhang_len == 0):
        return False
    
    # Check minimum overhang requirement
    if a.overhang_len < min_overhang or b.overhang_len < min_overhang:
        return False
    
    # Check length match
    if a.overhang_len != b.overhang_len:
        return False
    
    # Check overhang type match (both must be 5' or both must be 3')
    if a.overhang_type != b.overhang_type:
        return False
    
    # Check complementarity
    return a.sticky_seq.upper() == revcomp(b.sticky_seq).upper()


def is_directional(a: EndInfo, b: EndInfo) -> bool:
    """
    Check if a compatible pair enforces directionality.
    
    A pair is directional if the sticky sequence is not palindromic
    (i.e., sticky_seq != revcomp(sticky_seq)).
    
    Directional pairs prevent vector self-ligation and enforce insert orientation.
    
    Args:
        a: First end
        b: Second end
        
    Returns:
        True if directional, False if palindromic
    """
    s = a.sticky_seq.upper()
    return s != revcomp(s)


def calculate_compatibility(
    ends: List[EndInfo],
    include_blunt: bool = False,
    min_overhang: int = 1,
    require_directional: bool = False
) -> List[CompatibilityResult]:
    """
    Calculate compatibility between all pairs of fragment ends.
    
    Args:
        ends: List of fragment ends
        include_blunt: Include blunt-blunt compatibility
        min_overhang: Minimum overhang length
        require_directional: Only return directional pairs
        
    Returns:
        List of CompatibilityResult objects
    """
    results = []
    
    # Check all pairs
    for i in range(len(ends)):
        for j in range(i + 1, len(ends)):
            end_a = ends[i]
            end_b = ends[j]
            
            # Check compatibility
            compatible = are_compatible(end_a, end_b, include_blunt, min_overhang)
            
            if not compatible:
                continue
            
            # Check directionality
            directional = is_directional(end_a, end_b)
            
            # Skip if directionality is required but not met
            if require_directional and not directional:
                continue
            
            # Calculate heuristics
            gc_a = calculate_gc_percent(end_a.sticky_seq)
            gc_b = calculate_gc_percent(end_b.sticky_seq)
            tm_a = calculate_tm(end_a.sticky_seq)
            tm_b = calculate_tm(end_b.sticky_seq)
            
            # Generate note
            note = _generate_compatibility_note(end_a, end_b, directional)
            
            results.append(CompatibilityResult(
                end_a=end_a,
                end_b=end_b,
                compatible=True,
                directional=directional,
                gc_percent_a=gc_a,
                gc_percent_b=gc_b,
                tm_a=tm_a,
                tm_b=tm_b,
                note=note
            ))
    
    return results


def _generate_compatibility_note(end_a: EndInfo, end_b: EndInfo, directional: bool) -> str:
    """Generate a descriptive note about the compatibility."""
    if end_a.overhang_len == 0:
        return "Blunt-blunt ligation"
    
    dir_note = "directional" if directional else "non-directional (palindromic)"
    return f"{end_a.overhang_type}, {end_a.overhang_len} bp overhang, {dir_note}"


# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

def format_pairs_output(results: List[CompatibilityResult]) -> str:
    """
    Format compatibility results as list of compatible pairs.
    
    Args:
        results: List of CompatibilityResult objects
        
    Returns:
        Formatted string
    """
    if not results:
        return "No compatible pairs found.\n"
    
    lines = []
    lines.append("=" * 80)
    lines.append("LIGATION COMPATIBILITY ANALYSIS - COMPATIBLE PAIRS")
    lines.append("=" * 80)
    lines.append("")
    
    for i, result in enumerate(results, 1):
        a = result.end_a
        b = result.end_b
        
        # Format end labels
        a_label = f"[frag{a.fragment_id}:{a.polarity}] {a.enzyme}"
        b_label = f"[frag{b.fragment_id}:{b.polarity}] {b.enzyme}"
        
        if a.overhang_len > 0:
            lines.append(f"Compatible pair #{i} (k={a.overhang_len}):")
            lines.append(f"  {a_label} {a.overhang_type}: {a.sticky_seq}")
            lines.append(f"  ↔ {b_label} {b.overhang_type}: {b.sticky_seq} (revcomp match)")
            
            # Add heuristics
            avg_gc = (result.gc_percent_a + result.gc_percent_b) / 2
            avg_tm = (result.tm_a + result.tm_b) / 2
            dir_str = "YES" if result.directional else "NO (palindromic)"
            lines.append(f"  directionality: {dir_str}, GC%: {avg_gc:.1f}, Tm≈{avg_tm:.0f}°C")
        else:
            lines.append(f"Compatible pair #{i} (blunt ends):")
            lines.append(f"  {a_label} ↔ {b_label}")
            lines.append(f"  Note: Blunt-blunt ligation (requires ligase, not sticky)")
        
        lines.append("")
    
    lines.append(f"Total compatible pairs: {len(results)}")
    lines.append("")
    
    return "\n".join(lines)


def format_matrix_output(results: List[CompatibilityResult], all_ends: List[EndInfo]) -> str:
    """
    Format compatibility results as a matrix.
    
    Args:
        results: List of CompatibilityResult objects
        all_ends: All fragment ends for the matrix
        
    Returns:
        Formatted string
    """
    lines = []
    lines.append("=" * 80)
    lines.append("LIGATION COMPATIBILITY ANALYSIS - COMPATIBILITY MATRIX")
    lines.append("=" * 80)
    lines.append("")
    
    if not all_ends:
        lines.append("No ends to display.")
        return "\n".join(lines)
    
    # Build compatibility lookup
    compat_lookup = set()
    for r in results:
        # Create unique key for pair
        key1 = (r.end_a.fragment_id, r.end_a.polarity, r.end_b.fragment_id, r.end_b.polarity)
        key2 = (r.end_b.fragment_id, r.end_b.polarity, r.end_a.fragment_id, r.end_a.polarity)
        compat_lookup.add(key1)
        compat_lookup.add(key2)
    
    # Create end labels
    end_labels = []
    for end in all_ends:
        label = f"F{end.fragment_id}{end.polarity[0].upper()}"
        end_labels.append(label)
    
    # Print header
    header = "    " + " ".join(f"{label:>4}" for label in end_labels)
    lines.append(header)
    lines.append("    " + "-" * (5 * len(end_labels)))
    
    # Print matrix
    for i, end_i in enumerate(all_ends):
        row = f"{end_labels[i]:<4}"
        for j, end_j in enumerate(all_ends):
            if i == j:
                # Same end
                cell = "  · "
            else:
                key = (end_i.fragment_id, end_i.polarity, end_j.fragment_id, end_j.polarity)
                if key in compat_lookup:
                    # Check if it's blunt or sticky
                    is_blunt = end_i.overhang_len == 0 and end_j.overhang_len == 0
                    cell = "  • " if is_blunt else "  ✓ "
                else:
                    cell = "  . "
            row += cell
        lines.append(row)
    
    lines.append("")
    lines.append("Legend:")
    lines.append("  ✓  = compatible (sticky ends)")
    lines.append("  •  = compatible (blunt ends)")
    lines.append("  .  = incompatible")
    lines.append("  ·  = same end")
    lines.append("")
    lines.append(f"Total compatible pairs: {len(results)}")
    lines.append("")
    
    return "\n".join(lines)


def format_detailed_output(results: List[CompatibilityResult]) -> str:
    """
    Format compatibility results with detailed information.
    
    Args:
        results: List of CompatibilityResult objects
        
    Returns:
        Formatted string
    """
    if not results:
        return "No compatible pairs found.\n"
    
    lines = []
    lines.append("=" * 80)
    lines.append("LIGATION COMPATIBILITY ANALYSIS - DETAILED REPORT")
    lines.append("=" * 80)
    lines.append("")
    
    for i, result in enumerate(results, 1):
        a = result.end_a
        b = result.end_b
        
        lines.append(f"Pair #{i}:")
        lines.append(f"  End A: Fragment {a.fragment_id} ({a.polarity}), Enzyme: {a.enzyme}")
        lines.append(f"         Type: {a.overhang_type}, Length: {a.overhang_len} bp")
        if a.overhang_len > 0:
            lines.append(f"         Sequence: 5'-{a.sticky_seq}-3'")
            lines.append(f"         GC%: {result.gc_percent_a:.1f}%, Tm: {result.tm_a:.1f}°C")
        
        lines.append(f"  End B: Fragment {b.fragment_id} ({b.polarity}), Enzyme: {b.enzyme}")
        lines.append(f"         Type: {b.overhang_type}, Length: {b.overhang_len} bp")
        if b.overhang_len > 0:
            lines.append(f"         Sequence: 5'-{b.sticky_seq}-3'")
            lines.append(f"         GC%: {result.gc_percent_b:.1f}%, Tm: {result.tm_b:.1f}°C")
        
        lines.append(f"  Compatible: {result.compatible}")
        lines.append(f"  Directional: {result.directional}")
        lines.append(f"  Note: {result.note}")
        lines.append("")
    
    lines.append(f"Total compatible pairs: {len(results)}")
    lines.append("")
    
    return "\n".join(lines)


def export_to_json(results: List[CompatibilityResult], output_path: str) -> None:
    """
    Export compatibility results to JSON file.
    
    Args:
        results: List of CompatibilityResult objects
        output_path: Path to output JSON file
    """
    json_data = []
    
    for result in results:
        entry = {
            "end_a": {
                "fragment_id": result.end_a.fragment_id,
                "polarity": result.end_a.polarity,
                "enzyme": result.end_a.enzyme,
                "overhang_type": result.end_a.overhang_type,
                "overhang_len": result.end_a.overhang_len,
                "sticky_seq": result.end_a.sticky_seq,
                "gc_percent": result.gc_percent_a,
                "tm": result.tm_a,
                "position": result.end_a.position
            },
            "end_b": {
                "fragment_id": result.end_b.fragment_id,
                "polarity": result.end_b.polarity,
                "enzyme": result.end_b.enzyme,
                "overhang_type": result.end_b.overhang_type,
                "overhang_len": result.end_b.overhang_len,
                "sticky_seq": result.end_b.sticky_seq,
                "gc_percent": result.gc_percent_b,
                "tm": result.tm_b,
                "position": result.end_b.position
            },
            "compatible": result.compatible,
            "directional": result.directional,
            "note": result.note
        }
        json_data.append(entry)
    
    with open(output_path, 'w') as f:
        json.dump(json_data, f, indent=2)


# ============================================================================
# ENZYME PAIR ANALYSIS (THEORETICAL)
# ============================================================================

def analyze_enzyme_pair_theoretical(
    enzyme_a: Dict[str, any],
    enzyme_b: Dict[str, any]
) -> Dict[str, any]:
    """
    Analyze theoretical compatibility between two enzymes based on metadata only.
    
    This doesn't require a digest - it just checks if the enzymes would produce
    compatible ends based on their recognition sites and cut patterns.
    
    Args:
        enzyme_a: Dictionary with keys: name, site, cut_index, overhang_type
        enzyme_b: Dictionary with keys: name, site, cut_index, overhang_type
        
    Returns:
        Dictionary with compatibility information
    """
    # Calculate overhang lengths
    len_a = abs(2 * enzyme_a['cut_index'] - len(enzyme_a['site']))
    len_b = abs(2 * enzyme_b['cut_index'] - len(enzyme_b['site']))
    
    # Basic checks
    if len_a == 0 or len_b == 0:
        compatible = (len_a == 0 and len_b == 0)
        note = "Both blunt" if compatible else "One blunt, one sticky - incompatible"
        return {
            "compatible": compatible,
            "enzyme_a": enzyme_a['name'],
            "enzyme_b": enzyme_b['name'],
            "overhang_len_a": len_a,
            "overhang_len_b": len_b,
            "note": note
        }
    
    # Check overhang type match
    if enzyme_a['overhang_type'] != enzyme_b['overhang_type']:
        return {
            "compatible": False,
            "enzyme_a": enzyme_a['name'],
            "enzyme_b": enzyme_b['name'],
            "overhang_len_a": len_a,
            "overhang_len_b": len_b,
            "note": "Overhang types don't match (5' vs 3')"
        }
    
    # Check length match
    if len_a != len_b:
        return {
            "compatible": False,
            "enzyme_a": enzyme_a['name'],
            "enzyme_b": enzyme_b['name'],
            "overhang_len_a": len_a,
            "overhang_len_b": len_b,
            "note": f"Overhang lengths don't match ({len_a} vs {len_b})"
        }
    
    # For same enzyme, always compatible
    if enzyme_a['name'] == enzyme_b['name']:
        return {
            "compatible": True,
            "enzyme_a": enzyme_a['name'],
            "enzyme_b": enzyme_b['name'],
            "overhang_len_a": len_a,
            "overhang_len_b": len_b,
            "note": "Same enzyme - always compatible"
        }
    
    return {
        "compatible": True,
        "enzyme_a": enzyme_a['name'],
        "enzyme_b": enzyme_b['name'],
        "overhang_len_a": len_a,
        "overhang_len_b": len_b,
        "note": "Potentially compatible (length and type match, need sequence verification)"
    }

