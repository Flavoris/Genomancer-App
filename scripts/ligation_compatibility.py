#!/usr/bin/env python3
"""
Ligation Compatibility Analysis Module
Analyzes sticky-end compatibility between fragment ends produced by restriction enzymes.
Supports theoretical analysis based on enzyme metadata without sequence digests.
"""

from typing import List, Tuple, Dict, NamedTuple, Any
import json


# ============================================================================
# IUPAC SUPPORT
# ============================================================================

# IUPAC degenerate base code mapping to possible bases
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}


def revcomp_iupac(s: str) -> str:
    """
    Calculate reverse complement of a DNA sequence with IUPAC support.
    
    Args:
        s: DNA sequence with IUPAC codes
        
    Returns:
        Reverse complement sequence with IUPAC codes
    """
    comp_map = str.maketrans(
        "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
        "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"
    )
    return s.translate(comp_map)[::-1]


def iupac_compatible(a: str, b: str) -> bool:
    """
    Check if two IUPAC sequences are compatible for ligation.
    
    Two sequences are compatible if at every position, the set of possible bases
    in sequence A intersects with the complement set of possible bases in sequence B.
    
    Args:
        a: First IUPAC sequence (5'→3')
        b: Second IUPAC sequence (5'→3')
        
    Returns:
        True if compatible, False otherwise
    """
    if len(a) != len(b):
        return False
    
    b_revcomp = revcomp_iupac(b)
    
    for i, (x, yc) in enumerate(zip(a, b_revcomp)):
        x_upper = x.upper()
        yc_upper = yc.upper()
        
        # Check if the sets intersect
        if x_upper not in IUPAC or yc_upper not in IUPAC:
            return False
        
        if not (IUPAC[x_upper] & IUPAC[yc_upper]):
            return False
    
    return True


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


class TheoreticalEnd(NamedTuple):
    """Theoretical enzyme end information derived from metadata only."""
    enzyme: str
    overhang_type: str      # "5' overhang" | "3' overhang" | "Blunt"
    k: int                  # overhang length (0 if blunt)
    sticky_template: str    # 5'→3' IUPAC template of the single-stranded overhang
    is_palindromic: bool    # True if sticky_template == revcomp(sticky_template)


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
            lines.append("  Note: Blunt-blunt ligation (requires ligase, not sticky)")
        
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
# THEORETICAL ANALYSIS (NO DIGEST REQUIRED)
# ============================================================================

def derive_template(site: str, cut_index: int, overhang_type: str, k: int) -> str:
    """
    Derive the sticky end template from enzyme metadata.
    
    Returns the 5'→3' sequence of the single-stranded overhang that would be
    produced at the ligating end.
    
    Args:
        site: Recognition site (IUPAC)
        cut_index: Position where enzyme cuts on top strand
        overhang_type: "5' overhang" | "3' overhang" | "Blunt"
        k: Overhang length
        
    Returns:
        IUPAC template string (5'→3') of the overhang, or "" if blunt
    """
    if overhang_type == "Blunt" or k == 0:
        return ""
    
    site_len = len(site)
    
    if overhang_type == "5' overhang":
        # For 5' overhangs:
        # Top strand: 5'---NNNN^NNNN---3'
        # Bot strand: 3'---NNNN    ^NNNN---5'
        # The 5' overhang is the protruding part of the top strand
        # starting from the top cut position
        
        # The sticky end is from cut_index to cut_index + k on the top strand
        if cut_index + k <= site_len:
            # Overhang is within the recognition site
            template = site[cut_index:cut_index + k]
        else:
            # Overhang extends beyond the site (Type IIS)
            template = site[cut_index:] + 'N' * (k - (site_len - cut_index))
    
    elif overhang_type == "3' overhang":
        # For 3' overhangs:
        # Top strand: 5'---NNNN    ^NNNN---3'
        # Bot strand: 3'---NNNN^NNNN---5'
        # The 3' overhang is the protruding part of the bottom strand
        # We need to return the complement of what protrudes
        
        # Bottom strand cuts at position (site_len - cut_index)
        # The overhang on the bottom strand goes from (site_len - cut_index - k) to (site_len - cut_index)
        # But we want the 5'→3' sequence of that bottom strand overhang
        
        # For a 3' overhang, the single-stranded region is on the 3' end
        # The sequence we want is from (cut_index - k) to cut_index on top strand
        if cut_index - k >= 0:
            # Overhang is within the recognition site
            template = site[cut_index - k:cut_index]
        else:
            # Overhang extends beyond the site (Type IIS)
            template = 'N' * (k - cut_index) + site[:cut_index]
    
    else:
        template = ""
    
    return template.upper()


def theoretical_end_from_enzyme(name: str, db: Dict[str, Any]) -> TheoreticalEnd:
    """
    Create a TheoreticalEnd from enzyme metadata.
    
    Args:
        name: Enzyme name
        db: Enzyme database with metadata
        
    Returns:
        TheoreticalEnd object with computed template
        
    Raises:
        KeyError: If enzyme not found in database
        ValueError: If enzyme metadata is invalid
    """
    if name not in db:
        raise KeyError(f"Enzyme '{name}' not found in database")
    
    meta = db[name]
    
    # Get site and cut information
    site = meta.get("sequence", meta.get("site", ""))
    if not site:
        raise ValueError(f"Enzyme '{name}' missing recognition site")
    
    cut = meta.get("cut_index")
    if cut is None:
        raise ValueError(f"Enzyme '{name}' missing cut_index")
    
    typ = meta.get("overhang_type", "Unknown")
    if typ not in ["5' overhang", "3' overhang", "Blunt"]:
        raise ValueError(f"Enzyme '{name}' has invalid overhang_type: {typ}")
    
    # Calculate overhang length
    if typ == "Blunt":
        k = 0
    else:
        k = abs(len(site) - 2 * cut)
    
    # Derive template
    tpl = derive_template(site, cut, typ, k)
    
    # Check if palindromic
    pal = (k > 0) and (tpl.upper() == revcomp_iupac(tpl).upper())
    
    return TheoreticalEnd(name, typ, k, tpl.upper(), pal)


def are_theoretical_ends_compatible(
    a: TheoreticalEnd,
    b: TheoreticalEnd,
    include_blunt: bool = False,
    min_overhang: int = 1
) -> Tuple[bool, str]:
    """
    Check if two theoretical enzyme ends are compatible.
    
    Args:
        a: First TheoreticalEnd
        b: Second TheoreticalEnd
        include_blunt: Allow blunt-blunt compatibility
        min_overhang: Minimum overhang length
        
    Returns:
        Tuple of (compatible, reason)
    """
    # Check if both are blunt
    if a.k == 0 and b.k == 0:
        if include_blunt:
            return True, "Blunt-blunt (compatible with option)"
        else:
            return False, "Both blunt (excluded by default)"
    
    # Check if one is blunt and one is sticky
    if (a.k == 0) != (b.k == 0):
        return False, "One blunt, one sticky - incompatible"
    
    # Check minimum overhang
    if a.k < min_overhang or b.k < min_overhang:
        return False, f"Overhang shorter than minimum ({min_overhang} bp)"
    
    # Check length match
    if a.k != b.k:
        return False, f"Length mismatch: {a.k} vs {b.k} bp"
    
    # Check overhang type match
    if a.overhang_type != b.overhang_type:
        return False, f"Type mismatch: {a.overhang_type} vs {b.overhang_type}"
    
    # Check IUPAC compatibility
    if not iupac_compatible(a.sticky_template, b.sticky_template):
        return False, "IUPAC sequences not compatible"
    
    return True, "Compatible"


def calculate_theoretical_compatibility(
    ends: List[TheoreticalEnd],
    include_blunt: bool = False,
    min_overhang: int = 1,
    require_directional: bool = False
) -> List[Tuple[TheoreticalEnd, TheoreticalEnd, bool, str]]:
    """
    Calculate pairwise compatibility for theoretical enzyme ends.
    
    Args:
        ends: List of TheoreticalEnd objects
        include_blunt: Include blunt-blunt pairs
        min_overhang: Minimum overhang length
        require_directional: Only include non-palindromic pairs
        
    Returns:
        List of tuples: (end_a, end_b, directional, reason)
    """
    results = []
    
    for i in range(len(ends)):
        for j in range(i + 1, len(ends)):
            a = ends[i]
            b = ends[j]
            
            # Check compatibility
            compatible, reason = are_theoretical_ends_compatible(
                a, b, include_blunt, min_overhang
            )
            
            if not compatible:
                continue
            
            # Check directionality (non-palindromic)
            directional = not a.is_palindromic
            
            # Skip if directional required but not met
            if require_directional and not directional:
                continue
            
            results.append((a, b, directional, reason))
    
    return results


def format_theoretical_pairs(
    results: List[Tuple[TheoreticalEnd, TheoreticalEnd, bool, str]]
) -> str:
    """
    Format theoretical compatibility results as pairs.
    
    Args:
        results: List of compatibility tuples
        
    Returns:
        Formatted string
    """
    if not results:
        return "No compatible enzyme pairs found.\n"
    
    lines = []
    lines.append("Theoretical sticky-end compatibility (no digest)")
    lines.append("=" * 80)
    lines.append("")
    
    has_n_bases = False
    for a, b, directional, reason in results:
        lines.append(f"{a.enzyme}  | {a.overhang_type} k={a.k} | tpl={a.sticky_template} | palindromic: {'YES' if a.is_palindromic else 'NO'}")
        lines.append(f"{b.enzyme}  | {b.overhang_type} k={b.k} | tpl={b.sticky_template} | palindromic: {'YES' if b.is_palindromic else 'NO'}")
        
        if directional:
            lines.append("  ✔ Compatible (directional)")
        else:
            lines.append("  ✔ Compatible (non-directional)")
        
        lines.append("")
        
        # Check if any templates contain 'N'
        if 'N' in a.sticky_template or 'N' in b.sticky_template:
            has_n_bases = True
    
    # Add note about 'N' bases if present
    if has_n_bases:
        lines.append("")
        lines.append("Note: 'N' in templates indicates bases outside the recognition site")
        lines.append("      (typically from Type IIS enzymes that cut outside their recognition sequence)")
        lines.append("")
    
    lines.append(f"Total compatible pairs: {len(results)}")
    return "\n".join(lines)


def format_theoretical_matrix(
    results: List[Tuple[TheoreticalEnd, TheoreticalEnd, bool, str]],
    all_ends: List[TheoreticalEnd]
) -> str:
    """
    Format theoretical compatibility as a matrix.
    
    Args:
        results: List of compatibility tuples
        all_ends: All enzyme ends for the matrix
        
    Returns:
        Formatted string
    """
    if not all_ends:
        return "No enzymes to display.\n"
    
    lines = []
    lines.append("Theoretical compatibility matrix (no digest)")
    lines.append("=" * 80)
    lines.append("")
    
    # Build compatibility lookup
    compat_set = set()
    directional_set = set()
    
    for a, b, directional, reason in results:
        key1 = (a.enzyme, b.enzyme)
        key2 = (b.enzyme, a.enzyme)
        compat_set.add(key1)
        compat_set.add(key2)
        if directional:
            directional_set.add(key1)
            directional_set.add(key2)
    
    # Create labels
    labels = [e.enzyme for e in all_ends]
    max_label_len = max(len(label) for label in labels)
    col_width = max(4, max_label_len + 1)
    
    # Print header
    header = " " * (col_width + 1) + " ".join(f"{label[:col_width]:>{col_width}}" for label in labels)
    lines.append(header)
    lines.append(" " * (col_width + 1) + "-" * (col_width * len(labels) + len(labels) - 1))
    
    # Print matrix
    for i, end_i in enumerate(all_ends):
        row = f"{labels[i][:col_width]:<{col_width}} "
        for j, end_j in enumerate(all_ends):
            if i == j:
                cell = "."
            else:
                key = (end_i.enzyme, end_j.enzyme)
                if key in compat_set:
                    if key in directional_set:
                        cell = "▶"  # Directional
                    else:
                        cell = "✓"  # Non-directional
                else:
                    cell = "."
            row += f"{cell:^{col_width}} "
        lines.append(row)
    
    lines.append("")
    lines.append("Legend:")
    lines.append("  ✓ = compatible (non-directional/palindromic)")
    lines.append("  ▶ = compatible (directional)")
    lines.append("  . = incompatible or same enzyme")
    lines.append("")
    lines.append(f"Total compatible pairs: {len(results)}")
    
    return "\n".join(lines)


def format_theoretical_detailed(
    results: List[Tuple[TheoreticalEnd, TheoreticalEnd, bool, str]]
) -> str:
    """
    Format theoretical compatibility with detailed information.
    
    Args:
        results: List of compatibility tuples
        
    Returns:
        Formatted string (JSON-like)
    
    Note:
        'N' in templates indicates bases outside the recognition site,
        typically from Type IIS enzymes that cut outside their recognition sequence.
    """
    if not results:
        return "[]"
    
    output = []
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
        output.append(entry)
    
    return json.dumps(output, indent=2)


# ============================================================================
# ENZYME PAIR ANALYSIS (THEORETICAL) - LEGACY
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

