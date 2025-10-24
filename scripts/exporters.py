#!/usr/bin/env python3
"""
Exporters Module
Provides functions to export restriction digest results to GenBank and CSV formats.
"""

import datetime
import csv
import re
from typing import List, Dict
from fragment_calculator import compute_end_metadata


# ============================================================================
# DATA MODELS (for type hints)
# ============================================================================

class Cut:
    """Represents a cut site with enzyme metadata."""
    def __init__(self, pos: int, enzyme: str, recognition_site: str, 
                 cut_index: int, overhang_type: str, overhang_len: int):
        self.pos = pos
        self.enzyme = enzyme
        self.recognition_site = recognition_site
        self.cut_index = cut_index
        self.overhang_type = overhang_type
        self.overhang_len = overhang_len


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def gb_date_today() -> str:
    """
    Get current date formatted for GenBank (DD-MMM-YYYY).
    
    Returns:
        Date string in GenBank format
    """
    return datetime.datetime.now().strftime("%d-%b-%Y").upper()


def sanitize_locus_name(path: str) -> str:
    """
    Sanitize filename to valid LOCUS name (alphanumeric + underscore, max 16 chars).
    
    Args:
        path: File path or name
        
    Returns:
        Sanitized LOCUS name
    """
    # Extract basename without extension
    import os
    basename = os.path.basename(path)
    name = os.path.splitext(basename)[0]
    
    # Replace non-alphanumeric with underscore and truncate
    sanitized = re.sub(r'[^A-Za-z0-9_]', '_', name)
    return sanitized[:16]


def sanitize_genbank_string(s: str) -> str:
    """
    Sanitize string for GenBank format (escape quotes, remove newlines).
    
    Args:
        s: String to sanitize
        
    Returns:
        Sanitized string
    """
    # Replace quotes with single quotes, remove newlines
    return s.replace('"', "'").replace('\n', ' ').replace('\r', ' ')


def wrap_origin(seq: str) -> str:
    """
    Format sequence for GenBank ORIGIN section.
    60 nucleotides per line, grouped in 10-nt blocks, with 1-based index.
    
    Args:
        seq: DNA sequence
        
    Returns:
        Formatted ORIGIN section lines
    """
    out = []
    for i in range(0, len(seq), 60):
        line = seq[i:i+60].lower()
        idx = f"{i+1:>9}"
        # Group into 10-nt blocks
        blocks = " ".join(line[j:j+10] for j in range(0, len(line), 10))
        out.append(f"{idx} {blocks}")
    return "\n".join(out)


def gb_loc_linear(start0: int, end0: int, n: int) -> str:
    """
    Convert 0-based coordinates to GenBank 1-based location string (linear).
    
    Args:
        start0: Start position (0-based, inclusive)
        end0: End position (0-based, exclusive)
        n: Sequence length
        
    Returns:
        GenBank location string
    """
    # 0-based [start, end) -> 1-based [start+1, end]
    return f"{start0+1}..{end0}"


def gb_loc_wrap(start0: int, end0: int, n: int) -> str:
    """
    Convert 0-based wrapping coordinates to GenBank join() location.
    
    Args:
        start0: Start position (0-based, inclusive)
        end0: End position (0-based, exclusive, wraps if < start0)
        n: Sequence length
        
    Returns:
        GenBank location string with join()
    """
    # join(start+1..N, 1..end)
    return f"join({start0+1}..{n},1..{end0})"


def write_feature(f, key: str, location: str, qualifiers: Dict[str, str]):
    """
    Write a feature entry to GenBank file.
    
    Args:
        f: File handle
        key: Feature key (e.g., "source", "misc_feature")
        location: Location string
        qualifiers: Dictionary of qualifier key-value pairs
    """
    # Feature key and location (left-aligned key in 5-char column, then location)
    f.write(f"     {key:<16}{location}\n")
    
    # Qualifiers (21-space indent for /)
    for qkey, qval in qualifiers.items():
        if qval:
            # Escape quotes in value
            safe_val = sanitize_genbank_string(qval)
            f.write(f"                     /{qkey}=\"{safe_val}\"\n")
        else:
            # Boolean qualifiers (no value)
            f.write(f"                     /{qkey}\n")


# ============================================================================
# GENBANK EXPORTER
# ============================================================================

def export_genbank(sequence: str, cuts: List[Dict], fragments: List[Dict], *,
                   path: str, topology: str, definition: str, organism: str) -> None:
    """
    Export restriction digest to GenBank format.
    
    Args:
        sequence: Full DNA sequence
        cuts: List of cut dictionaries with keys:
              - pos (int): 0-based cut position
              - enzyme (str): Enzyme name
              - recognition_site (str): Recognition sequence
              - cut_index (int): Cut position within site
              - overhang_type (str): "5' overhang" | "3' overhang" | "Blunt"
              - overhang_len (int): Length of overhang
        fragments: List of fragment dictionaries from fragment_calculator
        path: Output file path
        topology: "circular" or "linear"
        definition: DEFINITION line text
        organism: SOURCE organism name
    """
    n = len(sequence)
    locus_name = sanitize_locus_name(path)
    date_str = gb_date_today()
    
    # Determine topology string
    topo_str = "circular" if topology == "circular" else "linear"
    
    with open(path, 'w') as f:
        # LOCUS line
        # Format: LOCUS       name    length bp    DNA     topology  date
        f.write(f"LOCUS       {locus_name:<16} {n:>11} bp    DNA     {topo_str:<8} {date_str}\n")
        
        # DEFINITION
        f.write(f"DEFINITION  {sanitize_genbank_string(definition)}\n")
        
        # ACCESSION and VERSION (blank)
        f.write("ACCESSION   \n")
        f.write("VERSION     \n")
        
        # SOURCE
        f.write(f"SOURCE      {sanitize_genbank_string(organism)}\n")
        f.write(f"  ORGANISM  {sanitize_genbank_string(organism)}\n")
        f.write("            synthetic construct.\n")
        
        # FEATURES
        f.write("FEATURES             Location/Qualifiers\n")
        
        # 1. Source feature (entire molecule)
        if topology == "circular":
            source_loc = f"1..{n}"
            source_quals = {
                "mol_type": "other DNA",
                "note": "circular"
            }
        else:
            source_loc = f"1..{n}"
            source_quals = {
                "mol_type": "other DNA"
            }
        
        write_feature(f, "source", source_loc, source_quals)
        
        # 2. Restriction site features (one per cut)
        for cut in cuts:
            pos = cut['pos']
            enzyme = cut['enzyme']
            site = cut.get('recognition_site', '')
            cut_idx = cut.get('cut_index', 0)
            overhang = cut.get('overhang_type', 'Unknown')
            overhang_len = cut.get('overhang_len', 0)
            
            # Try to find the recognition site in sequence for accurate annotation
            # For now, annotate as a single base at the cut position
            if topology == "circular":
                site_loc = f"{pos+1}"
            else:
                site_loc = f"{pos+1}"
            
            # Build note with cut details
            note_parts = []
            if site:
                note_parts.append(f"site={site}")
            note_parts.append(f"cut_index={cut_idx}")
            if overhang:
                note_parts.append(f"overhang={overhang}")
            if overhang_len > 0:
                note_parts.append(f"k={overhang_len}")
            
            note = "; ".join(note_parts)
            
            quals = {
                "label": enzyme,
                "note": note
            }
            
            write_feature(f, "misc_feature", site_loc, quals)
        
        # 3. Fragment features
        for frag in fragments:
            frag_idx = frag['index']
            start = frag['start']
            end = frag['end']
            length = frag['length']
            wraps = frag.get('wraps', False)
            
            # Determine location string
            if wraps and topology == "circular":
                frag_loc = gb_loc_wrap(start, end, n)
            else:
                frag_loc = gb_loc_linear(start, end, n)
            
            # Build boundary info for note
            left_cut = frag['boundaries'].get('left_cut')
            right_cut = frag['boundaries'].get('right_cut')
            
            left_str = "START"
            right_str = "END"
            
            if left_cut and left_cut.get('enzymes'):
                left_enzymes = [e['enzyme'] for e in left_cut['enzymes']]
                left_enz_name = left_enzymes[0] if left_enzymes else "?"
                left_oh = left_cut['enzymes'][0].get('overhang_type', '') if left_cut['enzymes'] else ''
                left_str = f"{left_enz_name}({left_oh})"
            
            if right_cut and right_cut.get('enzymes'):
                right_enzymes = [e['enzyme'] for e in right_cut['enzymes']]
                right_enz_name = right_enzymes[0] if right_enzymes else "?"
                right_oh = right_cut['enzymes'][0].get('overhang_type', '') if right_cut['enzymes'] else ''
                right_str = f"{right_enz_name}({right_oh})"
            
            note = f"length={length}bp; left={left_str}, right={right_str}"
            
            quals = {
                "label": f"fragment_{frag_idx}",
                "note": note
            }
            
            write_feature(f, "misc_feature", frag_loc, quals)
        
        # ORIGIN
        f.write("ORIGIN\n")
        f.write(wrap_origin(sequence))
        f.write("\n//\n")
    
    print(f"✓ GenBank file exported: {path}")


# ============================================================================
# CSV EXPORTERS
# ============================================================================

def export_csv(prefix: str, cuts: List[Dict], fragments: List[Dict], 
               topology: str, dna_sequence: str) -> None:
    """
    Export restriction digest to CSV files (fragments and cuts).
    
    Args:
        prefix: Output file prefix (will create prefix_fragments.csv and prefix_cuts.csv)
        cuts: List of cut dictionaries
        fragments: List of fragment dictionaries
        topology: "circular" or "linear"
        dna_sequence: Full DNA sequence to extract fragment sequences
    """
    # Export fragments CSV
    frag_path = f"{prefix}_fragments.csv"
    with open(frag_path, 'w', newline='') as csvfile:
        fieldnames = [
            'fragment_id', 'start_idx', 'end_idx', 'mode', 'length',
            'left_enzyme', 'left_overhang_type', 'left_overhang_len', 'left_end_bases',
            'right_enzyme', 'right_overhang_type', 'right_overhang_len', 'right_end_bases',
            'sequence'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for frag in fragments:
            frag_idx = frag['index']
            start = frag['start']
            end = frag['end']
            length = frag['length']
            wraps = frag.get('wraps', False)
            
            # Extract sequence
            if wraps and topology == "circular":
                seq = dna_sequence[start:] + dna_sequence[:end]
            else:
                seq = dna_sequence[start:end]
            
            # Get end information using centralized function
            left_cut = frag['boundaries'].get('left_cut')
            right_cut = frag['boundaries'].get('right_cut')
            
            left_enz = ""
            left_oh_type = ""
            left_oh_len = 0
            left_bases = ""
            
            if left_cut and left_cut.get('enzymes'):
                enz_meta = left_cut['enzymes'][0]
                left_enz = enz_meta['enzyme']
                left_oh_type = enz_meta.get('overhang_type', '')
                
                # Use centralized function to compute end metadata
                end_meta = compute_end_metadata(
                    dna=dna_sequence,
                    cut_pos=left_cut['pos'],
                    recognition_site=enz_meta.get('site', ''),
                    cut_index=enz_meta.get('cut_index', 0),
                    overhang_type=left_oh_type,
                    is_left_end=True,
                    circular=(topology == "circular")
                )
                left_oh_len = end_meta['overhang_len']
                left_bases = end_meta['end_bases']
            
            right_enz = ""
            right_oh_type = ""
            right_oh_len = 0
            right_bases = ""
            
            if right_cut and right_cut.get('enzymes'):
                enz_meta = right_cut['enzymes'][0]
                right_enz = enz_meta['enzyme']
                right_oh_type = enz_meta.get('overhang_type', '')
                
                # Use centralized function to compute end metadata
                end_meta = compute_end_metadata(
                    dna=dna_sequence,
                    cut_pos=right_cut['pos'],
                    recognition_site=enz_meta.get('site', ''),
                    cut_index=enz_meta.get('cut_index', 0),
                    overhang_type=right_oh_type,
                    is_left_end=False,
                    circular=(topology == "circular")
                )
                right_oh_len = end_meta['overhang_len']
                right_bases = end_meta['end_bases']
            
            writer.writerow({
                'fragment_id': frag_idx,
                'start_idx': start,
                'end_idx': end,
                'mode': topology,
                'length': length,
                'left_enzyme': left_enz,
                'left_overhang_type': left_oh_type,
                'left_overhang_len': left_oh_len,
                'left_end_bases': left_bases,
                'right_enzyme': right_enz,
                'right_overhang_type': right_oh_type,
                'right_overhang_len': right_oh_len,
                'right_end_bases': right_bases,
                'sequence': seq
            })
    
    print(f"✓ Fragments CSV exported: {frag_path}")
    
    # Export cuts CSV
    cuts_path = f"{prefix}_cuts.csv"
    with open(cuts_path, 'w', newline='') as csvfile:
        fieldnames = [
            'cut_id', 'pos', 'enzyme', 'recognition_site', 
            'cut_index', 'overhang_type', 'overhang_len'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for idx, cut in enumerate(cuts, 1):
            writer.writerow({
                'cut_id': idx,
                'pos': cut['pos'],
                'enzyme': cut['enzyme'],
                'recognition_site': cut.get('recognition_site', ''),
                'cut_index': cut.get('cut_index', 0),
                'overhang_type': cut.get('overhang_type', ''),
                'overhang_len': cut.get('overhang_len', 0)
            })
    
    print(f"✓ Cuts CSV exported: {cuts_path}")

