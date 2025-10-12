#!/usr/bin/env python3
"""
Planner Utilities Module
Helper functions for the cloning planner including spec loading, validation, and formatting.
"""

import json
import os
from typing import Dict, Any, Optional


# ============================================================================
# SPEC LOADING
# ============================================================================

def load_json_or_yaml(path: str) -> Dict[str, Any]:
    """
    Load a specification from JSON or YAML file.
    
    Args:
        path: Path to spec file (.json or .yaml/.yml)
        
    Returns:
        Specification dictionary
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Spec file not found: {path}")
    
    # Determine format from extension
    _, ext = os.path.splitext(path.lower())
    
    if ext == '.json':
        return load_json(path)
    elif ext in ['.yaml', '.yml']:
        return load_yaml(path)
    else:
        # Try JSON first, then YAML
        try:
            return load_json(path)
        except json.JSONDecodeError:
            try:
                return load_yaml(path)
            except Exception:
                raise ValueError(f"Could not parse file as JSON or YAML: {path}")


def load_json(path: str) -> Dict[str, Any]:
    """
    Load JSON specification file.
    
    Args:
        path: Path to JSON file
        
    Returns:
        Specification dictionary
    """
    with open(path, 'r') as f:
        try:
            spec = json.load(f)
            return spec
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {path}: {e}")


def load_yaml(path: str) -> Dict[str, Any]:
    """
    Load YAML specification file.
    
    Args:
        path: Path to YAML file
        
    Returns:
        Specification dictionary
        
    Raises:
        ImportError: If PyYAML is not installed
    """
    try:
        import yaml
    except ImportError:
        raise ImportError(
            "PyYAML is required for YAML support. Install with: pip install pyyaml"
        )
    
    with open(path, 'r') as f:
        try:
            spec = yaml.safe_load(f)
            return spec
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in {path}: {e}")


# ============================================================================
# SPEC VALIDATION
# ============================================================================

def validate_spec(spec: Dict[str, Any]) -> tuple[bool, str]:
    """
    Validate a cloning specification.
    
    Args:
        spec: Specification dictionary
        
    Returns:
        Tuple of (valid, error_message)
    """
    # Check required top-level keys
    required_keys = ['vector', 'target']
    for key in required_keys:
        if key not in spec:
            return False, f"Missing required key: {key}"
    
    # Validate vector
    vector = spec['vector']
    if not isinstance(vector, dict):
        return False, "vector must be a dictionary"
    
    if 'name' not in vector:
        return False, "vector must have a 'name' field"
    
    if 'fasta' not in vector:
        return False, "vector must have a 'fasta' field"
    
    # Validate inserts
    inserts = spec.get('inserts', [])
    if not isinstance(inserts, list):
        return False, "inserts must be a list"
    
    for i, insert in enumerate(inserts):
        if not isinstance(insert, dict):
            return False, f"insert {i} must be a dictionary"
        
        if 'name' not in insert:
            return False, f"insert {i} must have a 'name' field"
        
        if 'fasta' not in insert:
            return False, f"insert {i} must have a 'fasta' field"
    
    # Validate target
    target = spec['target']
    if not isinstance(target, dict):
        return False, "target must be a dictionary"
    
    if 'order' not in target:
        return False, "target must have an 'order' field"
    
    order = target['order']
    if not isinstance(order, list):
        return False, "target.order must be a list"
    
    # Check that all parts in order are defined
    vector_name = vector['name']
    insert_names = {ins['name'] for ins in inserts}
    all_parts = {vector_name} | insert_names
    
    for part_name in order:
        if part_name not in all_parts:
            return False, f"Unknown part in order: {part_name}"
    
    # Validate junctions if present
    junctions = target.get('junctions', [])
    if not isinstance(junctions, list):
        return False, "target.junctions must be a list"
    
    for i, junction in enumerate(junctions):
        if not isinstance(junction, dict):
            return False, f"junction {i} must be a dictionary"
        
        required_junction_keys = ['left', 'right']
        for key in required_junction_keys:
            if key not in junction:
                return False, f"junction {i} missing required key: {key}"
    
    # Validate constraints
    constraints = spec.get('constraints', {})
    if not isinstance(constraints, dict):
        return False, "constraints must be a dictionary"
    
    return True, ""


# ============================================================================
# SPEC UTILITIES
# ============================================================================

def get_part_sequence(spec: Dict[str, Any], part_name: str) -> Optional[str]:
    """
    Get the sequence of a named part from the spec.
    
    Args:
        spec: Specification dictionary
        part_name: Name of the part (vector or insert)
        
    Returns:
        DNA sequence string, or None if not found
    """
    from sim import read_dna_sequence
    
    # Check vector
    vector = spec.get('vector', {})
    if vector.get('name') == part_name:
        try:
            return read_dna_sequence(vector['fasta'])
        except (FileNotFoundError, ValueError):
            return None
    
    # Check inserts
    for insert in spec.get('inserts', []):
        if insert.get('name') == part_name:
            try:
                return read_dna_sequence(insert['fasta'])
            except (FileNotFoundError, ValueError):
                return None
    
    return None


def get_junction_spec(spec: Dict[str, Any], left_part: str, right_part: str) -> Optional[Dict]:
    """
    Get junction specification between two parts.
    
    Args:
        spec: Specification dictionary
        left_part: Name of left part
        right_part: Name of right part
        
    Returns:
        Junction specification dict, or None if not found
    """
    junctions = spec.get('target', {}).get('junctions', [])
    
    for junction in junctions:
        if junction.get('left') == left_part and junction.get('right') == right_part:
            return junction
    
    return None


# ============================================================================
# GEL PREDICTION
# ============================================================================

def predict_gel_bands(fragment_lengths: list[int]) -> str:
    """
    Predict gel bands for a set of fragments.
    
    Uses the gel simulation from fragment_calculator but returns a simple
    text representation.
    
    Args:
        fragment_lengths: List of fragment sizes in bp
        
    Returns:
        String representation of predicted bands
    """
    if not fragment_lengths:
        return "No fragments"
    
    # Sort by size (largest first for gel representation)
    sorted_frags = sorted(fragment_lengths, reverse=True)
    
    # Create simple text representation
    lines = []
    lines.append("Predicted bands:")
    for i, size in enumerate(sorted_frags, 1):
        # Categorize size
        if size >= 10000:
            category = "very large"
        elif size >= 5000:
            category = "large"
        elif size >= 1000:
            category = "medium"
        elif size >= 500:
            category = "small"
        else:
            category = "very small"
        
        lines.append(f"  {i}. {size:>6} bp  ({category})")
    
    return '\n'.join(lines)


# ============================================================================
# READING FRAME UTILITIES
# ============================================================================

def check_frame_preservation(left_seq: str, right_seq: str, scar: str, 
                             left_frame: int = 0, right_frame: int = 0) -> tuple[bool, str]:
    """
    Check if a junction preserves reading frame and provides details.
    
    Note: For accurate stop codon detection, context sequences should ideally be at least
    9 bases. If shorter sequences are provided, analysis will use available sequence.
    
    Args:
        left_seq: Last 12 bases of left sequence (minimum: 3 bases recommended)
        right_seq: First 12 bases of right sequence (minimum: 3 bases recommended)
        scar: Scar sequence at junction
        left_frame: Reading frame of left sequence (0, 1, or 2)
        right_frame: Reading frame of right sequence (0, 1, or 2)
        
    Returns:
        Tuple of (frame_ok, reason)
    """
    # Validate minimum context length (soft requirement)
    min_context_len = 3  # Minimum needed for one codon
    if len(left_seq) < min_context_len or len(right_seq) < min_context_len:
        return True, f"Insufficient context for frame analysis (left={len(left_seq)}bp, right={len(right_seq)}bp, need >={min_context_len}bp)"
    
    # Calculate frame offset from scar
    scar_len = len(scar)
    
    # Frame preservation requires: (left_frame + scar_len) % 3 == right_frame
    expected_right_frame = (left_frame + scar_len) % 3
    
    if expected_right_frame != right_frame:
        return False, f"Frame shift: scar adds {scar_len} bp, shifting frame from {left_frame} to {expected_right_frame}, but right expects {right_frame}"
    
    # Check for stop codons in junction (only check IN-FRAME codons)
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    # Build junction context - use available sequence length (up to 9 bases)
    context_len = min(9, len(left_seq), len(right_seq))
    context_left = left_seq[-context_len:] if len(left_seq) >= context_len else left_seq
    context_right = right_seq[:context_len] if len(right_seq) >= context_len else right_seq
    junction_context = context_left + scar + context_right
    
    # Determine the reading frame offset for the start of junction_context
    # We assume context_left is from a sequence with reading frame left_frame
    # and has been trimmed to align with codon boundaries
    # The start of junction_context is at frame offset based on context_left length
    frame_offset = (left_frame) % 3
    
    # Scan for stop codons (only at in-frame positions)
    stops_found = []
    # Start at frame_offset and check every 3rd position (in-frame codons)
    for i in range(frame_offset, len(junction_context) - 2, 3):
        codon = junction_context[i:i+3].upper()
        if len(codon) == 3 and codon in stop_codons:
            stops_found.append((i, codon))
    
    if stops_found:
        stop_positions = ', '.join(f"{codon}@{pos}" for pos, codon in stops_found)
        return False, f"Stop codons found in junction: {stop_positions}"
    
    return True, "Frame preserved, no stop codons"


# ============================================================================
# OVERHANG UTILITIES
# ============================================================================

def validate_golden_gate_overhangs(overhangs: list[str]) -> tuple[bool, str]:
    """
    Validate a set of Golden Gate overhangs.
    
    Checks for:
    - Duplicate overhangs (causes unwanted ligations)
    - Palindromic overhangs (non-directional)
    - Self-complementary pairs
    
    Args:
        overhangs: List of 4-nt overhang sequences
        
    Returns:
        Tuple of (valid, error_message)
    """
    from ligation_compatibility import revcomp
    
    # Check for duplicates
    if len(overhangs) != len(set(overhangs)):
        duplicates = [oh for oh in overhangs if overhangs.count(oh) > 1]
        return False, f"Duplicate overhangs: {', '.join(set(duplicates))}"
    
    # Check for palindromes
    palindromes = []
    for oh in overhangs:
        if oh == revcomp(oh):
            palindromes.append(oh)
    
    if palindromes:
        return False, f"Palindromic overhangs (non-directional): {', '.join(palindromes)}"
    
    # Check for self-complementary pairs
    complements = []
    for oh in overhangs:
        oh_rc = revcomp(oh)
        if oh_rc in overhangs and oh != oh_rc:
            pair = tuple(sorted([oh, oh_rc]))
            if pair not in complements:
                complements.append(pair)
    
    if complements:
        pairs_str = ', '.join(f"{a}↔{b}" for a, b in complements)
        return False, f"Complementary overhang pairs (causes unwanted ligation): {pairs_str}"
    
    return True, "Overhangs are valid"


def design_golden_gate_overhangs(num_junctions: int, avoid_palindromes: bool = True) -> list[str]:
    """
    Design a set of Golden Gate overhangs for multi-part assembly.
    
    Args:
        num_junctions: Number of junctions (need num_junctions + 1 overhangs)
        avoid_palindromes: If True, avoid palindromic overhangs
        
    Returns:
        List of 4-nt overhang sequences
    """
    # Common non-palindromic overhangs used in Golden Gate
    # (from NEB Golden Gate Assembly Protocol)
    common_overhangs = [
        'AATG', 'AAGC', 'ACAG', 'ACAT',
        'AGAT', 'AGGT', 'ATAG', 'ATCC',
        'CAAT', 'CAAG', 'CCAG', 'CGAT',
        'CTAG', 'CTCC', 'GAAT', 'GATG',
        'GCAG', 'GGAT', 'GTAG', 'TAAT'
    ]
    
    # Filter out palindromes if requested
    if avoid_palindromes:
        from ligation_compatibility import revcomp
        common_overhangs = [oh for oh in common_overhangs if oh != revcomp(oh)]
    
    # Return the requested number
    num_needed = num_junctions + 1
    if num_needed > len(common_overhangs):
        raise ValueError(f"Cannot design {num_needed} overhangs with current library")
    
    return common_overhangs[:num_needed]


# ============================================================================
# PLAN FORMATTING
# ============================================================================

def format_plan_detailed(plan, export_dir: Optional[str] = None) -> str:
    """
    Format a plan with detailed information suitable for protocol generation.
    
    Args:
        plan: Plan object
        export_dir: Optional directory for exported files
        
    Returns:
        Formatted detailed plan string
    """
    lines = []
    
    lines.append("=" * 80)
    lines.append("DETAILED CLONING PROTOCOL")
    lines.append("=" * 80)
    lines.append("")
    
    if not plan.feasible:
        lines.append(f"⚠ PLAN NOT FEASIBLE: {plan.reason}")
        return '\n'.join(lines)
    
    lines.append(f"Total steps: {len(plan.steps)}")
    lines.append(f"Complexity score: {plan.score:.2f}")
    lines.append("")
    
    for step_num, step in enumerate(plan.steps, 1):
        lines.append("=" * 80)
        lines.append(f"STEP {step_num}: {step.action.upper()}")
        lines.append("=" * 80)
        lines.append(f"Name: {step.name}")
        lines.append("")
        
        # Inputs
        lines.append("Inputs:")
        for inp in step.inputs:
            lines.append(f"  - {inp}")
        lines.append("")
        
        # Action-specific details
        if step.action == "digest":
            enzymes = step.params.get('enzymes', [])
            deP = step.params.get('dephosphorylate', False)
            
            lines.append("Protocol:")
            lines.append(f"  1. Set up restriction digest:")
            for enz in enzymes:
                lines.append(f"     - Add {enz}")
            lines.append(f"  2. Incubate at 37°C for 1-2 hours")
            if deP:
                lines.append(f"  3. Add Antarctic Phosphatase (dephosphorylation)")
                lines.append(f"  4. Incubate at 37°C for 30 minutes")
            lines.append(f"  {4 if deP else 3}. Purify by gel extraction or column")
            lines.append("")
            
            # Predicted fragments
            if step.predicted_fragments:
                lines.append("Expected fragments:")
                for i, frag in enumerate(step.predicted_fragments):
                    lines.append(f"  Fragment {i+1}: {frag['length']} bp")
                lines.append("")
        
        elif step.action == "ligate":
            directional = step.params.get('directional', False)
            scar = step.params.get('scar', '')
            
            lines.append("Protocol:")
            lines.append(f"  1. Set up ligation reaction:")
            lines.append(f"     - Add vector (50 ng)")
            lines.append(f"     - Add insert (3:1 molar ratio)")
            lines.append(f"     - Add T4 DNA Ligase")
            lines.append(f"  2. Incubate at 16°C overnight (or room temp 1 hour)")
            lines.append(f"  3. Transform into competent cells")
            lines.append("")
            
            lines.append(f"Directionality: {'Yes ✓' if directional else 'No (screen colonies)'}")
            if scar:
                lines.append(f"Junction scar: {scar}")
            lines.append("")
        
        elif step.action == "GG":
            enzyme = step.params.get('enzyme', 'BsaI')
            overhangs = step.params.get('overhangs', [])
            
            lines.append("Protocol (Golden Gate):")
            lines.append(f"  1. Set up one-pot reaction:")
            lines.append(f"     - Add all parts (equimolar, 50 ng each)")
            lines.append(f"     - Add {enzyme}")
            lines.append(f"     - Add T4 DNA Ligase")
            lines.append(f"  2. Run thermocycler protocol:")
            lines.append(f"     - 26 cycles: [37°C 2 min, 16°C 5 min]")
            lines.append(f"     - Final: 50°C 5 min, 80°C 10 min")
            lines.append(f"  3. Transform into competent cells")
            lines.append("")
            
            if overhangs:
                lines.append(f"Designed overhangs:")
                for i, oh in enumerate(overhangs):
                    lines.append(f"  Junction {i+1}: {oh}")
                lines.append("")
        
        elif step.action == "PCR":
            primers = step.params.get('primers', {})
            
            lines.append("Protocol:")
            lines.append(f"  1. Set up PCR reaction:")
            lines.append(f"     - Template DNA")
            lines.append(f"     - Forward primer (see below)")
            lines.append(f"     - Reverse primer (see below)")
            lines.append(f"     - High-fidelity polymerase")
            lines.append(f"  2. Run PCR with appropriate conditions")
            lines.append(f"  3. Purify PCR product")
            lines.append("")
            
            if primers:
                lines.append("Primers:")
                lines.append(f"  Forward: {primers.get('forward', 'N/A')}")
                lines.append(f"  Reverse: {primers.get('reverse', 'N/A')}")
                lines.append("")
        
        # Outputs
        lines.append("Expected outputs:")
        for out in step.outputs:
            lines.append(f"  - {out}")
        lines.append("")
        
        # Export information
        if export_dir:
            lines.append(f"Export files: {export_dir}/step_{step_num:02d}_*")
            lines.append("")
    
    # Final construct
    lines.append("=" * 80)
    lines.append("FINAL CONSTRUCT")
    lines.append("=" * 80)
    
    if plan.final:
        lines.append(f"Name: {plan.final.name}")
        lines.append(f"Length: {len(plan.final.seq)} bp")
        lines.append(f"Circular: {plan.final.circular}")
        lines.append("")
    
    lines.append("=" * 80)
    
    return '\n'.join(lines)


def format_plan_json(plan) -> str:
    """
    Format plan as JSON for programmatic use.
    
    Args:
        plan: Plan object
        
    Returns:
        JSON string
    """
    plan_dict = {
        'feasible': plan.feasible,
        'score': plan.score,
        'steps': [],
        'final': None
    }
    
    if not plan.feasible:
        plan_dict['reason'] = plan.reason
        return json.dumps(plan_dict, indent=2)
    
    for step in plan.steps:
        step_dict = {
            'name': step.name,
            'action': step.action,
            'inputs': step.inputs,
            'outputs': step.outputs,
            'params': step.params,
            'predicted_fragments': step.predicted_fragments,
            'cost': step.cost
        }
        plan_dict['steps'].append(step_dict)
    
    if plan.final:
        plan_dict['final'] = {
            'name': plan.final.name,
            'length': len(plan.final.seq),
            'circular': plan.final.circular,
            'notes': plan.final.notes
        }
    
    return json.dumps(plan_dict, indent=2)

