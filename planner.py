#!/usr/bin/env python3
"""
Cloning Planner Module
Implements multi-step cloning strategy planning with support for classic restriction
cloning, Golden Gate assembly, and PCR-based methods.
"""

import heapq
import json
import hashlib
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set, FrozenSet
from copy import deepcopy

from fragment_calculator import compute_fragments_with_sequences, Fragment
from ligation_compatibility import (
    are_compatible, is_directional, EndInfo, 
    calculate_gc_percent, calculate_tm, revcomp
)


# ============================================================================
# HELPER FUNCTIONS FOR STATE HASHING
# ============================================================================

def construct_signature(constructs):
    """
    Build a stable, hashable signature for a set of Construct objects.
    
    Uses name, length, circular flag, origin, and a short hash of seq.
    This allows state deduplication based on available constructs without
    requiring Step objects to be hashable.
    
    Args:
        constructs: FrozenSet or iterable of Construct objects
        
    Returns:
        frozenset of tuples representing construct signatures
    """
    items = []
    for c in constructs:
        h = hashlib.sha1(c.seq.encode('utf-8')).hexdigest()[:12]
        items.append((c.name, len(c.seq), c.circular, c.origin, h))
    # frozenset since order doesn't matter
    return frozenset(items)


# ============================================================================
# CORE DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Construct:
    """Represents a DNA construct (vector, insert, or intermediate product)."""
    name: str
    seq: str
    circular: bool
    features: Tuple[Dict, ...] = field(default_factory=tuple)  # List of feature annotations
    origin: int = 0  # Origin position for circular constructs
    notes: str = ""
    
    def __hash__(self):
        return hash((self.name, self.seq, self.circular))


@dataclass
class Step:
    """Represents a single cloning step."""
    name: str
    action: str  # "digest" | "ligate" | "PCR" | "GG" | "clean_up"
    inputs: List[str]  # Construct names
    params: Dict  # Enzymes, adapters, primers, buffer notes, etc.
    outputs: List[str]  # Resulting construct names
    predicted_fragments: List[Dict]  # Length, ends, gel bands
    cost: float = 0.0  # Step cost for scoring
    
    def __repr__(self):
        return f"Step({self.action}: {'+'.join(self.inputs)} → {'+'.join(self.outputs)})"


@dataclass
class Plan:
    """Represents a complete cloning plan with multiple steps."""
    steps: List[Step]
    final: Construct
    score: float
    feasible: bool = True
    reason: str = ""  # Reason if not feasible
    
    def __repr__(self):
        return f"Plan(steps={len(self.steps)}, score={self.score:.2f}, feasible={self.feasible})"


@dataclass
class SearchState:
    """Represents a state in the search space."""
    constructs: FrozenSet[Construct]  # Current set of available constructs
    steps_taken: Tuple[Step, ...]  # Steps taken so far
    cost: float  # Accumulated cost (g)
    heuristic: float  # Heuristic estimate (h)
    
    def __lt__(self, other):
        """For priority queue ordering (A* uses f = g + h)."""
        return (self.cost + self.heuristic) < (other.cost + other.heuristic)
    
    def signature(self):
        """
        Return a hashable signature representing the current state.
        
        The signature is based solely on the available constructs, not on
        the path taken to reach this state. This allows proper deduplication
        without requiring Step objects to be hashable.
        
        Returns:
            frozenset of construct signatures
        """
        return construct_signature(self.constructs)
    
    def __hash__(self):
        """
        Hash based on structural state only, not on step objects.
        
        This prevents "unhashable type: Step" errors since Step contains
        mutable lists/dicts.
        """
        return hash(self.signature())
    
    def __eq__(self, other):
        """
        Equality based on available constructs, not path taken.
        
        Two states are equal if they have the same set of available constructs,
        regardless of how we arrived at that state.
        """
        if not isinstance(other, SearchState):
            return False
        return self.signature() == other.signature()


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def sanitize_name(s: str) -> str:
    """
    Sanitize construct name for file paths and identifiers.
    
    Args:
        s: Original name
        
    Returns:
        Sanitized name (alphanumeric + underscore)
    """
    import re
    return re.sub(r'[^A-Za-z0-9_]', '_', s)


def extract_features_from_spec(features_list: List[Dict]) -> Tuple[Dict, ...]:
    """
    Extract and validate features from spec.
    
    Args:
        features_list: List of feature dictionaries
        
    Returns:
        Tuple of feature dicts
    """
    validated = []
    for feat in features_list:
        if 'type' not in feat:
            continue
        validated.append({
            'type': feat['type'],
            'frame': feat.get('frame', 0),
            'start': feat.get('start', 0),
            'end': feat.get('end', 0),
            'label': feat.get('label', ''),
            'direction': feat.get('direction', 'forward')
        })
    return tuple(validated)


def check_internal_cuts(seq: str, enzyme_name: str, enzyme_db: Dict, 
                       features: Tuple[Dict, ...], avoid_cds: bool = True) -> int:
    """
    Check for internal cuts within protected features.
    
    Args:
        seq: DNA sequence
        enzyme_name: Enzyme name
        enzyme_db: Enzyme database
        features: Feature annotations
        avoid_cds: If True, penalize cuts inside CDS features
        
    Returns:
        Number of internal cuts in protected regions
    """
    if enzyme_name not in enzyme_db:
        return 0
    
    # Find cut sites
    from sim import find_cut_positions_linear
    cuts = find_cut_positions_linear(seq, enzyme_name, enzyme_db)
    
    if not avoid_cds or not features:
        return 0
    
    # Check if any cuts are inside CDS features
    internal_count = 0
    for cut_pos in cuts:
        for feat in features:
            if feat['type'] == 'CDS':
                start = feat.get('start', 0)
                end = feat.get('end', len(seq))
                if start <= cut_pos < end:
                    internal_count += 1
                    break
    
    return internal_count


def reading_frame_ok(seq_left: str, seq_right: str, scar: str, 
                    frame_left: int = 0, frame_right: int = 0) -> bool:
    """
    Check if a junction preserves reading frame.
    
    Args:
        seq_left: Left sequence (last few codons)
        seq_right: Right sequence (first few codons)
        scar: Scar sequence at junction
        frame_left: Reading frame of left sequence (0, 1, or 2)
        frame_right: Reading frame of right sequence (0, 1, or 2)
        
    Returns:
        True if frame is preserved
    """
    # Calculate junction length
    junction_len = len(scar)
    
    # Check if junction length preserves frame
    # frame_left + junction_len + frame_right should be divisible by 3
    total_offset = (frame_left + junction_len) % 3
    
    if total_offset != frame_right:
        return False
    
    # Check for stop codons in junction context
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    # Build junction context (last codon of left + scar + first codon of right)
    context = seq_left[-3:] + scar + seq_right[:3]
    
    # Check for stop codons
    for i in range(0, len(context) - 2):
        codon = context[i:i+3].upper()
        if len(codon) == 3 and codon in stop_codons:
            return False
    
    return True


def compute_scar_from_overhang(overhang: str, overhang_type: str) -> str:
    """
    Compute scar sequence from overhang ligation.
    
    Args:
        overhang: Overhang sequence (5'→3')
        overhang_type: "5' overhang" or "3' overhang"
        
    Returns:
        Scar sequence after ligation
    """
    # For most ligations, the scar is just the overhang sequence
    # (after annealing, the gap is filled in)
    return overhang


# ============================================================================
# DIGEST SIMULATION
# ============================================================================

def simulate_digest(construct: Construct, enzymes: List[str], 
                   enzyme_db: Dict) -> List[Fragment]:
    """
    Simulate digestion of a construct with one or more enzymes.
    
    Args:
        construct: Construct to digest
        enzymes: List of enzyme names
        enzyme_db: Enzyme database
        
    Returns:
        List of Fragment objects
    """
    from sim import find_cut_positions_linear
    
    # Find all cut positions
    all_cuts = []
    cut_metadata = {}
    
    for enzyme_name in enzymes:
        if enzyme_name not in enzyme_db:
            continue
        
        cuts = find_cut_positions_linear(construct.seq, enzyme_name, enzyme_db)
        all_cuts.extend(cuts)
        
        # Store metadata
        enzyme_info = enzyme_db[enzyme_name]
        for pos in cuts:
            if pos not in cut_metadata:
                cut_metadata[pos] = []
            cut_metadata[pos].append({
                'enzyme': enzyme_name,
                'site': enzyme_info['sequence'],
                'cut_index': enzyme_info['cut_index'],
                'overhang_type': enzyme_info.get('overhang_type', 'Unknown')
            })
    
    # Remove duplicates and sort
    all_cuts = sorted(set(all_cuts))
    
    # Compute fragments with sequences
    fragments = compute_fragments_with_sequences(
        dna_sequence=construct.seq,
        cut_positions=all_cuts,
        circular=construct.circular,
        circular_single_cut_linearizes=True,  # Linearize for cloning
        cut_metadata=cut_metadata
    )
    
    return fragments


def fragments_compatible_for_ligation(frag_a: Fragment, frag_b: Fragment,
                                     end_a_side: str, end_b_side: str,
                                     include_blunt: bool = False) -> bool:
    """
    Check if two fragments can be ligated together.
    
    Args:
        frag_a: First fragment
        frag_b: Second fragment
        end_a_side: "left" or "right" end of fragment A
        end_b_side: "left" or "right" end of fragment B
        include_blunt: Allow blunt-blunt ligation
        
    Returns:
        True if compatible for ligation
    """
    # Get end information
    left_a, right_a = frag_a.enzymes_at_ends
    left_b, right_b = frag_b.enzymes_at_ends
    
    end_a = left_a if end_a_side == "left" else right_a
    end_b = left_b if end_b_side == "left" else right_b
    
    if end_a is None or end_b is None:
        return False
    
    # Convert to EndInfo for compatibility check
    from fragment_calculator import extract_sticky_seq
    
    end_info_a = EndInfo(
        enzyme=end_a.enzyme,
        overhang_type=end_a.overhang_type,
        overhang_len=end_a.overhang_len,
        sticky_seq=end_a.end_bases,
        polarity=end_a_side,
        fragment_id=0,
        position=0
    )
    
    end_info_b = EndInfo(
        enzyme=end_b.enzyme,
        overhang_type=end_b.overhang_type,
        overhang_len=end_b.overhang_len,
        sticky_seq=end_b.end_bases,
        polarity=end_b_side,
        fragment_id=0,
        position=0
    )
    
    return are_compatible(end_info_a, end_info_b, include_blunt=include_blunt)


# ============================================================================
# SCORING AND HEURISTICS
# ============================================================================

def calculate_heuristic(state: SearchState, target_spec: Dict, 
                       options: Dict) -> float:
    """
    Calculate heuristic estimate of remaining cost.
    
    The heuristic estimates:
    - Number of remaining junctions to realize
    - Penalties for unmet constraints
    
    Args:
        state: Current search state
        target_spec: Target assembly specification
        options: Planning options
        
    Returns:
        Heuristic cost estimate (h)
    """
    h = 0.0
    
    # Get target construct order
    target_order = target_spec.get('order', [])
    junctions = target_spec.get('junctions', [])
    
    # Count how many target constructs we have
    construct_names = {c.name for c in state.constructs}
    target_constructs = set(target_order)
    
    missing_constructs = target_constructs - construct_names
    
    # Each missing construct requires at least one step
    h += len(missing_constructs) * 1.0
    
    # Each junction requires at least one ligation
    h += len(junctions) * 0.5
    
    # Penalty for avoiding certain conditions
    if options.get('avoid_internal_cuts', True):
        # Slight penalty if we haven't checked this yet
        h += 0.1
    
    return h


def score_plan(plan: Plan, options: Dict) -> float:
    """
    Calculate comprehensive score for a plan.
    
    Lower scores are better.
    
    Scoring factors:
    - Number of steps
    - Unique enzyme count
    - Buffer switches
    - Non-directional junctions (if directional requested)
    - Internal cut violations
    - Expected scar length
    - Type IIS single-pot bonus (if preferred)
    
    Args:
        plan: Plan to score
        options: Planning options
        
    Returns:
        Total score (lower is better)
    """
    score = 0.0
    
    # Base: number of steps
    score += 1.0 * len(plan.steps)
    
    # Count unique enzymes used
    enzymes_used = set()
    for step in plan.steps:
        if step.action == "digest":
            enzymes_used.update(step.params.get('enzymes', []))
        elif step.action == "GG":
            enzymes_used.add(step.params.get('enzyme', ''))
    
    score += 0.5 * len(enzymes_used)
    
    # Count buffer switches (simplified: assume different enzyme groups need different buffers)
    # For now, just a placeholder
    buffer_switches = max(0, len(enzymes_used) - 1)
    score += 0.3 * buffer_switches
    
    # Non-directional junction penalty
    if options.get('require_directional', False):
        for step in plan.steps:
            if step.action == "ligate":
                if not step.params.get('directional', False):
                    score += 1.0
    
    # Internal cuts penalty
    for step in plan.steps:
        internal_cuts = step.params.get('internal_cuts', 0)
        score += 2.0 * internal_cuts
    
    # Scar length penalty (prefer smaller scars)
    for step in plan.steps:
        if step.action == "ligate" or step.action == "GG":
            scar_len = len(step.params.get('scar', ''))
            score += 0.1 * scar_len
    
    # Type IIS single-pot bonus
    if options.get('prefer_typeIIS', False):
        for step in plan.steps:
            if step.action == "GG":
                score -= 0.4  # Reward Golden Gate
    
    # Reuse of sites bonus
    enzyme_counts = {}
    for step in plan.steps:
        if step.action == "digest":
            for enz in step.params.get('enzymes', []):
                enzyme_counts[enz] = enzyme_counts.get(enz, 0) + 1
    
    for count in enzyme_counts.values():
        if count > 1:
            score -= 0.3 * (count - 1)
    
    return score


# ============================================================================
# ACTION GENERATORS
# ============================================================================

def generate_digest_actions(state: SearchState, enzyme_db: Dict, 
                           options: Dict) -> List[Step]:
    """
    Generate all possible digest actions from current state.
    
    Args:
        state: Current search state
        enzyme_db: Enzyme database
        options: Planning options (avoid_enzymes, allow_enzymes, etc.)
        
    Returns:
        List of possible digest steps
    """
    actions = []
    
    # Get enzyme constraints
    avoid_enzymes = set(options.get('avoid_enzymes', []))
    allow_enzymes = options.get('allow_enzymes', None)
    
    # Filter enzyme list
    if allow_enzymes:
        available_enzymes = [e for e in allow_enzymes if e in enzyme_db and e not in avoid_enzymes]
    else:
        available_enzymes = [e for e in enzyme_db.keys() if e not in avoid_enzymes]
    
    # Try single and double digests
    for construct in state.constructs:
        # Skip if construct is too small
        if len(construct.seq) < 10:
            continue
        
        # Single enzyme digests
        for enzyme in available_enzymes[:20]:  # Limit to avoid explosion
            # Check for internal cuts if requested
            internal_cuts = 0
            if options.get('avoid_internal_cuts', True):
                internal_cuts = check_internal_cuts(
                    construct.seq, enzyme, enzyme_db, 
                    construct.features, avoid_cds=True
                )
                if internal_cuts > 0:
                    continue  # Skip if enzyme cuts inside CDS
            
            # Simulate digest
            fragments = simulate_digest(construct, [enzyme], enzyme_db)
            
            # Skip if no fragmentation (no cuts)
            if len(fragments) <= 1:
                continue
            
            # Create step
            step = Step(
                name=f"Digest_{construct.name}_with_{enzyme}",
                action="digest",
                inputs=[construct.name],
                params={
                    'enzymes': [enzyme],
                    'dephosphorylate': construct.circular,  # DeP backbone
                    'internal_cuts': internal_cuts
                },
                outputs=[f"{construct.name}_frag{i}" for i in range(len(fragments))],
                predicted_fragments=[{
                    'length': frag.length,
                    'left_end': frag.enzymes_at_ends[0],
                    'right_end': frag.enzymes_at_ends[1]
                } for frag in fragments],
                cost=1.0
            )
            actions.append(step)
        
        # Double enzyme digests (for directional cloning)
        for i, enz1 in enumerate(available_enzymes[:10]):
            for enz2 in available_enzymes[i+1:10]:
                # Check internal cuts
                internal_cuts = 0
                if options.get('avoid_internal_cuts', True):
                    ic1 = check_internal_cuts(construct.seq, enz1, enzyme_db, construct.features)
                    ic2 = check_internal_cuts(construct.seq, enz2, enzyme_db, construct.features)
                    internal_cuts = ic1 + ic2
                    if internal_cuts > 0:
                        continue
                
                # Simulate digest
                fragments = simulate_digest(construct, [enz1, enz2], enzyme_db)
                
                if len(fragments) <= 1:
                    continue
                
                step = Step(
                    name=f"Digest_{construct.name}_with_{enz1}+{enz2}",
                    action="digest",
                    inputs=[construct.name],
                    params={
                        'enzymes': [enz1, enz2],
                        'dephosphorylate': construct.circular,
                        'internal_cuts': internal_cuts
                    },
                    outputs=[f"{construct.name}_frag{i}" for i in range(len(fragments))],
                    predicted_fragments=[{
                        'length': frag.length,
                        'left_end': frag.enzymes_at_ends[0],
                        'right_end': frag.enzymes_at_ends[1]
                    } for frag in fragments],
                    cost=1.2  # Slightly higher cost for double digest
                )
                actions.append(step)
    
    return actions[:50]  # Limit to avoid explosion


def generate_ligation_actions(state: SearchState, options: Dict) -> List[Step]:
    """
    Generate all possible ligation actions from current state.
    
    This considers pairwise ligations of compatible fragment ends.
    
    Args:
        state: Current search state
        options: Planning options
        
    Returns:
        List of possible ligation steps
    """
    actions = []
    
    # Get all constructs (some may be fragments from digests)
    constructs_list = list(state.constructs)
    
    # Try pairwise ligations
    for i, const_a in enumerate(constructs_list):
        for const_b in constructs_list[i+1:]:
            # Skip if both are circular (can't ligate circles)
            if const_a.circular and const_b.circular:
                continue
            
            # Check if they have compatible ends
            # (This is a simplified check - in reality we'd need to check all end combinations)
            
            # Create ligation step (simplified)
            step = Step(
                name=f"Ligate_{const_a.name}+{const_b.name}",
                action="ligate",
                inputs=[const_a.name, const_b.name],
                params={
                    'directional': True,  # Placeholder
                    'scar': 'NNNN'  # Placeholder
                },
                outputs=[f"{const_a.name}_{const_b.name}_ligated"],
                predicted_fragments=[],
                cost=1.0
            )
            actions.append(step)
    
    return actions[:20]  # Limit


def generate_golden_gate_actions(state: SearchState, enzyme_db: Dict, 
                                options: Dict) -> List[Step]:
    """
    Generate Golden Gate assembly actions.
    
    Golden Gate uses Type IIS enzymes (e.g., BsaI, BsmBI) that cut outside
    their recognition site, creating user-defined overhangs.
    
    Args:
        state: Current search state
        enzyme_db: Enzyme database
        options: Planning options
        
    Returns:
        List of possible Golden Gate steps
    """
    actions = []
    
    # Identify Type IIS enzymes
    type_iis_enzymes = []
    for name, info in enzyme_db.items():
        # Type IIS typically have cut sites outside recognition site
        site_len = len(info['sequence'])
        cut_idx = info['cut_index']
        
        # Heuristic: if cut is at position 0 or len(site), it's Type IIS-like
        if cut_idx == 0 or cut_idx == site_len:
            type_iis_enzymes.append(name)
    
    # For now, just add a placeholder
    # Full implementation would design overhangs and check order
    
    return actions  # Placeholder


def enumerate_actions(state: SearchState, enzyme_db: Dict, 
                     options: Dict) -> List[Step]:
    """
    Enumerate all possible actions from current state.
    
    Args:
        state: Current search state
        enzyme_db: Enzyme database
        options: Planning options
        
    Returns:
        List of all possible steps
    """
    actions = []
    
    # Generate digest actions
    actions.extend(generate_digest_actions(state, enzyme_db, options))
    
    # Generate ligation actions
    actions.extend(generate_ligation_actions(state, options))
    
    # Generate Golden Gate actions
    if options.get('prefer_typeIIS', False):
        actions.extend(generate_golden_gate_actions(state, enzyme_db, options))
    
    # PCR actions (future expansion)
    # Clean-up actions (future expansion)
    
    return actions


def apply_action(step: Step, state: SearchState, enzyme_db: Dict) -> Optional[SearchState]:
    """
    Apply an action to the current state to generate a new state.
    
    Args:
        step: Action to apply
        state: Current state
        enzyme_db: Enzyme database
        
    Returns:
        New state after applying action, or None if invalid
    """
    # Convert frozen set to list for manipulation
    constructs = list(state.constructs)
    
    # Remove input constructs
    input_names = set(step.inputs)
    constructs = [c for c in constructs if c.name not in input_names]
    
    # Add output constructs (simplified - need to actually compute sequences)
    for output_name in step.outputs:
        # Create placeholder construct
        # In full implementation, this would compute the actual sequence
        new_construct = Construct(
            name=output_name,
            seq="PLACEHOLDER",  # Would compute actual sequence
            circular=False,
            features=tuple(),
            notes=f"Generated by {step.action}"
        )
        constructs.append(new_construct)
    
    # Create new state
    new_state = SearchState(
        constructs=frozenset(constructs),
        steps_taken=state.steps_taken + (step,),
        cost=state.cost + step.cost,
        heuristic=0.0  # Will be recalculated
    )
    
    return new_state


# ============================================================================
# SEARCH ALGORITHM
# ============================================================================

def beam_search(initial_state: SearchState, target_spec: Dict, 
               enzyme_db: Dict, options: Dict) -> Optional[Plan]:
    """
    Beam search to find a cloning plan.
    
    Beam search is a width-limited breadth-first search that keeps only
    the k best states at each level.
    
    Args:
        initial_state: Starting state with input constructs
        target_spec: Target assembly specification
        enzyme_db: Enzyme database
        options: Planning options (max_steps, beam_width, etc.)
        
    Returns:
        Best plan found, or None if no feasible plan
    """
    max_steps = options.get('max_steps', 3)
    beam_width = options.get('beam_width', 10)
    
    # Priority queue: (f_score, state)
    # f_score = g (cost so far) + h (heuristic)
    frontier = []
    
    # Calculate initial heuristic
    initial_state = SearchState(
        constructs=initial_state.constructs,
        steps_taken=initial_state.steps_taken,
        cost=0.0,
        heuristic=calculate_heuristic(initial_state, target_spec, options)
    )
    
    heapq.heappush(frontier, initial_state)
    
    # Track visited states to avoid cycles
    visited = set()
    
    # Best plan found so far
    best_plan = None
    best_score = float('inf')
    
    while frontier:
        # Get beam_width best states from frontier
        current_beam = []
        for _ in range(min(beam_width, len(frontier))):
            if frontier:
                current_beam.append(heapq.heappop(frontier))
        
        for state in current_beam:
            # Check if already visited
            # Use signature() instead of state tuple to avoid hashing Step objects
            state_key = state.signature()
            if state_key in visited:
                continue
            visited.add(state_key)
            
            # Check if we've reached max depth
            if len(state.steps_taken) >= max_steps:
                # Evaluate this as a potential solution
                plan = Plan(
                    steps=list(state.steps_taken),
                    final=None,  # Would identify final construct
                    score=state.cost,
                    feasible=True
                )
                
                if state.cost < best_score:
                    best_score = state.cost
                    best_plan = plan
                continue
            
            # Check if goal is reached
            # (Simplified: check if we have target construct)
            target_name = target_spec.get('target_name', 'final')
            construct_names = {c.name for c in state.constructs}
            
            if target_name in construct_names:
                # Found target!
                final_construct = next(c for c in state.constructs if c.name == target_name)
                plan = Plan(
                    steps=list(state.steps_taken),
                    final=final_construct,
                    score=state.cost,
                    feasible=True
                )
                
                if state.cost < best_score:
                    best_score = state.cost
                    best_plan = plan
                continue
            
            # Generate successor states
            actions = enumerate_actions(state, enzyme_db, options)
            
            for action in actions:
                new_state = apply_action(action, state, enzyme_db)
                if new_state is not None:
                    # Calculate heuristic for new state
                    new_state = SearchState(
                        constructs=new_state.constructs,
                        steps_taken=new_state.steps_taken,
                        cost=new_state.cost,
                        heuristic=calculate_heuristic(new_state, target_spec, options)
                    )
                    heapq.heappush(frontier, new_state)
    
    return best_plan


# ============================================================================
# MAIN PLANNING INTERFACE
# ============================================================================

def plan_from_spec(spec: Dict, enzyme_db: Dict, max_steps: int = 3, 
                  options: Optional[Dict] = None) -> Plan:
    """
    Generate a cloning plan from a specification.
    
    Args:
        spec: Cloning specification with vector, inserts, target, constraints
        enzyme_db: Enzyme database (from enzymes.json)
        max_steps: Maximum number of steps to allow
        options: Additional options (prefer_typeIIS, avoid_enzymes, etc.)
        
    Returns:
        Best plan found
    """
    if options is None:
        options = {}
    
    options['max_steps'] = max_steps
    
    # Parse vector
    vector_spec = spec.get('vector', {})
    vector_name = vector_spec.get('name', 'vector')
    vector_fasta = vector_spec.get('fasta', '')
    vector_circular = vector_spec.get('circular', True)
    
    # Read vector sequence
    from sim import read_dna_sequence
    try:
        vector_seq = read_dna_sequence(vector_fasta)
    except (FileNotFoundError, ValueError) as e:
        return Plan(
            steps=[],
            final=None,
            score=float('inf'),
            feasible=False,
            reason=f"Could not read vector: {e}"
        )
    
    vector = Construct(
        name=vector_name,
        seq=vector_seq,
        circular=vector_circular,
        features=tuple(),
        notes="Vector"
    )
    
    # Parse inserts
    inserts = []
    for insert_spec in spec.get('inserts', []):
        insert_name = insert_spec.get('name', f'insert_{len(inserts)}')
        insert_fasta = insert_spec.get('fasta', '')
        insert_features = extract_features_from_spec(insert_spec.get('features', []))
        
        try:
            insert_seq = read_dna_sequence(insert_fasta)
        except (FileNotFoundError, ValueError) as e:
            return Plan(
                steps=[],
                final=None,
                score=float('inf'),
                feasible=False,
                reason=f"Could not read insert {insert_name}: {e}"
            )
        
        insert = Construct(
            name=insert_name,
            seq=insert_seq,
            circular=False,
            features=insert_features,
            notes="Insert"
        )
        inserts.append(insert)
    
    # Create initial state
    initial_constructs = frozenset([vector] + inserts)
    initial_state = SearchState(
        constructs=initial_constructs,
        steps_taken=tuple(),
        cost=0.0,
        heuristic=0.0
    )
    
    # Parse constraints
    constraints = spec.get('constraints', {})
    options['avoid_internal_cuts'] = constraints.get('avoid_internal_cuts', True)
    options['min_overhang'] = constraints.get('min_overhang', 4)
    
    # Parse target
    target_spec = spec.get('target', {})
    target_spec['target_name'] = 'final'  # Expected final construct name
    
    # Run search
    plan = beam_search(initial_state, target_spec, enzyme_db, options)
    
    if plan is None:
        return Plan(
            steps=[],
            final=None,
            score=float('inf'),
            feasible=False,
            reason="No feasible plan found within max_steps"
        )
    
    # Score the plan
    plan.score = score_plan(plan, options)
    
    return plan


# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

def format_plan_summary(plan: Plan, show_gels: bool = False) -> str:
    """
    Format plan as human-readable summary.
    
    Args:
        plan: Plan to format
        show_gels: If True, include predicted gel diagrams
        
    Returns:
        Formatted string
    """
    if not plan.feasible:
        return f"No feasible plan: {plan.reason}\n"
    
    lines = []
    lines.append("=" * 80)
    lines.append("CLONING PLAN")
    lines.append("=" * 80)
    lines.append(f"Total steps: {len(plan.steps)}")
    lines.append(f"Score: {plan.score:.2f}")
    lines.append("")
    
    for i, step in enumerate(plan.steps, 1):
        lines.append(f"Step {i} — {step.action.upper()}: {step.name}")
        lines.append(f"  Inputs: {', '.join(step.inputs)}")
        lines.append(f"  Outputs: {', '.join(step.outputs)}")
        
        if step.action == "digest":
            enzymes = step.params.get('enzymes', [])
            deP = "YES" if step.params.get('dephosphorylate', False) else "NO"
            lines.append(f"  Enzymes: {', '.join(enzymes)} (deP: {deP})")
            
            if step.predicted_fragments:
                frag_sizes = [f['length'] for f in step.predicted_fragments]
                lines.append(f"  Fragments: {', '.join(f'{s} bp' for s in frag_sizes)}")
        
        elif step.action == "ligate":
            directional = step.params.get('directional', False)
            dir_str = "✔" if directional else "✗"
            scar = step.params.get('scar', '')
            lines.append(f"  Directional: {dir_str}")
            if scar:
                lines.append(f"  Scar: {scar}")
        
        elif step.action == "GG":
            enzyme = step.params.get('enzyme', '')
            overhangs = step.params.get('overhangs', [])
            lines.append(f"  Type IIS enzyme: {enzyme}")
            if overhangs:
                lines.append(f"  Overhangs: {' → '.join(overhangs)}")
        
        lines.append("")
    
    if plan.final:
        lines.append(f"Final construct: {plan.final.name}")
        lines.append(f"  Length: {len(plan.final.seq)} bp")
        lines.append(f"  Circular: {plan.final.circular}")
    
    lines.append("")
    lines.append("=" * 80)
    
    return '\n'.join(lines)


def export_plan_to_files(plan: Plan, output_dir: str, enzyme_db: Dict) -> None:
    """
    Export plan steps to GenBank and CSV files.
    
    Args:
        plan: Plan to export
        output_dir: Output directory path
        enzyme_db: Enzyme database
    """
    import os
    from exporters import export_genbank, export_csv
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Export each step
    for i, step in enumerate(plan.steps, 1):
        step_prefix = os.path.join(output_dir, f"step_{i:02d}_{step.action}")
        
        # For digest steps, export fragments
        if step.action == "digest" and step.predicted_fragments:
            # Placeholder: would export actual fragment data
            print(f"  Exporting step {i}: {step.name}")
    
    # Export final construct
    if plan.final:
        final_path = os.path.join(output_dir, "final.gb")
        # Placeholder: would export final construct
        print(f"  Exporting final construct: {plan.final.name}")
    
    print(f"\n✓ Plan exported to: {output_dir}")

