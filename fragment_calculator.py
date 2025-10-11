#!/usr/bin/env python3
"""
Fragment Calculator Module
Handles fragment computation for both linear and circular DNA digestion.
"""

from typing import List, Dict, Tuple


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

