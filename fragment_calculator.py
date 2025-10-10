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

