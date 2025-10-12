#!/usr/bin/env python3
"""
Graphics Module for Restriction Enzyme Simulator
Generates SVG/PNG visualizations of plasmid maps, linear maps, and fragment diagrams.
"""

import math
import hashlib
from typing import List, Dict


def _hash_color(text: str) -> str:
    """
    Generate a deterministic color from a string using hashing.
    
    Args:
        text: String to hash (e.g., enzyme name)
        
    Returns:
        HSL color string
    """
    # Hash the text to get a deterministic hue
    hash_val = int(hashlib.md5(text.encode()).hexdigest()[:8], 16)
    hue = hash_val % 360
    # Keep saturation and lightness fixed for consistency
    return f"hsl({hue}, 65%, 50%)"


def _avoid_label_collision(
    labels: List[Dict],
    min_distance: float = 10.0,
    max_iterations: int = 10
) -> List[Dict]:
    """
    Adjust label positions to avoid collisions using an iterative approach.
    
    When multiple labels collide in sequence, this function uses an iterative
    algorithm to re-check all labels after each adjustment, preventing cascading
    collision issues.
    
    Args:
        labels: List of label dicts with 'x', 'y', 'text' keys
        min_distance: Minimum distance between labels
        max_iterations: Maximum number of adjustment passes
        
    Returns:
        List of adjusted labels with updated positions
    """
    if len(labels) <= 1:
        return labels
    
    # Copy labels to avoid modifying originals
    adjusted = [dict(label) for label in labels]
    
    # Iteratively resolve collisions
    for iteration in range(max_iterations):
        collision_found = False
        
        for i in range(len(adjusted)):
            # Check this label against ALL other labels
            for j in range(len(adjusted)):
                if i == j:
                    continue
                
                dx = adjusted[i]['x'] - adjusted[j]['x']
                dy = adjusted[i]['y'] - adjusted[j]['y']
                dist = math.sqrt(dx*dx + dy*dy)
                
                if dist < min_distance:
                    # Push this label outward radially
                    adjusted[i]['y'] += 12
                    adjusted[i]['needs_leader'] = True
                    collision_found = True
                    break  # Re-check in next iteration
        
        # If no collisions found, we're done
        if not collision_found:
            break
    
    # Fallback: If still colliding after max iterations, abbreviate overlapping labels
    if iteration == max_iterations - 1:
        for i in range(len(adjusted)):
            for j in range(i + 1, len(adjusted)):
                dx = adjusted[i]['x'] - adjusted[j]['x']
                dy = adjusted[i]['y'] - adjusted[j]['y']
                dist = math.sqrt(dx*dx + dy*dy)
                
                if dist < min_distance:
                    # Abbreviate labels as last resort
                    text = adjusted[i].get('text', '')
                    if len(text) > 5:
                        adjusted[i]['text'] = text[:4] + '…'
    
    return adjusted


def render_plasmid_map(
    L: int,
    cuts: List[Dict],
    *,
    title: str = "Plasmid",
    origin: int = 0,
    radius: int = 150,
    margin: int = 24,
    show_sites: bool = False,
    show_overhangs: bool = True,
    theme: str = "light"
) -> str:
    """
    Render a circular plasmid map with restriction sites.
    
    Args:
        L: Plasmid length in bp
        cuts: List of cut dictionaries with keys:
              - pos: int (cut position)
              - enzyme: str (enzyme name)
              - site: str (recognition sequence)
              - overhang_type: str (e.g., "5' overhang", "3' overhang", "Blunt")
        title: Map title
        origin: Origin position for angle calculation
        radius: Circle radius in pixels
        margin: Margin around the map
        show_sites: Include recognition sequences in labels
        show_overhangs: Include overhang type badges
        theme: "light" or "dark"
        
    Returns:
        SVG string
    """
    # Theme colors
    if theme == "dark":
        bg_color = "#1a1a1a"
        fg_color = "#e0e0e0"
        circle_color = "#555555"
        tick_color = "#888888"
    else:
        bg_color = "#ffffff"
        fg_color = "#000000"
        circle_color = "#cccccc"
        tick_color = "#999999"
    
    # Calculate SVG dimensions
    svg_size = 2 * (radius + margin + 80)  # Extra space for labels
    center = svg_size / 2
    
    # Start SVG
    svg_lines = []
    svg_lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                     f'width="{svg_size}" height="{svg_size}" '
                     f'viewBox="0 0 {svg_size} {svg_size}">')
    
    # Background
    svg_lines.append(f'  <rect width="{svg_size}" height="{svg_size}" fill="{bg_color}"/>')
    
    # Style
    svg_lines.append('  <style>')
    svg_lines.append(f'    text {{ font-family: Inter, "Segoe UI", Arial, sans-serif; fill: {fg_color}; }}')
    svg_lines.append('    .title { font-size: 16px; font-weight: bold; }')
    svg_lines.append('    .label { font-size: 11px; }')
    svg_lines.append(f'    .size {{ font-size: 10px; fill: {tick_color}; }}')
    svg_lines.append('    .badge { font-size: 9px; }')
    svg_lines.append('  </style>')
    
    # Title
    svg_lines.append(f'  <text x="{center}" y="{margin}" '
                     f'text-anchor="middle" class="title">{title}</text>')
    svg_lines.append(f'  <text x="{center}" y="{margin + 18}" '
                     f'text-anchor="middle" class="size">{L} bp</text>')
    
    # Draw main circle
    svg_lines.append(f'  <circle cx="{center}" cy="{center}" r="{radius}" '
                     f'fill="none" stroke="{circle_color}" stroke-width="2"/>')
    
    # Draw tick marks every 1000 bp or ceil(L/10)
    tick_interval = max(100, int(math.ceil(L / 10 / 100) * 100))
    for pos in range(0, L, tick_interval):
        angle_rad = 2 * math.pi * ((pos - origin) % L) / L - math.pi / 2
        x1 = center + radius * math.cos(angle_rad)
        y1 = center + radius * math.sin(angle_rad)
        
        # Small tick
        tick_len = 8 if pos % 1000 == 0 else 4
        x2 = center + (radius - tick_len) * math.cos(angle_rad)
        y2 = center + (radius - tick_len) * math.sin(angle_rad)
        
        svg_lines.append(f'  <line x1="{x1:.1f}" y1="{y1:.1f}" '
                        f'x2="{x2:.1f}" y2="{y2:.1f}" '
                        f'stroke="{tick_color}" stroke-width="1"/>')
        
        # Label major ticks
        if pos % 1000 == 0:
            label_radius = radius - 20
            lx = center + label_radius * math.cos(angle_rad)
            ly = center + label_radius * math.sin(angle_rad)
            svg_lines.append(f'  <text x="{lx:.1f}" y="{ly:.1f}" '
                           f'text-anchor="middle" dominant-baseline="middle" '
                           f'class="size">{pos}</text>')
    
    # Group cuts by position
    cuts_by_pos = {}
    for cut in cuts:
        pos = cut['pos']
        if pos not in cuts_by_pos:
            cuts_by_pos[pos] = []
        cuts_by_pos[pos].append(cut)
    
    # Draw cut markers and labels
    labels = []
    for pos, enzymes_at_pos in sorted(cuts_by_pos.items()):
        angle_rad = 2 * math.pi * ((pos - origin) % L) / L - math.pi / 2
        
        # Draw tick mark
        x1 = center + radius * math.cos(angle_rad)
        y1 = center + radius * math.sin(angle_rad)
        x2 = center + (radius + 15) * math.cos(angle_rad)
        y2 = center + (radius + 15) * math.sin(angle_rad)
        
        # Use color from first enzyme at this position
        color = _hash_color(enzymes_at_pos[0]['enzyme'])
        
        svg_lines.append(f'  <line x1="{x1:.1f}" y1="{y1:.1f}" '
                        f'x2="{x2:.1f}" y2="{y2:.1f}" '
                        f'stroke="{color}" stroke-width="2"/>')
        
        # Create label
        label_radius = radius + 30
        lx = center + label_radius * math.cos(angle_rad)
        ly = center + label_radius * math.sin(angle_rad)
        
        # Build label text
        if len(enzymes_at_pos) > 1:
            label_text = f"×{len(enzymes_at_pos)} @ {pos}"
        else:
            enzyme = enzymes_at_pos[0]
            parts = [enzyme['enzyme']]
            if show_sites:
                parts.append(f"({enzyme['site']})")
            label_text = f"{' '.join(parts)} @ {pos}"
        
        labels.append({
            'x': lx,
            'y': ly,
            'text': label_text,
            'color': color,
            'pos': pos,
            'angle': angle_rad,
            'enzymes': enzymes_at_pos
        })
    
    # Avoid label collisions
    labels = _avoid_label_collision(labels, min_distance=15.0)
    
    # Draw labels
    for label in labels:
        lx, ly = label['x'], label['y']
        
        # Draw leader line if needed
        if label.get('needs_leader', False):
            x_tick = center + (radius + 15) * math.cos(label['angle'])
            y_tick = center + (radius + 15) * math.sin(label['angle'])
            svg_lines.append(f'  <line x1="{x_tick:.1f}" y1="{y_tick:.1f}" '
                           f'x2="{lx:.1f}" y2="{ly:.1f}" '
                           f'stroke="{label["color"]}" stroke-width="0.5" '
                           f'stroke-dasharray="2,2"/>')
        
        svg_lines.append(f'  <text x="{lx:.1f}" y="{ly:.1f}" '
                        f'text-anchor="middle" class="label" '
                        f'fill="{label["color"]}">{label["text"]}</text>')
        
        # Add overhang badge if enabled
        if show_overhangs and len(label['enzymes']) == 1:
            overhang = label['enzymes'][0]['overhang_type']
            badge_map = {"5' overhang": "5'", "3' overhang": "3'", "Blunt": "B"}
            badge_text = badge_map.get(overhang, "?")
            svg_lines.append(f'  <text x="{lx:.1f}" y="{ly + 12:.1f}" '
                           f'text-anchor="middle" class="badge" '
                           f'fill="{label["color"]}">[{badge_text}]</text>')
    
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)


def render_linear_map(
    L: int,
    cuts: List[Dict],
    *,
    title: str = "Restriction Map",
    width: int = 900,
    height: int = 180,
    ticks: int = 10,
    show_sites: bool = False,
    show_overhangs: bool = True,
    theme: str = "light"
) -> str:
    """
    Render a linear restriction map with cut sites.
    
    Args:
        L: DNA length in bp
        cuts: List of cut dictionaries
        title: Map title
        width: SVG width in pixels
        height: SVG height in pixels
        ticks: Number of scale ticks
        show_sites: Include recognition sequences
        show_overhangs: Include overhang badges
        theme: "light" or "dark"
        
    Returns:
        SVG string
    """
    # Theme colors
    if theme == "dark":
        bg_color = "#1a1a1a"
        fg_color = "#e0e0e0"
        axis_color = "#555555"
    else:
        bg_color = "#ffffff"
        fg_color = "#000000"
        axis_color = "#333333"
    
    margin_x = 60
    ruler_y = height / 2
    ruler_width = width - 2 * margin_x
    
    # Start SVG
    svg_lines = []
    svg_lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                     f'width="{width}" height="{height}" '
                     f'viewBox="0 0 {width} {height}">')
    
    # Background
    svg_lines.append(f'  <rect width="{width}" height="{height}" fill="{bg_color}"/>')
    
    # Style
    svg_lines.append('  <style>')
    svg_lines.append(f'    text {{ font-family: Inter, "Segoe UI", Arial, sans-serif; fill: {fg_color}; }}')
    svg_lines.append('    .title { font-size: 14px; font-weight: bold; }')
    svg_lines.append('    .label { font-size: 10px; }')
    svg_lines.append('    .tick-label { font-size: 9px; }')
    svg_lines.append('  </style>')
    
    # Title
    svg_lines.append(f'  <text x="{width/2}" y="20" text-anchor="middle" class="title">{title}</text>')
    svg_lines.append(f'  <text x="{width/2}" y="35" text-anchor="middle" class="tick-label">{L} bp</text>')
    
    # Draw ruler
    svg_lines.append(f'  <line x1="{margin_x}" y1="{ruler_y}" '
                     f'x2="{width - margin_x}" y2="{ruler_y}" '
                     f'stroke="{axis_color}" stroke-width="2"/>')
    
    # Draw scale ticks
    for i in range(ticks + 1):
        pos = int(i * L / ticks) if ticks > 0 else 0
        x = margin_x + (i * ruler_width / ticks) if ticks > 0 else margin_x
        
        # Tick mark
        svg_lines.append(f'  <line x1="{x:.1f}" y1="{ruler_y - 5}" '
                        f'x2="{x:.1f}" y2="{ruler_y + 5}" '
                        f'stroke="{axis_color}" stroke-width="1"/>')
        
        # Tick label
        svg_lines.append(f'  <text x="{x:.1f}" y="{ruler_y + 20}" '
                        f'text-anchor="middle" class="tick-label">{pos}</text>')
    
    # Group cuts by position for collision detection
    cuts_by_pos = {}
    for cut in cuts:
        pos = cut['pos']
        if pos not in cuts_by_pos:
            cuts_by_pos[pos] = []
        cuts_by_pos[pos].append(cut)
    
    # Draw cut markers
    cut_positions = sorted(cuts_by_pos.keys())
    rows = [[], []]  # Two rows for collision avoidance (top and bottom)
    
    for pos in cut_positions:
        enzymes_at_pos = cuts_by_pos[pos]
        x = margin_x + (pos / L * ruler_width)
        color = _hash_color(enzymes_at_pos[0]['enzyme'])
        
        # Determine row (alternate top/bottom)
        row = 0
        if rows[0]:
            last_x = rows[0][-1]['x']
            if abs(x - last_x) < 80:  # Collision threshold
                row = 1
        
        # Draw marker triangle
        marker_y = ruler_y - 25 if row == 0 else ruler_y + 25
        points = f"{x:.1f},{marker_y:.1f} {x-4:.1f},{marker_y-8 if row == 0 else marker_y+8:.1f} {x+4:.1f},{marker_y-8 if row == 0 else marker_y+8:.1f}"
        svg_lines.append(f'  <polygon points="{points}" fill="{color}"/>')
        
        # Build label
        if len(enzymes_at_pos) > 1:
            label_text = " • ".join([e['enzyme'] for e in enzymes_at_pos])
        else:
            enzyme = enzymes_at_pos[0]
            parts = [enzyme['enzyme']]
            if show_sites:
                parts.append(f"({enzyme['site']})")
            if show_overhangs:
                overhang = enzyme['overhang_type']
                badge_map = {"5' overhang": "5'", "3' overhang": "3'", "Blunt": "B"}
                badge = badge_map.get(overhang, "?")
                parts.append(f"[{badge}]")
            label_text = " ".join(parts)
        
        # Draw label
        label_y = marker_y - 15 if row == 0 else marker_y + 20
        svg_lines.append(f'  <text x="{x:.1f}" y="{label_y:.1f}" '
                        f'text-anchor="middle" class="label" fill="{color}">{label_text}</text>')
        
        rows[row].append({'x': x, 'pos': pos})
    
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)


def render_fragment_diagram(
    fragments: List[Dict],
    L: int,
    *,
    title: str = "Fragments",
    width: int = 900,
    height: int = 140,
    theme: str = "light",
    annotate_sizes: bool = True
) -> str:
    """
    Render a fragment diagram showing linearized intervals.
    
    Args:
        fragments: List of fragment dicts with keys:
                   - start: int
                   - end: int
                   - length: int
                   - wraps: bool
                   - boundaries: dict
        L: Total DNA length
        title: Diagram title
        width: SVG width in pixels
        height: SVG height in pixels
        theme: "light" or "dark"
        annotate_sizes: Show size labels on fragments
        
    Returns:
        SVG string
    """
    # Theme colors
    if theme == "dark":
        bg_color = "#1a1a1a"
        fg_color = "#e0e0e0"
        frag_color = "#4a90e2"
        wrap_color = "#e24a4a"
    else:
        bg_color = "#ffffff"
        fg_color = "#000000"
        frag_color = "#6ba3e8"
        wrap_color = "#e87c7c"
    
    margin_x = 60
    bar_height = 30
    bar_y = height / 2 - bar_height / 2
    ruler_width = width - 2 * margin_x
    
    # Start SVG
    svg_lines = []
    svg_lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                     f'width="{width}" height="{height}" '
                     f'viewBox="0 0 {width} {height}">')
    
    # Background
    svg_lines.append(f'  <rect width="{width}" height="{height}" fill="{bg_color}"/>')
    
    # Style
    svg_lines.append('  <style>')
    svg_lines.append(f'    text {{ font-family: Inter, "Segoe UI", Arial, sans-serif; fill: {fg_color}; }}')
    svg_lines.append('    .title { font-size: 14px; font-weight: bold; }')
    svg_lines.append('    .size-label { font-size: 10px; }')
    svg_lines.append('  </style>')
    
    # Title
    svg_lines.append(f'  <text x="{width/2}" y="20" text-anchor="middle" class="title">{title}</text>')
    
    # Draw fragments
    for frag in fragments:
        start = frag['start']
        end = frag['end']
        length = frag['length']
        wraps = frag['wraps']
        
        color = wrap_color if wraps else frag_color
        
        if wraps:
            # Draw split block for wrap fragment
            # Part 1: from start to L
            x1 = margin_x + (start / L * ruler_width)
            w1 = (L - start) / L * ruler_width
            svg_lines.append(f'  <rect x="{x1:.1f}" y="{bar_y}" '
                           f'width="{w1:.1f}" height="{bar_height}" '
                           f'fill="{color}" opacity="0.7" '
                           f'stroke="{fg_color}" stroke-width="1"/>')
            
            # Part 2: from 0 to end
            x2 = margin_x
            w2 = (end / L * ruler_width)
            svg_lines.append(f'  <rect x="{x2:.1f}" y="{bar_y}" '
                           f'width="{w2:.1f}" height="{bar_height}" '
                           f'fill="{color}" opacity="0.7" '
                           f'stroke="{fg_color}" stroke-width="1"/>')
            
            # Dashed line connecting the two parts
            mid_y = bar_y + bar_height / 2
            svg_lines.append(f'  <line x1="{x1 + w1:.1f}" y1="{mid_y}" '
                           f'x2="{margin_x:.1f}" y2="{mid_y}" '
                           f'stroke="{color}" stroke-width="1" '
                           f'stroke-dasharray="4,4"/>')
            
            # Size label in the middle
            if annotate_sizes:
                label_x = width / 2
                label_y = bar_y - 10
                svg_lines.append(f'  <text x="{label_x:.1f}" y="{label_y:.1f}" '
                               f'text-anchor="middle" class="size-label" '
                               f'fill="{color}">{length:,} bp (wrap)</text>')
        else:
            # Normal fragment
            x = margin_x + (start / L * ruler_width)
            w = (length / L * ruler_width)
            
            svg_lines.append(f'  <rect x="{x:.1f}" y="{bar_y}" '
                           f'width="{w:.1f}" height="{bar_height}" '
                           f'fill="{color}" opacity="0.7" '
                           f'stroke="{fg_color}" stroke-width="1"/>')
            
            # Size label
            if annotate_sizes:
                label_x = x + w / 2
                label_y = bar_y + bar_height / 2 + 4
                svg_lines.append(f'  <text x="{label_x:.1f}" y="{label_y:.1f}" '
                               f'text-anchor="middle" class="size-label" '
                               f'dominant-baseline="middle">{length:,} bp</text>')
    
    # Draw scale ruler at bottom
    ruler_y = bar_y + bar_height + 30
    svg_lines.append(f'  <line x1="{margin_x}" y1="{ruler_y}" '
                     f'x2="{width - margin_x}" y2="{ruler_y}" '
                     f'stroke="{fg_color}" stroke-width="1" opacity="0.3"/>')
    
    # Ruler ticks
    for i in range(11):
        pos = int(i * L / 10)
        x = margin_x + (i * ruler_width / 10)
        svg_lines.append(f'  <line x1="{x:.1f}" y1="{ruler_y}" '
                        f'x2="{x:.1f}" y2="{ruler_y + 5}" '
                        f'stroke="{fg_color}" stroke-width="1" opacity="0.3"/>')
        svg_lines.append(f'  <text x="{x:.1f}" y="{ruler_y + 18}" '
                        f'text-anchor="middle" class="size-label" '
                        f'opacity="0.5">{pos}</text>')
    
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)


def svg_to_png(svg: str, out_path: str, *, scale: float = 2.0) -> None:
    """
    Convert SVG string to PNG file.
    
    Args:
        svg: SVG string content
        out_path: Output PNG file path
        scale: Scaling factor for resolution (default 2.0 for high-DPI)
        
    Raises:
        ImportError: If cairosvg is not installed
    """
    try:
        import cairosvg
    except ImportError:
        raise ImportError(
            "cairosvg is required for PNG export. Install it with:\n"
            "  pip install cairosvg\n"
            "Note: cairosvg requires cairo library to be installed on your system."
        )
    
    cairosvg.svg2png(
        bytestring=svg.encode('utf-8'),
        write_to=out_path,
        scale=scale
    )

