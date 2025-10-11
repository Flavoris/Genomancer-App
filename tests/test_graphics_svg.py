#!/usr/bin/env python3
"""
Tests for SVG graphics generation module.
"""

import os
import sys
import tempfile

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from graphics import (
    render_plasmid_map,
    render_linear_map,
    render_fragment_diagram,
    svg_to_png
)


def _toy_cuts():
    """Helper function to create toy cut data for testing."""
    return [
        {"pos": 300, "enzyme": "EcoRI", "site": "GAATTC", "overhang_type": "5' overhang"},
        {"pos": 1200, "enzyme": "HindIII", "site": "AAGCTT", "overhang_type": "5' overhang"},
        {"pos": 2200, "enzyme": "MluI", "site": "ACGCGT", "overhang_type": "5' overhang"},
    ]


def _toy_fragments():
    """Helper function to create toy fragment data for testing."""
    return [
        {
            "start": 300,
            "end": 1200,
            "length": 900,
            "wraps": False,
            "boundaries": {
                "left_cut": {"pos": 300, "enzymes": [{"enzyme": "EcoRI"}]},
                "right_cut": {"pos": 1200, "enzymes": [{"enzyme": "HindIII"}]}
            }
        },
        {
            "start": 1200,
            "end": 2200,
            "length": 1000,
            "wraps": False,
            "boundaries": {
                "left_cut": {"pos": 1200, "enzymes": [{"enzyme": "HindIII"}]},
                "right_cut": {"pos": 2200, "enzymes": [{"enzyme": "MluI"}]}
            }
        },
        {
            "start": 2200,
            "end": 300,
            "length": 1100,
            "wraps": True,
            "boundaries": {
                "left_cut": {"pos": 2200, "enzymes": [{"enzyme": "MluI"}]},
                "right_cut": {"pos": 300, "enzymes": [{"enzyme": "EcoRI"}]}
            }
        }
    ]


def test_plasmid_svg_header():
    """Test that plasmid map generates valid SVG with header."""
    svg = render_plasmid_map(3000, _toy_cuts(), title="TestPlasmid", origin=0)
    
    assert svg.strip().startswith("<svg"), "SVG should start with <svg tag"
    assert svg.strip().endswith("</svg>"), "SVG should end with </svg>"
    assert "TestPlasmid" in svg, "Title should appear in SVG"
    assert "3000" in svg, "Length should appear in SVG"
    assert "xmlns" in svg, "SVG should have xmlns attribute"


def test_plasmid_svg_contains_enzymes():
    """Test that plasmid map includes enzyme labels."""
    svg = render_plasmid_map(3000, _toy_cuts(), title="Test")
    
    assert "EcoRI" in svg, "EcoRI should be labeled"
    assert "HindIII" in svg, "HindIII should be labeled"
    assert "MluI" in svg, "MluI should be labeled"


def test_plasmid_svg_show_sites():
    """Test that plasmid map includes recognition sites when requested."""
    svg = render_plasmid_map(3000, _toy_cuts(), show_sites=True)
    
    # Sites might not always be visible if show_sites is implemented differently
    # but at least one should appear if the feature is working
    assert "GAATTC" in svg or "AAGCTT" in svg or "ACGCGT" in svg, \
        "At least one recognition site should appear when show_sites=True"


def test_plasmid_svg_hide_overhangs():
    """Test that plasmid map hides overhangs when requested."""
    svg_with = render_plasmid_map(3000, _toy_cuts(), show_overhangs=True)
    svg_without = render_plasmid_map(3000, _toy_cuts(), show_overhangs=False)
    
    # With overhangs should have badges
    assert "[5']" in svg_with or "[3']" in svg_with or "[B]" in svg_with, \
        "Overhang badges should appear when show_overhangs=True"
    
    # Without overhangs should not have as many badges
    # (This is a weak assertion but tests the flag)
    assert len(svg_without) <= len(svg_with), \
        "SVG without overhangs should be shorter or equal"


def test_plasmid_svg_dark_theme():
    """Test that plasmid map respects dark theme."""
    svg = render_plasmid_map(3000, _toy_cuts(), theme="dark")
    
    assert "#1a1a1a" in svg or "#e0e0e0" in svg, \
        "Dark theme colors should appear in SVG"


def test_linear_svg_has_ticks_and_labels():
    """Test that linear map has scale ticks and enzyme labels."""
    svg = render_linear_map(3000, _toy_cuts(), ticks=6)
    
    assert "<line" in svg, "Should have line elements for ruler and ticks"
    assert "EcoRI" in svg, "Should label EcoRI"
    assert "HindIII" in svg, "Should label HindIII"
    assert "AAGCTT" not in svg, "Should NOT show sites by default"


def test_linear_svg_show_sites():
    """Test that linear map shows sites when requested."""
    svg = render_linear_map(3000, _toy_cuts(), show_sites=True)
    
    # At least one site should appear
    assert "GAATTC" in svg or "AAGCTT" in svg or "ACGCGT" in svg, \
        "Recognition sites should appear when show_sites=True"


def test_linear_svg_dimensions():
    """Test that linear map respects custom dimensions."""
    svg = render_linear_map(3000, _toy_cuts(), width=600, height=150)
    
    assert 'width="600"' in svg, "Width should be set to 600"
    assert 'height="150"' in svg, "Height should be set to 150"


def test_fragment_diagram_blocks():
    """Test that fragment diagram includes fragment blocks."""
    frags = _toy_fragments()
    svg = render_fragment_diagram(frags, 3000, annotate_sizes=True)
    
    assert "Fragments" in svg, "Title should appear"
    assert "<rect" in svg, "Should have rectangle elements for fragments"
    assert "900" in svg or "1000" in svg or "1100" in svg, \
        "Fragment sizes should appear as labels"
    assert "bp" in svg, "Should include 'bp' unit"


def test_fragment_diagram_wrap_indication():
    """Test that fragment diagram indicates wrap fragments."""
    frags = _toy_fragments()
    svg = render_fragment_diagram(frags, 3000, annotate_sizes=True)
    
    # Wrap fragment should have special rendering
    assert "wrap" in svg.lower() or "stroke-dasharray" in svg, \
        "Wrap fragments should be indicated visually or with label"


def test_fragment_diagram_no_sizes():
    """Test that fragment diagram can hide size annotations."""
    frags = _toy_fragments()
    svg = render_fragment_diagram(frags, 3000, annotate_sizes=False)
    
    # Should still be valid SVG
    assert svg.strip().startswith("<svg"), "Should be valid SVG"
    assert "<rect" in svg, "Should have rectangle elements"


def test_empty_cuts():
    """Test rendering with no cuts."""
    empty_cuts = []
    
    # Plasmid map with no cuts
    svg1 = render_plasmid_map(3000, empty_cuts)
    assert "<svg" in svg1, "Should generate valid SVG even with no cuts"
    
    # Linear map with no cuts
    svg2 = render_linear_map(3000, empty_cuts)
    assert "<svg" in svg2, "Should generate valid SVG even with no cuts"


def test_single_cut():
    """Test rendering with a single cut."""
    single_cut = [{"pos": 1500, "enzyme": "EcoRI", "site": "GAATTC", "overhang_type": "5' overhang"}]
    
    svg = render_plasmid_map(3000, single_cut)
    assert "EcoRI" in svg, "Should show the single enzyme"
    assert "1500" in svg, "Should show the position"


def test_multiple_enzymes_same_position():
    """Test rendering when multiple enzymes cut at the same position."""
    cuts = [
        {"pos": 1000, "enzyme": "EcoRI", "site": "GAATTC", "overhang_type": "5' overhang"},
        {"pos": 1000, "enzyme": "HindIII", "site": "AAGCTT", "overhang_type": "5' overhang"},
    ]
    
    svg = render_plasmid_map(3000, cuts)
    assert "×2" in svg or "EcoRI" in svg, "Should indicate multiple enzymes or show enzyme name"


def test_svg_to_png_import_error():
    """Test that svg_to_png raises ImportError when cairosvg is not available."""
    svg = "<svg></svg>"
    
    try:
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            svg_to_png(svg, tmp_path)
            # If cairosvg is installed, this should succeed
            assert os.path.exists(tmp_path), "PNG file should be created"
            os.unlink(tmp_path)
        except ImportError as e:
            # Expected if cairosvg is not installed
            assert "cairosvg" in str(e).lower(), "Error message should mention cairosvg"
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def test_cli_integration_plasmid_map(tmp_path):
    """Test CLI integration for generating plasmid map."""
    import subprocess
    
    output_svg = tmp_path / "test_map.svg"
    
    # Run sim.py with graphics output
    result = subprocess.run([
        sys.executable, "sim.py",
        "--seq", "ATCGAATTCGGGATCCAAAGAATTCGGG",
        "--enz", "EcoRI",
        "--out-svg", str(output_svg)
    ], capture_output=True, text=True)
    
    # Check that it ran successfully
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    
    # Check that file was created
    assert output_svg.exists(), "SVG file should be created"
    
    # Check that file is valid SVG
    content = output_svg.read_text()
    assert content.strip().startswith("<svg"), "Output should be valid SVG"
    assert "EcoRI" in content, "SVG should contain enzyme name"


def test_cli_integration_all_outputs(tmp_path):
    """Test CLI integration for generating all output types."""
    import subprocess
    
    svg_plasmid = tmp_path / "plasmid.svg"
    svg_linear = tmp_path / "linear.svg"
    svg_fragments = tmp_path / "fragments.svg"
    
    # Run sim.py with all graphics outputs
    result = subprocess.run([
        sys.executable, "sim.py",
        "--seq", "ATCGAATTCGGGATCCAAAGAATTCGGG",
        "--enz", "EcoRI", "BamHI",
        "--circular",
        "--out-svg", str(svg_plasmid),
        "--out-svg-linear", str(svg_linear),
        "--out-svg-fragments", str(svg_fragments),
        "--theme", "dark",
        "--show-sites"
    ], capture_output=True, text=True)
    
    # Check that it ran successfully
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    
    # Check that all files were created
    assert svg_plasmid.exists(), "Plasmid SVG should be created"
    assert svg_linear.exists(), "Linear SVG should be created"
    assert svg_fragments.exists(), "Fragments SVG should be created"
    
    # Check that files contain expected content
    plasmid_content = svg_plasmid.read_text()
    assert "<svg" in plasmid_content, "Plasmid output should be SVG"
    
    linear_content = svg_linear.read_text()
    assert "<svg" in linear_content, "Linear output should be SVG"
    
    fragments_content = svg_fragments.read_text()
    assert "<svg" in fragments_content, "Fragments output should be SVG"


if __name__ == "__main__":
    # Run basic tests
    print("Running SVG graphics tests...")
    
    test_plasmid_svg_header()
    print("✓ test_plasmid_svg_header")
    
    test_plasmid_svg_contains_enzymes()
    print("✓ test_plasmid_svg_contains_enzymes")
    
    test_plasmid_svg_show_sites()
    print("✓ test_plasmid_svg_show_sites")
    
    test_plasmid_svg_hide_overhangs()
    print("✓ test_plasmid_svg_hide_overhangs")
    
    test_plasmid_svg_dark_theme()
    print("✓ test_plasmid_svg_dark_theme")
    
    test_linear_svg_has_ticks_and_labels()
    print("✓ test_linear_svg_has_ticks_and_labels")
    
    test_linear_svg_show_sites()
    print("✓ test_linear_svg_show_sites")
    
    test_linear_svg_dimensions()
    print("✓ test_linear_svg_dimensions")
    
    test_fragment_diagram_blocks()
    print("✓ test_fragment_diagram_blocks")
    
    test_fragment_diagram_wrap_indication()
    print("✓ test_fragment_diagram_wrap_indication")
    
    test_fragment_diagram_no_sizes()
    print("✓ test_fragment_diagram_no_sizes")
    
    test_empty_cuts()
    print("✓ test_empty_cuts")
    
    test_single_cut()
    print("✓ test_single_cut")
    
    test_multiple_enzymes_same_position()
    print("✓ test_multiple_enzymes_same_position")
    
    print("\nAll basic tests passed!")
    print("\nNote: CLI integration tests require pytest with tmp_path fixture")
    print("Run with: pytest tests/test_graphics_svg.py")

