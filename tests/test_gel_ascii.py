#!/usr/bin/env python3
"""
Tests for ASCII Agarose Gel Simulation
"""

import sys
import os
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fragment_calculator import (
    simulate_gel,
    gel_coefficients,
    calculate_migration_row,
    merge_bands,
    get_band_glyph,
    render_circular_topology_bands
)
from gel_ladders import get_ladder, get_available_ladders


class TestGelCoefficients:
    """Test gel migration coefficient calculation."""
    
    def test_gel_coefficients_returns_tuple(self):
        """Test that gel_coefficients returns a tuple of two floats."""
        result = gel_coefficients(1.0)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)
    
    def test_gel_coefficients_vary_by_percent(self):
        """Test that different percentages yield different coefficients."""
        coef_0_8 = gel_coefficients(0.8)
        coef_1_0 = gel_coefficients(1.0)
        coef_2_0 = gel_coefficients(2.0)
        
        # Coefficients should differ
        assert coef_0_8 != coef_1_0
        assert coef_1_0 != coef_2_0


class TestMigrationCalculation:
    """Test fragment migration calculation."""
    
    def test_calculate_migration_row_basic(self):
        """Test basic migration row calculation."""
        row = calculate_migration_row(500, 1.0, 24, 0.85)
        assert isinstance(row, int)
        assert 2 <= row < 24  # Should be below wells and within gel
    
    def test_small_fragments_migrate_further(self):
        """Test that smaller fragments migrate further (higher row number)."""
        row_small = calculate_migration_row(100, 1.0, 24, 0.85)
        row_large = calculate_migration_row(5000, 1.0, 24, 0.85)
        
        assert row_small > row_large  # Smaller fragments travel further
    
    def test_zero_bp_returns_zero(self):
        """Test that 0 bp returns row 0."""
        row = calculate_migration_row(0, 1.0, 24, 0.85)
        assert row == 0


class TestBandMerging:
    """Test band merging logic."""
    
    def test_merge_bands_basic(self):
        """Test basic band merging."""
        fragments = [500, 1000, 2000]
        bands = merge_bands(fragments, 20, 1.0, 24, 0.85)
        
        assert isinstance(bands, dict)
        # Should have 3 separate bands (no merging)
        total_fragments = sum(len(bp_list) for bp_list in bands.values())
        assert total_fragments == 3
    
    def test_merge_bands_collapses_close_fragments(self):
        """Test that close fragments get merged."""
        fragments = [498, 500, 502]  # All within 20bp
        bands = merge_bands(fragments, 20, 1.0, 24, 0.85)
        
        # Should merge into fewer bands
        assert len(bands) <= 2
    
    def test_merge_bands_empty_list(self):
        """Test that empty fragment list returns empty dict."""
        bands = merge_bands([], 20, 1.0, 24, 0.85)
        assert bands == {}


class TestBandGlyph:
    """Test band glyph selection."""
    
    def test_get_band_glyph_intensity_levels(self):
        """Test that different intensities yield different glyphs."""
        glyph1 = get_band_glyph(1)
        glyph2 = get_band_glyph(2)
        glyph3 = get_band_glyph(3)
        glyph4 = get_band_glyph(4)
        
        # All should be single characters
        assert len(glyph1) == 1
        assert len(glyph2) == 1
        assert len(glyph3) == 1
        assert len(glyph4) == 1
        
        # Higher intensity should be different
        assert glyph1 != glyph4


class TestCircularTopology:
    """Test circular DNA topology rendering."""
    
    def test_render_circular_topology_returns_bands(self):
        """Test that circular topology rendering returns band dict."""
        bands = render_circular_topology_bands(3000, "native", 1.0, 24, 0.85)
        
        assert isinstance(bands, dict)
        assert len(bands) >= 2  # Should have SC and OC at minimum
    
    def test_circular_topology_has_labeled_forms(self):
        """Test that circular topology includes SC and OC labels."""
        bands = render_circular_topology_bands(3000, "native", 1.0, 24, 0.85)
        
        # Check that we have tuples with labels
        all_bands = []
        for bp_list in bands.values():
            all_bands.extend(bp_list)
        
        labels = [item[0] for item in all_bands if isinstance(item, tuple)]
        assert 'SC' in labels or 'OC' in labels


class TestGelLadders:
    """Test ladder presets."""
    
    def test_get_available_ladders(self):
        """Test that we can get list of available ladders."""
        ladders = get_available_ladders()
        assert isinstance(ladders, list)
        assert len(ladders) > 0
        assert '1kb' in ladders or '100bp' in ladders
    
    def test_get_ladder_1kb(self):
        """Test loading 1kb ladder."""
        ladder = get_ladder('1kb')
        assert isinstance(ladder, list)
        assert len(ladder) > 0
        assert all(isinstance(x, int) for x in ladder)
        assert 1000 in ladder
    
    def test_get_ladder_100bp(self):
        """Test loading 100bp ladder."""
        ladder = get_ladder('100bp')
        assert isinstance(ladder, list)
        assert 100 in ladder
        assert 500 in ladder
    
    def test_get_ladder_invalid_raises_error(self):
        """Test that invalid ladder name raises error."""
        with pytest.raises(ValueError):
            get_ladder('invalid_ladder')


class TestGelSimulation:
    """Test complete gel simulation."""
    
    def test_basic_ladder_and_lane_render(self):
        """Test basic gel with ladder and one lane."""
        gel = simulate_gel(
            lanes=[{
                "label": "Demo",
                "fragments": [500, 1500, 3000],
                "topology": "linearized",
                "circular": False
            }],
            ladder_bp=[100, 200, 500, 1000, 1500, 2000],
            gel_percent=1.0,
            gel_length=24,
            gel_width=60
        )
        
        assert isinstance(gel, str)
        assert "Ladder" in gel
        assert "Demo" in gel
        
        # Must have some band glyphs
        glyph_count = gel.count("█") + gel.count("▮") + gel.count("•") + gel.count("·")
        assert glyph_count >= 3
    
    def test_merge_threshold_collapses_close_bands(self):
        """Test that close bands are merged and labeled with multiplicity."""
        gel = simulate_gel(
            lanes=[{
                "label": "Close",
                "fragments": [498, 510],
                "topology": "linearized",
                "circular": False
            }],
            ladder_bp=[500],
            merge_threshold=20,
            gel_length=24,
            gel_width=60
        )
        
        assert "Close" in gel
        # Should mention merged bands or show multiplicity marker
        assert "×" in gel or "2" in gel
    
    def test_circular_zero_cut_native_shows_topology_forms(self):
        """Test that circular DNA with 0 cuts shows SC/OC forms."""
        L = 3000
        lane = {
            "label": "Plasmid",
            "fragments": [L],
            "topology": "native",
            "circular": True
        }
        gel = simulate_gel(
            [lane],
            ladder_bp=[500, 1000, 2000, 3000],
            gel_length=24,
            gel_width=60
        )
        
        assert "Plasmid" in gel
        # Should mention SC or OC forms
        assert "SC" in gel or "OC" in gel
    
    def test_one_cut_auto_linearizes(self):
        """Test that one cut in auto mode creates linearized band."""
        L = 4000
        lane = {
            "label": "1cut",
            "fragments": [L],
            "topology": "auto",
            "circular": True
        }
        gel = simulate_gel(
            [lane],
            ladder_bp=[1000, 4000],
            gel_length=24,
            gel_width=60
        )
        
        assert "1cut" in gel
        assert str(L) in gel
    
    def test_multiple_lanes(self):
        """Test gel with multiple sample lanes."""
        gel = simulate_gel(
            lanes=[
                {"label": "Lane1", "fragments": [500, 1000], "topology": "linearized", "circular": False},
                {"label": "Lane2", "fragments": [750, 1500], "topology": "linearized", "circular": False},
            ],
            ladder_bp=[500, 1000, 1500],
            gel_length=24,
            gel_width=80
        )
        
        assert "Lane1" in gel
        assert "Lane2" in gel
        assert "Ladder" in gel
    
    def test_gel_with_smear(self):
        """Test that smear option doesn't crash."""
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_length=24,
            gel_width=60,
            smear="light"
        )
        
        assert isinstance(gel, str)
        assert len(gel) > 0
    
    def test_gel_percent_clamped(self):
        """Test that gel percent is clamped to valid range."""
        # Should not crash with out-of-range percent
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_percent=10.0,  # Too high
            gel_length=24,
            gel_width=60
        )
        
        assert isinstance(gel, str)
    
    def test_gel_contains_wells(self):
        """Test that gel renders wells at the top."""
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_length=24,
            gel_width=60
        )
        
        # Wells should be drawn with box characters
        assert "┏" in gel or "━" in gel or "┓" in gel
    
    def test_gel_contains_dye_front(self):
        """Test that gel renders dye front marker."""
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_length=24,
            gel_width=60,
            dye_front=0.85
        )
        
        # Dye front marked with ~
        assert "~" in gel
    
    def test_gel_legend_includes_agarose_percent(self):
        """Test that legend includes agarose percentage."""
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_percent=1.5,
            gel_length=24,
            gel_width=60
        )
        
        assert "Agarose" in gel
        assert "1.5" in gel


class TestGelSimulationEdgeCases:
    """Test edge cases for gel simulation."""
    
    def test_no_fragments(self):
        """Test gel with lane that has no fragments."""
        gel = simulate_gel(
            lanes=[{"label": "Empty", "fragments": [], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_length=24,
            gel_width=60
        )
        
        assert "Empty" in gel
    
    def test_very_small_gel(self):
        """Test that small gel dimensions are handled."""
        gel = simulate_gel(
            lanes=[{"label": "Test", "fragments": [1000], "topology": "linearized", "circular": False}],
            ladder_bp=[1000],
            gel_length=15,
            gel_width=40
        )
        
        assert isinstance(gel, str)
        assert len(gel) > 0
    
    def test_many_lanes(self):
        """Test gel with many lanes."""
        lanes = [
            {"label": f"Lane{i}", "fragments": [500 + i*100], "topology": "linearized", "circular": False}
            for i in range(5)
        ]
        
        gel = simulate_gel(
            lanes=lanes,
            ladder_bp=[500, 1000],
            gel_length=24,
            gel_width=120
        )
        
        assert isinstance(gel, str)
        for i in range(5):
            assert f"Lane{i}" in gel


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

