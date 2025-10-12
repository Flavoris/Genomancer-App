#!/usr/bin/env python3
"""
Test suite for planner utility functions
Tests spec loading, validation, formatting, and helper utilities.
"""

import pytest
import sys
import os
import json
import tempfile

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from planner_utils import (
    load_json, validate_spec, check_frame_preservation,
    validate_golden_gate_overhangs, design_golden_gate_overhangs,
    predict_gel_bands, format_plan_detailed, format_plan_json
)


# ============================================================================
# FIXTURES
# ============================================================================

@pytest.fixture
def valid_spec_dict():
    """Create a valid specification dictionary."""
    return {
        "vector": {
            "name": "pTest",
            "fasta": "test_vector.fasta",
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": "test_insert.fasta",
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": [
                {
                    "left": "pTest",
                    "right": "Insert1",
                    "directional": True,
                    "scar": "none"
                }
            ]
        },
        "constraints": {
            "avoid_internal_cuts": True,
            "min_overhang": 4
        }
    }


@pytest.fixture
def temp_json_file(valid_spec_dict):
    """Create a temporary JSON file with valid spec."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(valid_spec_dict, f)
        temp_path = f.name
    
    yield temp_path
    
    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


# ============================================================================
# SPEC LOADING TESTS
# ============================================================================

def test_load_valid_json(temp_json_file):
    """Test loading a valid JSON spec file."""
    spec = load_json(temp_json_file)
    
    assert 'vector' in spec
    assert 'inserts' in spec
    assert spec['vector']['name'] == 'pTest'


def test_load_invalid_json():
    """Test loading an invalid JSON file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        f.write("{ invalid json }")
        temp_path = f.name
    
    try:
        with pytest.raises(ValueError, match="Invalid JSON"):
            load_json(temp_path)
    finally:
        os.unlink(temp_path)


def test_load_nonexistent_file():
    """Test loading a file that doesn't exist."""
    from planner_utils import load_json_or_yaml
    
    with pytest.raises(FileNotFoundError):
        load_json_or_yaml("nonexistent_file.json")


# ============================================================================
# SPEC VALIDATION TESTS
# ============================================================================

def test_validate_valid_spec(valid_spec_dict):
    """Test validation of a valid specification."""
    valid, error_msg = validate_spec(valid_spec_dict)
    
    assert valid is True
    assert error_msg == ""


def test_validate_missing_vector():
    """Test validation catches missing vector."""
    invalid_spec = {
        "inserts": [],
        "target": {}
    }
    
    valid, error_msg = validate_spec(invalid_spec)
    
    assert valid is False
    assert "vector" in error_msg.lower()


def test_validate_missing_target():
    """Test validation catches missing target."""
    invalid_spec = {
        "vector": {"name": "pTest", "fasta": "test.fasta"},
        "inserts": []
    }
    
    valid, error_msg = validate_spec(invalid_spec)
    
    assert valid is False
    assert "target" in error_msg.lower()


def test_validate_invalid_vector_format():
    """Test validation catches invalid vector format."""
    invalid_spec = {
        "vector": {"name": "pTest"},  # Missing 'fasta'
        "inserts": [],
        "target": {"order": []}
    }
    
    valid, error_msg = validate_spec(invalid_spec)
    
    assert valid is False
    assert "fasta" in error_msg.lower()


def test_validate_unknown_part_in_order():
    """Test validation catches undefined parts in order."""
    invalid_spec = {
        "vector": {"name": "pTest", "fasta": "test.fasta"},
        "inserts": [{"name": "Insert1", "fasta": "i1.fasta"}],
        "target": {
            "order": ["pTest", "UnknownPart"],  # UnknownPart not defined
            "junctions": []
        }
    }
    
    valid, error_msg = validate_spec(invalid_spec)
    
    assert valid is False
    assert "unknown" in error_msg.lower()


def test_validate_invalid_junctions():
    """Test validation catches invalid junction specifications."""
    invalid_spec = {
        "vector": {"name": "pTest", "fasta": "test.fasta"},
        "inserts": [],
        "target": {
            "order": ["pTest"],
            "junctions": [
                {"left": "pTest"}  # Missing 'right'
            ]
        }
    }
    
    valid, error_msg = validate_spec(invalid_spec)
    
    assert valid is False
    assert "junction" in error_msg.lower()


# ============================================================================
# FRAME CHECKING TESTS
# ============================================================================

def test_frame_preservation_no_shift():
    """Test frame preservation with no shift."""
    left_seq = "ATGATGATG"  # ATG codons
    right_seq = "ATGATGATG"
    scar = ""  # No scar
    
    ok, reason = check_frame_preservation(left_seq, right_seq, scar, 0, 0)
    
    assert ok is True


def test_frame_preservation_with_shift():
    """Test frame preservation detects shift."""
    left_seq = "ATGATGATG"
    right_seq = "ATGATGATG"
    scar = "A"  # 1 bp scar causes frameshift
    
    ok, reason = check_frame_preservation(left_seq, right_seq, scar, 0, 0)
    
    assert ok is False
    assert "shift" in reason.lower()


def test_frame_preservation_3bp_scar():
    """Test frame preservation with 3 bp scar (no shift)."""
    left_seq = "ATGATGATG"
    right_seq = "ATGATGATG"
    scar = "GGC"  # 3 bp scar preserves frame
    
    ok, reason = check_frame_preservation(left_seq, right_seq, scar, 0, 0)
    
    assert ok is True


def test_frame_preservation_stop_codon():
    """Test detection of stop codons in junction."""
    left_seq = "ATGATGAT"
    right_seq = "AATGATG"
    scar = "TAG"  # Stop codon
    
    ok, reason = check_frame_preservation(left_seq, right_seq, scar, 0, 0)
    
    assert ok is False
    assert "stop" in reason.lower()


# ============================================================================
# GOLDEN GATE OVERHANG TESTS
# ============================================================================

def test_validate_valid_overhangs():
    """Test validation of valid Golden Gate overhangs."""
    overhangs = ["AATG", "AGGT", "GCTT"]
    
    valid, error_msg = validate_golden_gate_overhangs(overhangs)
    
    assert valid is True


def test_validate_duplicate_overhangs():
    """Test detection of duplicate overhangs."""
    overhangs = ["AATG", "AGGT", "AATG"]  # Duplicate
    
    valid, error_msg = validate_golden_gate_overhangs(overhangs)
    
    assert valid is False
    assert "duplicate" in error_msg.lower()


def test_validate_palindromic_overhangs():
    """Test detection of palindromic overhangs."""
    overhangs = ["AATT", "AGGT"]  # AATT is palindromic
    
    valid, error_msg = validate_golden_gate_overhangs(overhangs)
    
    assert valid is False
    assert "palindromic" in error_msg.lower()


def test_validate_complementary_overhangs():
    """Test detection of complementary overhang pairs."""
    overhangs = ["AATG", "CATT"]  # Complementary pair
    
    valid, error_msg = validate_golden_gate_overhangs(overhangs)
    
    assert valid is False
    assert "complementary" in error_msg.lower()


def test_design_golden_gate_overhangs():
    """Test designing Golden Gate overhangs."""
    overhangs = design_golden_gate_overhangs(num_junctions=3, avoid_palindromes=True)
    
    # Should return 4 overhangs for 3 junctions
    assert len(overhangs) == 4
    
    # All should be 4 bp
    assert all(len(oh) == 4 for oh in overhangs)
    
    # Should all be unique
    assert len(overhangs) == len(set(overhangs))


def test_design_too_many_overhangs():
    """Test requesting too many overhangs."""
    with pytest.raises(ValueError):
        design_golden_gate_overhangs(num_junctions=100, avoid_palindromes=True)


# ============================================================================
# GEL PREDICTION TESTS
# ============================================================================

def test_predict_gel_bands_empty():
    """Test gel prediction with no fragments."""
    result = predict_gel_bands([])
    
    assert "No fragments" in result


def test_predict_gel_bands_single():
    """Test gel prediction with single fragment."""
    result = predict_gel_bands([1000])
    
    assert "1000" in result
    assert "medium" in result.lower()


def test_predict_gel_bands_multiple():
    """Test gel prediction with multiple fragments."""
    result = predict_gel_bands([5000, 1000, 500])
    
    assert "5000" in result
    assert "1000" in result
    assert "500" in result


def test_predict_gel_bands_categorization():
    """Test gel band size categorization."""
    result = predict_gel_bands([15000, 7000, 2000, 750, 200])
    
    assert "very large" in result.lower()
    assert "large" in result.lower()
    assert "medium" in result.lower()
    assert "small" in result.lower()
    assert "very small" in result.lower()


# ============================================================================
# PLAN FORMATTING TESTS
# ============================================================================

def test_format_plan_json_infeasible():
    """Test JSON formatting of infeasible plan."""
    from planner import Plan
    
    plan = Plan(
        steps=[],
        final=None,
        score=float('inf'),
        feasible=False,
        reason="No solution found"
    )
    
    json_str = format_plan_json(plan)
    plan_dict = json.loads(json_str)
    
    assert plan_dict['feasible'] is False
    assert 'reason' in plan_dict


def test_format_plan_json_feasible():
    """Test JSON formatting of feasible plan."""
    from planner import Plan, Step, Construct
    
    step = Step(
        name="Test",
        action="digest",
        inputs=["vector"],
        params={'enzymes': ['EcoRI']},
        outputs=["frag1", "frag2"],
        predicted_fragments=[],
        cost=1.0
    )
    
    final = Construct("final", "ATCG", True)
    
    plan = Plan(
        steps=[step],
        final=final,
        score=1.5,
        feasible=True
    )
    
    json_str = format_plan_json(plan)
    plan_dict = json.loads(json_str)
    
    assert plan_dict['feasible'] is True
    assert len(plan_dict['steps']) == 1
    assert plan_dict['final']['name'] == 'final'


def test_format_plan_detailed():
    """Test detailed plan formatting."""
    from planner import Plan, Step, Construct
    
    step = Step(
        name="Digest",
        action="digest",
        inputs=["pTest"],
        params={'enzymes': ['EcoRI'], 'dephosphorylate': True},
        outputs=["frag1", "frag2"],
        predicted_fragments=[{'length': 3000}, {'length': 1000}],
        cost=1.0
    )
    
    final = Construct("final", "ATCG" * 100, True)
    
    plan = Plan(
        steps=[step],
        final=final,
        score=1.5,
        feasible=True
    )
    
    detailed = format_plan_detailed(plan)
    
    assert "DETAILED CLONING PROTOCOL" in detailed
    assert "STEP 1" in detailed
    assert "EcoRI" in detailed
    assert "3000 bp" in detailed
    assert "FINAL CONSTRUCT" in detailed


# ============================================================================
# EDGE CASE TESTS
# ============================================================================

def test_validate_spec_with_extra_fields(valid_spec_dict):
    """Test that validation allows extra fields."""
    valid_spec_dict['extra_field'] = 'extra_value'
    
    valid, error_msg = validate_spec(valid_spec_dict)
    
    # Should still be valid (extra fields are OK)
    assert valid is True


def test_empty_inserts_list():
    """Test validation with empty inserts list."""
    spec = {
        "vector": {"name": "pTest", "fasta": "test.fasta"},
        "inserts": [],  # Empty
        "target": {"order": ["pTest"], "junctions": []}
    }
    
    valid, error_msg = validate_spec(spec)
    
    # Should be valid (inserts are optional)
    assert valid is True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

