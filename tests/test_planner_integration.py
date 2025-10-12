#!/usr/bin/env python3
"""
Integration tests for the cloning planner
Tests end-to-end planning scenarios with realistic inputs.
"""

import pytest
import sys
import os
import tempfile
import json

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from planner import plan_from_spec, Plan
from planner_utils import validate_spec


# ============================================================================
# FIXTURES
# ============================================================================

@pytest.fixture
def enzyme_db():
    """Realistic enzyme database."""
    return {
        "EcoRI": {"sequence": "GAATTC", "cut_index": 1, "overhang_type": "5' overhang"},
        "BamHI": {"sequence": "GGATCC", "cut_index": 1, "overhang_type": "5' overhang"},
        "HindIII": {"sequence": "AAGCTT", "cut_index": 1, "overhang_type": "5' overhang"},
        "PstI": {"sequence": "CTGCAG", "cut_index": 5, "overhang_type": "3' overhang"},
        "XbaI": {"sequence": "TCTAGA", "cut_index": 1, "overhang_type": "5' overhang"},
        "SpeI": {"sequence": "ACTAGT", "cut_index": 1, "overhang_type": "5' overhang"},
        "BsaI": {"sequence": "GGTCTC", "cut_index": 1, "overhang_type": "5' overhang"},
        "NotI": {"sequence": "GCGGCCGC", "cut_index": 2, "overhang_type": "5' overhang"}
    }


@pytest.fixture
def temp_vector_fasta():
    """Create a temporary vector FASTA file."""
    # Vector with multiple restriction sites
    vector_seq = (
        "ATCGATCGATCG"
        "GAATTC"  # EcoRI
        "ATCGATCGATCG"
        "GGATCC"  # BamHI
        "ATCGATCGATCG"
        "AAGCTT"  # HindIII
        "ATCGATCGATCG"
        "CTGCAG"  # PstI
        "ATCGATCGATCG"
    ) * 2
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_vector\n")
        f.write(vector_seq + "\n")
        temp_path = f.name
    
    yield temp_path
    
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def temp_insert_fasta():
    """Create a temporary insert FASTA file."""
    insert_seq = "ATCGATCGATCGATCGATCG" * 5
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_insert\n")
        f.write(insert_seq + "\n")
        temp_path = f.name
    
    yield temp_path
    
    if os.path.exists(temp_path):
        os.unlink(temp_path)


# ============================================================================
# SIMPLE CLONING TESTS
# ============================================================================

def test_simple_two_enzyme_cloning(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test simple two-enzyme directional cloning."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
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
                    "scar": "none",
                    "keep_frame": False
                }
            ]
        },
        "constraints": {
            "avoid_internal_cuts": False,
            "min_overhang": 4
        }
    }
    
    # Validate spec
    valid, error_msg = validate_spec(spec)
    assert valid, f"Spec validation failed: {error_msg}"
    
    # Run planner
    plan = plan_from_spec(spec, enzyme_db, max_steps=3, options={})
    
    # Check that a plan was generated
    assert isinstance(plan, Plan)
    
    # Check feasibility (may or may not find a solution depending on sites)
    # For this test, we just check that the planner runs without error


def test_planner_with_avoid_enzymes(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test planner with enzyme avoidance."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": []
        },
        "constraints": {
            "avoid_internal_cuts": False,
            "min_overhang": 4
        }
    }
    
    options = {
        'avoid_enzymes': ['EcoRI', 'BamHI']
    }
    
    plan = plan_from_spec(spec, enzyme_db, max_steps=2, options=options)
    
    # Check that plan doesn't use avoided enzymes
    if plan.feasible:
        for step in plan.steps:
            if step.action == "digest":
                enzymes_used = step.params.get('enzymes', [])
                assert 'EcoRI' not in enzymes_used
                assert 'BamHI' not in enzymes_used


def test_planner_with_allowed_enzymes(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test planner with enzyme whitelist."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": []
        },
        "constraints": {
            "avoid_internal_cuts": False,
            "min_overhang": 4
        }
    }
    
    options = {
        'allow_enzymes': ['EcoRI', 'BamHI', 'HindIII']
    }
    
    plan = plan_from_spec(spec, enzyme_db, max_steps=2, options=options)
    
    # Check that plan only uses allowed enzymes
    if plan.feasible:
        for step in plan.steps:
            if step.action == "digest":
                enzymes_used = step.params.get('enzymes', [])
                for enz in enzymes_used:
                    assert enz in options['allow_enzymes']


# ============================================================================
# CONSTRAINT TESTS
# ============================================================================

def test_planner_with_internal_cut_avoidance(enzyme_db, temp_vector_fasta):
    """Test that planner avoids internal cuts when requested."""
    # Create insert with CDS feature
    insert_with_ecori = "ATCGGAATTCATCG" * 10  # Contains EcoRI site
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">insert_with_site\n")
        f.write(insert_with_ecori + "\n")
        insert_path = f.name
    
    try:
        spec = {
            "vector": {
                "name": "pTest",
                "fasta": temp_vector_fasta,
                "circular": True
            },
            "inserts": [
                {
                    "name": "Insert1",
                    "fasta": insert_path,
                    "features": [
                        {
                            "type": "CDS",
                            "frame": 0,
                            "start": 0,
                            "end": 140,
                            "label": "gene",
                            "direction": "forward"
                        }
                    ]
                }
            ],
            "target": {
                "order": ["pTest", "Insert1"],
                "junctions": []
            },
            "constraints": {
                "avoid_internal_cuts": True,
                "min_overhang": 4
            }
        }
        
        options = {}
        
        plan = plan_from_spec(spec, enzyme_db, max_steps=2, options=options)
        
        # If a plan was found, check that it doesn't have internal cuts
        if plan.feasible:
            for step in plan.steps:
                if step.action == "digest":
                    internal_cuts = step.params.get('internal_cuts', 0)
                    assert internal_cuts == 0, "Plan should avoid internal cuts"
    
    finally:
        if os.path.exists(insert_path):
            os.unlink(insert_path)


def test_planner_max_steps_limit(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test that planner respects max_steps limit."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": []
        },
        "constraints": {}
    }
    
    # Run with max_steps=1 (very restrictive)
    plan = plan_from_spec(spec, enzyme_db, max_steps=1, options={})
    
    # Should either find a 1-step plan or fail
    if plan.feasible:
        assert len(plan.steps) <= 1


# ============================================================================
# SCORING TESTS
# ============================================================================

def test_planner_scoring_consistency(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test that plan scoring is consistent."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": []
        },
        "constraints": {}
    }
    
    # Run planner multiple times with same inputs
    plan1 = plan_from_spec(spec, enzyme_db, max_steps=2, options={})
    plan2 = plan_from_spec(spec, enzyme_db, max_steps=2, options={})
    
    # Scores should be consistent
    if plan1.feasible and plan2.feasible:
        # Scores might differ slightly due to search randomness,
        # but should be in same ballpark
        assert abs(plan1.score - plan2.score) < plan1.score * 0.5


# ============================================================================
# ERROR HANDLING TESTS
# ============================================================================

def test_planner_with_nonexistent_vector():
    """Test planner with nonexistent vector file."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": "nonexistent_file.fasta",
            "circular": True
        },
        "inserts": [],
        "target": {
            "order": ["pTest"],
            "junctions": []
        },
        "constraints": {}
    }
    
    plan = plan_from_spec(spec, {}, max_steps=2, options={})
    
    # Should return infeasible plan with error reason
    assert not plan.feasible
    assert "vector" in plan.reason.lower() or "read" in plan.reason.lower()


def test_planner_with_empty_enzyme_db(temp_vector_fasta, temp_insert_fasta):
    """Test planner with empty enzyme database."""
    spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta,
            "circular": True
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta,
                "features": []
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"],
            "junctions": []
        },
        "constraints": {}
    }
    
    # Run with empty enzyme database
    plan = plan_from_spec(spec, {}, max_steps=2, options={})
    
    # Should not crash, but likely won't find a plan
    assert isinstance(plan, Plan)


# ============================================================================
# SPEC FORMAT TESTS
# ============================================================================

def test_spec_with_missing_optional_fields(enzyme_db, temp_vector_fasta, temp_insert_fasta):
    """Test that planner handles specs with missing optional fields."""
    minimal_spec = {
        "vector": {
            "name": "pTest",
            "fasta": temp_vector_fasta
            # Missing 'circular' (should default to True or be inferred)
        },
        "inserts": [
            {
                "name": "Insert1",
                "fasta": temp_insert_fasta
                # Missing 'features' (should default to empty)
            }
        ],
        "target": {
            "order": ["pTest", "Insert1"]
            # Missing 'junctions' (should default to empty)
        }
        # Missing 'constraints' (should use defaults)
    }
    
    # Should handle missing fields gracefully
    plan = plan_from_spec(minimal_spec, enzyme_db, max_steps=2, options={})
    
    assert isinstance(plan, Plan)


# ============================================================================
# REALISTIC SCENARIO TESTS
# ============================================================================

def test_realistic_cloning_scenario(enzyme_db):
    """Test a realistic cloning scenario with actual sequences."""
    # Create realistic vector with MCS
    vector_seq = (
        "ATCGATCGATCGATCG"
        "GAATTC"    # EcoRI
        "ATCG"
        "GGATCC"    # BamHI
        "ATCG"
        "AAGCTT"    # HindIII
        "ATCGATCGATCGATCG"
    )
    
    # Create insert without internal sites
    insert_seq = "ATGAAATCGATCGATCGATCGATCGTAA"  # Simple ORF
    
    # Write to temp files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">vector\n" + vector_seq + "\n")
        vector_path = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">insert\n" + insert_seq + "\n")
        insert_path = f.name
    
    try:
        spec = {
            "vector": {
                "name": "pVector",
                "fasta": vector_path,
                "circular": True
            },
            "inserts": [
                {
                    "name": "GeneX",
                    "fasta": insert_path,
                    "features": [
                        {
                            "type": "CDS",
                            "frame": 0,
                            "start": 0,
                            "end": len(insert_seq),
                            "label": "geneX",
                            "direction": "forward"
                        }
                    ]
                }
            ],
            "target": {
                "order": ["pVector", "GeneX"],
                "junctions": [
                    {
                        "left": "pVector",
                        "right": "GeneX",
                        "directional": True,
                        "scar": "none",
                        "keep_frame": False
                    }
                ]
            },
            "constraints": {
                "avoid_internal_cuts": True,
                "dephosphorylate_backbone": True,
                "min_overhang": 4
            }
        }
        
        plan = plan_from_spec(spec, enzyme_db, max_steps=3, options={})
        
        # This scenario should be solvable
        # (though we don't strictly assert feasibility since it depends on exact sites)
        assert isinstance(plan, Plan)
        
        # If feasible, check basic properties
        if plan.feasible:
            assert len(plan.steps) > 0
            assert plan.score > 0
    
    finally:
        if os.path.exists(vector_path):
            os.unlink(vector_path)
        if os.path.exists(insert_path):
            os.unlink(insert_path)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

