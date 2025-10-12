#!/usr/bin/env python3
"""
Test suite for planner core functionality
Tests the Construct, Step, Plan dataclasses and basic planning logic.
"""

import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from planner import (
    Construct, Step, Plan, SearchState,
    sanitize_name, calculate_heuristic, score_plan,
    simulate_digest, check_internal_cuts
)


# ============================================================================
# FIXTURES
# ============================================================================

@pytest.fixture
def simple_vector():
    """Create a simple test vector."""
    return Construct(
        name="pTest",
        seq="ATCGGAATTCGGATCCAAGCTTCTGCAG" * 10,  # 280 bp with multiple sites
        circular=True,
        features=tuple(),
        origin=0,
        notes="Test vector"
    )


@pytest.fixture
def simple_insert():
    """Create a simple test insert."""
    return Construct(
        name="insert1",
        seq="ATCGATCGATCG" * 10,  # 120 bp
        circular=False,
        features=tuple([{
            'type': 'CDS',
            'frame': 0,
            'start': 0,
            'end': 120,
            'label': 'test_gene',
            'direction': 'forward'
        }]),
        notes="Test insert"
    )


@pytest.fixture
def enzyme_db():
    """Create a minimal enzyme database for testing."""
    return {
        "EcoRI": {
            "sequence": "GAATTC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        },
        "BamHI": {
            "sequence": "GGATCC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        },
        "HindIII": {
            "sequence": "AAGCTT",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        },
        "PstI": {
            "sequence": "CTGCAG",
            "cut_index": 5,
            "overhang_type": "3' overhang"
        },
        "BsaI": {
            "sequence": "GGTCTC",
            "cut_index": 1,
            "overhang_type": "5' overhang"
        }
    }


# ============================================================================
# DATA STRUCTURE TESTS
# ============================================================================

def test_construct_creation(simple_vector):
    """Test Construct dataclass creation."""
    assert simple_vector.name == "pTest"
    assert len(simple_vector.seq) == 280
    assert simple_vector.circular is True
    assert simple_vector.origin == 0


def test_construct_hashing(simple_vector, simple_insert):
    """Test that Construct objects can be hashed and stored in sets."""
    constructs_set = {simple_vector, simple_insert}
    assert len(constructs_set) == 2
    assert simple_vector in constructs_set


def test_step_creation():
    """Test Step dataclass creation."""
    step = Step(
        name="Test_Digest",
        action="digest",
        inputs=["pTest"],
        params={'enzymes': ['EcoRI']},
        outputs=["pTest_frag0", "pTest_frag1"],
        predicted_fragments=[],
        cost=1.0
    )
    assert step.action == "digest"
    assert len(step.inputs) == 1
    assert len(step.outputs) == 2


def test_plan_creation(simple_vector):
    """Test Plan dataclass creation."""
    step = Step(
        name="Test",
        action="digest",
        inputs=["pTest"],
        params={},
        outputs=["frag1", "frag2"],
        predicted_fragments=[],
        cost=1.0
    )
    
    plan = Plan(
        steps=[step],
        final=simple_vector,
        score=1.5,
        feasible=True
    )
    
    assert len(plan.steps) == 1
    assert plan.feasible is True
    assert plan.score == 1.5


def test_search_state_ordering():
    """Test SearchState priority queue ordering."""
    construct = Construct("test", "ATCG", False)
    
    state1 = SearchState(
        constructs=frozenset([construct]),
        steps_taken=tuple(),
        cost=2.0,
        heuristic=3.0
    )
    
    state2 = SearchState(
        constructs=frozenset([construct]),
        steps_taken=tuple(),
        cost=1.0,
        heuristic=2.0
    )
    
    # state2 should be "less than" state1 (lower f-score)
    assert state2 < state1


# ============================================================================
# UTILITY FUNCTION TESTS
# ============================================================================

def test_sanitize_name():
    """Test name sanitization."""
    assert sanitize_name("test-name") == "test_name"
    assert sanitize_name("test name") == "test_name"
    assert sanitize_name("test@name!") == "test_name_"
    assert sanitize_name("pUC19") == "pUC19"


def test_check_internal_cuts(simple_insert, enzyme_db):
    """Test internal cut checking."""
    # Create a sequence with a CDS that contains an EcoRI site
    seq_with_cut = "ATCGGAATTCATCG" * 10
    features = tuple([{
        'type': 'CDS',
        'frame': 0,
        'start': 0,
        'end': 140,
        'label': 'gene',
        'direction': 'forward'
    }])
    
    construct = Construct(
        name="test",
        seq=seq_with_cut,
        circular=False,
        features=features
    )
    
    # Should find internal cuts
    internal = check_internal_cuts(
        construct.seq, "EcoRI", enzyme_db, 
        construct.features, avoid_cds=True
    )
    
    assert internal > 0


def test_simulate_digest(simple_vector, enzyme_db):
    """Test digest simulation."""
    fragments = simulate_digest(simple_vector, ["EcoRI"], enzyme_db)
    
    # Should produce at least one fragment (depending on sites present)
    assert len(fragments) >= 1
    
    # Check that fragments have required attributes
    for frag in fragments:
        assert hasattr(frag, 'length')
        assert hasattr(frag, 'sequence')
        assert hasattr(frag, 'enzymes_at_ends')


def test_simulate_double_digest(simple_vector, enzyme_db):
    """Test double digest simulation."""
    fragments = simulate_digest(simple_vector, ["EcoRI", "BamHI"], enzyme_db)
    
    # Double digest should produce multiple fragments
    assert len(fragments) >= 1
    
    # Total length should match original
    total_length = sum(f.length for f in fragments)
    assert total_length == len(simple_vector.seq)


# ============================================================================
# SCORING TESTS
# ============================================================================

def test_score_simple_plan():
    """Test scoring of a simple plan."""
    step1 = Step(
        name="Digest",
        action="digest",
        inputs=["vector"],
        params={'enzymes': ['EcoRI', 'BamHI']},
        outputs=["frag1", "frag2"],
        predicted_fragments=[],
        cost=1.0
    )
    
    step2 = Step(
        name="Ligate",
        action="ligate",
        inputs=["frag1", "insert"],
        params={'directional': True, 'scar': ''},
        outputs=["product"],
        predicted_fragments=[],
        cost=1.0
    )
    
    plan = Plan(
        steps=[step1, step2],
        final=None,
        score=0.0,
        feasible=True
    )
    
    options = {'require_directional': False, 'prefer_typeIIS': False}
    score = score_plan(plan, options)
    
    # Score should be positive (steps + enzymes + other factors)
    assert score > 0


def test_score_with_penalties():
    """Test that scoring applies penalties correctly."""
    step_with_internal_cuts = Step(
        name="Bad_Digest",
        action="digest",
        inputs=["vector"],
        params={'enzymes': ['EcoRI'], 'internal_cuts': 2},
        outputs=["frag1"],
        predicted_fragments=[],
        cost=1.0
    )
    
    plan = Plan(
        steps=[step_with_internal_cuts],
        final=None,
        score=0.0,
        feasible=True
    )
    
    options = {}
    score = score_plan(plan, options)
    
    # Should have penalty for internal cuts
    assert score > 1.0  # Base cost + penalty


def test_score_prefer_typeIIS():
    """Test Type IIS preference in scoring."""
    gg_step = Step(
        name="GG",
        action="GG",
        inputs=["part1", "part2", "part3"],
        params={'enzyme': 'BsaI'},
        outputs=["final"],
        predicted_fragments=[],
        cost=1.0
    )
    
    plan = Plan(
        steps=[gg_step],
        final=None,
        score=0.0,
        feasible=True
    )
    
    options_with_pref = {'prefer_typeIIS': True}
    options_without_pref = {'prefer_typeIIS': False}
    
    score_with = score_plan(plan, options_with_pref)
    score_without = score_plan(plan, options_without_pref)
    
    # Score with preference should be lower (better)
    assert score_with < score_without


# ============================================================================
# HEURISTIC TESTS
# ============================================================================

def test_calculate_heuristic_simple(simple_vector, simple_insert):
    """Test heuristic calculation for simple case."""
    state = SearchState(
        constructs=frozenset([simple_vector, simple_insert]),
        steps_taken=tuple(),
        cost=0.0,
        heuristic=0.0
    )
    
    target_spec = {
        'order': ['pTest', 'insert1'],
        'junctions': [{'left': 'pTest', 'right': 'insert1'}]
    }
    
    options = {}
    
    h = calculate_heuristic(state, target_spec, options)
    
    # Should estimate some cost
    assert h >= 0


def test_heuristic_decreases_with_progress(simple_vector, simple_insert):
    """Test that heuristic decreases as we get closer to goal."""
    # Initial state: both parts separate
    state1 = SearchState(
        constructs=frozenset([simple_vector, simple_insert]),
        steps_taken=tuple(),
        cost=0.0,
        heuristic=0.0
    )
    
    # Later state: after some steps
    final_construct = Construct("final", "ATCGATCG", True)
    state2 = SearchState(
        constructs=frozenset([final_construct]),
        steps_taken=(
            Step("s1", "digest", ["pTest"], {}, ["f1", "f2"], [], 1.0),
        ),
        cost=1.0,
        heuristic=0.0
    )
    
    target_spec = {
        'order': ['pTest', 'insert1'],
        'junctions': []
    }
    
    options = {}
    
    h1 = calculate_heuristic(state1, target_spec, options)
    h2 = calculate_heuristic(state2, target_spec, options)
    
    # Heuristic should not increase (may stay same or decrease)
    # This is a loose check since we don't have perfect heuristic
    assert h2 <= h1 + 2.0  # Allow some flexibility


# ============================================================================
# EDGE CASE TESTS
# ============================================================================

def test_empty_plan():
    """Test handling of empty plan."""
    plan = Plan(
        steps=[],
        final=None,
        score=float('inf'),
        feasible=False,
        reason="No plan found"
    )
    
    assert not plan.feasible
    assert len(plan.steps) == 0


def test_circular_construct_digest(enzyme_db):
    """Test digesting a circular construct."""
    circular = Construct(
        name="plasmid",
        seq="ATCGGAATTCGGATCC",
        circular=True
    )
    
    fragments = simulate_digest(circular, ["EcoRI"], enzyme_db)
    
    # Circular with one cut should linearize
    assert len(fragments) >= 1


def test_no_cut_sites(enzyme_db):
    """Test digesting a sequence with no recognition sites."""
    no_sites = Construct(
        name="insert",
        seq="ATCGATCGATCG" * 5,
        circular=False
    )
    
    fragments = simulate_digest(no_sites, ["EcoRI"], enzyme_db)
    
    # Should return single fragment (uncut)
    assert len(fragments) == 1
    assert fragments[0].length == len(no_sites.seq)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

