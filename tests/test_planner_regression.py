#!/usr/bin/env python3
"""
Regression test for the multi-step cloning planner.

This test verifies that the planner runs without crashing with a minimal spec,
particularly ensuring no "unhashable type: 'Step'" errors occur.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
from planner import plan_from_spec
from sim import load_enzyme_database


def test_planner_basic_run():
    """
    Test that the planner can run without crashing.
    
    This is a regression test for the fix to the "unhashable type: Step" error.
    The planner should be able to:
    1. Create SearchState objects
    2. Hash them for visited state tracking
    3. Run beam_search without TypeError
    """
    # Load enzyme database using the proper function
    enzyme_db = load_enzyme_database()
    
    # Load minimal spec
    spec_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "test_planner_spec.json"
    )
    
    with open(spec_path, 'r') as f:
        spec = json.load(f)
    
    # Run planner
    print("Testing planner with minimal spec...")
    print("=" * 60)
    
    try:
        plan = plan_from_spec(
            spec=spec,
            enzyme_db=enzyme_db,
            max_steps=2,
            options={'beam_width': 5}
        )
        
        print(f"✓ Planner completed without crash")
        print(f"  Feasible: {plan.feasible}")
        print(f"  Steps: {len(plan.steps)}")
        print(f"  Score: {plan.score:.2f}")
        
        if not plan.feasible:
            print(f"  Reason: {plan.reason}")
        
        print()
        print("=" * 60)
        print("SUCCESS: No TypeError occurred")
        print("The 'unhashable type: Step' issue has been fixed!")
        
    except TypeError as e:
        if "unhashable type" in str(e):
            print(f"✗ FAILED: {e}")
            print()
            print("The 'unhashable type: Step' error still occurs!")
            raise
        else:
            raise


def test_planner_state_hashing():
    """
    Test that SearchState objects can be hashed and used in sets.
    """
    from planner import SearchState, Construct
    
    print("\nTesting SearchState hashing...")
    
    # Create test constructs
    construct1 = Construct(
        name="test1",
        seq="ATCGATCG",
        circular=False
    )
    
    construct2 = Construct(
        name="test2",
        seq="GCTAGCTA",
        circular=True
    )
    
    # Create states
    state1 = SearchState(
        constructs=frozenset([construct1]),
        steps_taken=tuple(),
        cost=0.0,
        heuristic=0.0
    )
    
    state2 = SearchState(
        constructs=frozenset([construct1, construct2]),
        steps_taken=tuple(),
        cost=1.0,
        heuristic=0.5
    )
    
    # Test hashing
    try:
        hash1 = hash(state1)
        hash2 = hash(state2)
        print(f"  State 1 hash: {hash1}")
        print(f"  State 2 hash: {hash2}")
        
        # Test set membership
        visited = set()
        visited.add(state1)
        visited.add(state2)
        
        print(f"  ✓ Successfully added {len(visited)} states to visited set")
        
        # Test signature equality
        state1_copy = SearchState(
            constructs=frozenset([construct1]),
            steps_taken=tuple(),  # Different path
            cost=10.0,  # Different cost
            heuristic=5.0  # Different heuristic
        )
        
        assert state1.signature() == state1_copy.signature(), \
            "States with same constructs should have same signature"
        
        print(f"  ✓ Signature-based equality works correctly")
        
    except TypeError as e:
        print(f"  ✗ FAILED: {e}")
        raise


if __name__ == "__main__":
    print("=" * 60)
    print("PLANNER REGRESSION TEST")
    print("=" * 60)
    print()
    
    try:
        # Test state hashing
        test_planner_state_hashing()
        print()
        
        # Test basic planner run
        test_planner_basic_run()
        
        print()
        print("=" * 60)
        print("ALL TESTS PASSED ✓")
        print("=" * 60)
        sys.exit(0)
        
    except Exception as e:
        print()
        print("=" * 60)
        print(f"TEST FAILED: {e}")
        print("=" * 60)
        import traceback
        traceback.print_exc()
        sys.exit(1)

