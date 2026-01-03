"""
Test CKSNAP (Composition of K-Spaced Nucleotide Acid Pairs) implementation.
"""
import torch
import numpy as np
import sys
sys.path.insert(0, '.')


def test_cksnap_output_shape():
    """Test CKSNAP produces correct output dimensions."""
    from dataset import compute_cksnap

    test_seq = "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCG"

    # Default: max_gap=5 → 16 pairs × 6 gaps = 96 dimensions
    result = compute_cksnap(test_seq, max_gap=5)

    expected_dim = 16 * 6  # 96
    assert result.shape == (expected_dim,), f"Expected shape ({expected_dim},), got {result.shape}"

    print(f"CKSNAP output shape: {result.shape}")
    print(f"Expected: ({expected_dim},)")
    print("\n✓ TEST PASSED: CKSNAP output shape correct")

    return True


def test_cksnap_normalization():
    """Test that CKSNAP frequencies sum to ~1 for each gap."""
    from dataset import compute_cksnap

    test_seq = "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCGATCGATCG" * 3
    result = compute_cksnap(test_seq, max_gap=5)

    # Each gap should have frequencies summing to ~1
    for gap in range(6):
        start = gap * 16
        end = start + 16
        gap_sum = result[start:end].sum().item()

        # Allow small tolerance for floating point
        assert 0.99 < gap_sum < 1.01, f"Gap {gap} sum = {gap_sum}, expected ~1.0"
        print(f"Gap {gap} sum: {gap_sum:.6f}")

    print("\n✓ TEST PASSED: CKSNAP frequencies properly normalized")
    return True


def test_cksnap_different_sequences():
    """Test that different sequences produce different features."""
    from dataset import compute_cksnap

    seq1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # All A's
    seq2 = "ATCGATCGATCGATCGATCGATCGATCG"  # Alternating
    seq3 = "GGGGCCCCTTTTAAAACCCCGGGGTTTT"  # Blocks

    feat1 = compute_cksnap(seq1)
    feat2 = compute_cksnap(seq2)
    feat3 = compute_cksnap(seq3)

    # All should be different
    diff_12 = (feat1 - feat2).abs().sum().item()
    diff_13 = (feat1 - feat3).abs().sum().item()
    diff_23 = (feat2 - feat3).abs().sum().item()

    assert diff_12 > 0.1, f"Sequences 1 and 2 too similar: diff={diff_12}"
    assert diff_13 > 0.1, f"Sequences 1 and 3 too similar: diff={diff_13}"
    assert diff_23 > 0.1, f"Sequences 2 and 3 too similar: diff={diff_23}"

    print(f"Diff(seq1, seq2): {diff_12:.4f}")
    print(f"Diff(seq1, seq3): {diff_13:.4f}")
    print(f"Diff(seq2, seq3): {diff_23:.4f}")

    # All-A sequence should have AA = 1.0 for gap=0
    aa_idx = 0  # AA is first in sorted dinucleotide list
    assert feat1[aa_idx] > 0.9, f"All-A sequence should have high AA frequency"
    print(f"All-A sequence AA frequency (gap=0): {feat1[aa_idx]:.4f}")

    print("\n✓ TEST PASSED: CKSNAP discriminates between sequences")
    return True


def test_cksnap_gap_logic():
    """Test that gap logic is correct."""
    from dataset import compute_cksnap

    # Sequence: ACGT
    # Gap=0: pairs at (0,1), (1,2), (2,3) → AC, CG, GT
    # Gap=1: pairs at (0,2), (1,3) → AG, CT
    # Gap=2: pairs at (0,3) → AT

    seq = "ACGT"
    result = compute_cksnap(seq, max_gap=2)

    # Dinucleotide order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
    # Indices:            0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15

    # Gap 0 (indices 0-15): AC=1, CG=6, GT=11 → each should be 1/3
    gap0 = result[0:16]
    assert gap0[1] > 0.3, f"AC should be present at gap=0, got {gap0[1]}"   # AC
    assert gap0[6] > 0.3, f"CG should be present at gap=0, got {gap0[6]}"   # CG
    assert gap0[11] > 0.3, f"GT should be present at gap=0, got {gap0[11]}" # GT

    # Gap 1 (indices 16-31): AG=2, CT=7 → each should be 1/2
    gap1 = result[16:32]
    assert gap1[2] > 0.4, f"AG should be present at gap=1, got {gap1[2]}"   # AG
    assert gap1[7] > 0.4, f"CT should be present at gap=1, got {gap1[7]}"   # CT

    print("Gap 0 AC frequency:", gap0[1].item())
    print("Gap 1 AG frequency:", gap1[2].item())

    print("\n✓ TEST PASSED: CKSNAP gap logic correct")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT 2 VALIDATION: CKSNAP Feature Encoding")
    print("=" * 60)

    test_cksnap_output_shape()
    print()
    test_cksnap_normalization()
    print()
    test_cksnap_different_sequences()
    print()
    test_cksnap_gap_logic()

    print("\n" + "=" * 60)
    print("ALL PROMPT 2 CKSNAP TESTS PASSED ✓")
    print("=" * 60)
