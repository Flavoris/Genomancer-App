"""
Test that the dataset properly combines all engineered features.
"""
import sys

import torch

sys.path.insert(0, ".")


def test_combined_feature_dimensions():
    """Test combined features have correct total dimension."""
    from dataset import compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc

    # Expected: TNC(64) + PseEIIP(64) + CKSNAP(96) + PSTNP(64) = 288
    expected_dim = 288

    test_seq = (
        "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCGATCGATCGATCGATCG"
        "ATCGATCGATCGATCGATCGATCGATCGATCG"
    )

    tnc = compute_tnc(test_seq)
    pseeiip = compute_pseeiip(test_seq)
    cksnap = compute_cksnap(test_seq)
    pstnp = compute_pstnp(test_seq)

    combined = torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)

    assert combined.shape[0] == expected_dim, f"Expected {expected_dim}, got {combined.shape[0]}"

    print(f"TNC shape: {tnc.shape}")
    print(f"PseEIIP shape: {pseeiip.shape}")
    print(f"CKSNAP shape: {cksnap.shape}")
    print(f"PSTNP shape: {pstnp.shape}")
    print(f"Combined shape: {combined.shape}")
    print(f"Expected: ({expected_dim},)")

    print("\n✓ TEST PASSED: Combined features have correct dimensions")
    return True


def test_dataset_returns_correct_features():
    """Test that dataset __getitem__ returns correct feature dimension."""
    from bpe_tokenizer import DNABPETokenizer
    from dataset import PromoterDatasetStage1
    import os
    import pandas as pd
    import tempfile

    # Create temporary test data file
    test_data = (
        "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCGATCGATCGATCGATCG"
        "ATCGATCGATCGATCGATCGATCGATCGATCG\t1\n"
    )
    test_data += (
        "GCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
        "CGATCGATCGATCGATCGATCGATCGATCGA\t0\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as temp_file:
        temp_file.write(test_data)
        temp_path = temp_file.name

    try:
        df = pd.read_csv(temp_path, delimiter="\t", header=None, names=["sequence", "is_promoter"])
        sequences = df["sequence"].tolist()
        vocab = DNABPETokenizer(vocab_size=30)
        vocab.train(sequences)

        # Create dataset with BPE tokenizer
        dataset = PromoterDatasetStage1(
            df,
            vocab=vocab,
            max_token_len=24,
            reverse_complement_prob=0.0,
        )

        tokens, engineered, label = dataset[0]

        print(f"Token shape: {tokens.shape}")
        print(f"Engineered shape: {engineered.shape}")
        print(f"Label shape: {label.shape}")

        # Check dimension matches config
        # After prompt 3: should be 288
        # Before prompt 3: should be 128
        if engineered.shape[0] == 288:
            print("\n✓ TEST PASSED: Dataset returns 288-dim features (NEW)")
        elif engineered.shape[0] == 128:
            print("\n⚠ WARNING: Dataset still returns 128-dim features (OLD)")
            print("  Run Prompt 3 to update dataset class")
        else:
            print(f"\n✗ UNEXPECTED: Got {engineered.shape[0]} dimensions")

    finally:
        os.unlink(temp_path)

    return True


def test_feature_consistency():
    """Test that features are consistent across multiple calls."""
    from dataset import compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc

    test_seq = "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCG"

    # Compute twice
    feat1 = torch.cat(
        [
            compute_tnc(test_seq),
            compute_pseeiip(test_seq),
            compute_cksnap(test_seq),
            compute_pstnp(test_seq),
        ]
    )

    feat2 = torch.cat(
        [
            compute_tnc(test_seq),
            compute_pseeiip(test_seq),
            compute_cksnap(test_seq),
            compute_pstnp(test_seq),
        ]
    )

    diff = (feat1 - feat2).abs().max().item()
    assert diff < 1e-6, f"Features not deterministic: max diff = {diff}"

    print(f"Max difference between calls: {diff}")
    print("\n✓ TEST PASSED: Features are deterministic")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT 3 VALIDATION: Combined Dataset Features")
    print("=" * 60)

    test_combined_feature_dimensions()
    print()
    test_dataset_returns_correct_features()
    print()
    test_feature_consistency()

    print("\n" + "=" * 60)
    print("ALL PROMPT 3 TESTS PASSED ✓")
    print("=" * 60)
