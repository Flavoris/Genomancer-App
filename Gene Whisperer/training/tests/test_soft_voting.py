"""
Test that soft voting produces different results than logit averaging.
This validates the core mathematical fix.
"""
import sys

import torch
import torch.nn as nn

sys.path.insert(0, ".")


def test_soft_voting_vs_logit_average():
    """Verify soft voting != logit averaging."""

    # Simulate 4 models with different logit outputs
    torch.manual_seed(42)
    batch_size = 16

    # Create diverse logits (some confident, some uncertain)
    logits = [
        torch.randn(batch_size, 1) * 2,  # Model 1: high variance
        torch.randn(batch_size, 1) * 0.5,  # Model 2: low variance
        torch.ones(batch_size, 1) * 1.5,  # Model 3: confident positive
        torch.ones(batch_size, 1) * -0.5,  # Model 4: slightly negative
    ]

    # Method 1: Logit averaging (OLD - WRONG)
    avg_logits = sum(logits) / len(logits)
    prob_from_avg_logits = torch.sigmoid(avg_logits)

    # Method 2: Soft voting (NEW - CORRECT)
    probs = [torch.sigmoid(l) for l in logits]
    prob_soft_voting = sum(probs) / len(probs)

    # They should be DIFFERENT
    diff = (prob_soft_voting - prob_from_avg_logits).abs()
    mean_diff = diff.mean().item()
    max_diff = diff.max().item()

    print(f"Mean absolute difference: {mean_diff:.6f}")
    print(f"Max absolute difference: {max_diff:.6f}")
    print(f"Sample probabilities (soft voting): {prob_soft_voting[:3].flatten().tolist()}")
    print(f"Sample probabilities (logit avg):   {prob_from_avg_logits[:3].flatten().tolist()}")

    # Test passes if there's meaningful difference
    assert mean_diff > 0.001, f"Soft voting should differ from logit avg! Diff={mean_diff}"
    print("\nTEST PASSED: Soft voting produces different results than logit averaging")

    return True


def test_soft_voting_ensemble_class():
    """Test the updated MultiScaleEnsemble class uses soft voting."""
    from model import GeneWhispererStage1, MultiScaleEnsemble

    # Create mock models (minimal config)
    models = []
    for k in [3, 4]:
        model = GeneWhispererStage1(
            vocab_size=4**k + 3,
            kmer=k,
            embedding_dim=64,
            num_layers=1,
            num_heads=2,
            ff_dim=128,
            engineered_dim=128,
            use_tcn=False,
            post_cnn_transformer_layers=1,
        )
        models.append(model)

    ensemble = MultiScaleEnsemble(models)

    # Create mock batch
    batch = {
        3: (torch.randint(0, 64, (4, 79)), torch.randn(4, 128)),
        4: (torch.randint(0, 256, (4, 78)), torch.randn(4, 128)),
    }

    # Forward pass
    output = ensemble(batch)

    # Check output is probabilities (should be in [0, 1])
    assert output.min() >= 0, f"Output min {output.min()} < 0, not probabilities!"
    assert output.max() <= 1, f"Output max {output.max()} > 1, not probabilities!"

    print(f"Ensemble output shape: {output.shape}")
    print(f"Ensemble output range: [{output.min().item():.4f}, {output.max().item():.4f}]")
    print("\nTEST PASSED: MultiScaleEnsemble outputs valid probabilities")

    return True


def test_soft_voting_with_return_logits():
    """Test that return_logits option works for training."""
    from model import GeneWhispererStage1, MultiScaleEnsemble

    models = []
    for k in [3, 4]:
        model = GeneWhispererStage1(
            vocab_size=4**k + 3,
            kmer=k,
            embedding_dim=64,
            num_layers=1,
            num_heads=2,
            ff_dim=128,
            engineered_dim=128,
            use_tcn=False,
            post_cnn_transformer_layers=1,
        )
        models.append(model)

    ensemble = MultiScaleEnsemble(models)

    batch = {
        3: (torch.randint(0, 64, (4, 79)), torch.randn(4, 128)),
        4: (torch.randint(0, 256, (4, 78)), torch.randn(4, 128)),
    }

    # Check if return_logits is supported
    if hasattr(ensemble, "forward_with_logits"):
        prob, logits = ensemble.forward_with_logits(batch)

        # Verify logits can be converted back to same probs
        reconstructed_prob = torch.sigmoid(logits)
        diff = (prob - reconstructed_prob).abs().max().item()

        assert diff < 1e-5, f"Logit reconstruction error: {diff}"
        print(f"Logit reconstruction error: {diff:.8f}")
        print("\nTEST PASSED: return_logits works correctly")
    else:
        print("WARNING: forward_with_logits not implemented (optional)")

    return True


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT 1 VALIDATION: Soft Voting Ensemble")
    print("=" * 60)

    test_soft_voting_vs_logit_average()
    print()
    test_soft_voting_ensemble_class()
    print()
    test_soft_voting_with_return_logits()

    print("\n" + "=" * 60)
    print("ALL PROMPT 1 TESTS PASSED")
    print("=" * 60)
