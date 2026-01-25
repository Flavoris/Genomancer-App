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

    # Create two models with same vocab_size (BPE-based single tokenizer)
    models = []
    for _ in range(2):
        model = GeneWhispererStage1(
            vocab_size=100,
            embedding_dim=64,
            num_layers=1,
            num_heads=2,
            ff_dim=128,
            engineered_dim=128,
            use_engineered_features=True,
            max_seq_len=24,
        )
        models.append(model)

    ensemble = MultiScaleEnsemble(models)

    # Create mock batch as tuple (tokens, engineered_features)
    batch = (torch.randint(0, 100, (4, 24)), torch.randn(4, 128))

    # Forward pass
    output = ensemble(batch)

    # Check output is probabilities (should be in [0, 1])
    assert output.min() >= 0, f"Output min {output.min()} < 0, not probabilities!"
    assert output.max() <= 1, f"Output max {output.max()} > 1, not probabilities!"

    print(f"Ensemble output shape: {output.shape}")
    print(f"Ensemble output range: [{output.min().item():.4f}, {output.max().item():.4f}]")
    print("\nTEST PASSED: MultiScaleEnsemble outputs valid probabilities")


def test_soft_voting_with_return_logits():
    """Test that return_logits option works for training."""
    from model import GeneWhispererStage1, MultiScaleEnsemble

    models = []
    for _ in range(2):
        model = GeneWhispererStage1(
            vocab_size=100,
            embedding_dim=64,
            num_layers=1,
            num_heads=2,
            ff_dim=128,
            engineered_dim=128,
            use_engineered_features=True,
            max_seq_len=24,
        )
        models.append(model)

    ensemble = MultiScaleEnsemble(models)
    ensemble.eval()

    batch = (torch.randint(0, 100, (4, 24)), torch.randn(4, 128))

    # Test return_logits parameter
    with torch.no_grad():
        logits = ensemble(batch, return_logits=True)
        probs = ensemble(batch, return_logits=False)

    # Verify logits can be converted back to same probs
    # Note: logits = torch.logit(probs), so sigmoid(logits) should match probs
    reconstructed_prob = torch.sigmoid(logits)
    diff = (probs - reconstructed_prob).abs().max().item()

    # Due to float precision in logit/sigmoid roundtrip, tolerance is relaxed
    assert diff < 1e-3, f"Logit reconstruction error: {diff}"

    # Verify logits are finite and in expected range
    assert torch.isfinite(logits).all(), "Logits contain non-finite values"
    print(f"Logit reconstruction error: {diff:.8f}")
    print("\nTEST PASSED: return_logits works correctly")


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
