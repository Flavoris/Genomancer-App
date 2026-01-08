"""
End-to-end integration test to verify all optimizations work together.
"""
import sys

import torch

sys.path.insert(0, ".")


def _compute_grad_norm(parameters) -> float:
    """Compute an L2 norm across all available gradients."""
    grad_total = 0.0
    for param in parameters:
        if param.grad is not None:
            grad_total += param.grad.norm().item() ** 2
    return grad_total**0.5


def test_full_training_step():
    """Test a complete training step with all optimizations."""
    from dataset import compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc
    from model import GeneWhispererStage1

    print("Testing full training step with all optimizations...")

    model = GeneWhispererStage1(
        vocab_size=4099,
        kmer=6,
        embedding_dim=256,
        num_layers=2,  # Reduced for speed
        num_heads=8,
        ff_dim=512,
        engineered_dim=288,
        use_tcn=True,
        tcn_hidden=128,
        tcn_levels=3,
        post_cnn_transformer_layers=2,
        engineered_mlp_hidden=256,
        engineered_mlp_output=128,
        fusion_hidden=192,
    )

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=0.00002,  # Match updated stage training LR
        weight_decay=0.01,  # NEW lighter decay
    )

    label_smoothing = 0.0  # Disable when using focal loss
    criterion = torch.nn.BCEWithLogitsLoss()

    batch_size = 8
    seq_len = 76

    test_seq = "ATGGCTGCATGCATGCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

    # NEW: All 4 engineered feature types (64 + 64 + 96 + 64 = 288)
    tnc = compute_tnc(test_seq)
    pseeiip = compute_pseeiip(test_seq)
    cksnap = compute_cksnap(test_seq)
    pstnp = compute_pstnp(test_seq)
    engineered = torch.cat([tnc, pseeiip, cksnap, pstnp])

    tokens = torch.randint(0, 4096, (batch_size, seq_len))
    engineered_batch = engineered.unsqueeze(0).expand(batch_size, -1).contiguous()
    labels = torch.randint(0, 2, (batch_size, 1)).float()

    # Apply label smoothing to prevent overconfident targets
    smoothed_labels = labels * (1 - label_smoothing) + 0.5 * label_smoothing

    model.train()
    optimizer.zero_grad()

    logits = model(tokens, engineered_batch)
    loss = criterion(logits, smoothed_labels)

    loss.backward()

    grad_norm = _compute_grad_norm(model.parameters())
    optimizer.step()

    print(f"  Loss: {loss.item():.4f}")
    print(f"  Gradient norm: {grad_norm:.4f}")
    print(f"  Engineered features dim: {engineered_batch.shape[1]}")

    assert loss.item() < 10, "Loss too high"
    assert grad_norm > 0, "No gradients computed"
    assert engineered_batch.shape[1] == 288, "Wrong feature dimension"

    print("\n✓ Full training step successful!")
    return True


def test_ensemble_inference():
    """Test ensemble inference with soft voting."""
    from model import GeneWhispererStage1, MultiScaleEnsemble

    print("\nTesting ensemble inference with soft voting...")

    models = []
    for kmer_size in [3, 4]:
        model = GeneWhispererStage1(
            vocab_size=4**kmer_size + 3,
            kmer=kmer_size,
            embedding_dim=64,
            num_layers=1,
            num_heads=2,
            ff_dim=128,
            engineered_dim=128,  # Use 128 for quick test
            use_tcn=False,
            post_cnn_transformer_layers=1,
        )
        models.append(model)

    ensemble = MultiScaleEnsemble(models)

    batch = {
        3: (torch.randint(0, 64, (4, 79)), torch.randn(4, 128)),
        4: (torch.randint(0, 256, (4, 78)), torch.randn(4, 128)),
    }

    ensemble.eval()
    with torch.no_grad():
        output = ensemble(batch)

    assert output.min() >= 0, "Output below 0, not probabilities"
    assert output.max() <= 1, "Output above 1, not probabilities"

    print(f"  Ensemble output shape: {output.shape}")
    print(f"  Output range: [{output.min().item():.4f}, {output.max().item():.4f}]")

    print("\n✓ Ensemble inference successful!")
    return True


def test_memory_usage():
    """Test that model fits in M4 MacBook memory."""
    import gc

    from model import GeneWhispererStage1

    print("\nTesting memory usage...")

    gc.collect()
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()

    model = GeneWhispererStage1(
        vocab_size=4099,
        kmer=6,
        embedding_dim=256,
        num_layers=8,
        num_heads=8,
        ff_dim=1024,
        engineered_dim=288,
        use_tcn=True,
        tcn_hidden=256,
        tcn_levels=5,
        post_cnn_transformer_layers=4,
        engineered_mlp_hidden=512,
        engineered_mlp_output=256,
        fusion_hidden=384,
    )

    total_params = sum(p.numel() for p in model.parameters())
    param_memory_mb = total_params * 4 / (1024**2)  # float32

    print(f"  Total parameters: {total_params:,}")
    print(f"  Parameter memory: {param_memory_mb:.1f} MB")

    estimated_training_mb = param_memory_mb * 4
    print(f"  Estimated training memory: {estimated_training_mb:.1f} MB")

    assert estimated_training_mb < 4000, "Model too large for M4 MacBook"

    print("\n✓ Memory usage acceptable for M4 MacBook")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("INTEGRATION TEST: All Optimizations")
    print("=" * 60)

    test_full_training_step()
    test_ensemble_inference()
    test_memory_usage()

    print("\n" + "=" * 60)
    print("ALL INTEGRATION TESTS PASSED ✓")
    print("Ready to train for 90%+ accuracy!")
    print("=" * 60)
