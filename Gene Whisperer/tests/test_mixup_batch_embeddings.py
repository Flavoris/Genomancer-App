"""Unit test for mixup_batch_embeddings: shape, padding semantics, and forward pass."""

import sys
from pathlib import Path

import numpy as np
import torch

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model import GeneWhispererStage1Legacy
from train_stage1 import mixup_batch_embeddings


def test_mixup_batch_embeddings():
    """Test mixup_batch_embeddings output tensors and padding behavior."""
    # Seed for determinism
    np.random.seed(0)
    torch.manual_seed(0)

    # Build small Stage 1 model as in parity test (dropout=0.0)
    model = GeneWhispererStage1Legacy(
        vocab_size=50,
        kmer=6,
        embedding_dim=16,
        num_layers=1,
        num_heads=4,
        ff_dim=32,
        dropout=0.0,
        pad_token_id=0,
        engineered_dim=8,
        use_engineered_features=True,
        use_tcn=False,
        post_cnn_transformer_layers=1,
    )
    model.eval()

    # Create tokens with padding zeros in different places for different samples
    B, L = 4, 12
    tokens = torch.tensor([
        [1, 5, 10, 20, 30, 40, 15, 25, 35, 45, 8, 3],           # no padding
        [0, 0, 7, 14, 21, 28, 35, 42, 49, 12, 6, 2],            # 2 padding at start
        [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 0, 0],           # 2 padding at end
        [9, 18, 27, 0, 0, 0, 11, 22, 33, 44, 5, 6],             # 3 padding in middle
    ], dtype=torch.long)

    engineered = torch.randn(B, 8)
    labels = torch.tensor([[0.0], [1.0], [1.0], [0.0]])

    D = 16  # embedding_dim

    # Call mixup_batch_embeddings
    mixed_emb, mixed_engineered, labels_a, labels_b, lam, padding_mask = mixup_batch_embeddings(
        model, tokens, engineered, labels, alpha=0.4
    )

    # -------------------------------------------------------------------------
    # Assert shapes
    # -------------------------------------------------------------------------
    assert mixed_emb.shape == (B, L, D), (
        f"mixed_emb.shape: expected ({B}, {L}, {D}), got {mixed_emb.shape}"
    )
    print(f"[OK] mixed_emb.shape == ({B}, {L}, {D})")

    assert mixed_engineered.shape == engineered.shape, (
        f"mixed_engineered.shape: expected {engineered.shape}, got {mixed_engineered.shape}"
    )
    print(f"[OK] mixed_engineered.shape == {engineered.shape}")

    assert labels_a.shape == labels.shape, (
        f"labels_a.shape: expected {labels.shape}, got {labels_a.shape}"
    )
    print(f"[OK] labels_a.shape == {labels.shape}")

    assert labels_b.shape == labels.shape, (
        f"labels_b.shape: expected {labels.shape}, got {labels_b.shape}"
    )
    print(f"[OK] labels_b.shape == {labels.shape}")

    # -------------------------------------------------------------------------
    # Assert lam in (0, 1)
    # -------------------------------------------------------------------------
    assert 0.0 < lam < 1.0, f"lam should be in (0, 1), got {lam}"
    print(f"[OK] 0.0 < lam ({lam:.4f}) < 1.0")

    # -------------------------------------------------------------------------
    # Assert padding_mask semantics
    # -------------------------------------------------------------------------
    if padding_mask is not None:
        assert padding_mask.dtype == torch.bool, (
            f"padding_mask dtype: expected bool, got {padding_mask.dtype}"
        )
        print(f"[OK] padding_mask.dtype == torch.bool")

        # Check that mixed_emb at fully-padded positions is all zeros
        if padding_mask.any():
            padded_embs = mixed_emb[padding_mask]
            max_abs = padded_embs.abs().max().item()
            assert max_abs < 1e-6, (
                f"mixed_emb at padding positions should be ~0, max abs = {max_abs}"
            )
            print(f"[OK] mixed_emb[padding_mask] is all zeros (max abs = {max_abs:.2e})")
        else:
            print("[INFO] No fully-padded positions (no overlap of pad tokens)")
    else:
        print("[INFO] padding_mask is None (no pad_token_id in model)")

    # -------------------------------------------------------------------------
    # Forward pass through forward_from_embeds
    # -------------------------------------------------------------------------
    with torch.no_grad():
        probs, logits = model.forward_from_embeds(
            mixed_emb, padding_mask, mixed_engineered, return_logits=True
        )

    # Check output shape: model outputs (B, 1) for binary classification
    assert probs.shape == (B, 1), f"probs.shape: expected ({B}, 1), got {probs.shape}"
    print(f"[OK] probs.shape == ({B}, 1)")

    # Sanity check: logits are finite and probs are valid probabilities
    assert torch.isfinite(logits).all(), "Logits contain NaN or Inf values"
    assert (probs >= 0.0).all() and (probs <= 1.0).all(), (
        f"probs out of [0, 1] range: {probs.squeeze().tolist()}"
    )
    print(f"[OK] probs in [0, 1] range: {probs.squeeze().tolist()}")

    print("\n=== All assertions passed ===")
    print(f"lam = {lam:.4f}")
    print(f"probs = {probs.squeeze().tolist()}")


if __name__ == "__main__":
    test_mixup_batch_embeddings()
