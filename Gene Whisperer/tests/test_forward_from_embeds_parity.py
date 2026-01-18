"""Deterministic parity test: forward(tokens) == forward_from_embeds(embedding(tokens))."""

import torch

import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model import GeneWhispererStage1Legacy


def test_forward_from_embeds_parity():
    """Test that forward and forward_from_embeds produce identical outputs."""
    # Set seed for reproducibility
    torch.manual_seed(42)

    # Instantiate a tiny model
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

    # Set to eval mode (critical for determinism - disables dropout)
    model.eval()

    # Create test inputs
    # tokens: (B=4, L=12) ints in [0..49], include some zeros as pad
    tokens = torch.tensor([
        [1, 5, 10, 20, 30, 40, 15, 25, 35, 45, 8, 3],
        [0, 0, 7, 14, 21, 28, 35, 42, 49, 12, 6, 2],  # 2 padding tokens at start
        [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48],
        [9, 18, 27, 36, 45, 0, 0, 0, 11, 22, 33, 44],  # 3 padding tokens in middle
    ], dtype=torch.long)

    # engineered_features: (4, 8) floats
    engineered_features = torch.randn(4, 8)

    # Compute using forward()
    with torch.no_grad():
        probs1 = model(tokens, engineered_features)

    # Compute using forward_from_embeds()
    with torch.no_grad():
        embeds = model.embedding(tokens)
        padding_mask = tokens.eq(0)
        probs2 = model.forward_from_embeds(embeds, padding_mask, engineered_features)

    # Assert shapes match
    assert probs1.shape == probs2.shape, (
        f"Shape mismatch: forward={probs1.shape}, forward_from_embeds={probs2.shape}"
    )

    # Assert values match
    if not torch.allclose(probs1, probs2, atol=1e-6, rtol=1e-6):
        max_diff = (probs1 - probs2).abs().max().item()
        raise AssertionError(
            f"Output mismatch: max diff = {max_diff:.2e}\n"
            f"probs1 = {probs1.squeeze().tolist()}\n"
            f"probs2 = {probs2.squeeze().tolist()}"
        )

    print(f"Shape: {probs1.shape}")
    print(f"probs1: {probs1.squeeze().tolist()}")
    print(f"probs2: {probs2.squeeze().tolist()}")
    print(f"Max diff: {(probs1 - probs2).abs().max().item():.2e}")
    print("OK - forward and forward_from_embeds produce identical outputs")


if __name__ == "__main__":
    test_forward_from_embeds_parity()
