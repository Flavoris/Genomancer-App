"""Smoke tests for GeneWhispererStage2 forward pass.

Verifies that the Stage 2 model:
1. Runs forward pass without errors
2. Produces output shape (B, 1)
3. Produces finite values (no NaN/Inf)
4. Works with all head types: transformer, bilstm, both

Usage:
    python tests/test_stage2_forward.py

    # Or with pytest
    pytest tests/test_stage2_forward.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import torch
import torch.nn as nn

from model import (
    GeneWhispererStage2,
    StrengthTransformerHead,
    StrengthBiLSTMHead,
    CombinedStrengthHead,
)


def test_strength_transformer_head():
    """Test StrengthTransformerHead forward pass."""
    batch_size = 4
    seq_len = 79
    input_dim = 256

    head = StrengthTransformerHead(
        input_dim=input_dim,
        n_layers=2,
        n_heads=4,
        ff_dim=512,
        dropout=0.1,
    )

    # Random token-level features
    x = torch.randn(batch_size, seq_len, input_dim)

    # Forward pass
    output = head(x)

    # Check output shape
    assert output.shape == (batch_size, input_dim), f"Expected ({batch_size}, {input_dim}), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    print(f"StrengthTransformerHead: input {x.shape} -> output {output.shape}")


def test_strength_bilstm_head():
    """Test StrengthBiLSTMHead forward pass."""
    batch_size = 4
    seq_len = 79
    input_dim = 256
    hidden_dim = 128

    head = StrengthBiLSTMHead(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        num_layers=1,
        dropout=0.1,
    )

    # Random token-level features
    x = torch.randn(batch_size, seq_len, input_dim)

    # Forward pass
    output = head(x)

    # Check output shape (BiLSTM outputs 2 * hidden_dim)
    expected_dim = 2 * hidden_dim
    assert output.shape == (batch_size, expected_dim), f"Expected ({batch_size}, {expected_dim}), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    print(f"StrengthBiLSTMHead: input {x.shape} -> output {output.shape}")


def test_combined_strength_head_concat():
    """Test CombinedStrengthHead with concatenation mode."""
    batch_size = 4
    seq_len = 79
    input_dim = 256
    bilstm_hidden = 128

    head = CombinedStrengthHead(
        input_dim=input_dim,
        transformer_layers=2,
        transformer_heads=4,
        bilstm_hidden=bilstm_hidden,
        bilstm_layers=1,
        dropout=0.1,
        combine_mode="concat",
        order="transformer_first",
    )

    # Random token-level features
    x = torch.randn(batch_size, seq_len, input_dim)

    # Forward pass
    output = head(x)

    # Check output shape (transformer_dim + 2*bilstm_hidden)
    expected_dim = input_dim + 2 * bilstm_hidden
    assert output.shape == (batch_size, expected_dim), f"Expected ({batch_size}, {expected_dim}), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    print(f"CombinedStrengthHead (concat): input {x.shape} -> output {output.shape}")


def test_combined_strength_head_avg():
    """Test CombinedStrengthHead with averaging mode."""
    batch_size = 4
    seq_len = 79
    input_dim = 256

    head = CombinedStrengthHead(
        input_dim=input_dim,
        transformer_layers=2,
        transformer_heads=4,
        bilstm_hidden=128,
        bilstm_layers=1,
        dropout=0.1,
        combine_mode="avg",
        order="transformer_first",
    )

    # Random token-level features
    x = torch.randn(batch_size, seq_len, input_dim)

    # Forward pass
    output = head(x)

    # Check output shape (avg mode projects to input_dim)
    assert output.shape == (batch_size, input_dim), f"Expected ({batch_size}, {input_dim}), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    print(f"CombinedStrengthHead (avg): input {x.shape} -> output {output.shape}")


def test_stage2_transformer_head():
    """Test GeneWhispererStage2 with transformer head."""
    batch_size = 4
    seq_len = 79  # (81bp - kmer + 1) for k=3
    vocab_size = 67
    engineered_dim = 128

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=192,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=engineered_dim,
        use_engineered_features=True,
        stage2_head_type="transformer",
        stage2_transformer_layers=2,
        stage2_transformer_heads=4,
    )

    # Random inputs
    tokens = torch.randint(1, vocab_size, (batch_size, seq_len))
    engineered = torch.randn(batch_size, engineered_dim)

    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(tokens, engineered)

    # Check output shape
    assert output.shape == (batch_size, 1), f"Expected ({batch_size}, 1), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    # Check probability range (model returns logits)
    probs = torch.sigmoid(output)
    assert (probs >= 0).all() and (probs <= 1).all(), "Sigmoid output should be in [0, 1] range"

    print(f"GeneWhispererStage2 (transformer): tokens {tokens.shape}, eng {engineered.shape} -> output {output.shape}")


def test_stage2_bilstm_head():
    """Test GeneWhispererStage2 with BiLSTM head."""
    batch_size = 4
    seq_len = 79
    vocab_size = 67
    engineered_dim = 128

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=192,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=engineered_dim,
        use_engineered_features=True,
        stage2_head_type="bilstm",
        stage2_bilstm_hidden=128,
        stage2_bilstm_layers=1,
    )

    # Random inputs
    tokens = torch.randint(1, vocab_size, (batch_size, seq_len))
    engineered = torch.randn(batch_size, engineered_dim)

    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(tokens, engineered)

    # Check output shape
    assert output.shape == (batch_size, 1), f"Expected ({batch_size}, 1), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    # Check probability range (model returns logits)
    probs = torch.sigmoid(output)
    assert (probs >= 0).all() and (probs <= 1).all(), "Sigmoid output should be in [0, 1] range"

    print(f"GeneWhispererStage2 (bilstm): tokens {tokens.shape}, eng {engineered.shape} -> output {output.shape}")


def test_stage2_both_heads():
    """Test GeneWhispererStage2 with both transformer and BiLSTM heads."""
    batch_size = 4
    seq_len = 79
    vocab_size = 67
    engineered_dim = 128

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=192,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=engineered_dim,
        use_engineered_features=True,
        stage2_head_type="both",
        stage2_transformer_layers=2,
        stage2_transformer_heads=4,
        stage2_bilstm_hidden=128,
        stage2_bilstm_layers=1,
        stage2_combine_mode="concat",
    )

    # Random inputs
    tokens = torch.randint(1, vocab_size, (batch_size, seq_len))
    engineered = torch.randn(batch_size, engineered_dim)

    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(tokens, engineered)

    # Check output shape
    assert output.shape == (batch_size, 1), f"Expected ({batch_size}, 1), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    # Check probability range (model returns logits)
    probs = torch.sigmoid(output)
    assert (probs >= 0).all() and (probs <= 1).all(), "Sigmoid output should be in [0, 1] range"

    print(f"GeneWhispererStage2 (both): tokens {tokens.shape}, eng {engineered.shape} -> output {output.shape}")


def test_stage2_without_engineered_features():
    """Test GeneWhispererStage2 without engineered features."""
    batch_size = 4
    seq_len = 79
    vocab_size = 67

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=192,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=0,
        use_engineered_features=False,
        stage2_head_type="transformer",
    )

    # Random inputs (no engineered features)
    tokens = torch.randint(1, vocab_size, (batch_size, seq_len))

    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(tokens, None)

    # Check output shape
    assert output.shape == (batch_size, 1), f"Expected ({batch_size}, 1), got {output.shape}"

    # Check finite values
    assert torch.isfinite(output).all(), "Output contains NaN or Inf values"

    print(f"GeneWhispererStage2 (no eng): tokens {tokens.shape} -> output {output.shape}")


def test_stage2_gradient_flow():
    """Test that gradients flow correctly through Stage 2 model."""
    batch_size = 4
    seq_len = 79
    vocab_size = 67
    engineered_dim = 128

    model = GeneWhispererStage2(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=192,
        num_layers=4,
        num_heads=4,
        ff_dim=384,
        dropout=0.1,
        pad_token_id=0,
        engineered_dim=engineered_dim,
        use_engineered_features=True,
        stage2_head_type="transformer",
    )

    # Random inputs
    tokens = torch.randint(1, vocab_size, (batch_size, seq_len))
    engineered = torch.randn(batch_size, engineered_dim)
    targets = torch.randint(0, 2, (batch_size, 1)).float()

    # Forward pass
    model.train()
    output = model(tokens, engineered)

    # Compute loss
    loss = nn.BCEWithLogitsLoss()(output, targets)

    # Backward pass
    loss.backward()

    # Check that some gradients are non-zero
    has_gradients = False
    for name, param in model.named_parameters():
        if param.grad is not None and param.grad.abs().sum() > 0:
            has_gradients = True
            break

    assert has_gradients, "No gradients found in any parameter"

    # Check that loss is finite
    assert torch.isfinite(torch.tensor(loss.item())), "Loss is not finite"

    print(f"Gradient flow test passed. Loss: {loss.item():.4f}")


def run_all_tests():
    """Run all smoke tests."""
    print("=" * 60)
    print("Running GeneWhispererStage2 Smoke Tests")
    print("=" * 60)

    tests = [
        ("StrengthTransformerHead", test_strength_transformer_head),
        ("StrengthBiLSTMHead", test_strength_bilstm_head),
        ("CombinedStrengthHead (concat)", test_combined_strength_head_concat),
        ("CombinedStrengthHead (avg)", test_combined_strength_head_avg),
        ("Stage2 with Transformer Head", test_stage2_transformer_head),
        ("Stage2 with BiLSTM Head", test_stage2_bilstm_head),
        ("Stage2 with Both Heads", test_stage2_both_heads),
        ("Stage2 without Engineered Features", test_stage2_without_engineered_features),
        ("Gradient Flow", test_stage2_gradient_flow),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            print(f"\n[TEST] {name}")
            test_func()
            print(f"[PASS] {name}")
            passed += 1
        except Exception as e:
            print(f"[FAIL] {name}: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
