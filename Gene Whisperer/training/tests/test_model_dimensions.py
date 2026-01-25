"""
Test that model properly handles expanded feature dimensions.
"""
import sys

import torch

sys.path.insert(0, ".")


def test_model_with_288_dim_features():
    """Test model works with 288-dim engineered features."""
    from model_legacy import GeneWhispererStage1Legacy

    model = GeneWhispererStage1Legacy(
        vocab_size=4099,
        kmer=6,
        embedding_dim=256,
        num_layers=8,
        num_heads=8,
        ff_dim=1024,
        dropout=0.12,
        engineered_dim=288,  # NEW: 288 instead of 128
        use_engineered_features=True,
        use_tcn=True,
        tcn_hidden=256,
        tcn_levels=5,
        post_cnn_transformer_layers=4,
        engineered_mlp_hidden=512,  # NEW: 512 instead of 256
        engineered_mlp_output=256,  # NEW: 256 instead of 128
        fusion_hidden=384,  # NEW: 384 instead of 256
    )

    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    print(f"Total parameters: {total_params:,}")
    print(f"Trainable parameters: {trainable_params:,}")

    # Test forward pass
    batch_size = 4
    seq_len = 76  # 81 - 6 + 1 for 6-mer

    tokens = torch.randint(0, 4096, (batch_size, seq_len))
    engineered = torch.randn(batch_size, 288)  # NEW dimension

    model.eval()
    with torch.no_grad():
        output = model(tokens, engineered)

    assert output.shape == (batch_size, 1), f"Expected shape ({batch_size}, 1), got {output.shape}"

    print("\nInput shapes:")
    print(f"  Tokens: {tokens.shape}")
    print(f"  Engineered: {engineered.shape}")
    print(f"Output shape: {output.shape}")

    print("\n✓ TEST PASSED: Model works with 288-dim features")
    return True


def test_model_backward_compatible():
    """Test model still works with old 128-dim features."""
    from model_legacy import GeneWhispererStage1Legacy

    model = GeneWhispererStage1Legacy(
        vocab_size=67,
        kmer=3,
        embedding_dim=256,
        num_layers=8,
        num_heads=8,
        ff_dim=1024,
        engineered_dim=128,  # OLD dimension
        use_engineered_features=True,
        use_tcn=True,
        tcn_hidden=256,
        tcn_levels=4,
        post_cnn_transformer_layers=3,
        engineered_mlp_hidden=256,  # OLD
        engineered_mlp_output=128,  # OLD
        fusion_hidden=256,  # OLD
    )

    batch_size = 4
    seq_len = 79

    tokens = torch.randint(0, 64, (batch_size, seq_len))
    engineered = torch.randn(batch_size, 128)

    model.eval()
    with torch.no_grad():
        output = model(tokens, engineered)

    assert output.shape == (batch_size, 1)

    print("✓ TEST PASSED: Model backward compatible with 128-dim features")
    return True


def test_model_component_dimensions():
    """Test internal component dimensions are correct."""
    from model import EngineeredFeatureMLP
    from model_legacy import GatedAttentionFusion

    # Test EngineeredFeatureMLP with new dimensions
    mlp = EngineeredFeatureMLP(
        input_dim=288,
        hidden_dim=512,
        output_dim=256,
    )

    x = torch.randn(4, 288)
    out = mlp(x)
    assert out.shape == (4, 256), f"MLP output shape wrong: {out.shape}"
    print(f"EngineeredFeatureMLP: {x.shape} -> {out.shape} ✓")

    # Test GatedAttentionFusion with new dimensions
    fusion = GatedAttentionFusion(
        seq_dim=256,  # TCN hidden
        eng_dim=256,  # MLP output
        hidden_dim=384,
    )

    seq_feat = torch.randn(4, 256)
    eng_feat = torch.randn(4, 256)
    out = fusion(seq_feat, eng_feat)
    assert out.shape == (4, 384), f"Fusion output shape wrong: {out.shape}"
    print(f"GatedAttentionFusion: seq{seq_feat.shape} + eng{eng_feat.shape} -> {out.shape} ✓")

    print("\n✓ TEST PASSED: All component dimensions correct")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT 5 VALIDATION: Model Architecture Dimensions")
    print("=" * 60)

    test_model_with_288_dim_features()
    print()
    test_model_backward_compatible()
    print()
    test_model_component_dimensions()

    print("\n" + "=" * 60)
    print("ALL PROMPT 5 TESTS PASSED ✓")
    print("=" * 60)
