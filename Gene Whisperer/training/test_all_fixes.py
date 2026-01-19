#!/usr/bin/env python3
"""
Comprehensive integration test for Gene Whisperer fixes.
Run this after applying all fixes to verify everything works.
"""

import sys
import yaml
import numpy as np
from pathlib import Path

# Add training directory to path
script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))

def test_config():
    """Test config.yaml has correct values."""
    print("\n" + "=" * 60)
    print("TEST 1: Configuration Values")
    print("=" * 60)

    with open(script_dir / 'config.yaml', 'r') as f:
        cfg = yaml.safe_load(f)

    simplified = cfg.get('simplified_model', {}) if isinstance(cfg.get('simplified_model'), dict) else {}
    tests = [
        ('max_bp_len', 81, cfg.get('max_bp_len')),
        ('mlm_window_size', 234, cfg.get('mlm_window_size')),
        ('stage1_lr', 0.00002, cfg.get('stage1_lr')),
        ('stage2_lr', 0.00001, cfg.get('stage2_lr')),
        ('mlm_lr', 0.0002, cfg.get('mlm_lr')),
        ('stage1_reverse_complement_prob', 0.5, cfg.get('stage1_reverse_complement_prob')),
        ('simplified_pooling_type', 'attention', simplified.get('pooling_type')),
        ('simplified_classifier_hidden', 256, simplified.get('classifier_hidden')),
        ('simplified_classifier_dropout', 0.15, simplified.get('classifier_dropout')),
        ('simplified_fusion_method', 'concat', simplified.get('fusion_method')),
    ]

    all_passed = True
    for name, expected, actual in tests:
        status = "✓" if actual == expected else "✗"
        if actual != expected:
            all_passed = False
        print(f"  {status} {name}: {actual} (expected {expected})")

    assert all_passed, "Config validation failed!"
    print("\n✓ All config tests passed!")
    return True


def test_pseeiip():
    """Test PseEIIP produces correct number of non-zero elements."""
    print("\n" + "=" * 60)
    print("TEST 2: PseEIIP Feature Engineering")
    print("=" * 60)

    from dataset import compute_pseeiip, compute_tnc, compute_cksnap, compute_pstnp

    # Test sequence (81bp like benchmark)
    test_seq = "TTGACAATTTTTCTTGATAATGTAACTCACTTAATCTTGATAAATGCTATAATGTGTCGAAAAAAAAAAAAAAAAAAAA"

    pseeiip = compute_pseeiip(test_seq)
    tnc = compute_tnc(test_seq)
    cksnap = compute_cksnap(test_seq)
    pstnp = compute_pstnp(test_seq)

    pseeiip_nonzero = (pseeiip != 0).sum().item()
    tnc_nonzero = (tnc != 0).sum().item()
    cksnap_nonzero = (cksnap != 0).sum().item()
    pstnp_nonzero = (pstnp != 0).sum().item()

    print(f"  TNC:     {tnc.shape[0]}D, {tnc_nonzero} non-zero")
    print(f"  PseEIIP: {pseeiip.shape[0]}D, {pseeiip_nonzero} non-zero")
    print(f"  CKSNAP:  {cksnap.shape[0]}D, {cksnap_nonzero} non-zero")
    print(f"  PSTNP:   {pstnp.shape[0]}D, {pstnp_nonzero} non-zero")

    # PseEIIP should have many non-zero elements, not just 3
    assert pseeiip_nonzero > 20, f"PseEIIP BROKEN: only {pseeiip_nonzero} non-zero (should be >20)"
    assert pseeiip.shape[0] == 64, f"PseEIIP wrong shape: {pseeiip.shape}"

    # Verify total engineered dim
    total_dim = tnc.shape[0] + pseeiip.shape[0] + cksnap.shape[0] + pstnp.shape[0]
    assert total_dim == 288, f"Total engineered dim should be 288, got {total_dim}"

    print(f"\n  Total engineered features: {total_dim}D")
    print("✓ PseEIIP fix verified!")
    return True


def test_dataset_creation():
    """Test that datasets can be created with correct dimensions."""
    print("\n" + "=" * 60)
    print("TEST 3: Dataset Creation")
    print("=" * 60)

    from dataset import PromoterDatasetStage1, KmerVocabulary
    import pandas as pd

    with open(script_dir / 'config.yaml', 'r') as f:
        cfg = yaml.safe_load(f)

    max_bp_len = cfg.get('max_bp_len', 81)
    engineered_dim = cfg.get('engineered_dim', 288)

    # Create a minimal test dataframe with realistic sequence
    test_seq = "TTGACAATTTTTCTTGATAATGTAACTCACTTAATCTTGATAAATGCTATAATGTGTCGAAAAAAAAAAAAAAAAAAAA"
    test_df = pd.DataFrame({
        'sequence': [test_seq] * 10,
        'is_promoter': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    })

    # Build vocab
    vocab = KmerVocabulary.build_from_sequences(test_df['sequence'].tolist(), k=6)

    # Create dataset
    dataset = PromoterDatasetStage1(
        test_df,
        max_bp_len=max_bp_len,
        vocab=vocab,
        use_engineered_features=True,
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.5,  # Test with augmentation
    )

    # Get a sample
    tokens, engineered, label = dataset[0]

    print(f"  Max BP length: {max_bp_len}")
    print(f"  Token shape: {tokens.shape}")
    print(f"  Engineered shape: {engineered.shape}")
    print(f"  Engineered non-zero: {(engineered != 0).sum().item()}")

    # Verify dimensions
    expected_tokens = max_bp_len - vocab.k + 1
    assert tokens.shape[0] == expected_tokens, f"Token length wrong: {tokens.shape[0]} vs {expected_tokens}"
    assert engineered.shape[0] == engineered_dim, f"Engineered dim wrong: {engineered.shape[0]} vs {engineered_dim}"

    # Engineered features should have many non-zero (not sparse due to broken PseEIIP)
    nonzero_ratio = (engineered != 0).sum().item() / engineered_dim
    assert nonzero_ratio > 0.3, f"Too few non-zero engineered features: {nonzero_ratio:.1%}"

    print("✓ Dataset creation verified!")
    return True


def test_model_creation():
    """Test that model can be created with correct architecture."""
    print("\n" + "=" * 60)
    print("TEST 4: Model Architecture")
    print("=" * 60)

    import torch
    from model import GeneWhispererStage1
    from dataset import KmerVocabulary

    with open(script_dir / 'config.yaml', 'r') as f:
        cfg = yaml.safe_load(f)

    # Build a test vocab
    vocab = KmerVocabulary.build_from_sequences(['ATGC' * 20], k=6)

    simplified = cfg.get('simplified_model', {}) if isinstance(cfg.get('simplified_model'), dict) else {}
    pooling_type = simplified.get('pooling_type', 'attention')
    use_attention_pool = pooling_type == 'attention'
    classifier_hidden = simplified.get('classifier_hidden', 256)
    classifier_dropout = simplified.get('classifier_dropout', 0.15)

    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=6,
        embedding_dim=cfg.get('embedding_dim', 384),
        num_layers=cfg.get('transformer_layers', 12),
        num_heads=cfg.get('transformer_heads', 12),
        ff_dim=cfg.get('transformer_ff_dim', 1536),
        dropout=cfg.get('transformer_dropout', 0.12),
        pad_token_id=vocab.pad_id,
        engineered_dim=cfg.get('engineered_dim', 288),
        use_engineered_features=True,
        use_attention_pool=use_attention_pool,
        engineered_mlp_hidden=cfg.get('engineered_mlp_hidden', 512),
        engineered_mlp_output=cfg.get('engineered_mlp_output', 384),
        classifier_hidden=classifier_hidden,
        classifier_dropout=classifier_dropout,
    )

    if use_attention_pool:
        assert model.pool is not None, "Expected attention pooling to be enabled"
    else:
        assert model.pool is None, "Expected mean pooling to be used"

    # Test forward pass
    batch_size = 2
    max_bp_len = cfg.get('max_bp_len', 81)
    seq_len = max_bp_len - 6 + 1  # k=6

    tokens = torch.randint(0, len(vocab.itos), (batch_size, seq_len))
    engineered = torch.randn(batch_size, cfg.get('engineered_dim', 288))

    with torch.no_grad():
        output = model(tokens, engineered)

    print(f"  Input tokens shape: {tokens.shape}")
    print(f"  Input engineered shape: {engineered.shape}")
    print(f"  Output shape: {output.shape}")

    assert output.shape == (batch_size, 1), f"Output shape wrong: {output.shape}"

    print("✓ Model architecture verified!")
    return True


def main():
    print("=" * 60)
    print("GENE WHISPERER COMPREHENSIVE FIX VERIFICATION")
    print("=" * 60)

    all_passed = True

    try:
        test_config()
    except Exception as e:
        print(f"✗ Config test failed: {e}")
        all_passed = False

    try:
        test_pseeiip()
    except Exception as e:
        print(f"✗ PseEIIP test failed: {e}")
        all_passed = False

    try:
        test_dataset_creation()
    except Exception as e:
        print(f"✗ Dataset test failed: {e}")
        all_passed = False

    try:
        test_model_creation()
    except Exception as e:
        print(f"✗ Model test failed: {e}")
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("ALL TESTS PASSED! ✓")
        print("=" * 60)
        print("""
Ready for retraining. Expected improvements:
- PseEIIP fix: +3-5% accuracy
- Sequence length fix: +1-2% accuracy
- Reverse complement: +1-2% accuracy
- Architecture tuning: +1-2% accuracy

Total expected: 84-85% → 90-94% accuracy
""")
    else:
        print("SOME TESTS FAILED! ✗")
        print("=" * 60)
        print("Please fix the failing tests before retraining.")
        sys.exit(1)


if __name__ == "__main__":
    main()
