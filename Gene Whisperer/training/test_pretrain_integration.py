"""Test that pretraining integration is correct."""

import yaml
from pathlib import Path

def test_integration():
    print("Testing Pretraining Integration")
    print("=" * 50)
    
    # Test 1: Import early stopping
    try:
        from early_stopping import PretrainingEarlyStopping
        print("✓ early_stopping module imports correctly")
    except ImportError as e:
        print(f"✗ Failed to import early_stopping: {e}")
        return False
    
    # Test 2: Import pretrain_mlm functions
    try:
        from pretrain_mlm import get_kmer_pretrain_config
        print("✓ get_kmer_pretrain_config function exists")
    except ImportError as e:
        print(f"✗ Failed to import get_kmer_pretrain_config: {e}")
        return False
    
    # Test 3: Test config generation
    try:
        with open("config.yaml") as f:
            cfg = yaml.safe_load(f)
        
        for k in [3, 4, 5, 6]:
            kmer_cfg = get_kmer_pretrain_config(cfg, k)
            expected_vocab = (4 ** k) + 3
            actual_vocab = kmer_cfg.get("vocab_size")
            
            if actual_vocab != expected_vocab:
                print(f"✗ k={k}: Wrong vocab size (expected {expected_vocab}, got {actual_vocab})")
                return False
            
            print(f"✓ k={k}: vocab_size={actual_vocab}, path={kmer_cfg.get('mlm_vocab_path')}")
    except Exception as e:
        print(f"✗ Config generation failed: {e}")
        return False
    
    # Test 4: Check pretrain_all_kmers exists
    try:
        from pretrain_mlm import pretrain_all_kmers
        print("✓ pretrain_all_kmers function exists")
    except ImportError:
        print("⚠ pretrain_all_kmers not found (optional)")
    
    print("=" * 50)
    print("INTEGRATION TEST PASSED")
    return True

if __name__ == "__main__":
    test_integration()
