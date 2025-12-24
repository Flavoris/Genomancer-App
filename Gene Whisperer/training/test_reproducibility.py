"""Acceptance tests for deterministic reproducibility.

Tests verify:
1. Same seed → same golden batch loss within tolerance
2. Changing seed changes batch composition

Usage:
    # Run all tests
    python test_reproducibility.py
    
    # Run with pytest
    pytest test_reproducibility.py -v
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from typing import Callable, Tuple

import torch
import torch.nn as nn

# Add training directory to path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from seed_utils import (
    set_global_seed,
    get_deterministic_dataloader_kwargs,
    verify_reproducibility,
)


class TestSeedUtils(unittest.TestCase):
    """Tests for seed_utils module."""
    
    def test_verify_reproducibility(self):
        """Verify that seed utility correctly enables reproducibility."""
        # Test multiple seeds
        for seed in [42, 1337, 0, 999999]:
            self.assertTrue(
                verify_reproducibility(seed),
                f"Reproducibility check failed for seed {seed}"
            )
    
    def test_different_seeds_different_values(self):
        """Verify that different seeds produce different random values."""
        seed1, seed2 = 42, 43
        
        set_global_seed(seed1)
        vals1 = torch.rand(10).tolist()
        
        set_global_seed(seed2)
        vals2 = torch.rand(10).tolist()
        
        self.assertNotEqual(vals1, vals2, "Different seeds should produce different values")
    
    def test_same_seed_same_values(self):
        """Verify that same seed produces identical values."""
        seed = 42
        
        set_global_seed(seed)
        vals1 = torch.rand(100).tolist()
        
        set_global_seed(seed)
        vals2 = torch.rand(100).tolist()
        
        self.assertEqual(vals1, vals2, "Same seed should produce identical values")


class TestGoldenBatch(unittest.TestCase):
    """Tests for golden batch functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that require loading configs."""
        import yaml
        
        config_path = SCRIPT_DIR / "config.yaml"
        if not config_path.exists():
            raise unittest.SkipTest("config.yaml not found")
        
        with config_path.open("r", encoding="utf-8") as f:
            cls.cfg = yaml.safe_load(f) or {}
        
        # Check if required data files exist
        fasta_path = Path(cls.cfg.get("mlm_fasta_path", "")).expanduser()
        if not fasta_path.is_absolute():
            fasta_path = (SCRIPT_DIR / fasta_path).resolve()
        else:
            fasta_path = fasta_path.resolve()
        if not fasta_path.exists():
            raise unittest.SkipTest(f"FASTA file not found: {fasta_path}")
    
    def _create_model_and_loader(self, seed: int):
        """Helper to create model and dataloader with specific seed."""
        from pretrain_mlm import prepare_dataset, MLMDataset, MLMCollator, DNAMLM
        from model import DNAEncoder
        from torch.utils.data import DataLoader, random_split
        
        set_global_seed(seed)
        
        token_tensors, vocab = prepare_dataset(self.cfg)
        full_dataset = MLMDataset(token_tensors)
        
        val_ratio = float(self.cfg.get("mlm_val_ratio", 0.1))
        val_size = int(len(full_dataset) * val_ratio)
        train_size = len(full_dataset) - val_size
        train_dataset, _ = random_split(
            full_dataset,
            [train_size, val_size],
            generator=torch.Generator().manual_seed(seed)
        )
        
        batch_size = int(self.cfg.get("mlm_batch_size", 32))
        collator = MLMCollator(vocab, mask_prob=0.15, track_spans=True)
        deterministic_kwargs = get_deterministic_dataloader_kwargs(seed, num_workers=0)
        
        train_loader = DataLoader(
            train_dataset,
            batch_size=batch_size,
            shuffle=True,
            collate_fn=collator,
            num_workers=0,
            drop_last=True,
            **deterministic_kwargs,
        )
        
        # Build minimal model
        embedding_dim = int(self.cfg.get("mlm_embedding_dim", 64))
        encoder = DNAEncoder(
            vocab_size=len(vocab.itos),
            kmer=vocab.k,
            embedding_dim=embedding_dim,
            num_layers=2,  # Minimal for testing
            num_heads=4,
            ff_dim=embedding_dim * 2,
            dropout=0.0,  # No dropout for determinism
            pad_token_id=vocab.pad_id,
            drop_path_rate=0.0,
        )
        
        special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
        model = DNAMLM(
            encoder,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
            tie_weights=True,
            use_output_norm=True,
        )
        
        return model, train_loader, vocab
    
    def test_same_seed_same_batch_composition(self):
        """Test that same seed produces identical batch composition."""
        seed = 42
        
        # Get batch with first initialization
        model1, loader1, vocab1 = self._create_model_and_loader(seed)
        batch1 = next(iter(loader1))
        inputs1, labels1 = batch1[0], batch1[1]
        
        # Get batch with second initialization (same seed)
        model2, loader2, vocab2 = self._create_model_and_loader(seed)
        batch2 = next(iter(loader2))
        inputs2, labels2 = batch2[0], batch2[1]
        
        self.assertTrue(
            torch.equal(inputs1, inputs2),
            "Same seed should produce identical input tensors"
        )
        self.assertTrue(
            torch.equal(labels1, labels2),
            "Same seed should produce identical label tensors"
        )
    
    def test_same_seed_same_loss(self):
        """Test that same seed produces identical loss (within tolerance)."""
        seed = 42
        tolerance = 1e-5
        
        # First run
        model1, loader1, _ = self._create_model_and_loader(seed)
        model1.eval()  # Deterministic behavior
        batch1 = next(iter(loader1))
        inputs1, labels1 = batch1[0], batch1[1]
        
        with torch.no_grad():
            _, loss1 = model1(inputs1, labels1)
        
        # Second run (same seed, fresh model)
        model2, loader2, _ = self._create_model_and_loader(seed)
        model2.eval()
        batch2 = next(iter(loader2))
        inputs2, labels2 = batch2[0], batch2[1]
        
        with torch.no_grad():
            _, loss2 = model2(inputs2, labels2)
        
        diff = abs(loss1.item() - loss2.item())
        self.assertLessEqual(
            diff, tolerance,
            f"Same seed should produce identical loss (diff={diff:.2e} > tol={tolerance:.2e})"
        )
    
    def test_same_seed_same_backward_loss(self):
        """Test that same seed produces identical loss after forward/backward pass."""
        seed = 42
        tolerance = 1e-5
        
        # First run with backward pass
        model1, loader1, _ = self._create_model_and_loader(seed)
        model1.train()
        model1.zero_grad()
        batch1 = next(iter(loader1))
        inputs1, labels1 = batch1[0], batch1[1]
        
        _, loss1 = model1(inputs1, labels1)
        loss1.backward()
        loss1_val = loss1.item()
        
        # Second run with backward pass (same seed, fresh model)
        model2, loader2, _ = self._create_model_and_loader(seed)
        model2.train()
        model2.zero_grad()
        batch2 = next(iter(loader2))
        inputs2, labels2 = batch2[0], batch2[1]
        
        _, loss2 = model2(inputs2, labels2)
        loss2.backward()
        loss2_val = loss2.item()
        
        diff = abs(loss1_val - loss2_val)
        self.assertLessEqual(
            diff, tolerance,
            f"Same seed should produce identical loss after backward (diff={diff:.2e})"
        )
    
    def test_different_seeds_different_batches(self):
        """Test that different seeds produce different batch compositions."""
        seed1, seed2 = 42, 43
        
        # Get batch with seed1
        _, loader1, _ = self._create_model_and_loader(seed1)
        batch1 = next(iter(loader1))
        inputs1, labels1 = batch1[0], batch1[1]
        
        # Get batch with seed2
        _, loader2, _ = self._create_model_and_loader(seed2)
        batch2 = next(iter(loader2))
        inputs2, labels2 = batch2[0], batch2[1]
        
        # At least one should differ
        inputs_differ = not torch.equal(inputs1, inputs2)
        labels_differ = not torch.equal(labels1, labels2)
        
        self.assertTrue(
            inputs_differ or labels_differ,
            "Different seeds should produce different batch compositions"
        )


class TestDeterministicDataLoader(unittest.TestCase):
    """Tests for deterministic DataLoader behavior."""
    
    def test_deterministic_generator(self):
        """Test that deterministic generator produces reproducible shuffling."""
        from torch.utils.data import DataLoader, TensorDataset
        
        # Create simple dataset
        data = torch.arange(100).reshape(100, 1)
        dataset = TensorDataset(data)
        
        seed = 42
        
        # First loader
        kwargs1 = get_deterministic_dataloader_kwargs(seed, num_workers=0)
        loader1 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs1)
        order1 = torch.cat([b[0] for b in loader1]).flatten().tolist()
        
        # Second loader (same seed)
        kwargs2 = get_deterministic_dataloader_kwargs(seed, num_workers=0)
        loader2 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs2)
        order2 = torch.cat([b[0] for b in loader2]).flatten().tolist()
        
        self.assertEqual(order1, order2, "Same seed should produce same shuffle order")
    
    def test_different_seed_different_shuffle(self):
        """Test that different seeds produce different shuffling."""
        from torch.utils.data import DataLoader, TensorDataset
        
        data = torch.arange(100).reshape(100, 1)
        dataset = TensorDataset(data)
        
        # First loader with seed 42
        kwargs1 = get_deterministic_dataloader_kwargs(42, num_workers=0)
        loader1 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs1)
        order1 = torch.cat([b[0] for b in loader1]).flatten().tolist()
        
        # Second loader with seed 43
        kwargs2 = get_deterministic_dataloader_kwargs(43, num_workers=0)
        loader2 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs2)
        order2 = torch.cat([b[0] for b in loader2]).flatten().tolist()
        
        self.assertNotEqual(order1, order2, "Different seeds should produce different shuffle order")


def run_quick_reproducibility_test():
    """Quick standalone test that can be run without pytest."""
    print("=" * 70)
    print("QUICK REPRODUCIBILITY TEST")
    print("=" * 70)
    
    # Test 1: Seed utility works
    print("\n[Test 1] Seed utility reproducibility...")
    passed = verify_reproducibility(42)
    status = "✓ PASS" if passed else "✗ FAIL"
    print(f"  {status}")
    
    # Test 2: Same seed = same torch random
    print("\n[Test 2] Same seed produces identical torch values...")
    set_global_seed(42)
    vals1 = torch.rand(100)
    set_global_seed(42)
    vals2 = torch.rand(100)
    match = torch.equal(vals1, vals2)
    status = "✓ PASS" if match else "✗ FAIL"
    print(f"  {status}")
    
    # Test 3: Different seeds = different values
    print("\n[Test 3] Different seeds produce different values...")
    set_global_seed(42)
    vals1 = torch.rand(100)
    set_global_seed(43)
    vals2 = torch.rand(100)
    differ = not torch.equal(vals1, vals2)
    status = "✓ PASS" if differ else "✗ FAIL"
    print(f"  {status}")
    
    # Test 4: DataLoader determinism
    print("\n[Test 4] DataLoader shuffle is deterministic...")
    from torch.utils.data import DataLoader, TensorDataset
    
    data = torch.arange(50).reshape(50, 1)
    dataset = TensorDataset(data)
    
    kwargs1 = get_deterministic_dataloader_kwargs(42)
    loader1 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs1)
    order1 = torch.cat([b[0] for b in loader1]).flatten().tolist()
    
    kwargs2 = get_deterministic_dataloader_kwargs(42)
    loader2 = DataLoader(dataset, batch_size=10, shuffle=True, **kwargs2)
    order2 = torch.cat([b[0] for b in loader2]).flatten().tolist()
    
    match = order1 == order2
    status = "✓ PASS" if match else "✗ FAIL"
    print(f"  {status}")
    
    print("\n" + "=" * 70)
    all_passed = passed and match and differ
    if all_passed:
        print("ALL TESTS PASSED! ✓")
    else:
        print("SOME TESTS FAILED! ✗")
    print("=" * 70)
    
    return all_passed


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run reproducibility tests")
    parser.add_argument("--quick", action="store_true", help="Run quick standalone test")
    parser.add_argument("--full", action="store_true", help="Run full test suite with pytest")
    args = parser.parse_args()
    
    if args.quick or (not args.quick and not args.full):
        # Default to quick test
        success = run_quick_reproducibility_test()
        sys.exit(0 if success else 1)
    elif args.full:
        # Run full test suite
        unittest.main(argv=[''], exit=True, verbosity=2)
