"""Unit tests for progressive length curriculum.

Acceptance tests:
1. Length schedule matches config
2. No shape mismatches during ramps
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from typing import List

import torch

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent))

from length_curriculum import (
    LengthCurriculumConfig,
    LengthCurriculumScheduler,
    CurriculumMLMCollator,
    create_curriculum_collator,
)


class MockVocab:
    """Mock vocabulary for testing."""
    
    def __init__(self, k: int = 3, vocab_size: int = 67):
        self.k = k
        self.itos = [f"tok_{i}" for i in range(vocab_size)]
        self.stoi = {t: i for i, t in enumerate(self.itos)}
        self.mask_id = vocab_size - 3
        self.unk_id = vocab_size - 2
        self.pad_id = vocab_size - 1
        self._base_token_ids = list(range(vocab_size - 3))


class TestLengthScheduleMatchesConfig(unittest.TestCase):
    """Test: Length schedule matches config values."""
    
    def test_start_length_matches_config(self):
        """At step 0, length should equal start_tokens."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.0,
            end_pct=0.3,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        self.assertEqual(scheduler.get_current_length(0), 32)
    
    def test_end_length_matches_config(self):
        """After ramp, length should equal end_tokens."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.0,
            end_pct=0.3,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        # At 30% of training (end of ramp)
        self.assertEqual(scheduler.get_current_length(300), 79)
        # After ramp
        self.assertEqual(scheduler.get_current_length(500), 79)
        self.assertEqual(scheduler.get_current_length(1000), 79)
    
    def test_ramp_boundaries_match_config_percentages(self):
        """Ramp should start and end at configured percentages."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.1,  # 10%
            end_pct=0.4,    # 40%
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        # Before 10%: should be at start length
        self.assertEqual(scheduler.get_current_length(0), 32)
        self.assertEqual(scheduler.get_current_length(99), 32)
        
        # At 10%: ramp starts
        length_at_start = scheduler.get_current_length(100)
        self.assertEqual(length_at_start, 32)
        
        # At 40%: ramp ends, should be at end length
        self.assertEqual(scheduler.get_current_length(400), 79)
    
    def test_disabled_curriculum_returns_end_length(self):
        """When disabled, should always return end_tokens."""
        config = LengthCurriculumConfig(
            enabled=False,
            start_tokens=32,
            end_tokens=79,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        for step in [0, 100, 500, 1000]:
            self.assertEqual(
                scheduler.get_current_length(step), 79,
                f"Disabled curriculum should return 79 at step {step}"
            )
    
    def test_monotonic_increase_during_ramp(self):
        """Length should never decrease during ramp."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.0,
            end_pct=0.5,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        prev_length = 0
        for step in range(0, 1001, 10):
            current_length = scheduler.get_current_length(step)
            self.assertGreaterEqual(
                current_length, prev_length,
                f"Length decreased from {prev_length} to {current_length} at step {step}"
            )
            prev_length = current_length
    
    def test_config_validation(self):
        """Config should reject invalid values."""
        # start_tokens too small
        with self.assertRaises(ValueError):
            LengthCurriculumConfig(start_tokens=4, end_tokens=79)
        
        # end_tokens <= start_tokens
        with self.assertRaises(ValueError):
            LengthCurriculumConfig(start_tokens=80, end_tokens=79)
        
        # Invalid percentages
        with self.assertRaises(ValueError):
            LengthCurriculumConfig(start_pct=-0.1, end_pct=0.5)
        
        with self.assertRaises(ValueError):
            LengthCurriculumConfig(start_pct=0.5, end_pct=0.3)


class TestNoShapeMismatchDuringRamps(unittest.TestCase):
    """Test: No shape mismatches during length transitions."""
    
    def setUp(self):
        """Set up mock vocabulary and test data."""
        self.vocab = MockVocab(k=3, vocab_size=67)
        self.max_tokens = 79
        self.batch_size = 8
    
    def _create_test_batch(self, seq_len: int) -> List[torch.LongTensor]:
        """Create a batch of random token sequences."""
        return [
            torch.randint(0, 64, (seq_len,), dtype=torch.long)
            for _ in range(self.batch_size)
        ]
    
    def test_collator_output_shape_matches_curriculum_length(self):
        """Output shape should match current curriculum length."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.0,
            end_pct=0.5,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        collator = CurriculumMLMCollator(
            vocab=self.vocab,
            scheduler=scheduler,
            mask_prob=0.15,
        )
        
        # Test at various steps during ramp
        test_steps = [0, 100, 250, 500, 750, 1000]
        
        for step in test_steps:
            scheduler.step(step)
            expected_length = scheduler.current_length
            
            # Create batch with full-length sequences
            batch = self._create_test_batch(self.max_tokens)
            
            # Collate
            masked_inputs, labels = collator(batch)
            
            # Verify shapes
            self.assertEqual(
                masked_inputs.shape, (self.batch_size, expected_length),
                f"Step {step}: masked_inputs shape mismatch. "
                f"Expected (8, {expected_length}), got {masked_inputs.shape}"
            )
            self.assertEqual(
                labels.shape, (self.batch_size, expected_length),
                f"Step {step}: labels shape mismatch. "
                f"Expected (8, {expected_length}), got {labels.shape}"
            )
    
    def test_shape_consistency_through_entire_ramp(self):
        """Shapes should be consistent through every step of the ramp."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=64,
            start_pct=0.0,
            end_pct=1.0,
        )
        scheduler = LengthCurriculumScheduler(total_steps=100, config=config)
        collator = CurriculumMLMCollator(
            vocab=self.vocab,
            scheduler=scheduler,
            mask_prob=0.15,
        )
        
        # Test every step
        for step in range(101):
            scheduler.step(step)
            expected_length = scheduler.current_length
            
            batch = self._create_test_batch(self.max_tokens)
            masked_inputs, labels = collator(batch)
            
            self.assertEqual(masked_inputs.shape[1], expected_length)
            self.assertEqual(labels.shape[1], expected_length)
    
    def test_batch_dimension_preserved(self):
        """Batch dimension should be preserved regardless of length."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=16,
            end_tokens=64,
        )
        scheduler = LengthCurriculumScheduler(total_steps=100, config=config)
        collator = CurriculumMLMCollator(
            vocab=self.vocab,
            scheduler=scheduler,
            mask_prob=0.15,
        )
        
        for batch_size in [1, 4, 8, 16, 32]:
            batch = [
                torch.randint(0, 64, (self.max_tokens,), dtype=torch.long)
                for _ in range(batch_size)
            ]
            
            masked_inputs, labels = collator(batch)
            
            self.assertEqual(
                masked_inputs.shape[0], batch_size,
                f"Batch size mismatch: expected {batch_size}, got {masked_inputs.shape[0]}"
            )
    
    def test_truncation_when_curriculum_shorter_than_input(self):
        """When curriculum length < input length, should truncate."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=64,
        )
        scheduler = LengthCurriculumScheduler(total_steps=100, config=config)
        collator = CurriculumMLMCollator(
            vocab=self.vocab,
            scheduler=scheduler,
            mask_prob=0.15,
        )
        
        # At step 0, length should be 32
        scheduler.step(0)
        self.assertEqual(scheduler.current_length, 32)
        
        # Create batch with 79 tokens
        batch = self._create_test_batch(79)
        masked_inputs, labels = collator(batch)
        
        # Output should be truncated to 32
        self.assertEqual(masked_inputs.shape[1], 32)
        self.assertEqual(labels.shape[1], 32)
    
    def test_padding_when_curriculum_longer_than_input(self):
        """When curriculum length > input length, should pad."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=40,
            end_tokens=64,
        )
        scheduler = LengthCurriculumScheduler(total_steps=100, config=config)
        scheduler.step(100)  # At end, length = 64
        
        collator = CurriculumMLMCollator(
            vocab=self.vocab,
            scheduler=scheduler,
            mask_prob=0.15,
        )
        
        # Create batch with only 32 tokens
        batch = self._create_test_batch(32)
        masked_inputs, labels = collator(batch)
        
        # Output should be padded to 64
        self.assertEqual(masked_inputs.shape[1], 64)
        
        # Last 32 positions should be pad tokens
        pad_positions = masked_inputs[:, 32:]
        pad_count = (pad_positions == self.vocab.pad_id).sum().item()
        expected_pads = self.batch_size * 32
        self.assertEqual(pad_count, expected_pads)


class TestSchedulerStepMethod(unittest.TestCase):
    """Test the step() method of the scheduler."""
    
    def test_step_updates_internal_state(self):
        """step() should update current_step and current_length."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        self.assertEqual(scheduler.current_step, 0)
        self.assertEqual(scheduler.current_length, 32)
        
        scheduler.step(500)
        
        self.assertEqual(scheduler.current_step, 500)
        self.assertGreater(scheduler.current_length, 32)
    
    def test_step_without_argument_increments(self):
        """step() without argument should increment by 1."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        for expected_step in range(1, 11):
            scheduler.step()
            self.assertEqual(scheduler.current_step, expected_step)
    
    def test_get_progress_returns_correct_phase(self):
        """get_progress() should return correct phase information."""
        config = LengthCurriculumConfig(
            enabled=True,
            start_tokens=32,
            end_tokens=79,
            start_pct=0.1,
            end_pct=0.5,
        )
        scheduler = LengthCurriculumScheduler(total_steps=1000, config=config)
        
        # Before ramp starts
        scheduler.step(50)
        progress = scheduler.get_progress()
        self.assertEqual(progress["phase"], "warmup")
        
        # During ramp
        scheduler.step(300)
        progress = scheduler.get_progress()
        self.assertEqual(progress["phase"], "ramping")
        self.assertTrue(0 < progress["ramp_progress"] < 1)
        
        # After ramp
        scheduler.step(600)
        progress = scheduler.get_progress()
        self.assertEqual(progress["phase"], "full")
        self.assertEqual(progress["ramp_progress"], 1.0)


class TestFromConfigFactory(unittest.TestCase):
    """Test the from_config factory method."""
    
    def test_creates_config_from_dict(self):
        """from_config should create config from dictionary."""
        cfg = {
            "mlm_length_curriculum_enabled": True,
            "mlm_length_curriculum_start_tokens": 48,
            "mlm_length_curriculum_end_tokens": 96,
            "mlm_length_curriculum_start_pct": 0.05,
            "mlm_length_curriculum_end_pct": 0.4,
        }
        
        config = LengthCurriculumConfig.from_config(cfg)
        
        self.assertTrue(config.enabled)
        self.assertEqual(config.start_tokens, 48)
        self.assertEqual(config.end_tokens, 96)
        self.assertEqual(config.start_pct, 0.05)
        self.assertEqual(config.end_pct, 0.4)
    
    def test_uses_defaults_for_missing_keys(self):
        """from_config should use defaults for missing keys."""
        cfg = {}
        
        config = LengthCurriculumConfig.from_config(cfg)
        
        self.assertFalse(config.enabled)  # Default disabled
        self.assertEqual(config.start_tokens, 32)
        self.assertEqual(config.end_tokens, 79)


class TestCreateCurriculumCollator(unittest.TestCase):
    """Test the create_curriculum_collator factory function."""
    
    def test_returns_collator_and_scheduler(self):
        """Factory should return both collator and scheduler."""
        vocab = MockVocab()
        cfg = {
            "mlm_length_curriculum_enabled": True,
            "mlm_length_curriculum_start_tokens": 32,
            "mlm_length_curriculum_end_tokens": 64,
        }
        
        collator, scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=1000,
            cfg=cfg,
        )
        
        self.assertIsInstance(collator, CurriculumMLMCollator)
        self.assertIsInstance(scheduler, LengthCurriculumScheduler)
    
    def test_scheduler_is_shared_with_collator(self):
        """Collator and scheduler should share the same scheduler instance."""
        vocab = MockVocab()
        cfg = {
            "mlm_length_curriculum_enabled": True,
            "mlm_length_curriculum_start_tokens": 32,
            "mlm_length_curriculum_end_tokens": 64,
        }
        
        collator, scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=1000,
            cfg=cfg,
        )
        
        # When we update the scheduler, collator should see the change
        scheduler.step(500)
        self.assertEqual(collator.scheduler.current_step, 500)


def run_tests():
    """Run all tests and print summary."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestLengthScheduleMatchesConfig))
    suite.addTests(loader.loadTestsFromTestCase(TestNoShapeMismatchDuringRamps))
    suite.addTests(loader.loadTestsFromTestCase(TestSchedulerStepMethod))
    suite.addTests(loader.loadTestsFromTestCase(TestFromConfigFactory))
    suite.addTests(loader.loadTestsFromTestCase(TestCreateCurriculumCollator))
    
    # Run with verbosity
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 70)
    print("LENGTH CURRICULUM TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    print()
    
    if result.wasSuccessful():
        print("✓ ALL ACCEPTANCE TESTS PASSED")
        print("  - Length schedule matches config ✓")
        print("  - No shape mismatches during ramps ✓")
    else:
        print("✗ SOME TESTS FAILED")
        if result.failures:
            print("\nFailures:")
            for test, traceback in result.failures:
                print(f"  - {test}: {traceback.split(chr(10))[0]}")
        if result.errors:
            print("\nErrors:")
            for test, traceback in result.errors:
                print(f"  - {test}: {traceback.split(chr(10))[0]}")
    
    print("=" * 70)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)

