"""Test early stopping implementation with intelligent convergence detection."""

import torch
import torch.nn as nn
import tempfile
from pathlib import Path
from early_stopping import EarlyStopping, PretrainingEarlyStopping


def test_basic_early_stopping():
    """Test basic early stopping behavior."""
    print("Test 1: Basic Early Stopping")
    print("-" * 40)

    # Create simple model
    model = nn.Linear(10, 2)

    # Create early stopper with patience=3
    stopper = EarlyStopping(patience=3, min_delta=0.01, min_epochs=2, verbose=False)

    # Simulate training: loss decreases then plateaus
    losses = [1.0, 0.8, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5]

    stopped_at = None
    for epoch, loss in enumerate(losses):
        should_stop = stopper(loss, model, epoch)
        if should_stop:
            stopped_at = epoch
            break

    print(f"  Losses: {losses}")
    print(f"  Stopped at epoch: {stopped_at}")
    print(f"  Best score: {stopper.best_score}")
    print(f"  Best epoch: {stopper.best_epoch}")

    # Should stop at epoch 6 (patience=3 after epoch 3)
    assert stopped_at == 6, f"Expected stop at epoch 6, got {stopped_at}"
    assert stopper.best_score == 0.5, f"Expected best_score=0.5, got {stopper.best_score}"
    assert stopper.best_epoch == 3, f"Expected best_epoch=3, got {stopper.best_epoch}"

    print("  PASSED")
    print()


def test_min_delta():
    """Test that min_delta prevents stopping on tiny improvements."""
    print("Test 2: Min Delta")
    print("-" * 40)

    model = nn.Linear(10, 2)

    # min_delta=0.01 means improvements < 0.01 don't count
    stopper = EarlyStopping(patience=2, min_delta=0.01, min_epochs=0, verbose=False)

    # Tiny improvements that shouldn't count
    losses = [1.0, 0.999, 0.998, 0.997, 0.996]

    stopped_at = None
    for epoch, loss in enumerate(losses):
        should_stop = stopper(loss, model, epoch)
        if should_stop:
            stopped_at = epoch
            break

    print(f"  Losses: {losses}")
    print(f"  Min delta: 0.01")
    print(f"  Stopped at epoch: {stopped_at}")

    # Should stop at epoch 2 (0.001 improvement doesn't count)
    assert stopped_at == 2, f"Expected stop at epoch 2, got {stopped_at}"
    print("  PASSED")
    print()


def test_min_epochs():
    """Test that early stopping respects min_epochs."""
    print("Test 3: Min Epochs")
    print("-" * 40)

    model = nn.Linear(10, 2)

    # min_epochs=5 means can't stop before epoch 5
    stopper = EarlyStopping(patience=1, min_delta=0.001, min_epochs=5, verbose=False)

    # Constant loss from start
    losses = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    stopped_at = None
    for epoch, loss in enumerate(losses):
        should_stop = stopper(loss, model, epoch)
        if should_stop:
            stopped_at = epoch
            break

    print(f"  Losses: {losses}")
    print(f"  Min epochs: 5, Patience: 1")
    print(f"  Stopped at epoch: {stopped_at}")

    # Should stop at epoch 6 (min_epochs=5, then patience=1)
    assert stopped_at == 6, f"Expected stop at epoch 6, got {stopped_at}"
    print("  PASSED")
    print()


def test_checkpoint_saving():
    """Test that best checkpoint is saved and restored."""
    print("Test 4: Checkpoint Saving/Restoring")
    print("-" * 40)

    model = nn.Linear(10, 2)

    # Get initial weights
    initial_weight = model.weight.data.clone()

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = EarlyStopping(patience=2, min_delta=0.01, min_epochs=0, verbose=False)

        # Epoch 0: loss=1.0, save checkpoint
        stopper(1.0, model, 0, Path(tmpdir))

        # Modify model weights
        model.weight.data = torch.randn_like(model.weight.data)
        modified_weight = model.weight.data.clone()

        # Epoch 1: loss=0.5, better, save new checkpoint
        stopper(0.5, model, 1, Path(tmpdir))
        best_weight = model.weight.data.clone()

        # Epoch 2: loss=0.9, worse
        model.weight.data = torch.randn_like(model.weight.data)
        stopper(0.9, model, 2, Path(tmpdir))

        # Epoch 3: loss=0.9, worse again, should stop
        stopper(0.9, model, 3, Path(tmpdir))

        # Restore best weights
        stopper.restore_best_weights(model)

        restored_weight = model.weight.data

        print(f"  Best epoch: {stopper.best_epoch}")
        print(f"  Best score: {stopper.best_score}")

        # Restored weights should match epoch 1 weights
        assert torch.allclose(restored_weight, best_weight), "Restored weights don't match best"
        print("  PASSED")
    print()


def test_pretraining_early_stopping():
    """Test pretraining-specific early stopping."""
    print("Test 5: Pretraining Early Stopping (basic)")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=2,
            min_delta=0.01,
            min_epochs=0,
            kmer=6,
            verbose=False,
            use_kmer_defaults=False,  # Disable to use explicit values
        )

        # Simulate pretraining with loss and perplexity
        epochs_data = [
            (1.0, 10.0),  # (loss, perplexity)
            (0.5, 5.0),
            (0.4, 4.0),
            (0.4, 4.0),
            (0.4, 4.0),
        ]

        for epoch, (loss, ppl) in enumerate(epochs_data):
            should_stop = stopper(loss, model, epoch, Path(tmpdir), ppl)
            if should_stop:
                break

        summary = stopper.get_summary()
        print(f"  K-mer: {stopper.kmer}")
        print(f"  Best loss: {summary['best_score']:.4f}")
        print(f"  Best epoch: {summary['best_epoch']}")
        print(f"  Perplexity history: {stopper.perplexity_history}")

        # Check checkpoint file name includes k-mer
        assert stopper.best_checkpoint_path.name == f"mlm_encoder_k6_best.pt"
        print("  PASSED")
    print()


def test_kmer_defaults():
    """Test that k-mer specific defaults are applied correctly."""
    print("Test 6: K-mer Specific Defaults")
    print("-" * 40)

    # Test each k-mer size gets appropriate defaults
    expected_defaults = {
        3: {"min_epochs": 5, "patience": 8, "target_accuracy": 0.85},
        4: {"min_epochs": 8, "patience": 10, "target_accuracy": 0.82},
        5: {"min_epochs": 10, "patience": 12, "target_accuracy": 0.80},
        6: {"min_epochs": 12, "patience": 15, "target_accuracy": 0.78},
    }

    for kmer, expected in expected_defaults.items():
        stopper = PretrainingEarlyStopping(
            kmer=kmer,
            use_kmer_defaults=True,
            verbose=False,
        )

        print(f"  k={kmer}: min_epochs={stopper.min_epochs}, "
              f"patience={stopper.patience}, target_acc={stopper.target_accuracy}")

        assert stopper.min_epochs == expected["min_epochs"], \
            f"k={kmer}: expected min_epochs={expected['min_epochs']}, got {stopper.min_epochs}"
        assert stopper.patience == expected["patience"], \
            f"k={kmer}: expected patience={expected['patience']}, got {stopper.patience}"
        assert stopper.target_accuracy == expected["target_accuracy"], \
            f"k={kmer}: expected target_accuracy={expected['target_accuracy']}, got {stopper.target_accuracy}"

    print("  PASSED")
    print()


def test_target_accuracy_stop():
    """Test that training stops when target accuracy is reached."""
    print("Test 7: Target Accuracy Stop ('Good Enough')")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=100,  # High patience so we don't stop from that
            min_delta=0.001,
            min_epochs=2,  # Low min_epochs
            kmer=3,
            target_accuracy=0.80,  # Stop at 80% accuracy
            target_loss=0.1,  # Very low, won't trigger
            use_kmer_defaults=False,
            verbose=False,
        )

        # Simulate training where accuracy reaches target
        epochs_data = [
            (1.0, None, 0.50),  # (loss, perplexity, accuracy)
            (0.9, None, 0.60),
            (0.8, None, 0.70),
            (0.7, None, 0.75),
            (0.6, None, 0.82),  # Exceeds target - should stop here
            (0.5, None, 0.85),
        ]

        stopped_at = None
        for epoch, (loss, ppl, acc) in enumerate(epochs_data):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=ppl,
                val_accuracy=acc,
            )
            if should_stop:
                stopped_at = epoch
                break

        print(f"  Target accuracy: 80%")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        assert stopped_at == 4, f"Expected stop at epoch 4, got {stopped_at}"
        assert "target_accuracy_reached" in stopper.stop_reason
        print("  PASSED")
    print()


def test_target_loss_stop():
    """Test that training stops when target loss is reached."""
    print("Test 8: Target Loss Stop ('Good Enough')")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=100,
            min_delta=0.001,
            min_epochs=2,
            kmer=3,
            target_accuracy=0.99,  # Very high, won't trigger
            target_loss=0.4,  # Stop when loss drops below 0.4
            use_kmer_defaults=False,
            verbose=False,
        )

        # Simulate training where loss reaches target
        epochs_data = [
            (1.0, None, 0.50),
            (0.8, None, 0.55),
            (0.6, None, 0.60),
            (0.35, None, 0.65),  # Below target loss - should stop here
            (0.3, None, 0.70),
        ]

        stopped_at = None
        for epoch, (loss, ppl, acc) in enumerate(epochs_data):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=ppl,
                val_accuracy=acc,
            )
            if should_stop:
                stopped_at = epoch
                break

        print(f"  Target loss: 0.4")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        assert stopped_at == 3, f"Expected stop at epoch 3, got {stopped_at}"
        assert "target_loss_reached" in stopper.stop_reason
        print("  PASSED")
    print()


def test_accuracy_stagnation():
    """Test that training stops when accuracy stagnates."""
    print("Test 9: Accuracy Stagnation Detection")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=100,  # High patience
            min_delta=0.001,
            min_epochs=3,
            kmer=3,
            accuracy_stagnation_window=5,
            accuracy_stagnation_threshold=0.01,  # Need >1% improvement
            target_accuracy=0.99,  # Won't trigger
            target_loss=0.01,  # Won't trigger
            use_kmer_defaults=False,
            verbose=False,
        )

        # Simulate training where accuracy stagnates (very flat - less than 1% over 5 epochs)
        # Loss keeps improving but accuracy is essentially stuck
        epochs_data = [
            (1.00, None, 0.70),
            (0.95, None, 0.72),
            (0.90, None, 0.74),
            (0.85, None, 0.742),  # Stagnation starts (min_epochs passed)
            (0.80, None, 0.744),
            (0.75, None, 0.745),
            (0.70, None, 0.746),
            (0.65, None, 0.747),  # 5 epochs with <1% improvement - should stop
            (0.60, None, 0.75),
        ]

        stopped_at = None
        for epoch, (loss, ppl, acc) in enumerate(epochs_data):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=ppl,
                val_accuracy=acc,
            )
            if should_stop:
                stopped_at = epoch
                break

        if stopped_at is not None:
            print(f"  Accuracy over epochs: {[d[2] for d in epochs_data[:stopped_at+1]]}")
        else:
            print(f"  Accuracy over epochs: {[d[2] for d in epochs_data]}")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        assert stopped_at is not None, "Should have stopped due to accuracy stagnation"
        assert "accuracy_stagnation" in stopper.stop_reason
        print("  PASSED")
    print()


def test_diminishing_returns():
    """Test that training stops when improvement rate drops (diminishing returns)."""
    print("Test 10: Diminishing Returns Detection")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=100,  # High patience
            min_delta=0.0001,  # Small delta so small improvements count
            min_epochs=5,
            kmer=3,
            rate_of_change_window=10,
            rate_of_change_min_improvement=0.01,  # Need >0.01 improvement per epoch
            accuracy_stagnation_window=100,  # Disable accuracy check
            accuracy_stagnation_threshold=0.001,
            target_accuracy=0.99,  # Won't trigger
            target_loss=0.01,  # Won't trigger
            use_kmer_defaults=False,
            verbose=False,
        )

        # Simulate training where loss improvement rate drops significantly
        # First few epochs: big improvements, then tiny improvements
        losses = [
            1.0, 0.8, 0.6, 0.5, 0.45,  # Good improvement rate
            0.44, 0.435, 0.43, 0.428, 0.426,  # Slowing down
            0.425, 0.424, 0.4235, 0.423, 0.4228,  # Very slow - should trigger
            0.4227, 0.4226,
        ]

        stopped_at = None
        for epoch, loss in enumerate(losses):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=None,
                val_accuracy=0.5 + epoch * 0.01,  # Some accuracy improvement
            )
            if should_stop:
                stopped_at = epoch
                break

        print(f"  Losses: {losses[:stopped_at+1] if stopped_at is not None else losses}")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        assert stopped_at is not None, "Should have stopped due to diminishing returns"
        assert "diminishing_returns" in stopper.stop_reason
        print("  PASSED")
    print()


def test_realistic_k3_scenario():
    """Test a realistic k=3 training scenario that should stop early."""
    print("Test 11: Realistic k=3 Scenario (should stop quickly)")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Use k=3 defaults
        stopper = PretrainingEarlyStopping(
            kmer=3,
            use_kmer_defaults=True,
            verbose=False,
        )

        print(f"  K=3 settings: min_epochs={stopper.min_epochs}, "
              f"patience={stopper.patience}, target_acc={stopper.target_accuracy}")

        # Simulate realistic k=3 training: fast initial progress, then plateau
        # This mimics what the user experienced: 74% -> 76% stagnation
        epochs_data = [
            (2.5, None, 0.40),   # Epoch 0
            (1.8, None, 0.55),   # Epoch 1
            (1.2, None, 0.65),   # Epoch 2
            (0.9, None, 0.72),   # Epoch 3
            (0.7, None, 0.74),   # Epoch 4
            (0.65, None, 0.745), # Epoch 5 - min_epochs passed, stagnation starts
            (0.62, None, 0.75),  # Epoch 6
            (0.60, None, 0.752), # Epoch 7
            (0.58, None, 0.755), # Epoch 8
            (0.57, None, 0.758), # Epoch 9
            (0.56, None, 0.76),  # Epoch 10 - should stop by here due to stagnation
            (0.55, None, 0.762),
            (0.54, None, 0.764),
        ]

        stopped_at = None
        for epoch, (loss, ppl, acc) in enumerate(epochs_data):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=ppl,
                val_accuracy=acc,
            )
            if should_stop:
                stopped_at = epoch
                break

        print(f"  Final accuracy: {epochs_data[stopped_at][2] if stopped_at is not None else 'N/A'}")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        # Should stop reasonably early due to accuracy stagnation
        assert stopped_at is not None, "Should have stopped"
        assert stopped_at <= 12, f"Should stop within 12 epochs, stopped at {stopped_at}"
        print("  PASSED")
    print()


def test_no_false_positives_good_training():
    """Test that early stopping doesn't trigger during good training progress."""
    print("Test 12: No False Positives During Good Training")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            kmer=3,
            use_kmer_defaults=True,
            verbose=False,
        )

        # Simulate good training with steady improvement
        epochs_data = [
            (2.5, None, 0.40),
            (2.0, None, 0.50),
            (1.5, None, 0.60),
            (1.0, None, 0.70),
            (0.7, None, 0.78),
            (0.5, None, 0.83),
            (0.4, None, 0.86),  # Should stop here - reached target accuracy 85%
        ]

        stopped_at = None
        for epoch, (loss, ppl, acc) in enumerate(epochs_data):
            should_stop = stopper(
                val_loss=loss,
                model=model,
                epoch=epoch,
                checkpoint_dir=Path(tmpdir),
                val_perplexity=ppl,
                val_accuracy=acc,
            )
            if should_stop:
                stopped_at = epoch
                break

        print(f"  Accuracies: {[d[2] for d in epochs_data]}")
        print(f"  Stopped at epoch: {stopped_at}")
        print(f"  Stop reason: {stopper.stop_reason}")

        # Should stop due to reaching a "good enough" threshold, not stagnation
        assert stopped_at is not None
        assert "target" in stopper.stop_reason, \
            f"Expected 'good enough' stop reason, got: {stopper.stop_reason}"
        print("  PASSED")
    print()


def run_all_tests():
    print("=" * 60)
    print("EARLY STOPPING TESTS (with intelligent convergence detection)")
    print("=" * 60)
    print()

    test_basic_early_stopping()
    test_min_delta()
    test_min_epochs()
    test_checkpoint_saving()
    test_pretraining_early_stopping()
    test_kmer_defaults()
    test_target_accuracy_stop()
    test_target_loss_stop()
    test_accuracy_stagnation()
    test_diminishing_returns()
    test_realistic_k3_scenario()
    test_no_false_positives_good_training()

    print("=" * 60)
    print("ALL EARLY STOPPING TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
