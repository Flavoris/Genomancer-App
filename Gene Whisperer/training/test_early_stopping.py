"""Test early stopping implementation."""

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

    print("  ✓ Basic early stopping test PASSED")
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
    print("  ✓ Min delta test PASSED")
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
    print("  ✓ Min epochs test PASSED")
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
        print("  ✓ Checkpoint saving/restoring test PASSED")
    print()


def test_pretraining_early_stopping():
    """Test pretraining-specific early stopping."""
    print("Test 5: Pretraining Early Stopping")
    print("-" * 40)

    model = nn.Linear(10, 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        stopper = PretrainingEarlyStopping(
            patience=2,
            min_delta=0.01,
            min_epochs=0,
            kmer=6,
            verbose=False,
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
        print("  ✓ Pretraining early stopping test PASSED")
    print()


def run_all_tests():
    print("=" * 50)
    print("EARLY STOPPING TESTS")
    print("=" * 50)
    print()

    test_basic_early_stopping()
    test_min_delta()
    test_min_epochs()
    test_checkpoint_saving()
    test_pretraining_early_stopping()

    print("=" * 50)
    print("ALL EARLY STOPPING TESTS PASSED ✓")
    print("=" * 50)


if __name__ == "__main__":
    run_all_tests()
