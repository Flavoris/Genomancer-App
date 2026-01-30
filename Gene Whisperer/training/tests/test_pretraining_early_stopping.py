"""Tests for pretraining early stopping convergence checks."""
from __future__ import annotations

import torch

from early_stopping import PretrainingEarlyStopping


def _run_sequence(
    stopper: PretrainingEarlyStopping,
    losses: list[float],
    accuracies: list[float],
) -> bool:
    model = torch.nn.Linear(1, 1)
    should_stop = False
    for epoch, (loss, acc) in enumerate(zip(losses, accuracies), start=1):
        should_stop = stopper(
            val_loss=loss,
            model=model,
            epoch=epoch,
            checkpoint_dir=None,
            val_accuracy=acc,
        )
    return should_stop


def test_accuracy_stagnation_ignored_when_loss_improves() -> None:
    stopper = PretrainingEarlyStopping(
        patience=50,
        min_delta=0.001,
        min_epochs=1,
        accuracy_stagnation_window=5,
        accuracy_stagnation_threshold=0.01,
        rate_of_change_window=10,
        rate_of_change_min_improvement=0.0001,
        use_kmer_defaults=False,
    )

    losses = [7.0, 6.9, 6.8, 6.7, 6.6, 6.5]
    accuracies = [0.10, 0.101, 0.101, 0.101, 0.101, 0.101]

    should_stop = _run_sequence(stopper, losses, accuracies)
    assert not should_stop, "Accuracy stagnation should not stop when loss still improves"


def test_accuracy_stagnation_triggers_when_loss_plateaus() -> None:
    stopper = PretrainingEarlyStopping(
        patience=50,
        min_delta=0.001,
        min_epochs=1,
        accuracy_stagnation_window=5,
        accuracy_stagnation_threshold=0.01,
        rate_of_change_window=10,
        rate_of_change_min_improvement=0.0001,
        use_kmer_defaults=False,
    )

    losses = [7.0, 6.999, 6.998, 6.997, 6.996, 6.995]
    accuracies = [0.10, 0.101, 0.101, 0.101, 0.101, 0.101]

    should_stop = _run_sequence(stopper, losses, accuracies)
    assert should_stop, "Accuracy stagnation should stop when loss plateaus"
    assert stopper.stop_reason is not None
    assert "accuracy_stagnation" in stopper.stop_reason
