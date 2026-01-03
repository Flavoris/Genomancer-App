"""
Sanity checks for multi-scale loss composition.
"""
import sys

import torch

sys.path.insert(0, ".")

from train_stage1 import combine_multiscale_losses


def test_combine_multiscale_losses_basic():
    ensemble_loss = torch.tensor(1.5)
    individual_losses = [torch.tensor(2.0), torch.tensor(3.0)]
    total = combine_multiscale_losses(ensemble_loss, individual_losses, 0.2)
    expected = 1.5 + 0.2 * (2.0 + 3.0)
    assert torch.isclose(total, torch.tensor(expected))


def test_combine_multiscale_losses_empty():
    ensemble_loss = torch.tensor(0.7)
    total = combine_multiscale_losses(ensemble_loss, [], 0.2)
    assert torch.isclose(total, ensemble_loss)
