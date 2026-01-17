"""
Tests for quick_variant_sanity_check.py

These tests verify the sanity check utilities work correctly without
running full training loops.
"""

import sys

import torch
import torch.nn as nn

sys.path.insert(0, ".")

from quick_variant_sanity_check import (
    check_for_nan_inf,
    check_gradients,
    create_sanity_dataloader,
)


def test_check_for_nan_inf_clean_tensor():
    """Test NaN/Inf check on clean tensor."""
    tensor = torch.randn(10, 10)
    issues = check_for_nan_inf(tensor, "test_tensor")
    assert len(issues) == 0, f"Clean tensor should have no issues: {issues}"


def test_check_for_nan_inf_with_nan():
    """Test NaN/Inf check detects NaN values."""
    tensor = torch.randn(10, 10)
    tensor[0, 0] = float("nan")
    issues = check_for_nan_inf(tensor, "test_tensor")
    assert len(issues) == 1, "Should detect NaN"
    assert "NaN" in issues[0], "Should mention NaN"


def test_check_for_nan_inf_with_inf():
    """Test NaN/Inf check detects Inf values."""
    tensor = torch.randn(10, 10)
    tensor[0, 0] = float("inf")
    issues = check_for_nan_inf(tensor, "test_tensor")
    assert len(issues) == 1, "Should detect Inf"
    assert "Inf" in issues[0], "Should mention Inf"


def test_check_for_nan_inf_with_both():
    """Test NaN/Inf check detects both NaN and Inf."""
    tensor = torch.randn(10, 10)
    tensor[0, 0] = float("nan")
    tensor[1, 1] = float("inf")
    issues = check_for_nan_inf(tensor, "test_tensor")
    assert len(issues) == 2, "Should detect both NaN and Inf"


def test_check_gradients_with_gradients():
    """Test gradient check on model with gradients."""
    model = nn.Linear(10, 2)
    x = torch.randn(4, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()

    ok, issues = check_gradients(model)
    assert ok, f"Gradient check should pass: {issues}"
    assert len(issues) == 0


def test_check_gradients_without_gradients():
    """Test gradient check on model without backward pass."""
    model = nn.Linear(10, 2)
    # No backward call - no gradients computed

    ok, issues = check_gradients(model)
    assert not ok, "Should fail without gradients"
    assert any("No gradients" in issue for issue in issues)


def test_check_gradients_with_nan_gradient():
    """Test gradient check detects NaN in gradients."""
    model = nn.Linear(10, 2)
    x = torch.randn(4, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()

    # Manually inject NaN into gradient
    model.weight.grad[0, 0] = float("nan")

    ok, issues = check_gradients(model)
    assert not ok, "Should fail with NaN gradient"
    assert any("NaN" in issue for issue in issues)


class MockDataset(torch.utils.data.Dataset):
    """Mock dataset for testing dataloader creation."""

    def __init__(self, size: int = 1000):
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        return torch.zeros(10), torch.zeros(5), torch.tensor(0)


def test_create_sanity_dataloader_limits_samples():
    """Test that sanity dataloader respects max_samples."""
    full_dataset = MockDataset(size=1000)
    full_loader = torch.utils.data.DataLoader(full_dataset, batch_size=64)

    sanity_loader = create_sanity_dataloader(full_loader, max_samples=100, batch_size=32)

    assert len(sanity_loader.dataset) == 100, "Should limit to 100 samples"


def test_create_sanity_dataloader_small_dataset():
    """Test dataloader with dataset smaller than max_samples."""
    full_dataset = MockDataset(size=50)
    full_loader = torch.utils.data.DataLoader(full_dataset, batch_size=64)

    sanity_loader = create_sanity_dataloader(full_loader, max_samples=100, batch_size=32)

    assert len(sanity_loader.dataset) == 50, "Should use all 50 samples"


def test_create_sanity_dataloader_batch_size():
    """Test that sanity dataloader uses specified batch size."""
    full_dataset = MockDataset(size=1000)
    full_loader = torch.utils.data.DataLoader(full_dataset, batch_size=64)

    sanity_loader = create_sanity_dataloader(full_loader, max_samples=100, batch_size=16)

    assert sanity_loader.batch_size == 16, "Should use batch_size=16"


if __name__ == "__main__":
    # Run tests manually
    test_check_for_nan_inf_clean_tensor()
    test_check_for_nan_inf_with_nan()
    test_check_for_nan_inf_with_inf()
    test_check_for_nan_inf_with_both()
    test_check_gradients_with_gradients()
    test_check_gradients_without_gradients()
    test_check_gradients_with_nan_gradient()
    test_create_sanity_dataloader_limits_samples()
    test_create_sanity_dataloader_small_dataset()
    test_create_sanity_dataloader_batch_size()
    print("All tests passed!")
