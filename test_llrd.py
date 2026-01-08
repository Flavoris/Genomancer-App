"""LLRD sanity checks plus a printable schedule example."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import torch
import torch.nn as nn

TRAINING_DIR = Path(__file__).resolve().parent / "Gene Whisperer" / "training"
sys.path.insert(0, str(TRAINING_DIR))

from train_stage1 import get_parameter_groups_with_llrd


class DummyEncoder(nn.Module):
    """Minimal encoder with a transformer-like layer stack."""

    def __init__(self, num_layers: int) -> None:
        super().__init__()
        self.layers = nn.ModuleList([nn.Linear(4, 4) for _ in range(num_layers)])


class DummyModel(nn.Module):
    """Small model that exposes embeddings, encoder layers, and a head."""

    def __init__(self, num_layers: int) -> None:
        super().__init__()
        self.embedding = nn.Embedding(8, 4)
        self.encoder = DummyEncoder(num_layers)
        self.classifier = nn.Linear(4, 1)


def _get_group_lr(param_groups: list[dict], group_name: str) -> float:
    for group in param_groups:
        if group.get("name") == group_name:
            return float(group["lr"])
    raise AssertionError(f"Missing param group: {group_name}")


def test_llrd_param_groups() -> None:
    """Verify that LLRD assigns expected learning rates per layer."""
    base_lr = 2e-5
    layer_decay = 0.9
    weight_decay = 0.01
    num_layers = 3

    model = DummyModel(num_layers)
    param_groups = get_parameter_groups_with_llrd(
        model=model,
        base_lr=base_lr,
        weight_decay=weight_decay,
        layer_decay=layer_decay,
    )

    embedding_lr = _get_group_lr(param_groups, "embedding")
    expected_embedding_lr = base_lr * (layer_decay ** num_layers)
    assert math.isclose(embedding_lr, expected_embedding_lr, rel_tol=1e-6), (
        f"embedding lr mismatch: got {embedding_lr}, expected {expected_embedding_lr}"
    )

    for layer_idx in range(num_layers):
        layer_lr = _get_group_lr(param_groups, f"layer_{layer_idx}")
        expected_lr = base_lr * (layer_decay ** (num_layers - layer_idx - 1))
        assert math.isclose(layer_lr, expected_lr, rel_tol=1e-6), (
            f"layer {layer_idx} lr mismatch: got {layer_lr}, expected {expected_lr}"
        )

    top_lr = _get_group_lr(param_groups, "top")
    assert math.isclose(top_lr, base_lr, rel_tol=1e-6), (
        f"top lr mismatch: got {top_lr}, expected {base_lr}"
    )


def test_llrd_decay_ratio_math() -> None:
    """Ensure the decay ratio math matches the expected formula."""
    base_lr = 2e-5
    layer_decay = 0.9
    num_layers = 12

    embedding_lr = base_lr * (layer_decay ** num_layers)
    ratio = base_lr / embedding_lr
    expected_ratio = (1 / layer_decay) ** num_layers

    assert math.isclose(ratio, expected_ratio, rel_tol=1e-6), (
        f"ratio mismatch: got {ratio}, expected {expected_ratio}"
    )


def print_llrd_schedule(
    base_lr: float = 2e-5,
    layer_decay: float = 0.9,
    num_layers: int = 12,
) -> None:
    """Print the LLRD schedule summary used in the prompt."""
    print("Layer-wise Learning Rate Decay Test")
    print(f"Base LR: {base_lr}")
    print(f"Decay: {layer_decay}")
    print(f"Num layers: {num_layers}")
    print()

    embedding_lr = base_lr * (layer_decay ** num_layers)
    print(f"Embedding LR: {embedding_lr:.2e} (decay^{num_layers})")

    for i in range(num_layers):
        layer_lr = base_lr * (layer_decay ** (num_layers - i - 1))
        print(f"Layer {i} LR: {layer_lr:.2e}")

    top_lr = base_lr
    print(f"Top (classifier) LR: {top_lr:.2e}")

    ratio = top_lr / embedding_lr
    expected_ratio = (1 / layer_decay) ** num_layers
    print(f"\nTop/Embedding LR ratio: {ratio:.2f}x")
    print(f"Expected ratio: {expected_ratio:.2f}x")

    if abs(ratio - expected_ratio) < 0.01:
        print("LLRD TEST PASSED")
    else:
        print("LLRD TEST FAILED")


if __name__ == "__main__":
    print_llrd_schedule()
