"""Tests for Stage 1 cosine warmup scheduler."""
from __future__ import annotations

from pathlib import Path
import sys

import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).parent.parent))

from train_stage1 import get_cosine_schedule_with_warmup


def test_cosine_schedule_warmup_then_decay() -> None:
    model = nn.Linear(2, 2)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
    warmup_steps = 4
    total_steps = 12

    scheduler = get_cosine_schedule_with_warmup(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=0.01,
    )

    lrs: list[float] = []
    for _ in range(total_steps):
        optimizer.step()
        scheduler.step()
        lrs.append(optimizer.param_groups[0]["lr"])

    for step in range(warmup_steps):
        expected = 0.1 * ((step + 1) / warmup_steps)
        assert abs(lrs[step] - expected) < 1e-6

    assert abs(lrs[warmup_steps - 1] - 0.1) < 1e-6
    assert lrs[warmup_steps] < lrs[warmup_steps - 1]


def test_cosine_schedule_respects_min_lr_ratio_after_warmup() -> None:
    model = nn.Linear(2, 2)
    base_lr = 0.2
    optimizer = torch.optim.SGD(model.parameters(), lr=base_lr)
    warmup_steps = 3
    total_steps = 30
    min_lr_ratio = 0.2

    scheduler = get_cosine_schedule_with_warmup(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=min_lr_ratio,
    )

    lrs: list[float] = []
    for _ in range(total_steps + 1):
        optimizer.step()
        scheduler.step()
        lrs.append(optimizer.param_groups[0]["lr"])

    min_after_warmup = min(lrs[warmup_steps:])
    assert min_after_warmup >= base_lr * min_lr_ratio - 1e-6
