from pathlib import Path
import sys

import torch

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.training.pretrain_optim import (
    build_grad_scaler,
    build_warmup_cosine_scheduler,
    is_amp_enabled,
)


def test_warmup_cosine_scheduler_decreases_lr() -> None:
    model = torch.nn.Linear(4, 4)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    scheduler = build_warmup_cosine_scheduler(
        optimizer=optimizer,
        total_steps=10,
        warmup_ratio=0.2,
        min_lr_ratio=0.1,
    )

    lrs = []
    for _ in range(10):
        optimizer.step()
        scheduler.step()
        lrs.append(optimizer.param_groups[0]["lr"])

    assert lrs[0] > 0
    assert min(lrs) >= 1e-4
    assert lrs[-1] <= lrs[2]


def test_amp_helpers_cpu_disabled() -> None:
    device = torch.device("cpu")
    assert not is_amp_enabled(use_amp=True, device=device)
    scaler = build_grad_scaler(use_amp=True, device=device)
    assert scaler.is_enabled() is False
