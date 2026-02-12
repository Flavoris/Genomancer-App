"""Optimization utilities for MLM pretraining."""
from __future__ import annotations

import math

import torch


def build_warmup_cosine_scheduler(
    optimizer: torch.optim.Optimizer,
    total_steps: int,
    warmup_ratio: float,
    min_lr_ratio: float,
) -> torch.optim.lr_scheduler.LambdaLR:
    warmup_steps = max(1, int(total_steps * warmup_ratio))

    def lr_lambda(current_step: int) -> float:
        if current_step < warmup_steps:
            return float(current_step + 1) / float(warmup_steps)

        progress_numerator = current_step - warmup_steps
        progress_denominator = max(1, total_steps - warmup_steps)
        progress = min(1.0, progress_numerator / progress_denominator)

        cosine = 0.5 * (1.0 + math.cos(math.pi * progress))
        return min_lr_ratio + (1.0 - min_lr_ratio) * cosine

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lr_lambda)


def build_grad_scaler(use_amp: bool, device: torch.device) -> torch.amp.GradScaler:
    return torch.amp.GradScaler(enabled=use_amp and device.type == "cuda")


def is_amp_enabled(use_amp: bool, device: torch.device) -> bool:
    return use_amp and device.type == "cuda"
