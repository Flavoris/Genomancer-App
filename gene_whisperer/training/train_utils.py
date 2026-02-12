"""Training utilities shared across pretraining and fine-tuning."""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Iterable, List

import numpy as np
import torch
from torch import nn


@dataclass
class TrainingConfig:
    lr: float = 5e-4
    weight_decay: float = 0.01
    layer_decay: float = 0.9
    use_llrd: bool = False


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def _collect_params(module: nn.Module) -> List[nn.Parameter]:
    return [param for param in module.parameters() if param.requires_grad]


def _build_weight_decay_groups(model: nn.Module, weight_decay: float) -> List[dict]:
    decay_params: List[nn.Parameter] = []
    no_decay_params: List[nn.Parameter] = []

    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue

        is_bias = name.endswith(".bias")
        is_norm = "norm" in name.lower()
        is_vector = param.ndim == 1
        if is_bias or is_norm or is_vector:
            no_decay_params.append(param)
        else:
            decay_params.append(param)

    groups: List[dict] = []
    if decay_params:
        groups.append({"params": decay_params, "weight_decay": weight_decay})
    if no_decay_params:
        groups.append({"params": no_decay_params, "weight_decay": 0.0})
    return groups


def get_parameter_groups_with_llrd(
    model: nn.Module,
    base_lr: float,
    weight_decay: float,
    layer_decay: float,
) -> List[dict]:
    """Create parameter groups with layer-wise learning rate decay."""
    param_groups: List[dict] = []
    param_set = set()

    encoder = getattr(model, "encoder", None)
    if encoder is None:
        raise ValueError("Model must expose an encoder attribute for LLRD")

    embedding_modules = []
    for name in ("token_embed", "pos_embed"):
        module = getattr(encoder, name, None)
        if module is not None:
            embedding_modules.append(module)

    layer_modules = None
    if hasattr(encoder, "encoder") and hasattr(encoder.encoder, "layers"):
        layer_modules = list(encoder.encoder.layers)
    elif hasattr(encoder, "layers"):
        layer_modules = list(encoder.layers)
    else:
        layer_modules = []

    num_layers = len(layer_modules)

    if embedding_modules:
        embed_params: List[nn.Parameter] = []
        for module in embedding_modules:
            embed_params.extend(_collect_params(module))
        param_set.update(embed_params)
        param_groups.append(
            {
                "params": embed_params,
                "lr": base_lr * (layer_decay ** num_layers),
                "weight_decay": weight_decay,
                "name": "embedding",
            }
        )

    for idx, layer in enumerate(layer_modules):
        layer_params = _collect_params(layer)
        param_set.update(layer_params)
        lr = base_lr * (layer_decay ** (num_layers - idx - 1))
        param_groups.append(
            {
                "params": layer_params,
                "lr": lr,
                "weight_decay": weight_decay,
                "name": f"layer_{idx}",
            }
        )

    remaining_params = [
        param for param in model.parameters() if param.requires_grad and param not in param_set
    ]
    if remaining_params:
        param_groups.append(
            {
                "params": remaining_params,
                "lr": base_lr,
                "weight_decay": weight_decay,
                "name": "top",
            }
        )
    return param_groups


def build_optimizer(model: nn.Module, config: TrainingConfig) -> torch.optim.Optimizer:
    if config.use_llrd:
        groups = get_parameter_groups_with_llrd(
            model, base_lr=config.lr, weight_decay=config.weight_decay, layer_decay=config.layer_decay
        )
        return torch.optim.AdamW(groups)
    param_groups = _build_weight_decay_groups(model, config.weight_decay)
    return torch.optim.AdamW(param_groups, lr=config.lr)


def cycle_batches(loader: Iterable):
    while True:
        for batch in loader:
            yield batch
