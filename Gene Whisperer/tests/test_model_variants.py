"""Verification tests for Stage 1 model variants.

Run this file directly to execute all checks and print a comparison summary.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import torch
import torch.nn as nn

TRAINING_DIR = Path(__file__).resolve().parent.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

from model import (  # noqa: E402
    GeneWhispererCombined,
    GeneWhispererFeaturesOnly,
    GeneWhispererTransformerOnly,
    create_model_variant,
)

VARIANT_NAMES = ("transformer_only", "features_only", "combined")


@dataclass(frozen=True)
class VariantSummary:
    name: str
    parameter_count: int
    inputs_required: str
    output_shape: tuple[int, int]


def _small_config() -> dict:
    return {
        "vocab_size": 16,
        "kmer": 3,
        "embedding_dim": 16,
        "transformer_layers": 1,
        "transformer_heads": 2,
        "transformer_ff_dim": 32,
        "transformer_dropout": 0.0,
        "pad_token_id": 0,
        "engineered_dim": 6,
        "engineered_mlp_hidden": 8,
        "engineered_mlp_output": 6,
        "max_seq_len": 8,
        "drop_path_rate": 0.0,
        "use_relative_position_bias": False,
        "use_glu_ffn": False,
    }


def _set_seed(seed: int) -> None:
    torch.manual_seed(seed)


def _count_parameters(model: nn.Module) -> int:
    return sum(param.numel() for param in model.parameters())


def _assert_tensor_finite(tensor: torch.Tensor, label: str) -> None:
    if not torch.isfinite(tensor).all():
        raise AssertionError(f"{label} contains NaN or Inf values.")


def _assert_sequential_layout(
    module: nn.Sequential, expected_types: Sequence[type], label: str
) -> None:
    layers = list(module)
    assert len(layers) == len(expected_types), (
        f"{label} expected {len(expected_types)} layers, got {len(layers)}."
    )
    for idx, (layer, expected) in enumerate(zip(layers, expected_types)):
        assert isinstance(layer, expected), (
            f"{label} layer {idx} expected {expected.__name__}, got {type(layer).__name__}."
        )


def _assert_gradients_finite(model: nn.Module, label: str) -> None:
    saw_grad = False
    for name, param in model.named_parameters():
        if param.grad is None:
            continue
        saw_grad = True
        if not torch.isfinite(param.grad).all():
            raise AssertionError(f"{label} gradient has NaN/Inf values in {name}.")
    if not saw_grad:
        raise AssertionError(f"{label} produced no gradients.")


def _build_forward_inputs(config: dict, batch_size: int) -> tuple[torch.Tensor, torch.Tensor]:
    seq_len = config["max_seq_len"]
    tokens = torch.randint(1, config["vocab_size"], (batch_size, seq_len))
    engineered = torch.randn(batch_size, config["engineered_dim"])
    return tokens, engineered


def _build_training_batch(
    config: dict, batch_size: int
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    seq_len = config["max_seq_len"]
    split = max(2, config["vocab_size"] // 2)
    if split >= config["vocab_size"]:
        raise ValueError("vocab_size must be at least 3 for training batch generation.")

    half = batch_size // 2
    # Shift token ranges and feature means to keep the signal easy to learn.
    high_tokens = torch.randint(split, config["vocab_size"], (half, seq_len))
    low_tokens = torch.randint(1, split, (batch_size - half, seq_len))
    tokens = torch.cat([high_tokens, low_tokens], dim=0)

    engineered_pos = torch.randn(half, config["engineered_dim"]) + 1.0
    engineered_neg = torch.randn(batch_size - half, config["engineered_dim"]) - 1.0
    engineered = torch.cat([engineered_pos, engineered_neg], dim=0)

    labels = torch.zeros(batch_size, 1)
    labels[:half] = 1.0
    return tokens, engineered, labels

def _verify_transformer_architecture(model: GeneWhispererTransformerOnly, config: dict) -> None:
    assert isinstance(model, GeneWhispererTransformerOnly)
    assert model.encoder is not None
    assert model.embedding is not None
    assert model.classifier[0].in_features == config["embedding_dim"]
    _assert_sequential_layout(
        model.classifier,
        (nn.Linear, nn.GELU, nn.Dropout, nn.Linear),
        "transformer_only classifier",
    )


def _verify_features_architecture(model: GeneWhispererFeaturesOnly, config: dict) -> None:
    assert isinstance(model, GeneWhispererFeaturesOnly)
    assert model.encoder is None
    assert model.engineered_dim == config["engineered_dim"]
    assert model.classifier[0].in_features == config["engineered_dim"]
    _assert_sequential_layout(
        model.classifier,
        (nn.Linear, nn.GELU, nn.Dropout, nn.Linear, nn.GELU, nn.Dropout, nn.Linear),
        "features_only classifier",
    )


def _verify_combined_architecture(model: GeneWhispererCombined, config: dict) -> None:
    assert isinstance(model, GeneWhispererCombined)
    assert model.encoder is not None
    assert model.embedding is not None
    assert model.engineered_dim == config["engineered_dim"]
    assert model.engineered_mlp[0].in_features == config["engineered_dim"]
    assert model.classifier[0].in_features == (
        config["embedding_dim"] + config["engineered_mlp_output"]
    )
    _assert_sequential_layout(
        model.engineered_mlp,
        (nn.Linear, nn.GELU, nn.Dropout, nn.Linear, nn.GELU),
        "combined engineered_mlp",
    )
    _assert_sequential_layout(
        model.classifier,
        (nn.Linear, nn.GELU, nn.Dropout, nn.Linear),
        "combined classifier",
    )


def _assert_raises(expected: type[Exception] | tuple[type[Exception], ...], fn, label: str) -> None:
    try:
        fn()
    except expected:
        return
    except Exception as exc:
        raise AssertionError(f"{label} raised unexpected {type(exc).__name__}.") from exc
    raise AssertionError(f"{label} did not raise.")


def _verify_forward_pass(
    variant_name: str, model: nn.Module, config: dict
) -> tuple[int, int]:
    batch_size = 4
    tokens, engineered = _build_forward_inputs(config, batch_size)

    if variant_name == "transformer_only":
        logits = model(tokens, None)
        _assert_tensor_finite(logits, "transformer_only logits")
        assert logits.shape == (batch_size, 1)
        logits_with_features = model(tokens, engineered)
        _assert_tensor_finite(logits_with_features, "transformer_only logits with features")
    elif variant_name == "features_only":
        logits = model(None, engineered)
        _assert_tensor_finite(logits, "features_only logits")
        assert logits.shape == (batch_size, 1)
        _assert_raises(ValueError, lambda: model(tokens, None), "features_only missing features")
    elif variant_name == "combined":
        logits = model(tokens, engineered)
        _assert_tensor_finite(logits, "combined logits")
        assert logits.shape == (batch_size, 1)
        _assert_raises(ValueError, lambda: model(tokens, None), "combined missing features")
        _assert_raises(Exception, lambda: model(None, engineered), "combined missing tokens")
    else:
        raise ValueError(f"Unknown variant: {variant_name}")

    return tuple(logits.shape)


def _run_training_simulation(variant_name: str, model: nn.Module, config: dict) -> None:
    model.train()
    batch_size = 16
    tokens, engineered, labels = _build_training_batch(config, batch_size)
    loss_fn = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    losses: list[float] = []
    for _ in range(10):
        optimizer.zero_grad()
        if variant_name == "transformer_only":
            logits = model(tokens, None)
        elif variant_name == "features_only":
            logits = model(None, engineered)
        elif variant_name == "combined":
            logits = model(tokens, engineered)
        else:
            raise ValueError(f"Unknown variant: {variant_name}")
        loss = loss_fn(logits, labels)
        _assert_tensor_finite(loss, f"{variant_name} loss")
        loss.backward()
        _assert_gradients_finite(model, variant_name)
        optimizer.step()
        losses.append(loss.item())

    if losses[-1] >= losses[0]:
        raise AssertionError(
            f"{variant_name} loss did not decrease: start={losses[0]:.4f}, end={losses[-1]:.4f}"
        )


def _verify_variant(
    variant_name: str, config: dict, seed: int, verbose: bool
) -> VariantSummary:
    _set_seed(seed)
    model = create_model_variant(variant_name, config)
    param_count = _count_parameters(model)

    if variant_name == "transformer_only":
        _verify_transformer_architecture(model, config)
        inputs_required = "tokens"
    elif variant_name == "features_only":
        _verify_features_architecture(model, config)
        inputs_required = "engineered"
    elif variant_name == "combined":
        _verify_combined_architecture(model, config)
        inputs_required = "tokens + engineered"
    else:
        raise ValueError(f"Unknown variant: {variant_name}")

    output_shape = _verify_forward_pass(variant_name, model, config)
    _run_training_simulation(variant_name, model, config)

    if verbose:
        param_millions = param_count / 1e6
        print(f"{variant_name} parameters: {param_millions:.2f}M ({param_count})")
    return VariantSummary(
        name=variant_name,
        parameter_count=param_count,
        inputs_required=inputs_required,
        output_shape=output_shape,
    )


def _run_all_variants(verbose: bool) -> list[VariantSummary]:
    config = _small_config()
    seeds = {"transformer_only": 11, "features_only": 17, "combined": 23}
    summaries = []
    for variant in VARIANT_NAMES:
        summaries.append(_verify_variant(variant, config, seeds[variant], verbose))
    return summaries


def _format_summary_table(summaries: Iterable[VariantSummary]) -> str:
    header = "| Variant | Parameters | Inputs Required | Output Shape |"
    divider = "|---------|------------|-----------------|--------------|"
    rows = []
    for summary in summaries:
        param_millions = summary.parameter_count / 1e6
        shape_label = f"(B, {summary.output_shape[1]})"
        rows.append(
            f"| {summary.name} | {param_millions:.2f}M | {summary.inputs_required} | "
            f"{shape_label} |"
        )
    return "\n".join([header, divider, *rows])

def test_model_variants_end_to_end() -> None:
    _run_all_variants(verbose=False)


if __name__ == "__main__":
    summary = _run_all_variants(verbose=True)
    print()
    print(_format_summary_table(summary))
