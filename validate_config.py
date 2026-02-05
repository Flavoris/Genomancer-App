#!/usr/bin/env python3
"""
Validate gene_whisperer fine-tuning config defaults.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    import yaml
except ImportError as exc:  # pragma: no cover - exercised in CLI usage
    raise SystemExit("PyYAML is required (pip install pyyaml).") from exc


CHECK_MARK = "\u2713"
DEFAULT_CONFIG_PATH = (
    Path(__file__).resolve().parent / "gene_whisperer" / "configs" / "finetune.yaml"
)
FLOAT_TOLERANCE = 1e-12


@dataclass(frozen=True)
class ConfigCheck:
    key: str
    expected: Any
    sections: Sequence[str] = ()


CONFIG_CHECKS: Sequence[ConfigCheck] = (
    ConfigCheck("lr", 0.00002, ("training",)),
    ConfigCheck("batch_size", 64, ("training",)),
    ConfigCheck("grad_accum_steps", 1, ("training",)),
    ConfigCheck("use_llrd", True, ("training",)),
    ConfigCheck("layer_decay", 0.9, ("training",)),
    ConfigCheck("embedding_dim", 256, ("model",)),
    ConfigCheck("num_layers", 6, ("model",)),
    ConfigCheck("max_length", 128, ("model",)),
    ConfigCheck("use_pstnp", True, ("stage1",)),
)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load YAML config and return a mapping."""
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError("config.yaml must parse to a mapping")
    return config


def _normalize_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _values_match(actual: Any, expected: Any) -> bool:
    if isinstance(expected, bool):
        return isinstance(actual, bool) and actual is expected
    if isinstance(expected, int) and not isinstance(expected, bool):
        actual_num = _normalize_number(actual)
        return actual_num is not None and abs(actual_num - expected) <= FLOAT_TOLERANCE
    if isinstance(expected, float):
        actual_num = _normalize_number(actual)
        return actual_num is not None and abs(actual_num - expected) <= FLOAT_TOLERANCE
    return actual == expected


def _format_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        formatted = f"{value:.5f}".rstrip("0").rstrip(".")
        return formatted or "0"
    return str(value)


def get_config_value(
    config: dict[str, Any],
    key: str,
    sections: Iterable[str] = (),
) -> Any:
    """Return config value from root or fallback sections."""
    if key in config:
        return config[key]
    for section in sections:
        section_config = config.get(section)
        if isinstance(section_config, dict) and key in section_config:
            return section_config[key]
    return None


def validate_config_values(config: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    """Validate required config values and return errors plus extracted values."""
    errors: list[str] = []
    extracted: dict[str, Any] = {}

    for check in CONFIG_CHECKS:
        actual = get_config_value(config, check.key, check.sections)
        extracted[check.key] = actual
        if actual is None:
            errors.append(f"{check.key}: missing (expected {check.expected})")
            continue
        if not _values_match(actual, check.expected):
            errors.append(f"{check.key}: expected {check.expected}, got {actual}")

    batch_size = _normalize_number(extracted.get("batch_size"))
    grad_accum_steps = _normalize_number(extracted.get("grad_accum_steps"))
    if batch_size is None or grad_accum_steps is None:
        errors.append("effective_batch_size: batch_size/grad_accum_steps invalid")
        effective_batch_size = None
    else:
        effective_batch_size = int(batch_size * grad_accum_steps)
        if effective_batch_size != 64:
            errors.append(
                f"effective_batch_size: expected 64, got {effective_batch_size}"
            )

    extracted["effective_batch_size"] = effective_batch_size
    return errors, extracted


def print_success(values: dict[str, Any]) -> None:
    print("CONFIG VALIDATION PASSED")
    print()
    print(f"lr: {_format_scalar(values.get('lr'))} {CHECK_MARK}")
    print(f"batch_size: {_format_scalar(values.get('batch_size'))} {CHECK_MARK}")
    print(f"grad_accum_steps: {_format_scalar(values.get('grad_accum_steps'))} {CHECK_MARK}")
    print(
        "effective_batch_size: "
        f"{_format_scalar(values.get('effective_batch_size'))} {CHECK_MARK}"
    )
    print(f"use_llrd: {_format_scalar(values.get('use_llrd'))} {CHECK_MARK}")
    print(f"layer_decay: {_format_scalar(values.get('layer_decay'))} {CHECK_MARK}")


def print_errors(errors: Sequence[str]) -> None:
    print("CONFIG VALIDATION FAILED")
    for error in errors:
        print(f"- {error}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate finetune.yaml defaults for gene_whisperer"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to config.yaml",
    )
    args = parser.parse_args()

    try:
        config = load_config(args.config)
    except (FileNotFoundError, ValueError, yaml.YAMLError) as exc:
        print(f"CONFIG VALIDATION FAILED\n- {exc}")
        raise SystemExit(1) from exc

    errors, values = validate_config_values(config)
    if errors:
        print_errors(errors)
        raise SystemExit(1)
    print_success(values)


if __name__ == "__main__":
    main()
