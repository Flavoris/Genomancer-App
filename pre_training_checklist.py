#!/usr/bin/env python3
"""
Pre-training checklist for gene_whisperer configuration values.
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


DEFAULT_CONFIG_PATH = (
    Path(__file__).resolve().parent / "gene_whisperer" / "configs" / "pretrain.yaml"
)
FLOAT_TOLERANCE = 1e-12


@dataclass(frozen=True)
class ConfigCheck:
    label: str
    keys: Sequence[str]
    expected: Any
    operator: str
    sections: Sequence[str] = ()


CHECKS: Sequence[ConfigCheck] = (
    ConfigCheck(
        label="mlm_lr",
        keys=("lr",),
        expected=0.0003,
        operator="==",
        sections=("training",),
    ),
    ConfigCheck(
        label="batch_size",
        keys=("batch_size",),
        expected=64,
        operator="==",
        sections=("training",),
    ),
    ConfigCheck(
        label="grad_accum_steps",
        keys=("grad_accum_steps",),
        expected=2,
        operator="==",
        sections=("training",),
    ),
    ConfigCheck(
        label="mask_prob",
        keys=("mask_prob",),
        expected=0.15,
        operator=">=",
        sections=("mlm",),
    ),
    ConfigCheck(
        label="vocab_size",
        keys=("vocab_size",),
        expected=4096,
        operator="==",
        sections=("mlm",),
    ),
    ConfigCheck(
        label="embedding_dim",
        keys=("embedding_dim",),
        expected=320,
        operator="==",
        sections=("model",),
    ),
    ConfigCheck(
        label="num_layers",
        keys=("num_layers",),
        expected=8,
        operator="==",
        sections=("model",),
    ),
)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load YAML config and return the parsed mapping."""
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


def _values_equal(actual: Any, expected: Any) -> bool:
    if isinstance(expected, bool):
        return isinstance(actual, bool) and actual is expected
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        actual_num = _normalize_number(actual)
        if actual_num is None:
            return False
        return abs(actual_num - float(expected)) <= FLOAT_TOLERANCE
    return actual == expected


def _compare_values(actual: Any, expected: Any, operator: str) -> bool:
    if operator == "==":
        return _values_equal(actual, expected)
    actual_num = _normalize_number(actual)
    expected_num = _normalize_number(expected)
    if actual_num is None or expected_num is None:
        return False
    if operator == ">=":
        return actual_num + FLOAT_TOLERANCE >= expected_num
    if operator == "<=":
        return actual_num <= expected_num + FLOAT_TOLERANCE
    if operator == "<":
        return actual_num < expected_num - FLOAT_TOLERANCE
    raise ValueError(f"Unsupported operator: {operator}")


def _format_scalar(value: Any) -> str:
    if value is None:
        return "missing"
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
    keys: Iterable[str],
    sections: Iterable[str] = (),
) -> Any:
    """Return the first matching key from root or fallback sections."""
    for key in keys:
        if key in config:
            return config[key]
    for section in sections:
        section_config = config.get(section)
        if isinstance(section_config, dict):
            for key in keys:
                if key in section_config:
                    return section_config[key]
    return None


def run_checks(
    config: dict[str, Any],
    checks: Sequence[ConfigCheck] = CHECKS,
) -> tuple[list[tuple[ConfigCheck, Any, bool]], bool]:
    """Evaluate checks and return results plus overall status."""
    results: list[tuple[ConfigCheck, Any, bool]] = []
    all_passed = True
    for check in checks:
        actual = get_config_value(config, check.keys, check.sections)
        passed = _compare_values(actual, check.expected, check.operator)
        results.append((check, actual, passed))
        all_passed = all_passed and passed
    return results, all_passed


def print_report(results: Sequence[tuple[ConfigCheck, Any, bool]], all_passed: bool) -> None:
    print("Pre-Training Configuration Checklist")
    print("=" * 50)
    for check, actual, passed in results:
        status = "OK" if passed else "FAIL"
        actual_display = _format_scalar(actual)
        expected_display = _format_scalar(check.expected)
        print(
            f"{status} {check.label}: {actual_display} "
            f"(expected {check.operator} {expected_display})"
        )
    print("=" * 50)
    if all_passed:
        print("ALL CHECKS PASSED - Ready for training!")
    else:
        print("SOME CHECKS FAILED - Fix before training!")


def check_config(config_path: Path | str = DEFAULT_CONFIG_PATH) -> bool:
    """Load config, run checks, and print a summary."""
    path = Path(config_path)
    config = load_config(path)
    results, all_passed = run_checks(config)
    print_report(results, all_passed)
    return all_passed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a pre-training checklist for config.yaml settings."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to config.yaml",
    )
    args = parser.parse_args()

    try:
        all_passed = check_config(args.config)
    except (FileNotFoundError, ValueError, yaml.YAMLError) as exc:
        print(f"SOME CHECKS FAILED - Fix before training!\n- {exc}")
        raise SystemExit(1) from exc

    raise SystemExit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
