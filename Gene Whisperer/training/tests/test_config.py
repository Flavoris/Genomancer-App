"""
Test that config.yaml has the optimized hyperparameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

CONFIG_PATH = Path(__file__).resolve().parents[1] / "config.yaml"
FLOAT_TOLERANCE = 1e-4


@dataclass(frozen=True)
class HyperparameterCheck:
    key: str
    expected: float | int
    min_val: float | int
    max_val: float | int
    label: str


def load_config(config_path: Path) -> dict[str, Any]:
    """Load config.yaml and return the parsed mapping."""
    assert config_path.exists(), f"config.yaml not found at {config_path}"
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    assert isinstance(cfg, dict), "config.yaml should parse to a mapping"
    return cfg


def is_value_optimal(actual: Any, expected: float | int) -> bool:
    """Check closeness to expected value with float tolerance."""
    if isinstance(expected, float):
        if not isinstance(actual, (int, float)):
            return False
        return abs(float(actual) - expected) <= FLOAT_TOLERANCE
    return actual == expected


def is_value_in_range(actual: Any, min_val: float | int, max_val: float | int) -> bool:
    """Validate numeric ranges safely for YAML-loaded values."""
    if not isinstance(actual, (int, float)):
        return False
    return min_val <= float(actual) <= max_val


def test_config_hyperparameters() -> None:
    """Verify critical hyperparameters are updated."""
    cfg = load_config(CONFIG_PATH)

    checks = [
        HyperparameterCheck("max_bp_len", 81, 64, 128, "Stage max_bp_len"),
        HyperparameterCheck("lr", 0.00002, 0.00001, 0.00005, "Learning rate"),
        HyperparameterCheck("stage1_lr", 0.00002, 0.00001, 0.00005, "Stage 1 learning rate"),
        HyperparameterCheck("stage2_lr", 0.00001, 0.000005, 0.00002, "Stage 2 learning rate"),
        HyperparameterCheck("label_smoothing", 0.0, 0.0, 0.01, "Label smoothing"),
        HyperparameterCheck("mixup_alpha", 0.1, 0.05, 0.2, "Mixup alpha"),
        HyperparameterCheck("warmup_ratio", 0.10, 0.05, 0.15, "Warmup ratio"),
        HyperparameterCheck("weight_decay", 0.01, 0.005, 0.02, "Weight decay"),
        HyperparameterCheck("grad_accum_steps", 4, 2, 8, "Gradient accumulation"),
        HyperparameterCheck("epochs", 150, 100, 200, "Epochs"),
        HyperparameterCheck("early_stopping_patience", 15, 10, 25, "Early stopping patience"),
        HyperparameterCheck("engineered_dim", 288, 128, 320, "Engineered feature dimension"),
        HyperparameterCheck("mlm_window_size", 81, 64, 128, "MLM window size"),
        HyperparameterCheck(
            "post_cnn_transformer_layers",
            3,
            3,
            3,
            "Post-CNN transformer layers",
        ),
    ]

    failures: list[str] = []

    print("Hyperparameter Validation:")
    print("-" * 60)

    for check in checks:
        actual = cfg.get(check.key)

        if actual is None:
            print(f"FAIL {check.label}: NOT FOUND (expected {check.expected})")
            failures.append(f"{check.key}: missing")
            continue

        in_range = is_value_in_range(actual, check.min_val, check.max_val)
        is_optimal = is_value_optimal(actual, check.expected)

        status = "OK" if in_range else "FAIL"
        optimal_note = " (OPTIMAL)" if is_optimal else ""

        print(f"{status} {check.label}: {actual}{optimal_note}")
        print(f"   Expected: {check.expected}, Range: [{check.min_val}, {check.max_val}]")

        if not in_range:
            failures.append(
                f"{check.key}: {actual} not in range [{check.min_val}, {check.max_val}]"
            )

    print("-" * 60)

    all_passed = not failures
    if all_passed:
        print("\nOK ALL CONFIG CHECKS PASSED")
    else:
        print("\nFAIL SOME CONFIG VALUES OUT OF RANGE")

    assert all_passed, "Hyperparameter checks failed: " + ", ".join(failures)


def test_new_config_sections() -> None:
    """Verify new config sections are present."""
    cfg = load_config(CONFIG_PATH)

    new_sections = [
        ("ensemble_method", "Ensemble method configuration"),
        ("features", "Feature flags"),
        ("mlm_data", "MLM data configuration"),
    ]

    print("\nNew Config Sections:")
    print("-" * 60)

    for key, desc in new_sections:
        if key in cfg:
            print(f"OK {desc}: PRESENT")
            if isinstance(cfg[key], dict):
                for setting_key, setting_value in cfg[key].items():
                    print(f"   - {setting_key}: {setting_value}")
        else:
            print(f"WARN {desc}: NOT FOUND (optional)")



def test_separate_training_phase_settings() -> None:
    """Ensure stage/MLM settings are separated with backward-compatible root keys."""
    cfg = load_config(CONFIG_PATH)

    assert cfg.get("max_bp_len") == 81, f"max_bp_len should be 81, got {cfg.get('max_bp_len')}"
    assert cfg.get("mlm_window_size") == 81, (
        f"mlm_window_size should be 81, got {cfg.get('mlm_window_size')}"
    )
    assert cfg.get("lr") == 0.00002, f"lr should be 0.00002, got {cfg.get('lr')}"
    assert cfg.get("mlm_lr") == 0.0002, f"mlm_lr should be 0.0002, got {cfg.get('mlm_lr')}"
    assert cfg.get("stage1_lr") == 0.00002, "stage1_lr should be 0.00002"
    assert cfg.get("stage2_lr") == 0.00001, "stage2_lr should be 0.00001"


def test_mps_optimized_settings() -> None:
    """Verify MPS-optimized settings for M4 MacBook."""
    cfg = load_config(CONFIG_PATH)

    print("\nMPS Optimization Settings:")
    print("-" * 60)

    checks = [
        ("num_workers", 0, "Should be 0 for MPS"),
        ("amp_enabled", True, "Mixed precision should be enabled"),
        ("amp_dtype", "bfloat16", "bfloat16 recommended for MPS"),
    ]

    failures: list[str] = []

    for key, expected, reason in checks:
        actual = cfg.get(key)
        match = actual == expected
        status = "OK" if match else "FAIL"
        print(f"{status} {key}: {actual} (expected: {expected})")
        print(f"   {reason}")
        if not match:
            failures.append(f"{key}: {actual} != {expected}")

    assert not failures, "MPS config checks failed: " + ", ".join(failures)


def _run_check(check_fn) -> bool:
    """Run a check and capture failures for CLI usage."""
    try:
        check_fn()
        return True
    except AssertionError as exc:
        print(f"\nFAIL {check_fn.__name__}: {exc}")
        return False


if __name__ == "__main__":
    print("=" * 60)
    print("PROMPT 4 VALIDATION: Hyperparameter Configuration")
    print("=" * 60)

    all_ok = True
    all_ok &= _run_check(test_config_hyperparameters)
    print()
    all_ok &= _run_check(test_new_config_sections)
    print()
    all_ok &= _run_check(test_separate_training_phase_settings)
    print()
    all_ok &= _run_check(test_mps_optimized_settings)

    print("\n" + "=" * 60)
    print("CONFIG VALIDATION COMPLETE")
    print("=" * 60)

    raise SystemExit(0 if all_ok else 1)
