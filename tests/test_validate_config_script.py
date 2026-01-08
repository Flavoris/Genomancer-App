from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from validate_config import DEFAULT_CONFIG_PATH, _format_scalar, load_config, validate_config_values


def test_validate_config_values_pass() -> None:
    config = load_config(DEFAULT_CONFIG_PATH)
    errors, values = validate_config_values(config)
    assert errors == []
    assert values["effective_batch_size"] == 64


def test_format_scalar_handles_small_floats() -> None:
    assert _format_scalar(0.00002) == "0.00002"
    assert _format_scalar(0.1) == "0.1"
