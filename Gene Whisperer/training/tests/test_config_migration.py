"""Tests for simplified config defaults and legacy migration warnings."""

from pathlib import Path
import sys
import warnings

import pytest
import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from train_stage1 import load_config


LEGACY_WARNING = "legacy CNN/TCN options"


def _write_yaml(path: Path, data: dict) -> None:
    path.write_text(yaml.safe_dump(data, sort_keys=True), encoding="utf-8")


def test_load_config_applies_simplified_defaults(tmp_path: Path) -> None:
    config_path = tmp_path / "simplified.yaml"
    _write_yaml(config_path, {"simplified_model": {"pooling_type": "attention"}})

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        loaded = load_config(config_path)
    assert captured == []
    simplified = loaded["simplified_model"]
    assert simplified["pooling_type"] == "attention"
    assert simplified["classifier_hidden"] == 256
    assert simplified["classifier_dropout"] == 0.15
    assert simplified["fusion_method"] == "concat"
    assert loaded["use_simplified_architecture"] is True


def test_load_config_warns_on_legacy_tcn_options(tmp_path: Path) -> None:
    config_path = tmp_path / "legacy.yaml"
    _write_yaml(config_path, {"use_tcn": True, "tcn_hidden": 384})

    with pytest.warns(UserWarning, match=LEGACY_WARNING):
        loaded = load_config(config_path)

    simplified = loaded["simplified_model"]
    assert simplified["pooling_type"] == "attention"
    assert simplified["classifier_hidden"] == 256
