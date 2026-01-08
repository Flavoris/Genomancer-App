from pathlib import Path
import sys

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from pre_training_checklist import DEFAULT_CONFIG_PATH, check_config


def test_pre_training_checklist_passes_default_config() -> None:
    assert check_config(DEFAULT_CONFIG_PATH) is True


def test_pre_training_checklist_fails_for_bad_config(tmp_path: Path) -> None:
    bad_config = {
        "training": {
            "stage1_lr": 0.00003,
            "grad_accum_steps": 2,
            "warmup_ratio": 0.05,
            "early_stopping_patience": 30,
            "layer_lr_decay": 1.0,
        },
        "data": {"batch_size": 8},
        "loss": {"stage1_use_focal_loss": False},
    }
    config_path = tmp_path / "config.yaml"
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(bad_config, handle, sort_keys=True)

    assert check_config(config_path) is False
