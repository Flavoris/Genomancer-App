from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from validate_config import DEFAULT_CONFIG_PATH, load_config


REMOVED_MODEL_KEYS = (
    "use_tcn",
    "tcn_hidden",
    "tcn_levels",
    "tcn_kernel",
    "multiscale_channels",
    "multiscale_kernels",
    "post_cnn_transformer_layers",
)


def test_simplified_model_block_present() -> None:
    config = load_config(DEFAULT_CONFIG_PATH)
    simplified = config.get("simplified_model")
    assert isinstance(simplified, dict)
    assert simplified["pooling_type"] == "attention"
    assert simplified["classifier_hidden"] == 256
    assert simplified["classifier_dropout"] == 0.15
    assert simplified["fusion_method"] == "concat"


def test_removed_model_options_absent() -> None:
    config = load_config(DEFAULT_CONFIG_PATH)
    model = config.get("model")
    for key in REMOVED_MODEL_KEYS:
        assert key not in config
        if isinstance(model, dict):
            assert key not in model
