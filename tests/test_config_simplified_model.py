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


def test_model_block_present() -> None:
    config = load_config(DEFAULT_CONFIG_PATH)
    model = config.get("model")
    assert isinstance(model, dict)
    assert model["embedding_dim"] == 256
    assert model["num_layers"] == 6
    assert model["max_length"] == 128


def test_removed_model_options_absent() -> None:
    config = load_config(DEFAULT_CONFIG_PATH)
    model = config.get("model")
    for key in REMOVED_MODEL_KEYS:
        assert key not in config
        if isinstance(model, dict):
            assert key not in model
