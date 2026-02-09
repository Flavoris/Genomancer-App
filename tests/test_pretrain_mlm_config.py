from pathlib import Path
import sys

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.training.pretrain_mlm import _load_config


def test_load_config_defaults_and_relative_paths(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    config_dir.mkdir(parents=True)
    config_path = config_dir / "pretrain.yaml"

    config_data = {
        "mlm": {
            "fasta_paths": ["../data/test.fna"],
            "tokenizer_path": "../artifacts/tokenizer.json",
            "vocab_size": 128,
            "max_length": 64,
            "window_size": 32,
            "mask_prob": 0.2,
        },
        "model": {
            "embedding_dim": 32,
            "num_layers": 2,
            "num_heads": 4,
            "ff_dim": 64,
            "dropout": 0.1,
        },
        "training": {
            "batch_size": 8,
            "epochs": 2,
            "lr": 1e-4,
            "weight_decay": 0.01,
            "grad_accum_steps": 1,
            "seed": 7,
            "output_dir": "../artifacts/mlm",
        },
    }
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config_data, handle, sort_keys=False)

    config = _load_config(config_path)

    assert config.fasta_paths[0] == (tmp_path / "data" / "test.fna").resolve()
    assert config.tokenizer_path == (tmp_path / "artifacts" / "tokenizer.json").resolve()
    assert config.output_dir == (tmp_path / "artifacts" / "mlm").resolve()
    assert config.samples_per_epoch == 4096
    assert config.log_interval == 50
    assert config.num_workers == 0
    assert config.max_bases_per_file is None


def test_load_config_with_optional_fields(tmp_path: Path) -> None:
    config_path = tmp_path / "pretrain.yaml"
    config_data = {
        "mlm": {
            "fasta_paths": ["/tmp/a.fna"],
            "tokenizer_path": "/tmp/tokenizer.json",
            "vocab_size": 256,
            "max_length": 128,
            "window_size": 64,
            "mask_prob": 0.15,
            "max_bases_per_file": 1000,
        },
        "model": {
            "embedding_dim": 64,
            "num_layers": 4,
            "num_heads": 4,
            "ff_dim": 128,
            "dropout": 0.1,
        },
        "training": {
            "batch_size": 16,
            "epochs": 3,
            "samples_per_epoch": 500,
            "log_interval": 10,
            "num_workers": 2,
            "lr": 2e-4,
            "weight_decay": 0.01,
            "grad_accum_steps": 2,
            "seed": 9,
            "output_dir": "/tmp/out",
        },
    }
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config_data, handle, sort_keys=False)

    config = _load_config(config_path)

    assert config.max_bases_per_file == 1000
    assert config.samples_per_epoch == 500
    assert config.log_interval == 10
    assert config.num_workers == 2
