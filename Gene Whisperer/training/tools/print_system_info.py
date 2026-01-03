#!/usr/bin/env python3
"""Print system information for Gene Whisperer training reproducibility.

Prints:
- Python version
- PyTorch version
- Device selected (cuda/mps/cpu)
- Total model parameter count for Stage 1 model instance
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Add parent directory to path so we can import model and config
script_dir = Path(__file__).resolve().parent
training_dir = script_dir.parent
sys.path.insert(0, str(training_dir))

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import torch
import yaml

from model import GeneWhispererStage1


def get_device() -> torch.device:
    """Select the best available device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    else:
        return torch.device("cpu")


def load_config(config_path: Path) -> dict:
    """Load YAML configuration."""
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def build_stage1_model(cfg: dict) -> GeneWhispererStage1:
    """Build a Stage 1 model instance from config."""
    return GeneWhispererStage1(
        vocab_size=int(cfg.get("vocab_size", 67)),
        kmer=int(cfg.get("kmer", 3)),
        embedding_dim=int(cfg.get("embedding_dim", 256)),
        num_layers=int(cfg.get("transformer_layers", 6)),
        num_heads=int(cfg.get("transformer_heads", 8)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1024)),
        dropout=float(cfg.get("transformer_dropout", 0.15)),
        pad_token_id=int(cfg.get("pad_token_id", 66)),
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 256)),
        tcn_levels=int(cfg.get("tcn_levels", 4)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        fusion_hidden=int(cfg.get("fusion_hidden", 256)),
    )


def get_system_info(cfg: dict) -> dict:
    """Collect all system information."""
    device = get_device()
    model = build_stage1_model(cfg)
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    return {
        "python_version": sys.version,
        "torch_version": torch.__version__,
        "device": str(device),
        "cuda_available": torch.cuda.is_available(),
        "mps_available": hasattr(torch.backends, "mps") and torch.backends.mps.is_available(),
        "stage1_total_params": total_params,
        "stage1_trainable_params": trainable_params,
        "stage1_total_params_millions": round(total_params / 1e6, 2),
    }


def print_system_info(info: dict, output_format: str = "text") -> None:
    """Print system information in the specified format."""
    if output_format == "json":
        print(json.dumps(info, indent=2))
    else:
        print("=" * 60)
        print("GENE WHISPERER SYSTEM INFO")
        print("=" * 60)
        print(f"Python version:      {info['python_version'].split()[0]}")
        print(f"PyTorch version:     {info['torch_version']}")
        print(f"Device selected:     {info['device']}")
        print(f"  CUDA available:    {info['cuda_available']}")
        print(f"  MPS available:     {info['mps_available']}")
        print("-" * 60)
        print(f"Stage 1 Model Parameters:")
        print(f"  Total params:      {info['stage1_total_params']:,} ({info['stage1_total_params_millions']}M)")
        print(f"  Trainable params:  {info['stage1_trainable_params']:,}")
        print("=" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(description="Print Gene Whisperer system info")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output in JSON format",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (training_dir / config_path).resolve()

    if not config_path.exists():
        fallback = training_dir / "config.yaml"
        if fallback.exists():
            config_path = fallback
        else:
            print(f"Error: Config file not found: {config_path}", file=sys.stderr)
            sys.exit(1)

    cfg = load_config(config_path)
    info = get_system_info(cfg)
    print_system_info(info, output_format="json" if args.json else "text")


if __name__ == "__main__":
    main()
