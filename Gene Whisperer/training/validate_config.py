#!/usr/bin/env python3
"""
Validate config.yaml against a strict schema.
Fails hard on parsing errors or unknown/malformed keys.
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Type, Union

import yaml


# =============================================================================
# SCHEMA DEFINITION
# =============================================================================
# Each key maps to its expected type(s). Use tuple for multiple allowed types.
# None means the value can be null.

CONFIG_SCHEMA: Dict[str, Union[Type, tuple]] = {
    # Random seed
    "seed": int,
    
    # Data settings
    "max_bp_len": int,
    "batch_size": int,
    "delimiter": str,
    "train_val_split": float,
    "stage1_train": str,
    "stage1_val": str,
    "stage2_train": str,
    "stage2_val": str,
    "stage2_test": str,
    "stage2_strength_positive": int,
    "stage2_strength_negative": int,
    
    # Model architecture
    "embedding_dim": int,
    "transformer_layers": int,
    "transformer_heads": int,
    "transformer_ff_dim": int,
    "transformer_dropout": float,
    "use_alibi": bool,
    "use_positional_bias": bool,
    "use_attention_pool": bool,
    
    # K-mer tokenization
    "kmer": int,
    "vocab_size": int,
    "pad_token_id": int,
    "vocab_cache_dir": str,
    
    # TCN + Multi-scale CNN
    "use_tcn": bool,
    "tcn_hidden": int,
    "tcn_levels": int,
    "tcn_kernel": int,
    "multiscale_channels": int,
    "multiscale_kernels": list,
    "lstm_hidden": int,
    
    # Post-CNN transformer
    "post_cnn_transformer_layers": int,
    
    # Engineered features
    "stage1_use_engineered_features": bool,
    "stage1_engineered_dim": int,
    "engineered_mlp_hidden": int,
    "engineered_mlp_output": int,
    
    # Fusion
    "fusion_hidden": int,
    
    # Augmentation
    "stage1_reverse_complement_prob": float,
    
    # Training settings
    "lr": float,
    "weight_decay": float,
    "warmup_ratio": float,
    "epochs": int,
    "early_stopping_patience": int,
    "grad_accum_steps": int,
    "grad_clip_norm": float,
    "layer_lr_decay": float,
    "drop_path_rate": float,
    
    # Loss function
    "stage1_use_focal_loss": bool,
    "stage1_focal_alpha": float,
    "stage1_focal_gamma": float,
    "label_smoothing": float,
    
    # Mixup augmentation
    "use_mixup": bool,
    "mixup_alpha": float,
    
    # SWA
    "use_swa": bool,
    "swa_start_epoch": int,
    "swa_lr": float,
    "swa_anneal_epochs": int,
    
    # MLM pretraining
    "stage1_load_mlm_weights": bool,
    "stage2_load_mlm_weights": bool,
    "stage1_freeze_layers": int,
    "stage1_unfreeze_epoch": (int, type(None)),  # Can be null
    "mlm_encoder_checkpoint": str,
    "include_reverse_complements": bool,
    "mlm_fasta_path": str,
    "mlm_window_size": int,
    "mlm_stride": int,
    "mlm_batch_size": int,
    "mlm_epochs": int,
    "mlm_lr": float,
    "mlm_weight_decay": float,
    "mlm_kmer": int,
    "mlm_vocab_path": str,
    "mlm_encoder_path": str,
    "mlm_embedding_dim": int,
    "mlm_transformer_layers": int,
    "mlm_transformer_heads": int,
    "mlm_transformer_ff_dim": int,
    "mlm_transformer_dropout": float,
    "mlm_tie_weights": bool,
    "mlm_use_output_norm": bool,
    "mlm_use_layer_lr_decay": bool,
    "mlm_layer_lr_decay": float,
    "mlm_val_ratio": float,
    "mlm_grad_accum_steps": int,
    "mlm_warmup_ratio": float,
    "mlm_patience": int,
    
    # Checkpoints
    "stage1_checkpoint_name": (str, type(None)),  # Can be null
    
    # Visualization
    "stage1_embed_plot_epoch": (int, type(None)),  # Can be null
    "stage1_embed_plot_method": str,
    "stage1_embed_plot_samples": int,
    
    # Metrics
    "metrics_stage1": list,
    "metrics_stage2": list,
    
    # Multi-scale ensemble
    "use_multi_scale": bool,
    "multi_scale_kmers": list,
}


def get_type_name(expected_type: Union[Type, tuple]) -> str:
    """Get human-readable type name."""
    if isinstance(expected_type, tuple):
        names = []
        for t in expected_type:
            if t is type(None):
                names.append("null")
            else:
                names.append(t.__name__)
        return " | ".join(names)
    return expected_type.__name__


def validate_type(key: str, value: Any, expected_type: Union[Type, tuple]) -> bool:
    """Validate that value matches expected type(s)."""
    if isinstance(expected_type, tuple):
        return isinstance(value, expected_type)
    return isinstance(value, expected_type)


def validate_config(config_path: Path) -> Dict[str, Any]:
    """
    Validate config file against schema.
    
    Args:
        config_path: Path to YAML config file
        
    Returns:
        Parsed config dict if valid
        
    Raises:
        SystemExit: On any validation error
    """
    # Check file exists
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    
    # Parse YAML
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        print(f"ERROR: Failed to parse YAML: {e}", file=sys.stderr)
        sys.exit(1)
    
    if config is None:
        print("ERROR: Config file is empty", file=sys.stderr)
        sys.exit(1)
    
    if not isinstance(config, dict):
        print(f"ERROR: Config must be a dict, got {type(config).__name__}", file=sys.stderr)
        sys.exit(1)
    
    errors: List[str] = []
    
    # Check for unknown keys
    unknown_keys = set(config.keys()) - set(CONFIG_SCHEMA.keys())
    if unknown_keys:
        for key in sorted(unknown_keys):
            errors.append(f"Unknown key: '{key}'")
    
    # Check for missing keys
    missing_keys = set(CONFIG_SCHEMA.keys()) - set(config.keys())
    if missing_keys:
        for key in sorted(missing_keys):
            errors.append(f"Missing required key: '{key}'")
    
    # Validate types for existing keys
    for key, value in config.items():
        if key not in CONFIG_SCHEMA:
            continue  # Already reported as unknown
        
        expected_type = CONFIG_SCHEMA[key]
        if not validate_type(key, value, expected_type):
            actual_type = type(value).__name__ if value is not None else "null"
            expected_name = get_type_name(expected_type)
            errors.append(f"Type mismatch for '{key}': expected {expected_name}, got {actual_type}")
    
    # Check for duplicate keys by re-parsing with custom loader
    try:
        with open(config_path, 'r') as f:
            content = f.read()
        
        # Check for obvious concatenation errors in raw content
        lines = content.split('\n')
        for i, line in enumerate(lines, 1):
            # Skip comments and empty lines
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            
            # Check for multiple colons that might indicate concatenated keys
            if stripped.count(':') > 2:
                errors.append(f"Line {i}: Possible concatenated keys: '{stripped[:60]}...'")
            
            # Check for value immediately followed by another key (no newline)
            # Pattern like "true" or number followed directly by text
            import re
            concat_pattern = r':\s*(true|false|\d+\.?\d*)\s*[a-z_]+:'
            if re.search(concat_pattern, stripped, re.IGNORECASE):
                errors.append(f"Line {i}: Possible concatenated key-value: '{stripped[:60]}'")
    except Exception as e:
        errors.append(f"Failed to check for duplicates: {e}")
    
    if errors:
        print("=" * 60, file=sys.stderr)
        print("CONFIG VALIDATION FAILED", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        for error in sorted(errors):
            print(f"  ✗ {error}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        sys.exit(1)
    
    return config


def print_config(config: Dict[str, Any]) -> None:
    """Print config with sorted keys and type annotations."""
    print("=" * 60)
    print("CONFIG VALIDATION PASSED")
    print("=" * 60)
    print()
    
    # Find max key length for alignment
    max_key_len = max(len(k) for k in config.keys())
    
    for key in sorted(config.keys()):
        value = config[key]
        value_type = type(value).__name__ if value is not None else "NoneType"
        
        # Format value for display
        if isinstance(value, str):
            display_value = f'"{value}"'
        elif isinstance(value, list):
            display_value = str(value)
        elif value is None:
            display_value = "null"
        else:
            display_value = str(value)
        
        print(f"  {key:<{max_key_len}} : {display_value:<40} ({value_type})")
    
    print()
    print("=" * 60)
    print(f"Total keys: {len(config)}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Validate config.yaml against strict schema"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to config file (default: config.yaml)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only output errors, no success message"
    )
    
    args = parser.parse_args()
    
    # Validate config
    config = validate_config(args.config)
    
    # Print validated config
    if not args.quiet:
        print_config(config)
    else:
        print("✓ Config is valid")


if __name__ == "__main__":
    main()
