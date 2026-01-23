"""Core experiment logic for transformer contribution diagnosis.

Provides utilities to build models, load checkpoints, and run the three
ablation experiments (full model, transformer-only, features-only).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

from dataset import KmerVocabulary, PromoterDatasetStage1
from model import create_model_variant
from hooks import TransformerOutputRandomizer

LOGGER = logging.getLogger(__name__)


def build_model_from_config(config: dict) -> nn.Module:
    """Instantiate the correct model variant from config."""
    variant = str(config.get("model_variant", "v6"))
    kmer = int(config.get("kmer", 6))
    max_bp_len = int(config.get("max_bp_len", 81))

    model_config = {
        "vocab_size": int(config.get("vocab_size", 4099)),
        "embedding_dim": int(config.get("embedding_dim", 384)),
        "transformer_layers": int(config.get("transformer_layers", 12)),
        "transformer_heads": int(config.get("transformer_heads", 12)),
        "transformer_ff_dim": int(config.get("transformer_ff_dim", 1536)),
        "transformer_dropout": float(config.get("transformer_dropout", 0.1)),
        "pad_token_id": config.get("pad_token_id"),
        "kmer": kmer,
        "max_bp_len": max_bp_len,
        "max_seq_len": max_bp_len - kmer + 1,
        "engineered_dim": int(config.get("engineered_dim", 288)),
        "stage1_use_engineered_features": bool(
            config.get("stage1_use_engineered_features", True)
        ),
        # Simplified model settings
        "classifier_hidden": config.get("simplified_model", {}).get(
            "classifier_hidden", 256
        ),
        "classifier_dropout": config.get("simplified_model", {}).get(
            "classifier_dropout", 0.2
        ),
        "drop_path_rate": float(config.get("drop_path_rate", 0.1)),
        # Modern architecture settings
        "use_rope": bool(config.get("use_rope", True)),
        "rope_base": float(config.get("rope_base", 10000.0)),
        "ffn_type": config.get("ffn_type", "swiglu"),
        "norm_type": config.get("norm_type", "rmsnorm"),
        "ffn_mult": config.get("ffn_mult", 2.67),
        "use_positional_motif_bias": bool(
            config.get("use_positional_motif_bias", False)
        ),
        "motif_regions": config.get("motif_regions"),
        # V6 specific
        "use_swiglu": bool(config.get("use_swiglu", True)),
        "use_rmsnorm": bool(config.get("use_rmsnorm", True)),
        "use_positional_bias": bool(config.get("use_positional_bias", True)),
        "use_gqa": bool(config.get("use_gqa", False)),
        "num_kv_heads": config.get("num_kv_heads"),
        # Engineered feature MLP settings
        "engineered_mlp_hidden": int(
            config.get("engineered_mlp_hidden", 512)
        ),
        "engineered_mlp_output": int(
            config.get("engineered_mlp_output", 384)
        ),
    }
    LOGGER.info("Building model variant: %s", variant)
    return create_model_variant(variant, model_config)


def load_checkpoint(
    model: nn.Module,
    checkpoint_path: Path,
    device: torch.device,
    threshold_override: Optional[float] = None,
) -> float:
    """Load checkpoint weights and return the decision threshold."""
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)

    state_dict = checkpoint.get("model_state", checkpoint)
    model.load_state_dict(state_dict, strict=False)
    LOGGER.info("Loaded model weights from checkpoint")

    if threshold_override is not None:
        return threshold_override

    threshold = checkpoint.get("best_threshold", 0.5)
    return float(threshold)


def build_val_dataloader(
    config: dict,
    max_samples: Optional[int] = None,
) -> DataLoader:
    """Build the validation DataLoader from config."""
    val_path = Path(config.get("stage1_val", "../data/stage1_val.tsv"))
    if not val_path.exists():
        data_cfg = config.get("data", {})
        val_path = Path(data_cfg.get("stage1_val", "../data/stage1_val.tsv"))

    delimiter = config.get("delimiter", "\t")
    val_df = pd.read_csv(val_path, sep=delimiter)
    LOGGER.info("Loaded %d validation samples from %s", len(val_df), val_path)

    kmer = int(config.get("kmer", 6))
    vocab_path = _resolve_vocab_path(config, kmer)
    vocab = KmerVocabulary.load(Path(vocab_path))
    LOGGER.info("Loaded vocabulary (k=%d, size=%d)", kmer, len(vocab))

    max_bp_len = int(config.get("max_bp_len", 81))
    engineered_dim = int(config.get("engineered_dim", 288))

    dataset = PromoterDatasetStage1(
        df=val_df,
        max_bp_len=max_bp_len,
        vocab=vocab,
        use_engineered_features=True,
        feature_enable_tnc=bool(config.get("stage1_feature_enable_tnc", True)),
        feature_enable_pseeiip=bool(config.get("stage1_feature_enable_pseeiip", True)),
        feature_enable_cksnap=bool(config.get("stage1_feature_enable_cksnap", True)),
        feature_enable_pstnp=bool(config.get("stage1_feature_enable_pstnp", True)),
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.0,
        cache_engineered_features=True,
        cache_tokens=True,
    )

    if max_samples and max_samples < len(dataset):
        indices = list(range(max_samples))
        dataset = Subset(dataset, indices)

    batch_size = int(config.get("batch_size", 16))
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        pin_memory=False,
    )


def run_experiment_a(
    model: nn.Module,
    val_loader: DataLoader,
    device: torch.device,
    threshold: float,
) -> float:
    """Experiment A: Full model (transformer + engineered features)."""
    correct, total = _evaluate(model, val_loader, device, threshold, mode="full")
    accuracy = correct / total if total > 0 else 0.0
    print(f"Experiment A (Full model):        {accuracy * 100:.1f}% accuracy")
    return accuracy


def run_experiment_b(
    model: nn.Module,
    val_loader: DataLoader,
    device: torch.device,
    threshold: float,
) -> float:
    """Experiment B: Transformer-only (zero out engineered features)."""
    correct, total = _evaluate(
        model, val_loader, device, threshold, mode="zero_features"
    )
    accuracy = correct / total if total > 0 else 0.0
    print(f"Experiment B (Transformer only):  {accuracy * 100:.1f}% accuracy")
    return accuracy


def run_experiment_c(
    model: nn.Module,
    val_loader: DataLoader,
    device: torch.device,
    threshold: float,
) -> float:
    """Experiment C: Features-only (randomize transformer output)."""
    correct, total = _evaluate(
        model, val_loader, device, threshold, mode="random_transformer"
    )
    accuracy = correct / total if total > 0 else 0.0
    print(f"Experiment C (Features only):     {accuracy * 100:.1f}% accuracy")
    return accuracy


def print_analysis(acc_a: float, acc_b: float, acc_c: float) -> None:
    """Print interpretation of experiment results."""
    transformer_contribution = acc_a - acc_c
    features_contribution = acc_a - acc_b

    print("\nAnalysis:")
    print(f"  Transformer contributes: +{transformer_contribution * 100:.1f}% above features alone")
    print(f"  Engineered features contribute: +{features_contribution * 100:.1f}% above transformer alone")
    print()

    if acc_b <= 0.55:
        print("  ** The transformer is NOT learning useful patterns **")
        print("  -> Pre-training needs improvement")
    elif acc_b <= 0.65:
        print("  The transformer is learning WEAK patterns")
        print("  -> Pre-training may benefit from more epochs or data")
    else:
        print("  The transformer IS learning useful patterns")

    if acc_c >= acc_a * 0.98:
        print("  ** Engineered features do nearly ALL the work **")
    elif acc_c >= acc_a * 0.95:
        print("  Engineered features do MOST of the work")

    print()


def _evaluate(
    model: nn.Module,
    val_loader: DataLoader,
    device: torch.device,
    threshold: float,
    mode: str,
) -> Tuple[int, int]:
    """Run evaluation in a given experiment mode.

    Args:
        mode: One of "full", "zero_features", "random_transformer"

    Returns:
        (correct_count, total_count)
    """
    hook_handle = None
    if mode == "random_transformer":
        randomizer = TransformerOutputRandomizer(model)
        hook_handle = randomizer.attach()

    correct = 0
    total = 0

    with torch.no_grad():
        for tokens, engineered, labels in val_loader:
            tokens = tokens.to(device)
            engineered = engineered.to(device)
            labels = labels.to(device)

            if mode == "zero_features":
                engineered = torch.zeros_like(engineered)

            logits = model(tokens, engineered)
            probs = torch.sigmoid(logits).squeeze(-1)
            preds = (probs >= threshold).float()
            correct += (preds == labels).sum().item()
            total += labels.size(0)

    if hook_handle is not None:
        hook_handle.remove()

    return int(correct), int(total)


def _resolve_vocab_path(config: dict, kmer: int) -> str:
    """Find vocab path from config."""
    vocab_by_k = config.get("mlm_vocab_by_k", {})
    if not vocab_by_k:
        vocab_by_k = config.get("vocab_paths", {}).get("mlm_vocab_by_k", {})
    if kmer in vocab_by_k:
        return vocab_by_k[kmer]
    cache_dir = config.get("vocab_cache_dir", "../artifacts/vocabs")
    return f"{cache_dir}/mlm_k{kmer}_vocab.json"
