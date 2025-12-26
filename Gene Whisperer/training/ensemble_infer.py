"""Promoter prediction via Stage-1 ensemble over multiple k-mer models."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

import torch
import yaml

from dataset import KmerVocabulary, compute_tnc, compute_pseeiip
from model import GeneWhispererStage1, MultiScaleEnsemble

LOGGER = logging.getLogger("gene_whisperer.ensemble")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def load_checkpoint(model: GeneWhispererStage1, path: Path, device: torch.device) -> None:
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    state_dict = checkpoint.get("model_state", checkpoint)

    # Use strict=True to ensure all weights load correctly
    # Log any issues for debugging
    try:
        incompatible = model.load_state_dict(state_dict, strict=True)
        if incompatible.missing_keys:
            LOGGER.warning("Missing keys: %s", incompatible.missing_keys)
        if incompatible.unexpected_keys:
            LOGGER.warning("Unexpected keys: %s", incompatible.unexpected_keys)
    except RuntimeError as e:
        LOGGER.error("Failed to load checkpoint with strict=True: %s", e)
        LOGGER.info("Attempting load with strict=False...")
        incompatible = model.load_state_dict(state_dict, strict=False)
        if incompatible.missing_keys:
            LOGGER.warning("Missing keys (not loaded): %s", incompatible.missing_keys[:10])
        if incompatible.unexpected_keys:
            LOGGER.warning("Unexpected keys (ignored): %s", incompatible.unexpected_keys[:10])

    model.to(device)
    model.eval()
    LOGGER.info("Loaded model from %s", path)


def load_vocab(vocab_path: Path) -> KmerVocabulary:
    return KmerVocabulary.load(vocab_path)


def build_model(cfg: Dict, vocab: KmerVocabulary, device: torch.device) -> GeneWhispererStage1:
    """Build GeneWhispererStage1 model with V3 architecture parameters."""
    embedding_dim = int(cfg.get("embedding_dim", 192))
    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.2))
    
    # TCN parameters
    use_tcn = bool(cfg.get("use_tcn", True))
    tcn_hidden = int(cfg.get("tcn_hidden", 256))
    tcn_levels = int(cfg.get("tcn_levels", 4))
    tcn_kernel = int(cfg.get("tcn_kernel", 3))
    multiscale_channels = int(cfg.get("multiscale_channels", 64))
    multiscale_kernels = tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15]))
    lstm_hidden = int(cfg.get("lstm_hidden", 192))
    
    # New V3 architecture parameters
    post_cnn_transformer_layers = int(cfg.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg.get("fusion_hidden", 256))
    
    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        engineered_dim=int(cfg.get("engineered_dim", 128)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        # TCN parameters
        use_tcn=use_tcn,
        tcn_hidden=tcn_hidden,
        tcn_levels=tcn_levels,
        tcn_kernel=tcn_kernel,
        multiscale_channels=multiscale_channels,
        multiscale_kernels=multiscale_kernels,
        lstm_hidden=lstm_hidden,
        # V3 architecture parameters
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        fusion_hidden=fusion_hidden,
    ).to(device)
    return model


def compute_engineered_features(sequence: str, cfg: Optional[Dict] = None) -> torch.Tensor:
    """Compute engineered features: TNC(64) + PseEIIP(64) = 128."""
    tnc = compute_tnc(sequence)          # 64-dim
    pseeiip = compute_pseeiip(sequence)  # 64-dim
    if cfg is not None:
        if not bool(cfg.get("stage1_feature_enable_tnc", True)):
            tnc = torch.zeros_like(tnc)
        if not bool(cfg.get("stage1_feature_enable_pseeiip", True)):
            pseeiip = torch.zeros_like(pseeiip)
    return torch.cat([tnc, pseeiip], dim=0)


def prepare_inputs(
    sequence: str, vocab: KmerVocabulary, max_bp_len: int, cfg: Optional[Dict] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    tokens = vocab.tokenize(sequence, max_bp_len)
    engineered = compute_engineered_features(sequence, cfg)
    return tokens.unsqueeze(0), engineered.unsqueeze(0)


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage-1 ensemble inference")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--sequence", required=True, help="Raw DNA sequence (A/C/G/T)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Promoter decision threshold")
    parser.add_argument("--kmer_configs", type=str, help="JSON mapping of k to checkpoint/vocab paths")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    sequence = args.sequence.strip().upper()
    max_bp_len = int(cfg.get("max_bp_len", 81))
    device = select_device()

    if args.kmer_configs:
        kmer_config = json.loads(args.kmer_configs)
    else:
        # fallback to default naming convention
        kmer_config = {
            "3": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k3.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k3_vocab.json").resolve()),
            },
            "4": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k4.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k4_vocab.json").resolve()),
            },
            "6": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k6.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k6_vocab.json").resolve()),
            },
        }

    models = []
    batch_inputs = {}
    for k_str, paths in kmer_config.items():
        k = int(k_str)
        vocab_path = Path(paths["vocab"]).resolve()
        checkpoint_path = Path(paths["checkpoint"]).resolve()
        if not vocab_path.exists() or not checkpoint_path.exists():
            LOGGER.warning("Skipping k=%s (missing vocab or checkpoint)", k_str)
            continue
        vocab = load_vocab(vocab_path)
        model = build_model(cfg, vocab, device)
        load_checkpoint(model, checkpoint_path, device)
        tokens, engineered = prepare_inputs(sequence, vocab, max_bp_len, cfg)
        tokens = tokens.to(device)
        engineered = engineered.to(device)
        models.append(model)
        batch_inputs[k] = (tokens, engineered)

    if not models:
        raise RuntimeError("No models loaded for ensemble inference")

    ensemble = MultiScaleEnsemble(models).to(device)
    with torch.no_grad():
        outputs = ensemble(batch_inputs)
    prob = outputs.squeeze().item()
    label = "promoter" if prob >= args.threshold else "non-promoter"
    print(json.dumps({"probability": prob, "label": label}))


if __name__ == "__main__":
    main()
