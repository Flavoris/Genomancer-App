"""Promoter prediction via Stage-1 ensemble over multiple k-mer models."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Tuple

import torch
import yaml

from dataset import KmerVocabulary, compute_eiip, compute_pstnp, compute_tnc
from model import GeneWhispererStage1, Stage1Ensemble

LOGGER = logging.getLogger("gene_whisperer.ensemble")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def load_checkpoint(model: GeneWhispererStage1, path: Path, device: torch.device) -> None:
    checkpoint = torch.load(path, map_location=device)
    state_dict = checkpoint.get("model_state", checkpoint)
    model.load_state_dict(state_dict, strict=False)
    model.to(device)
    model.eval()
    LOGGER.info("Loaded model from %s", path)


def load_vocab(vocab_path: Path) -> KmerVocabulary:
    return KmerVocabulary.load(vocab_path)


def build_model(cfg: Dict, vocab: KmerVocabulary, device: torch.device) -> GeneWhispererStage1:
    embedding_dim = int(cfg.get("embedding_dim", 160))
    transformer_layers = int(cfg.get("transformer_layers", 4))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 256))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.1))
    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        use_positional_bias=bool(cfg.get("use_positional_bias", True)),
        pad_token_id=vocab.pad_id,
        engineered_dim=int(cfg.get("stage1_engineered_dim", 192)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
    ).to(device)
    return model


def compute_engineered_features(sequence: str) -> torch.Tensor:
    tnc = compute_tnc(sequence)
    pstnp = compute_pstnp(sequence)
    eiip = compute_eiip(sequence)
    return torch.cat([tnc, pstnp, eiip], dim=0)


def prepare_inputs(sequence: str, vocab: KmerVocabulary, max_bp_len: int) -> Tuple[torch.Tensor, torch.Tensor]:
    tokens = vocab.tokenize(sequence, max_bp_len)
    engineered = compute_engineered_features(sequence)
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
        tokens, engineered = prepare_inputs(sequence, vocab, max_bp_len)
        tokens = tokens.to(device)
        engineered = engineered.to(device)
        models.append(model)
        batch_inputs[k] = (tokens, engineered)

    if not models:
        raise RuntimeError("No models loaded for ensemble inference")

    ensemble = Stage1Ensemble(models).to(device)
    with torch.no_grad():
        outputs = ensemble(batch_inputs)
    prob = outputs.squeeze().item()
    label = "promoter" if prob >= args.threshold else "non-promoter"
    print(json.dumps({"probability": prob, "label": label}))


if __name__ == "__main__":
    main()
