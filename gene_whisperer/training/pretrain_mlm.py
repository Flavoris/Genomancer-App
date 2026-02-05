"""Pretrain a DNA MLM model using genome sequences."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required (pip install pyyaml).") from exc

from gene_whisperer.datasets.fasta import load_fasta_sequences
from gene_whisperer.datasets.mlm_dataset import MLMDataset
from gene_whisperer.models.mlm_model import DNAMLMModel, MLMConfig
from gene_whisperer.models.transformer import TransformerConfig
from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.train_utils import TrainingConfig, build_optimizer, set_seed


@dataclass
class PretrainConfig:
    fasta_paths: List[Path]
    tokenizer_path: Path
    vocab_size: int
    max_length: int
    window_size: int
    mask_prob: float
    batch_size: int
    epochs: int
    lr: float
    weight_decay: float
    grad_accum_steps: int
    embedding_dim: int
    num_layers: int
    num_heads: int
    ff_dim: int
    dropout: float
    seed: int
    output_dir: Path


def _load_config(path: Path) -> PretrainConfig:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    mlm = data["mlm"]
    model = data["model"]
    training = data["training"]
    return PretrainConfig(
        fasta_paths=[Path(p) for p in mlm["fasta_paths"]],
        tokenizer_path=Path(mlm["tokenizer_path"]),
        vocab_size=int(mlm["vocab_size"]),
        max_length=int(mlm["max_length"]),
        window_size=int(mlm["window_size"]),
        mask_prob=float(mlm.get("mask_prob", 0.15)),
        batch_size=int(training["batch_size"]),
        epochs=int(training["epochs"]),
        lr=float(training["lr"]),
        weight_decay=float(training["weight_decay"]),
        grad_accum_steps=int(training.get("grad_accum_steps", 1)),
        embedding_dim=int(model["embedding_dim"]),
        num_layers=int(model["num_layers"]),
        num_heads=int(model["num_heads"]),
        ff_dim=int(model["ff_dim"]),
        dropout=float(model.get("dropout", 0.1)),
        seed=int(training.get("seed", 1337)),
        output_dir=Path(training["output_dir"]),
    )


def _maybe_train_tokenizer(config: PretrainConfig, sequences: List[str]) -> BPETokenizer:
    if config.tokenizer_path.exists():
        return BPETokenizer.load(config.tokenizer_path)
    tokenizer = BPETokenizer.train(sequences, vocab_size=config.vocab_size)
    tokenizer.save(config.tokenizer_path)
    return tokenizer


def train(config: PretrainConfig) -> None:
    set_seed(config.seed)
    sequences = load_fasta_sequences(config.fasta_paths)
    if not sequences:
        raise ValueError("No genome sequences loaded for MLM")

    tokenizer = _maybe_train_tokenizer(config, sequences)

    dataset = MLMDataset(
        sequences=sequences,
        tokenizer=tokenizer,
        window_size=config.window_size,
        max_length=config.max_length,
        mask_prob=config.mask_prob,
        seed=config.seed,
    )
    loader = DataLoader(dataset, batch_size=config.batch_size, shuffle=True)

    transformer_cfg = TransformerConfig(
        vocab_size=len(tokenizer.vocab),
        embedding_dim=config.embedding_dim,
        num_layers=config.num_layers,
        num_heads=config.num_heads,
        ff_dim=config.ff_dim,
        max_seq_len=config.max_length,
        dropout=config.dropout,
        pad_token_id=tokenizer.pad_token_id,
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = DNAMLMModel(MLMConfig(transformer=transformer_cfg)).to(device)

    optimizer = build_optimizer(
        model, TrainingConfig(lr=config.lr, weight_decay=config.weight_decay)
    )
    loss_fn = nn.CrossEntropyLoss(ignore_index=-100)

    config.output_dir.mkdir(parents=True, exist_ok=True)

    step = 0
    for epoch in range(config.epochs):
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()
        for batch in loader:
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)
            labels = batch["labels"].to(device)
            logits = model(input_ids, attention_mask)
            loss = loss_fn(
                logits.view(-1, logits.size(-1)), labels.view(-1)
            )
            loss = loss / config.grad_accum_steps
            loss.backward()
            total_loss += loss.item()
            if (step + 1) % config.grad_accum_steps == 0:
                optimizer.step()
                optimizer.zero_grad()
            step += 1
        avg_loss = total_loss / max(1, len(loader))
        print(f"Epoch {epoch + 1}/{config.epochs} - loss: {avg_loss:.4f}")

        checkpoint = {
            "model_state": model.state_dict(),
            "epoch": epoch + 1,
            "loss": avg_loss,
            "config": config.__dict__,
        }
        torch.save(checkpoint, config.output_dir / f"mlm_epoch_{epoch + 1}.pt")


def main() -> None:
    parser = argparse.ArgumentParser(description="Pretrain gene_whisperer MLM")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs" / "pretrain.yaml",
        help="Path to pretraining config",
    )
    args = parser.parse_args()
    config = _load_config(args.config)
    train(config)


if __name__ == "__main__":
    main()
