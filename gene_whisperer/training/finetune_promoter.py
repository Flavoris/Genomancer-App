"""Fine-tune promoter classifier with engineered features."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required (pip install pyyaml).") from exc

from gene_whisperer.datasets.promoter_dataset import PromoterDataset
from gene_whisperer.features.engineered import FeatureExtractor
from gene_whisperer.models.promoter_model import PromoterConfig, PromoterModel
from gene_whisperer.models.transformer import TransformerConfig
from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.train_utils import TrainingConfig, build_optimizer, set_seed
from gene_whisperer.utils.metrics import compute_binary_metrics


@dataclass
class FinetuneConfig:
    train_path: Path
    val_path: Path
    tokenizer_path: Path
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
    max_length: int
    use_llrd: bool
    layer_decay: float
    seed: int
    output_dir: Path
    use_pstnp: bool


def _load_config(path: Path, stage: str) -> FinetuneConfig:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    stage_cfg = data[stage]
    model = data["model"]
    training = data["training"]
    return FinetuneConfig(
        train_path=Path(stage_cfg["train_path"]),
        val_path=Path(stage_cfg["val_path"]),
        tokenizer_path=Path(stage_cfg["tokenizer_path"]),
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
        max_length=int(model.get("max_length", 128)),
        use_llrd=bool(training.get("use_llrd", False)),
        layer_decay=float(training.get("layer_decay", 0.9)),
        seed=int(training.get("seed", 1337)),
        output_dir=Path(training["output_dir"]),
        use_pstnp=bool(stage_cfg.get("use_pstnp", True)),
    )


def _load_sequences(path: Path, task: str) -> Tuple[List[str], List[int]]:
    sequences: List[str] = []
    labels: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            seq = row.get("sequence", "")
            if task == "stage1":
                label = int(row.get("is_promoter", 0))
            else:
                strength = row.get("strength", "weak").strip().lower()
                label = 1 if strength == "strong" else 0
            sequences.append(seq)
            labels.append(label)
    return sequences, labels


def _build_feature_extractor(config: FinetuneConfig, task: str) -> FeatureExtractor:
    extractor = FeatureExtractor(use_pstnp=config.use_pstnp)
    if config.use_pstnp:
        seqs, labels = _load_sequences(config.train_path, task)
        positives = [seq for seq, label in zip(seqs, labels) if label == 1]
        negatives = [seq for seq, label in zip(seqs, labels) if label == 0]
        extractor.fit_pstnp(positives, negatives)
    return extractor


def _save_pstnp(extractor: FeatureExtractor, output_dir: Path, task: str) -> None:
    if extractor.pstnp_matrix is None:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    data = extractor.pstnp_matrix.to_dict()
    with (output_dir / f"pstnp_{task}.json").open("w", encoding="utf-8") as handle:
        import json

        json.dump(data, handle, indent=2, sort_keys=True)


def train(config: FinetuneConfig, task: str) -> None:
    set_seed(config.seed)
    tokenizer = BPETokenizer.load(config.tokenizer_path)

    feature_extractor = _build_feature_extractor(config, task)
    _save_pstnp(feature_extractor, config.output_dir, task)

    train_dataset = PromoterDataset(
        config.train_path, tokenizer, feature_extractor, max_length=config.max_length, task=task
    )
    val_dataset = PromoterDataset(
        config.val_path, tokenizer, feature_extractor, max_length=config.max_length, task=task
    )

    train_loader = DataLoader(train_dataset, batch_size=config.batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=config.batch_size, shuffle=False)

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

    model = PromoterModel(
        PromoterConfig(
            transformer=transformer_cfg,
            engineered_dim=feature_extractor.feature_dim(),
            engineered_hidden=max(64, config.embedding_dim // 2),
            fusion_hidden=config.embedding_dim,
            dropout=config.dropout,
            use_engineered_features=True,
        )
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    optimizer = build_optimizer(
        model,
        TrainingConfig(
            lr=config.lr,
            weight_decay=config.weight_decay,
            layer_decay=config.layer_decay,
            use_llrd=config.use_llrd,
        ),
    )
    loss_fn = nn.BCEWithLogitsLoss()

    config.output_dir.mkdir(parents=True, exist_ok=True)

    step = 0
    for epoch in range(config.epochs):
        model.train()
        optimizer.zero_grad()
        running_loss = 0.0
        for batch in train_loader:
            logits = model(
                batch["input_ids"].to(device),
                attention_mask=batch["attention_mask"].to(device),
                engineered_features=batch["engineered"].to(device),
                task=task,
            ).squeeze(-1)
            labels = batch["label"].to(device)
            loss = loss_fn(logits, labels)
            loss = loss / config.grad_accum_steps
            loss.backward()
            running_loss += loss.item()
            if (step + 1) % config.grad_accum_steps == 0:
                optimizer.step()
                optimizer.zero_grad()
            step += 1

        avg_loss = running_loss / max(1, len(train_loader))
        metrics = evaluate(model, val_loader, device, task)
        print(
            f"Epoch {epoch + 1}/{config.epochs} - loss: {avg_loss:.4f} "
            f"val_acc: {metrics['accuracy']:.3f}"
        )

        checkpoint = {
            "model_state": model.state_dict(),
            "epoch": epoch + 1,
            "metrics": metrics,
            "config": config.__dict__,
        }
        torch.save(checkpoint, config.output_dir / f"{task}_epoch_{epoch + 1}.pt")


def evaluate(
    model: PromoterModel,
    loader: DataLoader,
    device: torch.device,
    task: str,
) -> dict:
    model.eval()
    all_logits: List[float] = []
    all_labels: List[float] = []
    with torch.no_grad():
        for batch in loader:
            logits = model(
                batch["input_ids"].to(device),
                attention_mask=batch["attention_mask"].to(device),
                engineered_features=batch["engineered"].to(device),
                task=task,
            ).squeeze(-1)
            all_logits.append(logits.cpu().numpy())
            all_labels.append(batch["label"].cpu().numpy())

    logits = np.concatenate(all_logits)
    labels = np.concatenate(all_labels)
    return compute_binary_metrics(logits, labels)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fine-tune promoter classifier")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs" / "finetune.yaml",
        help="Path to fine-tuning config",
    )
    parser.add_argument(
        "--task",
        choices=["stage1", "stage2"],
        default="stage1",
        help="Which classification stage to train",
    )
    args = parser.parse_args()
    config = _load_config(args.config, args.task)
    train(config, args.task)


if __name__ == "__main__":
    main()
