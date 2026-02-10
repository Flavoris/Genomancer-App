"""Pretrain a DNA MLM model using genome sequences."""
from __future__ import annotations

import argparse
import random
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List

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
    max_bases_per_file: int | None
    tokenizer_max_bases: int | None
    tokenizer_max_sequences: int | None
    tokenizer_window_size: int
    batch_size: int
    epochs: int
    samples_per_epoch: int
    log_interval: int
    num_workers: int
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


def _resolve_path(path_str: str, base_dir: Path) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _load_config(path: Path) -> PretrainConfig:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    config_dir = path.parent
    mlm = data["mlm"]
    model = data["model"]
    training = data["training"]

    max_bases_per_file = mlm.get("max_bases_per_file")
    if max_bases_per_file is not None:
        max_bases_per_file = int(max_bases_per_file)
    tokenizer_max_bases = mlm.get("tokenizer_max_bases")
    if tokenizer_max_bases is not None:
        tokenizer_max_bases = int(tokenizer_max_bases)
    tokenizer_max_sequences = mlm.get("tokenizer_max_sequences")
    if tokenizer_max_sequences is not None:
        tokenizer_max_sequences = int(tokenizer_max_sequences)

    return PretrainConfig(
        fasta_paths=[_resolve_path(str(p), config_dir) for p in mlm["fasta_paths"]],
        tokenizer_path=_resolve_path(str(mlm["tokenizer_path"]), config_dir),
        vocab_size=int(mlm["vocab_size"]),
        max_length=int(mlm["max_length"]),
        window_size=int(mlm["window_size"]),
        mask_prob=float(mlm.get("mask_prob", 0.15)),
        max_bases_per_file=max_bases_per_file,
        tokenizer_max_bases=tokenizer_max_bases,
        tokenizer_max_sequences=tokenizer_max_sequences,
        tokenizer_window_size=max(64, int(mlm.get("tokenizer_window_size", 2048))),
        batch_size=int(training["batch_size"]),
        epochs=int(training["epochs"]),
        samples_per_epoch=int(training.get("samples_per_epoch", 4096)),
        log_interval=max(1, int(training.get("log_interval", 50))),
        num_workers=max(0, int(training.get("num_workers", 0))),
        lr=float(training["lr"]),
        weight_decay=float(training["weight_decay"]),
        grad_accum_steps=int(training.get("grad_accum_steps", 1)),
        embedding_dim=int(model["embedding_dim"]),
        num_layers=int(model["num_layers"]),
        num_heads=int(model["num_heads"]),
        ff_dim=int(model["ff_dim"]),
        dropout=float(model.get("dropout", 0.1)),
        seed=int(training.get("seed", 1337)),
        output_dir=_resolve_path(str(training["output_dir"]), config_dir),
    )


def _sample_tokenizer_corpus(
    sequences: List[str],
    seed: int,
    max_bases: int | None,
    max_sequences: int | None,
    window_size: int,
) -> List[str]:
    if max_bases is None and max_sequences is None:
        return sequences

    if not sequences:
        return []

    rng = random.Random(seed)
    indices = list(range(len(sequences)))
    rng.shuffle(indices)

    sampled: List[str] = []
    consumed_bases = 0
    sampled_sequences = 0

    for idx in indices:
        if max_sequences is not None and sampled_sequences >= max_sequences:
            break
        if max_bases is not None and consumed_bases >= max_bases:
            break

        seq = sequences[idx]
        if not seq:
            continue

        span = min(window_size, len(seq))
        if span <= 0:
            continue
        if len(seq) == span:
            window = seq
        else:
            start = rng.randint(0, len(seq) - span)
            window = seq[start : start + span]

        if max_bases is not None:
            remaining = max_bases - consumed_bases
            if remaining <= 0:
                break
            if len(window) > remaining:
                window = window[:remaining]

        if not window:
            continue

        sampled.append(window)
        consumed_bases += len(window)
        sampled_sequences += 1

    return sampled


def _maybe_train_tokenizer(config: PretrainConfig, sequences: List[str]) -> BPETokenizer:
    if config.tokenizer_path.exists():
        print(f"Loading tokenizer: {config.tokenizer_path}", flush=True)
        return BPETokenizer.load(config.tokenizer_path)

    tokenizer_sequences = _sample_tokenizer_corpus(
        sequences=sequences,
        seed=config.seed,
        max_bases=config.tokenizer_max_bases,
        max_sequences=config.tokenizer_max_sequences,
        window_size=config.tokenizer_window_size,
    )
    if not tokenizer_sequences:
        raise ValueError("Tokenizer corpus sampling produced zero sequences.")

    tokenizer_bases = sum(len(seq) for seq in tokenizer_sequences)
    print(
        "Training tokenizer with "
        f"vocab_size={config.vocab_size} "
        f"on {len(tokenizer_sequences)} sampled sequences ({tokenizer_bases:,} bases)",
        flush=True,
    )
    tokenizer = BPETokenizer.train(
        tokenizer_sequences,
        vocab_size=config.vocab_size,
        verbose=True,
        log_interval=100,
    )
    tokenizer.save(config.tokenizer_path)
    print(f"Saved tokenizer: {config.tokenizer_path}", flush=True)
    return tokenizer


def _print_device_info(device: torch.device) -> None:
    print(f"Using device: {device}", flush=True)
    if device.type == "cuda":
        gpu_name = torch.cuda.get_device_name(0)
        print(f"GPU: {gpu_name}", flush=True)


def train(config: PretrainConfig) -> None:
    set_seed(config.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _print_device_info(device)

    print("Loading genome sequences...", flush=True)
    load_start = time.time()
    sequences = load_fasta_sequences(
        config.fasta_paths,
        max_bases_per_file=config.max_bases_per_file,
        verbose=True,
    )
    if not sequences:
        raise ValueError("No genome sequences loaded for MLM")

    total_bases = sum(len(sequence) for sequence in sequences)
    print(
        f"Loaded {len(sequences)} sequences ({total_bases:,} bases) in {time.time() - load_start:.1f}s",
        flush=True,
    )

    tokenizer = _maybe_train_tokenizer(config, sequences)

    dataset = MLMDataset(
        sequences=sequences,
        tokenizer=tokenizer,
        window_size=config.window_size,
        max_length=config.max_length,
        mask_prob=config.mask_prob,
        num_samples=config.samples_per_epoch,
        seed=config.seed,
    )

    pin_memory = device.type == "cuda"
    loader = DataLoader(
        dataset,
        batch_size=config.batch_size,
        shuffle=True,
        num_workers=config.num_workers,
        pin_memory=pin_memory,
        persistent_workers=config.num_workers > 0,
    )

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
    model = DNAMLMModel(MLMConfig(transformer=transformer_cfg)).to(device)

    optimizer = build_optimizer(
        model,
        TrainingConfig(lr=config.lr, weight_decay=config.weight_decay),
    )
    loss_fn = nn.CrossEntropyLoss(ignore_index=-100)

    config.output_dir.mkdir(parents=True, exist_ok=True)

    batches_per_epoch = len(loader)
    print(
        "Training setup: "
        f"epochs={config.epochs} batches_per_epoch={batches_per_epoch} "
        f"batch_size={config.batch_size} grad_accum={config.grad_accum_steps}",
        flush=True,
    )

    global_step = 0
    for epoch in range(1, config.epochs + 1):
        epoch_start = time.time()
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()

        for batch_idx, batch in enumerate(loader, start=1):
            input_ids = batch["input_ids"].to(device, non_blocking=pin_memory)
            attention_mask = batch["attention_mask"].to(device, non_blocking=pin_memory)
            labels = batch["labels"].to(device, non_blocking=pin_memory)

            logits = model(input_ids, attention_mask)
            loss = loss_fn(logits.view(-1, logits.size(-1)), labels.view(-1))
            scaled_loss = loss / config.grad_accum_steps
            scaled_loss.backward()
            total_loss += loss.item()

            if batch_idx % config.grad_accum_steps == 0:
                optimizer.step()
                optimizer.zero_grad()
                global_step += 1

            if batch_idx % config.log_interval == 0 or batch_idx == batches_per_epoch:
                print(
                    f"Epoch {epoch}/{config.epochs} "
                    f"batch {batch_idx}/{batches_per_epoch} "
                    f"loss={loss.item():.4f}",
                    flush=True,
                )

        if batches_per_epoch % config.grad_accum_steps != 0:
            optimizer.step()
            optimizer.zero_grad()
            global_step += 1

        avg_loss = total_loss / max(1, batches_per_epoch)
        elapsed = time.time() - epoch_start
        print(
            f"Epoch {epoch}/{config.epochs} complete "
            f"avg_loss={avg_loss:.4f} epoch_time={elapsed:.1f}s global_step={global_step}",
            flush=True,
        )

        checkpoint = {
            "model_state": model.state_dict(),
            "epoch": epoch,
            "loss": avg_loss,
            "global_step": global_step,
            "config": config.__dict__,
        }
        checkpoint_path = config.output_dir / f"mlm_epoch_{epoch}.pt"
        torch.save(checkpoint, checkpoint_path)
        print(f"Saved checkpoint: {checkpoint_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pretrain gene_whisperer MLM")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs" / "pretrain.yaml",
        help="Path to pretraining config",
    )
    args = parser.parse_args()
    config = _load_config(args.config.resolve())
    train(config)


if __name__ == "__main__":
    main()
