"""Pretrain a DNA MLM model using genome sequences."""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import torch
from torch import nn
from torch.utils.data import DataLoader

from gene_whisperer.datasets.fasta import load_fasta_sequences
from gene_whisperer.datasets.mlm_dataset import MLMDataset
from gene_whisperer.models.mlm_model import DNAMLMModel, MLMConfig
from gene_whisperer.models.transformer import TransformerConfig
from gene_whisperer.training.pretrain_config import (
    PretrainConfig,
    load_pretrain_config,
)
from gene_whisperer.training.pretrain_optim import (
    build_grad_scaler,
    build_warmup_cosine_scheduler,
    is_amp_enabled,
)
from gene_whisperer.training.pretrain_runtime import (
    build_checkpoint,
    maybe_train_tokenizer,
    print_device_info,
)
from gene_whisperer.training.train_utils import TrainingConfig, build_optimizer, set_seed


def _is_improvement(current_loss: float, best_loss: float, min_delta: float) -> bool:
    return current_loss < (best_loss - min_delta)


def train(config: PretrainConfig) -> None:
    set_seed(config.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print_device_info(device)

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
        f"Loaded {len(sequences)} sequences ({total_bases:,} bases) "
        f"in {time.time() - load_start:.1f}s",
        flush=True,
    )

    tokenizer = maybe_train_tokenizer(config, sequences)

    dataset = MLMDataset(
        sequences=sequences,
        tokenizer=tokenizer,
        window_size=config.window_size,
        max_length=config.max_length,
        mask_prob=config.mask_prob,
        num_samples=config.samples_per_epoch,
        seed=config.seed,
        sample_by_length=config.sample_by_length,
        mask_ambiguous_tokens=config.mask_ambiguous_tokens,
        min_masked_tokens=config.min_masked_tokens,
        min_maskable_tokens=config.min_maskable_tokens,
        min_tokenized_tokens=config.min_tokenized_tokens,
        resample_attempts=config.resample_attempts,
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
    model = DNAMLMModel(
        MLMConfig(transformer=transformer_cfg, head_dropout=config.dropout)
    ).to(device)

    optimizer = build_optimizer(
        model,
        TrainingConfig(lr=config.lr, weight_decay=config.weight_decay),
    )
    steps_per_epoch = math.ceil(len(loader) / config.grad_accum_steps)
    total_steps = max(1, steps_per_epoch * config.epochs)
    scheduler = build_warmup_cosine_scheduler(
        optimizer=optimizer,
        total_steps=total_steps,
        warmup_ratio=config.warmup_ratio,
        min_lr_ratio=config.min_lr_ratio,
    )
    amp_enabled = is_amp_enabled(config.use_amp, device)
    scaler = build_grad_scaler(use_amp=config.use_amp, device=device)
    loss_fn = nn.CrossEntropyLoss(ignore_index=-100)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    best_checkpoint_path = config.output_dir / "mlm_best.pt"

    batches_per_epoch = len(loader)
    print(
        "Training setup: "
        f"max_epochs={config.epochs} batches_per_epoch={batches_per_epoch} "
        f"batch_size={config.batch_size} grad_accum={config.grad_accum_steps} "
        f"lr={config.lr:.2e} warmup_ratio={config.warmup_ratio:.3f} "
        f"min_lr_ratio={config.min_lr_ratio:.3f} amp={amp_enabled} "
        f"mask_ambiguous={config.mask_ambiguous_tokens} "
        f"min_masked_tokens={config.min_masked_tokens} "
        f"min_maskable_tokens={config.min_maskable_tokens} "
        f"min_tokenized_tokens={config.min_tokenized_tokens} "
        f"resample_attempts={config.resample_attempts} "
        f"min_epochs={config.min_epochs} "
        f"patience={config.early_stopping_patience} "
        f"min_delta={config.early_stopping_min_delta}",
        flush=True,
    )

    global_step = 0
    best_loss = math.inf
    best_epoch = 0
    no_improve_epochs = 0

    for epoch in range(1, config.epochs + 1):
        epoch_start = time.time()
        model.train()
        total_loss = 0.0
        trained_batches = 0
        skipped_batches = 0
        supervised_tokens = 0
        accum_batches = 0
        optimizer.zero_grad()

        for batch_idx, batch in enumerate(loader, start=1):
            input_ids = batch["input_ids"].to(device, non_blocking=pin_memory)
            attention_mask = batch["attention_mask"].to(device, non_blocking=pin_memory)
            labels = batch["labels"].to(device, non_blocking=pin_memory)
            valid_targets = int((labels != -100).sum().item())
            if valid_targets == 0:
                skipped_batches += 1
                continue

            with torch.autocast(
                device_type=device.type,
                dtype=torch.float16,
                enabled=amp_enabled,
            ):
                logits = model(input_ids, attention_mask)
                loss = loss_fn(logits.view(-1, logits.size(-1)), labels.view(-1))
            scaled_loss = loss / config.grad_accum_steps
            scaler.scale(scaled_loss).backward()
            total_loss += loss.item() * valid_targets
            trained_batches += 1
            supervised_tokens += valid_targets
            accum_batches += 1

            if accum_batches >= config.grad_accum_steps:
                if config.max_grad_norm > 0:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(
                        model.parameters(),
                        config.max_grad_norm,
                    )
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()
                scheduler.step()
                global_step += 1
                accum_batches = 0

            if batch_idx % config.log_interval == 0 or batch_idx == batches_per_epoch:
                current_lr = optimizer.param_groups[0]["lr"]
                print(
                    f"Epoch {epoch}/{config.epochs} "
                    f"batch {batch_idx}/{batches_per_epoch} "
                    f"loss={loss.item():.4f} targets={valid_targets} lr={current_lr:.2e}",
                    flush=True,
                )

        if accum_batches > 0:
            if config.max_grad_norm > 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(
                    model.parameters(),
                    config.max_grad_norm,
                )
            scaler.step(optimizer)
            scaler.update()
            optimizer.zero_grad()
            scheduler.step()
            global_step += 1

        if trained_batches == 0:
            raise RuntimeError(
                "No supervised MLM targets were found in this epoch. "
                "Lower min_maskable_tokens, lower min_tokenized_tokens, "
                "or enable mask_ambiguous_tokens."
            )
        avg_loss = total_loss / supervised_tokens
        elapsed = time.time() - epoch_start
        print(
            f"Epoch {epoch}/{config.epochs} complete "
            f"avg_loss={avg_loss:.4f} epoch_time={elapsed:.1f}s "
            f"global_step={global_step} trained_batches={trained_batches} "
            f"skipped_batches={skipped_batches} supervised_tokens={supervised_tokens}",
            flush=True,
        )

        if _is_improvement(avg_loss, best_loss, config.early_stopping_min_delta):
            best_loss = avg_loss
            best_epoch = epoch
            no_improve_epochs = 0
            checkpoint = build_checkpoint(
                model=model,
                config=config,
                epoch=epoch,
                loss=avg_loss,
                global_step=global_step,
                best_loss=best_loss,
                best_epoch=best_epoch,
            )
            torch.save(checkpoint, best_checkpoint_path)
            print(
                f"Saved best checkpoint: {best_checkpoint_path} (loss={best_loss:.4f})",
                flush=True,
            )
        else:
            no_improve_epochs += 1
            print(
                f"No improvement for {no_improve_epochs}/{config.early_stopping_patience} "
                f"epochs (best_loss={best_loss:.4f} at epoch {best_epoch})",
                flush=True,
            )

        if not config.save_best_only:
            checkpoint = build_checkpoint(
                model=model,
                config=config,
                epoch=epoch,
                loss=avg_loss,
                global_step=global_step,
                best_loss=best_loss,
                best_epoch=best_epoch,
            )
            epoch_checkpoint = config.output_dir / f"mlm_epoch_{epoch}.pt"
            torch.save(checkpoint, epoch_checkpoint)
            print(f"Saved epoch checkpoint: {epoch_checkpoint}", flush=True)

        if (
            epoch >= config.min_epochs
            and no_improve_epochs >= config.early_stopping_patience
        ):
            print(
                f"Early stopping triggered at epoch {epoch}. "
                f"Best epoch={best_epoch} best_loss={best_loss:.4f}",
                flush=True,
            )
            break

    print(
        f"Training complete. Best epoch={best_epoch} best_loss={best_loss:.4f} "
        f"checkpoint={best_checkpoint_path}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Pretrain gene_whisperer MLM")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs" / "pretrain.yaml",
        help="Path to pretraining config",
    )
    args = parser.parse_args()
    config = load_pretrain_config(args.config.resolve())
    train(config)


if __name__ == "__main__":
    main()
