"""Runtime helpers for MLM pretraining."""
from __future__ import annotations

from typing import List

import torch

from gene_whisperer.models.mlm_model import DNAMLMModel
from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.pretrain_config import (
    PretrainConfig,
    sample_tokenizer_corpus,
)


def build_checkpoint(
    model: DNAMLMModel,
    config: PretrainConfig,
    epoch: int,
    loss: float,
    global_step: int,
    best_loss: float,
    best_epoch: int,
) -> dict:
    return {
        "model_state": model.state_dict(),
        "epoch": epoch,
        "loss": loss,
        "global_step": global_step,
        "best_loss": best_loss,
        "best_epoch": best_epoch,
        "config": config.__dict__,
    }


def maybe_train_tokenizer(config: PretrainConfig, sequences: List[str]) -> BPETokenizer:
    if config.tokenizer_path.exists():
        print(f"Loading tokenizer: {config.tokenizer_path}", flush=True)
        return BPETokenizer.load(config.tokenizer_path)

    tokenizer_sequences = sample_tokenizer_corpus(
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


def print_device_info(device: torch.device) -> None:
    print(f"Using device: {device}", flush=True)
    if device.type == "cuda":
        gpu_name = torch.cuda.get_device_name(0)
        print(f"GPU: {gpu_name}", flush=True)
