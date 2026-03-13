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


def _tokenizer_matches_config(
    tokenizer: BPETokenizer,
    config: PretrainConfig,
) -> bool:
    metadata = tokenizer.metadata
    if not metadata:
        return False

    expected_vocab_size = int(config.vocab_size)
    actual_vocab_size = int(metadata.get("vocab_size", -1))
    expected_min_freq = int(config.tokenizer_min_freq)
    actual_min_freq = int(metadata.get("min_freq", -1))
    expected_max_token_length = int(config.tokenizer_max_token_length)
    actual_max_token_length = int(metadata.get("max_token_length", -1))

    return (
        actual_vocab_size == expected_vocab_size
        and actual_min_freq == expected_min_freq
        and actual_max_token_length == expected_max_token_length
    )


def _describe_tokenizer(tokenizer: BPETokenizer) -> str:
    metadata = tokenizer.metadata or {}
    vocab_size = metadata.get("vocab_size", len(tokenizer.vocab))
    min_freq = metadata.get("min_freq", "unknown")
    max_token_length = metadata.get("max_token_length", "unknown")
    return (
        f"vocab_size={vocab_size} "
        f"actual_vocab={len(tokenizer.vocab)} "
        f"min_freq={min_freq} "
        f"max_token_length={max_token_length}"
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
        tokenizer = BPETokenizer.load(config.tokenizer_path)
        if (
            not config.tokenizer_retrain_if_mismatch
            or _tokenizer_matches_config(tokenizer, config)
        ):
            print(
                f"Using existing tokenizer ({_describe_tokenizer(tokenizer)})",
                flush=True,
            )
            return tokenizer

        print(
            "Tokenizer metadata is missing or mismatched; retraining tokenizer "
            f"with vocab_size={config.vocab_size} min_freq={config.tokenizer_min_freq} "
            f"max_token_length={config.tokenizer_max_token_length}",
            flush=True,
        )

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
        f"min_freq={config.tokenizer_min_freq} "
        f"max_token_length={config.tokenizer_max_token_length} "
        f"on {len(tokenizer_sequences)} sampled sequences ({tokenizer_bases:,} bases)",
        flush=True,
    )
    tokenizer = BPETokenizer.train(
        tokenizer_sequences,
        vocab_size=config.vocab_size,
        min_freq=config.tokenizer_min_freq,
        max_token_length=config.tokenizer_max_token_length,
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
