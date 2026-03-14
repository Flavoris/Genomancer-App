from __future__ import annotations

import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.pretrain_config import PretrainConfig
from gene_whisperer.training.pretrain_runtime import maybe_train_tokenizer


def _build_config(tokenizer_path: Path) -> PretrainConfig:
    return PretrainConfig(
        fasta_paths=[],
        tokenizer_path=tokenizer_path,
        vocab_size=64,
        max_length=64,
        window_size=32,
        mask_prob=0.15,
        max_bases_per_file=None,
        tokenizer_max_bases=None,
        tokenizer_max_sequences=None,
        tokenizer_window_size=128,
        tokenizer_min_freq=8,
        tokenizer_max_token_length=8,
        tokenizer_retrain_if_mismatch=True,
        sample_by_length=True,
        mask_ambiguous_tokens=False,
        min_masked_tokens=2,
        min_maskable_tokens=4,
        min_tokenized_tokens=8,
        resample_attempts=4,
        batch_size=8,
        epochs=2,
        min_epochs=1,
        early_stopping_patience=2,
        early_stopping_min_delta=0.0,
        save_best_only=True,
        samples_per_epoch=16,
        log_interval=10,
        num_workers=0,
        warmup_ratio=0.03,
        min_lr_ratio=0.05,
        max_grad_norm=1.0,
        use_amp=False,
        lr=1e-4,
        weight_decay=0.01,
        grad_accum_steps=1,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        dropout=0.1,
        seed=7,
        output_dir=tokenizer_path.parent / "out",
    )


def test_maybe_train_tokenizer_retrains_on_metadata_mismatch(tmp_path: Path) -> None:
    tokenizer_path = tmp_path / "tokenizer.json"
    old_tokenizer = BPETokenizer.train(
        ["ACGTACGTACGT"] * 4,
        vocab_size=64,
        min_freq=2,
        max_token_length=32,
    )
    old_tokenizer.save(tokenizer_path)

    config = _build_config(tokenizer_path)
    tokenizer = maybe_train_tokenizer(
        config=config,
        sequences=["ACGTACGTACGT"] * 8,
    )

    assert tokenizer.metadata["min_freq"] == config.tokenizer_min_freq
    assert tokenizer.metadata["max_token_length"] == config.tokenizer_max_token_length


def test_maybe_train_tokenizer_reuses_matching_metadata(tmp_path: Path) -> None:
    tokenizer_path = tmp_path / "tokenizer.json"
    config = _build_config(tokenizer_path)
    tokenizer = BPETokenizer.train(
        ["ACGTACGTACGT"] * 8,
        vocab_size=config.vocab_size,
        min_freq=config.tokenizer_min_freq,
        max_token_length=config.tokenizer_max_token_length,
    )
    tokenizer.save(tokenizer_path)

    loaded = maybe_train_tokenizer(
        config=config,
        sequences=["TTTTCCCCAAAA"] * 8,
    )

    assert loaded.metadata == tokenizer.metadata
    assert loaded.vocab == tokenizer.vocab
