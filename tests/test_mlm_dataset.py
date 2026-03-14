from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.datasets.mlm_dataset import (
    MLMDataset,
    _count_maskable_tokens,
    _build_maskable_token_ids,
    _build_replacement_token_ids,
    _mask_tokens,
)
from gene_whisperer.tokenization.bpe import BPETokenizer


def _build_tokenizer() -> BPETokenizer:
    return BPETokenizer.train(
        sequences=["ACGTACGTACGT", "TTTTCCCCAAAAGGGG"],
        vocab_size=64,
    )


def _build_near_character_tokenizer() -> BPETokenizer:
    return BPETokenizer.train(
        sequences=["ACGTACGT", "TGCATGCA"],
        vocab_size=8,
        min_freq=100,
    )


def test_mask_tokens_forces_at_least_one_prediction() -> None:
    tokenizer = _build_tokenizer()
    token_ids, _ = tokenizer.encode("ACGTACGT", add_special_tokens=True, max_length=32, pad_to_max=False)
    original = list(token_ids)
    rng = random.Random(42)

    labels = _mask_tokens(
        token_ids=token_ids,
        tokenizer=tokenizer,
        mask_prob=0.0,
        rng=rng,
        replacement_token_ids=_build_replacement_token_ids(
            tokenizer,
            _build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=True),
        ),
        maskable_token_ids=_build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=True),
        min_masked_tokens=1,
    )

    assert any(label != -100 for label in labels)
    assert token_ids != original
    assert tokenizer.mask_token_id in token_ids


def test_mask_tokens_keeps_at_least_one_context_token() -> None:
    tokenizer = _build_tokenizer()
    token_ids, _ = tokenizer.encode("ACGT", add_special_tokens=True, max_length=16, pad_to_max=False)
    original_token_ids = list(token_ids)
    maskable_ids = _build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=True)
    labels = _mask_tokens(
        token_ids=token_ids,
        tokenizer=tokenizer,
        mask_prob=1.0,
        rng=random.Random(3),
        replacement_token_ids=_build_replacement_token_ids(tokenizer, maskable_ids),
        maskable_token_ids=maskable_ids,
        min_masked_tokens=3,
    )
    masked_count = int(np.sum(np.array(labels) != -100))
    candidate_count = _count_maskable_tokens(original_token_ids, maskable_ids)
    assert masked_count <= max(1, candidate_count - 1)


def test_replacement_ids_exclude_reserved_tokens() -> None:
    tokenizer = _build_tokenizer()
    maskable_ids = _build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=True)
    replacement_ids = _build_replacement_token_ids(tokenizer, maskable_ids)
    reserved_ids = {tokenizer.vocab[token] for token in tokenizer.reserved_tokens}
    assert replacement_ids
    assert not any(token_id in reserved_ids for token_id in replacement_ids)


def test_maskable_ids_exclude_ambiguous_tokens_when_disabled() -> None:
    tokenizer = _build_tokenizer()
    maskable_ids = _build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=False)
    assert tokenizer.vocab["N"] not in maskable_ids


def test_mask_tokens_skips_ambiguous_token_targets() -> None:
    tokenizer = _build_tokenizer()
    token_ids, _ = tokenizer.encode("NNNN", add_special_tokens=True, max_length=16, pad_to_max=False)
    labels = _mask_tokens(
        token_ids=token_ids,
        tokenizer=tokenizer,
        mask_prob=1.0,
        rng=random.Random(1),
        replacement_token_ids=[tokenizer.vocab["A"]],
        maskable_token_ids=_build_maskable_token_ids(tokenizer, mask_ambiguous_tokens=False),
        min_masked_tokens=1,
    )
    assert all(label == -100 for label in labels)


def test_length_weighted_sampling_prefers_longer_sequences() -> None:
    tokenizer = _build_tokenizer()
    short = "A" * 32
    long = "C" * 1024
    dataset = MLMDataset(
        sequences=[short, long],
        tokenizer=tokenizer,
        window_size=32,
        max_length=16,
        num_samples=32,
        seed=7,
        sample_by_length=True,
    )
    rng = random.Random(7)
    picks = [dataset._sample_sequence(rng=rng, idx=i) for i in range(500)]
    short_count = np.sum([len(seq) == len(short) for seq in picks])
    long_count = np.sum([len(seq) == len(long) for seq in picks])

    assert long_count > (short_count * 10)


def test_min_masked_tokens_enforces_multiple_targets() -> None:
    tokenizer = _build_tokenizer()
    dataset = MLMDataset(
        sequences=["ACGT" * 32],
        tokenizer=tokenizer,
        window_size=64,
        max_length=64,
        mask_prob=0.0,
        num_samples=1,
        seed=11,
        min_masked_tokens=2,
    )
    sample = dataset[0]
    labels = sample["labels"].numpy()
    assert int(np.sum(labels != -100)) >= 2
    assert int(sample["masked_count"].item()) >= 2
    assert int(sample["maskable_count"].item()) >= 2
    assert int(sample["tokenized_count"].item()) >= 2


def test_dataset_resamples_when_window_has_too_few_maskable_tokens() -> None:
    tokenizer = _build_tokenizer()
    dataset = MLMDataset(
        sequences=["N" * 512, "ACGT" * 128],
        tokenizer=tokenizer,
        window_size=128,
        max_length=64,
        mask_prob=0.15,
        num_samples=1,
        seed=5,
        sample_by_length=False,
        mask_ambiguous_tokens=False,
        min_masked_tokens=2,
        min_maskable_tokens=4,
        resample_attempts=4,
    )
    sample = dataset[0]
    labels = sample["labels"].numpy()
    assert int(np.sum(labels != -100)) >= 2


def test_dataset_resamples_when_tokenized_window_is_too_short() -> None:
    tokenizer = _build_near_character_tokenizer()
    dataset = MLMDataset(
        sequences=["A" * 16, "ACGT" * 64],
        tokenizer=tokenizer,
        window_size=128,
        max_length=64,
        mask_prob=0.15,
        num_samples=1,
        seed=9,
        sample_by_length=False,
        min_masked_tokens=2,
        min_maskable_tokens=6,
        min_tokenized_tokens=32,
        resample_attempts=4,
    )
    sample = dataset[0]
    attention_mask = sample["attention_mask"].numpy()
    assert int(np.sum(attention_mask)) >= 32
