"""Dataset for masked language modeling on genome sequences."""
from __future__ import annotations

import random
from bisect import bisect_right
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np
import torch
from torch.utils.data import Dataset, get_worker_info

from gene_whisperer.tokenization.bpe import BPETokenizer



def _select_window(sequence: str, window_size: int, rng: random.Random) -> str:
    if len(sequence) <= window_size:
        return sequence
    start = rng.randint(0, len(sequence) - window_size)
    return sequence[start : start + window_size]


def _mask_tokens(
    token_ids: List[int],
    tokenizer: BPETokenizer,
    mask_prob: float,
    rng: random.Random,
    replacement_token_ids: Sequence[int],
    maskable_token_ids: Set[int],
    min_masked_tokens: int,
) -> List[int]:
    labels = np.full(len(token_ids), -100, dtype=np.int64)
    candidate_indices = [
        idx
        for idx, token_id in enumerate(token_ids)
        if token_id in maskable_token_ids
    ]
    if not candidate_indices:
        return labels.tolist()

    selected_indices = [idx for idx in candidate_indices if rng.random() < mask_prob]
    if len(candidate_indices) > 1:
        max_masked_count = len(candidate_indices) - 1
    else:
        max_masked_count = 1
    target_masked_count = max(min_masked_tokens, len(selected_indices))
    target_masked_count = min(target_masked_count, max_masked_count)

    if len(selected_indices) < target_masked_count:
        selected_set = set(selected_indices)
        unselected = [idx for idx in candidate_indices if idx not in selected_set]
        rng.shuffle(unselected)
        selected_indices.extend(unselected[: target_masked_count - len(selected_indices)])
    elif len(selected_indices) > target_masked_count:
        rng.shuffle(selected_indices)
        selected_indices = selected_indices[:target_masked_count]

    for idx in selected_indices:
        labels[idx] = token_ids[idx]
        roll = rng.random()
        if roll < 0.8:
            token_ids[idx] = tokenizer.mask_token_id
        elif roll < 0.9:
            token_ids[idx] = rng.choice(replacement_token_ids)
        else:
            continue
    return labels.tolist()


def _build_maskable_token_ids(
    tokenizer: BPETokenizer,
    mask_ambiguous_tokens: bool,
) -> Set[int]:
    maskable_ids: Set[int] = set()
    for token, token_id in tokenizer.vocab.items():
        if token in tokenizer.reserved_tokens:
            continue
        if not mask_ambiguous_tokens and "N" in token:
            continue
        maskable_ids.add(token_id)
    return maskable_ids


def _build_replacement_token_ids(
    tokenizer: BPETokenizer,
    maskable_token_ids: Set[int],
) -> List[int]:
    replacement_ids = [
        token_id
        for token, token_id in tokenizer.vocab.items()
        if token_id in maskable_token_ids
    ]
    if replacement_ids:
        return replacement_ids
    return [tokenizer.unk_token_id]


def _count_maskable_tokens(token_ids: Sequence[int], maskable_token_ids: Set[int]) -> int:
    return sum(1 for token_id in token_ids if token_id in maskable_token_ids)


class MLMDataset(Dataset[Dict[str, torch.Tensor]]):
    """Random window sampling over genome sequences for MLM."""

    def __init__(
        self,
        sequences: Sequence[str],
        tokenizer: BPETokenizer,
        window_size: int = 234,
        max_length: int = 256,
        mask_prob: float = 0.15,
        num_samples: int | None = None,
        seed: int = 1337,
        sample_by_length: bool = True,
        mask_ambiguous_tokens: bool = False,
        min_masked_tokens: int = 1,
        min_maskable_tokens: int | None = None,
        resample_attempts: int = 8,
    ) -> None:
        cleaned_sequences = [seq for seq in sequences if seq]
        if not cleaned_sequences:
            raise ValueError("MLMDataset requires at least one sequence")
        self.sequences = cleaned_sequences
        self.sequence_lengths = [len(seq) for seq in self.sequences]
        self.cumulative_lengths = np.cumsum(self.sequence_lengths).tolist()
        self.total_bases = self.cumulative_lengths[-1]
        self.tokenizer = tokenizer
        self.window_size = window_size
        self.max_length = max_length
        self.mask_prob = mask_prob
        self.num_samples = num_samples or len(self.sequences)
        self.seed = seed
        self.sample_by_length = sample_by_length
        self.min_masked_tokens = max(1, min_masked_tokens)
        self.min_maskable_tokens = max(1, min_maskable_tokens or self.min_masked_tokens)
        self.resample_attempts = max(1, resample_attempts)
        self.maskable_token_ids = _build_maskable_token_ids(
            tokenizer=tokenizer,
            mask_ambiguous_tokens=mask_ambiguous_tokens,
        )
        self.replacement_token_ids = _build_replacement_token_ids(
            tokenizer=tokenizer,
            maskable_token_ids=self.maskable_token_ids,
        )
        self._rng: random.Random | None = None
        self._rng_worker_id: int | None = None

    def __len__(self) -> int:
        return self.num_samples

    def _get_rng(self) -> random.Random:
        worker_info = get_worker_info()
        worker_id = worker_info.id if worker_info is not None else 0
        if self._rng is None or self._rng_worker_id != worker_id:
            self._rng_worker_id = worker_id
            self._rng = random.Random(self.seed + (worker_id * 10_007))
        return self._rng

    def _sample_sequence(self, rng: random.Random, idx: int) -> str:
        if not self.sample_by_length:
            return self.sequences[idx % len(self.sequences)]
        target = rng.randrange(self.total_bases)
        seq_idx = bisect_right(self.cumulative_lengths, target)
        return self.sequences[min(seq_idx, len(self.sequences) - 1)]

    def _sample_encoded_window(
        self,
        idx: int,
        rng: random.Random,
    ) -> Tuple[List[int], List[int]]:
        best_token_ids: List[int] | None = None
        best_attention_mask: List[int] | None = None
        best_maskable_count = -1

        for attempt in range(self.resample_attempts):
            seq = self._sample_sequence(rng=rng, idx=idx + attempt)
            window = _select_window(seq, self.window_size, rng)
            token_ids, attention_mask = self.tokenizer.encode(
                window,
                add_special_tokens=True,
                max_length=self.max_length,
                pad_to_max=True,
            )
            maskable_count = _count_maskable_tokens(token_ids, self.maskable_token_ids)
            if maskable_count >= self.min_maskable_tokens:
                return token_ids, attention_mask
            if maskable_count > best_maskable_count:
                best_maskable_count = maskable_count
                best_token_ids = token_ids
                best_attention_mask = attention_mask

        if best_token_ids is None or best_attention_mask is None:
            raise RuntimeError("Failed to sample MLM window with tokenization output.")
        return best_token_ids, best_attention_mask

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        rng = self._get_rng()
        token_ids, attention_mask = self._sample_encoded_window(idx=idx, rng=rng)
        labels = _mask_tokens(
            token_ids=token_ids,
            tokenizer=self.tokenizer,
            mask_prob=self.mask_prob,
            rng=rng,
            replacement_token_ids=self.replacement_token_ids,
            maskable_token_ids=self.maskable_token_ids,
            min_masked_tokens=self.min_masked_tokens,
        )
        return {
            "input_ids": torch.tensor(token_ids, dtype=torch.long),
            "attention_mask": torch.tensor(attention_mask, dtype=torch.long),
            "labels": torch.tensor(labels, dtype=torch.long),
        }
