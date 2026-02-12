"""Dataset for masked language modeling on genome sequences."""
from __future__ import annotations

from bisect import bisect_right
import random
from typing import Dict, List, Sequence

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
) -> List[int]:
    labels = np.full(len(token_ids), -100, dtype=np.int64)
    candidate_indices: List[int] = []
    for idx, token_id in enumerate(token_ids):
        if token_id in (tokenizer.pad_token_id, tokenizer.cls_token_id, tokenizer.sep_token_id):
            continue
        candidate_indices.append(idx)
        if rng.random() < mask_prob:
            labels[idx] = token_id
            roll = rng.random()
            if roll < 0.8:
                token_ids[idx] = tokenizer.mask_token_id
            elif roll < 0.9:
                token_ids[idx] = rng.choice(replacement_token_ids)
            else:
                token_ids[idx] = token_id

    if candidate_indices and np.all(labels == -100):
        forced_idx = rng.choice(candidate_indices)
        labels[forced_idx] = token_ids[forced_idx]
        token_ids[forced_idx] = tokenizer.mask_token_id
    return labels.tolist()


def _build_replacement_token_ids(tokenizer: BPETokenizer) -> List[int]:
    replacement_ids = [
        token_id
        for token, token_id in tokenizer.vocab.items()
        if token not in tokenizer.reserved_tokens
    ]
    if replacement_ids:
        return replacement_ids
    return [tokenizer.unk_token_id]


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
        self.replacement_token_ids = _build_replacement_token_ids(tokenizer)
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

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        rng = self._get_rng()
        seq = self._sample_sequence(rng=rng, idx=idx)
        window = _select_window(seq, self.window_size, rng)
        token_ids, attention_mask = self.tokenizer.encode(
            window, add_special_tokens=True, max_length=self.max_length, pad_to_max=True
        )
        labels = _mask_tokens(
            token_ids=token_ids,
            tokenizer=self.tokenizer,
            mask_prob=self.mask_prob,
            rng=rng,
            replacement_token_ids=self.replacement_token_ids,
        )
        return {
            "input_ids": torch.tensor(token_ids, dtype=torch.long),
            "attention_mask": torch.tensor(attention_mask, dtype=torch.long),
            "labels": torch.tensor(labels, dtype=torch.long),
        }
