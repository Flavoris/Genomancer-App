"""Dataset for masked language modeling on genome sequences."""
from __future__ import annotations

import random
from typing import Dict, List, Sequence

import numpy as np
import torch
from torch.utils.data import Dataset

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
) -> List[int]:
    labels = np.full(len(token_ids), -100, dtype=np.int64)
    for idx, token_id in enumerate(token_ids):
        if token_id in (tokenizer.pad_token_id, tokenizer.cls_token_id, tokenizer.sep_token_id):
            continue
        if rng.random() < mask_prob:
            labels[idx] = token_id
            roll = rng.random()
            if roll < 0.8:
                token_ids[idx] = tokenizer.mask_token_id
            elif roll < 0.9:
                token_ids[idx] = rng.randint(0, len(tokenizer.vocab) - 1)
            else:
                token_ids[idx] = token_id
    return labels.tolist()


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
    ) -> None:
        if not sequences:
            raise ValueError("MLMDataset requires at least one sequence")
        self.sequences = list(sequences)
        self.tokenizer = tokenizer
        self.window_size = window_size
        self.max_length = max_length
        self.mask_prob = mask_prob
        self.num_samples = num_samples or len(self.sequences)
        self.rng = random.Random(seed)

    def __len__(self) -> int:
        return self.num_samples

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        seq = self.sequences[idx % len(self.sequences)]
        window = _select_window(seq, self.window_size, self.rng)
        token_ids, attention_mask = self.tokenizer.encode(
            window, add_special_tokens=True, max_length=self.max_length, pad_to_max=True
        )
        labels = _mask_tokens(token_ids, self.tokenizer, self.mask_prob, self.rng)
        return {
            "input_ids": torch.tensor(token_ids, dtype=torch.long),
            "attention_mask": torch.tensor(attention_mask, dtype=torch.long),
            "labels": torch.tensor(labels, dtype=torch.long),
        }
