"""Minimal BPE tokenizer inspired by DNABERT2."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

DNA_ALPHABET = set("ACGTN")
DEFAULT_RESERVED_TOKENS = ["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"]


@dataclass(frozen=True)
class EncodedBatch:
    input_ids: List[List[int]]
    attention_mask: List[List[int]]


def normalize_sequence(sequence: str) -> str:
    """Uppercase and strip to canonical A/C/G/T/N characters."""
    if not isinstance(sequence, str):
        return ""
    return "".join(base for base in sequence.upper().strip() if base in DNA_ALPHABET)


def _merge_pair(tokens: List[str], pair: Tuple[str, str]) -> List[str]:
    merged: List[str] = []
    i = 0
    while i < len(tokens):
        if i < len(tokens) - 1 and tokens[i] == pair[0] and tokens[i + 1] == pair[1]:
            merged.append(tokens[i] + tokens[i + 1])
            i += 2
        else:
            merged.append(tokens[i])
            i += 1
    return merged


def _get_pair_counts(seqs: Sequence[List[str]]) -> Dict[Tuple[str, str], int]:
    counts: Dict[Tuple[str, str], int] = {}
    for tokens in seqs:
        for left, right in zip(tokens, tokens[1:]):
            pair = (left, right)
            counts[pair] = counts.get(pair, 0) + 1
    return counts


class BPETokenizer:
    """Simple BPE tokenizer for DNA sequences."""

    def __init__(
        self,
        vocab: Dict[str, int],
        merges: Sequence[Tuple[str, str]],
        reserved_tokens: Sequence[str] = DEFAULT_RESERVED_TOKENS,
    ) -> None:
        self.vocab = dict(vocab)
        self.id_to_token = {idx: tok for tok, idx in self.vocab.items()}
        self.merges = list(merges)
        self.bpe_ranks = {pair: idx for idx, pair in enumerate(self.merges)}
        self.reserved_tokens = list(reserved_tokens)
        self.pad_token_id = self.vocab[self.reserved_tokens[0]]
        self.unk_token_id = self.vocab[self.reserved_tokens[1]]
        self.cls_token_id = self.vocab[self.reserved_tokens[2]]
        self.sep_token_id = self.vocab[self.reserved_tokens[3]]
        self.mask_token_id = self.vocab[self.reserved_tokens[4]]
        self._cache: Dict[str, List[str]] = {}

    @classmethod
    def train(
        cls,
        sequences: Iterable[str],
        vocab_size: int = 4096,
        min_freq: int = 2,
        reserved_tokens: Sequence[str] = DEFAULT_RESERVED_TOKENS,
    ) -> "BPETokenizer":
        if vocab_size <= len(reserved_tokens):
            raise ValueError("vocab_size must exceed number of reserved tokens")

        normalized = [normalize_sequence(seq) for seq in sequences]
        token_seqs = [list(seq) for seq in normalized if seq]
        if not token_seqs:
            raise ValueError("No sequences provided for BPE training")

        vocab: Dict[str, int] = {tok: idx for idx, tok in enumerate(reserved_tokens)}
        for base in sorted(DNA_ALPHABET):
            if base not in vocab:
                vocab[base] = len(vocab)

        merges: List[Tuple[str, str]] = []
        while len(vocab) < vocab_size:
            pair_counts = _get_pair_counts(token_seqs)
            if not pair_counts:
                break
            best_pair, best_count = max(pair_counts.items(), key=lambda item: item[1])
            if best_count < min_freq:
                break

            merged_token = "".join(best_pair)
            if merged_token in vocab:
                break

            token_seqs = [_merge_pair(tokens, best_pair) for tokens in token_seqs]
            vocab[merged_token] = len(vocab)
            merges.append(best_pair)

        return cls(vocab=vocab, merges=merges, reserved_tokens=reserved_tokens)

    def _apply_bpe(self, tokens: List[str]) -> List[str]:
        if not tokens:
            return []
        cache_key = "".join(tokens)
        if cache_key in self._cache:
            return list(self._cache[cache_key])

        while True:
            pairs = set(zip(tokens, tokens[1:]))
            if not pairs:
                break
            best_pair = min(pairs, key=lambda pair: self.bpe_ranks.get(pair, 1e12))
            if best_pair not in self.bpe_ranks:
                break
            tokens = _merge_pair(tokens, best_pair)
        self._cache[cache_key] = list(tokens)
        return tokens

    def encode(
        self,
        sequence: str,
        add_special_tokens: bool = True,
        max_length: int | None = None,
        pad_to_max: bool = False,
    ) -> Tuple[List[int], List[int]]:
        normalized = normalize_sequence(sequence)
        tokens = self._apply_bpe(list(normalized))
        token_ids = [self.vocab.get(tok, self.unk_token_id) for tok in tokens]

        if add_special_tokens:
            token_ids = [self.cls_token_id] + token_ids + [self.sep_token_id]

        if max_length is not None:
            token_ids = token_ids[:max_length]

        attention_mask = [1] * len(token_ids)
        if pad_to_max and max_length is not None and len(token_ids) < max_length:
            pad_len = max_length - len(token_ids)
            token_ids = token_ids + [self.pad_token_id] * pad_len
            attention_mask = attention_mask + [0] * pad_len

        return token_ids, attention_mask

    def encode_batch(
        self,
        sequences: Sequence[str],
        add_special_tokens: bool = True,
        max_length: int | None = None,
        pad_to_max: bool = True,
    ) -> EncodedBatch:
        input_ids: List[List[int]] = []
        attention_masks: List[List[int]] = []
        for seq in sequences:
            ids, mask = self.encode(
                seq,
                add_special_tokens=add_special_tokens,
                max_length=max_length,
                pad_to_max=pad_to_max,
            )
            input_ids.append(ids)
            attention_masks.append(mask)
        return EncodedBatch(input_ids=input_ids, attention_mask=attention_masks)

    def decode(self, token_ids: Sequence[int], skip_special_tokens: bool = True) -> str:
        tokens = [self.id_to_token.get(idx, "") for idx in token_ids]
        if skip_special_tokens:
            tokens = [tok for tok in tokens if tok not in self.reserved_tokens]
        return "".join(tokens)

    def save(self, path: Path) -> None:
        data = {
            "vocab": self.vocab,
            "merges": [list(pair) for pair in self.merges],
            "reserved_tokens": self.reserved_tokens,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, sort_keys=True)

    @classmethod
    def load(cls, path: Path) -> "BPETokenizer":
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        merges = [tuple(pair) for pair in data["merges"]]
        return cls(
            vocab={str(k): int(v) for k, v in data["vocab"].items()},
            merges=merges,
            reserved_tokens=data.get("reserved_tokens", DEFAULT_RESERVED_TOKENS),
        )
