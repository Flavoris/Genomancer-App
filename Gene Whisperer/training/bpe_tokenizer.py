"""
DNA-specific Byte Pair Encoding (BPE) tokenizer following DNABERT-2.

BPE replaces k-mer tokenization by iteratively merging the most frequent
co-occurring genome segments, producing variable-length tokens that:
- Prevent information leakage (unlike overlapping k-mer)
- Maintain robust representations (unlike non-overlapping k-mer)
- Reduce sequence length by ~5x (improving computational efficiency)
- Capture biologically meaningful motifs (TATA, CAAT, etc.)

Reference: Zhou et al., "DNABERT-2: Efficient Foundation Model and
Benchmark for Multi-Species Genomes", ICLR 2024.
"""

import json
import random
import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import torch


# Special tokens matching common genome language model conventions
SPECIAL_TOKENS = {
    "[PAD]": 0,
    "[UNK]": 1,
    "[CLS]": 2,
    "[SEP]": 3,
    "[MASK]": 4,
}

# DNA nucleotide bases
DNA_BASES = ["A", "T", "C", "G"]


class DNABPETokenizer:
    """BPE tokenizer optimized for DNA sequences following DNABERT-2.

    Learns a fixed-size vocabulary of variable-length tokens from DNA
    sequences by iteratively merging the most frequent adjacent pairs.
    Target vocabulary size of 4096 balances model performance with
    computational efficiency (per DNABERT-2 empirical analysis).
    """

    def __init__(self, vocab_size: int = 4096):
        self.vocab_size = vocab_size
        self.vocab: Dict[str, int] = {}
        self.id_to_token: Dict[int, str] = {}
        self.merges: List[Tuple[str, str]] = []

    @property
    def pad_id(self) -> int:
        return SPECIAL_TOKENS["[PAD]"]

    @property
    def unk_id(self) -> int:
        return SPECIAL_TOKENS["[UNK]"]

    @property
    def mask_id(self) -> int:
        return SPECIAL_TOKENS["[MASK]"]

    def train(
        self,
        sequences: List[str],
        min_frequency: int = 2,
        verbose: bool = True,
    ) -> None:
        """Train BPE vocabulary from DNA sequences.

        Algorithm (following DNABERT-2 / SentencePiece BPE):
        1. Initialize vocab with A, T, C, G + special tokens
        2. Split all sequences into individual characters
        3. Count all adjacent token pairs across corpus
        4. Merge the most frequent pair into a new token
        5. Repeat until vocab_size is reached

        Uses optimized incremental pair counting to avoid O(n) rescans.

        Args:
            sequences: List of DNA sequence strings (A/T/C/G only).
            min_frequency: Minimum pair frequency to consider for merging.
            verbose: Print progress every 500 merges.
        """
        self._initialize_vocab()

        # Represent each sequence as a list of character tokens
        corpus = self._prepare_corpus(sequences)

        # Initial pair count (only done once)
        pair_counts = self._count_pairs(corpus)

        num_merges = self.vocab_size - len(self.vocab)
        for merge_idx in range(num_merges):
            if not pair_counts:
                break

            best_pair = max(pair_counts, key=pair_counts.get)
            if pair_counts[best_pair] < min_frequency:
                break

            merged_token = best_pair[0] + best_pair[1]
            self.merges.append(best_pair)
            token_id = len(self.vocab)
            self.vocab[merged_token] = token_id
            self.id_to_token[token_id] = merged_token

            # Apply merge and update counts incrementally
            corpus, pair_counts = self._apply_merge_incremental(
                corpus, best_pair, merged_token, pair_counts
            )

            if verbose and (merge_idx + 1) % 500 == 0:
                print(
                    f"  Merge {merge_idx + 1}/{num_merges}: "
                    f"'{best_pair[0]}' + '{best_pair[1]}' -> '{merged_token}'"
                )

    def tokenize(self, sequence: str) -> List[int]:
        """Tokenize a DNA sequence using learned BPE vocabulary.

        Uses greedy longest-match strategy: applies merges in the order
        they were learned during training (priority-based).

        Args:
            sequence: DNA sequence string (A/T/C/G characters).

        Returns:
            List of token IDs.
        """
        sequence = self._sanitize(sequence)
        tokens = list(sequence)
        # Apply merges in learned order (most frequent first)
        for left, right in self.merges:
            tokens = self._apply_merge_to_tokens(tokens, left, right)

        return [self.vocab.get(t, self.unk_id) for t in tokens]

    def decode(self, token_ids: List[int]) -> str:
        """Decode token IDs back to a DNA sequence string.

        Skips special tokens (PAD, CLS, SEP, UNK, MASK) and unknown IDs,
        returning only the reconstructed DNA bases.

        Args:
            token_ids: List of integer token IDs.

        Returns:
            Reconstructed DNA sequence string.
        """
        special_ids = set(SPECIAL_TOKENS.values())
        tokens = []
        for tid in token_ids:
            if tid in special_ids:
                continue
            if tid in self.id_to_token:
                tokens.append(self.id_to_token[tid])
        return "".join(tokens)

    def encode_with_special(self, sequence: str) -> List[int]:
        """Tokenize with [CLS] prepended and [SEP] appended.

        Args:
            sequence: DNA sequence string.

        Returns:
            Token IDs with special tokens wrapping the sequence.
        """
        ids = self.tokenize(sequence)
        return [SPECIAL_TOKENS["[CLS]"]] + ids + [SPECIAL_TOKENS["[SEP]"]]

    def tokenize_and_pad(
        self, sequence: str, max_token_len: int
    ) -> torch.LongTensor:
        """Tokenize and pad/truncate to fixed length for batching.

        Args:
            sequence: DNA sequence string.
            max_token_len: Target output length (pad or truncate to this).

        Returns:
            LongTensor of shape (max_token_len,) with token IDs.
        """
        ids = self.tokenize(sequence)
        # Truncate if too long
        if len(ids) > max_token_len:
            ids = ids[:max_token_len]
        # Pad if too short
        while len(ids) < max_token_len:
            ids.append(self.pad_id)
        return torch.LongTensor(ids)

    def random_token_id(self) -> int:
        """Return a random non-special token ID for MLM random replacement."""
        base_start = len(SPECIAL_TOKENS)
        return random.randint(base_start, len(self.vocab) - 1)

    @property
    def itos(self) -> List[str]:
        """Index-to-string list for compatibility with pipeline code."""
        result = [""] * len(self.vocab)
        for token, idx in self.vocab.items():
            result[idx] = token
        return result

    def save(self, path: str) -> None:
        """Save vocabulary and merges to a JSON file.

        Args:
            path: File path to save the tokenizer state.
        """
        state = {
            "vocab_size": self.vocab_size,
            "vocab": self.vocab,
            "merges": self.merges,
        }
        with open(path, "w") as f:
            json.dump(state, f, indent=2)

    @classmethod
    def load(cls, path: str) -> "DNABPETokenizer":
        """Load a pre-trained tokenizer from a JSON file.

        Args:
            path: File path to load the tokenizer state from.

        Returns:
            Initialized DNABPETokenizer with loaded vocabulary.
        """
        with open(path, "r") as f:
            state = json.load(f)

        tokenizer = cls(vocab_size=state["vocab_size"])
        tokenizer.vocab = state["vocab"]
        tokenizer.merges = [tuple(m) for m in state["merges"]]
        tokenizer.id_to_token = {v: k for k, v in tokenizer.vocab.items()}
        return tokenizer

    def _initialize_vocab(self) -> None:
        """Initialize vocabulary with special tokens and single nucleotides."""
        self.vocab = dict(SPECIAL_TOKENS)
        for i, base in enumerate(DNA_BASES):
            idx = len(SPECIAL_TOKENS) + i
            self.vocab[base] = idx
        self.id_to_token = {v: k for k, v in self.vocab.items()}

    def _sanitize(self, sequence: str) -> str:
        """Normalize a DNA sequence: uppercase, keep only A/T/C/G."""
        sequence = sequence.upper()
        return re.sub(r"[^ATCG]", "", sequence)

    def _prepare_corpus(self, sequences: List[str]) -> List[List[str]]:
        """Sanitize and split sequences into character-level token lists."""
        corpus = []
        for seq in sequences:
            cleaned = self._sanitize(seq)
            if cleaned:
                corpus.append(list(cleaned))
        return corpus

    def _count_pairs(self, corpus: List[List[str]]) -> Counter:
        """Count frequencies of all adjacent token pairs in corpus."""
        pair_counts: Counter = Counter()
        for tokens in corpus:
            for i in range(len(tokens) - 1):
                pair_counts[(tokens[i], tokens[i + 1])] += 1
        return pair_counts

    def _apply_merge(
        self,
        corpus: List[List[str]],
        pair: Tuple[str, str],
        merged: str,
    ) -> List[List[str]]:
        """Replace all occurrences of a pair in corpus with merged token."""
        new_corpus = []
        for tokens in corpus:
            new_corpus.append(
                self._apply_merge_to_tokens(tokens, pair[0], pair[1])
            )
        return new_corpus

    def _apply_merge_incremental(
        self,
        corpus: List[List[str]],
        pair: Tuple[str, str],
        merged: str,
        pair_counts: Counter,
    ) -> Tuple[List[List[str]], Counter]:
        """Apply merge and update pair counts incrementally.

        This is O(occurrences) instead of O(corpus_size), making it much faster.

        When we merge (A, B) -> AB at position i:
        1. Remove (A, B) count
        2. If there was X before A: remove (X, A), add (X, AB)
        3. If there was Y after B: remove (B, Y), add (AB, Y)
        """
        left, right = pair
        new_corpus = []

        # Remove the merged pair count entirely
        if pair in pair_counts:
            del pair_counts[pair]

        for tokens in corpus:
            if len(tokens) < 2:
                new_corpus.append(tokens)
                continue

            # First pass: find all merge positions (non-overlapping, left to right)
            merge_positions = set()
            i = 0
            while i < len(tokens) - 1:
                if tokens[i] == left and tokens[i + 1] == right:
                    merge_positions.add(i)
                    i += 2  # Skip to avoid overlapping merges
                else:
                    i += 1

            if not merge_positions:
                new_corpus.append(tokens)
                continue

            # Second pass: build new tokens and update pair counts
            new_tokens = []
            i = 0
            while i < len(tokens):
                if i in merge_positions:
                    # Merging tokens[i] and tokens[i+1]
                    # Update left neighbor pair
                    if new_tokens:
                        prev_token = new_tokens[-1]
                        # Remove old pair (prev, left)
                        old_pair = (prev_token, left)
                        if old_pair in pair_counts:
                            pair_counts[old_pair] -= 1
                            if pair_counts[old_pair] <= 0:
                                del pair_counts[old_pair]
                        # Add new pair (prev, merged)
                        pair_counts[(prev_token, merged)] += 1

                    # Update right neighbor pair
                    if i + 2 < len(tokens) and (i + 2) not in merge_positions:
                        next_token = tokens[i + 2]
                        # Remove old pair (right, next)
                        old_pair = (right, next_token)
                        if old_pair in pair_counts:
                            pair_counts[old_pair] -= 1
                            if pair_counts[old_pair] <= 0:
                                del pair_counts[old_pair]
                        # Add new pair (merged, next)
                        pair_counts[(merged, next_token)] += 1

                    new_tokens.append(merged)
                    i += 2
                else:
                    new_tokens.append(tokens[i])
                    i += 1

            new_corpus.append(new_tokens)

        return new_corpus, pair_counts

    def _apply_merge_to_tokens(
        self, tokens: List[str], left: str, right: str
    ) -> List[str]:
        """Merge adjacent left+right tokens into a single token."""
        merged = left + right
        result = []
        i = 0
        while i < len(tokens):
            if (
                i < len(tokens) - 1
                and tokens[i] == left
                and tokens[i + 1] == right
            ):
                result.append(merged)
                i += 2
            else:
                result.append(tokens[i])
                i += 1
        return result

    def get_token_lengths(self) -> Dict[int, int]:
        """Return distribution of token lengths in vocabulary.

        Returns:
            Dictionary mapping token length to count of tokens with that length.
        """
        length_dist: Dict[int, int] = {}
        for token in self.vocab:
            if token not in SPECIAL_TOKENS:
                length = len(token)
                length_dist[length] = length_dist.get(length, 0) + 1
        return length_dist

    def compression_ratio(self, sequence: str) -> float:
        """Compute the compression ratio for a given sequence.

        Returns the ratio of original length to tokenized length,
        indicating how much shorter the tokenized representation is.

        Args:
            sequence: DNA sequence string.

        Returns:
            Float ratio (higher means more compression).
        """
        cleaned = self._sanitize(sequence)
        if not cleaned:
            return 0.0
        token_ids = self.tokenize(cleaned)
        if not token_ids:
            return 0.0
        return len(cleaned) / len(token_ids)

    def __len__(self) -> int:
        """Return current vocabulary size."""
        return len(self.vocab)
