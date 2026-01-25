"""Dataset utilities for Gene Whisperer promoter models."""
from __future__ import annotations

import itertools
import json
import logging
import math
import os
import random
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler

from bpe_tokenizer import DNABPETokenizer

LOGGER = logging.getLogger(__name__)
if not LOGGER.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    LOGGER.addHandler(handler)
LOGGER.setLevel(logging.INFO)

# Note: Random seeding is now handled by seed_utils.set_global_seed()
# The GENE_WHISPERER_SEED environment variable is set by seed_utils when called
# For backwards compatibility, we still check the env var and apply it if present
RANDOM_SEED = int(os.environ.get("GENE_WHISPERER_SEED", 1337))
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)

BASES: Tuple[str, ...] = ("A", "C", "G", "T")
BASES_SET: set = {"A", "C", "G", "T"}  # For fast membership checks
TRI_KMER_VOCAB: List[str] = ["".join(kmer) for kmer in itertools.product(BASES, repeat=3)]
TRI_KMER_TO_IDX: Dict[str, int] = {kmer: idx for idx, kmer in enumerate(TRI_KMER_VOCAB)}
DI_KMER_VOCAB: List[str] = ["".join(kmer) for kmer in itertools.product(BASES, repeat=2)]
DI_KMER_TO_IDX: Dict[str, int] = {kmer: idx for idx, kmer in enumerate(DI_KMER_VOCAB)}
BASE_TO_IDX: Dict[str, int] = {base: idx for idx, base in enumerate(BASES)}
PSTNP_CORE_WEIGHT = 2.0
EIIP_VALUES: Dict[str, float] = {
    "A": 0.1260,
    "C": 0.1340,
    "G": 0.0806,
    "T": 0.1335,
}
# Complement map: A<->T, C<->G
REV_COMP_MAP = str.maketrans("ACGT", "TGCA")


class KmerVocabulary:
    """Vocabulary built from observed k-mers with common special tokens."""

    def __init__(
        self,
        k: int,
        tokens: Sequence[str],
        unk_token: str = "<UNK>",
        pad_token: str = "<PAD>",
        mask_token: str = "[MASK]",
    ):
        self.k = k
        self.unk_token = unk_token
        self.pad_token = pad_token
        self.mask_token = mask_token
        unique_tokens: List[str] = []
        seen = set()
        for token in tokens:
            if token not in seen and token not in (unk_token, pad_token, mask_token):
                unique_tokens.append(token)
                seen.add(token)
        for special in (mask_token, unk_token, pad_token):
            if special not in seen:
                unique_tokens.append(special)
                seen.add(special)
        self.itos = unique_tokens
        self.stoi = {token: idx for idx, token in enumerate(self.itos)}
        self.mask_id = self.stoi[mask_token]
        self.unk_id = self.stoi[unk_token]
        self.pad_id = self.stoi[pad_token]
        self._base_token_ids = [
            idx for idx, token in enumerate(self.itos) if token not in {mask_token, unk_token, pad_token}
        ]

    def __len__(self) -> int:
        return len(self.itos)

    @classmethod
    def build_from_sequences(
        cls,
        sequences: Sequence[str],
        k: int,
        unk_token: str = "<UNK>",
        pad_token: str = "<PAD>",
        mask_token: str = "[MASK]",
    ):
        observed = set()
        for seq in sequences:
            cleaned = _sanitize_sequence(seq)
            limit = len(cleaned) - k + 1
            for i in range(max(0, limit)):
                observed.add(cleaned[i : i + k])
        if not observed:
            observed.add("A" * k)
        tokens = sorted(observed)
        return cls(k=k, tokens=tokens, unk_token=unk_token, pad_token=pad_token, mask_token=mask_token)

    @classmethod
    def load(cls, path: Path) -> "KmerVocabulary":
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        if "mask_token" not in data:
            data["mask_token"] = "[MASK]"
        vocab = cls(
            k=data["k"],
            tokens=data["itos"],
            unk_token=data.get("unk_token", "<UNK>"),
            pad_token=data.get("pad_token", "<PAD>"),
            mask_token=data.get("mask_token", "[MASK]"),
        )
        return vocab

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "k": self.k,
                    "itos": self.itos,
                    "unk_token": self.unk_token,
                    "pad_token": self.pad_token,
                    "mask_token": self.mask_token,
                },
                handle,
                indent=2,
            )

    def tokenize(self, sequence: str, max_bp_len: int) -> torch.LongTensor:
        cleaned = _sanitize_sequence(sequence)
        if len(cleaned) > max_bp_len:
            cleaned = cleaned[:max_bp_len]
        else:
            cleaned = cleaned.ljust(max_bp_len, "A")
        max_tokens = max(max_bp_len - self.k + 1, 0)
        token_ids: List[int] = []
        for i in range(max_tokens):
            kmer = cleaned[i : i + self.k]
            if len(kmer) < self.k:
                token_ids.append(self.pad_id)
            else:
                token_ids.append(self.stoi.get(kmer, self.unk_id))
        return torch.tensor(token_ids, dtype=torch.long)

    def random_token_id(self) -> int:
        if not self._base_token_ids:
            return 0
        return random.choice(self._base_token_ids)

    def decode(self, token_ids: List[int]) -> str:
        """Decode token IDs back to a DNA sequence.

        For k-mer tokenization, consecutive k-mers overlap by k-1 bases.
        Reconstruction: first k-mer provides k bases, each subsequent k-mer adds 1 base.
        """
        if not token_ids:
            return ""

        # Filter out padding tokens
        filtered_ids = [tid for tid in token_ids if tid != self.pad_id]
        if not filtered_ids:
            return ""

        # Get first k-mer
        first_kmer = self.itos[filtered_ids[0]] if filtered_ids[0] < len(self.itos) else "N" * self.k
        if first_kmer in (self.unk_token, self.mask_token, self.pad_token):
            first_kmer = "N" * self.k

        result = list(first_kmer)

        # Each subsequent k-mer adds one base (the last character)
        for tid in filtered_ids[1:]:
            if tid < len(self.itos):
                kmer = self.itos[tid]
                if kmer in (self.unk_token, self.mask_token, self.pad_token):
                    result.append("N")
                else:
                    result.append(kmer[-1])  # Add last base of k-mer
            else:
                result.append("N")

        return "".join(result)


def _get_or_create_vocab(sequences: Sequence[str], k: int, vocab_dir: Path) -> KmerVocabulary:
    vocab_dir.mkdir(parents=True, exist_ok=True)
    cache_path = vocab_dir / f"k{k}_vocab.json"
    if cache_path.exists():
        vocab = KmerVocabulary.load(cache_path)
        if vocab.k != k:
            raise ValueError(f"Vocabulary at {cache_path} was built for k={vocab.k}, expected {k}")
        return vocab
    vocab = KmerVocabulary.build_from_sequences(sequences, k)
    vocab.save(cache_path)
    LOGGER.info("Built k=%d vocabulary with %d entries", k, len(vocab.itos))
    return vocab


def _cfg_value(cfg: Any, key: str, default: Any = None) -> Any:
    if isinstance(cfg, dict):
        return cfg.get(key, default)
    return getattr(cfg, key, default)


def _sanitize_sequence(seq: str) -> str:
    """
    Normalize DNA sequence for downstream encoding.

    - Uppercase
    - Strip whitespace
    - KEEP only A/C/G/T
    - Drop everything else (spaces, N, etc.)
    """
    seq = (seq or "").strip().upper()
    return "".join(base for base in seq if base in BASES_SET)


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    
    Sanitizes the input first (uppercase, keeps only ACGT),
    then complements (A<->T, C<->G) and reverses.
    """
    sanitized = _sanitize_sequence(seq)
    return sanitized.translate(REV_COMP_MAP)[::-1]


# Alias for backwards compatibility
_reverse_complement = reverse_complement


def _iter_kmers(sequence: str, k: int = 3) -> Iterable[str]:
    cleaned = _sanitize_sequence(sequence)
    limit = len(cleaned) - k + 1
    for i in range(max(0, limit)):
        yield cleaned[i : i + k]


def compute_tnc(sequence: str) -> torch.FloatTensor:
    counts = np.zeros(len(TRI_KMER_VOCAB), dtype=np.float32)
    total = 0
    for kmer in _iter_kmers(sequence, 3):
        counts[TRI_KMER_TO_IDX[kmer]] += 1
        total += 1
    if total > 0:
        counts /= float(total)
    return torch.from_numpy(counts)


def compute_eiip(sequence: str) -> torch.FloatTensor:
    accum = np.zeros(len(TRI_KMER_VOCAB), dtype=np.float32)
    kmers = list(_iter_kmers(sequence, 3))
    total = len(kmers)
    if total == 0:
        return torch.from_numpy(accum)
    for kmer in kmers:
        kmer_idx = TRI_KMER_TO_IDX[kmer]
        value = sum(EIIP_VALUES.get(base, 0.0) for base in kmer) / 3.0
        accum[kmer_idx] += value
    accum /= float(total)
    return torch.from_numpy(accum.astype(np.float32))


def compute_pseeiip(sequence: str) -> torch.FloatTensor:
    """
    Compute PseEIIP (Pseudo Electron-Ion Interaction Potential) features.

    Following PROCABLES paper:
    V[i] = EIIP_trinuc[i] * frequency[i]

    Where EIIP_trinuc is the SUM of three nucleotide EIIP values.
    """
    seq = _sanitize_sequence(sequence)
    L = len(seq)

    if L < 3:
        return torch.zeros(64, dtype=torch.float32)

    # Count trinucleotides.
    counts = np.zeros(64, dtype=np.float32)
    for i in range(L - 2):
        trinuc = seq[i : i + 3]
        if all(b in BASES_SET for b in trinuc):
            # Trinucleotide index (AAA=0, AAC=1, ..., TTT=63)
            idx = (
                BASE_TO_IDX[trinuc[0]] * 16
                + BASE_TO_IDX[trinuc[1]] * 4
                + BASE_TO_IDX[trinuc[2]]
            )
            counts[idx] += 1.0

    # Compute frequencies.
    total = L - 2
    frequencies = counts / total if total > 0 else counts

    features = np.zeros(64, dtype=np.float32)
    for i, trinuc in enumerate(TRI_KMER_VOCAB):
        eiip_sum = EIIP_VALUES[trinuc[0]] + EIIP_VALUES[trinuc[1]] + EIIP_VALUES[trinuc[2]]
        features[i] = eiip_sum * frequencies[i]

    return torch.from_numpy(features)


def compute_cksnap(sequence: str, max_gap: int = 5) -> torch.FloatTensor:
    """
    Compute CKSNAP (Composition of K-Spaced Nucleotide Acid Pairs).

    For each gap size k in [0, max_gap], count all 16 dinucleotide pairs
    separated by k bases. Counts are normalized to frequencies per gap.
    """
    cleaned = _sanitize_sequence(sequence)
    gap_count = max(0, int(max_gap)) + 1
    counts = np.zeros((gap_count, len(DI_KMER_VOCAB)), dtype=np.float32)
    seq_len = len(cleaned)
    if seq_len < 2:
        return torch.from_numpy(counts.reshape(-1))

    for gap in range(gap_count):
        total_pairs = seq_len - gap - 1
        if total_pairs <= 0:
            continue
        for i in range(total_pairs):
            left = cleaned[i]
            right = cleaned[i + gap + 1]
            pair_idx = BASE_TO_IDX[left] * 4 + BASE_TO_IDX[right]
            counts[gap, pair_idx] += 1.0
        counts[gap] /= float(total_pairs)

    return torch.from_numpy(counts.reshape(-1))


def compute_pstnp(sequence: str) -> torch.FloatTensor:
    """
    Compute PSTNP (Position-Specific Trinucleotide Propensity).

    Uses position-weighted trinucleotide frequencies, emphasizing promoter core
    regions (-10 and -35) relative to the sequence end (TSS).
    """
    cleaned = _sanitize_sequence(sequence)
    counts = np.zeros(len(TRI_KMER_VOCAB), dtype=np.float32)
    seq_len = len(cleaned)
    if seq_len < 3:
        return torch.from_numpy(counts)

    # Core regions for 81bp sequences: -10 ~ 66-76, -35 ~ 41-51 (1-based).
    # Generalize relative to sequence end (TSS at seq_len).
    start_10 = max(0, seq_len - 16)
    end_10 = min(seq_len - 1, seq_len - 6)
    start_35 = max(0, seq_len - 41)
    end_35 = min(seq_len - 1, seq_len - 31)

    total_weight = 0.0
    for i in range(seq_len - 2):
        kmer = cleaned[i : i + 3]
        weight = 1.0
        if (start_10 <= i <= end_10) or (start_35 <= i <= end_35):
            weight = PSTNP_CORE_WEIGHT
        counts[TRI_KMER_TO_IDX[kmer]] += weight
        total_weight += weight

    if total_weight > 0:
        counts /= float(total_weight)
    return torch.from_numpy(counts)


def _load_dataframe(path: Optional[str], delimiter: str) -> Optional[pd.DataFrame]:
    if not path:
        return None
    csv_path = Path(path)
    if not csv_path.exists():
        LOGGER.warning("Dataset file %s not found", path)
        return None
    df = pd.read_csv(csv_path, delimiter=delimiter)
    LOGGER.info("Loaded %d rows from %s", len(df), path)
    return df


def _split_dataframe(df: pd.DataFrame, train_ratio: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not 0.0 < train_ratio < 1.0:
        raise ValueError("train_ratio must be between 0 and 1")
    indices = df.index.to_list()
    rng = np.random.default_rng(RANDOM_SEED)
    rng.shuffle(indices)
    train_cut = max(1, int(round(train_ratio * len(indices))))
    train_indices = indices[:train_cut]
    val_indices = indices[train_cut:]
    if not val_indices:
        val_indices = train_indices[-1:]
        train_indices = train_indices[:-1] or train_indices
    train_df = df.loc[train_indices].reset_index(drop=True)
    val_df = df.loc[val_indices].reset_index(drop=True)
    LOGGER.info("Performed random split: %d train / %d val", len(train_df), len(val_df))
    return train_df, val_df


def _maybe_ablation_subsample(
    df: pd.DataFrame,
    max_samples: Any,
    seed: int,
    label: str,
) -> pd.DataFrame:
    if max_samples is None:
        return df
    try:
        limit = int(max_samples)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} ablation limit must be an int, got {max_samples!r}") from exc
    if limit < 0:
        raise ValueError(f"{label} ablation limit must be >= 0, got {limit}")
    if len(df) <= limit:
        return df
    original_len = len(df)
    df = df.sample(n=limit, random_state=seed).reset_index(drop=True)
    LOGGER.info("Ablation subsample %s: %d -> %d (seed=%d)", label, original_len, len(df), seed)
    return df


def _maybe_sampler(labels: Sequence[float], imbalance_threshold: float = 0.2) -> Optional[WeightedRandomSampler]:
    total = len(labels)
    if total == 0:
        return None
    pos = sum(1 for label in labels if label >= 0.5)
    neg = total - pos
    if pos == 0 or neg == 0:
        return None
    pos_fraction = pos / total
    if abs(pos_fraction - 0.5) <= imbalance_threshold:
        return None
    weights = torch.zeros(total, dtype=torch.double)
    neg_weight = total / (2.0 * neg)
    pos_weight = total / (2.0 * pos)
    for idx, label in enumerate(labels):
        weights[idx] = pos_weight if label >= 0.5 else neg_weight
    LOGGER.info(
        "Using WeightedRandomSampler (pos_frac=%.3f, threshold=%.2f)",
        pos_fraction,
        imbalance_threshold,
    )
    return WeightedRandomSampler(weights, num_samples=total, replacement=True)


def random_base_substitution(sequence: str, prob: float = 0.05) -> str:
    """
    Randomly substitute bases with probability `prob`.
    Simulates sequencing errors and improves robustness.
    """
    if prob <= 0:
        return sequence
    bases = list(sequence)
    for i in range(len(bases)):
        if bases[i] in BASES_SET and random.random() < prob:
            bases[i] = random.choice([b for b in BASES if b != bases[i]])
    return "".join(bases)


def random_base_deletion(sequence: str, prob: float = 0.02) -> str:
    """
    Randomly delete bases with probability `prob`.
    Simulates indel errors.
    """
    if prob <= 0:
        return sequence
    return "".join(b for b in sequence if random.random() >= prob)


def random_base_insertion(sequence: str, prob: float = 0.02) -> str:
    """
    Randomly insert bases with probability `prob`.
    Simulates indel errors.
    """
    if prob <= 0:
        return sequence
    result = []
    for b in sequence:
        if random.random() < prob:
            result.append(random.choice(BASES))
        result.append(b)
    return "".join(result)


class PromoterDatasetStage1(Dataset):
    """
    Stage 1 dataset with enhanced augmentation for 90%+ accuracy.

    Augmentation strategies:
    1. Reverse complement (50% prob) - DNA is double-stranded
    2. Random base substitution (5% per base) - Sequencing error simulation
    3. Random indels (2% per base) - Optional, more aggressive augmentation

    Features: TNC(64) + PseEIIP(64) + CKSNAP(96) + PSTNP(64) = 288

    Performance optimization:
    - ALWAYS pre-caches tokens for both forward and reverse complement sequences
    - Pre-caches engineered features for both orientations
    - At __getitem__ time, we only do random selection and tensor copying (no computation)
    - This makes data loading ~10x faster than on-the-fly tokenization
    """
    ENGINEERED_DIM = 64 + 64 + 96 + 64

    def __init__(
        self,
        df: pd.DataFrame,
        max_token_len: int,
        vocab: DNABPETokenizer,
        use_engineered_features: bool = True,
        feature_enable_tnc: bool = True,
        feature_enable_pseeiip: bool = True,
        feature_enable_cksnap: bool = True,
        feature_enable_pstnp: bool = True,
        engineered_dim: int = ENGINEERED_DIM,
        reverse_complement_prob: float = 0.5,
        base_substitution_prob: float = 0.0,
        base_indel_prob: float = 0.0,
        cache_engineered_features: bool = False,
        cache_tokens: bool = False,
    ):
        if "sequence" not in df or "is_promoter" not in df:
            raise ValueError("Stage 1 data requires 'sequence' and 'is_promoter' columns")
        self.sequences = df["sequence"].astype(str).tolist()
        self.labels = df["is_promoter"].astype(float).tolist()
        self.max_token_len = max_token_len
        self.vocab = vocab
        self.use_engineered_features = use_engineered_features
        self.feature_enable_tnc = bool(feature_enable_tnc)
        self.feature_enable_pseeiip = bool(feature_enable_pseeiip)
        self.feature_enable_cksnap = bool(feature_enable_cksnap)
        self.feature_enable_pstnp = bool(feature_enable_pstnp)
        self.engineered_dim = int(engineered_dim)
        self.reverse_complement_prob = max(0.0, min(1.0, reverse_complement_prob))
        self.base_substitution_prob = max(0.0, min(0.2, base_substitution_prob))
        self.base_indel_prob = max(0.0, min(0.1, base_indel_prob))
        self._zero_engineered = torch.zeros(self.engineered_dim, dtype=torch.float32)

        # Legacy parameters (ignored, kept for API compatibility)
        self.cache_engineered_features = cache_engineered_features
        self.cache_tokens = cache_tokens

        # ALWAYS pre-cache for performance (regardless of cache_* flags)
        # We cache both forward and reverse complement versions
        self._cached_tokens_fwd: Optional[torch.Tensor] = None
        self._cached_tokens_rev: Optional[torch.Tensor] = None
        self._cached_engineered_fwd: Optional[torch.Tensor] = None
        self._cached_engineered_rev: Optional[torch.Tensor] = None

        # Pre-compute everything at initialization
        self._precompute_all()

    def _precompute_all(self) -> None:
        """Pre-compute tokens and engineered features for all sequences.

        Caches both forward and reverse complement versions so that
        __getitem__ only needs to do random selection and tensor indexing.
        """
        n = len(self.sequences)
        LOGGER.info(
            "Pre-computing tokens and features for %d sequences "
            "(forward + reverse complement)...",
            n,
        )

        # Pre-allocate tensors for efficiency
        tokens_fwd = torch.zeros((n, self.max_token_len), dtype=torch.long)
        tokens_rev = torch.zeros((n, self.max_token_len), dtype=torch.long)
        engineered_fwd = torch.zeros((n, self.engineered_dim), dtype=torch.float32)
        engineered_rev = torch.zeros((n, self.engineered_dim), dtype=torch.float32)

        for i, seq in enumerate(self.sequences):
            # Forward sequence
            seq_clean = _sanitize_sequence(seq)
            tokens_fwd[i] = self.vocab.tokenize_and_pad(seq_clean, self.max_token_len)
            engineered_fwd[i] = self._compute_engineered_features(seq_clean)

            # Reverse complement
            seq_rev = reverse_complement(seq_clean)
            tokens_rev[i] = self.vocab.tokenize_and_pad(seq_rev, self.max_token_len)
            engineered_rev[i] = self._compute_engineered_features(seq_rev)

            if (i + 1) % 10000 == 0:
                LOGGER.info("  Pre-computed %d/%d sequences...", i + 1, n)

        self._cached_tokens_fwd = tokens_fwd
        self._cached_tokens_rev = tokens_rev
        self._cached_engineered_fwd = engineered_fwd
        self._cached_engineered_rev = engineered_rev

        LOGGER.info(
            "Pre-computation complete. Tokens: %s, Features: %s",
            tokens_fwd.shape,
            engineered_fwd.shape,
        )

    def _precompute_engineered_features(self) -> None:
        """Legacy method - no longer used, pre-computation happens in _precompute_all."""
        pass

    def _precompute_tokens(self) -> None:
        """Legacy method - no longer used, pre-computation happens in _precompute_all."""
        pass

    def __len__(self) -> int:
        return len(self.sequences)

    def _build_engineered_features(self, sequence: str) -> torch.FloatTensor:
        if not self.use_engineered_features:
            engineered = self._zero_engineered.clone()
            assert engineered.shape[-1] == self.engineered_dim, (
                f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
                f"expected {self.engineered_dim}"
            )
            return engineered
        engineered = self._compute_engineered_features(sequence)
        assert engineered.shape[-1] == self.engineered_dim, (
            f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
            f"expected {self.engineered_dim}"
        )
        return engineered

    def _compute_engineered_features(self, sequence: str) -> torch.FloatTensor:
        """Compute all engineered features for a sequence.

        Note: Feature masking is applied at __getitem__ time, not here,
        so we always compute all features during pre-computation.
        """
        tnc = compute_tnc(sequence)           # 64-dim
        pseeiip = compute_pseeiip(sequence)   # 64-dim
        cksnap = compute_cksnap(sequence)     # 96-dim
        pstnp = compute_pstnp(sequence)       # 64-dim
        return torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)

    def __getitem__(self, idx: int):
        # Fast path: use pre-cached tensors
        # Randomly select forward or reverse complement based on probability
        use_reverse = (
            self.reverse_complement_prob > 0
            and random.random() < self.reverse_complement_prob
        )

        if use_reverse:
            tokens = self._cached_tokens_rev[idx]
            engineered = self._cached_engineered_rev[idx]
        else:
            tokens = self._cached_tokens_fwd[idx]
            engineered = self._cached_engineered_fwd[idx]

        # Apply base substitution if enabled (rare, usually disabled)
        # This is the only case where we need on-the-fly computation
        if self.base_substitution_prob > 0 or self.base_indel_prob > 0:
            sequence = self.sequences[idx]
            if use_reverse:
                sequence = reverse_complement(sequence)
            sequence = self._apply_base_mutations(sequence)
            tokens = self.vocab.tokenize_and_pad(sequence, self.max_token_len)
            engineered = self._compute_engineered_features(sequence)

        # Apply feature masking if any features are disabled
        if not self.use_engineered_features:
            engineered = self._zero_engineered
        elif not (
            self.feature_enable_tnc
            and self.feature_enable_pseeiip
            and self.feature_enable_cksnap
            and self.feature_enable_pstnp
        ):
            # Need to apply feature masking
            engineered = self._apply_feature_mask(engineered)

        label = torch.tensor(float(self.labels[idx]), dtype=torch.float32)
        return tokens, engineered, label

    def _apply_base_mutations(self, sequence: str) -> str:
        """Apply base substitution and indel mutations (rarely used)."""
        if self.base_substitution_prob > 0:
            sequence = random_base_substitution(sequence, self.base_substitution_prob)
        if self.base_indel_prob > 0:
            if random.random() < 0.5:
                sequence = random_base_deletion(sequence, self.base_indel_prob)
            else:
                sequence = random_base_insertion(sequence, self.base_indel_prob)
        return sequence

    def _apply_feature_mask(self, engineered: torch.Tensor) -> torch.Tensor:
        """Apply feature masking for disabled feature types."""
        result = engineered.clone()
        if not self.feature_enable_tnc:
            result[0:64] = 0.0
        if not self.feature_enable_pseeiip:
            result[64:128] = 0.0
        if not self.feature_enable_cksnap:
            result[128:224] = 0.0
        if not self.feature_enable_pstnp:
            result[224:288] = 0.0
        return result


class PromoterDatasetStage2(Dataset):
    """
    Stage 2 dataset for strong vs weak promoter classification.

    Performance optimization:
    - ALWAYS pre-caches tokens and engineered features at initialization
    - No augmentation in Stage 2, so caching is fully effective
    - At __getitem__ time, only tensor indexing is needed (no computation)
    """

    def __init__(
        self,
        df: pd.DataFrame,
        max_token_len: int,
        positive_strength: float,
        negative_strength: float,
        vocab: DNABPETokenizer,
        use_engineered_features: bool = True,
        feature_enable_tnc: bool = True,
        feature_enable_pseeiip: bool = True,
        feature_enable_cksnap: bool = True,
        feature_enable_pstnp: bool = True,
        engineered_dim: int = PromoterDatasetStage1.ENGINEERED_DIM,
        cache_engineered_features: bool = False,
        cache_tokens: bool = False,
    ):
        if "sequence" not in df or "strength" not in df:
            raise ValueError("Stage 2 data requires 'sequence' and 'strength' columns")
        self.sequences = df["sequence"].astype(str).tolist()
        self.labels = [self._to_label(val, positive_strength, negative_strength) for val in df["strength"].tolist()]
        self.max_token_len = max_token_len
        self.vocab = vocab
        self.engineered_dim = int(engineered_dim)
        self.use_engineered_features = bool(use_engineered_features)
        self.feature_enable_tnc = bool(feature_enable_tnc)
        self.feature_enable_pseeiip = bool(feature_enable_pseeiip)
        self.feature_enable_cksnap = bool(feature_enable_cksnap)
        self.feature_enable_pstnp = bool(feature_enable_pstnp)
        self._zero_engineered = torch.zeros(self.engineered_dim, dtype=torch.float32)

        # Legacy parameters (ignored, kept for API compatibility)
        self.cache_engineered_features = cache_engineered_features
        self.cache_tokens = cache_tokens

        # ALWAYS pre-cache for performance
        self._cached_tokens: Optional[torch.Tensor] = None
        self._cached_engineered: Optional[torch.Tensor] = None

        # Pre-compute everything at initialization
        self._precompute_all()

    def _precompute_all(self) -> None:
        """Pre-compute tokens and engineered features for all sequences."""
        n = len(self.sequences)
        LOGGER.info("Pre-computing Stage2 tokens and features for %d sequences...", n)

        # Pre-allocate tensors for efficiency
        tokens = torch.zeros((n, self.max_token_len), dtype=torch.long)
        engineered = torch.zeros((n, self.engineered_dim), dtype=torch.float32)

        for i, seq in enumerate(self.sequences):
            seq_clean = _sanitize_sequence(seq)
            tokens[i] = self.vocab.tokenize_and_pad(seq_clean, self.max_token_len)
            engineered[i] = self._compute_engineered_features(seq_clean)

            if (i + 1) % 5000 == 0:
                LOGGER.info("  Pre-computed %d/%d Stage2 sequences...", i + 1, n)

        self._cached_tokens = tokens
        self._cached_engineered = engineered

        LOGGER.info(
            "Stage2 pre-computation complete. Tokens: %s, Features: %s",
            tokens.shape,
            engineered.shape,
        )

    def _precompute_engineered_features(self) -> None:
        """Legacy method - no longer used, pre-computation happens in _precompute_all."""
        pass

    def _precompute_tokens(self) -> None:
        """Legacy method - no longer used, pre-computation happens in _precompute_all."""
        pass

    def _compute_engineered_features(self, sequence: str) -> torch.FloatTensor:
        """Compute all engineered features for a sequence."""
        tnc = compute_tnc(sequence)           # 64-dim
        pseeiip = compute_pseeiip(sequence)   # 64-dim
        cksnap = compute_cksnap(sequence)     # 96-dim
        pstnp = compute_pstnp(sequence)       # 64-dim
        return torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)

    @staticmethod
    def _to_label(value, positive_strength: float, negative_strength: float) -> float:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            as_text = str(value).strip().lower()
            if as_text.startswith("strong"):
                numeric = float(positive_strength)
            elif as_text.startswith("weak"):
                numeric = float(negative_strength)
            else:
                raise ValueError(f"Unrecognized strength label: {value}")
        if math.isclose(numeric, float(positive_strength), rel_tol=1e-5):
            return 1.0
        if math.isclose(numeric, float(negative_strength), rel_tol=1e-5):
            return 0.0
        return float(numeric)

    def __len__(self) -> int:
        return len(self.sequences)

    def _apply_feature_mask(self, engineered: torch.Tensor) -> torch.Tensor:
        """Apply feature masking for disabled feature types."""
        result = engineered.clone()
        if not self.feature_enable_tnc:
            result[0:64] = 0.0
        if not self.feature_enable_pseeiip:
            result[64:128] = 0.0
        if not self.feature_enable_cksnap:
            result[128:224] = 0.0
        if not self.feature_enable_pstnp:
            result[224:288] = 0.0
        return result

    def __getitem__(self, idx: int):
        # Fast path: use pre-cached tensors
        tokens = self._cached_tokens[idx]
        engineered = self._cached_engineered[idx]

        # Apply feature masking if needed
        if not self.use_engineered_features:
            engineered = self._zero_engineered
        elif not (
            self.feature_enable_tnc
            and self.feature_enable_pseeiip
            and self.feature_enable_cksnap
            and self.feature_enable_pstnp
        ):
            engineered = self._apply_feature_mask(engineered)

        label = torch.tensor(float(self.labels[idx]), dtype=torch.float32)
        return tokens, engineered, label


def build_dataloaders(cfg) -> Dict[str, Dict[str, Optional[DataLoader]]]:
    max_token_len = int(_cfg_value(cfg, "max_token_len", 24))
    batch_size = int(_cfg_value(cfg, "batch_size", 64))
    delimiter = _cfg_value(cfg, "delimiter", "\t")
    train_val_split = float(_cfg_value(cfg, "train_val_split", 0.9))
    # DataLoader performance options
    num_workers = int(_cfg_value(cfg, "num_workers", 0))
    persistent_workers = bool(_cfg_value(cfg, "persistent_workers", False)) and num_workers > 0
    prefetch_factor_raw = _cfg_value(cfg, "prefetch_factor", None)
    prefetch_factor = int(prefetch_factor_raw) if prefetch_factor_raw is not None and num_workers > 0 else None
    pin_memory = bool(_cfg_value(cfg, "pin_memory", False)) and torch.cuda.is_available()

    if num_workers > 0:
        LOGGER.info(
            "DataLoader config: num_workers=%d, persistent_workers=%s, prefetch_factor=%s, pin_memory=%s",
            num_workers,
            persistent_workers,
            prefetch_factor,
            pin_memory,
        )

    # Caching options for performance optimization
    cache_engineered_features = bool(_cfg_value(cfg, "cache_engineered_features", False))
    cache_tokens = bool(_cfg_value(cfg, "cache_tokens", False))
    if cache_engineered_features or cache_tokens:
        LOGGER.info(
            "Caching config: cache_engineered_features=%s, cache_tokens=%s",
            cache_engineered_features,
            cache_tokens,
        )

    stage1_use_engineered = bool(_cfg_value(cfg, "stage1_use_engineered_features", True))
    stage1_feature_enable_tnc = bool(_cfg_value(cfg, "stage1_feature_enable_tnc", True))
    stage1_feature_enable_pseeiip = bool(_cfg_value(cfg, "stage1_feature_enable_pseeiip", True))
    stage1_feature_enable_cksnap = bool(_cfg_value(cfg, "stage1_feature_enable_cksnap", True))
    stage1_feature_enable_pstnp = bool(_cfg_value(cfg, "stage1_feature_enable_pstnp", True))
    stage2_use_engineered = bool(_cfg_value(cfg, "stage2_use_engineered_features", True))
    stage2_feature_enable_tnc = bool(_cfg_value(cfg, "stage2_feature_enable_tnc", stage1_feature_enable_tnc))
    stage2_feature_enable_pseeiip = bool(_cfg_value(cfg, "stage2_feature_enable_pseeiip", stage1_feature_enable_pseeiip))
    stage2_feature_enable_cksnap = bool(_cfg_value(cfg, "stage2_feature_enable_cksnap", stage1_feature_enable_cksnap))
    stage2_feature_enable_pstnp = bool(_cfg_value(cfg, "stage2_feature_enable_pstnp", stage1_feature_enable_pstnp))
    engineered_dim = int(_cfg_value(cfg, "engineered_dim", PromoterDatasetStage1.ENGINEERED_DIM))
    stage1_rc_prob = float(_cfg_value(cfg, "stage1_reverse_complement_prob", 0.5))

    # Load BPE vocabulary
    bpe_vocab_path = _cfg_value(cfg, "bpe_vocab_path", "../artifacts/vocabs/bpe_vocab.json")
    bpe_path = Path(bpe_vocab_path)
    if not bpe_path.exists():
        raise FileNotFoundError(
            f"BPE vocabulary not found at {bpe_path}. "
            "Run train_bpe_vocab.py first to generate it."
        )
    stage1_vocab = DNABPETokenizer.load(str(bpe_path))
    LOGGER.info("Loaded BPE vocab from %s (vocab_size=%d)", bpe_path, len(stage1_vocab))

    ablation_seed = int(_cfg_value(cfg, "ablation_seed", 1337))
    ablation_cfg = _cfg_value(cfg, "ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_enabled = bool(ablation_cfg.get("enabled", False))
    stage1_train_limit = ablation_cfg.get("stage1_max_train_samples") if ablation_enabled else None
    stage1_val_limit = ablation_cfg.get("stage1_max_val_samples") if ablation_enabled else None
    if stage1_train_limit is None:
        stage1_train_limit = _cfg_value(cfg, "ablation_max_train_samples_stage1")
    if stage1_val_limit is None:
        stage1_val_limit = _cfg_value(cfg, "ablation_max_val_samples_stage1")
    stage2_train_limit = ablation_cfg.get("stage2_max_train_samples") if ablation_enabled else None
    stage2_val_limit = ablation_cfg.get("stage2_max_val_samples") if ablation_enabled else None
    if stage2_train_limit is None:
        stage2_train_limit = _cfg_value(cfg, "ablation_max_train_samples_stage2")
    if stage2_val_limit is None:
        stage2_val_limit = _cfg_value(cfg, "ablation_max_val_samples_stage2")

    stage1_train_df = _load_dataframe(_cfg_value(cfg, "stage1_train"), delimiter)
    if stage1_train_df is None:
        raise FileNotFoundError("Stage 1 training data is required")
    stage1_val_df = _load_dataframe(_cfg_value(cfg, "stage1_val"), delimiter)
    if stage1_val_df is None:
        stage1_train_df, stage1_val_df = _split_dataframe(stage1_train_df, train_val_split)
    stage1_train_df = _maybe_ablation_subsample(
        stage1_train_df,
        stage1_train_limit,
        ablation_seed,
        "stage1 train",
    )
    stage1_val_df = _maybe_ablation_subsample(
        stage1_val_df,
        stage1_val_limit,
        ablation_seed,
        "stage1 val",
    )

    stage1_train_ds = PromoterDatasetStage1(
        stage1_train_df,
        max_token_len,
        vocab=stage1_vocab,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        feature_enable_cksnap=stage1_feature_enable_cksnap,
        feature_enable_pstnp=stage1_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        reverse_complement_prob=stage1_rc_prob,
        cache_engineered_features=False,
        cache_tokens=False,
    )
    stage1_val_ds = PromoterDatasetStage1(
        stage1_val_df,
        max_token_len,
        vocab=stage1_vocab,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        feature_enable_cksnap=stage1_feature_enable_cksnap,
        feature_enable_pstnp=stage1_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.0,
        cache_engineered_features=cache_engineered_features,
        cache_tokens=cache_tokens,
    )
    stage1_sampler = _maybe_sampler(stage1_train_ds.labels)
    # Build DataLoader kwargs with optional performance parameters
    loader_kwargs: Dict[str, Any] = {
        "num_workers": num_workers,
        "pin_memory": pin_memory,
    }
    if persistent_workers:
        loader_kwargs["persistent_workers"] = True
    if prefetch_factor is not None:
        loader_kwargs["prefetch_factor"] = prefetch_factor

    stage1_train_loader = DataLoader(
        stage1_train_ds,
        batch_size=batch_size,
        shuffle=stage1_sampler is None,
        sampler=stage1_sampler,
        **loader_kwargs,
    )
    stage1_val_loader = DataLoader(
        stage1_val_ds,
        batch_size=batch_size,
        shuffle=False,
        **loader_kwargs,
    )

    stage2_train_df = _load_dataframe(_cfg_value(cfg, "stage2_train"), delimiter)
    if stage2_train_df is None:
        raise FileNotFoundError("Stage 2 training data is required")
    stage2_val_df = _load_dataframe(_cfg_value(cfg, "stage2_val"), delimiter)
    if stage2_val_df is None:
        stage2_train_df, stage2_val_df = _split_dataframe(stage2_train_df, train_val_split)
    stage2_train_df = _maybe_ablation_subsample(
        stage2_train_df,
        stage2_train_limit,
        ablation_seed,
        "stage2 train",
    )
    stage2_val_df = _maybe_ablation_subsample(
        stage2_val_df,
        stage2_val_limit,
        ablation_seed,
        "stage2 val",
    )
    stage2_test_path = _cfg_value(cfg, "stage2_test")
    stage2_test_df = None
    if stage2_test_path:
        if Path(stage2_test_path).exists():
            stage2_test_df = _load_dataframe(stage2_test_path, delimiter)
        else:
            LOGGER.info(
                "Stage 2 test set not found at %s; reusing stage2_val for test loader",
                stage2_test_path,
            )
            stage2_test_df = stage2_val_df
    else:
        LOGGER.info("Stage 2 test set not configured; reusing stage2_val for test loader")
        stage2_test_df = stage2_val_df

    pos_strength = float(_cfg_value(cfg, "stage2_strength_positive", 1))
    neg_strength = float(_cfg_value(cfg, "stage2_strength_negative", 0))

    stage2_train_ds = PromoterDatasetStage2(
        stage2_train_df,
        max_token_len,
        pos_strength,
        neg_strength,
        stage1_vocab,
        use_engineered_features=stage2_use_engineered,
        feature_enable_tnc=stage2_feature_enable_tnc,
        feature_enable_pseeiip=stage2_feature_enable_pseeiip,
        feature_enable_cksnap=stage2_feature_enable_cksnap,
        feature_enable_pstnp=stage2_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        cache_engineered_features=cache_engineered_features,
        cache_tokens=cache_tokens,
    )
    stage2_val_ds = PromoterDatasetStage2(
        stage2_val_df,
        max_token_len,
        pos_strength,
        neg_strength,
        stage1_vocab,
        use_engineered_features=stage2_use_engineered,
        feature_enable_tnc=stage2_feature_enable_tnc,
        feature_enable_pseeiip=stage2_feature_enable_pseeiip,
        feature_enable_cksnap=stage2_feature_enable_cksnap,
        feature_enable_pstnp=stage2_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        cache_engineered_features=cache_engineered_features,
        cache_tokens=cache_tokens,
    )
    stage2_test_ds = (
        PromoterDatasetStage2(
            stage2_test_df,
            max_token_len,
            pos_strength,
            neg_strength,
            stage1_vocab,
            use_engineered_features=stage2_use_engineered,
            feature_enable_tnc=stage2_feature_enable_tnc,
            feature_enable_pseeiip=stage2_feature_enable_pseeiip,
            feature_enable_cksnap=stage2_feature_enable_cksnap,
            feature_enable_pstnp=stage2_feature_enable_pstnp,
            engineered_dim=engineered_dim,
            cache_engineered_features=cache_engineered_features,
            cache_tokens=cache_tokens,
        )
        if stage2_test_df is not None
        else None
    )
    stage2_sampler = _maybe_sampler(stage2_train_ds.labels)

    stage2_train_loader = DataLoader(
        stage2_train_ds,
        batch_size=batch_size,
        shuffle=stage2_sampler is None,
        sampler=stage2_sampler,
        **loader_kwargs,
    )
    stage2_val_loader = DataLoader(
        stage2_val_ds,
        batch_size=batch_size,
        shuffle=False,
        **loader_kwargs,
    )
    stage2_test_loader = (
        DataLoader(
            stage2_test_ds,
            batch_size=batch_size,
            shuffle=False,
            **loader_kwargs,
        )
        if stage2_test_ds is not None
        else None
    )

    return {
        "stage1": {"train": stage1_train_loader, "val": stage1_val_loader},
        "stage2": {
            "train": stage2_train_loader,
            "val": stage2_val_loader,
            "test": stage2_test_loader,
        },
    }


if __name__ == "__main__":
    # Test PseEIIP fix
    test_seq = "TTGACAATTTTTCTTGATAATGTAACTCACTTAATCTTGATAAATGCTATAATGTGTCGAAAAAAAAAAAAAAAAAAAA"  # 81bp
    pseeiip = compute_pseeiip(test_seq)
    non_zero = (pseeiip != 0).sum().item()
    print(f"PseEIIP test: {non_zero} non-zero elements (should be 25-40, NOT 3)")
    assert non_zero > 20, f"FAILED: PseEIIP still broken, only {non_zero} non-zero elements"
    assert pseeiip.shape[0] == 64, f"FAILED: Wrong shape {pseeiip.shape}"
    print("✓ PseEIIP fix validated successfully!")

    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise SystemExit("PyYAML is required for this sanity check (pip install pyyaml)") from exc

    here = Path(__file__).resolve().parent
    config_path = here / "config.yaml"
    if not config_path.exists():
        raise SystemExit(f"Config file not found: {config_path}")

    os.chdir(here)
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    limit_keys = (
        "ablation_max_train_samples_stage1",
        "ablation_max_val_samples_stage1",
        "ablation_max_train_samples_stage2",
        "ablation_max_val_samples_stage2",
    )
    if not any(_cfg_value(cfg, key) is not None for key in limit_keys):
        print("No ablation_max_* limits set; add them to config.yaml to sanity-check subsampling.")
        raise SystemExit(0)

    cfg_no_limits = dict(cfg)
    for key in limit_keys:
        cfg_no_limits.pop(key, None)

    original = build_dataloaders(cfg_no_limits)
    limited = build_dataloaders(cfg)

    for stage in ("stage1", "stage2"):
        orig_train = len(original[stage]["train"].dataset)
        orig_val = len(original[stage]["val"].dataset)
        lim_train = len(limited[stage]["train"].dataset)
        lim_val = len(limited[stage]["val"].dataset)
        print(f"{stage}: train {orig_train} -> {lim_train}, val {orig_val} -> {lim_val}")
