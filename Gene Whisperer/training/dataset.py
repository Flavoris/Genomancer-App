"""Dataset utilities for Gene Whisperer promoter models."""
from __future__ import annotations

import json
import logging
import math
import os
import random
from itertools import product
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler

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
TRI_KMER_VOCAB: List[str] = ["".join(kmer) for kmer in product(BASES, repeat=3)]
TRI_KMER_TO_IDX: Dict[str, int] = {kmer: idx for idx, kmer in enumerate(TRI_KMER_VOCAB)}
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


def compute_pseeiip(sequence: str, lag: int = 3) -> torch.FloatTensor:
    """
    Compute PseEIIP (Pseudo Electron-Ion Interaction Potential) features.
    
    Uses EIIP values per base and computes auto-correlation of EIIP signal
    with specified lag. Returns a fixed 64-dim vector.
    """
    seq = _sanitize_sequence(sequence)
    L = len(seq)
    if L == 0:
        return torch.zeros(64, dtype=torch.float32)

    # Basic EIIP sequence
    eiip_seq = np.array([EIIP_VALUES.get(b, 0.0) for b in seq], dtype=np.float32)

    # PseEIIP: auto-correlation of EIIP signal with a lag
    features = []
    for d in range(1, lag + 1):
        if L > d:
            corr = (eiip_seq[:-d] * eiip_seq[d:]).mean()
        else:
            corr = 0.0
        features.append(corr)
    
    # Pad to 64 dims
    vec = np.zeros(64, dtype=np.float32)
    vec[:len(features)] = features
    return torch.from_numpy(vec)


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
    
    Features: TNC(64) + PseEIIP(64) = 128
    """
    ENGINEERED_DIM = 64 + 64

    def __init__(
        self,
        df: pd.DataFrame,
        max_bp_len: int,
        vocab: KmerVocabulary,
        use_engineered_features: bool = True,
        feature_enable_tnc: bool = True,
        feature_enable_pseeiip: bool = True,
        engineered_dim: int = ENGINEERED_DIM,
        reverse_complement_prob: float = 0.5,
        base_substitution_prob: float = 0.0,  # Set > 0 for training augmentation
        base_indel_prob: float = 0.0,  # Usually 0 for promoter prediction
    ):
        if "sequence" not in df or "is_promoter" not in df:
            raise ValueError("Stage 1 data requires 'sequence' and 'is_promoter' columns")
        self.sequences = df["sequence"].astype(str).tolist()
        self.labels = df["is_promoter"].astype(float).tolist()
        self.max_bp_len = max_bp_len
        self.vocab = vocab
        self.use_engineered_features = use_engineered_features
        self.feature_enable_tnc = bool(feature_enable_tnc)
        self.feature_enable_pseeiip = bool(feature_enable_pseeiip)
        self.engineered_dim = int(engineered_dim)
        self.reverse_complement_prob = max(0.0, min(1.0, reverse_complement_prob))
        self.base_substitution_prob = max(0.0, min(0.2, base_substitution_prob))
        self.base_indel_prob = max(0.0, min(0.1, base_indel_prob))
        self._zero_engineered = torch.zeros(self.engineered_dim, dtype=torch.float32)

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
        tnc = compute_tnc(sequence)           # 64-dim
        pseeiip = compute_pseeiip(sequence)   # 64-dim
        if not self.feature_enable_tnc:
            tnc = torch.zeros_like(tnc)
        if not self.feature_enable_pseeiip:
            pseeiip = torch.zeros_like(pseeiip)
        engineered = torch.cat([tnc, pseeiip], dim=0)
        assert engineered.shape[-1] == self.engineered_dim, (
            f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
            f"expected {self.engineered_dim}"
        )
        return engineered

    def _augment_sequence(self, sequence: str) -> str:
        """Apply augmentation chain."""
        # Reverse complement
        if self.reverse_complement_prob > 0 and random.random() < self.reverse_complement_prob:
            sequence = _reverse_complement(sequence)
        
        # Base substitution
        if self.base_substitution_prob > 0:
            sequence = random_base_substitution(sequence, self.base_substitution_prob)
        
        # Indels (optional, usually disabled)
        if self.base_indel_prob > 0:
            if random.random() < 0.5:
                sequence = random_base_deletion(sequence, self.base_indel_prob)
            else:
                sequence = random_base_insertion(sequence, self.base_indel_prob)
        
        return sequence

    def __getitem__(self, idx: int):
        sequence = self.sequences[idx]
        sequence = self._augment_sequence(sequence)
        tokens = self.vocab.tokenize(sequence, self.max_bp_len)
        engineered = self._build_engineered_features(sequence)
        label = torch.tensor(float(self.labels[idx]), dtype=torch.float32)
        return tokens, engineered, label


class PromoterDatasetStage2(Dataset):
    def __init__(
        self,
        df: pd.DataFrame,
        max_bp_len: int,
        positive_strength: float,
        negative_strength: float,
        vocab: KmerVocabulary,
        engineered_dim: int = PromoterDatasetStage1.ENGINEERED_DIM,
    ):
        if "sequence" not in df or "strength" not in df:
            raise ValueError("Stage 2 data requires 'sequence' and 'strength' columns")
        self.sequences = df["sequence"].astype(str).tolist()
        self.labels = [self._to_label(val, positive_strength, negative_strength) for val in df["strength"].tolist()]
        self.max_bp_len = max_bp_len
        self.vocab = vocab
        self.engineered_dim = int(engineered_dim)

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

    def __getitem__(self, idx: int):
        sequence = self.sequences[idx]
        tokens = self.vocab.tokenize(sequence, self.max_bp_len)
        tnc = compute_tnc(sequence)
        pseeiip = compute_pseeiip(sequence)
        label = torch.tensor(float(self.labels[idx]), dtype=torch.float32)
        engineered = torch.cat([tnc, pseeiip], dim=0)
        assert engineered.shape[-1] == self.engineered_dim, (
            f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
            f"expected {self.engineered_dim}"
        )
        return tokens, engineered, label


def build_dataloaders(cfg) -> Dict[str, Dict[str, Optional[DataLoader]]]:
    max_bp_len = int(_cfg_value(cfg, "max_bp_len", 81))
    batch_size = int(_cfg_value(cfg, "batch_size", 64))
    delimiter = _cfg_value(cfg, "delimiter", "\t")
    train_val_split = float(_cfg_value(cfg, "train_val_split", 0.9))
    num_workers = int(_cfg_value(cfg, "num_workers", 0))
    # pin_memory only works on CUDA, not MPS
    pin_memory = bool(_cfg_value(cfg, "pin_memory", False)) and torch.cuda.is_available()
    stage1_use_engineered = bool(_cfg_value(cfg, "stage1_use_engineered_features", True))
    stage1_feature_enable_tnc = bool(_cfg_value(cfg, "stage1_feature_enable_tnc", True))
    stage1_feature_enable_pseeiip = bool(_cfg_value(cfg, "stage1_feature_enable_pseeiip", True))
    engineered_dim = int(_cfg_value(cfg, "engineered_dim", PromoterDatasetStage1.ENGINEERED_DIM))
    stage1_rc_prob = float(_cfg_value(cfg, "stage1_reverse_complement_prob", 0.5))
    kmer = int(_cfg_value(cfg, "kmer", 3))
    vocab_cache_dir = Path(_cfg_value(cfg, "vocab_cache_dir", "../artifacts/vocabs")).resolve()
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
    stage1_vocab = _get_or_create_vocab(stage1_train_df["sequence"].astype(str).tolist(), kmer, vocab_cache_dir)

    stage1_train_ds = PromoterDatasetStage1(
        stage1_train_df,
        max_bp_len,
        vocab=stage1_vocab,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        engineered_dim=engineered_dim,
        reverse_complement_prob=stage1_rc_prob,
    )
    stage1_val_ds = PromoterDatasetStage1(
        stage1_val_df,
        max_bp_len,
        vocab=stage1_vocab,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.0,
    )
    stage1_sampler = _maybe_sampler(stage1_train_ds.labels)
    stage1_train_loader = DataLoader(
        stage1_train_ds,
        batch_size=batch_size,
        shuffle=stage1_sampler is None,
        sampler=stage1_sampler,
        num_workers=num_workers,
        pin_memory=pin_memory,
    )
    stage1_val_loader = DataLoader(
        stage1_val_ds,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=pin_memory,
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
        max_bp_len,
        pos_strength,
        neg_strength,
        stage1_vocab,
        engineered_dim=engineered_dim,
    )
    stage2_val_ds = PromoterDatasetStage2(
        stage2_val_df,
        max_bp_len,
        pos_strength,
        neg_strength,
        stage1_vocab,
        engineered_dim=engineered_dim,
    )
    stage2_test_ds = (
        PromoterDatasetStage2(
            stage2_test_df,
            max_bp_len,
            pos_strength,
            neg_strength,
            stage1_vocab,
            engineered_dim=engineered_dim,
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
        num_workers=num_workers,
        pin_memory=pin_memory,
    )
    stage2_val_loader = DataLoader(
        stage2_val_ds,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=pin_memory,
    )
    stage2_test_loader = (
        DataLoader(
            stage2_test_ds,
            batch_size=batch_size,
            shuffle=False,
            num_workers=num_workers,
            pin_memory=pin_memory,
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
