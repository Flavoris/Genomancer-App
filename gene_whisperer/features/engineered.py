"""PROCABLES-inspired engineered feature extraction."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence

import numpy as np

DNA_BASES = "ACGT"
TRIMERS = [a + b + c for a in DNA_BASES for b in DNA_BASES for c in DNA_BASES]
TRIMER_INDEX = {tri: idx for idx, tri in enumerate(TRIMERS)}
PAIR_INDEX = {a + b: idx for idx, (a, b) in enumerate((x, y) for x in DNA_BASES for y in DNA_BASES)}

EIIP_VALUES = {
    "A": 0.1260,
    "C": 0.1340,
    "G": 0.0806,
    "T": 0.1335,
}


def _normalize_sequence(sequence: str) -> str:
    return "".join(base for base in sequence.upper().strip() if base in (DNA_BASES + "N"))


def compute_cksnap(sequence: str, kmax: int = 2) -> np.ndarray:
    """Composition of k-spaced nucleotide pairs (CKSNAP)."""
    seq = _normalize_sequence(sequence)
    features: List[float] = []
    for k in range(kmax + 1):
        counts = np.zeros(len(PAIR_INDEX), dtype=np.float32)
        total = 0
        for i in range(len(seq) - k - 1):
            pair = seq[i] + seq[i + k + 1]
            if "N" in pair or pair not in PAIR_INDEX:
                continue
            counts[PAIR_INDEX[pair]] += 1
            total += 1
        if total == 0:
            features.extend(counts.tolist())
        else:
            features.extend((counts / total).tolist())
    return np.asarray(features, dtype=np.float32)


def compute_tnc(sequence: str) -> np.ndarray:
    """Trinucleotide composition frequencies (64D)."""
    seq = _normalize_sequence(sequence)
    counts = np.zeros(len(TRIMERS), dtype=np.float32)
    total = 0
    for i in range(len(seq) - 2):
        tri = seq[i : i + 3]
        if "N" in tri:
            continue
        idx = TRIMER_INDEX.get(tri)
        if idx is None:
            continue
        counts[idx] += 1
        total += 1
    if total == 0:
        return counts
    return counts / total


def compute_pseeiip(sequence: str) -> np.ndarray:
    """PseEIIP descriptor: EIIP-weighted trinucleotide frequencies (64D)."""
    seq = _normalize_sequence(sequence)
    counts = np.zeros(len(TRIMERS), dtype=np.float32)
    total = 0
    for i in range(len(seq) - 2):
        tri = seq[i : i + 3]
        if "N" in tri:
            continue
        idx = TRIMER_INDEX.get(tri)
        if idx is None:
            continue
        eiip = EIIP_VALUES[tri[0]] + EIIP_VALUES[tri[1]] + EIIP_VALUES[tri[2]]
        counts[idx] += eiip
        total += 1
    if total == 0:
        return counts
    return counts / total


@dataclass
class PSTNPMatrix:
    """Position-specific trinucleotide propensity (PSTNP) matrix."""

    log_odds: np.ndarray
    sequence_length: int

    @classmethod
    def from_sequences(
        cls,
        positives: Sequence[str],
        negatives: Sequence[str],
        smoothing: float = 1.0,
    ) -> "PSTNPMatrix":
        if not positives or not negatives:
            raise ValueError("PSTNP requires both positive and negative sequences")

        length = len(_normalize_sequence(positives[0]))
        if length < 3:
            raise ValueError("Sequences must be at least length 3 for PSTNP")

        def _check_lengths(seqs: Sequence[str]) -> None:
            for seq in seqs:
                if len(_normalize_sequence(seq)) != length:
                    raise ValueError("All sequences must share the same length")

        _check_lengths(positives)
        _check_lengths(negatives)

        positions = length - 2
        pos_counts = np.zeros((positions, len(TRIMERS)), dtype=np.float64)
        neg_counts = np.zeros((positions, len(TRIMERS)), dtype=np.float64)

        for seq in positives:
            seq = _normalize_sequence(seq)
            for pos in range(positions):
                tri = seq[pos : pos + 3]
                if "N" in tri:
                    continue
                pos_counts[pos, TRIMER_INDEX[tri]] += 1

        for seq in negatives:
            seq = _normalize_sequence(seq)
            for pos in range(positions):
                tri = seq[pos : pos + 3]
                if "N" in tri:
                    continue
                neg_counts[pos, TRIMER_INDEX[tri]] += 1

        pos_totals = pos_counts.sum(axis=1, keepdims=True)
        neg_totals = neg_counts.sum(axis=1, keepdims=True)

        pos_probs = (pos_counts + smoothing) / (pos_totals + smoothing * len(TRIMERS))
        neg_probs = (neg_counts + smoothing) / (neg_totals + smoothing * len(TRIMERS))

        log_odds = np.log(pos_probs) - np.log(neg_probs)
        return cls(log_odds=log_odds.astype(np.float32), sequence_length=length)

    @classmethod
    def from_dict(cls, data: dict) -> "PSTNPMatrix":
        log_odds = np.asarray(data["log_odds"], dtype=np.float32)
        return cls(log_odds=log_odds, sequence_length=int(data["sequence_length"]))

    def to_dict(self) -> dict:
        return {
            "sequence_length": self.sequence_length,
            "log_odds": self.log_odds.tolist(),
        }

    def transform(self, sequence: str) -> np.ndarray:
        seq = _normalize_sequence(sequence)
        positions = self.sequence_length - 2
        if len(seq) < self.sequence_length:
            seq = seq.ljust(self.sequence_length, "N")
        if len(seq) > self.sequence_length:
            seq = seq[: self.sequence_length]

        values = np.zeros(positions, dtype=np.float32)
        for pos in range(positions):
            tri = seq[pos : pos + 3]
            if "N" in tri:
                continue
            values[pos] = self.log_odds[pos, TRIMER_INDEX[tri]]
        return values


@dataclass
class FeatureExtractor:
    """Bundle engineered features used alongside the transformer."""

    use_cksnap: bool = True
    use_tnc: bool = True
    use_pseeiip: bool = True
    use_pstnp: bool = True
    use_gc: bool = True
    kmax: int = 2
    pstnp_matrix: PSTNPMatrix | None = None

    def fit_pstnp(self, positives: Sequence[str], negatives: Sequence[str]) -> None:
        if not self.use_pstnp:
            return
        self.pstnp_matrix = PSTNPMatrix.from_sequences(positives, negatives)

    def feature_dim(self) -> int:
        dim = 0
        if self.use_cksnap:
            dim += (self.kmax + 1) * len(PAIR_INDEX)
        if self.use_tnc:
            dim += len(TRIMERS)
        if self.use_pseeiip:
            dim += len(TRIMERS)
        if self.use_pstnp:
            dim += (self.pstnp_matrix.sequence_length - 2) if self.pstnp_matrix else 0
        if self.use_gc:
            dim += 2
        return dim

    def transform(self, sequence: str) -> np.ndarray:
        parts: List[np.ndarray] = []
        if self.use_cksnap:
            parts.append(compute_cksnap(sequence, kmax=self.kmax))
        if self.use_tnc:
            parts.append(compute_tnc(sequence))
        if self.use_pseeiip:
            parts.append(compute_pseeiip(sequence))
        if self.use_pstnp:
            if self.pstnp_matrix is None:
                raise ValueError("PSTNP matrix not initialized. Call fit_pstnp first.")
            parts.append(self.pstnp_matrix.transform(sequence))
        if self.use_gc:
            seq = _normalize_sequence(sequence)
            if seq:
                gc = (seq.count("G") + seq.count("C")) / len(seq)
                at = (seq.count("A") + seq.count("T")) / len(seq)
            else:
                gc = 0.0
                at = 0.0
            parts.append(np.asarray([gc, at], dtype=np.float32))
        if not parts:
            return np.zeros(0, dtype=np.float32)
        return np.concatenate(parts).astype(np.float32)

    def transform_batch(self, sequences: Iterable[str]) -> np.ndarray:
        features = [self.transform(seq) for seq in sequences]
        return np.stack(features, axis=0)
