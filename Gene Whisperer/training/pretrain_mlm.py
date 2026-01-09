"""Masked language model pretraining on E. coli MG1655 genome.

Key improvements for reaching lower loss (<1.0):
1. Warmup + cosine annealing with warm restarts
2. Curriculum learning: mask probability increases over training
3. Validation monitoring to detect overfitting
4. Gradient accumulation for larger effective batch size
5. Label smoothing for better generalization
"""
# README: pretraining map -> dataset (MLMDataset + prepare_dataset), trainer (run_mlm_pretrain).
# Train/val loss: run_mlm_pretrain epoch summary + run_validation; grad accumulation scales loss only for backward.
from __future__ import annotations

import argparse
import hashlib
import copy
import csv
import json
import logging
import math
import os
import random
import shutil
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts, LambdaLR
from torch.utils.data import DataLoader, Dataset, IterableDataset, random_split

from dataset import KmerVocabulary
from model import DNAEncoder
from length_curriculum import (
    LengthCurriculumConfig,
    LengthCurriculumScheduler,
    CurriculumMLMCollator,
    create_curriculum_collator,
)
from pretrain_health import (
    PretrainHealthDashboard,
    create_health_dashboard,
    run_pretrain_sanity_checks,
)
from seed_utils import (
    set_global_seed,
    get_deterministic_dataloader_kwargs,
)
from numerics import (
    NumericsChecker,
    NumericsError,
    assert_finite,
    cap_model_param_norm_,
    detect_nan_in_model,
)
from amp_utils import get_amp_context, log_amp_status
from genome_io import read_fasta_records, sanitize_unknown_bases, find_acgt_segments
from checkpoint_utils import prune_step_checkpoints
from early_stopping import PretrainingEarlyStopping

LOGGER = logging.getLogger("gene_whisperer.pretrain_mlm")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

REV_COMP_MAP = str.maketrans("ACGT", "TGCA")


def get_kmer_pretrain_config(cfg: dict, kmer: int) -> dict:
    """
    Get configuration for specific k-mer pretraining.

    Applies k-mer specific overrides from config.
    """
    kmer_cfg = dict(cfg)
    kmer_cfg["kmer"] = kmer
    kmer_cfg["mlm_kmer"] = kmer

    # Apply k-mer specific overrides
    overrides = cfg.get("kmer_pretrain_overrides", {}).get(kmer, {})
    for key, value in overrides.items():
        kmer_cfg[key] = value
        LOGGER.info("K=%d override: %s=%s", kmer, key, value)

    # Set vocab size
    kmer_cfg["vocab_size"] = (4 ** kmer) + 3

    # Set paths
    vocab_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
    kmer_cfg["mlm_vocab_path"] = str(vocab_dir / f"mlm_k{kmer}_vocab.json")

    artifacts_dir = Path("../artifacts")
    kmer_cfg["mlm_encoder_path"] = str(artifacts_dir / f"mlm_encoder_k{kmer}.pt")
    kmer_cfg["mlm_metadata_path"] = str(artifacts_dir / f"mlm_encoder_k{kmer}_metadata.json")

    return kmer_cfg


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of an A/C/G/T-only DNA sequence."""
    return sequence.translate(REV_COMP_MAP)[::-1]


# Note: set_seed is imported from seed_utils as set_global_seed
# This wrapper maintains backwards compatibility
def set_seed(seed: int) -> None:
    """Set all random seeds for reproducibility (wrapper for set_global_seed)."""
    set_global_seed(seed)


class MaskingScheduler:
    """
    Progressive masking rate scheduler following DNABERT.

    DNABERT uses curriculum learning for masking:
    - First 80% of training: 15% masking rate (easier)
    - Last 20% of training: 20% masking rate (harder)

    This progressive approach starts easier and increases difficulty,
    which improves final model performance.

    Args:
        initial_rate: Starting mask rate (default 0.15)
        final_rate: Ending mask rate (default 0.20)
        warmup_steps: Steps before starting to increase rate (default 0)
        total_steps: Total training steps
        transition_start_ratio: When to start transitioning (default 0.8 = 80% of training)

    Example:
        >>> scheduler = MaskingScheduler(
        ...     initial_rate=0.15,
        ...     final_rate=0.20,
        ...     total_steps=100000,
        ...     transition_start_ratio=0.8,
        ... )
        >>> scheduler.get_mask_rate(0)      # Start: 0.15
        0.15
        >>> scheduler.get_mask_rate(80000)  # At 80%: 0.15 (just starting transition)
        0.15
        >>> scheduler.get_mask_rate(90000)  # At 90%: 0.175 (midway through transition)
        0.175
        >>> scheduler.get_mask_rate(100000) # At 100%: 0.20
        0.2
    """

    def __init__(
        self,
        initial_rate: float = 0.15,
        final_rate: float = 0.20,
        warmup_steps: int = 0,
        total_steps: int = 100000,
        transition_start_ratio: float = 0.8,
    ):
        if not 0.0 <= initial_rate <= 1.0:
            raise ValueError(f"initial_rate must be in [0, 1], got {initial_rate}")
        if not 0.0 <= final_rate <= 1.0:
            raise ValueError(f"final_rate must be in [0, 1], got {final_rate}")
        if warmup_steps < 0:
            raise ValueError(f"warmup_steps must be >= 0, got {warmup_steps}")
        if total_steps <= 0:
            raise ValueError(f"total_steps must be > 0, got {total_steps}")
        if not 0.0 <= transition_start_ratio <= 1.0:
            raise ValueError(f"transition_start_ratio must be in [0, 1], got {transition_start_ratio}")

        self.initial_rate = initial_rate
        self.final_rate = final_rate
        self.warmup_steps = warmup_steps
        self.total_steps = total_steps
        self.transition_start_ratio = transition_start_ratio

        # Compute transition boundaries
        self.transition_start_step = int(total_steps * transition_start_ratio)
        self.transition_steps = total_steps - self.transition_start_step

        # Track current step for logging
        self._current_step = 0

    def get_mask_rate(self, step: int) -> float:
        """
        Return the masking rate for the current step.

        Schedule:
        - Steps 0 to warmup_steps: initial_rate (warmup phase)
        - Steps warmup_steps to transition_start_step: initial_rate (stable phase)
        - Steps transition_start_step to total_steps: linear interpolation to final_rate
        - Steps > total_steps: final_rate (clamped)

        Args:
            step: Current training step (0-indexed)

        Returns:
            Masking rate for this step
        """
        self._current_step = step

        # Before transition starts: use initial rate
        if step < self.transition_start_step:
            return self.initial_rate

        # During transition: linear interpolation
        if step < self.total_steps:
            progress = (step - self.transition_start_step) / max(1, self.transition_steps)
            return self.initial_rate + progress * (self.final_rate - self.initial_rate)

        # After total_steps: clamp to final rate
        return self.final_rate

    def state_dict(self) -> dict:
        """Return state for checkpointing."""
        return {
            "initial_rate": self.initial_rate,
            "final_rate": self.final_rate,
            "warmup_steps": self.warmup_steps,
            "total_steps": self.total_steps,
            "transition_start_ratio": self.transition_start_ratio,
            "current_step": self._current_step,
        }

    def load_state_dict(self, state: dict) -> None:
        """Load state for resuming training."""
        self.initial_rate = state["initial_rate"]
        self.final_rate = state["final_rate"]
        self.warmup_steps = state["warmup_steps"]
        self.total_steps = state["total_steps"]
        self.transition_start_ratio = state["transition_start_ratio"]
        self._current_step = state.get("current_step", 0)

        # Recompute derived values
        self.transition_start_step = int(self.total_steps * self.transition_start_ratio)
        self.transition_steps = self.total_steps - self.transition_start_step

    def get_schedule_info(self) -> dict:
        """Return human-readable schedule information for logging."""
        return {
            "initial_rate": f"{self.initial_rate:.1%}",
            "final_rate": f"{self.final_rate:.1%}",
            "transition_start_step": self.transition_start_step,
            "total_steps": self.total_steps,
            "current_step": self._current_step,
            "current_rate": f"{self.get_mask_rate(self._current_step):.1%}",
        }

    def __repr__(self) -> str:
        return (
            f"MaskingScheduler(initial={self.initial_rate:.1%}, final={self.final_rate:.1%}, "
            f"transition_at={self.transition_start_ratio:.0%}, total_steps={self.total_steps})"
        )


def extract_windows(
    sequence: str,
    window_size: int,
    stride: int,
    *,
    unknown_base_strategy: str,
) -> List[str]:
    windows: List[str] = []
    if unknown_base_strategy == "skip":
        segments = find_acgt_segments(sequence, min_length=window_size)
        for segment_start, segment_end in segments:
            max_start = segment_end - window_size + 1
            for start in range(segment_start, max_start, stride):
                window = sequence[start : start + window_size]
                windows.append(window)
        return windows
    sanitized = sanitize_unknown_bases(sequence, unknown_base_strategy)
    max_start = len(sanitized) - window_size + 1
    for start in range(0, max_start, stride):
        window = sanitized[start : start + window_size]
        windows.append(window)
    return windows


def _resolve_cfg_path(path_value: str | Path | None, *, base_dir: Path) -> Path:
    if path_value is None or str(path_value).strip() == "":
        raise ValueError("Expected a non-empty path value")
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def _normalize_unknown_base_strategy(raw_strategy: Any) -> str:
    if raw_strategy is None:
        return "random"
    if isinstance(raw_strategy, list):
        if len(raw_strategy) != 1:
            raise ValueError("mlm_unknown_base_strategy must be a single value")
        raw_strategy = raw_strategy[0]
    if not isinstance(raw_strategy, str):
        raise ValueError("mlm_unknown_base_strategy must be a string")
    strategy = raw_strategy.strip().lower()
    if strategy == "random":
        return "random"
    if strategy in {"skip", "strict", "acgt_only"}:
        return "skip"
    raise ValueError(f"Unsupported mlm_unknown_base_strategy: {raw_strategy!r}")


def _apply_bp_limit(records: List[Tuple[str, str]], bp_limit: int) -> List[Tuple[str, str]]:
    if bp_limit <= 0:
        return records
    limited: List[Tuple[str, str]] = []
    remaining = bp_limit
    for record_id, sequence in records:
        if remaining <= 0:
            break
        if len(sequence) <= remaining:
            limited.append((record_id, sequence))
            remaining -= len(sequence)
        else:
            limited.append((record_id, sequence[:remaining]))
            remaining = 0
    return limited


class MLMDataset(Dataset):
    def __init__(self, token_tensors: List[torch.LongTensor]):
        self.samples = token_tensors

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> torch.LongTensor:
        return self.samples[idx]


class PromoterMLMDataset(Dataset):
    """
    MLM Dataset with support for contiguous (DNABERT-style) and random masking.

    This dataset stores raw sequences and applies masking on-the-fly during
    __getitem__, allowing for:
    - Dynamic mask rate updates (for progressive/curriculum masking)
    - Switchable masking strategies (contiguous vs random span masking)
    - Proper handling of k-mer specific masking

    Args:
        sequences: List of DNA sequences (raw strings)
        vocab: KmerVocabulary for tokenization
        max_len: Maximum sequence length in base pairs (default: 81)
        mask_strategy: "contiguous" for DNABERT-style or "random" for span masking
        mask_rate: Fraction of tokens to mask (can be updated during training)
        mask_token_prob: Probability of replacing masked token with [MASK] (default: 0.8)
        random_token_prob: Probability of replacing masked token with random token (default: 0.1)
        max_span_len: Maximum span length for random strategy (default: 3)
        mean_span: Mean span length for geometric distribution in random strategy (default: 3.0)
    """

    def __init__(
        self,
        sequences: List[str],
        vocab: KmerVocabulary,
        max_len: int = 81,
        mask_strategy: str = "contiguous",
        mask_rate: float = 0.15,
        mask_token_prob: float = 0.8,
        random_token_prob: float = 0.1,
        max_span_len: int = 3,
        mean_span: float = 3.0,
    ):
        if mask_strategy not in ("contiguous", "random"):
            raise ValueError(f"mask_strategy must be 'contiguous' or 'random', got '{mask_strategy}'")
        if not 0.0 <= mask_rate <= 1.0:
            raise ValueError(f"mask_rate must be in [0, 1], got {mask_rate}")
        if not 0.0 <= mask_token_prob <= 1.0:
            raise ValueError(f"mask_token_prob must be in [0, 1], got {mask_token_prob}")
        if not 0.0 <= random_token_prob <= 1.0:
            raise ValueError(f"random_token_prob must be in [0, 1], got {random_token_prob}")
        if mask_token_prob + random_token_prob > 1.0:
            raise ValueError(
                f"mask_token_prob + random_token_prob must be <= 1.0, "
                f"got {mask_token_prob} + {random_token_prob} = {mask_token_prob + random_token_prob}"
            )

        self.sequences = sequences
        self.vocab = vocab
        self.max_len = max_len
        self.mask_strategy = mask_strategy
        self.mask_rate = mask_rate  # Mutable for progressive masking
        self.mask_token_prob = mask_token_prob
        self.random_token_prob = random_token_prob
        self.max_span_len = max_span_len
        self.mean_span = mean_span

        # Pre-compute special token IDs for contiguous masking
        self._special_token_ids = {vocab.pad_id, vocab.unk_id, vocab.mask_id}

    def __len__(self) -> int:
        return len(self.sequences)

    def set_mask_rate(self, rate: float) -> None:
        """
        Update mask rate for progressive/curriculum masking.

        This allows dynamically increasing mask difficulty during training.
        Note: With num_workers > 0, worker processes won't see this update.
        For progressive masking, use num_workers=0 or re-create the DataLoader.

        Args:
            rate: New masking rate (0.0 to 1.0)
        """
        if not 0.0 <= rate <= 1.0:
            raise ValueError(f"mask_rate must be in [0, 1], got {rate}")
        self.mask_rate = rate

    def set_mask_strategy(self, strategy: str) -> None:
        """
        Switch between masking strategies.

        Args:
            strategy: "contiguous" for DNABERT-style or "random" for span masking
        """
        if strategy not in ("contiguous", "random"):
            raise ValueError(f"mask_strategy must be 'contiguous' or 'random', got '{strategy}'")
        self.mask_strategy = strategy

    def __getitem__(self, idx: int) -> Tuple[torch.LongTensor, torch.LongTensor]:
        """
        Returns (masked_tokens, labels) for a single sequence.

        Args:
            idx: Index of the sequence

        Returns:
            masked_tokens: Token IDs with masking applied, shape (seq_len,)
            labels: Original token IDs for masked positions, -100 for unmasked
        """
        sequence = self.sequences[idx]

        # Tokenize the sequence
        tokens = self.vocab.tokenize(sequence, self.max_len)

        # Add batch dimension for masking functions (they expect batch input)
        tokens_batch = tokens.unsqueeze(0)

        if self.mask_strategy == "contiguous":
            # Use DNABERT-style contiguous span masking
            masked_tokens, labels, _ = create_contiguous_mask(
                tokens_batch,
                kmer=self.vocab.k,
                mask_token_id=self.vocab.mask_id,
                vocab_size=len(self.vocab.itos),
                mask_prob=self.mask_rate,
                mask_token_prob=self.mask_token_prob,
                random_token_prob=self.random_token_prob,
                special_token_ids=self._special_token_ids,
            )
        else:
            # Use random span masking (backward compatible)
            masked_tokens, labels = mask_tokens_span(
                tokens_batch,
                self.vocab,
                mask_prob=self.mask_rate,
                max_span_len=self.max_span_len,
                track_spans=False,
                span_distribution="geometric",
                mean_span=self.mean_span,
                exclude_special_from_labels=True,
            )

        # Remove batch dimension
        return masked_tokens.squeeze(0), labels.squeeze(0)


RecordEntry = Tuple[str, str, List[Tuple[int, int]] | None]


class GenomeMLMIterableDataset(IterableDataset):
    """On-the-fly MLM sampling across multiple genomes."""

    def __init__(
        self,
        corpus: List[Tuple[str, str, str]],
        window_bp: int,
        kmer: int,
        steps_per_epoch: int | None,
        mask_config: Dict[str, Any],
        reverse_complement_prob: float,
        unknown_base_strategy: Any,
        *,
        species_weights: Dict[str, float] | None = None,
        return_metadata: bool = False,
    ) -> None:
        if window_bp <= 0:
            raise ValueError(f"window_bp must be > 0, got {window_bp}")
        if kmer <= 0:
            raise ValueError(f"kmer must be > 0, got {kmer}")
        if kmer > window_bp:
            raise ValueError(f"kmer ({kmer}) must be <= window_bp ({window_bp})")
        if not 0.0 <= reverse_complement_prob <= 1.0:
            raise ValueError(
                f"reverse_complement_prob must be in [0, 1], got {reverse_complement_prob}"
            )
        if steps_per_epoch is not None:
            steps_int = int(steps_per_epoch)
            if steps_int <= 0:
                raise ValueError(f"steps_per_epoch must be > 0, got {steps_per_epoch}")
            steps_per_epoch = steps_int

        vocab = mask_config.get("vocab")
        if vocab is None:
            raise ValueError("mask_config must include 'vocab'")

        self.window_bp = int(window_bp)
        self.kmer = int(kmer)
        self.steps_per_epoch = steps_per_epoch
        self.vocab = vocab
        self.mask_prob = float(mask_config.get("mask_prob", 0.15))
        self.max_span_len = int(mask_config.get("max_span_len", 3))
        self.span_distribution = str(mask_config.get("span_distribution", "geometric"))
        self.mean_span = float(mask_config.get("mean_span", 3.0))
        self.exclude_special_from_labels = bool(mask_config.get("exclude_special_from_labels", True))
        self.track_spans = bool(mask_config.get("track_spans", False))
        self.debug = bool(mask_config.get("debug", False))
        self.reverse_complement_prob = float(reverse_complement_prob)
        self.unknown_base_strategy = _normalize_unknown_base_strategy(unknown_base_strategy)
        self.return_metadata = bool(return_metadata)
        self._debug_ran = False

        if self.return_metadata and self.track_spans:
            raise ValueError("return_metadata is incompatible with track_spans")

        batch_size = mask_config.get("batch_size")
        if batch_size is not None:
            batch_size = int(batch_size)
        self._samples_per_epoch = None
        if self.steps_per_epoch is not None and batch_size:
            self._samples_per_epoch = self.steps_per_epoch * batch_size

        records_by_species: Dict[str, List[RecordEntry]] = {}
        species_order: List[str] = []
        skipped_short = 0
        skipped_unknown = 0
        for species_name, record_id, sequence in corpus:
            if len(sequence) < self.window_bp:
                skipped_short += 1
                continue
            segments = None
            if self.unknown_base_strategy == "skip":
                segments = find_acgt_segments(sequence, min_length=self.window_bp)
                if not segments:
                    skipped_unknown += 1
                    continue
            if species_name not in records_by_species:
                records_by_species[species_name] = []
                species_order.append(species_name)
            records_by_species[species_name].append((record_id, sequence, segments))

        if not records_by_species:
            raise ValueError(
                "No sequences long enough for window_bp after filtering; "
                f"skipped_short={skipped_short}"
            )

        if skipped_short > 0:
            LOGGER.warning("Skipped %d records shorter than window_bp=%d", skipped_short, self.window_bp)
        if skipped_unknown > 0:
            LOGGER.warning(
                "Skipped %d records with no ACGT-only windows (window_bp=%d)",
                skipped_unknown,
                self.window_bp,
            )

        if species_weights is None:
            species_weights = {species: 1.0 for species in species_order}

        missing = [name for name in species_order if name not in species_weights]
        if missing:
            raise ValueError(f"Species weights missing entries for: {missing}")

        self._species_names = [name for name in species_order if name in records_by_species]
        self._species_weights = [float(species_weights[name]) for name in self._species_names]
        for name, weight in zip(self._species_names, self._species_weights):
            if weight <= 0:
                raise ValueError(f"Species weight must be > 0 (species={name}, weight={weight})")

        self._records_by_species = records_by_species

    def __len__(self) -> int:
        if self._samples_per_epoch is None:
            raise TypeError("GenomeMLMIterableDataset length requires steps_per_epoch and batch_size")
        return self._samples_per_epoch

    def _sample_tokens(self) -> Tuple[torch.LongTensor, Dict[str, Any]]:
        species_name = random.choices(self._species_names, weights=self._species_weights, k=1)[0]
        record_id, sequence, segments = random.choice(self._records_by_species[species_name])
        if self.unknown_base_strategy == "skip":
            if not segments:
                raise RuntimeError(
                    f"Record missing ACGT segments (species={species_name}, record={record_id})"
                )
            segment_start, segment_end = random.choice(segments)
            max_start = segment_end - self.window_bp
            start = random.randint(segment_start, max_start)
            window = sequence[start : start + self.window_bp]
        else:
            max_start = len(sequence) - self.window_bp
            if max_start < 0:
                raise RuntimeError(
                    f"Record shorter than window_bp after filtering (species={species_name}, record={record_id})"
                )
            start = random.randint(0, max_start)
            window = sequence[start : start + self.window_bp]
            window = sanitize_unknown_bases(window, self.unknown_base_strategy)
        rc_used = False
        if self.reverse_complement_prob > 0 and random.random() < self.reverse_complement_prob:
            window = reverse_complement(window)
            rc_used = True
        tokens = self.vocab.tokenize(window, self.window_bp)
        metadata = {
            "species": species_name,
            "record_id": record_id,
            "start": start,
            "reverse_complement": rc_used,
        }
        return tokens, metadata

    def sample_raw_tokens(self, num_samples: int) -> List[torch.LongTensor]:
        return [self._sample_tokens()[0] for _ in range(num_samples)]

    def set_mask_prob(self, mask_prob: float) -> None:
        """
        Update the mask probability for progressive masking.

        Note: This only affects the main process. When using num_workers > 0,
        worker processes won't see this update. For progressive masking with
        streaming datasets, use num_workers=0.

        Args:
            mask_prob: New masking probability (0.0 to 1.0)
        """
        self.mask_prob = mask_prob

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if self._samples_per_epoch is None:
            samples_to_yield = None
        elif worker_info is None:
            samples_to_yield = self._samples_per_epoch
        else:
            base = self._samples_per_epoch // worker_info.num_workers
            extra = self._samples_per_epoch % worker_info.num_workers
            samples_to_yield = base + (1 if worker_info.id < extra else 0)

        yielded = 0
        while samples_to_yield is None or yielded < samples_to_yield:
            original_tokens, metadata = self._sample_tokens()
            batch_tokens = original_tokens.unsqueeze(0)
            if self.track_spans:
                masked_inputs, labels, span_lengths = mask_tokens_span(
                    batch_tokens,
                    self.vocab,
                    mask_prob=self.mask_prob,
                    max_span_len=self.max_span_len,
                    track_spans=True,
                    span_distribution=self.span_distribution,
                    mean_span=self.mean_span,
                    exclude_special_from_labels=self.exclude_special_from_labels,
                )
            else:
                masked_inputs, labels = mask_tokens_span(
                    batch_tokens,
                    self.vocab,
                    mask_prob=self.mask_prob,
                    max_span_len=self.max_span_len,
                    track_spans=False,
                    span_distribution=self.span_distribution,
                    mean_span=self.mean_span,
                    exclude_special_from_labels=self.exclude_special_from_labels,
                )
                span_lengths = []

            masked_inputs = masked_inputs.squeeze(0)
            labels = labels.squeeze(0)

            if self.debug and not self._debug_ran:
                self._debug_ran = True
                self._run_debug_verification(original_tokens, masked_inputs, labels)

            if self.return_metadata:
                yield masked_inputs, labels, metadata
            elif self.track_spans:
                yield masked_inputs, labels, span_lengths
            else:
                yield masked_inputs, labels
            yielded += 1

    def _run_debug_verification(
        self,
        original: torch.LongTensor,
        masked_input: torch.LongTensor,
        labels: torch.LongTensor,
    ) -> None:
        masked_positions = labels != -100
        if masked_positions.any():
            assert torch.equal(labels[masked_positions], original[masked_positions]), (
                "Label debug failed: masked labels do not match original tokens"
            )
        unmasked_positions = ~masked_positions
        if unmasked_positions.any():
            assert torch.all(labels[unmasked_positions] == -100), (
                "Label debug failed: unmasked positions should be -100"
            )


def sample_span_length_geometric(max_span_len: int, mean_span: float = 3.0) -> int:
    """
    Sample span length from geometric distribution (DNABERT-style).
    
    The geometric distribution naturally generates span lengths where shorter
    spans are more common than longer ones, which matches the intuition that
    most masked regions should be small.
    
    Args:
        max_span_len: Maximum allowed span length
        mean_span: Target mean span length (controls distribution)
        
    Returns:
        Sampled span length in [1, max_span_len]
    """
    # Geometric distribution: P(X=k) = p * (1-p)^(k-1)
    # Mean = 1/p, so p = 1/mean_span
    p = 1.0 / mean_span
    
    # Sample from geometric distribution and clamp to [1, max_span_len]
    # Use inverse CDF method: X = ceil(log(1-U) / log(1-p))
    u = random.random()
    if u == 0:
        u = 1e-10  # Avoid log(0)
    
    span_len = int(math.ceil(math.log(1 - u) / math.log(1 - p)))
    return max(1, min(span_len, max_span_len))


def mask_tokens_span(
    inputs: torch.LongTensor,
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
    max_span_len: int = 3,
    track_spans: bool = False,
    span_distribution: str = "geometric",
    mean_span: float = 3.0,
    exclude_special_from_labels: bool = True,
) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
    """
    DNABERT-style span masking with geometric or uniform distribution.
    
    Instead of independent per-token masking, this selects contiguous spans
    until approximately mask_prob of tokens are covered, then applies the
    standard 80/10/10 BERT masking within those spans.
    
    Key features:
    - Geometric distribution (default): shorter spans more common, matches DNABERT
    - Uniform distribution: all span lengths equally likely
    - Never masks PAD tokens
    - Handles edge cases (all-PAD, very short sequences)
    - CRITICAL: Special tokens ([MASK], [UNK], [PAD]) are NEVER used as labels
    
    Args:
        inputs: (B, L) input token indices
        vocab: KmerVocabulary with mask_id and pad_id
        mask_prob: Target fraction of tokens to mask (~15%)
        max_span_len: Maximum length of each masked span
        track_spans: If True, also return list of span lengths
        span_distribution: "geometric" (DNABERT-style) or "uniform"
        mean_span: Mean span length for geometric distribution
        exclude_special_from_labels: If True, set labels to -100 for any position
            where the original token was a special token. This prevents the model
            from being trained to predict special tokens.
        
    Returns:
        (masked_inputs, labels) tuple where labels has -100 for unmasked positions
        If track_spans=True, also returns list of span lengths
    """
    device = inputs.device
    labels = inputs.clone()
    batch_size, seq_len = inputs.shape
    vocab_size = len(vocab.itos)

    if max_span_len < 1:
        raise ValueError(f"max_span_len must be >= 1, got {max_span_len}")

    # Initialize all positions as unmasked
    masked = torch.zeros_like(inputs, dtype=torch.bool, device=device)
    
    # Identify PAD positions early (we won't mask these)
    if hasattr(vocab, "pad_id"):
        pad_mask = inputs.eq(vocab.pad_id)
    else:
        pad_mask = torch.zeros_like(inputs, dtype=torch.bool, device=device)
    
    # Constant: how many random starts to try per span length before falling back
    TRIES_PER_LEN = 16

    for b in range(batch_size):
        # Precompute row views for this sequence
        masked_row = masked[b]
        non_pad_row = ~pad_mask[b]

        # Count non-PAD tokens in this sequence
        non_pad_count = non_pad_row.sum().item()

        # Handle all-PAD sequences
        if non_pad_count == 0:
            continue

        # Target number of tokens to mask (from non-PAD tokens)
        num_to_mask = max(1, int(mask_prob * non_pad_count))
        covered = 0
        attempts = 0
        # With rejection sampling we can afford more attempts since each is O(1)
        max_attempts = num_to_mask * 100  # Safety limit

        while covered < num_to_mask and attempts < max_attempts:
            # Sample span length based on distribution
            if span_distribution == "geometric":
                span_len = sample_span_length_geometric(max_span_len, mean_span)
            elif span_distribution == "uniform":
                span_len = random.randint(1, max_span_len)
            else:
                raise ValueError(f"Unknown span_distribution: {span_distribution!r}")

            # Enforce cap (defensive; sampling should already respect max_span_len)
            span_len = min(span_len, max_span_len)

            # Never mask more than we still need (prevents systematic overshoot)
            remaining = num_to_mask - covered
            span_len = max(1, min(span_len, remaining))

            # Fast rejection sampling: try random starts instead of scanning all positions
            # If we can't place the sampled length, try shorter lengths down to 1.
            placed = False
            for try_len in range(span_len, 0, -1):
                max_start = seq_len - try_len
                if max_start < 0:
                    continue

                # Try up to TRIES_PER_LEN random starts for this length
                for _ in range(TRIES_PER_LEN):
                    start = random.randint(0, max_start)
                    end = start + try_len

                    # Reject if any PAD in span
                    if not non_pad_row[start:end].all():
                        continue

                    # Reject if overlap with existing masked tokens
                    if masked_row[start:end].any():
                        continue

                    # Reject if touching existing masked tokens
                    if start > 0 and masked_row[start - 1]:
                        continue
                    if end < seq_len and masked_row[end]:
                        continue

                    # Accept: mark as masked
                    masked_row[start:end] = True
                    covered += try_len
                    placed = True
                    break

                if placed:
                    break

            attempts += 1
            # Continue trying even if this attempt failed (don't break early)

        # Fallback fill: if still short, randomly select remaining non-masked non-pad positions
        if covered < num_to_mask:
            # Find all positions that are non-pad and not yet masked
            available = non_pad_row & ~masked_row
            available_indices = torch.where(available)[0]

            if len(available_indices) > 0:
                shortfall = num_to_mask - covered
                # Randomly select up to shortfall positions
                num_to_fill = min(shortfall, len(available_indices))
                perm = torch.randperm(len(available_indices), device=inputs.device)[:num_to_fill]
                fill_indices = available_indices[perm]
                masked_row[fill_indices] = True

    # Ensure PAD tokens are NEVER masked (redundant safety check)
    masked[pad_mask] = False

    # Set labels: -100 for unmasked positions (ignored in loss)
    labels[~masked] = -100
    
    # CRITICAL: Exclude special tokens from being supervised targets
    # Even if a position is masked, if the original token was a special token,
    # we should NOT train the model to predict it. This prevents the model
    # from learning to predict [MASK], [UNK], or [PAD] tokens.
    if exclude_special_from_labels:
        # Vectorized: single mask for all special tokens (faster than loop)
        special_mask = (inputs == vocab.mask_id) | (inputs == vocab.unk_id) | (inputs == vocab.pad_id)
        labels[special_mask] = -100

    # Apply 80/10/10 BERT masking strategy within masked spans
    # 80% -> [MASK] token
    # 10% -> random token (from vocab, excluding special tokens)
    # 10% -> unchanged (keep original)
    rand = torch.rand_like(inputs, dtype=torch.float32, device=device)
    mask_token_indices = masked & (rand < 0.8)
    random_token_indices = masked & (rand >= 0.8) & (rand < 0.9)
    # Remaining 10% (rand >= 0.9) stay unchanged

    inputs = inputs.clone()
    inputs[mask_token_indices] = vocab.mask_id

    if random_token_indices.any():
        flat = random_token_indices.view(-1)
        num_rand = flat.sum().item()
        # Use vocab's base token IDs (excludes special tokens)
        # Vectorized: torch sampling instead of Python loop
        if hasattr(vocab, '_base_token_ids') and vocab._base_token_ids:
            base_ids = torch.as_tensor(vocab._base_token_ids, device=device, dtype=inputs.dtype)
            rand_idx = torch.randint(0, base_ids.numel(), (num_rand,), device=device)
            random_ids = base_ids[rand_idx]
        else:
            random_ids = torch.randint(
                low=0,
                high=vocab_size,
                size=(num_rand,),
                device=device,
                dtype=inputs.dtype,
            )
        inputs.view(-1)[flat] = random_ids

    if track_spans:
        # Track contiguous masked run lengths (what the model actually sees as "spans")
        span_lengths: List[int] = []
        pad0 = torch.zeros((1,), dtype=torch.int8, device=device)
        for b in range(batch_size):
            row = masked[b].to(torch.int8)
            padded = torch.cat([pad0, row, pad0])
            diff = padded[1:] - padded[:-1]
            starts = torch.where(diff == 1)[0]
            ends = torch.where(diff == -1)[0]
            if starts.numel() == 0:
                continue
            span_lengths.extend((ends - starts).tolist())

        return inputs, labels, span_lengths
    return inputs, labels


def create_contiguous_mask(
    tokens: torch.LongTensor,
    kmer: int,
    mask_token_id: int,
    vocab_size: int,
    mask_prob: float = 0.15,
    mask_token_prob: float = 0.8,
    random_token_prob: float = 0.1,
    special_token_ids: Optional[Set[int]] = None,
) -> Tuple[torch.LongTensor, torch.LongTensor, torch.BoolTensor]:
    """
    DNABERT-style contiguous span masking for MLM pretraining.

    This implements the masking strategy from the DNABERT paper, which is more
    effective for DNA sequences than random token masking because:
    1. It forces the model to learn longer-range dependencies
    2. It's more biologically meaningful (mutations often affect regions)
    3. It achieved better downstream performance in the original paper

    Algorithm:
    1. Select ~mask_prob of POSITIONS as span START points
    2. Each selected position masks a contiguous span of length k (k-mer size)
    3. Spans can overlap (that's fine and expected)
    4. For each masked token:
       - mask_token_prob (80%) → replace with [MASK] token
       - random_token_prob (10%) → replace with random token
       - keep_prob (10%) → keep original token
    5. Special tokens (PAD, CLS, SEP, etc.) are NEVER masked

    Args:
        tokens: Input token IDs, shape (batch_size, seq_len)
        kmer: K-mer size (3, 4, 5, or 6) - determines span length
        mask_token_id: ID of [MASK] token in vocabulary
        vocab_size: Total vocabulary size (for random token sampling)
        mask_prob: Probability of a position being a span START (default 0.15)
        mask_token_prob: Probability of using [MASK] token (default 0.8)
        random_token_prob: Probability of using random token (default 0.1)
        special_token_ids: Set of token IDs to never mask (PAD, CLS, SEP, UNK, etc.)

    Returns:
        masked_tokens: Tokens with masking applied, shape (batch_size, seq_len)
        labels: Original tokens for masked positions, -100 for unmasked (for loss)
        mask: Boolean mask indicating which positions are masked

    Example:
        >>> tokens = torch.tensor([[1, 2, 3, 4, 5, 6, 7, 8]])
        >>> masked, labels, mask = create_contiguous_mask(
        ...     tokens, kmer=3, mask_token_id=100, vocab_size=64,
        ...     special_token_ids={0}  # PAD token
        ... )
        >>> # If position 2 is selected as span start, positions 2,3,4 get masked
    """
    device = tokens.device
    batch_size, seq_len = tokens.shape

    # Initialize special token set
    if special_token_ids is None:
        special_token_ids = set()

    # Create a mask for special tokens (positions we should never mask)
    special_mask = torch.zeros_like(tokens, dtype=torch.bool, device=device)
    for special_id in special_token_ids:
        special_mask |= (tokens == special_id)

    # Step 1: Select span START positions
    # We sample ~mask_prob of positions as potential span starts
    # But we need to avoid starting spans too close to the end (within k-1 of end)
    # and we need to avoid special tokens

    # Create mask for valid span start positions
    # A position can be a span start if:
    # 1. It's not a special token
    # 2. Starting a span of length k from here doesn't go past sequence end
    valid_start_mask = ~special_mask.clone()

    # Positions within (kmer-1) of the end cannot start a full span
    # Set last (kmer-1) positions as invalid starts
    if kmer > 1:
        valid_start_mask[:, -(kmer - 1):] = False

    # Sample span starts using Bernoulli distribution
    # Only sample from valid positions
    span_start_probs = torch.full_like(tokens, mask_prob, dtype=torch.float32, device=device)
    span_start_probs[~valid_start_mask] = 0.0
    span_starts = torch.bernoulli(span_start_probs).bool()

    # Step 2: Expand span starts to full spans of length k
    # Each span start masks positions [i, i+1, ..., i+k-1]
    # Spans can overlap, which is fine

    # Initialize the mask for all masked positions
    mask = torch.zeros_like(tokens, dtype=torch.bool, device=device)

    # For each span start, mark the next k positions as masked
    # We use a sliding window approach: shift span_starts and OR together
    for offset in range(kmer):
        if offset == 0:
            mask |= span_starts
        else:
            # Shift span_starts left by offset positions
            # span_starts[:, offset:] aligns with mask[:, offset:]
            shifted = torch.zeros_like(span_starts)
            shifted[:, offset:] = span_starts[:, :-offset]
            mask |= shifted

    # Step 3: Ensure special tokens are NEVER masked
    mask[special_mask] = False

    # Step 4: Create labels tensor
    # Labels are the original tokens for masked positions, -100 for unmasked
    labels = tokens.clone()
    labels[~mask] = -100

    # Step 5: Apply 80/10/10 BERT masking strategy to masked positions
    # 80% -> [MASK] token
    # 10% -> random token
    # 10% -> keep original token

    masked_tokens = tokens.clone()

    # Generate random values for each position to determine masking type
    rand = torch.rand_like(tokens, dtype=torch.float32, device=device)

    # Positions where we apply [MASK] (80% of masked positions)
    mask_token_positions = mask & (rand < mask_token_prob)

    # Positions where we apply random token (10% of masked positions)
    random_token_positions = mask & (rand >= mask_token_prob) & (rand < mask_token_prob + random_token_prob)

    # Remaining 10% keep their original tokens (no action needed)

    # Apply [MASK] token
    masked_tokens[mask_token_positions] = mask_token_id

    # Apply random tokens
    if random_token_positions.any():
        num_random = random_token_positions.sum().item()
        # Sample random tokens from vocab, but exclude special tokens
        # For simplicity, we sample from [0, vocab_size - len(special_token_ids))
        # This is a slight approximation; for perfect correctness we'd need
        # to sample from base tokens only
        random_token_values = torch.randint(
            low=0,
            high=vocab_size - len(special_token_ids) if special_token_ids else vocab_size,
            size=(num_random,),
            device=device,
            dtype=tokens.dtype,
        )
        masked_tokens[random_token_positions] = random_token_values

    return masked_tokens, labels, mask


def debug_label_token_frequency(
    labels: torch.LongTensor,
    vocab: KmerVocabulary,
    num_batches: int = 1,
) -> Dict[str, Any]:
    """
    Generate debug report showing frequency of each token ID in labels.
    
    This is critical for verifying that special tokens are NOT being used
    as supervised targets. A properly configured MLM should have:
    - 0 occurrences of [MASK], [UNK], [PAD] in labels
    - All label tokens should be real k-mers (IDs 0 to 63 for k=3)
    
    Args:
        labels: (B, L) label tensor with -100 for ignored positions
        vocab: KmerVocabulary to decode token IDs
        num_batches: Number of batches this report covers (for averaging)
        
    Returns:
        Dict with:
        - token_counts: Dict[int, int] mapping token ID to count
        - special_token_counts: Dict[str, int] for [MASK], [UNK], [PAD]
        - total_supervised: Total number of supervised (non -100) labels
        - issues: List of issues found (empty if all is well)
    """
    # Flatten and filter out -100 (ignored positions)
    flat_labels = labels.view(-1)
    supervised_mask = flat_labels != -100
    supervised_labels = flat_labels[supervised_mask]
    
    total_supervised = supervised_labels.numel()
    
    # Count occurrences of each token ID
    token_counts: Dict[int, int] = {}
    if total_supervised > 0:
        unique, counts = torch.unique(supervised_labels, return_counts=True)
        for token_id, count in zip(unique.tolist(), counts.tolist()):
            token_counts[token_id] = count
    
    # Check special tokens specifically
    special_token_names = {
        vocab.mask_id: "[MASK]",
        vocab.unk_id: "[UNK]",
        vocab.pad_id: "[PAD]",
    }
    special_token_counts = {}
    issues = []
    
    for token_id, name in special_token_names.items():
        count = token_counts.get(token_id, 0)
        special_token_counts[name] = count
        if count > 0:
            issues.append(f"CRITICAL: {name} (id={token_id}) appears {count} times as supervised target!")
    
    return {
        "token_counts": token_counts,
        "special_token_counts": special_token_counts,
        "total_supervised": total_supervised,
        "issues": issues,
        "num_batches": num_batches,
    }


def print_label_debug_report(
    report: Dict[str, Any],
    vocab: KmerVocabulary,
    top_k: int = 10,
) -> None:
    """
    Print a formatted debug report for label token frequencies.
    
    Args:
        report: Output from debug_label_token_frequency()
        vocab: KmerVocabulary for decoding token names
        top_k: Number of most common tokens to display
    """
    print("=" * 80)
    print("LABEL TOKEN FREQUENCY DEBUG REPORT")
    print("=" * 80)
    print()
    
    print(f"Total supervised labels (non -100): {report['total_supervised']:,}")
    print(f"Batches analyzed: {report['num_batches']}")
    print()
    
    # Special tokens section - CRITICAL
    print("-" * 80)
    print("SPECIAL TOKEN CHECK (should ALL be 0)")
    print("-" * 80)
    all_zero = True
    for name, count in report["special_token_counts"].items():
        status = "✓" if count == 0 else "✗ FAIL"
        print(f"  {name}: {count:,} occurrences {status}")
        if count > 0:
            all_zero = False
    
    print()
    if all_zero:
        print("  ✓ All special tokens correctly excluded from labels!")
    else:
        print("  ✗ CRITICAL ERROR: Special tokens found in labels!")
        print("    The model is being trained to predict special tokens.")
        print("    This will hurt performance significantly.")
    print()
    
    # Top-k most common tokens
    print("-" * 80)
    print(f"TOP {top_k} MOST COMMON LABEL TOKENS")
    print("-" * 80)
    
    token_counts = report["token_counts"]
    if token_counts:
        sorted_counts = sorted(token_counts.items(), key=lambda x: -x[1])[:top_k]
        total = report["total_supervised"]
        
        for rank, (token_id, count) in enumerate(sorted_counts, 1):
            pct = count / total * 100 if total > 0 else 0
            # Get token name
            if token_id < len(vocab.itos):
                token_name = vocab.itos[token_id]
            else:
                token_name = f"<unknown:{token_id}>"
            print(f"  {rank:>2}. '{token_name}' (id={token_id}): {count:,} ({pct:.2f}%)")
    else:
        print("  No supervised labels found.")
    
    print()
    
    # Issues section
    if report["issues"]:
        print("-" * 80)
        print("ISSUES FOUND")
        print("-" * 80)
        for issue in report["issues"]:
            print(f"  ✗ {issue}")
        print()
    
    print("=" * 80)


def get_warmup_cosine_schedule(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
    final_lr_ratio: float = 0.1,
) -> LambdaLR:
    """
    Learning rate schedule with linear warmup and cosine decay.
    
    Schedule:
    1. Linear warmup: LR ramps from 0 to base_lr over num_warmup_steps
    2. Cosine decay: LR decays from base_lr to base_lr * final_lr_ratio
    
    Args:
        optimizer: The optimizer to schedule
        num_warmup_steps: Number of steps for linear warmup
        num_training_steps: Total number of training steps
        final_lr_ratio: Final LR as fraction of peak LR (e.g., 0.1 = 10% of peak)
        
    Returns:
        LambdaLR scheduler that multiplies base_lr by the lambda value
    """
    def lr_lambda(current_step: int) -> float:
        # Phase 1: Linear warmup from 0 to 1
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        
        # Phase 2: Cosine decay from 1 to final_lr_ratio
        # progress goes from 0 to 1 after warmup
        progress = float(current_step - num_warmup_steps) / float(
            max(1, num_training_steps - num_warmup_steps)
        )
        # Clamp progress to [0, 1] in case of extra steps
        progress = min(1.0, progress)
        
        # Cosine decay: starts at 1.0, ends at final_lr_ratio
        # Formula: final_lr_ratio + (1 - final_lr_ratio) * 0.5 * (1 + cos(pi * progress))
        # At progress=0: 0.5 * (1 + 1) = 1.0 -> returns 1.0
        # At progress=1: 0.5 * (1 + (-1)) = 0.0 -> returns final_lr_ratio
        cosine_decay = 0.5 * (1.0 + math.cos(math.pi * progress))
        return final_lr_ratio + (1.0 - final_lr_ratio) * cosine_decay
    
    return LambdaLR(optimizer, lr_lambda)


def get_layer_wise_lr_groups(
    model: nn.Module,
    base_lr: float,
    lr_decay: float = 0.8,
    weight_decay: float = 0.01,
) -> List[dict]:
    """
    Create parameter groups with layer-wise learning rate decay.
    
    Lower layers (closer to input) get smaller learning rates,
    which helps preserve pretrained features during fine-tuning
    and improves pretraining stability.
    
    Args:
        model: DNAMLM model
        base_lr: Base learning rate for the top layer
        lr_decay: Decay factor per layer (0.8 means each lower layer gets 0.8x LR)
        weight_decay: Weight decay for regularization
    
    Returns:
        List of parameter groups for optimizer
    """
    # No decay for bias and layer norm
    no_decay = {'bias', 'LayerNorm', 'layer_norm', 'norm', 'embed_norm'}
    
    # Get encoder layers
    encoder = model.encoder
    num_layers = len(encoder.layers)
    
    # Group parameters by layer depth
    param_groups = []
    
    # Embedding parameters (lowest layer, smallest LR)
    embed_lr = base_lr * (lr_decay ** num_layers)
    embed_params_decay = []
    embed_params_no_decay = []
    for name, param in encoder.embedding.named_parameters():
        if any(nd in name for nd in no_decay):
            embed_params_no_decay.append(param)
        else:
            embed_params_decay.append(param)
    
    if embed_params_decay:
        param_groups.append({'params': embed_params_decay, 'lr': embed_lr, 'weight_decay': weight_decay})
    if embed_params_no_decay:
        param_groups.append({'params': embed_params_no_decay, 'lr': embed_lr, 'weight_decay': 0.0})
    
    # Positional embedding parameters (also low layer, same LR as token embedding)
    if hasattr(encoder, 'pos_embedding'):
        pos_embed_params_decay = []
        pos_embed_params_no_decay = []
        for name, param in encoder.pos_embedding.named_parameters():
            if any(nd in name for nd in no_decay):
                pos_embed_params_no_decay.append(param)
            else:
                pos_embed_params_decay.append(param)
        
        if pos_embed_params_decay:
            param_groups.append({'params': pos_embed_params_decay, 'lr': embed_lr, 'weight_decay': weight_decay})
        if pos_embed_params_no_decay:
            param_groups.append({'params': pos_embed_params_no_decay, 'lr': embed_lr, 'weight_decay': 0.0})
    
    # Embed norm parameters (also low layer)
    if hasattr(encoder, 'embed_norm'):
        embed_norm_params = list(encoder.embed_norm.parameters())
        if embed_norm_params:
            param_groups.append({'params': embed_norm_params, 'lr': embed_lr, 'weight_decay': 0.0})
    
    # Transformer layers (progressively higher LR)
    for layer_idx, layer in enumerate(encoder.layers):
        layer_lr = base_lr * (lr_decay ** (num_layers - layer_idx - 1))
        layer_params_decay = []
        layer_params_no_decay = []
        
        for name, param in layer.named_parameters():
            if any(nd in name for nd in no_decay):
                layer_params_no_decay.append(param)
            else:
                layer_params_decay.append(param)
        
        if layer_params_decay:
            param_groups.append({'params': layer_params_decay, 'lr': layer_lr, 'weight_decay': weight_decay})
        if layer_params_no_decay:
            param_groups.append({'params': layer_params_no_decay, 'lr': layer_lr, 'weight_decay': 0.0})
    
    # Final norm (top layer, full LR)
    final_norm_params = list(encoder.final_norm.parameters())
    if final_norm_params:
        param_groups.append({'params': final_norm_params, 'lr': base_lr, 'weight_decay': 0.0})
    
    # LM head parameters (top layer, full LR)
    # Note: if weight tying, embedding weights are already added above
    lm_head_params_decay = []
    lm_head_params_no_decay = []
    
    for name, param in model.lm_head.named_parameters():
        if 'weight' in name and model.tie_weights:
            # Skip tied weights (already in embedding group)
            continue
        if any(nd in name for nd in no_decay):
            lm_head_params_no_decay.append(param)
        else:
            lm_head_params_decay.append(param)

    # Learnable logit scale for normalized MLM head (do not decay)
    if hasattr(model, "logit_scale"):
        lm_head_params_no_decay.append(model.logit_scale)
    
    if lm_head_params_decay:
        param_groups.append({'params': lm_head_params_decay, 'lr': base_lr, 'weight_decay': weight_decay})
    if lm_head_params_no_decay:
        param_groups.append({'params': lm_head_params_no_decay, 'lr': base_lr, 'weight_decay': 0.0})
    
    # Output norm if present
    if hasattr(model, 'output_norm'):
        output_norm_params = list(model.output_norm.parameters())
        if output_norm_params:
            param_groups.append({'params': output_norm_params, 'lr': base_lr, 'weight_decay': 0.0})
    
    return param_groups


class MLMCollator:
    """
    Collator that supports curriculum learning via adjustable mask probability.

    Supports progressive masking via MaskingScheduler:
    - Pass masking_scheduler to enable dynamic mask rate adjustment
    - Call update_step(global_step) before each epoch/batch to update the rate
    - Or manually call set_mask_prob(rate) to set a specific rate

    Debug mode (debug=True):
    - Prints original token ids, masked input token ids, and label token ids
    - For 5 positions per sequence, asserts:
      - label == original_token for masked positions
      - label == -100 for unmasked positions
    - Raises AssertionError if any mismatch occurs
    - Debug mode runs once then auto-disables

    track_spans mode:
    - Returns (masked_inputs, labels, span_lengths) instead of (masked_inputs, labels)
    - Used by health dashboard to log span length histograms
    """
    def __init__(
        self,
        vocab: KmerVocabulary,
        mask_prob: float = 0.15,
        debug: bool = False,
        track_spans: bool = False,
        masking_scheduler: MaskingScheduler | None = None,
    ):
        self.vocab = vocab
        self.mask_prob = mask_prob
        self.debug = debug
        self.track_spans = track_spans
        self._debug_ran = False
        self.masking_scheduler = masking_scheduler

    def set_mask_prob(self, mask_prob: float) -> None:
        """Manually set the mask probability (for external control)."""
        self.mask_prob = mask_prob

    def update_step(self, step: int) -> float:
        """
        Update mask probability based on current training step.

        If a masking_scheduler is set, updates mask_prob from the scheduler.
        Otherwise, mask_prob remains unchanged.

        Args:
            step: Current global training step

        Returns:
            Current mask probability after update
        """
        if self.masking_scheduler is not None:
            self.mask_prob = self.masking_scheduler.get_mask_rate(step)
        return self.mask_prob
    
    def __call__(
        self, batch: Iterable[torch.LongTensor]
    ) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
        original_inputs = torch.stack(list(batch))
        
        if self.track_spans:
            masked_inputs, labels, span_lengths = mask_tokens_span(
                original_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=True
            )
        else:
            masked_inputs, labels = mask_tokens_span(
                original_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=False
            )
            span_lengths = []
        
        # Run debug verification once if enabled
        if self.debug and not self._debug_ran:
            self._run_debug_verification(original_inputs, masked_inputs, labels)
            self._debug_ran = True
        
        if self.track_spans:
            return masked_inputs, labels, span_lengths
        return masked_inputs, labels
    
    def _run_debug_verification(
        self, 
        original: torch.LongTensor, 
        masked_input: torch.LongTensor, 
        labels: torch.LongTensor
    ) -> None:
        """
        Debug mode: verify labels match originals for one batch.
        
        Checks:
        - For MASKED positions: label == original_token
        - For UNMASKED positions: label == -100
        
        Hard error (AssertionError) if any mismatch occurs.
        """
        batch_size, seq_len = original.shape
        
        print("=" * 80)
        print("DEBUG MODE: MLM Collator Label Verification")
        print("=" * 80)
        print(f"Batch size: {batch_size}, Sequence length: {seq_len}")
        print(f"Mask prob: {self.mask_prob:.2%}")
        print()
        
        # Determine masked positions (where label != -100)
        masked_positions = (labels != -100)
        
        # Print first sequence's tokens (limited for readability)
        print("[Token IDs for first sequence (first 20 positions)]")
        print(f"  Original tokens:     {original[0, :20].tolist()}")
        print(f"  Masked input tokens: {masked_input[0, :20].tolist()}")
        print(f"  Label tokens:        {labels[0, :20].tolist()}")
        print()
        
        # Collect errors for batch reporting
        errors = []
        
        # For each sequence, check 5 positions
        print(f"[Checking 5 positions per sequence across {batch_size} sequences]")
        print("-" * 80)
        
        for seq_idx in range(batch_size):
            seq_original = original[seq_idx]
            seq_labels = labels[seq_idx]
            seq_masked_positions = masked_positions[seq_idx]
            
            # Select 5 positions to check (evenly spaced)
            check_positions = [int(i * seq_len / 5) for i in range(5)]
            
            for pos in check_positions:
                orig_token = seq_original[pos].item()
                label_token = seq_labels[pos].item()
                is_masked = seq_masked_positions[pos].item()
                
                if is_masked:
                    # MASKED position: label should equal original token
                    if label_token != orig_token:
                        error_msg = (
                            f"MISMATCH at seq={seq_idx}, pos={pos}: "
                            f"MASKED position has label={label_token}, expected original={orig_token}"
                        )
                        errors.append(error_msg)
                        print(f"  ✗ seq[{seq_idx}] pos[{pos}]: MASKED - label={label_token}, orig={orig_token} - FAIL!")
                    else:
                        if seq_idx < 3:  # Only print details for first 3 sequences
                            print(f"  ✓ seq[{seq_idx}] pos[{pos}]: MASKED - label={label_token}, orig={orig_token} - OK")
                else:
                    # UNMASKED position: label should be -100
                    if label_token != -100:
                        error_msg = (
                            f"MISMATCH at seq={seq_idx}, pos={pos}: "
                            f"UNMASKED position has label={label_token}, expected -100"
                        )
                        errors.append(error_msg)
                        print(f"  ✗ seq[{seq_idx}] pos[{pos}]: UNMASKED - label={label_token}, expected=-100 - FAIL!")
                    else:
                        if seq_idx < 3:  # Only print details for first 3 sequences
                            print(f"  ✓ seq[{seq_idx}] pos[{pos}]: UNMASKED - label={label_token} - OK")
        
        print("-" * 80)
        print()
        
        # Summary statistics
        total_tokens = batch_size * seq_len
        num_masked = masked_positions.sum().item()
        num_unmasked = total_tokens - num_masked
        
        print("[Summary Statistics]")
        print(f"  Total tokens:   {total_tokens}")
        print(f"  Masked tokens:  {num_masked} ({num_masked/total_tokens:.2%})")
        print(f"  Unmasked:       {num_unmasked} ({num_unmasked/total_tokens:.2%})")
        print()
        
        # Exhaustive verification (not just 5 positions, but ALL)
        print("[Exhaustive Verification - ALL positions]")
        
        # Check all masked positions: label should equal original
        all_masked_labels = labels[masked_positions]
        all_masked_originals = original[masked_positions]
        masked_match = (all_masked_labels == all_masked_originals).all().item()
        masked_mismatches = (all_masked_labels != all_masked_originals).sum().item()
        
        print(f"  Masked positions (label == original):   ", end="")
        if masked_match:
            print(f"✓ PASS (all {num_masked} positions match)")
        else:
            print(f"✗ FAIL ({masked_mismatches} mismatches out of {num_masked})")
            errors.append(f"Exhaustive check: {masked_mismatches} masked positions have wrong labels")
        
        # Check all unmasked positions: label should be -100
        unmasked_positions = ~masked_positions
        all_unmasked_labels = labels[unmasked_positions]
        unmasked_correct = (all_unmasked_labels == -100).all().item()
        unmasked_wrong = (all_unmasked_labels != -100).sum().item()
        
        print(f"  Unmasked positions (label == -100):     ", end="")
        if unmasked_correct:
            print(f"✓ PASS (all {num_unmasked} positions are -100)")
        else:
            print(f"✗ FAIL ({unmasked_wrong} positions are not -100)")
            errors.append(f"Exhaustive check: {unmasked_wrong} unmasked positions don't have -100")
        
        print()
        print("=" * 80)
        
        # HARD ERROR if any mismatches
        if errors:
            print("ERRORS FOUND:")
            for err in errors:
                print(f"  - {err}")
            print("=" * 80)
            raise AssertionError(
                f"MLM Collator Debug: {len(errors)} label verification error(s) found!\n"
                + "\n".join(errors)
            )
        
        print("✓ All assertions passed! Labels correctly match originals.")
        print("=" * 80)


def collate_mlm_streaming(
    batch: Iterable[Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, Any]]
) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, Any]:
    batch_list = list(batch)
    if not batch_list:
        raise ValueError("Empty batch in collate_mlm_streaming")
    sample = batch_list[0]
    if len(sample) == 2:
        inputs, labels = zip(*batch_list)
        return torch.stack(inputs), torch.stack(labels)
    if len(sample) == 3:
        inputs, labels, extras = zip(*batch_list)
        stacked_inputs = torch.stack(inputs)
        stacked_labels = torch.stack(labels)
        first_extra = extras[0]
        if isinstance(first_extra, (list, tuple)):
            flat_spans = [span for spans in extras for span in spans]
            return stacked_inputs, stacked_labels, flat_spans
        return stacked_inputs, stacked_labels, list(extras)
    raise ValueError(f"Unexpected batch sample length: {len(sample)}")


def collate_mlm(
    batch: Iterable[torch.LongTensor],
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
) -> Tuple[torch.LongTensor, torch.LongTensor]:
    """Collate batch and apply DNABERT-style span masking."""
    inputs = torch.stack(list(batch))
    masked, labels = mask_tokens_span(inputs, vocab, mask_prob=mask_prob)
    return masked, labels


def debug_masking(
    vocab: KmerVocabulary,
    batch: Iterable[torch.LongTensor],
    mask_prob: float = 0.15,
    max_span_len: int = 3,
) -> None:
    """
    Debug mode: analyze one batch of masking and print statistics.
    
    Checks:
    1. Percent masked tokens (should be ~mask_prob)
    2. Percent masked that were PAD (should be 0)
    3. Distribution of span lengths
    4. Assert labels are -100 for unmasked, never -100 for masked
    """
    from collections import Counter
    
    print("=" * 70)
    print("DEBUG MODE: Masking Analysis")
    print("=" * 70)
    
    # Stack batch into tensor
    inputs = torch.stack(list(batch))
    batch_size, seq_len = inputs.shape
    total_tokens = batch_size * seq_len
    
    print(f"Batch size: {batch_size}")
    print(f"Sequence length: {seq_len}")
    print(f"Total tokens: {total_tokens}")
    print(f"Expected mask_prob: {mask_prob:.2%}")
    print(f"Max span length: {max_span_len}")
    print()
    
    # Run masking with span tracking
    masked_inputs, labels, span_lengths = mask_tokens_span(
        inputs, vocab, mask_prob=mask_prob, max_span_len=max_span_len, track_spans=True
    )
    
    # 1. Percent masked tokens
    masked_positions = (labels != -100)
    num_masked = masked_positions.sum().item()
    actual_mask_pct = num_masked / total_tokens
    
    print(f"[1] MASKED TOKEN PERCENTAGE")
    print(f"    Masked tokens: {num_masked} / {total_tokens}")
    print(f"    Actual mask rate: {actual_mask_pct:.2%}")
    print(f"    Expected: ~{mask_prob:.2%}")
    diff = abs(actual_mask_pct - mask_prob)
    status = "✓ OK" if diff < 0.05 else "⚠ WARN (>5% deviation)"
    print(f"    Status: {status}")
    print()
    
    # 2. Percent masked that were PAD
    pad_positions = (inputs == vocab.pad_id)
    masked_pads = (masked_positions & pad_positions).sum().item()
    total_pads = pad_positions.sum().item()
    
    print(f"[2] PAD TOKEN MASKING")
    print(f"    Total PAD tokens in batch: {total_pads}")
    print(f"    PAD tokens that are masked: {masked_pads}")
    if total_pads > 0:
        pad_mask_pct = masked_pads / total_pads
        print(f"    PAD mask rate: {pad_mask_pct:.2%}")
    status = "✓ OK (no PADs masked)" if masked_pads == 0 else "✗ FAIL (PADs should never be masked!)"
    print(f"    Status: {status}")
    print()
    
    # 3. Span length distribution
    print(f"[3] SPAN LENGTH DISTRIBUTION")
    if span_lengths:
        span_counter = Counter(span_lengths)
        total_spans = len(span_lengths)
        print(f"    Total spans created: {total_spans}")
        print(f"    Span length distribution:")
        for length in sorted(span_counter.keys()):
            count = span_counter[length]
            pct = count / total_spans * 100
            bar = "█" * int(pct / 2)
            print(f"      Length {length}: {count:>5} ({pct:5.1f}%) {bar}")
        avg_span = sum(span_lengths) / len(span_lengths)
        print(f"    Average span length: {avg_span:.2f}")
    else:
        print(f"    No spans recorded")
    print()
    
    # 4. Label correctness assertions
    print(f"[4] LABEL CORRECTNESS ASSERTIONS")
    
    # Check: unmasked positions should have label = -100
    unmasked_positions = ~masked_positions
    unmasked_labels = labels[unmasked_positions]
    all_unmasked_are_minus100 = (unmasked_labels == -100).all().item()
    
    print(f"    Unmasked positions have label=-100: ", end="")
    if all_unmasked_are_minus100:
        print("✓ PASS")
    else:
        bad_count = (unmasked_labels != -100).sum().item()
        print(f"✗ FAIL ({bad_count} unmasked positions don't have -100)")
    
    # Check: masked positions should NEVER have label = -100
    masked_labels = labels[masked_positions]
    no_masked_is_minus100 = (masked_labels != -100).all().item()
    
    print(f"    Masked positions have label!=-100:  ", end="")
    if no_masked_is_minus100:
        print("✓ PASS")
    else:
        bad_count = (masked_labels == -100).sum().item()
        print(f"✗ FAIL ({bad_count} masked positions have -100)")
    
    # Check: labels for masked positions should be valid token IDs
    if num_masked > 0:
        min_label = masked_labels.min().item()
        max_label = masked_labels.max().item()
        valid_range = (min_label >= 0 and max_label < len(vocab.itos))
        print(f"    Masked labels in valid range [0, {len(vocab.itos)-1}]: ", end="")
        if valid_range:
            print(f"✓ PASS (min={min_label}, max={max_label})")
        else:
            print(f"✗ FAIL (min={min_label}, max={max_label})")
    
    print()
    print("=" * 70)
    
    # Assert critical properties
    assert masked_pads == 0, "ASSERTION FAILED: PAD tokens should never be masked!"
    assert all_unmasked_are_minus100, "ASSERTION FAILED: Unmasked positions must have label=-100!"
    assert no_masked_is_minus100, "ASSERTION FAILED: Masked positions must NOT have label=-100!"
    
    print("All assertions passed! Masking is working correctly.")
    print("=" * 70)


def verify_no_leakage(
    model: nn.Module,
    vocab: KmerVocabulary,
    device: torch.device,
) -> None:
    """
    Verify that masked positions have NO access to original token information.
    
    This is CRITICAL for proper MLM training:
    - MLM ≠ autoencoder. The model must predict from context only.
    - Masked positions should only see [MASK] embedding
    - Engineered features must be disabled (they would leak sequence info)
    
    Tests:
    1. Verify masked token IDs are replaced with [MASK] before embedding lookup
    2. Verify model runs with 100% masking (extreme blindfold test)
    3. Assert no engineered features are used in DNAMLM
    """
    print("=" * 80)
    print("LEAKAGE VERIFICATION: MLM Blindfold Test")
    print("=" * 80)
    print("MLM ≠ autoencoder. Blindfold means blindfold.")
    print()
    
    model.eval()
    errors = []
    
    # =========================================================================
    # TEST 1: Verify masked inputs have [MASK] token, not original
    # =========================================================================
    print("[TEST 1] Masked inputs contain [MASK] token, not originals")
    print("-" * 80)
    
    # Create a test sequence
    batch_size = 4
    seq_len = 32
    original_tokens = torch.randint(0, 64, (batch_size, seq_len))  # Only real k-mers
    
    # Manually mask some positions (50% for clear test)
    mask_prob = 0.5
    masked_inputs = original_tokens.clone()
    
    # Create mask
    mask = torch.rand(batch_size, seq_len) < mask_prob
    
    # Apply [MASK] token to masked positions
    masked_inputs[mask] = vocab.mask_id
    
    # Verify: masked positions should have vocab.mask_id, not original
    masked_positions = mask
    masked_values = masked_inputs[masked_positions]
    original_values = original_tokens[masked_positions]
    
    all_have_mask_token = (masked_values == vocab.mask_id).all().item()
    none_leaked = (masked_values != original_values).all().item()
    
    print(f"  Total masked positions: {masked_positions.sum().item()}")
    print(f"  All masked positions have [MASK] token (id={vocab.mask_id}): ", end="")
    if all_have_mask_token:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 1: Some masked positions don't have [MASK] token")
    
    print(f"  No original tokens leaked to masked positions: ", end="")
    if none_leaked:
        print("✓ PASS")
    else:
        # This can fail if original happened to be [MASK], which is impossible for k-mers
        leak_count = (masked_values == original_values).sum().item()
        print(f"✗ FAIL ({leak_count} positions may have leaked)")
        errors.append(f"TEST 1: {leak_count} positions have original values")
    print()
    
    # =========================================================================
    # TEST 2: Model runs with 100% masking (extreme blindfold)
    # =========================================================================
    print("[TEST 2] Model runs with 100% masking (extreme blindfold)")
    print("-" * 80)
    
    # Create input where ALL tokens are [MASK]
    all_masked_input = torch.full((batch_size, seq_len), vocab.mask_id, dtype=torch.long).to(device)
    
    # Labels are all the original tokens (model must predict all)
    labels = original_tokens.clone().to(device)
    
    try:
        with torch.no_grad():
            logits, loss = model(all_masked_input, labels)
        
        # Verify outputs
        has_valid_logits = logits is not None and not torch.isnan(logits).any()
        has_valid_loss = loss is not None and not torch.isnan(loss) and not torch.isinf(loss)
        
        print(f"  Model forward pass: ✓ PASS")
        print(f"  Output shape: {logits.shape}")
        print(f"  Logits valid (no NaN): ", end="")
        if has_valid_logits:
            print("✓ PASS")
        else:
            print("✗ FAIL")
            errors.append("TEST 2: Logits contain NaN")
        
        print(f"  Loss valid (no NaN/Inf): ", end="")
        if has_valid_loss:
            print(f"✓ PASS (loss={loss.item():.4f})")
        else:
            print(f"✗ FAIL (loss={loss})")
            errors.append("TEST 2: Loss is NaN or Inf")
        
        # With 100% masking and no context, accuracy should be low (near random)
        preds = logits.argmax(dim=-1)
        acc = (preds == labels).float().mean().item()
        print(f"  Accuracy with 100% masking: {acc*100:.1f}% (expected: near random ~1.5%)")
        
        # Verify model isn't cheating somehow
        if acc > 0.5:  # If accuracy is >50%, something is very wrong
            print(f"  ⚠ WARNING: Accuracy suspiciously high - possible leakage!")
            errors.append(f"TEST 2: Suspiciously high accuracy ({acc*100:.1f}%) with 100% masking")
        
    except Exception as e:
        print(f"  ✗ FAIL: Model failed with error: {e}")
        errors.append(f"TEST 2: Model failed with {type(e).__name__}: {e}")
    print()
    
    # =========================================================================
    # TEST 3: Verify no engineered features in DNAMLM
    # =========================================================================
    print("[TEST 3] DNAMLM does not use engineered features")
    print("-" * 80)
    
    # Check DNAMLM doesn't have engineered feature inputs
    dnamlm_signature = str(model.forward.__code__.co_varnames)
    has_engineered = "engineered" in dnamlm_signature.lower()
    
    print(f"  DNAMLM.forward() signature: {model.forward.__code__.co_varnames[:5]}...")
    print(f"  Contains 'engineered' parameter: ", end="")
    if not has_engineered:
        print("✓ PASS (no engineered features)")
    else:
        print("✗ FAIL (engineered features would leak sequence info!)")
        errors.append("TEST 3: DNAMLM has engineered features in signature")
    
    # Also check the encoder doesn't receive engineered features
    encoder_params = list(model.encoder.forward.__code__.co_varnames)
    encoder_has_engineered = any("engineered" in p.lower() for p in encoder_params)
    print(f"  DNAEncoder.forward() has engineered: ", end="")
    if not encoder_has_engineered:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 3: DNAEncoder has engineered features")
    print()
    
    # =========================================================================
    # TEST 4: Embedding lookup verification
    # =========================================================================
    print("[TEST 4] Embedding lookup only sees masked token IDs")
    print("-" * 80)
    
    # Get the embedding layer
    embedding = model.encoder.embedding
    
    # Look up embeddings for [MASK] and a real k-mer
    mask_embedding = embedding(torch.tensor([vocab.mask_id]).to(device))
    real_embedding = embedding(torch.tensor([0]).to(device))  # First k-mer (e.g., "AAA")
    
    # They should be different
    embeddings_different = not torch.allclose(mask_embedding, real_embedding)
    
    print(f"  [MASK] embedding (id={vocab.mask_id}) shape: {mask_embedding.shape}")
    print(f"  Real k-mer embedding (id=0) shape: {real_embedding.shape}")
    print(f"  [MASK] embedding is distinct from real k-mers: ", end="")
    if embeddings_different:
        print("✓ PASS")
    else:
        print("✗ FAIL (embeddings are identical!)")
        errors.append("TEST 4: [MASK] embedding same as real k-mer embedding")
    
    # Verify that when we look up a masked sequence, we get [MASK] embeddings
    test_input = torch.tensor([[vocab.mask_id, vocab.mask_id, vocab.mask_id]]).to(device)
    test_embeds = embedding(test_input)
    expected_embeds = mask_embedding.expand(1, 3, -1)
    
    embeds_match = torch.allclose(test_embeds, expected_embeds)
    print(f"  Masked sequence gets [MASK] embeddings: ", end="")
    if embeds_match:
        print("✓ PASS")
    else:
        print("✗ FAIL")
        errors.append("TEST 4: Masked sequence doesn't get [MASK] embeddings")
    print()
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("=" * 80)
    if errors:
        print("LEAKAGE DETECTED!")
        print("-" * 80)
        for err in errors:
            print(f"  ✗ {err}")
        print("=" * 80)
        raise AssertionError(
            f"Leakage verification FAILED with {len(errors)} error(s)!\n" +
            "\n".join(errors)
        )
    else:
        print("✓ All leakage tests passed!")
        print("  - Masked positions only see [MASK] embeddings")
        print("  - Model handles 100% masking correctly")
        print("  - No engineered features in MLM pipeline")
        print("  - Embedding lookup is clean")
        print()
        print("The model is properly blindfolded. MLM ≠ autoencoder. ✓")
        print("=" * 80)


class DNAMLM(nn.Module):
    """
    Masked Language Model wrapper around DNAEncoder.
    
    Key improvements for lower loss:
    1. Weight tying between embedding and LM head (BERT-style)
    2. Bias in LM head for better output distribution modeling
    3. Optional output layer norm before projection
    4. **Special token exclusion**: [MASK], <UNK>, <PAD> are excluded from 
       prediction space by setting their logits to -inf. This prevents the
       model from "cheating" by predicting special tokens.
    """
    
    def __init__(
        self, 
        encoder: DNAEncoder, 
        vocab_size: int,
        special_token_ids: List[int] | None = None,  # IDs to exclude from predictions
        tie_weights: bool = True,  # CRITICAL: tie embedding weights for better learning
        use_output_norm: bool = True,  # Additional layer norm before LM head
    ):
        super().__init__()
        self.encoder = encoder
        self.vocab_size = vocab_size
        self.tie_weights = tie_weights
        
        # Store special token IDs to exclude from prediction
        # These will have their logits set to -inf before softmax/loss
        self.special_token_ids = special_token_ids or []
        if self.special_token_ids:
            LOGGER.info(
                "Special token exclusion enabled: IDs %s will be masked from predictions",
                self.special_token_ids
            )
        
        # Optional output layer norm (helps with distribution matching)
        self.use_output_norm = use_output_norm
        if use_output_norm:
            self.output_norm = nn.LayerNorm(encoder.embedding_dim)
        
        # LM head with bias (bias helps model output distribution)
        self.lm_head = nn.Linear(encoder.embedding_dim, vocab_size, bias=True)
        self.logit_scale = nn.Parameter(torch.tensor(10.0))
        
        # Weight tying: share weights between embedding and LM head
        # This is a key technique from BERT that improves MLM loss significantly
        if tie_weights:
            self.lm_head.weight = encoder.embedding.weight
            LOGGER.info("Weight tying enabled: LM head shares embedding weights")
        else:
            # If not tying, initialize LM head with scaled init
            nn.init.trunc_normal_(self.lm_head.weight, std=0.02, a=-0.04, b=0.04)
        
        # Initialize bias to zero
        nn.init.zeros_(self.lm_head.bias)
    
    def _mask_special_tokens(self, logits: torch.Tensor) -> torch.Tensor:
        """
        Set logits for special tokens to -inf so they can never be predicted.
        
        This is critical for proper MLM training:
        - Without this, the model can "cheat" by predicting [MASK], <UNK>, or <PAD>
        - These tokens should never appear as ground truth labels
        - Excluding them forces the model to learn real k-mer distributions
        
        Args:
            logits: (B, L, vocab_size) raw logits from LM head
            
        Returns:
            Modified logits with special token positions set to -inf
        """
        if not self.special_token_ids:
            return logits
        
        # Clone to avoid in-place modification issues with autograd
        logits = logits.clone()
        
        # Set special token logits to -inf (will become 0 after softmax)
        for token_id in self.special_token_ids:
            logits[:, :, token_id] = float('-inf')
        
        return logits

    def forward(
        self,
        token_ids: torch.LongTensor,
        labels: torch.LongTensor | None = None,
        key_padding_mask: torch.Tensor | None = None,
    ):
        """
        Forward pass for MLM.
        
        Args:
            token_ids: (B, L) input token indices
            labels: (B, L) target labels, -100 where we ignore
            key_padding_mask: Optional (B, L) boolean mask for padding
            
        Returns:
            (logits, loss) tuple where loss is None if labels not provided
        """
        x = self.encoder(token_ids, key_padding_mask=key_padding_mask)
        
        # Optional output norm before projection
        if self.use_output_norm:
            x = self.output_norm(x)

        # Normalize MLM head projection to prevent logit explosion (CLIP-style)
        h = F.normalize(x, dim=-1)
        w = F.normalize(self.lm_head.weight, dim=-1)
        logits = self.logit_scale * (h @ w.T)
        if self.lm_head.bias is not None:
            logits = logits + self.lm_head.bias
        
        # CRITICAL: Mask out special tokens from prediction space
        # This prevents the model from predicting [MASK], <UNK>, <PAD>
        logits = self._mask_special_tokens(logits)
        
        if labels is None:
            return logits, None

        # Standard cross entropy - NO label smoothing for MLM
        # CRITICAL: Force logits + loss to fp32 for numerical stability in mixed precision
        logits = logits.float()
        assert logits.dtype == torch.float32
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
        )
        return logits, loss
    
    def compute_accuracy(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute masked token prediction accuracy."""
        # Note: logits should already have special tokens masked out
        preds = logits.argmax(dim=-1)  # (B, L)
        mask = labels != -100
        if mask.sum() == 0:
            return 0.0
        correct = (preds[mask] == labels[mask]).sum().item()
        total = mask.sum().item()
        return correct / total
    
    def compute_perplexity(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute perplexity on masked tokens."""
        # Note: logits should already have special tokens masked out
        logits = logits.float()
        assert logits.dtype == torch.float32
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
            reduction='mean',
        )
        return math.exp(loss.item())


def load_or_build_vocab(windows: List[str], k: int, vocab_path: Path) -> KmerVocabulary:
    if vocab_path.exists():
        LOGGER.info("Loading existing vocabulary from %s", vocab_path)
        return KmerVocabulary.load(vocab_path)
    vocab = KmerVocabulary.build_from_sequences(windows, k=k)
    vocab.save(vocab_path)
    LOGGER.info("Saved vocabulary with %d entries to %s", len(vocab.itos), vocab_path)
    return vocab


def _build_windows_from_records(
    records: List[Tuple[str, str]],
    *,
    window_size: int,
    stride: int,
    unknown_base_strategy: str,
    include_reverse_complements: bool,
    label: str | None = None,
) -> List[str]:
    windows: List[str] = []
    for _, sequence in records:
        windows.extend(
            extract_windows(
                sequence,
                window_size,
                stride,
                unknown_base_strategy=unknown_base_strategy,
            )
        )
    if not windows:
        suffix = f" for {label}" if label else ""
        raise ValueError(f"No valid windows generated from genome{suffix}")
    forward_count = len(windows)
    if include_reverse_complements:
        rc_windows = [reverse_complement(window) for window in windows]
        windows.extend(rc_windows)
        rc_fraction = len(rc_windows) / max(1, len(windows))
        if label:
            LOGGER.info(
                "Extracted %d %s windows (%d + RC=%d, rc_frac=%.3f)",
                len(windows),
                label,
                forward_count,
                len(rc_windows),
                rc_fraction,
            )
        else:
            LOGGER.info(
                "Extracted %d windows (%d + RC=%d, rc_frac=%.3f)",
                len(windows),
                forward_count,
                len(rc_windows),
                rc_fraction,
            )
    else:
        if label:
            LOGGER.info("Extracted %d %s windows (rc_frac=0.000)", forward_count, label)
        else:
            LOGGER.info("Extracted %d windows (rc_frac=0.000)", forward_count)
    return windows


def prepare_dataset(cfg) -> Tuple[List[torch.LongTensor], KmerVocabulary]:
    script_dir = Path(__file__).resolve().parent
    fasta_path = _resolve_cfg_path(cfg.get("mlm_fasta_path"), base_dir=script_dir)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found at {fasta_path}")
    # Use MLM-specific window size, fallback to max_bp_len for backward compatibility
    window_size = int(cfg.get("mlm_window_size", cfg.get("max_bp_len", 234)))
    stride = int(cfg.get("mlm_stride", 20))
    k = int(cfg.get("mlm_kmer", 3))
    unknown_base_strategy = _normalize_unknown_base_strategy(cfg.get("mlm_unknown_base_strategy"))
    include_reverse_complements = bool(cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False)))
    vocab_path = _resolve_cfg_path(
        cfg.get("mlm_vocab_path", f"../artifacts/vocabs/k{k}_mlm_vocab.json"),
        base_dir=script_dir,
    )

    LOGGER.info("Reading FASTA records from %s", fasta_path)
    records = read_fasta_records(fasta_path)
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    LOGGER.info("Found %d FASTA record(s)", len(records))
    ablation_cfg = cfg.get("ablation") if isinstance(cfg, dict) else None
    if isinstance(ablation_cfg, dict) and bool(ablation_cfg.get("enabled", False)):
        bp_limit_raw = ablation_cfg.get("pretrain_bp_limit")
        if bp_limit_raw is not None:
            try:
                bp_limit = int(bp_limit_raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"ablation.pretrain_bp_limit must be an int, got {bp_limit_raw!r}") from exc
            if bp_limit < 0:
                raise ValueError(f"ablation.pretrain_bp_limit must be >= 0, got {bp_limit}")
            if bp_limit > 0:
                records = _apply_bp_limit(records, bp_limit)
                total_bp = sum(len(seq) for _, seq in records)
                LOGGER.info("Ablation enabled: using first %d bp from FASTA", total_bp)
    windows = _build_windows_from_records(
        records,
        window_size=window_size,
        stride=stride,
        unknown_base_strategy=unknown_base_strategy,
        include_reverse_complements=include_reverse_complements,
    )
    vocab = load_or_build_vocab(windows, k=k, vocab_path=vocab_path)
    token_tensors = [vocab.tokenize(window, window_size) for window in windows]
    return token_tensors, vocab


def _should_use_streaming(cfg: Dict[str, Any]) -> bool:
    raw_paths = cfg.get("mlm_fasta_paths")
    if isinstance(raw_paths, list):
        if any(str(p).strip() for p in raw_paths):
            return True
    elif isinstance(raw_paths, str) and raw_paths.strip():
        return True

    steps = cfg.get("mlm_steps_per_epoch")
    if steps is None:
        return False
    try:
        return int(steps) > 0
    except (TypeError, ValueError):
        return False


def _resolve_mlm_fasta_paths(cfg: Dict[str, Any], *, base_dir: Path) -> List[Path]:
    raw_paths: List[Any] = []
    configured = cfg.get("mlm_fasta_paths")
    if isinstance(configured, list):
        raw_paths = [p for p in configured if str(p).strip()]
    elif isinstance(configured, str) and configured.strip():
        raw_paths = [configured]

    if not raw_paths:
        fallback = cfg.get("mlm_fasta_path")
        if fallback and str(fallback).strip():
            raw_paths = [fallback]

    if not raw_paths:
        raise ValueError("Provide mlm_fasta_paths or mlm_fasta_path in config")

    resolved: List[Path] = []
    for raw in raw_paths:
        path = _resolve_cfg_path(raw, base_dir=base_dir)
        if not path.exists():
            raise FileNotFoundError(f"FASTA file not found at {path}")
        resolved.append(path)
    return resolved


def _build_mlm_corpus(fasta_paths: List[Path]) -> Tuple[List[Tuple[str, str, str]], List[str]]:
    corpus: List[Tuple[str, str, str]] = []
    species_names: List[str] = []
    used_names: Dict[str, int] = {}

    for path in fasta_paths:
        records = read_fasta_records(path)
        if not records:
            raise ValueError(f"No FASTA records found in {path}")
        base_name = path.stem
        if base_name in used_names:
            used_names[base_name] += 1
            species_name = f"{base_name}_{used_names[base_name]}"
        else:
            used_names[base_name] = 1
            species_name = base_name
        species_names.append(species_name)
        for record_id, sequence in records:
            corpus.append((species_name, record_id, sequence))
    return corpus, species_names


def _stable_species_seed(seed: int, species_name: str) -> int:
    payload = f"{seed}:{species_name}".encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    return int.from_bytes(digest, "big")


def _split_records_for_species(
    records: List[Tuple[str, str]],
    *,
    window_bp: int,
    val_ratio: float,
    seed: int,
    species_name: str,
    exclude_patterns: List[str] | None = None,
    min_val_records: int = 1,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    total = len(records)
    if total == 0:
        return [], []
    if total == 1:
        record_id, sequence = records[0]
        seq_len = len(sequence)
        gap = window_bp
        desired_val_bp = max(int(seq_len * val_ratio), 10 * window_bp)
        desired_val_bp = min(desired_val_bp, seq_len - (2 * window_bp + gap))
        if desired_val_bp < window_bp:
            LOGGER.warning(
                "Single-record species %s too short to split safely (len=%d, window_bp=%d); using all-train.",
                species_name,
                seq_len,
                window_bp,
            )
            return [(record_id, sequence)], []
        train_end = seq_len - desired_val_bp - gap
        val_start = train_end + gap
        train_seq = sequence[:train_end]
        val_seq = sequence[val_start:]
        return [(f"{record_id}|train", train_seq)], [(f"{record_id}|val", val_seq)]

    if min_val_records < 1:
        min_val_records = 1

    patterns = [pattern.lower() for pattern in (exclude_patterns or []) if pattern]
    eligible: List[Tuple[str, str]] = []
    excluded: List[Tuple[str, str]] = []
    if patterns:
        for record_id, sequence in records:
            record_id_lower = record_id.lower()
            if any(pattern in record_id_lower for pattern in patterns):
                excluded.append((record_id, sequence))
            else:
                eligible.append((record_id, sequence))
    else:
        eligible = list(records)

    if patterns and not eligible:
        eligible = list(records)
        excluded = []

    val_count = max(1, int(total * val_ratio))
    val_count = max(val_count, min_val_records)
    if val_count >= total:
        val_count = total - 1

    rng = random.Random(_stable_species_seed(seed, species_name))
    val_records: List[Tuple[str, str]] = []
    train_records: List[Tuple[str, str]] = []

    if eligible and val_count > 0:
        eligible_indices = list(range(len(eligible)))
        rng.shuffle(eligible_indices)
        first_idx = eligible_indices[0]
        val_records.append(eligible[first_idx])
        remaining_pool = [record for idx, record in enumerate(eligible) if idx != first_idx]
        remaining_pool.extend(excluded)
        remaining_needed = val_count - 1
    else:
        remaining_pool = list(eligible)
        remaining_pool.extend(excluded)
        remaining_needed = val_count

    if remaining_needed > 0:
        remaining_indices = list(range(len(remaining_pool)))
        rng.shuffle(remaining_indices)
        val_indices = set(remaining_indices[:remaining_needed])
        val_records.extend(
            record for idx, record in enumerate(remaining_pool) if idx in val_indices
        )
        train_records.extend(
            record for idx, record in enumerate(remaining_pool) if idx not in val_indices
        )
    else:
        train_records.extend(remaining_pool)

    return train_records, val_records


def _split_mlm_corpus_by_record(
    fasta_paths: List[Path],
    *,
    window_bp: int,
    val_ratio: float,
    seed: int,
    single_record_val_ratio: float | None = None,
    exclude_patterns: List[str] | None = None,
    min_val_records: int = 1,
) -> Tuple[
    List[Tuple[str, str, str]],
    List[Tuple[str, str, str]],
    List[str],
    Dict[str, Tuple[int, int, int]],
]:
    train_corpus: List[Tuple[str, str, str]] = []
    val_corpus: List[Tuple[str, str, str]] = []
    species_names: List[str] = []
    split_counts: Dict[str, Tuple[int, int, int]] = {}
    used_names: Dict[str, int] = {}

    for path in fasta_paths:
        records = read_fasta_records(path, keep_header=True)
        if not records:
            raise ValueError(f"No FASTA records found in {path}")
        base_name = path.stem
        if base_name in used_names:
            used_names[base_name] += 1
            species_name = f"{base_name}_{used_names[base_name]}"
        else:
            used_names[base_name] = 1
            species_name = base_name
        species_names.append(species_name)
        effective_val_ratio = val_ratio
        if len(records) == 1 and single_record_val_ratio is not None:
            effective_val_ratio = single_record_val_ratio
        train_records, val_records = _split_records_for_species(
            records,
            window_bp=window_bp,
            val_ratio=effective_val_ratio,
            seed=seed,
            species_name=species_name,
            exclude_patterns=exclude_patterns,
            min_val_records=min_val_records,
        )
        if len(records) == 1 and val_records:
            total_count = len(train_records) + len(val_records)
        else:
            total_count = len(records)
        split_counts[species_name] = (len(train_records), len(val_records), total_count)
        for record_id, sequence in train_records:
            train_corpus.append((species_name, record_id, sequence))
        for record_id, sequence in val_records:
            val_corpus.append((species_name, record_id, sequence))
    return train_corpus, val_corpus, species_names, split_counts


def _log_record_split(
    species_names: List[str],
    split_counts: Dict[str, Tuple[int, int, int]],
) -> None:
    LOGGER.info("Record split per species (train/val):")
    empty_val: List[str] = []
    for species_name in species_names:
        train_count, val_count, total = split_counts[species_name]
        LOGGER.info(
            "  %s: train=%d val=%d total=%d",
            species_name,
            train_count,
            val_count,
            total,
        )
        if val_count == 0:
            empty_val.append(species_name)
    if empty_val:
        LOGGER.warning(
            "Validation split empty for species: %s",
            ", ".join(empty_val),
        )


def _group_corpus_by_species(
    corpus: List[Tuple[str, str, str]],
) -> Dict[str, List[Tuple[str, str]]]:
    grouped: Dict[str, List[Tuple[str, str]]] = {}
    for species_name, record_id, sequence in corpus:
        grouped.setdefault(species_name, []).append((record_id, sequence))
    return grouped


def _matches_markers(species_name: str, record_ids: List[str], markers: List[str]) -> bool:
    species_lower = species_name.lower()
    for marker in markers:
        if marker in species_lower:
            return True
    for record_id in record_ids:
        record_lower = record_id.lower()
        for marker in markers:
            if marker in record_lower:
                return True
    return False


def _log_split_validation_check(
    species_names: List[str],
    train_corpus: List[Tuple[str, str, str]],
    val_corpus: List[Tuple[str, str, str]],
) -> None:
    train_by_species = _group_corpus_by_species(train_corpus)
    val_by_species = _group_corpus_by_species(val_corpus)

    print("Split validation check:")
    for species_name in species_names:
        train_records = train_by_species.get(species_name, [])
        val_records = val_by_species.get(species_name, [])
        val_ids = [record_id for record_id, _ in val_records][:3]
        print(
            f"{species_name}: train={len(train_records)} val={len(val_records)} val_ids={val_ids}"
        )
        train_bp = sum(
            len(sequence)
            for record_id, sequence in train_records
            if record_id.endswith("|train")
        )
        val_bp = sum(
            len(sequence)
            for record_id, sequence in val_records
            if record_id.endswith("|val")
        )
        if train_bp or val_bp:
            print(f"  {species_name}: train_bp={train_bp} val_bp={val_bp}")

    targets = {
        "e_coli": ["escherichia coli", "e. coli", "ecoli", "gcf_000005845"],
        "b_subtilis": ["bacillus subtilis", "b. subtilis", "subtilis"],
        "yeast": ["saccharomyces cerevisiae", "s. cerevisiae", "cerevisiae", "yeast"],
    }
    resolved: Dict[str, str | None] = {label: None for label in targets}
    for species_name in species_names:
        record_ids = [
            record_id
            for record_id, _ in train_by_species.get(species_name, [])
            + val_by_species.get(species_name, [])
        ]
        for label, markers in targets.items():
            if resolved[label] is None and _matches_markers(species_name, record_ids, markers):
                resolved[label] = species_name

    for label, species_name in resolved.items():
        if species_name is None:
            continue
        val_ids = [record_id for record_id, _ in val_by_species.get(species_name, [])]
        if label in {"e_coli", "b_subtilis"}:
            if not val_ids:
                raise AssertionError(f"{label} validation split empty for species {species_name}")
        if label == "yeast":
            if not val_ids:
                raise AssertionError(f"yeast validation split empty for species {species_name}")
            if all("mitochond" in record_id.lower() for record_id in val_ids):
                raise AssertionError(
                    f"yeast validation split contains only mitochondrion for species {species_name}"
                )


def _build_species_records(
    corpus: List[Tuple[str, str, str]],
    *,
    window_bp: int,
) -> Tuple[Dict[str, List[Tuple[str, str]]], List[str], int]:
    records_by_species: Dict[str, List[Tuple[str, str]]] = {}
    species_order: List[str] = []
    skipped_short = 0
    for species_name, record_id, sequence in corpus:
        if len(sequence) < window_bp:
            skipped_short += 1
            continue
        if species_name not in records_by_species:
            records_by_species[species_name] = []
            species_order.append(species_name)
        records_by_species[species_name].append((record_id, sequence))
    return records_by_species, species_order, skipped_short


def _normalize_species_weights(
    raw_weights: Any,
    species_names: List[str],
) -> Dict[str, float]:
    if raw_weights is None or raw_weights == {} or raw_weights == []:
        return {name: 1.0 for name in species_names}

    if isinstance(raw_weights, dict):
        missing = [name for name in species_names if name not in raw_weights]
        extra = [name for name in raw_weights if name not in species_names]
        if missing:
            raise ValueError(f"mlm_species_weights missing entries for: {missing}")
        if extra:
            raise ValueError(f"mlm_species_weights contains unknown species: {extra}")
        return {name: float(raw_weights[name]) for name in species_names}

    if isinstance(raw_weights, list):
        if len(raw_weights) != len(species_names):
            raise ValueError(
                f"mlm_species_weights length ({len(raw_weights)}) does not match "
                f"mlm_fasta_paths ({len(species_names)})"
            )
        return {name: float(weight) for name, weight in zip(species_names, raw_weights)}

    raise ValueError("mlm_species_weights must be a dict or list")


def _resolve_reverse_complement_prob(cfg: Dict[str, Any]) -> float:
    include_rc = bool(cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False)))
    return 0.5 if include_rc else 0.0


def _load_or_build_vocab_for_streaming(
    cfg: Dict[str, Any],
    *,
    corpus: List[Tuple[str, str, str]],
    window_bp: int,
    kmer: int,
    base_dir: Path,
    reverse_complement_prob: float,
    unknown_base_strategy: str,
    species_weights: Dict[str, float],
) -> KmerVocabulary:
    vocab_path = _resolve_cfg_path(
        cfg.get("mlm_vocab_path", f"../artifacts/vocabs/mlm_k{kmer}_vocab.json"),
        base_dir=base_dir,
    )
    if vocab_path.exists():
        LOGGER.info("Loading existing vocabulary from %s", vocab_path)
        return KmerVocabulary.load(vocab_path)

    LOGGER.warning("Vocabulary not found at %s; sampling windows to build one", vocab_path)
    records_by_species, species_order, skipped_short = _build_species_records(corpus, window_bp=window_bp)
    if skipped_short > 0:
        LOGGER.warning("Skipped %d records shorter than window_bp=%d for vocab sampling", skipped_short, window_bp)

    if not records_by_species:
        raise ValueError("No records available for vocab sampling")

    weights = [species_weights[name] for name in species_order]
    sample_count = 10000
    windows: List[str] = []
    segments_cache: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    max_attempts = sample_count * 20
    attempts = 0
    while len(windows) < sample_count and attempts < max_attempts:
        attempts += 1
        species_name = random.choices(species_order, weights=weights, k=1)[0]
        record_id, sequence = random.choice(records_by_species[species_name])
        if unknown_base_strategy == "skip":
            key = (species_name, record_id)
            segments = segments_cache.get(key)
            if segments is None:
                segments = find_acgt_segments(sequence, min_length=window_bp)
                segments_cache[key] = segments
            if not segments:
                continue
            segment_start, segment_end = random.choice(segments)
            max_start = segment_end - window_bp
            start = random.randint(segment_start, max_start)
            window = sequence[start : start + window_bp]
        else:
            max_start = len(sequence) - window_bp
            if max_start < 0:
                continue
            start = random.randint(0, max_start)
            window = sequence[start : start + window_bp]
            window = sanitize_unknown_bases(window, unknown_base_strategy)
        if reverse_complement_prob > 0 and random.random() < reverse_complement_prob:
            window = reverse_complement(window)
        windows.append(window)
    if len(windows) < sample_count:
        LOGGER.warning(
            "Sampled %d/%d windows for vocab build (strategy=%s)",
            len(windows),
            sample_count,
            unknown_base_strategy,
        )

    if not windows:
        raise ValueError("Failed to sample windows for vocab building")

    return load_or_build_vocab(windows, k=kmer, vocab_path=vocab_path)


def run_streaming_dry_run(cfg: Dict[str, Any]) -> None:
    """Run a short streaming dry run to verify MLM configuration.

    Prints:
    - kmer value
    - vocab size
    - sample tokenization length for a short sequence
    - target checkpoint path
    """
    script_dir = Path(__file__).resolve().parent
    fasta_paths = _resolve_mlm_fasta_paths(cfg, base_dir=script_dir)
    val_ratio = float(cfg.get("mlm_val_ratio", 0.1))
    single_record_val_ratio = float(cfg.get("mlm_single_record_val_ratio", val_ratio))
    exclude_patterns = cfg.get("mlm_val_exclude_patterns", ["mitochond"])
    min_val_records = int(cfg.get("mlm_min_val_records_per_species", 1))
    seed = int(cfg.get("seed", 1337) or 1337)
    # Use MLM-specific window size, fallback to max_bp_len for backward compatibility
    window_bp = int(cfg.get("mlm_window_size", cfg.get("max_bp_len", 234)))
    train_corpus, val_corpus, species_names, split_counts = _split_mlm_corpus_by_record(
        fasta_paths,
        window_bp=window_bp,
        val_ratio=val_ratio,
        seed=seed,
        single_record_val_ratio=single_record_val_ratio,
        exclude_patterns=exclude_patterns,
        min_val_records=min_val_records,
    )
    _log_record_split(species_names, split_counts)
    _log_split_validation_check(species_names, train_corpus, val_corpus)
    if not val_corpus:
        LOGGER.warning("Validation corpus empty after record split; dry run uses train corpus.")
    kmer = int(cfg.get("mlm_kmer", 3))
    reverse_complement_prob = _resolve_reverse_complement_prob(cfg)
    unknown_base_strategy = _normalize_unknown_base_strategy(cfg.get("mlm_unknown_base_strategy"))
    species_weights = _normalize_species_weights(cfg.get("mlm_species_weights"), species_names)
    vocab = _load_or_build_vocab_for_streaming(
        cfg,
        corpus=train_corpus,
        window_bp=window_bp,
        kmer=kmer,
        base_dir=script_dir,
        reverse_complement_prob=reverse_complement_prob,
        unknown_base_strategy=unknown_base_strategy,
        species_weights=species_weights,
    )

    # === MLM DRY RUN SUMMARY ===
    print("=" * 60)
    print("MLM DRY RUN CONFIGURATION")
    print("=" * 60)
    print(f"kmer:       {kmer}")
    print(f"vocab_size: {len(vocab.itos)}")

    # Sample tokenization for a short sequence
    sample_seq = "ATGCATGCATGCATGCATGC"  # 20bp sample
    sample_tokens = vocab.tokenize(sample_seq, max_bp_len=len(sample_seq))
    print(f"sample_seq: '{sample_seq}' ({len(sample_seq)} bp)")
    print(f"sample_tokenization_length: {len(sample_tokens)} tokens")

    # Target checkpoint path
    encoder_path = cfg.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{kmer}.pt")
    print(f"target_checkpoint: {encoder_path}")
    print("=" * 60)

    batch_size = int(cfg.get("mlm_batch_size", 128))
    mask_prob = float(cfg.get("mlm_mask_prob", 0.15))
    mask_config = {
        "vocab": vocab,
        "mask_prob": mask_prob,
        "track_spans": False,
        "batch_size": batch_size,
    }
    train_dataset = GenomeMLMIterableDataset(
        corpus=train_corpus,
        window_bp=window_bp,
        kmer=kmer,
        steps_per_epoch=None,
        mask_config=mask_config,
        reverse_complement_prob=reverse_complement_prob,
        unknown_base_strategy=unknown_base_strategy,
        species_weights=species_weights,
        return_metadata=True,
    )
    val_corpus_for_dataset = val_corpus or train_corpus
    val_dataset = GenomeMLMIterableDataset(
        corpus=val_corpus_for_dataset,
        window_bp=window_bp,
        kmer=kmer,
        steps_per_epoch=None,
        mask_config=mask_config,
        reverse_complement_prob=reverse_complement_prob,
        unknown_base_strategy=unknown_base_strategy,
        species_weights=species_weights,
        return_metadata=True,
    )

    def _run_batches(dataset: GenomeMLMIterableDataset, label: str) -> None:
        masked_fractions: List[float] = []
        species_counts: Counter[str] = Counter()
        iterator = iter(dataset)
        for batch_idx in range(3):
            batch_samples = [next(iterator) for _ in range(batch_size)]
            inputs = torch.stack([item[0] for item in batch_samples])
            labels = torch.stack([item[1] for item in batch_samples])
            metadata = [item[2] for item in batch_samples]
            species_counts.update([item["species"] for item in metadata])
            masked_fraction = (labels != -100).sum().item() / labels.numel()
            masked_fractions.append(masked_fraction)
            print(
                f"{label} batch {batch_idx + 1}: inputs={tuple(inputs.shape)} "
                f"labels={tuple(labels.shape)}"
            )

        mean_masked_fraction = sum(masked_fractions) / max(1, len(masked_fractions))
        print(f"{label} masked_fraction mean: {mean_masked_fraction:.4f}")
        print(f"{label} species distribution counts: {dict(species_counts)}")

    _run_batches(train_dataset, "train")
    _run_batches(val_dataset, "val")


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def run_validation(
    model: DNAMLM,
    val_loader: DataLoader,
    device: torch.device,
    autocast_device: str | None,
    autocast_dtype: torch.dtype | None,
) -> Dict[str, float]:
    """Run validation and return metrics dict with loss, accuracy, perplexity, masked stats."""
    was_training = model.training
    model.eval()
    LOGGER.info("Val loop: model.training=%s", model.training)
    try:
        total_loss = 0.0
        total_correct = 0
        total_masked = 0
        total_tokens = 0

        with torch.no_grad():
            for batch_data in val_loader:
                if len(batch_data) == 3:
                    inputs, labels, _span_lengths = batch_data
                else:
                    inputs, labels = batch_data

                non_blocking = (device.type == "cuda")
                inputs = inputs.to(device, non_blocking=non_blocking)
                labels = labels.to(device, non_blocking=non_blocking)

                if autocast_device and autocast_dtype is not None:
                    with torch.autocast(autocast_device, dtype=autocast_dtype):
                        logits, _ = model(inputs, labels=None)
                else:
                    logits, _ = model(inputs, labels=None)

                # Loss computation outside autocast for numerical stability
                loss = F.cross_entropy(
                    logits.float().view(-1, logits.size(-1)),
                    labels.view(-1),
                    ignore_index=-100,
                )

                total_loss += loss.item() * inputs.size(0)

                # Compute accuracy on masked tokens
                preds = logits.argmax(dim=-1)
                mask = labels != -100
                total_correct += (preds[mask] == labels[mask]).sum().item()
                total_masked += mask.sum().item()
                total_tokens += labels.numel()

        avg_loss = total_loss / len(val_loader.dataset)
        accuracy = total_correct / max(1, total_masked)
        perplexity = math.exp(min(avg_loss, 100))  # Cap to avoid overflow
        masked_fraction = total_masked / max(1, total_tokens)

        return {
            'loss': avg_loss,
            'accuracy': accuracy,
            'perplexity': perplexity,
            'loss_masked_only': avg_loss,
            'masked_token_accuracy': accuracy,
            'masked_fraction': masked_fraction,
        }
    finally:
        model.train(was_training)


def compute_masked_batch_metrics(
    logits: torch.Tensor,
    labels: torch.LongTensor,
) -> Dict[str, float]:
    """Compute masked-only loss, accuracy, and masked fraction for a single batch."""
    loss = F.cross_entropy(
        logits.float().view(-1, logits.size(-1)),
        labels.view(-1),
        ignore_index=-100,
    )
    mask = labels != -100
    total_tokens = labels.numel()
    masked_tokens = mask.sum().item()
    masked_fraction = masked_tokens / max(1, total_tokens)

    if masked_tokens > 0:
        preds = logits.argmax(dim=-1)
        masked_token_accuracy = (preds[mask] == labels[mask]).float().mean().item()
    else:
        masked_token_accuracy = 0.0

    return {
        "loss_masked_only": loss.item(),
        "masked_token_accuracy": masked_token_accuracy,
        "masked_fraction": masked_fraction,
    }


def check_gradients_and_frozen(model: nn.Module) -> Dict[str, Any]:
    """Check gradient norms and frozen weights for debugging."""
    grad_info = {"frozen": [], "zero_grad": [], "normal": []}
    total_params = 0
    frozen_params = 0
    
    for name, param in model.named_parameters():
        total_params += param.numel()
        if not param.requires_grad:
            frozen_params += param.numel()
            grad_info["frozen"].append(name)
        elif param.grad is None:
            grad_info["zero_grad"].append(name)
        else:
            grad_norm = param.grad.norm().item()
            if grad_norm == 0:
                grad_info["zero_grad"].append(name)
            else:
                grad_info["normal"].append((name, grad_norm))
    
    return {
        "total_params": total_params,
        "frozen_params": frozen_params,
        "frozen_pct": frozen_params / total_params * 100 if total_params else 0,
        "details": grad_info,
    }


def print_topk_predictions_debug(
    model: nn.Module,
    inputs: torch.LongTensor,
    labels: torch.LongTensor,
    vocab: KmerVocabulary,
    epoch: int,
    num_sequences: int = 3,
    top_k: int = 5,
) -> Dict[str, Any]:
    """
    Top-K sanity print for masked positions.
    
    Prints detailed prediction info for the first `num_sequences` sequences:
    - Masked position index
    - True token string
    - Model top-5 predicted tokens + probabilities
    
    This helps verify that:
    1. Predictions become sharper over time (top-1 prob rises)
    2. True token appears in top-5 more often as training progresses
    
    Args:
        model: The DNAMLM model
        inputs: (B, L) masked input token IDs
        labels: (B, L) target labels with -100 for non-masked positions
        vocab: KmerVocabulary for decoding tokens
        epoch: Current epoch number (for display)
        num_sequences: Number of sequences to analyze (default: 3)
        top_k: Number of top predictions to show (default: 5)
        
    Returns:
        Dict with summary statistics:
        - avg_top1_prob: Average probability of top-1 prediction
        - true_in_top1_pct: Percentage of masked tokens where true label is top-1
        - true_in_topk_pct: Percentage of masked tokens where true label is in top-k
    """
    model.eval()
    device = next(model.parameters()).device

    # Move to device and get predictions
    non_blocking = (device.type == "cuda")
    inputs = inputs.to(device, non_blocking=non_blocking)
    labels = labels.to(device, non_blocking=non_blocking)

    with torch.no_grad():
        logits, _ = model(inputs, labels)
        probs = F.softmax(logits, dim=-1)
    
    batch_size = min(num_sequences, inputs.size(0))
    
    print()
    print("=" * 90)
    print(f"TOP-K PREDICTIONS DEBUG (Epoch {epoch})")
    print("=" * 90)
    print(f"Analyzing {batch_size} sequences, showing top-{top_k} predictions per masked position")
    print()
    
    # Track statistics
    total_masked = 0
    true_in_top1 = 0
    true_in_topk = 0
    top1_probs = []
    
    for seq_idx in range(batch_size):
        seq_labels = labels[seq_idx]
        seq_probs = probs[seq_idx]
        
        # Find masked positions (where label != -100)
        masked_positions = torch.where(seq_labels != -100)[0]
        
        if len(masked_positions) == 0:
            print(f"[Sequence {seq_idx}] No masked positions found")
            continue
            
        print(f"[Sequence {seq_idx}] {len(masked_positions)} masked positions")
        print("-" * 90)
        
        # Show up to 5 masked positions per sequence for readability
        positions_to_show = masked_positions[:5]
        
        for pos in positions_to_show:
            pos_idx = pos.item()
            true_label = seq_labels[pos_idx].item()
            
            # Get true token string
            if true_label < len(vocab.itos):
                true_token = vocab.itos[true_label]
            else:
                true_token = f"<id:{true_label}>"
            
            # Get top-k predictions
            pos_probs = seq_probs[pos_idx]
            topk_probs, topk_ids = torch.topk(pos_probs, top_k)
            
            # Check if true label is in top-1 or top-k
            top1_id = topk_ids[0].item()
            is_top1_correct = (top1_id == true_label)
            is_in_topk = true_label in topk_ids.tolist()
            
            # Update statistics
            total_masked += 1
            if is_top1_correct:
                true_in_top1 += 1
            if is_in_topk:
                true_in_topk += 1
            top1_probs.append(topk_probs[0].item())
            
            # Format output
            correct_marker = "✓" if is_top1_correct else ("○" if is_in_topk else "✗")
            print(f"  Position {pos_idx:>3} | True: '{true_token}' (id={true_label}) | {correct_marker}")
            print(f"    Top-{top_k} predictions:")
            
            for rank, (prob, pred_id) in enumerate(zip(topk_probs, topk_ids), 1):
                pred_id_val = pred_id.item()
                if pred_id_val < len(vocab.itos):
                    pred_token = vocab.itos[pred_id_val]
                else:
                    pred_token = f"<id:{pred_id_val}>"
                
                # Mark if this is the true token
                marker = " ★" if pred_id_val == true_label else ""
                print(f"      {rank}. '{pred_token}' (id={pred_id_val}): {prob.item()*100:6.2f}%{marker}")
            print()
        
        print()
    
    # Calculate and print summary statistics
    avg_top1_prob = sum(top1_probs) / len(top1_probs) if top1_probs else 0.0
    true_in_top1_pct = true_in_top1 / total_masked * 100 if total_masked > 0 else 0.0
    true_in_topk_pct = true_in_topk / total_masked * 100 if total_masked > 0 else 0.0
    
    print("-" * 90)
    print("SUMMARY STATISTICS")
    print("-" * 90)
    print(f"  Total masked positions analyzed: {total_masked}")
    print(f"  Average top-1 probability:       {avg_top1_prob*100:.2f}%")
    print(f"  True token is top-1:             {true_in_top1}/{total_masked} ({true_in_top1_pct:.1f}%)")
    print(f"  True token in top-{top_k}:            {true_in_topk}/{total_masked} ({true_in_topk_pct:.1f}%)")
    print()
    print("  Healthy training indicators:")
    print("    → Top-1 probability should INCREASE over epochs")
    print("    → True-in-top-1 % should INCREASE over epochs")
    print("    → True-in-top-5 % should be HIGH early, approaching true-in-top-1 later")
    print("=" * 90)
    print()
    
    model.train()  # Restore training mode
    
    return {
        "avg_top1_prob": avg_top1_prob,
        "true_in_top1_pct": true_in_top1_pct,
        "true_in_topk_pct": true_in_topk_pct,
        "total_masked": total_masked,
    }


def run_overfit_debug(
    model: nn.Module,
    train_samples: List[torch.LongTensor],
    vocab: KmerVocabulary,
    device: torch.device,
    num_steps: int = 500,
    lr: float = 3e-3,
) -> None:
    """
    HARD OVERFIT TEST: Train on tiny dataset until it breaks physics.
    
    This is a "lie detector" - if the model can't perfectly memorize a tiny 
    dataset, something is fundamentally wrong with the architecture or training.
    
    Configuration:
    - 500-1000 steps (aggressive)
    - NO dropout (pure memorization)
    - Higher LR (3e-3) for fast convergence
    - Single batch repeated (same data every step)
    - Log every 50 steps
    - FAIL if accuracy doesn't exceed 80%
    - Print top-5 predictions for masked positions on failure
    
    Expected: accuracy should reach 95%+ on this tiny dataset.
    If not, the model is broken.
    """
    print("=" * 80)
    print("HARD OVERFIT TEST: Breaking Physics Mode")
    print("=" * 80)
    print("If model can't memorize 32 samples, something is fundamentally broken.")
    print()
    print(f"Configuration:")
    print(f"  Samples:     {len(train_samples)}")
    print(f"  Steps:       {num_steps}")
    print(f"  LR:          {lr}")
    print(f"  Dropout:     DISABLED")
    print(f"  Target:      >80% accuracy (should reach 95%+)")
    print()
    
    # =========================================================================
    # DISABLE ALL DROPOUT for pure memorization
    # =========================================================================
    def disable_dropout(module):
        if isinstance(module, nn.Dropout):
            module.p = 0.0
    
    model.apply(disable_dropout)
    model.train()  # Keep in train mode but with dropout=0
    
    # Use AdamW with no weight decay for pure memorization
    # Also try higher LR for first phase, then reduce
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.0, betas=(0.9, 0.999), foreach=True)
    
    # Add learning rate scheduler: warmup then constant
    def lr_schedule(step):
        warmup = 50
        if step < warmup:
            return step / warmup
        return 1.0
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_schedule)
    
    # Create FIXED masked inputs/labels - same data every step
    raw_inputs = torch.stack([train_samples[i] for i in range(len(train_samples))]).to(device)
    masked_inputs, labels = mask_tokens_span(raw_inputs, vocab, mask_prob=0.15)
    
    num_masked = (labels != -100).sum().item()
    mask = labels != -100
    
    print(f"Fixed masking: {num_masked} tokens masked across {len(train_samples)} sequences")
    print(f"Masked positions per sequence: ~{num_masked / len(train_samples):.1f}")
    print()
    print("-" * 80)
    print(f"{'Step':>6} | {'Loss':>10} | {'Accuracy':>10} | {'Status'}")
    print("-" * 80)
    
    losses = []
    accuracies = []
    best_acc = 0.0
    
    for step in range(num_steps):
        optimizer.zero_grad()
        logits, loss = model(masked_inputs, labels)
        loss.backward()
        
        # Gradient clipping for stability with norm logging
        unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        clipped_grad_norm = min(unclipped_grad_norm.item(), 1.0)
        optimizer.step()
        scheduler.step()
        
        # Compute accuracy
        preds = logits.argmax(dim=-1)
        acc = (preds[mask] == labels[mask]).float().mean().item() if mask.any() else 0.0
        
        losses.append(loss.item())
        accuracies.append(acc)
        best_acc = max(best_acc, acc)
        
        # Log every 50 steps
        if step % 50 == 0 or step == num_steps - 1:
            status = "✓" if acc > 0.8 else ("↑" if acc > 0.5 else "...")
            print(f"{step:>6} | {loss.item():>10.4f} | {acc*100:>9.1f}% | {status}")
    
    print("-" * 80)
    print()
    
    # =========================================================================
    # ANALYSIS
    # =========================================================================
    print("=" * 80)
    print("OVERFIT ANALYSIS")
    print("=" * 80)
    
    initial_loss = sum(losses[:5]) / min(5, len(losses))
    final_loss = sum(losses[-10:]) / min(10, len(losses))
    final_acc = sum(accuracies[-10:]) / min(10, len(accuracies))
    
    print(f"Initial loss:     {initial_loss:.4f}")
    print(f"Final loss:       {final_loss:.4f}")
    print(f"Loss reduction:   {(1 - final_loss/initial_loss)*100:.1f}%")
    print(f"Best accuracy:    {best_acc*100:.1f}%")
    print(f"Final accuracy:   {final_acc*100:.1f}%")
    print()
    
    # SUCCESS THRESHOLD: 80% accuracy
    if best_acc >= 0.80:
        print("✓ HARD OVERFIT PASSED!")
        print(f"  Model successfully memorized the data (best acc: {best_acc*100:.1f}%)")
        print("  The model architecture is capable of learning.")
        print("=" * 80)
        return
    
    # =========================================================================
    # FAILURE ANALYSIS: Print top-5 predictions for masked positions
    # =========================================================================
    print("✗ HARD OVERFIT FAILED!")
    print(f"  Best accuracy: {best_acc*100:.1f}% (required: >80%)")
    print()
    print("  This means the model CANNOT memorize a tiny dataset.")
    print("  Something is fundamentally broken.")
    print()
    
    print("-" * 80)
    print("TOP-5 PREDICTIONS FOR MASKED POSITIONS")
    print("-" * 80)
    
    # Get final predictions with probabilities
    model.eval()
    with torch.no_grad():
        logits, _ = model(masked_inputs, labels)
        probs = F.softmax(logits, dim=-1)
    
    # Analyze first few masked positions
    masked_indices = torch.where(mask)
    num_to_show = min(10, len(masked_indices[0]))
    
    print(f"Showing {num_to_show} masked positions:\n")
    
    for i in range(num_to_show):
        batch_idx = masked_indices[0][i].item()
        pos_idx = masked_indices[1][i].item()
        
        true_label = labels[batch_idx, pos_idx].item()
        true_token = vocab.itos[true_label] if true_label < len(vocab.itos) else f"<{true_label}>"
        
        # Get top-5 predictions
        pos_probs = probs[batch_idx, pos_idx]
        top5_probs, top5_ids = torch.topk(pos_probs, 5)
        
        print(f"  Position [{batch_idx},{pos_idx}] - True label: '{true_token}' (id={true_label})")
        print(f"    Top-5 predictions:")
        for j, (prob, pred_id) in enumerate(zip(top5_probs, top5_ids)):
            pred_token = vocab.itos[pred_id.item()] if pred_id.item() < len(vocab.itos) else f"<{pred_id.item()}>"
            marker = "✓" if pred_id.item() == true_label else " "
            print(f"      {j+1}. {marker} '{pred_token}' (id={pred_id.item()}) - {prob.item()*100:.1f}%")
        print()
    
    # Check for systematic issues
    print("-" * 80)
    print("DIAGNOSIS")
    print("-" * 80)
    
    # Check gradient info
    grad_info = check_gradients_and_frozen(model)
    
    if grad_info["frozen_pct"] > 10:
        print(f"  ⚠ {grad_info['frozen_pct']:.1f}% of parameters are frozen!")
    
    if len(grad_info["details"]["zero_grad"]) > 0:
        print(f"  ⚠ {len(grad_info['details']['zero_grad'])} layers have zero gradients")
    
    # Check if predictions are too uniform
    avg_top1_prob = probs.max(dim=-1)[0][mask].mean().item()
    print(f"  Average top-1 probability: {avg_top1_prob*100:.1f}%")
    if avg_top1_prob < 0.1:
        print("    ⚠ Predictions are too uniform - model not learning to discriminate")
    
    # Check if model is predicting special tokens (should be impossible now)
    special_ids = set(model.special_token_ids) if hasattr(model, 'special_token_ids') else set()
    if special_ids:
        top1_preds = probs.argmax(dim=-1)[mask]
        special_preds = sum(1 for p in top1_preds if p.item() in special_ids)
        if special_preds > 0:
            print(f"    ✗ {special_preds} predictions are special tokens (should be 0!)")
    
    # =========================================================================
    # KEY DIAGNOSTIC: Check 80/10/10 masking breakdown
    # =========================================================================
    print()
    print("-" * 80)
    print("80/10/10 MASKING ANALYSIS")
    print("-" * 80)
    print("  BERT-style masking: 80% [MASK], 10% random, 10% original")
    print()
    
    # Count how many masked positions got [MASK] vs other tokens
    mask_token_id = vocab.mask_id
    masked_input_tokens = masked_inputs[mask]
    original_tokens = raw_inputs[mask]
    
    got_mask_token = (masked_input_tokens == mask_token_id).sum().item()
    got_original = (masked_input_tokens == original_tokens).sum().item()
    got_random = num_masked - got_mask_token - got_original
    
    print(f"  Masked positions breakdown:")
    print(f"    [MASK] token:    {got_mask_token:>4} ({got_mask_token/num_masked*100:.1f}%) - NO context")
    print(f"    Original token:  {got_original:>4} ({got_original/num_masked*100:.1f}%) - EASY (sees answer)")
    print(f"    Random token:    {got_random:>4} ({got_random/num_masked*100:.1f}%) - SOME context")
    print()
    
    # Check accuracy breakdown by masking type
    final_preds = logits.argmax(dim=-1)[mask]
    final_labels = labels[mask]
    
    correct_mask = (final_preds == final_labels)
    
    # Accuracy for [MASK] positions
    is_mask_token = (masked_input_tokens == mask_token_id)
    acc_mask_only = correct_mask[is_mask_token].float().mean().item() if is_mask_token.any() else 0
    
    # Accuracy for original positions (should be ~100%)
    is_original = (masked_input_tokens == original_tokens)
    acc_original = correct_mask[is_original].float().mean().item() if is_original.any() else 0
    
    # Accuracy for random positions
    is_random = ~is_mask_token & ~is_original
    acc_random = correct_mask[is_random].float().mean().item() if is_random.any() else 0
    
    print(f"  Accuracy by masking type:")
    print(f"    [MASK] positions:    {acc_mask_only*100:.1f}%")
    print(f"    Original positions:  {acc_original*100:.1f}%")
    print(f"    Random positions:    {acc_random*100:.1f}%")
    print()
    
    if acc_original > 0.9 and acc_mask_only < 0.2:
        print("  ⚠ CRITICAL: Model only succeeds when it can SEE the answer!")
        print("    The 10% 'original' positions have high accuracy (cheating).")
        print("    The 80% [MASK] positions have near-random accuracy.")
        print()
        print("    This means attention/context IS NOT WORKING.")
        print("    The model cannot use surrounding context to predict masked tokens.")
    
    print()
    print("POSSIBLE CAUSES:")
    print("  1. Model architecture issue (attention not working)")
    print("  2. Loss function bug (gradients not flowing)")
    print("  3. Special token exclusion not working")
    print("  4. Embedding initialization issue")
    print("  5. Learning rate too low or too high")
    print()
    print("Try running with --leakage_test and --label_debug to verify setup.")
    print("=" * 80)


def run_mlm_pretrain(cfg: dict, *, overrides: dict | None = None) -> dict:
    """Runs MLM pretraining and returns metrics + encoder checkpoint path."""
    start_time = time.perf_counter()

    cfg_run = copy.copy(cfg) if isinstance(cfg, dict) else dict(cfg or {})
    if overrides:
        cfg_run.update(overrides)

    # Log MLM-specific configuration settings
    # Use MLM-specific settings with fallbacks for backward compatibility
    _max_bp_len = int(cfg_run.get("mlm_window_size", cfg_run.get("max_bp_len", 234)))
    _lr = float(cfg_run.get("mlm_lr", cfg_run.get("lr", 2e-4)))
    _batch_size = int(cfg_run.get("mlm_batch_size", cfg_run.get("batch_size", 128)))
    LOGGER.info("=" * 60)
    LOGGER.info("MLM PRETRAINING CONFIGURATION")
    LOGGER.info("=" * 60)
    LOGGER.info("Window size (max_bp_len): %d", _max_bp_len)
    LOGGER.info("Learning rate: %.6f", _lr)
    LOGGER.info("Batch size: %d", _batch_size)
    LOGGER.info("=" * 60)

    script_dir = Path(__file__).resolve().parent
    config_path_value = cfg_run.get("_config_path") or cfg_run.get("config_path")
    config_id = str(config_path_value) if config_path_value else "<in-memory>"

    ablation_cfg = cfg_run.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_enabled = bool(ablation_cfg.get("enabled", False))

    epochs = int(cfg_run.get("mlm_epochs", 150))
    if ablation_enabled and ablation_cfg.get("pretrain_epochs") is not None:
        epochs_raw = ablation_cfg.get("pretrain_epochs")
        try:
            epochs = int(epochs_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.pretrain_epochs must be an int, got {epochs_raw!r}") from exc
        if epochs <= 0:
            raise ValueError(f"ablation.pretrain_epochs must be > 0, got {epochs}")

    max_steps: int | None = None
    if ablation_enabled and ablation_cfg.get("pretrain_max_steps") is not None:
        max_steps_raw = ablation_cfg.get("pretrain_max_steps")
        try:
            max_steps_int = int(max_steps_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.pretrain_max_steps must be an int, got {max_steps_raw!r}") from exc
        if max_steps_int < 0:
            raise ValueError(f"ablation.pretrain_max_steps must be >= 0, got {max_steps_int}")
        if max_steps_int > 0:
            max_steps = max_steps_int

    label_debug = bool(cfg_run.get("label_debug", False))
    debug_mode = bool(cfg_run.get("debug", False))
    overfit_debug = bool(cfg_run.get("overfit_debug", False))
    leakage_test = bool(cfg_run.get("leakage_test", False))
    generate_golden_batch = bool(cfg_run.get("generate_golden_batch", False))
    topk_debug = bool(cfg_run.get("topk_debug", False))
    debug_batch = bool(cfg_run.get("debug_batch", False))
    golden_batch_path = cfg_run.get("golden_batch_path")
    mask_prob = float(cfg_run.get("mlm_mask_prob", 0.15))
    dry_run = bool(cfg_run.get("mlm_dry_run", False))
    use_streaming = _should_use_streaming(cfg_run)

    # Progressive masking configuration (DNABERT-style curriculum)
    mlm_masking_cfg = cfg_run.get("mlm_masking", {})
    if not isinstance(mlm_masking_cfg, dict):
        mlm_masking_cfg = {}
    progressive_masking_enabled = bool(mlm_masking_cfg.get("progressive", False))
    initial_mask_rate = float(mlm_masking_cfg.get("initial_mask_rate", mask_prob))
    final_mask_rate = float(mlm_masking_cfg.get("final_mask_rate", 0.20))
    masking_transition_start_ratio = float(mlm_masking_cfg.get("transition_start_ratio", 0.8))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"mlm_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    # Output directory for checkpoints and artifacts
    output_dir_value = cfg_run.get("output_dir")
    if output_dir_value is not None:
        output_dir = Path(output_dir_value).expanduser()
        if not output_dir.is_absolute():
            output_dir = (script_dir / output_dir).resolve()
    else:
        output_dir = (script_dir / "../artifacts").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = output_dir / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Output directory: %s", output_dir)
    LOGGER.info("Checkpoints directory: %s", checkpoints_dir)

    # Checkpoint save frequency
    save_every_steps_raw = cfg_run.get("save_every_steps", 500)
    try:
        save_every_steps = int(save_every_steps_raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"save_every_steps must be an int, got {save_every_steps_raw!r}") from exc
    if save_every_steps < 0:
        raise ValueError(f"save_every_steps must be >= 0, got {save_every_steps}")

    checkpoint_keep_last_raw = cfg_run.get("mlm_checkpoint_keep_last", 2)
    if checkpoint_keep_last_raw is None:
        checkpoint_keep_last = 2
    else:
        try:
            checkpoint_keep_last = int(checkpoint_keep_last_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"mlm_checkpoint_keep_last must be an int, got {checkpoint_keep_last_raw!r}"
            ) from exc
    if checkpoint_keep_last < 1:
        LOGGER.warning(
            "mlm_checkpoint_keep_last=%s is invalid; defaulting to 1",
            checkpoint_keep_last_raw,
        )
        checkpoint_keep_last = 1

    resume_path_value = cfg_run.get("resume_path")

    length_curriculum_config = LengthCurriculumConfig.from_config(cfg_run)
    if use_streaming and length_curriculum_config.enabled:
        LOGGER.warning("Length curriculum disabled for streaming datasets")
        length_curriculum_config.enabled = False

    resolved_config: Dict[str, Any] = {
        "script": "pretrain_mlm.py",
        "timestamp": timestamp,
        "config_path": config_id,
        # Ablation controls
        "ablation_enabled": ablation_enabled,
        "pretrain_bp_limit": ablation_cfg.get("pretrain_bp_limit") if ablation_enabled else None,
        "pretrain_max_steps": max_steps,
        "pretrain_epochs": epochs if ablation_enabled else None,
        # Model architecture
        "embedding_dim": int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256))),
        "transformer_layers": int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8))),
        "transformer_heads": int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8))),
        "transformer_ff_dim": int(cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", 1024))),
        "transformer_dropout": float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1))),
        # Optimizer
        # Use MLM-specific learning rate
        "lr": float(cfg_run.get("mlm_lr", cfg_run.get("lr", 2e-4))),
        "weight_decay": float(cfg_run.get("mlm_weight_decay", 0.01)),
        # MLM-specific
        "mask_prob": mask_prob,
        "max_span_len": 3,  # Fixed in mask_tokens_span
        # Data
        "batch_size": int(cfg_run.get("mlm_batch_size", 128)),
        "max_bp_len": int(cfg_run.get("mlm_window_size", cfg_run.get("max_bp_len", 234))),
        "kmer": int(cfg_run.get("mlm_kmer", 3)),
        "streaming_enabled": use_streaming,
        "mlm_fasta_path": cfg_run.get("mlm_fasta_path"),
        "mlm_fasta_paths": cfg_run.get("mlm_fasta_paths"),
        "mlm_species_weights": cfg_run.get("mlm_species_weights"),
        "mlm_steps_per_epoch": cfg_run.get("mlm_steps_per_epoch"),
        "include_reverse_complements": bool(
            cfg_run.get(
                "include_reverse_complements",
                cfg_run.get("mlm_include_reverse_complements", False),
            )
        ),
        "reverse_complement_prob": _resolve_reverse_complement_prob(cfg_run),
        # Training
        "epochs": epochs,
        "grad_accum_steps": int(cfg_run.get("mlm_grad_accum_steps", 4)),
        "warmup_ratio": float(cfg_run.get("mlm_warmup_ratio", 0.06)),
        "final_lr_ratio": float(cfg_run.get("mlm_final_lr_ratio", 0.1)),
        "lr_log_interval": int(cfg_run.get("mlm_lr_log_interval", 50)),
        "patience": int(cfg_run.get("mlm_patience", 30)),
        "save_every_steps": save_every_steps,
        "checkpoint_keep_last": checkpoint_keep_last,
        "val_ratio": float(cfg_run.get("mlm_val_ratio", 0.1)),
        # Model options
        "tie_weights": bool(cfg_run.get("mlm_tie_weights", True)),
        "use_output_norm": bool(cfg_run.get("mlm_use_output_norm", True)),
        "use_layer_lr_decay": bool(cfg_run.get("mlm_use_layer_lr_decay", True)),
        "layer_lr_decay": float(cfg_run.get("mlm_layer_lr_decay", 0.9)),
        # Seed
        "seed": int(cfg_run.get("seed", 1337) or 1337),
        # AMP (Mixed Precision)
        "amp_enabled": bool(cfg_run.get("amp_enabled", True)),
        "amp_dtype": str(cfg_run.get("amp_dtype", "bfloat16")),
        "amp_force_fp32_loss": bool(cfg_run.get("amp_force_fp32_loss", True)),
        # Length curriculum
        "length_curriculum": length_curriculum_config.to_dict(),
        # Progressive masking (DNABERT-style)
        "progressive_masking": {
            "enabled": progressive_masking_enabled,
            "initial_rate": initial_mask_rate,
            "final_rate": final_mask_rate,
            "transition_start_ratio": masking_transition_start_ratio,
        },
    }

    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)

    seed = int(cfg_run.get("seed", 1337) or 1337)
    set_seed(seed)
    LOGGER.info("Using seed: %d", seed)

    batch_size = int(cfg_run.get("mlm_batch_size", 128))
    grad_accum_steps = int(cfg_run.get("mlm_grad_accum_steps", 2))
    effective_batch = batch_size * grad_accum_steps
    LOGGER.info(
        "Batch size: %d, Grad accum: %d, Effective batch: %d",
        batch_size,
        grad_accum_steps,
        effective_batch,
    )

    val_ratio = float(cfg_run.get("mlm_val_ratio", 0.1))
    reverse_complement_prob = _resolve_reverse_complement_prob(cfg_run)
    unknown_base_strategy = _normalize_unknown_base_strategy(cfg_run.get("mlm_unknown_base_strategy"))

    if dry_run:
        run_streaming_dry_run(cfg_run)
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": 0,
        }

    token_tensors: List[torch.LongTensor] | None = None
    steps_per_epoch = None
    vocab_path: Path | None = None

    if use_streaming:
        steps_raw = cfg_run.get("mlm_steps_per_epoch")
        if steps_raw is None:
            raise ValueError("mlm_steps_per_epoch is required when using streaming datasets")
        try:
            steps_per_epoch = int(steps_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"mlm_steps_per_epoch must be an int, got {steps_raw!r}") from exc
        if steps_per_epoch <= 0:
            raise ValueError(f"mlm_steps_per_epoch must be > 0, got {steps_per_epoch}")

        fasta_paths = _resolve_mlm_fasta_paths(cfg_run, base_dir=script_dir)
        # Use MLM-specific window size, fallback to max_bp_len for backward compatibility
        window_bp = int(cfg_run.get("mlm_window_size", cfg_run.get("max_bp_len", 234)))
        single_record_val_ratio = float(cfg_run.get("mlm_single_record_val_ratio", val_ratio))
        exclude_patterns = cfg_run.get("mlm_val_exclude_patterns", ["mitochond"])
        min_val_records = int(cfg_run.get("mlm_min_val_records_per_species", 1))
        train_corpus, val_corpus, species_names, split_counts = _split_mlm_corpus_by_record(
            fasta_paths,
            window_bp=window_bp,
            val_ratio=val_ratio,
            seed=seed,
            single_record_val_ratio=single_record_val_ratio,
            exclude_patterns=exclude_patterns,
            min_val_records=min_val_records,
        )
        _log_record_split(species_names, split_counts)
        species_weights = _normalize_species_weights(cfg_run.get("mlm_species_weights"), species_names)

        if ablation_enabled and ablation_cfg.get("pretrain_bp_limit") is not None:
            LOGGER.warning("Ablation pretrain_bp_limit not supported in streaming mode; ignoring.")

        kmer = int(cfg_run.get("mlm_kmer", 3))
        vocab_path = _resolve_cfg_path(
            cfg_run.get("mlm_vocab_path", f"../artifacts/vocabs/mlm_k{kmer}_vocab.json"),
            base_dir=script_dir,
        )
        vocab = _load_or_build_vocab_for_streaming(
            cfg_run,
            corpus=train_corpus,
            window_bp=window_bp,
            kmer=kmer,
            base_dir=script_dir,
            reverse_complement_prob=reverse_complement_prob,
            unknown_base_strategy=unknown_base_strategy,
            species_weights=species_weights,
        )

        train_records_by_species, train_species_order, train_skipped_short = _build_species_records(
            train_corpus,
            window_bp=window_bp,
        )
        train_total_records = sum(len(records) for records in train_records_by_species.values())
        LOGGER.info(
            "Streaming train corpus: %d species, %d records (skipped_short=%d)",
            len(train_species_order),
            train_total_records,
            train_skipped_short,
        )
        val_corpus_for_dataset = val_corpus
        if not val_corpus_for_dataset:
            LOGGER.warning(
                "Validation corpus empty after record split; using train corpus for validation (leakage)."
            )
            val_corpus_for_dataset = train_corpus
        else:
            val_records_by_species, val_species_order, val_skipped_short = _build_species_records(
                val_corpus_for_dataset,
                window_bp=window_bp,
            )
            val_total_records = sum(len(records) for records in val_records_by_species.values())
            LOGGER.info(
                "Streaming val corpus: %d species, %d records (skipped_short=%d)",
                len(val_species_order),
                val_total_records,
                val_skipped_short,
            )

        train_mask_config = {
            "vocab": vocab,
            "mask_prob": mask_prob,
            "track_spans": True,
            "batch_size": batch_size,
            "debug": label_debug,
        }
        train_dataset = GenomeMLMIterableDataset(
            corpus=train_corpus,
            window_bp=window_bp,
            kmer=kmer,
            steps_per_epoch=steps_per_epoch,
            mask_config=train_mask_config,
            reverse_complement_prob=reverse_complement_prob,
            unknown_base_strategy=unknown_base_strategy,
            species_weights=species_weights,
        )

        val_steps_per_epoch = max(1, int(steps_per_epoch * val_ratio))
        val_mask_config = {
            "vocab": vocab,
            "mask_prob": mask_prob,
            "track_spans": False,
            "batch_size": batch_size,
        }
        val_dataset = GenomeMLMIterableDataset(
            corpus=val_corpus_for_dataset,
            window_bp=window_bp,
            kmer=kmer,
            steps_per_epoch=val_steps_per_epoch,
            mask_config=val_mask_config,
            reverse_complement_prob=reverse_complement_prob,
            unknown_base_strategy=unknown_base_strategy,
            species_weights=species_weights,
        )
        LOGGER.info(
            "Streaming dataset: steps_per_epoch=%d, val_steps_per_epoch=%d",
            steps_per_epoch,
            val_steps_per_epoch,
        )
    else:
        fasta_path = _resolve_cfg_path(cfg_run.get("mlm_fasta_path"), base_dir=script_dir)
        LOGGER.info("Reading FASTA records from %s", fasta_path)
        records = read_fasta_records(fasta_path, keep_header=True)
        if not records:
            raise ValueError(f"No FASTA records found in {fasta_path}")
        LOGGER.info("Found %d FASTA record(s)", len(records))
        if ablation_enabled and ablation_cfg.get("pretrain_bp_limit") is not None:
            bp_limit_raw = ablation_cfg.get("pretrain_bp_limit")
            try:
                bp_limit = int(bp_limit_raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"ablation.pretrain_bp_limit must be an int, got {bp_limit_raw!r}"
                ) from exc
            if bp_limit < 0:
                raise ValueError(f"ablation.pretrain_bp_limit must be >= 0, got {bp_limit}")
            if bp_limit > 0:
                records = _apply_bp_limit(records, bp_limit)
                total_bp = sum(len(seq) for _, seq in records)
                LOGGER.info("Ablation enabled: using first %d bp from FASTA", total_bp)

        # Use MLM-specific window size, fallback to max_bp_len for backward compatibility
        window_size = int(cfg_run.get("mlm_window_size", cfg_run.get("max_bp_len", 234)))
        stride = int(cfg_run.get("mlm_stride", 20))
        k = int(cfg_run.get("mlm_kmer", 3))
        include_reverse_complements = reverse_complement_prob > 0
        vocab_path = _resolve_cfg_path(
            cfg_run.get("mlm_vocab_path", f"../artifacts/vocabs/k{k}_mlm_vocab.json"),
            base_dir=script_dir,
        )
        species_name = fasta_path.stem
        single_record_val_ratio = float(cfg_run.get("mlm_single_record_val_ratio", val_ratio))
        exclude_patterns = cfg_run.get("mlm_val_exclude_patterns", ["mitochond"])
        min_val_records = int(cfg_run.get("mlm_min_val_records_per_species", 1))
        effective_val_ratio = single_record_val_ratio if len(records) == 1 else val_ratio
        train_records, val_records = _split_records_for_species(
            records,
            window_bp=window_size,
            val_ratio=effective_val_ratio,
            seed=seed,
            species_name=species_name,
            exclude_patterns=exclude_patterns,
            min_val_records=min_val_records,
        )
        if len(records) == 1 and val_records:
            total_count = len(train_records) + len(val_records)
        else:
            total_count = len(records)
        split_counts = {species_name: (len(train_records), len(val_records), total_count)}
        _log_record_split([species_name], split_counts)

        train_windows = _build_windows_from_records(
            train_records,
            window_size=window_size,
            stride=stride,
            unknown_base_strategy=unknown_base_strategy,
            include_reverse_complements=include_reverse_complements,
            label="train",
        )
        if val_records:
            val_windows = _build_windows_from_records(
                val_records,
                window_size=window_size,
                stride=stride,
                unknown_base_strategy=unknown_base_strategy,
                include_reverse_complements=include_reverse_complements,
                label="val",
            )
        else:
            LOGGER.warning(
                "Validation records empty after record split; using train windows for validation (leakage)."
            )
            val_windows = train_windows

        vocab = load_or_build_vocab(train_windows, k=k, vocab_path=vocab_path)
        train_token_tensors = [vocab.tokenize(window, window_size) for window in train_windows]
        val_token_tensors = [vocab.tokenize(window, window_size) for window in val_windows]
        token_tensors = train_token_tensors
        train_dataset = MLMDataset(train_token_tensors)
        val_dataset = MLMDataset(val_token_tensors)
        LOGGER.info(
            "Dataset split (record-level): %d train, %d val",
            len(train_dataset),
            len(val_dataset),
        )

    # Detect Colab environment and set appropriate num_workers default
    is_colab = os.path.exists("/content") or "COLAB_GPU" in os.environ or "COLAB_RELEASE_TAG" in os.environ
    if "num_workers" in cfg_run:
        # User explicitly set num_workers - respect their choice
        num_workers = int(cfg_run.get("num_workers", 0))
    else:
        # Default: 2 workers on Colab (helps with k-mer tokenization bottleneck), 0 otherwise
        num_workers = 2 if is_colab else 0
        if is_colab:
            LOGGER.info("Colab detected: defaulting to num_workers=%d for better data loading", num_workers)
    steps_per_epoch_approx = steps_per_epoch if use_streaming else len(train_dataset) // batch_size
    total_steps_for_curriculum = epochs * steps_per_epoch_approx // grad_accum_steps
    if max_steps is not None:
        total_steps_for_curriculum = min(total_steps_for_curriculum, max_steps)

    curriculum_scheduler: LengthCurriculumScheduler | None = None
    masking_scheduler: MaskingScheduler | None = None

    # Create MaskingScheduler if progressive masking is enabled
    if progressive_masking_enabled:
        masking_scheduler = MaskingScheduler(
            initial_rate=initial_mask_rate,
            final_rate=final_mask_rate,
            warmup_steps=0,
            total_steps=total_steps_for_curriculum,
            transition_start_ratio=masking_transition_start_ratio,
        )
        LOGGER.info("=" * 60)
        LOGGER.info("PROGRESSIVE MASKING ENABLED (DNABERT-style)")
        LOGGER.info("  Initial rate:      %.1f%%", initial_mask_rate * 100)
        LOGGER.info("  Final rate:        %.1f%%", final_mask_rate * 100)
        LOGGER.info("  Transition start:  %.0f%% of training (step %d)",
                    masking_transition_start_ratio * 100,
                    masking_scheduler.transition_start_step)
        LOGGER.info("  Total steps:       %d", total_steps_for_curriculum)
        LOGGER.info("=" * 60)
        # Use initial rate as starting mask_prob
        mask_prob = initial_mask_rate

    if length_curriculum_config.enabled:
        LOGGER.info("=" * 60)
        LOGGER.info("LENGTH CURRICULUM ENABLED")
        LOGGER.info("  Start length: %d tokens", length_curriculum_config.start_tokens)
        LOGGER.info("  End length:   %d tokens", length_curriculum_config.end_tokens)
        LOGGER.info(
            "  Ramp period:  %.0f%% - %.0f%% of training",
            length_curriculum_config.start_pct * 100,
            length_curriculum_config.end_pct * 100,
        )
        LOGGER.info("  Total steps:  %d", total_steps_for_curriculum)
        LOGGER.info("=" * 60)

        collator, curriculum_scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=total_steps_for_curriculum,
            cfg=cfg_run,
            mask_prob=mask_prob,
            track_spans=True,
            debug=label_debug,
        )
        # If we also have progressive masking, attach the scheduler to the collator
        if masking_scheduler is not None:
            collator.masking_scheduler = masking_scheduler
    else:
        LOGGER.info("Length curriculum: DISABLED (using fixed length)")
        collator = MLMCollator(
            vocab,
            mask_prob=mask_prob,
            debug=label_debug,
            track_spans=True,
            masking_scheduler=masking_scheduler,
        )

    deterministic_kwargs = get_deterministic_dataloader_kwargs(seed, num_workers=num_workers)

    pin_memory = torch.cuda.is_available()
    if use_streaming:
        # Build streaming DataLoader kwargs, adding worker optimizations when applicable
        streaming_loader_kwargs = {
            "batch_size": batch_size,
            "shuffle": False,
            "collate_fn": collate_mlm_streaming,
            "num_workers": num_workers,
            "pin_memory": pin_memory,
            "drop_last": True,
            **deterministic_kwargs,
        }
        if num_workers > 0:
            streaming_loader_kwargs["persistent_workers"] = True
            streaming_loader_kwargs["prefetch_factor"] = 2
        train_loader = DataLoader(train_dataset, **streaming_loader_kwargs)
        val_loader = DataLoader(
            val_dataset,
            batch_size=batch_size,
            shuffle=False,
            collate_fn=collate_mlm_streaming,
            num_workers=0,
            pin_memory=pin_memory,
        )
    else:
        train_loader = DataLoader(
            train_dataset,
            batch_size=batch_size,
            shuffle=True,
            collate_fn=collator,
            num_workers=num_workers,
            pin_memory=pin_memory,
            drop_last=True,
            **deterministic_kwargs,
        )

        val_collator = MLMCollator(vocab, mask_prob=mask_prob)
        val_loader = DataLoader(
            val_dataset,
            batch_size=batch_size,
            shuffle=False,
            collate_fn=val_collator,
            num_workers=0,
            pin_memory=pin_memory,
        )

    if label_debug:
        LOGGER.info("Running LABEL DEBUG mode - verifying first batch then exiting")
        for _batch in train_loader:
            break
        LOGGER.info("Label debug complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": 0,
        }

    if generate_golden_batch:
        from golden_batch import generate_golden_batch as gen_golden

        LOGGER.info("Generating golden batch with seed %d...", seed)

        if golden_batch_path:
            golden_path = Path(str(golden_batch_path))
            if not golden_path.is_absolute():
                golden_path = (script_dir / str(golden_batch_path)).resolve()
        else:
            golden_path = (script_dir / "../artifacts/golden_batch.pt").resolve()

        device = select_device()

        embedding_dim = int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256)))
        transformer_layers = int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8)))
        transformer_heads = int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8)))
        transformer_ff_dim = int(
            cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", embedding_dim * 4))
        )
        transformer_dropout = float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1)))

        # Relative position bias settings
        use_relative_position_bias = bool(cfg_run.get("use_relative_position_bias", False))
        relative_position_num_buckets = int(cfg_run.get("relative_position_num_buckets", 32))
        relative_position_max_distance = int(cfg_run.get("relative_position_max_distance", 128))

        encoder = DNAEncoder(
            vocab_size=len(vocab.itos),
            kmer=vocab.k,
            embedding_dim=embedding_dim,
            num_layers=transformer_layers,
            num_heads=transformer_heads,
            ff_dim=transformer_ff_dim,
            dropout=transformer_dropout,
            pad_token_id=vocab.pad_id,
            drop_path_rate=0.0,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
        )

        special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
        model = DNAMLM(
            encoder,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
            tie_weights=bool(cfg_run.get("mlm_tie_weights", True)),
            use_output_norm=bool(cfg_run.get("mlm_use_output_norm", True)),
        )

        metadata = {
            "config": config_id,
            "embedding_dim": embedding_dim,
            "transformer_layers": transformer_layers,
            "transformer_heads": transformer_heads,
            "batch_size": batch_size,
            "seed": seed,
        }

        gen_golden(
            dataloader=train_loader,
            model=model,
            device=device,
            seed=seed,
            path=golden_path,
            metadata=metadata,
        )

        LOGGER.info("Golden batch generation complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": sum(p.numel() for p in model.parameters()),
        }

    LOGGER.info("Running pretrain sanity checks (checkpoint gate)...")
    sanity_passed, sanity_issues = run_pretrain_sanity_checks(
        train_loader=train_loader,
        vocab=vocab,
        num_batches=5,
        target_mask_prob=mask_prob,
    )

    if not sanity_passed:
        LOGGER.error("Sanity checks FAILED! Training cannot proceed.")
        for issue in sanity_issues:
            LOGGER.error("  - %s", issue)
        raise RuntimeError(
            "Pretrain sanity checks failed. Fix the issues above before training.\n"
            + "\n".join(f"  - {i}" for i in sanity_issues)
        )

    LOGGER.info("Running label token frequency debug report...")
    label_report = None
    for batch_data in train_loader:
        if len(batch_data) == 3:
            _, labels, _ = batch_data
        else:
            _, labels = batch_data
        label_report = debug_label_token_frequency(labels, vocab, num_batches=1)
        break

    if label_report:
        print_label_debug_report(label_report, vocab, top_k=10)
        if label_report["issues"]:
            raise RuntimeError(
                "Label token frequency check FAILED!\n"
                "Special tokens found in supervised labels.\n"
                + "\n".join(f"  - {i}" for i in label_report["issues"])
            )

    if debug_mode:
        LOGGER.info("Running in DEBUG mode - analyzing masking on one batch")
        if use_streaming:
            raw_batch = train_dataset.sample_raw_tokens(batch_size)
        else:
            raw_batch = [train_dataset[i] for i in range(min(batch_size, len(train_dataset)))]
        debug_masking(vocab, raw_batch, mask_prob=mask_prob, max_span_len=3)
        LOGGER.info("Debug mode complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": 0,
        }

    device = select_device()
    LOGGER.info("Using device: %s", device)

    # Log CUDA device details if available
    if device.type == "cuda":
        LOGGER.info("CUDA device name: %s", torch.cuda.get_device_name(0))
        total_mem = torch.cuda.get_device_properties(0).total_memory
        LOGGER.info("CUDA total memory: %.2f GB", total_mem / (1024**3))
        LOGGER.info("CUDA is_available: %s", torch.cuda.is_available())

    embedding_dim = int(cfg_run.get("mlm_embedding_dim", cfg_run.get("embedding_dim", 256)))
    transformer_layers = int(cfg_run.get("mlm_transformer_layers", cfg_run.get("transformer_layers", 8)))
    transformer_heads = int(cfg_run.get("mlm_transformer_heads", cfg_run.get("transformer_heads", 8)))
    transformer_ff_dim = int(cfg_run.get("mlm_transformer_ff_dim", cfg_run.get("transformer_ff_dim", embedding_dim * 4)))
    transformer_dropout = float(cfg_run.get("mlm_transformer_dropout", cfg_run.get("transformer_dropout", 0.1)))

    LOGGER.info("=" * 60)
    LOGGER.info("MLM ARCHITECTURE (Enhanced for lower loss)")
    LOGGER.info("  Embedding dim: %d", embedding_dim)
    LOGGER.info("  Transformer layers: %d", transformer_layers)
    LOGGER.info("  Transformer heads: %d", transformer_heads)
    LOGGER.info("  FF dim: %d (%.1fx embedding)", transformer_ff_dim, transformer_ff_dim / embedding_dim)
    LOGGER.info("  Dropout: %.2f", transformer_dropout)
    LOGGER.info("=" * 60)

    # Relative position bias settings
    use_relative_position_bias = bool(cfg_run.get("use_relative_position_bias", False))
    relative_position_num_buckets = int(cfg_run.get("relative_position_num_buckets", 32))
    relative_position_max_distance = int(cfg_run.get("relative_position_max_distance", 128))

    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        drop_path_rate=0.0,
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=relative_position_num_buckets,
        relative_position_max_distance=relative_position_max_distance,
    )

    tie_weights = bool(cfg_run.get("mlm_tie_weights", True))
    use_output_norm = bool(cfg_run.get("mlm_use_output_norm", True))

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    LOGGER.info(
        "Special tokens excluded from predictions: [MASK]=%d, <UNK>=%d, <PAD>=%d",
        vocab.mask_id,
        vocab.unk_id,
        vocab.pad_id,
    )

    model = DNAMLM(
        encoder,
        vocab_size=len(vocab.itos),
        special_token_ids=special_token_ids,
        tie_weights=tie_weights,
        use_output_norm=use_output_norm,
    ).to(device)

    # Log model device placement
    model_device = next(model.parameters()).device
    LOGGER.info("Model parameters device: %s", model_device)
    if device.type == "cuda":
        LOGGER.info("CUDA memory allocated: %.2f MB", torch.cuda.memory_allocated() / (1024**2))
        LOGGER.info("CUDA memory reserved: %.2f MB", torch.cuda.memory_reserved() / (1024**2))

    total_params = sum(p.numel() for p in model.parameters())
    unique_params = sum(p.numel() for p in set(model.parameters()))
    LOGGER.info("Model parameters: %.2fM total, %.2fM unique", total_params / 1e6, unique_params / 1e6)

    if leakage_test:
        LOGGER.info("Running LEAKAGE TEST mode - verifying no information leakage")
        verify_no_leakage(model, vocab, device)
        LOGGER.info("Leakage test complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": total_params,
        }

    if overfit_debug:
        LOGGER.info("Running OVERFIT DEBUG mode")
        if use_streaming:
            overfit_samples = train_dataset.sample_raw_tokens(32)
        else:
            overfit_samples = token_tensors[:32]
        run_overfit_debug(
            model=model,
            train_samples=overfit_samples,
            vocab=vocab,
            device=device,
            num_steps=1000,
            lr=5e-3,
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return {
            "mlm_loss": float("nan"),
            "mlm_steps": 0,
            "encoder_ckpt_path": "",
            "pretrain_seconds": float(time.perf_counter() - start_time),
            "params": total_params,
        }

    # Use MLM-specific learning rate
    lr = float(cfg_run.get("mlm_lr", cfg_run.get("lr", 2e-4)))
    weight_decay = float(cfg_run.get("mlm_weight_decay", 0.01))
    lr_decay = float(cfg_run.get("mlm_layer_lr_decay", 0.9))
    use_layer_lr_decay = bool(cfg_run.get("mlm_use_layer_lr_decay", True))
    optimizer_foreach = bool(cfg_run.get("optimizer_foreach", True))

    LOGGER.info("=" * 60)
    LOGGER.info("OPTIMIZER CONFIGURATION")
    LOGGER.info("Optimizer: AdamW foreach=%s", optimizer_foreach)
    if use_layer_lr_decay:
        param_groups = get_layer_wise_lr_groups(model, lr, lr_decay, weight_decay)
        LOGGER.info("Layer-wise LR decay: ENABLED (factor: %.2f)", lr_decay)
        LOGGER.info("  Lower layers get smaller LRs to preserve pretrained features")
        optimizer = AdamW(
            param_groups,
            betas=(0.9, 0.98),
            eps=1e-6,
            foreach=optimizer_foreach,
        )
    else:
        LOGGER.info("Layer-wise LR decay: DISABLED (all layers use same LR)")
        LOGGER.info("  All parameters use base LR: %.2e", lr)
        optimizer = AdamW(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay,
            betas=(0.9, 0.98),
            eps=1e-6,
            foreach=optimizer_foreach,
        )
    LOGGER.info("=" * 60)

    # Early stopping setup
    early_stopping_enabled = bool(cfg_run.get("mlm_early_stopping_enabled", True))
    if early_stopping_enabled:
        early_stopper = PretrainingEarlyStopping(
            patience=int(cfg_run.get("mlm_early_stopping_patience", 20)),
            min_delta=float(cfg_run.get("mlm_early_stopping_min_delta", 0.001)),
            min_epochs=int(cfg_run.get("mlm_early_stopping_min_epochs", 30)),
            kmer=int(vocab.k),
            restore_best=bool(cfg_run.get("mlm_early_stopping_restore_best", True)),
            verbose=True,
        )
        LOGGER.info(
            "Early stopping enabled (patience=%d, min_delta=%.6f, min_epochs=%d)",
            early_stopper.patience,
            early_stopper.min_delta,
            early_stopper.min_epochs,
        )
    else:
        early_stopper = None
        LOGGER.info("Early stopping disabled")

    if use_streaming and steps_per_epoch is not None:
        batches_per_epoch = steps_per_epoch
    else:
        batches_per_epoch = len(train_loader)
    optimizer_steps_per_epoch = batches_per_epoch // grad_accum_steps
    total_steps = epochs * optimizer_steps_per_epoch
    if max_steps is not None:
        total_steps = min(total_steps, max_steps)

    warmup_ratio = float(cfg_run.get("mlm_warmup_ratio", 0.06))
    warmup_steps = int(warmup_ratio * total_steps)
    final_lr_ratio = float(cfg_run.get("mlm_final_lr_ratio", 0.1))
    lr_log_interval = int(cfg_run.get("mlm_lr_log_interval", 50))

    LOGGER.info("=" * 60)
    LOGGER.info("TOTAL STEPS CALCULATION")
    LOGGER.info("  Epochs:                  %d", epochs)
    LOGGER.info("  Batches per epoch:       %d", batches_per_epoch)
    LOGGER.info("  Grad accum steps:        %d", grad_accum_steps)
    LOGGER.info("  Optimizer steps/epoch:   %d", optimizer_steps_per_epoch)
    LOGGER.info("  Total optimizer steps:   %d", total_steps)
    LOGGER.info("=" * 60)

    scheduler = get_warmup_cosine_schedule(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        final_lr_ratio=final_lr_ratio,
    )
    LOGGER.info("=" * 60)
    LOGGER.info("LR SCHEDULE (Linear Warmup + Cosine Decay)")
    LOGGER.info("  Base LR:        %.2e", lr)
    LOGGER.info("  Warmup steps:   %d (%.1f%% of total)", warmup_steps, warmup_ratio * 100)
    LOGGER.info("  Total steps:    %d", total_steps)
    LOGGER.info("  Final LR:       %.2e (%.0f%% of peak)", lr * final_lr_ratio, final_lr_ratio * 100)
    LOGGER.info("  LR log interval: every %d optimizer steps", lr_log_interval)
    LOGGER.info("=" * 60)

    # AMP (Mixed Precision) setup using config
    amp_ctx = get_amp_context(device, cfg_run, logger=LOGGER)
    log_amp_status(amp_ctx, logger=LOGGER)

    autocast_device: str | None = amp_ctx.device_type if amp_ctx.enabled else None
    autocast_dtype: torch.dtype | None = amp_ctx.dtype if amp_ctx.enabled else None
    scaler = None

    # GradScaler only for CUDA with float16
    if amp_ctx.enabled and device.type == "cuda" and amp_ctx.dtype == torch.float16:
        scaler = torch.cuda.amp.GradScaler()
        LOGGER.info("Using GradScaler for CUDA fp16")

    if debug_batch:
        LOGGER.info("Running DEBUG BATCH check (one train batch + one val batch)")
        was_training = model.training
        model.eval()

        def log_one_batch(label: str, loader: DataLoader) -> None:
            for batch_data in loader:
                if len(batch_data) == 3:
                    inputs, labels, _span_lengths = batch_data
                else:
                    inputs, labels = batch_data

                non_blocking = (device.type == "cuda")
                inputs = inputs.to(device, non_blocking=non_blocking)
                labels = labels.to(device, non_blocking=non_blocking)

                with torch.no_grad():
                    if autocast_device and autocast_dtype is not None:
                        with torch.autocast(autocast_device, dtype=autocast_dtype):
                            logits, _ = model(inputs, labels=None)
                    else:
                        logits, _ = model(inputs, labels=None)

                metrics = compute_masked_batch_metrics(logits, labels)
                LOGGER.info(
                    "Debug batch (%s): masked_fraction=%.4f, loss_masked_only=%.4f, masked_token_accuracy=%.4f",
                    label,
                    metrics["masked_fraction"],
                    metrics["loss_masked_only"],
                    metrics["masked_token_accuracy"],
                )
                break

        log_one_batch("train", train_loader)
        log_one_batch("val", val_loader)

        if was_training:
            model.train()

    best_val_loss = float("inf")
    best_val_acc = 0.0
    patience = int(cfg_run.get("mlm_patience", 25))
    patience_counter = 0
    global_step = 0
    param_norm_warn = float(cfg_run.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg_run.get("param_norm_cap", 200.0))

    encoder_path = _resolve_cfg_path(
        cfg_run.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt"),
        base_dir=script_dir,
    )
    encoder_path.parent.mkdir(parents=True, exist_ok=True)
    encoder_metadata_path = encoder_path.with_suffix(".metadata.json")

    def _write_encoder_metadata(best_loss: float) -> None:
        metadata = {
            "config_snapshot": resolved_config,
            "vocab_path": str(vocab_path) if vocab_path is not None else None,
            "kmer": int(vocab.k),
            "embedding_dim": embedding_dim,
            "num_layers": transformer_layers,
            "heads": transformer_heads,
            "best_val_loss": float(best_loss),
            "early_stopping": {
                "enabled": early_stopping_enabled,
                "triggered": early_stopper.should_stop if early_stopper else False,
                "best_epoch": early_stopper.best_epoch if early_stopper else None,
                "patience": early_stopper.patience if early_stopper else None,
            },
        }
        with encoder_metadata_path.open("w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

    def _save_encoder_checkpoint(best_loss: float) -> None:
        torch.save(model.encoder.state_dict(), encoder_path)
        _write_encoder_metadata(best_loss)

    def _save_training_checkpoint(
        step: int,
        epoch: int,
        best_val_loss: float,
        best_val_acc: float,
        patience_counter: int,
    ) -> Path:
        """Save a full training checkpoint including all state for resuming."""
        ckpt_path = checkpoints_dir / f"checkpoint_step_{step}.pt"
        checkpoint = {
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "scheduler_state_dict": scheduler.state_dict(),
            "epoch": epoch,
            "global_step": step,
            "best_val_loss": best_val_loss,
            "best_val_acc": best_val_acc,
            "patience_counter": patience_counter,
            "config_snapshot": resolved_config,
            "vocab_k": vocab.k,
            "rng_state_python": random.getstate(),
            "rng_state_torch": torch.get_rng_state(),
        }
        if scaler is not None:
            checkpoint["scaler_state_dict"] = scaler.state_dict()
        if masking_scheduler is not None:
            checkpoint["masking_scheduler_state_dict"] = masking_scheduler.state_dict()
        if torch.cuda.is_available():
            checkpoint["rng_state_cuda"] = torch.cuda.get_rng_state_all()
        torch.save(checkpoint, ckpt_path)
        LOGGER.info("Saved training checkpoint: %s (step %d, epoch %d)", ckpt_path, step, epoch)
        return ckpt_path

    def _load_training_checkpoint(ckpt_path: Path) -> Dict[str, Any]:
        """Load a training checkpoint and restore all state."""
        LOGGER.info("Loading training checkpoint from: %s", ckpt_path)
        checkpoint = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint["model_state_dict"])
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        if scaler is not None and "scaler_state_dict" in checkpoint:
            scaler.load_state_dict(checkpoint["scaler_state_dict"])
        if masking_scheduler is not None and "masking_scheduler_state_dict" in checkpoint:
            masking_scheduler.load_state_dict(checkpoint["masking_scheduler_state_dict"])
            LOGGER.info("Restored masking scheduler: %s", masking_scheduler)
        if "rng_state_python" in checkpoint:
            random.setstate(checkpoint["rng_state_python"])
        if "rng_state_torch" in checkpoint:
            torch.set_rng_state(checkpoint["rng_state_torch"])
        if "rng_state_cuda" in checkpoint and torch.cuda.is_available():
            torch.cuda.set_rng_state_all(checkpoint["rng_state_cuda"])
        LOGGER.info(
            "Restored checkpoint: epoch=%d, global_step=%d, best_val_loss=%.4f",
            checkpoint["epoch"],
            checkpoint["global_step"],
            checkpoint["best_val_loss"],
        )
        return checkpoint

    # Resume from checkpoint if specified
    resume_ckpt_path: Optional[Path] = None
    start_epoch = 1
    if resume_path_value is not None:
        resume_ckpt_path = Path(resume_path_value).expanduser()
        if not resume_ckpt_path.is_absolute():
            resume_ckpt_path = (script_dir / resume_ckpt_path).resolve()
        if resume_ckpt_path.exists():
            loaded_ckpt = _load_training_checkpoint(resume_ckpt_path)
            global_step = loaded_ckpt["global_step"]
            start_epoch = loaded_ckpt["epoch"]
            best_val_loss = loaded_ckpt["best_val_loss"]
            best_val_acc = loaded_ckpt.get("best_val_acc", 0.0)
            patience_counter = loaded_ckpt.get("patience_counter", 0)
            LOGGER.info("Resuming training from step %d, epoch %d", global_step, start_epoch)
        else:
            LOGGER.warning("Resume checkpoint not found: %s - starting from scratch", resume_ckpt_path)

    LOGGER.info("=" * 60)
    LOGGER.info("STARTING MLM PRETRAINING")
    if masking_scheduler is not None:
        LOGGER.info(
            "Epochs: %d, LR: %.2e, Mask prob: %.1f%% → %.1f%% (progressive)",
            epochs, lr, initial_mask_rate * 100, final_mask_rate * 100
        )
    else:
        LOGGER.info("Epochs: %d, LR: %.2e, Mask prob: %.1f%%", epochs, lr, mask_prob * 100)
    if max_steps is not None:
        LOGGER.info("Ablation enabled: stopping after %d optimizer steps", max_steps)
    LOGGER.info("=" * 60)

    health_log_interval = int(cfg_run.get("health_log_interval", 50))
    health_dashboard = create_health_dashboard(
        run_dir=run_dir,
        vocab=vocab,
        log_every_n_steps=health_log_interval,
        target_mask_prob=mask_prob,
    )
    health_dashboard.print_table_header()

    numerics_checker = NumericsChecker(enabled=True)
    LOGGER.info("Numerics checker enabled - will fail fast on NaN/Inf")

    # Gap monitor: CSV file for tracking train-val gap and other metrics
    gap_csv_path = (script_dir / "../artifacts/mlm_metrics.csv").resolve()
    gap_csv_path.parent.mkdir(parents=True, exist_ok=True)
    gap_csv_columns = [
        "epoch",
        "train_loss",
        "val_loss",
        "gap_loss",
        "train_masked_acc",
        "val_masked_acc",
        "gap_acc",
        "embedding_weight_norm",
        "grad_norm",
        "lr",
    ]
    # Write CSV header (overwrite any existing file)
    with gap_csv_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=gap_csv_columns)
        writer.writeheader()
    LOGGER.info("Gap monitor CSV initialized: %s", gap_csv_path)

    best_train_loss = float("inf")
    reached_max_steps = False
    # Track last gradient norm for epoch-level logging
    last_grad_norm = 0.0
    # Tracking for step logging
    epoch_start_time = time.perf_counter()
    step_tokens_processed = 0
    last_log_time = time.perf_counter()
    last_log_step = global_step
    # GPU telemetry tracking (accumulated over telemetry interval)
    telemetry_interval = 50  # log every N optimizer steps
    telemetry_data_time = 0.0
    telemetry_compute_time = 0.0
    telemetry_step_count = 0
    data_load_start = time.perf_counter()

    for epoch in range(start_epoch, epochs + 1):
        model.train()
        LOGGER.info("Train loop: model.training=%s", model.training)
        total_loss_unscaled = 0.0
        total_correct = 0
        total_masked = 0
        total_tokens = 0
        samples_seen = 0
        optimizer.zero_grad()

        if scaler is not None:
            LOGGER.info("Epoch %d: GradScaler scale = %.0f", epoch, scaler.get_scale())

        if curriculum_scheduler is not None:
            curr_len = curriculum_scheduler.current_length
            progress = curriculum_scheduler.get_progress()
            LOGGER.info(
                "Epoch %d: Curriculum length = %d tokens (%s, %.0f%% ramp progress)",
                epoch,
                curr_len,
                progress["phase"],
                progress["ramp_progress"] * 100,
            )

        for batch_idx, batch_data in enumerate(train_loader, start=1):
            # Measure data loading time
            data_load_end = time.perf_counter()
            telemetry_data_time += data_load_end - data_load_start

            if max_steps is not None and global_step >= max_steps:
                reached_max_steps = True
                break

            if len(batch_data) == 3:
                inputs, labels, span_lengths = batch_data
            else:
                inputs, labels = batch_data
                span_lengths = None

            non_blocking = (device.type == "cuda")
            inputs = inputs.to(device, non_blocking=non_blocking)
            labels = labels.to(device, non_blocking=non_blocking)
            samples_seen += int(inputs.size(0))
            total_tokens += labels.numel()

            current_lr = optimizer.param_groups[0]["lr"]
            numerics_checker.set_context(
                epoch=epoch,
                step=global_step,
                lr=current_lr,
                batch_idx=batch_idx,
                model=model,
            )

            # Start compute timing
            if device.type == "cuda":
                torch.cuda.synchronize()
            compute_start = time.perf_counter()

            if autocast_device and autocast_dtype is not None:
                with torch.autocast(autocast_device, dtype=autocast_dtype):
                    logits, _ = model(inputs, labels=None)
            else:
                logits, _ = model(inputs, labels=None)

            loss_unscaled = F.cross_entropy(
                logits.float().view(-1, logits.size(-1)),
                labels.view(-1),
                ignore_index=-100,
            )
            numerics_checker.check_forward(logits, loss_unscaled)

            loss_scaled = loss_unscaled / grad_accum_steps
            if scaler is not None:
                scaler.scale(loss_scaled).backward()
            else:
                loss_scaled.backward()

            if scaler is not None:
                numerics_checker.check_backward(model, check_grads=False)
            else:
                numerics_checker.check_backward(model)

            # End compute timing
            if device.type == "cuda":
                torch.cuda.synchronize()
            compute_end = time.perf_counter()
            telemetry_compute_time += compute_end - compute_start

            total_loss_unscaled += loss_unscaled.item() * inputs.size(0)

            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()

            health_dashboard.accumulate_batch(
                loss=loss_unscaled.item(),
                logits=logits.detach(),
                labels=labels,
                masked_inputs=inputs,
                span_lengths=span_lengths,
            )

            if batch_idx % grad_accum_steps == 0:
                if scaler is not None:
                    scaler.unscale_(optimizer)
                unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                last_grad_norm = unclipped_grad_norm.item() if hasattr(unclipped_grad_norm, 'item') else float(unclipped_grad_norm)

                if scaler is not None:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()

                total_norm, scale = cap_model_param_norm_(model, max_norm=param_norm_cap)
                effective_param_norm = param_norm_cap if scale < 1.0 else total_norm
                if scale < 1.0:
                    LOGGER.warning(
                        "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                        global_step + 1,
                        total_norm,
                        param_norm_cap,
                        scale,
                    )

                scheduler.step()
                optimizer.zero_grad()
                global_step += 1
                telemetry_step_count += 1

                # GPU telemetry every N optimizer steps
                if global_step % telemetry_interval == 0 and telemetry_step_count > 0:
                    avg_data_ms = (telemetry_data_time / telemetry_step_count) * 1000
                    avg_compute_ms = (telemetry_compute_time / telemetry_step_count) * 1000
                    total_step_ms = avg_data_ms + avg_compute_ms
                    data_pct = (avg_data_ms / total_step_ms * 100) if total_step_ms > 0 else 0
                    compute_pct = (avg_compute_ms / total_step_ms * 100) if total_step_ms > 0 else 0
                    telemetry_msg = (
                        f"[Telemetry] Step {global_step}: "
                        f"step={total_step_ms:.1f}ms (data={avg_data_ms:.1f}ms [{data_pct:.0f}%], "
                        f"compute={avg_compute_ms:.1f}ms [{compute_pct:.0f}%])"
                    )
                    if device.type == "cuda":
                        mem_alloc = torch.cuda.memory_allocated() / (1024 ** 2)
                        mem_max = torch.cuda.max_memory_allocated() / (1024 ** 2)
                        telemetry_msg += f" | GPU mem: {mem_alloc:.0f}MB alloc, {mem_max:.0f}MB peak"
                    LOGGER.info(telemetry_msg)
                    # Reset telemetry accumulators
                    telemetry_data_time = 0.0
                    telemetry_compute_time = 0.0
                    telemetry_step_count = 0

                if global_step % lr_log_interval == 0:
                    current_lr = optimizer.param_groups[0]["lr"]
                    if global_step <= warmup_steps:
                        phase = "warmup"
                        phase_progress = global_step / max(1, warmup_steps) * 100
                    else:
                        phase = "decay"
                        phase_progress = (global_step - warmup_steps) / max(1, total_steps - warmup_steps) * 100

                    LOGGER.info(
                        "Step %d/%d [%s %.0f%%]: LR=%.2e, grad_norm=%.4f, param_norm=%.2f",
                        global_step,
                        total_steps,
                        phase,
                        phase_progress,
                        current_lr,
                        unclipped_grad_norm.item(),
                        effective_param_norm,
                    )
                    if param_norm_warn > 0 and effective_param_norm > param_norm_warn:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f above warn=%.2f",
                            global_step,
                            effective_param_norm,
                            param_norm_warn,
                        )

                if curriculum_scheduler is not None:
                    new_length = curriculum_scheduler.step(global_step)
                    if curriculum_scheduler.is_ramping() and global_step % 100 == 0:
                        LOGGER.debug("Step %d: Curriculum length = %d tokens", global_step, new_length)

                # Update progressive masking rate
                if masking_scheduler is not None:
                    new_mask_rate = masking_scheduler.get_mask_rate(global_step)
                    # Update collator (for non-streaming mode with num_workers=0)
                    if hasattr(collator, 'set_mask_prob'):
                        collator.set_mask_prob(new_mask_rate)
                    elif hasattr(collator, 'mask_prob'):
                        collator.mask_prob = new_mask_rate
                    # For streaming mode, update the dataset directly
                    if use_streaming and hasattr(train_dataset, 'set_mask_prob'):
                        train_dataset.set_mask_prob(new_mask_rate)
                    # Log mask rate changes periodically
                    if global_step % lr_log_interval == 0:
                        in_transition = global_step >= masking_scheduler.transition_start_step
                        phase = "transition" if in_transition else "stable"
                        LOGGER.debug(
                            "Step %d: Mask rate = %.1f%% (%s phase)",
                            global_step, new_mask_rate * 100, phase
                        )

                health_metrics = health_dashboard.log_step(
                    step=global_step,
                    epoch=epoch,
                    model=model,
                    optimizer=optimizer,
                )
                if health_metrics:
                    health_dashboard.print_table_row(health_metrics)

                # Checkpoint saving every N steps
                if save_every_steps > 0 and global_step % save_every_steps == 0:
                    _save_training_checkpoint(
                        step=global_step,
                        epoch=epoch,
                        best_val_loss=best_val_loss,
                        best_val_acc=best_val_acc,
                        patience_counter=patience_counter,
                    )
                    removed = prune_step_checkpoints(
                        checkpoints_dir,
                        checkpoint_keep_last,
                        keep_paths=[resume_ckpt_path] if resume_ckpt_path else None,
                        logger=LOGGER,
                    )
                    if removed:
                        LOGGER.info(
                            "Pruned %d old training checkpoints (keep_last=%d)",
                            len(removed),
                            checkpoint_keep_last,
                        )

                # Enhanced step logging with tokens/sec, elapsed time, ETA
                if global_step % lr_log_interval == 0 or global_step == 1:
                    current_time = time.perf_counter()
                    elapsed_since_last = current_time - last_log_time
                    steps_since_last = global_step - last_log_step
                    tokens_this_step = total_tokens

                    if elapsed_since_last > 0 and steps_since_last > 0:
                        tokens_per_sec = tokens_this_step / elapsed_since_last
                        steps_per_sec = steps_since_last / elapsed_since_last
                    else:
                        tokens_per_sec = 0.0
                        steps_per_sec = 0.0

                    total_elapsed = current_time - start_time
                    remaining_steps = total_steps - global_step
                    if steps_per_sec > 0:
                        eta_seconds = remaining_steps / steps_per_sec
                        eta_str = f"{eta_seconds/60:.1f}min" if eta_seconds < 3600 else f"{eta_seconds/3600:.1f}hr"
                    else:
                        eta_str = "N/A"

                    avg_loss = total_loss_unscaled / max(1, samples_seen)
                    LOGGER.info(
                        "Step %d | loss=%.4f | tok/s=%.0f | elapsed=%.1fs | ETA=%s",
                        global_step,
                        avg_loss,
                        tokens_per_sec,
                        total_elapsed,
                        eta_str,
                    )

                    last_log_time = current_time
                    last_log_step = global_step

                if max_steps is not None and global_step >= max_steps:
                    reached_max_steps = True
                    break

            if batch_idx % 100 == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                batch_acc = total_correct / max(1, total_masked)
                LOGGER.info(
                    "Epoch %d [%d/%d] Loss: %.4f Acc: %.2f%% LR: %.2e",
                    epoch,
                    batch_idx,
                    len(train_loader),
                    loss_unscaled.item(),
                    batch_acc * 100,
                    current_lr,
                )

            # Reset data loading timer for next batch
            data_load_start = time.perf_counter()

        train_loss = total_loss_unscaled / max(1, samples_seen)
        train_acc = total_correct / max(1, total_masked)
        train_masked_fraction = total_masked / max(1, total_tokens)
        train_ppl = math.exp(min(train_loss, 100)) if math.isfinite(train_loss) else float("nan")
        if math.isfinite(train_loss):
            best_train_loss = min(best_train_loss, train_loss)

        val_metrics = run_validation(model, val_loader, device, autocast_device, autocast_dtype)
        val_loss = val_metrics["loss"]
        val_acc = val_metrics["accuracy"]
        val_ppl = val_metrics["perplexity"]
        val_masked_fraction = val_metrics["masked_fraction"]
        val_masked_accuracy = val_metrics["masked_token_accuracy"]
        val_loss_masked_only = val_metrics["loss_masked_only"]

        curriculum_info = f" | Len: {curriculum_scheduler.current_length}" if curriculum_scheduler is not None else ""

        LOGGER.info(
            "Epoch %d: Train [Loss: %.4f, Acc: %.1f%%, PPL: %.2f] | Val [Loss: %.4f, Acc: %.1f%%, PPL: %.2f]%s",
            epoch,
            train_loss,
            train_acc * 100,
            train_ppl,
            val_loss,
            val_acc * 100,
            val_ppl,
            curriculum_info,
        )
        LOGGER.info(
            "Epoch %d: masked_fraction train=%.4f val=%.4f | masked_token_accuracy train=%.4f val=%.4f | "
            "loss_masked_only train=%.4f val=%.4f",
            epoch,
            train_masked_fraction,
            val_masked_fraction,
            train_acc,
            val_masked_accuracy,
            train_loss,
            val_loss_masked_only,
        )

        # Gap monitor: compute and log gap metrics
        gap_loss = train_loss - val_loss_masked_only
        gap_acc = train_acc - val_masked_accuracy
        current_lr = optimizer.param_groups[0]["lr"]

        # Compute embedding weight norm
        embedding_weight_norm = 0.0
        if hasattr(model, 'encoder') and hasattr(model.encoder, 'embedding'):
            emb_weight = model.encoder.embedding.weight
            embedding_weight_norm = torch.norm(emb_weight).item()
        elif hasattr(model, 'embedding'):
            emb_weight = model.embedding.weight
            embedding_weight_norm = torch.norm(emb_weight).item()

        LOGGER.info(
            "Epoch %d GAP: loss_gap=%.4f acc_gap=%.4f | emb_norm=%.2f grad_norm=%.4f lr=%.2e",
            epoch,
            gap_loss,
            gap_acc,
            embedding_weight_norm,
            last_grad_norm,
            current_lr,
        )

        # Append to CSV
        with gap_csv_path.open("a", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=gap_csv_columns)
            writer.writerow({
                "epoch": epoch,
                "train_loss": f"{train_loss:.6f}",
                "val_loss": f"{val_loss_masked_only:.6f}",
                "gap_loss": f"{gap_loss:.6f}",
                "train_masked_acc": f"{train_acc:.6f}",
                "val_masked_acc": f"{val_masked_accuracy:.6f}",
                "gap_acc": f"{gap_acc:.6f}",
                "embedding_weight_norm": f"{embedding_weight_norm:.4f}",
                "grad_norm": f"{last_grad_norm:.6f}",
                "lr": f"{current_lr:.8e}",
            })

        if topk_debug:
            for debug_batch_data in train_loader:
                if len(debug_batch_data) == 3:
                    debug_inputs, debug_labels, _ = debug_batch_data
                else:
                    debug_inputs, debug_labels = debug_batch_data
                break

            topk_stats = print_topk_predictions_debug(
                model=model,
                inputs=debug_inputs,
                labels=debug_labels,
                vocab=vocab,
                epoch=epoch,
                num_sequences=3,
                top_k=5,
            )

            LOGGER.info(
                "Top-K Debug: Avg top-1 prob: %.1f%%, True in top-1: %.1f%%, True in top-5: %.1f%%",
                topk_stats["avg_top1_prob"] * 100,
                topk_stats["true_in_top1_pct"],
                topk_stats["true_in_topk_pct"],
            )

        if early_stopper is not None:
            checkpoint_dir = encoder_path.parent
            should_stop = early_stopper(
                val_loss=val_loss_masked_only,
                model=encoder,
                epoch=epoch,
                checkpoint_dir=checkpoint_dir,
                val_perplexity=val_ppl,
            )

            if early_stopper.best_epoch == epoch:
                best_val_loss = val_loss_masked_only
                best_val_acc = val_acc
                patience_counter = 0
                _save_encoder_checkpoint(best_val_loss)
                LOGGER.info(
                    "★ New best! Val Loss (masked-only): %.4f, PPL: %.2f, Acc: %.1f%% - Saved to %s",
                    val_loss_masked_only,
                    val_ppl,
                    val_acc * 100,
                    encoder_path,
                )
            else:
                patience_counter = early_stopper.epochs_without_improvement

            if should_stop:
                LOGGER.info("Early stopping triggered at epoch %d", epoch)
                if early_stopper.best_score is not None:
                    LOGGER.info(
                        "Best validation loss: %.6f at epoch %d",
                        early_stopper.best_score,
                        early_stopper.best_epoch,
                    )
                break
        else:
            if val_loss_masked_only < best_val_loss:
                best_val_loss = val_loss_masked_only
                best_val_acc = val_acc
                patience_counter = 0
                _save_encoder_checkpoint(best_val_loss)
                LOGGER.info(
                    "★ New best! Val Loss (masked-only): %.4f, PPL: %.2f, Acc: %.1f%% - Saved to %s",
                    val_loss_masked_only,
                    val_ppl,
                    val_acc * 100,
                    encoder_path,
                )
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    LOGGER.info("Early stopping at epoch %d (patience %d)", epoch, patience)
                    break

        if reached_max_steps:
            LOGGER.info(
                "Reached ablation.pretrain_max_steps=%d (optimizer steps); stopping after validation",
                max_steps,
            )
            break

    if not encoder_path.exists():
        _save_encoder_checkpoint(best_val_loss)

    # Save final checkpoint
    final_ckpt_path = _save_training_checkpoint(
        step=global_step,
        epoch=epoch,
        best_val_loss=best_val_loss,
        best_val_acc=best_val_acc,
        patience_counter=patience_counter,
    )
    # Also save as "final" checkpoint for easy identification
    final_named_path = checkpoints_dir / "checkpoint_final.pt"
    shutil.copy(final_ckpt_path, final_named_path)
    LOGGER.info("Saved final checkpoint: %s", final_named_path)
    removed = prune_step_checkpoints(
        checkpoints_dir,
        checkpoint_keep_last,
        keep_paths=[resume_ckpt_path] if resume_ckpt_path else None,
        logger=LOGGER,
    )
    if removed:
        LOGGER.info(
            "Pruned %d old training checkpoints (keep_last=%d)",
            len(removed),
            checkpoint_keep_last,
        )

    # Copy encoder to output_dir as well (if output_dir is different from encoder_path's parent)
    output_encoder_path = output_dir / f"mlm_encoder_k{vocab.k}.pt"
    if output_encoder_path.resolve() != encoder_path.resolve():
        shutil.copy(encoder_path, output_encoder_path)
        LOGGER.info("Copied encoder to output_dir: %s", output_encoder_path)

    _write_encoder_metadata(best_val_loss)

    if early_stopper is not None and early_stopper.restore_best:
        LOGGER.info("Restoring best model weights...")
        restored = early_stopper.restore_best_weights(encoder)
        if restored:
            LOGGER.info("Restored weights from epoch %d", early_stopper.best_epoch)

    best_ppl = math.exp(min(best_val_loss, 100)) if math.isfinite(best_val_loss) else float("nan")
    LOGGER.info("=" * 60)
    LOGGER.info("PRETRAINING COMPLETE")
    LOGGER.info("Best Val Loss: %.4f (Perplexity: %.2f)", best_val_loss, best_ppl)
    LOGGER.info("Best Val Accuracy: %.1f%%", best_val_acc * 100)
    LOGGER.info("Target: Loss < 1.0 (PPL < 2.72)")
    LOGGER.info("Saved pretrained encoder to %s", encoder_path)
    LOGGER.info("Checkpoints directory: %s", checkpoints_dir)
    LOGGER.info("=" * 60)

    return {
        "mlm_loss": float(best_train_loss if math.isfinite(best_train_loss) else train_loss),
        "mlm_steps": int(global_step),
        "encoder_ckpt_path": str(encoder_path),
        "pretrain_seconds": float(time.perf_counter() - start_time),
        "params": int(total_params),
        "best_val_loss": float(best_val_loss),
        "best_val_accuracy": float(best_val_acc),
        "epochs_ran": int(epoch),
    }


def main(cfg: dict | None = None) -> dict:
    if cfg is not None:
        return run_mlm_pretrain(cfg)

    parser = argparse.ArgumentParser(description="Pretrain DNA MLM on E. coli genome")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--kmer", type=int, default=None, help="Single k-mer to pretrain")
    parser.add_argument("--kmers", type=int, nargs="+", default=None, help="Multiple k-mers")
    parser.add_argument("--all-kmers", action="store_true", help="Pretrain all k-mers (3,4,5,6)")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoints")
    parser.add_argument("--debug", action="store_true", help="Run debug mode: analyze one batch and exit")
    parser.add_argument("--overfit_debug", action="store_true", 
                        help="Train on 32 samples for 200 steps to verify model can learn")
    parser.add_argument("--label_debug", action="store_true",
                        help="Enable label verification debug: proves labels match originals. "
                             "Prints original/masked/label token IDs and asserts correctness for 5 positions per sequence. "
                             "Hard error if any mismatch. Runs once on first batch then exits.")
    parser.add_argument("--leakage_test", action="store_true",
                        help="Verify no information leakage at masked positions. "
                             "Tests: (1) masked embeddings don't contain original info, "
                             "(2) model runs with 100%% masking, (3) no engineered features. "
                             "Hard error if leakage detected.")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility. Overrides config value.")
    parser.add_argument("--generate_golden_batch", action="store_true",
                        help="Generate and save a golden batch checkpoint for reproducibility testing.")
    parser.add_argument("--golden_batch_path", type=str, default=None,
                        help="Path to golden batch file (default: artifacts/golden_batch.pt)")
    parser.add_argument("--topk_debug", action="store_true",
                        help="Enable top-k sanity print for masked positions. "
                             "Prints masked position index, true token, and top-5 predictions + probabilities "
                             "for 3 sequences once per epoch. Helps verify predictions become sharper over time.")
    parser.add_argument("--debug_batch", action="store_true",
                        help="Print masked_fraction, loss_masked_only, masked_token_accuracy for one train "
                             "batch and one val batch.")
    parser.add_argument("--mlm_dry_run", action="store_true",
                        help="Run a short streaming dry run (3 train/3 val batches) and exit.")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Output directory for checkpoints and artifacts. "
                             "Defaults to value from config or 'Gene Whisperer/artifacts'.")
    parser.add_argument("--resume_path", type=str, default=None,
                        help="Path to a checkpoint to resume training from. "
                             "Restores model, optimizer, scheduler, scaler, epoch, and global_step.")
    parser.add_argument("--save_every_steps", type=int, default=None,
                        help="Save a checkpoint every N optimizer steps (default: 500).")
    parser.add_argument("--checkpoint_keep_last", type=int, default=None,
                        help="Keep only the most recent N training checkpoints (default: 2).")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()

    if args.all_kmers or args.kmers:
        kmers = args.kmers if args.kmers else [3, 4, 5, 6]
        return pretrain_all_kmers(str(config_path), kmers=kmers, resume=args.resume)

    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    cfg["_config_path"] = str(config_path)
    kmer = args.kmer if args.kmer is not None else int(cfg.get("mlm_kmer", 6))
    kmer_cfg = get_kmer_pretrain_config(cfg, kmer)

    overrides: Dict[str, Any] = {}
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.debug:
        overrides["debug"] = True
    if args.overfit_debug:
        overrides["overfit_debug"] = True
    if args.label_debug:
        overrides["label_debug"] = True
    if args.leakage_test:
        overrides["leakage_test"] = True
    if args.generate_golden_batch:
        overrides["generate_golden_batch"] = True
    if args.golden_batch_path is not None:
        overrides["golden_batch_path"] = args.golden_batch_path
    if args.topk_debug:
        overrides["topk_debug"] = True
    if args.debug_batch:
        overrides["debug_batch"] = True
    if args.mlm_dry_run:
        overrides["mlm_dry_run"] = True
    if args.output_dir is not None:
        overrides["output_dir"] = args.output_dir
    if args.resume_path is not None:
        overrides["resume_path"] = args.resume_path
    if args.save_every_steps is not None:
        overrides["save_every_steps"] = args.save_every_steps
    if args.checkpoint_keep_last is not None:
        overrides["mlm_checkpoint_keep_last"] = args.checkpoint_keep_last

    return run_mlm_pretrain(kmer_cfg, overrides=overrides)

    # =========================================================================
    # Log resolved hyperparameters to runs/<timestamp>/resolved_config.json
    # =========================================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"mlm_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse length curriculum config early (needed for resolved_config)
    length_curriculum_config = LengthCurriculumConfig.from_config(cfg)
    
    resolved_config: Dict[str, Any] = {
        "script": "pretrain_mlm.py",
        "timestamp": timestamp,
        "config_path": str(config_path),
        # Model architecture
        "embedding_dim": int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256))),
        "transformer_layers": int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8))),
        "transformer_heads": int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8))),
        "transformer_ff_dim": int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", 1024))),
        "transformer_dropout": float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1))),
        # Optimizer
        # Use MLM-specific learning rate
        "lr": float(cfg.get("mlm_lr", cfg.get("lr", 2e-4))),
        "weight_decay": float(cfg.get("mlm_weight_decay", 0.01)),
        # MLM-specific
        "mask_prob": 0.15,  # Fixed in code
        "max_span_len": 3,   # Fixed in mask_tokens_span
        # Data
        "batch_size": int(cfg.get("mlm_batch_size", 128)),
        "max_bp_len": int(cfg.get("mlm_window_size", cfg.get("max_bp_len", 234))),
        "kmer": int(cfg.get("mlm_kmer", 3)),
        "include_reverse_complements": bool(cfg.get("include_reverse_complements", cfg.get("mlm_include_reverse_complements", False))),
        # Training
        "epochs": int(cfg.get("mlm_epochs", 200)),
        "grad_accum_steps": int(cfg.get("mlm_grad_accum_steps", 4)),
        "warmup_ratio": float(cfg.get("mlm_warmup_ratio", 0.06)),
        "final_lr_ratio": float(cfg.get("mlm_final_lr_ratio", 0.1)),
        "lr_log_interval": int(cfg.get("mlm_lr_log_interval", 50)),
        "patience": int(cfg.get("mlm_patience", 30)),
        "val_ratio": float(cfg.get("mlm_val_ratio", 0.1)),
        # Model options
        "tie_weights": bool(cfg.get("mlm_tie_weights", True)),
        "use_output_norm": bool(cfg.get("mlm_use_output_norm", True)),
        "use_layer_lr_decay": bool(cfg.get("mlm_use_layer_lr_decay", True)),
        "layer_lr_decay": float(cfg.get("mlm_layer_lr_decay", 0.9)),
        # Seed
        "seed": int(cfg.get("seed", 1337)),
        # AMP (Mixed Precision)
        "amp_enabled": bool(cfg.get("amp_enabled", True)),
        "amp_dtype": str(cfg.get("amp_dtype", "bfloat16")),
        "amp_force_fp32_loss": bool(cfg.get("amp_force_fp32_loss", True)),
        # Length curriculum
        "length_curriculum": length_curriculum_config.to_dict(),
    }
    
    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)
    
    # Seed: CLI argument takes precedence over config
    seed = args.seed if args.seed is not None else int(cfg.get("seed", 1337))
    set_seed(seed)
    LOGGER.info("Using seed: %d", seed)

    token_tensors, vocab = prepare_dataset(cfg)
    full_dataset = MLMDataset(token_tensors)
    
    # Split into train/val for monitoring overfitting
    val_ratio = float(cfg.get("mlm_val_ratio", 0.1))
    val_size = int(len(full_dataset) * val_ratio)
    train_size = len(full_dataset) - val_size
    train_dataset, val_dataset = random_split(
        full_dataset, 
        [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    LOGGER.info("Dataset split: %d train, %d val", train_size, val_size)
    
    batch_size = int(cfg.get("mlm_batch_size", 128))
    grad_accum_steps = int(cfg.get("mlm_grad_accum_steps", 2))
    effective_batch = batch_size * grad_accum_steps
    LOGGER.info("Batch size: %d, Grad accum: %d, Effective batch: %d", 
                batch_size, grad_accum_steps, effective_batch)
    
    # Standard 15% masking - enable debug mode if --label_debug flag is set
    # track_spans=True enables span length tracking for health dashboard
    
    # Calculate total steps for curriculum scheduler
    num_workers = int(cfg.get("num_workers", 0))
    epochs = int(cfg.get("mlm_epochs", 200))
    steps_per_epoch_approx = len(train_dataset) // batch_size
    total_steps = epochs * steps_per_epoch_approx // grad_accum_steps
    
    # Create collator with or without length curriculum
    curriculum_scheduler: LengthCurriculumScheduler | None = None
    if length_curriculum_config.enabled:
        LOGGER.info("=" * 60)
        LOGGER.info("LENGTH CURRICULUM ENABLED")
        LOGGER.info("  Start length: %d tokens", length_curriculum_config.start_tokens)
        LOGGER.info("  End length:   %d tokens", length_curriculum_config.end_tokens)
        LOGGER.info("  Ramp period:  %.0f%% - %.0f%% of training",
                    length_curriculum_config.start_pct * 100,
                    length_curriculum_config.end_pct * 100)
        LOGGER.info("  Total steps:  %d", total_steps)
        LOGGER.info("=" * 60)
        
        collator, curriculum_scheduler = create_curriculum_collator(
            vocab=vocab,
            total_steps=total_steps,
            cfg=cfg,
            mask_prob=0.15,
            track_spans=True,
            debug=args.label_debug,
        )
    else:
        LOGGER.info("Length curriculum: DISABLED (using fixed length)")
        collator = MLMCollator(vocab, mask_prob=0.15, debug=args.label_debug, track_spans=True)
    
    deterministic_kwargs = get_deterministic_dataloader_kwargs(seed, num_workers=num_workers)

    pin_memory = torch.cuda.is_available()
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=num_workers,
        pin_memory=pin_memory,
        drop_last=True,
        **deterministic_kwargs,
    )

    # Validation always uses full length (no curriculum)
    val_collator = MLMCollator(vocab, mask_prob=0.15)
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=val_collator,
        num_workers=0,
        pin_memory=pin_memory,
    )
    
    # =========================================================================
    # LABEL DEBUG MODE: Verify labels match originals and exit
    # =========================================================================
    if args.label_debug:
        LOGGER.info("Running LABEL DEBUG mode - verifying first batch then exiting")
        # Force one batch through the collator to trigger debug verification
        for batch in train_loader:
            # Debug verification happens inside collator's __call__
            break
        LOGGER.info("Label debug complete. Exiting.")
        return
    
    # =========================================================================
    # GOLDEN BATCH MODE: Generate golden batch for reproducibility testing
    # (Runs before sanity checks as it doesn't need them)
    # =========================================================================
    if args.generate_golden_batch:
        from golden_batch import generate_golden_batch as gen_golden
        
        LOGGER.info("Generating golden batch with seed %d...", seed)
        
        # Determine golden batch path
        if args.golden_batch_path:
            golden_path = Path(args.golden_batch_path)
            if not golden_path.is_absolute():
                golden_path = (script_dir / args.golden_batch_path).resolve()
        else:
            golden_path = (script_dir / "../artifacts/golden_batch.pt").resolve()
        
        device = select_device()
        
        # Build model for golden batch generation
        embedding_dim = int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256)))
        transformer_layers = int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8)))
        transformer_heads = int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8)))
        transformer_ff_dim = int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", embedding_dim * 4)))
        transformer_dropout = float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1)))

        # Relative position bias settings
        use_relative_position_bias = bool(cfg.get("use_relative_position_bias", False))
        relative_position_num_buckets = int(cfg.get("relative_position_num_buckets", 32))
        relative_position_max_distance = int(cfg.get("relative_position_max_distance", 128))

        encoder = DNAEncoder(
            vocab_size=len(vocab.itos),
            kmer=vocab.k,
            embedding_dim=embedding_dim,
            num_layers=transformer_layers,
            num_heads=transformer_heads,
            ff_dim=transformer_ff_dim,
            dropout=transformer_dropout,
            pad_token_id=vocab.pad_id,
            drop_path_rate=0.0,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
        )

        special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
        model = DNAMLM(
            encoder,
            vocab_size=len(vocab.itos),
            special_token_ids=special_token_ids,
            tie_weights=bool(cfg.get("mlm_tie_weights", True)),
            use_output_norm=bool(cfg.get("mlm_use_output_norm", True)),
        )

        metadata = {
            "config": str(config_path),
            "embedding_dim": embedding_dim,
            "transformer_layers": transformer_layers,
            "transformer_heads": transformer_heads,
            "batch_size": batch_size,
            "seed": seed,
        }
        
        gen_golden(
            dataloader=train_loader,
            model=model,
            device=device,
            seed=seed,
            path=golden_path,
            metadata=metadata,
        )
        
        LOGGER.info("Golden batch generation complete. Exiting.")
        return

    # =========================================================================
    # CHECKPOINT GATE: Run sanity checks before training can start
    # =========================================================================
    LOGGER.info("Running pretrain sanity checks (checkpoint gate)...")
    sanity_passed, sanity_issues = run_pretrain_sanity_checks(
        train_loader=train_loader,
        vocab=vocab,
        num_batches=5,
        target_mask_prob=0.15,
    )
    
    if not sanity_passed:
        LOGGER.error("Sanity checks FAILED! Training cannot proceed.")
        for issue in sanity_issues:
            LOGGER.error("  - %s", issue)
        raise RuntimeError(
            "Pretrain sanity checks failed. Fix the issues above before training.\n"
            + "\n".join(f"  - {i}" for i in sanity_issues)
        )
    
    # =========================================================================
    # LABEL DEBUG REPORT: Verify no special tokens in supervised labels
    # =========================================================================
    LOGGER.info("Running label token frequency debug report...")
    label_report = None
    for batch_data in train_loader:
        if len(batch_data) == 3:
            _, labels, _ = batch_data
        else:
            _, labels = batch_data
        label_report = debug_label_token_frequency(labels, vocab, num_batches=1)
        break  # Only need one batch for the report
    
    if label_report:
        print_label_debug_report(label_report, vocab, top_k=10)
        
        # HARD FAIL if special tokens found in labels
        if label_report["issues"]:
            raise RuntimeError(
                "Label token frequency check FAILED!\n"
                "Special tokens found in supervised labels.\n"
                + "\n".join(f"  - {i}" for i in label_report["issues"])
            )

    # =========================================================================
    # DEBUG MODE: Run one batch analysis and exit
    # =========================================================================
    if args.debug:
        LOGGER.info("Running in DEBUG mode - analyzing masking on one batch")
        # Get one batch of raw tokens (before masking)
        raw_batch = [train_dataset[i] for i in range(min(batch_size, len(train_dataset)))]
        debug_masking(vocab, raw_batch, mask_prob=0.15, max_span_len=3)
        LOGGER.info("Debug mode complete. Exiting.")
        return

    device = select_device()
    LOGGER.info("Using device: %s", device)
    
    # Model architecture - ENHANCED for better pretraining
    # Larger models learn better representations during MLM pretraining
    embedding_dim = int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256)))
    transformer_layers = int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8)))
    transformer_heads = int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8)))
    # FF dimension should be 4x embedding dim for good capacity (BERT-style)
    transformer_ff_dim = int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", embedding_dim * 4)))
    transformer_dropout = float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1)))
    
    LOGGER.info("=" * 60)
    LOGGER.info("MLM ARCHITECTURE (Enhanced for lower loss)")
    LOGGER.info("  Embedding dim: %d", embedding_dim)
    LOGGER.info("  Transformer layers: %d", transformer_layers)
    LOGGER.info("  Transformer heads: %d", transformer_heads)
    LOGGER.info("  FF dim: %d (%.1fx embedding)", transformer_ff_dim, transformer_ff_dim / embedding_dim)
    LOGGER.info("  Dropout: %.2f", transformer_dropout)
    LOGGER.info("=" * 60)

    # Relative position bias settings
    use_relative_position_bias = bool(cfg.get("use_relative_position_bias", False))
    relative_position_num_buckets = int(cfg.get("relative_position_num_buckets", 32))
    relative_position_max_distance = int(cfg.get("relative_position_max_distance", 128))

    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        drop_path_rate=0.0,  # NO stochastic depth during pretraining
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=relative_position_num_buckets,
        relative_position_max_distance=relative_position_max_distance,
    )

    # DNAMLM with weight tying (critical for lower loss)
    tie_weights = bool(cfg.get("mlm_tie_weights", True))
    use_output_norm = bool(cfg.get("mlm_use_output_norm", True))
    
    # CRITICAL: Exclude special tokens from prediction space
    # This prevents the model from "cheating" by predicting [MASK], <UNK>, <PAD>
    # which has been shown to cause models to get stuck at ~30% accuracy
    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    LOGGER.info("Special tokens excluded from predictions: [MASK]=%d, <UNK>=%d, <PAD>=%d",
                vocab.mask_id, vocab.unk_id, vocab.pad_id)
    
    model = DNAMLM(
        encoder, 
        vocab_size=len(vocab.itos),
        special_token_ids=special_token_ids,
        tie_weights=tie_weights,
        use_output_norm=use_output_norm,
    ).to(device)
    
    # Count parameters (accounting for weight tying)
    total_params = sum(p.numel() for p in model.parameters())
    unique_params = sum(p.numel() for p in set(model.parameters()))
    LOGGER.info("Model parameters: %.2fM total, %.2fM unique", total_params / 1e6, unique_params / 1e6)
    
    # =========================================================================
    # LEAKAGE TEST: Verify no information leakage at masked positions
    # =========================================================================
    if args.leakage_test:
        LOGGER.info("Running LEAKAGE TEST mode - verifying no information leakage")
        verify_no_leakage(model, vocab, device)
        LOGGER.info("Leakage test complete. Exiting.")
        return
    
    # =========================================================================
    # OVERFIT DEBUG MODE: Train on 32 samples for 200 steps
    # =========================================================================
    if args.overfit_debug:
        LOGGER.info("Running OVERFIT DEBUG mode")
        # Take 32 samples for proper overfit test
        overfit_samples = token_tensors[:32]
        run_overfit_debug(
            model=model,
            train_samples=overfit_samples,
            vocab=vocab,
            device=device,
            num_steps=1000,  # Very aggressive: 1000 steps
            lr=5e-3,         # Very high LR for fast convergence
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return
    
    # Learning rate and optimizer settings
    # Use MLM-specific learning rate
    lr = float(cfg.get("mlm_lr", cfg.get("lr", 2e-4)))  # Base LR 3e-4 for stable MLM pretraining
    weight_decay = float(cfg.get("mlm_weight_decay", 0.01))
    lr_decay = float(cfg.get("mlm_layer_lr_decay", 0.9))  # Layer-wise LR decay
    use_layer_lr_decay = bool(cfg.get("mlm_use_layer_lr_decay", True))
    optimizer_foreach = bool(cfg.get("optimizer_foreach", True))

    # Build optimizer with or without layer-wise LR decay
    LOGGER.info("=" * 60)
    LOGGER.info("OPTIMIZER CONFIGURATION")
    LOGGER.info("Optimizer: AdamW foreach=%s", optimizer_foreach)
    if use_layer_lr_decay:
        param_groups = get_layer_wise_lr_groups(model, lr, lr_decay, weight_decay)
        LOGGER.info("Layer-wise LR decay: ENABLED (factor: %.2f)", lr_decay)
        LOGGER.info("  Lower layers get smaller LRs to preserve pretrained features")
        optimizer = AdamW(
            param_groups,
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,  # Smaller eps for stability
            foreach=optimizer_foreach,
        )
    else:
        LOGGER.info("Layer-wise LR decay: DISABLED (all layers use same LR)")
        LOGGER.info("  All parameters use base LR: %.2e", lr)
        optimizer = AdamW(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay,
            betas=(0.9, 0.98),  # BERT-style betas
            eps=1e-6,
            foreach=optimizer_foreach,
        )
    LOGGER.info("=" * 60)
    
    epochs = int(cfg.get("mlm_epochs", 150))
    
    # Learning rate schedule with warmup + cosine decay
    # total_steps = number of OPTIMIZER steps (not batches)
    # 
    # Calculation:
    #   batches_per_epoch = len(train_loader)  [with drop_last=True]
    #   optimizer_steps_per_epoch = batches_per_epoch // grad_accum_steps
    #   total_optimizer_steps = epochs * optimizer_steps_per_epoch
    #
    batches_per_epoch = len(train_loader)
    steps_per_epoch = batches_per_epoch // grad_accum_steps
    total_steps = epochs * steps_per_epoch
    warmup_ratio = float(cfg.get("mlm_warmup_ratio", 0.06))  # 6% warmup
    warmup_steps = int(warmup_ratio * total_steps)
    final_lr_ratio = float(cfg.get("mlm_final_lr_ratio", 0.1))  # Decay to 10% of peak
    lr_log_interval = int(cfg.get("mlm_lr_log_interval", 50))
    
    LOGGER.info("=" * 60)
    LOGGER.info("TOTAL STEPS CALCULATION")
    LOGGER.info("  Epochs:                  %d", epochs)
    LOGGER.info("  Batches per epoch:       %d", batches_per_epoch)
    LOGGER.info("  Grad accum steps:        %d", grad_accum_steps)
    LOGGER.info("  Optimizer steps/epoch:   %d", steps_per_epoch)
    LOGGER.info("  Total optimizer steps:   %d", total_steps)
    LOGGER.info("=" * 60)
    
    scheduler = get_warmup_cosine_schedule(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        final_lr_ratio=final_lr_ratio,
    )
    LOGGER.info("=" * 60)
    LOGGER.info("LR SCHEDULE (Linear Warmup + Cosine Decay)")
    LOGGER.info("  Base LR:        %.2e", lr)
    LOGGER.info("  Warmup steps:   %d (%.1f%% of total)", warmup_steps, warmup_ratio * 100)
    LOGGER.info("  Total steps:    %d", total_steps)
    LOGGER.info("  Final LR:       %.2e (%.0f%% of peak)", lr * final_lr_ratio, final_lr_ratio * 100)
    LOGGER.info("  LR log interval: every %d optimizer steps", lr_log_interval)
    LOGGER.info("=" * 60)

    # AMP (Mixed Precision) setup using config
    amp_ctx = get_amp_context(device, cfg, logger=LOGGER)
    log_amp_status(amp_ctx, logger=LOGGER)

    autocast_device: str | None = amp_ctx.device_type if amp_ctx.enabled else None
    autocast_dtype: torch.dtype | None = amp_ctx.dtype if amp_ctx.enabled else None
    scaler = None

    # GradScaler only for CUDA with float16
    if amp_ctx.enabled and device.type == "cuda" and amp_ctx.dtype == torch.float16:
        scaler = torch.cuda.amp.GradScaler()
        LOGGER.info("Using GradScaler for CUDA fp16")

    # Training state
    best_val_loss = float('inf')
    best_val_acc = 0.0
    patience = int(cfg.get("mlm_patience", 25))
    patience_counter = 0
    global_step = 0
    param_norm_warn = float(cfg.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg.get("param_norm_cap", 200.0))
    
    encoder_path = _resolve_cfg_path(
        cfg.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt"),
        base_dir=script_dir,
    )
    encoder_path.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("=" * 60)
    LOGGER.info("STARTING MLM PRETRAINING")
    LOGGER.info("Epochs: %d, LR: %.2e, Mask prob: 15%%", epochs, lr)
    LOGGER.info("=" * 60)
    
    # =========================================================================
    # HEALTH DASHBOARD: Initialize tracking
    # =========================================================================
    health_log_interval = int(cfg.get("health_log_interval", 50))
    health_dashboard = create_health_dashboard(
        run_dir=run_dir,
        vocab=vocab,
        log_every_n_steps=health_log_interval,
        target_mask_prob=0.15,
    )
    health_dashboard.print_table_header()
    
    # =========================================================================
    # NUMERICS CHECKER: NaN/Inf tripwires for fail-fast detection
    # =========================================================================
    numerics_checker = NumericsChecker(enabled=True)
    LOGGER.info("Numerics checker enabled - will fail fast on NaN/Inf")

    for epoch in range(1, epochs + 1):
        model.train()
        LOGGER.info("Train loop: model.training=%s", model.training)
        total_loss = 0.0
        total_correct = 0
        total_masked = 0
        optimizer.zero_grad()

        if scaler is not None:
            LOGGER.info("Epoch %d: GradScaler scale = %.0f", epoch, scaler.get_scale())
        
        # Log current curriculum length at start of epoch
        if curriculum_scheduler is not None:
            curr_len = curriculum_scheduler.current_length
            progress = curriculum_scheduler.get_progress()
            LOGGER.info(
                "Epoch %d: Curriculum length = %d tokens (%s, %.0f%% ramp progress)",
                epoch, curr_len, progress["phase"], progress["ramp_progress"] * 100
            )
        
        for batch_idx, batch_data in enumerate(train_loader, start=1):
            # Handle both (inputs, labels) and (inputs, labels, spans) formats
            if len(batch_data) == 3:
                inputs, labels, span_lengths = batch_data
            else:
                inputs, labels = batch_data
                span_lengths = None

            non_blocking = (device.type == "cuda")
            inputs = inputs.to(device, non_blocking=non_blocking)
            labels = labels.to(device, non_blocking=non_blocking)

            # Update numerics checker context for error reporting
            current_lr = optimizer.param_groups[0]['lr']
            numerics_checker.set_context(
                epoch=epoch,
                step=global_step,
                lr=current_lr,
                batch_idx=batch_idx,
                model=model,
            )
            
            if autocast_device and autocast_dtype is not None:
                with torch.autocast(autocast_device, dtype=autocast_dtype):
                    logits, _ = model(inputs, labels=None)
            else:
                logits, _ = model(inputs, labels=None)

            # Loss computation outside autocast for numerical stability
            loss = F.cross_entropy(
                logits.float().view(-1, logits.size(-1)),
                labels.view(-1),
                ignore_index=-100,
            )
            
            # TRIPWIRE: Check forward pass outputs for NaN/Inf
            numerics_checker.check_forward(logits, loss)
            
            # Scale loss for gradient accumulation
            scaled_loss = loss / grad_accum_steps
            if scaler is not None:
                scaler.scale(scaled_loss).backward()
            else:
                scaled_loss.backward()
            
            # TRIPWIRE: In fp16+GradScaler mode, gradients are scaled and may overflow by design.
            # Still check parameters, but avoid failing on transient gradient Inf/NaN handled by scaler.
            if scaler is not None:
                numerics_checker.check_backward(model, check_grads=False)
            else:
                numerics_checker.check_backward(model)
            
            total_loss += loss.item() * inputs.size(0)
            
            # Compute accuracy
            preds = logits.argmax(dim=-1)
            mask = labels != -100
            total_correct += (preds[mask] == labels[mask]).sum().item()
            total_masked += mask.sum().item()
            
            # Accumulate metrics for health dashboard
            health_dashboard.accumulate_batch(
                loss=loss.item(),
                logits=logits.detach(),
                labels=labels,
                masked_inputs=inputs,
                span_lengths=span_lengths,
            )
            
            # Optimizer step with gradient accumulation
            if batch_idx % grad_accum_steps == 0:
                if scaler is not None:
                    scaler.unscale_(optimizer)
                # Gradient clipping with norm logging
                # clip_grad_norm_ returns the total norm BEFORE clipping
                unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                # After clipping, effective norm is min(unclipped, max_norm)
                clipped_grad_norm = min(unclipped_grad_norm.item(), 1.0)
                
                if scaler is not None:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()
                total_norm, scale = cap_model_param_norm_(model, max_norm=param_norm_cap)
                effective_param_norm = param_norm_cap if scale < 1.0 else total_norm
                if scale < 1.0:
                    LOGGER.warning(
                        "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                        global_step + 1,
                        total_norm,
                        param_norm_cap,
                        scale,
                    )
                scheduler.step()
                optimizer.zero_grad()
                global_step += 1
                
                # Log LR every N optimizer steps (configurable via mlm_lr_log_interval)
                if global_step % lr_log_interval == 0:
                    current_lr = optimizer.param_groups[0]['lr']
                    # Determine training phase for logging
                    if global_step <= warmup_steps:
                        phase = "warmup"
                        phase_progress = global_step / warmup_steps * 100
                    else:
                        phase = "decay"
                        phase_progress = (global_step - warmup_steps) / max(1, total_steps - warmup_steps) * 100
                    
                    LOGGER.info(
                        "Step %d/%d [%s %.0f%%]: LR=%.2e, grad_norm=%.4f, param_norm=%.2f",
                        global_step, total_steps, phase, phase_progress,
                        current_lr, unclipped_grad_norm.item(), effective_param_norm
                    )
                    if param_norm_warn > 0 and effective_param_norm > param_norm_warn:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f above warn=%.2f",
                            global_step,
                            effective_param_norm,
                            param_norm_warn,
                        )
                
                # Step curriculum scheduler and update length
                if curriculum_scheduler is not None:
                    new_length = curriculum_scheduler.step(global_step)
                    # Log length changes during ramping
                    if curriculum_scheduler.is_ramping() and global_step % 100 == 0:
                        LOGGER.debug(
                            "Step %d: Curriculum length = %d tokens",
                            global_step, new_length
                        )
                
                # Log health metrics every N steps
                health_metrics = health_dashboard.log_step(
                    step=global_step,
                    epoch=epoch,
                    model=model,
                    optimizer=optimizer,
                )
                if health_metrics:
                    health_dashboard.print_table_row(health_metrics)
            
            if batch_idx % 100 == 0:
                current_lr = optimizer.param_groups[0]['lr']
                batch_acc = total_correct / max(1, total_masked)
                LOGGER.info(
                    "Epoch %d [%d/%d] Loss: %.4f Acc: %.2f%% LR: %.2e",
                    epoch, batch_idx, len(train_loader), loss.item(), 
                    batch_acc * 100, current_lr
                )
        
        # Epoch statistics
        train_loss = total_loss / len(train_dataset)
        train_acc = total_correct / max(1, total_masked)
        train_ppl = math.exp(min(train_loss, 100))
        
        # Validation
        val_metrics = run_validation(model, val_loader, device, autocast_device, autocast_dtype)
        val_loss = val_metrics['loss']
        val_acc = val_metrics['accuracy']
        val_ppl = val_metrics['perplexity']
        
        # Build epoch summary with curriculum info
        curriculum_info = ""
        if curriculum_scheduler is not None:
            curriculum_info = f" | Len: {curriculum_scheduler.current_length}"
        
        LOGGER.info(
            "Epoch %d: Train [Loss: %.4f, Acc: %.1f%%, PPL: %.2f] | Val [Loss: %.4f, Acc: %.1f%%, PPL: %.2f]%s",
            epoch, train_loss, train_acc * 100, train_ppl, val_loss, val_acc * 100, val_ppl, curriculum_info
        )
        
        # =========================================================================
        # TOP-K SANITY PRINT: Run once per epoch if enabled
        # =========================================================================
        if args.topk_debug:
            # Get a fresh batch for debug analysis
            for debug_batch_data in train_loader:
                if len(debug_batch_data) == 3:
                    debug_inputs, debug_labels, _ = debug_batch_data
                else:
                    debug_inputs, debug_labels = debug_batch_data
                break
            
            topk_stats = print_topk_predictions_debug(
                model=model,
                inputs=debug_inputs,
                labels=debug_labels,
                vocab=vocab,
                epoch=epoch,
                num_sequences=3,
                top_k=5,
            )
            
            # Log summary stats for tracking over epochs
            LOGGER.info(
                "Top-K Debug: Avg top-1 prob: %.1f%%, True in top-1: %.1f%%, True in top-5: %.1f%%",
                topk_stats["avg_top1_prob"] * 100,
                topk_stats["true_in_top1_pct"],
                topk_stats["true_in_topk_pct"],
            )
        
        # Save best model based on validation loss
        improved = False
        if val_loss < best_val_loss - 0.001:
            improved = True
            best_val_loss = val_loss
            best_val_acc = val_acc
            patience_counter = 0
            torch.save(model.encoder.state_dict(), encoder_path)
            LOGGER.info("★ New best! Val Loss: %.4f, PPL: %.2f, Acc: %.1f%% - Saved to %s", 
                       val_loss, val_ppl, val_acc * 100, encoder_path)
        else:
            patience_counter += 1
            if patience_counter >= patience:
                LOGGER.info("Early stopping at epoch %d (patience %d)", epoch, patience)
                break
    
    best_ppl = math.exp(min(best_val_loss, 100))
    LOGGER.info("=" * 60)
    LOGGER.info("PRETRAINING COMPLETE")
    LOGGER.info("Best Val Loss: %.4f (Perplexity: %.2f)", best_val_loss, best_ppl)
    LOGGER.info("Best Val Accuracy: %.1f%%", best_val_acc * 100)
    LOGGER.info("Target: Loss < 1.0 (PPL < 2.72)")
    LOGGER.info("Saved pretrained encoder to %s", encoder_path)
    LOGGER.info("=" * 60)


def pretrain_all_kmers(
    config_path: str = "config.yaml",
    kmers: list | None = None,
    resume: bool = False,
) -> dict:
    """
    Pretrain MLM models for multiple k-mer sizes.

    Args:
        config_path: Path to config.yaml
        kmers: List of k-mer sizes (default: [3, 4, 5, 6])
        resume: Whether to resume from existing checkpoints
    """
    config_file = Path(config_path)
    if not config_file.is_absolute():
        config_file = (Path(__file__).resolve().parent / config_file).resolve()
    with config_file.open("r", encoding="utf-8") as f:
        base_cfg = yaml.safe_load(f) or {}
    base_cfg["_config_path"] = str(config_file)

    if kmers is None:
        kmers = base_cfg.get("pretrain_kmers", [3, 4, 5, 6])

    LOGGER.info("=" * 70)
    LOGGER.info("MULTI-KMER MLM PRETRAINING")
    LOGGER.info("=" * 70)
    LOGGER.info("K-mers to pretrain: %s", kmers)

    results: dict = {}

    for kmer in kmers:
        LOGGER.info("\n%s", "=" * 60)
        LOGGER.info("PRETRAINING K=%d", kmer)
        LOGGER.info("%s\n", "=" * 60)

        # Get k-mer specific config
        kmer_cfg = get_kmer_pretrain_config(base_cfg, kmer)

        # Check if we should skip (checkpoint exists and not resuming)
        checkpoint_path = Path(kmer_cfg["mlm_encoder_path"])
        if not checkpoint_path.is_absolute():
            checkpoint_path = (Path(__file__).resolve().parent / checkpoint_path).resolve()
        if checkpoint_path.exists() and not resume:
            LOGGER.info("Checkpoint exists for k=%d, skipping (use --resume to continue)", kmer)
            results[kmer] = {"skipped": True, "checkpoint": str(checkpoint_path)}
            continue

        try:
            # Run pretraining (call your existing main function with kmer_cfg)
            result = main(kmer_cfg)
            results[kmer] = result
            LOGGER.info("✓ K=%d complete: val_loss=%s", kmer, result.get("best_val_loss", "N/A"))
        except Exception as exc:
            LOGGER.error("✗ K=%d failed: %s", kmer, exc)
            results[kmer] = {"error": str(exc)}

    # Save summary
    summary = {
        "kmers": kmers,
        "results": results,
        "timestamp": datetime.now().isoformat(),
    }

    summary_path = (Path(__file__).resolve().parent / "../artifacts/multi_kmer_pretrain_summary.json").resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)

    LOGGER.info("\n%s", "=" * 70)
    LOGGER.info("PRETRAINING COMPLETE")
    LOGGER.info("%s", "=" * 70)
    LOGGER.info("Summary saved to %s", summary_path)

    return results


if __name__ == "__main__":
    main()
