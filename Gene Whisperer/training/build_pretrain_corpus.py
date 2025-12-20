"""Build pretraining corpus from genomic FASTA files.

This script generates unlabeled DNA sequences for masked language model pretraining.
It supports:
- Sliding (overlapping) windows across a genome (recommended for large corpora)
- Random window sampling from genomic sequences (legacy / debugging)
- Promoter-enriched sampling (biases towards upstream regions)
- Shift augmentation (random +/- offsets around reference positions)
- N-content filtering and sanitization

Output: Line-delimited sequences in .txt or .tsv format with a summary report.

Usage:
    python build_pretrain_corpus.py --config config.yaml
    python build_pretrain_corpus.py --fasta ../data/genome.fna --seq_len 81 --num_samples 100000
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import random
import re
import statistics
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import yaml

from seed_utils import set_global_seed

LOGGER = logging.getLogger("gene_whisperer.corpus_builder")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Valid DNA bases
ALLOWED_BASES: Set[str] = {"A", "C", "G", "T"}
REV_COMP_MAP = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of an A/C/G/T-only DNA sequence."""
    return sequence.translate(REV_COMP_MAP)[::-1]


@dataclass
class CorpusStats:
    """Statistics about the generated corpus."""
    
    total_sequences: int = 0
    length_min: int = 0
    length_max: int = 0
    length_mean: float = 0.0
    length_std: float = 0.0
    
    # N-content stats (before sanitization)
    sequences_with_n: int = 0
    total_n_bases: int = 0
    n_filtered_sequences: int = 0

    # Windowing stats (primarily for sliding-window mode)
    total_windows_considered: int = 0
    skipped_windows: int = 0
    invalid_filtered_windows: int = 0
    
    # Sampling stats
    sliding_samples: int = 0
    random_samples: int = 0
    promoter_enriched_samples: int = 0
    shift_augmented_samples: int = 0

    # Reverse-complement augmentation stats
    reverse_complement_samples: int = 0
    reverse_complement_fraction: float = 0.0
    
    # Base composition
    base_counts: Dict[str, int] = field(default_factory=dict)
    
    def compute_base_percentages(self) -> Dict[str, float]:
        """Compute percentage of each base in corpus."""
        total = sum(self.base_counts.values())
        if total == 0:
            return {}
        return {base: count / total * 100 for base, count in self.base_counts.items()}


def read_fasta_sequences(path: Path) -> Dict[str, str]:
    """
    Read all sequences from a FASTA file.
    
    Returns a dictionary mapping sequence ID to sequence string.
    Handles multi-line sequences and multiple records.
    """
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_parts: List[str] = []
    
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = "".join(current_parts).upper()
                
                # Start new sequence
                # Extract ID (first word after >)
                header = line[1:].strip()
                current_id = header.split()[0] if header else f"seq_{len(sequences)}"
                current_parts = []
            else:
                current_parts.append(line.upper())
    
    # Save last sequence
    if current_id is not None:
        sequences[current_id] = "".join(current_parts).upper()
    
    return sequences


def count_n_content(sequence: str) -> Tuple[int, int]:
    """
    Count N bases and total bases in a sequence.
    
    Returns: (n_count, total_length)
    """
    n_count = sequence.upper().count("N")
    return n_count, len(sequence)


def sanitize_sequence(
    sequence: str,
    n_handling: str = "filter",
    n_replacement: str = "random",
) -> Optional[str]:
    """
    Sanitize a DNA sequence for pretraining.
    
    Args:
        sequence: Raw DNA sequence
        n_handling: How to handle N bases:
            - "filter": Return None if sequence contains N
            - "replace": Replace N with specified strategy
            - "skip": Remove N bases from sequence
        n_replacement: When n_handling="replace", how to replace:
            - "random": Random ACGT base
            - "A", "C", "G", "T": Specific base
            
    Returns:
        Sanitized sequence (A/C/G/T only) or None if filtered
    """
    seq = sequence.upper().strip()
    
    # Handle N bases
    has_n = "N" in seq
    
    if has_n:
        if n_handling == "filter":
            return None
        elif n_handling == "replace":
            if n_replacement == "random":
                # Replace each N with random base
                result = []
                for base in seq:
                    if base == "N":
                        result.append(random.choice(["A", "C", "G", "T"]))
                    else:
                        result.append(base)
                seq = "".join(result)
            else:
                # Replace with specific base
                seq = seq.replace("N", n_replacement.upper())
        elif n_handling == "skip":
            seq = seq.replace("N", "")
    
    # Keep only valid bases
    sanitized = "".join(base for base in seq if base in ALLOWED_BASES)
    
    return sanitized if sanitized else None


def sample_random_windows(
    sequence: str,
    window_size: int,
    num_samples: int,
    n_handling: str = "filter",
    n_replacement: str = "random",
    max_n_ratio: float = 0.1,
) -> Tuple[List[str], int, int, int]:
    """
    Sample random windows from a genomic sequence.
    
    Args:
        sequence: Full genomic sequence
        window_size: Length of each sampled window
        num_samples: Target number of samples
        n_handling: How to handle N bases
        n_replacement: Replacement strategy for N
        max_n_ratio: Maximum ratio of N bases to attempt replacement
        
    Returns:
        (windows, n_filtered_count, n_containing_count, total_n_bases)
    """
    windows: List[str] = []
    seq_len = len(sequence)
    n_filtered = 0
    n_containing = 0
    total_n = 0
    
    if seq_len < window_size:
        LOGGER.warning("Sequence length %d < window size %d", seq_len, window_size)
        return [], 0, 0, 0
    
    max_attempts = num_samples * 3  # Allow some failures
    attempts = 0
    
    while len(windows) < num_samples and attempts < max_attempts:
        attempts += 1
        
        # Random start position
        start = random.randint(0, seq_len - window_size)
        raw_window = sequence[start:start + window_size]
        
        # Count N content
        n_count, _ = count_n_content(raw_window)
        if n_count > 0:
            n_containing += 1
            total_n += n_count
            
            # Check if too many Ns
            if n_count / window_size > max_n_ratio and n_handling != "filter":
                n_filtered += 1
                continue
        
        # Sanitize
        clean_window = sanitize_sequence(raw_window, n_handling, n_replacement)
        
        if clean_window is None:
            n_filtered += 1
            continue
        
        # Ensure correct length after sanitization (skip mode can change length)
        if len(clean_window) != window_size:
            continue
            
        windows.append(clean_window)
    
    return windows, n_filtered, n_containing, total_n


def extract_overlapping_windows(
    sequence: str,
    window_size: int,
    stride: int,
    n_handling: str = "filter",
    n_replacement: str = "random",
    max_n_ratio: float = 0.1,
) -> Tuple[List[str], int, int, int, int, int]:
    """
    Extract overlapping (sliding) windows from a genomic sequence.

    Returns:
        (windows, n_filtered_count, n_containing_count, total_n_bases,
         total_windows_considered, invalid_filtered_windows)
    """
    if stride <= 0:
        raise ValueError("stride must be > 0")

    windows: List[str] = []
    seq_len = len(sequence)
    n_filtered = 0
    n_containing = 0
    total_n = 0
    total_windows = 0
    invalid_filtered = 0

    if seq_len < window_size:
        LOGGER.warning("Sequence length %d < window size %d", seq_len, window_size)
        return [], 0, 0, 0, 0, 0

    max_start = seq_len - window_size + 1
    for start in range(0, max_start, stride):
        total_windows += 1
        raw_window = sequence[start : start + window_size]

        n_count, _ = count_n_content(raw_window)
        if n_count > 0:
            n_containing += 1
            total_n += n_count

            if n_handling != "filter" and (n_count / window_size) > max_n_ratio:
                n_filtered += 1
                continue

        clean_window = sanitize_sequence(raw_window, n_handling, n_replacement)
        if clean_window is None:
            n_filtered += 1
            continue

        if len(clean_window) != window_size:
            invalid_filtered += 1
            continue

        if not set(clean_window) <= ALLOWED_BASES:
            invalid_filtered += 1
            continue

        windows.append(clean_window)

    return windows, n_filtered, n_containing, total_n, total_windows, invalid_filtered


def load_promoter_positions(
    tsv_path: Path,
    sequence_col: str = "sequence",
    label_col: str = "is_promoter",
) -> List[str]:
    """
    Load promoter sequences from a TSV file.
    
    Returns list of promoter sequences (where label == 1).
    """
    promoters: List[str] = []
    
    if not tsv_path.exists():
        LOGGER.warning("Promoter file not found: %s", tsv_path)
        return promoters
    
    with tsv_path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        
        try:
            seq_idx = header.index(sequence_col)
            label_idx = header.index(label_col)
        except ValueError as e:
            LOGGER.warning("Could not find columns in %s: %s", tsv_path, e)
            return promoters
        
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > max(seq_idx, label_idx):
                try:
                    if int(parts[label_idx]) == 1:
                        promoters.append(parts[seq_idx])
                except (ValueError, IndexError):
                    continue
    
    LOGGER.info("Loaded %d promoter sequences from %s", len(promoters), tsv_path.name)
    return promoters


def sample_promoter_enriched(
    sequence: str,
    window_size: int,
    num_samples: int,
    upstream_bias: int = 200,
    n_handling: str = "filter",
    n_replacement: str = "random",
) -> List[str]:
    """
    Sample windows with bias towards promoter-like upstream regions.
    
    This simulates promoter-enriched sampling by preferring positions
    upstream of potential gene starts (roughly every ~1kb in bacteria).
    
    Args:
        sequence: Full genomic sequence
        window_size: Length of each sampled window
        num_samples: Target number of samples
        upstream_bias: How many bp upstream of gene starts to prefer
        n_handling: How to handle N bases
        n_replacement: Replacement strategy
        
    Returns:
        List of sampled windows
    """
    windows: List[str] = []
    seq_len = len(sequence)
    
    if seq_len < window_size:
        return []
    
    # Estimate gene starts (roughly every 1kb for bacteria)
    avg_gene_spacing = 1000
    potential_gene_starts = list(range(
        upstream_bias + window_size,
        seq_len,
        avg_gene_spacing
    ))
    
    if not potential_gene_starts:
        # Fall back to random sampling
        return sample_random_windows(
            sequence, window_size, num_samples, n_handling, n_replacement
        )[0]
    
    max_attempts = num_samples * 3
    attempts = 0
    
    while len(windows) < num_samples and attempts < max_attempts:
        attempts += 1
        
        # Pick a gene start and sample upstream
        gene_start = random.choice(potential_gene_starts)
        
        # Random position in upstream region
        upstream_start = gene_start - upstream_bias - window_size
        if upstream_start < 0:
            upstream_start = 0
        
        jitter = random.randint(0, upstream_bias)
        start = max(0, upstream_start + jitter)
        
        if start + window_size > seq_len:
            continue
            
        raw_window = sequence[start:start + window_size]
        clean_window = sanitize_sequence(raw_window, n_handling, n_replacement)
        
        if clean_window and len(clean_window) == window_size:
            windows.append(clean_window)
    
    return windows


def apply_shift_augmentation(
    sequences: List[str],
    reference_sequence: str,
    window_size: int,
    shift_range: int = 10,
    augmentation_factor: float = 0.5,
    n_handling: str = "filter",
    n_replacement: str = "random",
) -> List[str]:
    """
    Apply shift augmentation to create additional samples.
    
    For each input sequence, attempts to find it in the reference and
    create shifted versions with random +/- offsets.
    
    Args:
        sequences: Original sequences to augment
        reference_sequence: Full reference genome for finding positions
        window_size: Target window size
        shift_range: Maximum shift in either direction (+/-)
        augmentation_factor: Fraction of original sequences to augment
        n_handling: How to handle N bases
        n_replacement: Replacement strategy
        
    Returns:
        List of shift-augmented sequences
    """
    augmented: List[str] = []
    
    num_to_augment = int(len(sequences) * augmentation_factor)
    to_augment = random.sample(sequences, min(num_to_augment, len(sequences)))
    
    ref_len = len(reference_sequence)
    
    for seq in to_augment:
        # Try to find this sequence in reference
        # (for efficiency, we'll just create shifted versions around a random position
        # since finding exact matches is expensive for large genomes)
        
        # Random shift
        shift = random.randint(-shift_range, shift_range)
        if shift == 0:
            continue
        
        # Create a shifted version by sampling adjacent region
        # We'll simulate this by taking a random position
        pos = random.randint(abs(shift), ref_len - window_size - abs(shift))
        shifted_start = pos + shift
        
        if shifted_start < 0 or shifted_start + window_size > ref_len:
            continue
            
        raw_window = reference_sequence[shifted_start:shifted_start + window_size]
        clean_window = sanitize_sequence(raw_window, n_handling, n_replacement)
        
        if clean_window and len(clean_window) == window_size:
            augmented.append(clean_window)
    
    return augmented


def compute_corpus_stats(
    sequences: List[str],
    n_filtered: int,
    n_containing: int,
    total_n: int,
    random_count: int,
    promoter_enriched_count: int,
    shift_augmented_count: int,
) -> CorpusStats:
    """Compute comprehensive statistics about the corpus."""
    stats = CorpusStats()
    
    if not sequences:
        return stats
    
    lengths = [len(s) for s in sequences]
    
    stats.total_sequences = len(sequences)
    stats.length_min = min(lengths)
    stats.length_max = max(lengths)
    stats.length_mean = statistics.mean(lengths)
    stats.length_std = statistics.stdev(lengths) if len(lengths) > 1 else 0.0
    
    stats.sequences_with_n = n_containing
    stats.total_n_bases = total_n
    stats.n_filtered_sequences = n_filtered
    
    stats.random_samples = random_count
    stats.promoter_enriched_samples = promoter_enriched_count
    stats.shift_augmented_samples = shift_augmented_count
    
    # Base composition
    stats.base_counts = Counter()
    for seq in sequences:
        stats.base_counts.update(seq)
    
    return stats


def write_corpus(
    sequences: List[str],
    output_path: Path,
    format: str = "txt",
    shuffle: bool = True,
) -> None:
    """
    Write sequences to output file.
    
    Args:
        sequences: List of DNA sequences
        output_path: Path to output file
        format: Output format ("txt" or "tsv")
        shuffle: Whether to shuffle sequences before writing
    """
    if shuffle:
        random.shuffle(sequences)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with output_path.open("w", encoding="utf-8") as f:
        if format == "tsv":
            f.write("sequence\n")
        
        for seq in sequences:
            f.write(seq + "\n")
    
    LOGGER.info("Wrote %d sequences to %s", len(sequences), output_path)


def write_report(
    stats: CorpusStats,
    output_path: Path,
    config: Dict[str, Any],
) -> None:
    """Write a detailed report about the corpus generation."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    report = {
        "generated_at": datetime.now().isoformat(),
        "config": config,
        "statistics": {
            "total_sequences": stats.total_sequences,
            "length_stats": {
                "min": stats.length_min,
                "max": stats.length_max,
                "mean": round(stats.length_mean, 2),
                "std": round(stats.length_std, 2),
            },
            "n_content": {
                "sequences_with_n_before_sanitization": stats.sequences_with_n,
                "total_n_bases_encountered": stats.total_n_bases,
                "sequences_filtered_due_to_n": stats.n_filtered_sequences,
            },
            "windowing": {
                "total_windows_considered": stats.total_windows_considered,
                "skipped_windows": stats.skipped_windows,
                "invalid_filtered_windows": stats.invalid_filtered_windows,
            },
            "sampling_breakdown": {
                "sliding_samples": stats.sliding_samples,
                "random_samples": stats.random_samples,
                "promoter_enriched_samples": stats.promoter_enriched_samples,
                "shift_augmented_samples": stats.shift_augmented_samples,
            },
            "reverse_complements": {
                "reverse_complement_samples": stats.reverse_complement_samples,
                "reverse_complement_fraction": round(stats.reverse_complement_fraction, 6),
            },
            "base_composition": stats.base_counts,
            "base_percentages": {
                k: round(v, 2) for k, v in stats.compute_base_percentages().items()
            },
        },
    }
    
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    
    LOGGER.info("Wrote report to %s", output_path)


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_pretrain_corpus(
    fasta_path: Path,
    output_path: Path,
    seq_len: int,
    num_samples: Optional[int],
    seed: int = 1337,
    n_handling: str = "filter",
    n_replacement: str = "random",
    max_n_ratio: float = 0.1,
    sampling: str = "sliding",
    stride: int = 20,
    promoter_enriched: bool = False,
    promoter_enriched_ratio: float = 0.3,
    promoter_tsv_path: Optional[Path] = None,
    shift_augmentation: bool = False,
    shift_range: int = 10,
    shift_augmentation_factor: float = 0.2,
    include_reverse_complements: bool = False,
    output_format: str = "txt",
    shuffle: bool = True,
) -> CorpusStats:
    """
    Build a pretraining corpus from genomic FASTA file.
    
    Args:
        fasta_path: Path to input FASTA file
        output_path: Path to output corpus file
        seq_len: Sequence length (window size)
        num_samples: Total target number of samples (random mode only)
        seed: Random seed for reproducibility
        n_handling: How to handle N bases ("filter", "replace", "skip")
        n_replacement: Replacement strategy when n_handling="replace"
        max_n_ratio: Maximum N ratio before filtering (when not using "filter")
        sampling: Sampling strategy ("sliding" or "random")
        stride: Step size between consecutive windows (sliding mode only)
        promoter_enriched: Enable promoter-enriched sampling
        promoter_enriched_ratio: Fraction of samples from promoter-enriched sampling
        promoter_tsv_path: Path to promoter TSV for reference
        shift_augmentation: Enable shift augmentation
        shift_range: Maximum shift range for augmentation
        shift_augmentation_factor: Fraction of samples to augment
        include_reverse_complements: If True, emit each window and its reverse
            complement as separate samples (approximately doubles the corpus)
        output_format: Output format ("txt" or "tsv")
        shuffle: Shuffle sequences before writing
        
    Returns:
        CorpusStats with generation statistics
    """
    set_global_seed(seed)
    LOGGER.info("Building pretrain corpus with seed %d", seed)
    
    # Read FASTA
    LOGGER.info("Reading FASTA file: %s", fasta_path)
    sequences_dict = read_fasta_sequences(fasta_path)
    
    if not sequences_dict:
        raise ValueError(f"No sequences found in {fasta_path}")
    
    LOGGER.info("Found %d sequence(s) in FASTA", len(sequences_dict))
    
    # Concatenate all sequences (for genomes with multiple chromosomes)
    contigs = list(sequences_dict.values())
    total_genome_len = sum(len(s) for s in contigs)
    LOGGER.info("Total genome length: %d bp", total_genome_len)
    
    all_windows: List[str] = []
    total_n_filtered = 0
    total_n_containing = 0
    total_n_bases = 0
    total_windows_considered = 0
    total_invalid_filtered = 0
    
    sliding_count = 0
    random_count = 0
    promoter_enriched_count = 0
    shift_augmented_count = 0
    reverse_complement_count = 0

    if sampling not in {"sliding", "random"}:
        raise ValueError(f"Unknown sampling mode: {sampling}")

    if sampling == "sliding":
        if promoter_enriched or shift_augmentation:
            LOGGER.warning(
                "promoter_enriched and shift_augmentation are ignored in sliding mode"
            )

        LOGGER.info("Extracting overlapping windows (window=%d, stride=%d)...", seq_len, stride)
        for contig in contigs:
            windows, n_filt, n_cont, n_total, total_win, invalid_filt = extract_overlapping_windows(
                contig,
                seq_len,
                stride,
                n_handling=n_handling,
                n_replacement=n_replacement,
                max_n_ratio=max_n_ratio,
            )
            all_windows.extend(windows)
            total_n_filtered += n_filt
            total_n_containing += n_cont
            total_n_bases += n_total
            total_windows_considered += total_win
            total_invalid_filtered += invalid_filt

        sliding_count = len(all_windows)
        LOGGER.info(
            "Generated %d overlapping windows (skipped %d / %d)",
            sliding_count,
            total_windows_considered - sliding_count,
            total_windows_considered,
        )

    else:
        # Random sampling
        if num_samples is None:
            num_samples = 100000

        # Calculate sample allocation
        if promoter_enriched:
            promoter_samples = int(num_samples * promoter_enriched_ratio)
            random_samples = num_samples - promoter_samples
        else:
            promoter_samples = 0
            random_samples = num_samples

        full_sequence = "".join(contigs)

        LOGGER.info("Sampling %d random windows...", random_samples)
        random_windows, n_filt, n_cont, n_total = sample_random_windows(
            full_sequence,
            seq_len,
            random_samples,
            n_handling=n_handling,
            n_replacement=n_replacement,
            max_n_ratio=max_n_ratio,
        )
        all_windows.extend(random_windows)
        random_count = len(random_windows)
        total_n_filtered += n_filt
        total_n_containing += n_cont
        total_n_bases += n_total
        LOGGER.info("Generated %d random samples", random_count)

        # Promoter-enriched sampling
        if promoter_enriched and promoter_samples > 0:
            LOGGER.info("Sampling %d promoter-enriched windows...", promoter_samples)
            promoter_windows = sample_promoter_enriched(
                full_sequence,
                seq_len,
                promoter_samples,
                n_handling=n_handling,
                n_replacement=n_replacement,
            )
            all_windows.extend(promoter_windows)
            promoter_enriched_count = len(promoter_windows)
            LOGGER.info("Generated %d promoter-enriched samples", promoter_enriched_count)

        # Shift augmentation
        if shift_augmentation:
            LOGGER.info("Applying shift augmentation (factor %.2f)...", shift_augmentation_factor)
            augmented = apply_shift_augmentation(
                all_windows,
                full_sequence,
                seq_len,
                shift_range=shift_range,
                augmentation_factor=shift_augmentation_factor,
                n_handling=n_handling,
                n_replacement=n_replacement,
            )
            all_windows.extend(augmented)
            shift_augmented_count = len(augmented)
            LOGGER.info("Generated %d shift-augmented samples", shift_augmented_count)

    if include_reverse_complements and all_windows:
        original_count = len(all_windows)
        rc_windows = [reverse_complement(window) for window in all_windows]
        reverse_complement_count = len(rc_windows)
        all_windows.extend(rc_windows)
        rc_fraction = reverse_complement_count / max(1, len(all_windows))
        LOGGER.info(
            "Reverse-complement augmentation: %d → %d windows (rc_frac=%.3f)",
            original_count,
            len(all_windows),
            rc_fraction,
        )
    elif include_reverse_complements:
        LOGGER.info("Reverse-complement augmentation: enabled (0 windows, rc_frac=0.000)")
    else:
        LOGGER.info("Reverse-complement augmentation: disabled (rc_frac=0.000)")
    
    # Compute stats
    stats = compute_corpus_stats(
        all_windows,
        total_n_filtered,
        total_n_containing,
        total_n_bases,
        random_count,
        promoter_enriched_count,
        shift_augmented_count,
    )
    stats.sliding_samples = sliding_count
    stats.total_windows_considered = total_windows_considered
    stats.invalid_filtered_windows = total_invalid_filtered
    stats.skipped_windows = max(0, total_windows_considered - sliding_count)
    stats.reverse_complement_samples = reverse_complement_count
    stats.reverse_complement_fraction = reverse_complement_count / max(1, stats.total_sequences)
    
    # Write output
    write_corpus(all_windows, output_path, format=output_format, shuffle=shuffle)
    
    return stats


def main() -> None:
    """Main entry point for corpus building."""
    parser = argparse.ArgumentParser(
        description="Build pretraining corpus from genomic FASTA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with config file
  python build_pretrain_corpus.py --config config.yaml

  # Sliding windows (large corpus)
  python build_pretrain_corpus.py --fasta ../data/genome.fna --sampling sliding --seq_len 81 --stride 20

  # Random sampling (legacy)
  python build_pretrain_corpus.py --fasta ../data/genome.fna --sampling random --seq_len 81 --num_samples 100000

  # With promoter enrichment
  python build_pretrain_corpus.py --fasta ../data/genome.fna --promoter_enriched --promoter_ratio 0.3

  # With shift augmentation
  python build_pretrain_corpus.py --fasta ../data/genome.fna --shift_augment --shift_range 15
        """,
    )
    
    # Input/output
    parser.add_argument("--config", type=str, default=None,
                        help="Path to YAML config file")
    parser.add_argument("--fasta", type=str, default=None,
                        help="Path to input FASTA file")
    parser.add_argument("--output", type=str, default=None,
                        help="Path to output corpus file")
    parser.add_argument("--output_format", type=str, default="txt",
                        choices=["txt", "tsv"],
                        help="Output format (default: txt)")
    
    # Sampling parameters
    parser.add_argument("--seq_len", type=int, default=None,
                        help="Sequence length / window size (default: from config)")
    parser.add_argument("--sampling", type=str, default="sliding",
                        choices=["sliding", "random"],
                        help="Sampling strategy (default: sliding)")
    parser.add_argument("--stride", type=int, default=None,
                        help="Stride in bp for sliding windows (default: from config or 20)")
    parser.add_argument("--num_samples", type=int, default=None,
                        help="Number of random samples to generate (default: 100000, random mode only)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed (default: from config or 1337)")
    
    # N-content handling
    parser.add_argument("--n_handling", type=str, default="filter",
                        choices=["filter", "replace", "skip"],
                        help="How to handle N bases (default: filter)")
    parser.add_argument("--n_replacement", type=str, default="random",
                        choices=["random", "A", "C", "G", "T"],
                        help="Replacement strategy for N (default: random)")
    parser.add_argument("--max_n_ratio", type=float, default=0.1,
                        help="Max N ratio before filtering (default: 0.1)")
    
    # Promoter enrichment
    parser.add_argument("--promoter_enriched", action="store_true",
                        help="Enable promoter-enriched sampling")
    parser.add_argument("--promoter_ratio", type=float, default=0.3,
                        help="Fraction from promoter-enriched sampling (default: 0.3)")
    parser.add_argument("--promoter_tsv", type=str, default=None,
                        help="Path to promoter TSV file")
    
    # Shift augmentation
    parser.add_argument("--shift_augment", action="store_true",
                        help="Enable shift augmentation")
    parser.add_argument("--shift_range", type=int, default=10,
                        help="Shift range in bp (default: 10)")
    parser.add_argument("--shift_factor", type=float, default=0.2,
                        help="Shift augmentation factor (default: 0.2)")
    
    # Other options
    parser.add_argument(
        "--include_reverse_complements",
        action="store_true",
        help="Duplicate corpus by adding reverse-complement window for every window.",
    )
    parser.add_argument("--no_shuffle", action="store_true",
                        help="Don't shuffle sequences before writing")
    
    args = parser.parse_args()
    
    script_dir = Path(__file__).resolve().parent
    
    # Load config if provided
    cfg: Dict[str, Any] = {}
    if args.config:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = script_dir / args.config
        cfg = load_config(config_path)
        LOGGER.info("Loaded config from %s", config_path)
    
    # Resolve parameters (CLI overrides config)
    fasta_path = args.fasta or cfg.get("mlm_fasta_path", "../data/GCF_000005845.2_ASM584v2_genomic.fna")
    if not Path(fasta_path).is_absolute():
        fasta_path = script_dir / fasta_path
    fasta_path = Path(fasta_path).resolve()
    
    seq_len = args.seq_len or cfg.get("mlm_window_size", cfg.get("max_bp_len", 81))
    seed = args.seed or cfg.get("seed", 1337)
    stride = args.stride if args.stride is not None else cfg.get("mlm_stride", 20)
    include_reverse_complements = bool(cfg.get("include_reverse_complements", False)) or bool(
        args.include_reverse_complements
    )
    
    # Output path
    if args.output:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            output_path = script_dir / args.output
    else:
        output_dir = script_dir.parent / "data"
        if args.sampling == "sliding":
            output_path = output_dir / f"pretrain_corpus_{seq_len}bp_stride{stride}.{args.output_format}"
        else:
            random_target = args.num_samples if args.num_samples is not None else 100000
            output_path = output_dir / f"pretrain_corpus_{seq_len}bp_random{random_target}.{args.output_format}"
    output_path = output_path.resolve()
    
    # Report path
    report_path = output_path.with_suffix(".report.json")
    
    # Promoter TSV
    promoter_tsv = None
    if args.promoter_tsv:
        promoter_tsv = Path(args.promoter_tsv)
    elif cfg.get("stage1_train"):
        promoter_tsv = script_dir / cfg["stage1_train"]
    
    LOGGER.info("=" * 60)
    LOGGER.info("PRETRAIN CORPUS BUILDER")
    LOGGER.info("=" * 60)
    LOGGER.info("FASTA: %s", fasta_path)
    LOGGER.info("Output: %s", output_path)
    LOGGER.info("Sequence length: %d bp", seq_len)
    LOGGER.info("Sampling: %s", args.sampling)
    if args.sampling == "sliding":
        LOGGER.info("Stride: %d bp", stride)
    else:
        random_target = args.num_samples if args.num_samples is not None else 100000
        LOGGER.info("Target samples: %d", random_target)
    LOGGER.info("Seed: %d", seed)
    LOGGER.info("N-handling: %s", args.n_handling)
    LOGGER.info("Promoter-enriched: %s", args.promoter_enriched)
    LOGGER.info("Shift augmentation: %s", args.shift_augment)
    LOGGER.info("Include reverse complements: %s", include_reverse_complements)
    LOGGER.info("=" * 60)
    
    # Build corpus
    stats = build_pretrain_corpus(
        fasta_path=fasta_path,
        output_path=output_path,
        seq_len=seq_len,
        num_samples=args.num_samples,
        seed=seed,
        n_handling=args.n_handling,
        n_replacement=args.n_replacement,
        max_n_ratio=args.max_n_ratio,
        sampling=args.sampling,
        stride=stride,
        promoter_enriched=args.promoter_enriched,
        promoter_enriched_ratio=args.promoter_ratio,
        promoter_tsv_path=promoter_tsv,
        shift_augmentation=args.shift_augment,
        shift_range=args.shift_range,
        shift_augmentation_factor=args.shift_factor,
        include_reverse_complements=include_reverse_complements,
        output_format=args.output_format,
        shuffle=not args.no_shuffle,
    )
    
    # Write report
    config_used = {
        "fasta_path": str(fasta_path),
        "output_path": str(output_path),
        "seq_len": seq_len,
        "sampling": args.sampling,
        "stride": stride,
        "num_samples": args.num_samples,
        "seed": seed,
        "n_handling": args.n_handling,
        "n_replacement": args.n_replacement,
        "max_n_ratio": args.max_n_ratio,
        "promoter_enriched": args.promoter_enriched,
        "promoter_enriched_ratio": args.promoter_ratio,
        "shift_augmentation": args.shift_augment,
        "shift_range": args.shift_range,
        "shift_augmentation_factor": args.shift_factor,
        "include_reverse_complements": include_reverse_complements,
    }
    write_report(stats, report_path, config_used)
    
    # Print summary
    print("\n" + "=" * 60)
    print("CORPUS GENERATION COMPLETE")
    print("=" * 60)
    print(f"Total sequences: {stats.total_sequences:,}")
    print(
        f"Length: {stats.length_min} - {stats.length_max} bp "
        f"(mean: {stats.length_mean:.1f}, std: {stats.length_std:.2f})"
    )
    if stats.total_windows_considered > 0:
        skipped_pct = (stats.skipped_windows / stats.total_windows_considered * 100) if stats.total_windows_considered else 0.0
        print(
            f"\nWindowing:"
            f"\n  Total windows considered: {stats.total_windows_considered:,}"
            f"\n  Skipped windows: {stats.skipped_windows:,} ({skipped_pct:.2f}%)"
        )
    print(f"\nSampling breakdown:")
    if stats.sliding_samples:
        print(f"  Sliding windows: {stats.sliding_samples:,}")
    print(f"  Random samples: {stats.random_samples:,}")
    print(f"  Promoter-enriched: {stats.promoter_enriched_samples:,}")
    print(f"  Shift-augmented: {stats.shift_augmented_samples:,}")
    print(f"\nN-content stats:")
    print(f"  Sequences with N (before sanitization): {stats.sequences_with_n:,}")
    print(f"  Total N bases encountered: {stats.total_n_bases:,}")
    print(f"  Sequences filtered due to N: {stats.n_filtered_sequences:,}")
    print(f"\nBase composition:")
    for base, pct in sorted(stats.compute_base_percentages().items()):
        print(f"  {base}: {pct:.2f}%")
    print(f"\nOutput: {output_path}")
    print(f"Report: {report_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
