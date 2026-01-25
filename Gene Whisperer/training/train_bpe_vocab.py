"""Train BPE vocabulary on FASTA genomic data.

Standalone script that trains a DNA-specific BPE vocabulary (DNABERT-2 style)
from genomic FASTA files. Must be run before any model training.

Usage:
    python train_bpe_vocab.py --config config.yaml
    python train_bpe_vocab.py --fasta data.fna --output vocab.json
"""

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import List, Optional

from bpe_tokenizer import DNABPETokenizer

LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)


def parse_fasta(path: Path) -> List[str]:
    """Parse FASTA file(s) and return list of sequence strings.

    Handles both single FASTA files and directories of FASTA files.
    """
    sequences: List[str] = []
    fasta_extensions = {".fna", ".fa", ".fasta", ".fsa"}

    if path.is_dir():
        fasta_files = [
            f for f in sorted(path.iterdir())
            if f.suffix.lower() in fasta_extensions
        ]
        if not fasta_files:
            LOGGER.warning("No FASTA files found in directory: %s", path)
            return sequences
        for fasta_file in fasta_files:
            sequences.extend(_parse_single_fasta(fasta_file))
    else:
        sequences.extend(_parse_single_fasta(path))

    return sequences


def _parse_single_fasta(path: Path) -> List[str]:
    """Parse a single FASTA file into sequence records."""
    sequences: List[str] = []
    current_seq: List[str] = []

    LOGGER.info("Reading FASTA: %s", path)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.upper())

    if current_seq:
        sequences.append("".join(current_seq))

    LOGGER.info("  Found %d sequence record(s)", len(sequences))
    return sequences


def extract_windows(
    sequences: List[str],
    window_size: int = 234,
    stride: int = 20,
    max_windows: Optional[int] = None,
) -> List[str]:
    """Extract fixed-size windows from sequences for BPE training.

    Args:
        sequences: List of full genomic sequences.
        window_size: Size of each window in base pairs.
        stride: Step between windows.
        max_windows: Maximum number of windows to extract (None for all).

    Returns:
        List of DNA window strings.
    """
    windows: List[str] = []
    dna_pattern = re.compile(r"[^ATCG]")

    for seq in sequences:
        cleaned = dna_pattern.sub("", seq)
        if len(cleaned) < window_size:
            if len(cleaned) > 0:
                windows.append(cleaned)
            continue

        for start in range(0, len(cleaned) - window_size + 1, stride):
            window = cleaned[start:start + window_size]
            windows.append(window)

            if max_windows and len(windows) >= max_windows:
                return windows

    return windows


def train_bpe_vocab(
    fasta_paths: List[str],
    output_path: str,
    vocab_size: int = 4096,
    window_size: int = 234,
    stride: int = 20,
    min_frequency: int = 2,
    max_windows: Optional[int] = None,
) -> DNABPETokenizer:
    """Train BPE vocabulary from FASTA files.

    Args:
        fasta_paths: List of paths to FASTA files or directories.
        output_path: Where to save the trained vocabulary JSON.
        vocab_size: Target vocabulary size (default 4096 per DNABERT-2).
        window_size: Window size for extracting training sequences.
        stride: Stride between windows.
        min_frequency: Minimum pair frequency for BPE merging.
        max_windows: Maximum training windows (None for unlimited).

    Returns:
        Trained DNABPETokenizer instance.
    """
    # Collect sequences from all FASTA sources
    all_sequences: List[str] = []
    for fasta_path in fasta_paths:
        path = Path(fasta_path)
        if not path.exists():
            LOGGER.warning("FASTA path not found, skipping: %s", fasta_path)
            continue
        seqs = parse_fasta(path)
        all_sequences.extend(seqs)

    if not all_sequences:
        raise ValueError("No sequences found in any FASTA path")

    LOGGER.info("Total genomic records: %d", len(all_sequences))
    total_bp = sum(len(s) for s in all_sequences)
    LOGGER.info("Total base pairs: %s", f"{total_bp:,}")

    # Extract training windows
    LOGGER.info(
        "Extracting windows (size=%d, stride=%d)...",
        window_size, stride,
    )
    windows = extract_windows(
        all_sequences,
        window_size=window_size,
        stride=stride,
        max_windows=max_windows,
    )
    LOGGER.info("Training windows: %d", len(windows))

    # Train BPE tokenizer
    LOGGER.info("Training BPE tokenizer (vocab_size=%d)...", vocab_size)
    tokenizer = DNABPETokenizer(vocab_size=vocab_size)
    tokenizer.train(windows, min_frequency=min_frequency)
    LOGGER.info("Final vocabulary size: %d", len(tokenizer))
    LOGGER.info("Merge operations learned: %d", len(tokenizer.merges))

    # Save vocabulary
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tokenizer.save(str(output))
    LOGGER.info("Saved BPE vocabulary to: %s", output)

    # Print compression statistics
    _print_compression_stats(tokenizer, windows)

    return tokenizer


def _print_compression_stats(
    tokenizer: DNABPETokenizer,
    windows: List[str],
    sample_size: int = 1000,
) -> None:
    """Print compression and token length statistics."""
    import random as rng

    sample = windows[:sample_size] if len(windows) <= sample_size else rng.sample(windows, sample_size)

    ratios = [tokenizer.compression_ratio(w) for w in sample]
    avg_ratio = sum(ratios) / len(ratios) if ratios else 0.0
    max_tokens = max(len(tokenizer.tokenize(w)) for w in sample)
    avg_tokens = sum(len(tokenizer.tokenize(w)) for w in sample) / len(sample)

    LOGGER.info("--- Compression Statistics (sample=%d) ---", len(sample))
    LOGGER.info("  Average compression ratio: %.2fx", avg_ratio)
    LOGGER.info("  Average tokens per window: %.1f", avg_tokens)
    LOGGER.info("  Max tokens per window: %d", max_tokens)

    # Token length distribution
    lengths = tokenizer.get_token_lengths()
    LOGGER.info("--- Token Length Distribution ---")
    for length in sorted(lengths.keys()):
        LOGGER.info("  %d-mer tokens: %d", length, lengths[length])

    # Recommendations for max_token_len
    LOGGER.info("--- Recommended max_token_len Settings ---")
    window_size = len(sample[0]) if sample else 234
    LOGGER.info(
        "  For %dbp windows (MLM): %d (max_tokens + 10%% headroom)",
        window_size,
        int(max_tokens * 1.1),
    )
    # For 81bp fine-tuning sequences
    short_sample = [w[:81] for w in sample[:200]]
    short_max = max(len(tokenizer.tokenize(s)) for s in short_sample)
    LOGGER.info(
        "  For 81bp sequences (fine-tuning): %d (max_tokens + 10%% headroom)",
        int(short_max * 1.1),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train BPE vocabulary on FASTA genomic data"
    )
    parser.add_argument(
        "--config", type=str, default=None,
        help="Path to config.yaml (uses mlm_fasta_paths, bpe_vocab_path, etc.)",
    )
    parser.add_argument(
        "--fasta", type=str, nargs="+", default=None,
        help="FASTA file(s) or directory(ies) to train on",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output path for BPE vocabulary JSON",
    )
    parser.add_argument(
        "--vocab-size", type=int, default=4096,
        help="Target vocabulary size (default: 4096)",
    )
    parser.add_argument(
        "--window-size", type=int, default=234,
        help="Window size in bp (default: 234)",
    )
    parser.add_argument(
        "--stride", type=int, default=20,
        help="Stride between windows (default: 20)",
    )
    parser.add_argument(
        "--min-frequency", type=int, default=2,
        help="Minimum pair frequency for merging (default: 2)",
    )
    parser.add_argument(
        "--max-windows", type=int, default=None,
        help="Max number of training windows (default: unlimited)",
    )
    args = parser.parse_args()

    # Determine FASTA paths and output
    fasta_paths: List[str] = []
    output_path: str = "../artifacts/vocabs/bpe_vocab.json"

    if args.config:
        try:
            import yaml
        except ImportError:
            sys.exit("PyYAML required for --config mode (pip install pyyaml)")

        config_path = Path(args.config)
        if not config_path.exists():
            sys.exit(f"Config file not found: {args.config}")

        # Change to config directory for relative paths
        os.chdir(config_path.parent)

        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f) or {}

        # Get FASTA paths from config
        mlm_fasta_paths = cfg.get("mlm_fasta_paths", [])
        if mlm_fasta_paths:
            fasta_paths = mlm_fasta_paths
        elif cfg.get("mlm_fasta_path"):
            fasta_paths = [cfg["mlm_fasta_path"]]

        # Get output path from config
        output_path = cfg.get("bpe_vocab_path", output_path)

        # Override defaults from config
        if args.vocab_size == 4096:
            args.vocab_size = int(cfg.get("vocab_size", 4096))
        if args.window_size == 234:
            args.window_size = int(cfg.get("mlm_window_size", 234))
        if args.stride == 20:
            args.stride = int(cfg.get("mlm_stride", 20))

    if args.fasta:
        fasta_paths = args.fasta
    if args.output:
        output_path = args.output

    if not fasta_paths:
        sys.exit(
            "No FASTA paths specified. Use --fasta or --config with mlm_fasta_paths."
        )

    LOGGER.info("BPE Vocabulary Training")
    LOGGER.info("  Vocab size: %d", args.vocab_size)
    LOGGER.info("  Window size: %d bp", args.window_size)
    LOGGER.info("  Stride: %d", args.stride)
    LOGGER.info("  Min frequency: %d", args.min_frequency)
    LOGGER.info("  FASTA paths: %s", fasta_paths)
    LOGGER.info("  Output: %s", output_path)

    train_bpe_vocab(
        fasta_paths=fasta_paths,
        output_path=output_path,
        vocab_size=args.vocab_size,
        window_size=args.window_size,
        stride=args.stride,
        min_frequency=args.min_frequency,
        max_windows=args.max_windows,
    )

    LOGGER.info("Done.")


if __name__ == "__main__":
    main()
