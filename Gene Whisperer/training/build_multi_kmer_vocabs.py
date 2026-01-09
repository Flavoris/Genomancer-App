"""
Build vocabularies for all k-mer sizes.

Each k-mer vocabulary contains:
- All possible k-mers (4^k)
- [MASK] token for MLM
- <UNK> token for unknown k-mers
- <PAD> token for padding
"""

from __future__ import annotations

import json
import logging
from itertools import product
from pathlib import Path
from typing import List

LOGGER = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

BASES = ["A", "C", "G", "T"]


def build_kmer_vocabulary(k: int) -> dict:
    """
    Build complete k-mer vocabulary.

    Args:
        k: K-mer size

    Returns:
        Vocabulary dict with itos, stoi, and special token info
    """
    # Generate all possible k-mers
    all_kmers = ["".join(kmer) for kmer in product(BASES, repeat=k)]

    # Sort for consistency
    all_kmers = sorted(all_kmers)

    # Add special tokens at the end
    special_tokens = ["[MASK]", "<UNK>", "<PAD>"]

    itos = all_kmers + special_tokens
    stoi = {token: idx for idx, token in enumerate(itos)}

    vocab = {
        "k": k,
        "itos": itos,
        "unk_token": "<UNK>",
        "pad_token": "<PAD>",
        "mask_token": "[MASK]",
        "vocab_size": len(itos),
        "num_kmers": len(all_kmers),
        "num_special_tokens": len(special_tokens),
    }

    return vocab


def save_vocabulary(vocab: dict, path: Path) -> None:
    """Save vocabulary to JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(vocab, f, indent=2)
    LOGGER.info(f"Saved k={vocab['k']} vocabulary to {path}")


def build_all_vocabularies(
    kmers: List[int] = [3, 4, 5, 6],
    output_dir: Path = Path("../artifacts/vocabs"),
) -> dict:
    """
    Build vocabularies for all k-mer sizes.

    Args:
        kmers: List of k-mer sizes
        output_dir: Directory to save vocabularies

    Returns:
        Dict mapping k to vocab path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vocab_paths = {}

    print("=" * 60)
    print("BUILDING K-MER VOCABULARIES")
    print("=" * 60)

    for k in kmers:
        vocab = build_kmer_vocabulary(k)

        # Save as MLM vocab
        vocab_path = output_dir / f"mlm_k{k}_vocab.json"
        save_vocabulary(vocab, vocab_path)

        # Also save as general vocab for stage training
        general_vocab_path = output_dir / f"k{k}_vocab.json"
        save_vocabulary(vocab, general_vocab_path)

        vocab_paths[k] = str(vocab_path)

        print(f"k={k}: {vocab['vocab_size']} tokens ({vocab['num_kmers']} k-mers + {vocab['num_special_tokens']} special)")

    print("=" * 60)
    print("ALL VOCABULARIES BUILT")
    print("=" * 60)

    return vocab_paths


def verify_vocabularies(
    kmers: List[int] = [3, 4, 5, 6],
    vocab_dir: Path = Path("../artifacts/vocabs"),
) -> bool:
    """
    Verify all vocabularies are correctly built.

    Returns:
        True if all verifications pass
    """
    vocab_dir = Path(vocab_dir)
    all_passed = True

    print("\nVERIFYING VOCABULARIES")
    print("-" * 40)

    for k in kmers:
        vocab_path = vocab_dir / f"mlm_k{k}_vocab.json"

        if not vocab_path.exists():
            print(f"✗ k={k}: Vocabulary file not found")
            all_passed = False
            continue

        with vocab_path.open() as f:
            vocab = json.load(f)

        # Check vocab size
        expected_size = (4 ** k) + 3  # 4^k k-mers + 3 special
        actual_size = vocab["vocab_size"]

        if actual_size != expected_size:
            print(f"✗ k={k}: Wrong vocab size (expected {expected_size}, got {actual_size})")
            all_passed = False
            continue

        # Check special tokens
        itos = vocab["itos"]
        has_mask = "[MASK]" in itos
        has_unk = "<UNK>" in itos
        has_pad = "<PAD>" in itos

        if not (has_mask and has_unk and has_pad):
            print(f"✗ k={k}: Missing special tokens")
            all_passed = False
            continue

        # Check some k-mers exist
        expected_kmers = ["".join(["A"] * k), "".join(["T"] * k), "".join(["C"] * k)]
        all_exist = all(km in itos for km in expected_kmers)

        if not all_exist:
            print(f"✗ k={k}: Some expected k-mers missing")
            all_passed = False
            continue

        print(f"✓ k={k}: {actual_size} tokens (correct)")

    print("-" * 40)
    if all_passed:
        print("ALL VOCABULARIES VERIFIED ✓")
    else:
        print("SOME VERIFICATIONS FAILED ✗")

    return all_passed


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build multi-k-mer vocabularies")
    parser.add_argument("--output-dir", type=str, default="../artifacts/vocabs")
    parser.add_argument("--kmers", type=int, nargs="+", default=[3, 4, 5, 6])

    args = parser.parse_args()

    build_all_vocabularies(
        kmers=args.kmers,
        output_dir=Path(args.output_dir),
    )

    verify_vocabularies(
        kmers=args.kmers,
        vocab_dir=Path(args.output_dir),
    )
