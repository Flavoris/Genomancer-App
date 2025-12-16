#!/usr/bin/env python3
"""
Data sanity checker for Gene Whisperer training data.

Validates stage1_train.tsv and stage1_val.tsv for:
1. Required columns exist
2. Sequences contain only A/C/G/T after sanitization
3. Class balance
4. Duplicate sequence rate
5. Train/val leakage (identical sequences in both)
6. Spot-check 50 random rows
7. Token coverage stats (UNK rate)
"""

import argparse
import json
import random
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd
import yaml

# Constants
BASES_SET = {"A", "C", "G", "T"}
REQUIRED_COLUMNS = {"sequence", "is_promoter"}


def sanitize_sequence(seq: str) -> str:
    """Normalize DNA sequence: uppercase, keep only A/C/G/T."""
    seq = (seq or "").strip().upper()
    return "".join(base for base in seq if base in BASES_SET)


def load_tsv(path: Path) -> pd.DataFrame:
    """Load TSV file with error handling."""
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    df = pd.read_csv(path, sep="\t")
    return df


def check_required_columns(df: pd.DataFrame, name: str) -> List[str]:
    """Check that required columns exist."""
    issues = []
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        issues.append(f"{name}: Missing columns: {missing}")
    else:
        print(f"  ✓ {name}: All required columns present: {list(REQUIRED_COLUMNS)}")
    return issues


def check_sequences_valid(df: pd.DataFrame, name: str) -> Tuple[List[str], Dict[str, Any]]:
    """Check sequences contain only valid bases after sanitization."""
    issues = []
    stats = {
        "total": len(df),
        "empty_raw": 0,
        "empty_sanitized": 0,
        "has_non_acgt": 0,
        "lengths": [],
    }
    
    for idx, row in df.iterrows():
        raw = str(row.get("sequence", ""))
        sanitized = sanitize_sequence(raw)
        
        if not raw.strip():
            stats["empty_raw"] += 1
        
        if not sanitized:
            stats["empty_sanitized"] += 1
        
        # Check if raw had non-ACGT chars
        raw_upper = raw.strip().upper()
        if raw_upper and set(raw_upper) - BASES_SET:
            stats["has_non_acgt"] += 1
        
        stats["lengths"].append(len(sanitized))
    
    if stats["empty_raw"] > 0:
        issues.append(f"{name}: {stats['empty_raw']} empty raw sequences")
    
    if stats["empty_sanitized"] > 0:
        issues.append(f"{name}: {stats['empty_sanitized']} empty sequences after sanitization")
    
    # Length stats
    lengths = stats["lengths"]
    if lengths:
        min_len = min(lengths)
        max_len = max(lengths)
        avg_len = sum(lengths) / len(lengths)
        print(f"  ✓ {name}: Sequence lengths - min: {min_len}, max: {max_len}, avg: {avg_len:.1f}")
    
    if stats["has_non_acgt"] > 0:
        pct = stats["has_non_acgt"] / stats["total"] * 100
        print(f"  ⚠ {name}: {stats['has_non_acgt']} sequences ({pct:.1f}%) had non-ACGT chars (sanitized)")
    else:
        print(f"  ✓ {name}: All sequences contain only A/C/G/T")
    
    return issues, stats


def check_class_balance(df: pd.DataFrame, name: str) -> Dict[str, Any]:
    """Check class balance."""
    if "is_promoter" not in df.columns:
        return {"error": "is_promoter column missing"}
    
    labels = df["is_promoter"].astype(int)
    counts = Counter(labels)
    total = len(labels)
    
    positive = counts.get(1, 0)
    negative = counts.get(0, 0)
    
    pos_pct = positive / total * 100 if total else 0
    neg_pct = negative / total * 100 if total else 0
    
    balance_ratio = min(positive, negative) / max(positive, negative) if max(positive, negative) > 0 else 0
    
    print(f"  ✓ {name}: Class balance - Positive: {positive} ({pos_pct:.1f}%), Negative: {negative} ({neg_pct:.1f}%)")
    
    if balance_ratio < 0.5:
        print(f"  ⚠ {name}: Classes are imbalanced (ratio: {balance_ratio:.2f})")
    
    return {
        "positive": positive,
        "negative": negative,
        "positive_pct": pos_pct,
        "negative_pct": neg_pct,
        "balance_ratio": balance_ratio,
    }


def check_duplicates(df: pd.DataFrame, name: str) -> Tuple[List[str], Dict[str, Any]]:
    """Check for duplicate sequences."""
    issues = []
    
    # Use sanitized sequences for duplicate check
    sanitized = df["sequence"].apply(sanitize_sequence)
    
    total = len(sanitized)
    unique = sanitized.nunique()
    duplicates = total - unique
    dup_rate = duplicates / total * 100 if total else 0
    
    print(f"  {'⚠' if dup_rate > 5 else '✓'} {name}: {unique} unique / {total} total sequences ({duplicates} duplicates, {dup_rate:.2f}%)")
    
    if dup_rate > 10:
        issues.append(f"{name}: High duplicate rate ({dup_rate:.1f}%)")
    
    return issues, {"total": total, "unique": unique, "duplicates": duplicates, "duplicate_rate": dup_rate}


def check_leakage(train_df: pd.DataFrame, val_df: pd.DataFrame) -> Tuple[List[str], Dict[str, Any]]:
    """Check for train/val data leakage (identical sequences in both)."""
    issues = []
    
    train_seqs = set(train_df["sequence"].apply(sanitize_sequence))
    val_seqs = set(val_df["sequence"].apply(sanitize_sequence))
    
    overlap = train_seqs & val_seqs
    overlap_count = len(overlap)
    
    if overlap_count > 0:
        pct_of_val = overlap_count / len(val_seqs) * 100 if val_seqs else 0
        issues.append(f"DATA LEAKAGE: {overlap_count} sequences appear in both train and val ({pct_of_val:.1f}% of val)")
        print(f"  ✗ DATA LEAKAGE: {overlap_count} identical sequences in train and val ({pct_of_val:.1f}% of val)")
    else:
        print(f"  ✓ No data leakage: train and val sets are disjoint")
    
    return issues, {
        "train_unique": len(train_seqs),
        "val_unique": len(val_seqs),
        "overlap": overlap_count,
    }


def spot_check_samples(df: pd.DataFrame, name: str, n: int = 50) -> None:
    """Print spot-check of random samples."""
    print(f"\n{'='*70}")
    print(f"SPOT CHECK: {n} random samples from {name}")
    print(f"{'='*70}")
    
    sample_indices = random.sample(range(len(df)), min(n, len(df)))
    
    print(f"{'idx':>6} | {'label':>5} | {'len':>4} | {'raw_seq (first 40)':<42} | {'sanitized (first 40)':<40}")
    print("-" * 110)
    
    for idx in sample_indices:
        row = df.iloc[idx]
        raw_seq = str(row.get("sequence", ""))
        sanitized = sanitize_sequence(raw_seq)
        label = int(row.get("is_promoter", -1))
        length = len(sanitized)
        
        raw_display = raw_seq[:40] + "..." if len(raw_seq) > 40 else raw_seq
        san_display = sanitized[:40] + "..." if len(sanitized) > 40 else sanitized
        
        print(f"{idx:>6} | {label:>5} | {length:>4} | {raw_display:<42} | {san_display:<40}")


def compute_unk_rate(
    df: pd.DataFrame,
    name: str,
    k: int = 3,
    vocab_path: Optional[Path] = None,
) -> Tuple[List[str], Dict[str, Any]]:
    """Compute token coverage and UNK rate."""
    issues = []
    
    # Build vocab from data or load existing
    all_kmers: Set[str] = set()
    for seq in df["sequence"]:
        sanitized = sanitize_sequence(str(seq))
        for i in range(max(0, len(sanitized) - k + 1)):
            kmer = sanitized[i:i + k]
            if len(kmer) == k:
                all_kmers.add(kmer)
    
    # If vocab_path provided, load and check coverage
    vocab_kmers: Set[str] = set()
    if vocab_path and vocab_path.exists():
        with vocab_path.open("r") as f:
            vocab_data = json.load(f)
            # Vocab format: {"k": k, "itos": [...]} or {"k": k, "tokens": [...]}
            tokens = vocab_data.get("itos") or vocab_data.get("tokens") or []
            vocab_kmers = set(tokens)
            # Remove special tokens
            vocab_kmers -= {"<UNK>", "<PAD>", "[MASK]"}
        
        missing = all_kmers - vocab_kmers
        unk_count = len(missing)
        unk_rate = unk_count / len(all_kmers) * 100 if all_kmers else 0
        
        print(f"  {'✗' if unk_rate > 1 else '✓'} {name}: UNK rate for k={k}: {unk_count}/{len(all_kmers)} k-mers ({unk_rate:.2f}%) not in vocab")
        
        if unk_rate > 1:
            issues.append(f"{name}: High UNK rate ({unk_rate:.1f}%) for k={k}")
            # Show some missing k-mers
            sample_missing = list(missing)[:10]
            print(f"      Sample missing k-mers: {sample_missing}")
    else:
        # Just report observed k-mers
        print(f"  ℹ {name}: Observed {len(all_kmers)} unique {k}-mers")
        
        # Check theoretical coverage
        max_possible = 4 ** k
        coverage = len(all_kmers) / max_possible * 100
        print(f"  ℹ {name}: Coverage of theoretical {k}-mer space: {len(all_kmers)}/{max_possible} ({coverage:.1f}%)")
    
    return issues, {
        "observed_kmers": len(all_kmers),
        "vocab_kmers": len(vocab_kmers) if vocab_kmers else None,
        "missing_kmers": len(all_kmers - vocab_kmers) if vocab_kmers else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Data sanity checker for Gene Whisperer")
    parser.add_argument("--config", type=Path, default=Path("config.yaml"), help="Config file path")
    parser.add_argument("--train", type=Path, help="Override train TSV path")
    parser.add_argument("--val", type=Path, help="Override val TSV path")
    parser.add_argument("--k", type=int, default=3, help="K-mer size for UNK rate check")
    parser.add_argument("--vocab", type=Path, help="Vocab JSON path for UNK rate check")
    parser.add_argument("--samples", type=int, default=50, help="Number of spot-check samples")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
    args = parser.parse_args()
    
    random.seed(args.seed)
    
    script_dir = Path(__file__).resolve().parent
    
    # Load config
    config_path = args.config
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    
    cfg = {}
    if config_path.exists():
        with config_path.open("r") as f:
            cfg = yaml.safe_load(f) or {}
    
    # Determine data paths
    train_path = args.train
    val_path = args.val
    
    if not train_path:
        train_rel = cfg.get("stage1_train", "../data/stage1_train.tsv")
        train_path = (script_dir / train_rel).resolve()
    
    if not val_path:
        val_rel = cfg.get("stage1_val", "../data/stage1_val.tsv")
        val_path = (script_dir / val_rel).resolve()
    
    # Vocab path
    vocab_path = args.vocab
    if not vocab_path:
        vocab_dir = cfg.get("vocab_cache_dir", "../artifacts/vocabs")
        vocab_path = (script_dir / vocab_dir / f"k{args.k}_vocab.json").resolve()
    
    print("=" * 70)
    print("GENE WHISPERER DATA SANITY CHECK")
    print("=" * 70)
    print(f"Train file: {train_path}")
    print(f"Val file:   {val_path}")
    print(f"K-mer size: {args.k}")
    print(f"Vocab file: {vocab_path}")
    print("=" * 70)
    
    all_issues: List[str] = []
    
    # Load data
    print("\n[1/7] Loading data files...")
    try:
        train_df = load_tsv(train_path)
        print(f"  ✓ Loaded train: {len(train_df)} rows")
    except Exception as e:
        print(f"  ✗ Failed to load train: {e}")
        return 1
    
    try:
        val_df = load_tsv(val_path)
        print(f"  ✓ Loaded val: {len(val_df)} rows")
    except Exception as e:
        print(f"  ✗ Failed to load val: {e}")
        return 1
    
    # Check 1: Required columns
    print("\n[2/7] Checking required columns...")
    all_issues.extend(check_required_columns(train_df, "train"))
    all_issues.extend(check_required_columns(val_df, "val"))
    
    # Check 2: Valid sequences
    print("\n[3/7] Checking sequence validity...")
    issues, _ = check_sequences_valid(train_df, "train")
    all_issues.extend(issues)
    issues, _ = check_sequences_valid(val_df, "val")
    all_issues.extend(issues)
    
    # Check 3: Class balance
    print("\n[4/7] Checking class balance...")
    check_class_balance(train_df, "train")
    check_class_balance(val_df, "val")
    
    # Check 4: Duplicates
    print("\n[5/7] Checking for duplicates...")
    issues, _ = check_duplicates(train_df, "train")
    all_issues.extend(issues)
    issues, _ = check_duplicates(val_df, "val")
    all_issues.extend(issues)
    
    # Check 5: Train/val leakage
    print("\n[6/7] Checking for train/val leakage...")
    issues, _ = check_leakage(train_df, val_df)
    all_issues.extend(issues)
    
    # Check 6: Token coverage / UNK rate
    print("\n[7/7] Checking token coverage (UNK rate)...")
    issues, _ = compute_unk_rate(train_df, "train", k=args.k, vocab_path=vocab_path)
    all_issues.extend(issues)
    issues, _ = compute_unk_rate(val_df, "val", k=args.k, vocab_path=vocab_path)
    all_issues.extend(issues)
    
    # Spot check
    spot_check_samples(train_df, "train", n=args.samples)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    if all_issues:
        print(f"\n⚠ Found {len(all_issues)} issue(s):\n")
        for i, issue in enumerate(all_issues, 1):
            print(f"  {i}. {issue}")
        print()
        return 1
    else:
        print("\n✓ All checks passed! Data looks healthy.\n")
        return 0


if __name__ == "__main__":
    sys.exit(main())

