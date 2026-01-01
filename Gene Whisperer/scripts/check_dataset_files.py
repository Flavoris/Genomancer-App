#!/usr/bin/env python3
"""
Dataset Sanity Check Script for Gene Whisperer.

Validates the 4 TSV files used for training:
- Stage 1: sequence + is_promoter (binary 0/1)
- Stage 2: sequence + strength (strong/weak)

Prints statistics and exits non-zero if any validation fails.
"""

import sys
from pathlib import Path

import yaml
import pandas as pd


def load_config() -> dict:
    """Load config.yaml from the training directory."""
    config_path = Path(__file__).parent.parent / "training" / "config.yaml"
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        return yaml.safe_load(f)


def resolve_path(relative_path: str, base_dir: Path) -> Path:
    """Resolve a path relative to the training directory."""
    return (base_dir / relative_path).resolve()


def check_stage1_file(filepath: Path, delimiter: str, split_name: str) -> tuple[bool, dict]:
    """
    Validate a Stage 1 TSV file.

    Returns:
        (is_valid, stats_dict)
    """
    errors = []
    stats = {"split": split_name, "file": filepath.name}

    if not filepath.exists():
        print(f"  ERROR: File not found: {filepath}")
        return False, stats

    # Load file
    try:
        df = pd.read_csv(filepath, sep=delimiter)
    except Exception as e:
        print(f"  ERROR: Failed to load {filepath}: {e}")
        return False, stats

    stats["rows"] = len(df)

    # Check required columns
    required_cols = {"sequence", "is_promoter"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        errors.append(f"Missing columns: {missing}")

    if "is_promoter" in df.columns:
        # Check label values are only 0 or 1
        unique_labels = set(df["is_promoter"].unique())
        valid_labels = {0, 1}
        if not unique_labels.issubset(valid_labels):
            invalid = unique_labels - valid_labels
            errors.append(f"Invalid is_promoter values: {invalid} (expected only 0, 1)")

        # Class balance
        label_counts = df["is_promoter"].value_counts().to_dict()
        stats["class_0"] = label_counts.get(0, 0)
        stats["class_1"] = label_counts.get(1, 0)
        total = stats["class_0"] + stats["class_1"]
        if total > 0:
            stats["balance"] = f"{stats['class_1']}/{total} = {stats['class_1']/total:.1%} promoter"

    if "sequence" in df.columns:
        # Sequence length stats
        seq_lens = df["sequence"].str.len()
        stats["min_len"] = int(seq_lens.min())
        stats["mean_len"] = float(seq_lens.mean())
        stats["max_len"] = int(seq_lens.max())

        # Check for empty sequences
        empty_count = (seq_lens == 0).sum()
        if empty_count > 0:
            errors.append(f"{empty_count} empty sequences found")

    # Report errors
    for err in errors:
        print(f"  ERROR: {err}")

    return len(errors) == 0, stats


def check_stage2_file(filepath: Path, delimiter: str, split_name: str) -> tuple[bool, dict]:
    """
    Validate a Stage 2 TSV file.

    Returns:
        (is_valid, stats_dict)
    """
    errors = []
    stats = {"split": split_name, "file": filepath.name}

    if not filepath.exists():
        print(f"  ERROR: File not found: {filepath}")
        return False, stats

    # Load file
    try:
        df = pd.read_csv(filepath, sep=delimiter)
    except Exception as e:
        print(f"  ERROR: Failed to load {filepath}: {e}")
        return False, stats

    stats["rows"] = len(df)

    # Check required columns
    required_cols = {"sequence", "strength"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        errors.append(f"Missing columns: {missing}")

    if "strength" in df.columns:
        # Check strength values map to strong/weak
        unique_strengths = set(df["strength"].unique())
        valid_strengths = {"strong", "weak"}
        if not unique_strengths.issubset(valid_strengths):
            invalid = unique_strengths - valid_strengths
            errors.append(f"Invalid strength values: {invalid} (expected only 'strong', 'weak')")

        # Class balance
        strength_counts = df["strength"].value_counts().to_dict()
        stats["class_weak"] = strength_counts.get("weak", 0)
        stats["class_strong"] = strength_counts.get("strong", 0)
        total = stats["class_weak"] + stats["class_strong"]
        if total > 0:
            stats["balance"] = f"{stats['class_strong']}/{total} = {stats['class_strong']/total:.1%} strong"

    if "sequence" in df.columns:
        # Sequence length stats
        seq_lens = df["sequence"].str.len()
        stats["min_len"] = int(seq_lens.min())
        stats["mean_len"] = float(seq_lens.mean())
        stats["max_len"] = int(seq_lens.max())

        # Check for empty sequences
        empty_count = (seq_lens == 0).sum()
        if empty_count > 0:
            errors.append(f"{empty_count} empty sequences found")

    # Report errors
    for err in errors:
        print(f"  ERROR: {err}")

    return len(errors) == 0, stats


def print_stats_table(all_stats: list[dict], stage: int) -> None:
    """Print a formatted table of statistics."""
    print(f"\n{'='*70}")
    print(f"Stage {stage} Dataset Statistics")
    print(f"{'='*70}")

    for stats in all_stats:
        print(f"\n  {stats.get('split', 'Unknown')} ({stats.get('file', 'N/A')})")
        print(f"    Rows:           {stats.get('rows', 'N/A'):,}")
        print(f"    Class balance:  {stats.get('balance', 'N/A')}")
        print(f"    Sequence length:")
        print(f"      Min:  {stats.get('min_len', 'N/A')}")
        print(f"      Mean: {stats.get('mean_len', 'N/A'):.1f}" if isinstance(stats.get('mean_len'), float) else f"      Mean: N/A")
        print(f"      Max:  {stats.get('max_len', 'N/A')}")


def main() -> int:
    """Run all dataset checks. Returns exit code (0=success, 1=failure)."""
    print("Gene Whisperer Dataset Sanity Check")
    print("="*40)

    config = load_config()
    delimiter = config.get("delimiter", "\t")
    training_dir = Path(__file__).parent.parent / "training"

    all_valid = True
    stage1_stats = []
    stage2_stats = []

    # Stage 1 files
    print("\n[Stage 1] Checking binary classification files...")
    stage1_files = [
        ("stage1_train", config.get("stage1_train", "../data/stage1_train.tsv")),
        ("stage1_val", config.get("stage1_val", "../data/stage1_val.tsv")),
    ]

    for split_name, rel_path in stage1_files:
        filepath = resolve_path(rel_path, training_dir)
        print(f"  Checking {filepath.name}...")
        valid, stats = check_stage1_file(filepath, delimiter, split_name)
        stage1_stats.append(stats)
        if not valid:
            all_valid = False

    # Stage 2 files
    print("\n[Stage 2] Checking strength classification files...")
    stage2_files = [
        ("stage2_train", config.get("stage2_train", "../data/stage2_train.tsv")),
        ("stage2_val", config.get("stage2_val", "../data/stage2_val.tsv")),
    ]

    for split_name, rel_path in stage2_files:
        filepath = resolve_path(rel_path, training_dir)
        print(f"  Checking {filepath.name}...")
        valid, stats = check_stage2_file(filepath, delimiter, split_name)
        stage2_stats.append(stats)
        if not valid:
            all_valid = False

    # Print summary tables
    print_stats_table(stage1_stats, stage=1)
    print_stats_table(stage2_stats, stage=2)

    # Final result
    print(f"\n{'='*70}")
    if all_valid:
        print("All checks PASSED")
        return 0
    else:
        print("Some checks FAILED - see errors above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
