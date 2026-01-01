#!/usr/bin/env python3
"""
Convert PROCABLES raw files into TSV format for Gene Whisperer training.

Input files:
  - P-Non.txt → Stage 1 (binary promoter classification)
  - SW.txt → Stage 2 (strong/weak promoter strength)

PROCABLES record format:
  - Header line: >|<0 or 1>|training <metadata...>
  - Next line: 81bp DNA sequence

Output: 4 TSV files with stratified train/val splits.
"""

import argparse
import random
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_procables_file(filepath: Path, stage: str) -> list[dict]:
    """
    Parse a PROCABLES format file.

    Args:
        filepath: Path to the input file
        stage: 'stage1' or 'stage2'

    Returns:
        List of parsed records
    """
    records = []
    invalid_chars_count = 0
    length_warnings = 0

    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        header_line = lines[i].strip()

        if not header_line.startswith('>|'):
            i += 1
            continue

        # Parse header: >|<label>|training <metadata...>
        # Format: >|0|training ... or >|1|training ...
        header_match = re.match(r'^>\|([01])\|training\s*(.*)', header_line)
        if not header_match:
            print(f"Warning: Could not parse header at line {i+1}: {header_line}")
            i += 1
            continue

        label = int(header_match.group(1))
        metadata = header_match.group(2).strip()

        # Get sequence from next line
        if i + 1 >= len(lines):
            print(f"Warning: Missing sequence after header at line {i+1}")
            break

        sequence = lines[i + 1].strip().upper()

        # Validate sequence characters
        invalid_chars = set(sequence) - {'A', 'C', 'G', 'T'}
        if invalid_chars:
            invalid_chars_count += 1
            if invalid_chars_count <= 5:
                print(f"Warning: Invalid characters {invalid_chars} in sequence at line {i+2}")

        # Validate length
        if len(sequence) != 81:
            length_warnings += 1
            if length_warnings <= 5:
                print(f"Warning: Sequence length {len(sequence)} != 81 at line {i+2}")

        record = {
            'sequence': sequence,
            'metadata': metadata
        }

        if stage == 'stage1':
            record['is_promoter'] = label
        elif stage == 'stage2':
            # Extract strength from metadata
            # Patterns: "... SIGMA24 STRONG", "... SIGMA70STRONG", "... STRONG", etc.
            metadata_upper = metadata.upper()

            # Check for STRONG/WEAK anywhere in metadata (handles SIGMA70STRONG etc.)
            if 'STRONG' in metadata_upper:
                strength = 'strong'
                expected_label = 1
            elif 'WEAK' in metadata_upper:
                strength = 'weak'
                expected_label = 0
            else:
                # Fall back to numeric label (e.g., SIGMA38 without STRONG/WEAK)
                strength = 'strong' if label == 1 else 'weak'
                expected_label = label

            # Verify label matches strength (only warn if metadata had explicit strength)
            if ('STRONG' in metadata_upper or 'WEAK' in metadata_upper) and label != expected_label:
                print(f"Warning: Label {label} doesn't match strength in metadata at line {i+1}: {metadata}")

            record['strength'] = strength

        records.append(record)
        i += 2

    if invalid_chars_count > 5:
        print(f"... and {invalid_chars_count - 5} more sequences with invalid characters")
    if length_warnings > 5:
        print(f"... and {length_warnings - 5} more sequences with non-81bp length")

    return records


def stratified_split(records: list[dict], label_key: str, train_ratio: float, seed: int) -> tuple[list[dict], list[dict]]:
    """
    Perform stratified split to preserve class balance.

    Args:
        records: List of records to split
        label_key: Key to use for stratification ('is_promoter' or 'strength')
        train_ratio: Fraction for training set
        seed: Random seed for reproducibility

    Returns:
        Tuple of (train_records, val_records)
    """
    # Group by class
    by_class = defaultdict(list)
    for record in records:
        by_class[record[label_key]].append(record)

    train_records = []
    val_records = []

    rng = random.Random(seed)

    for class_label, class_records in by_class.items():
        # Shuffle deterministically
        shuffled = class_records.copy()
        rng.shuffle(shuffled)

        # Split
        n_train = int(len(shuffled) * train_ratio)
        train_records.extend(shuffled[:n_train])
        val_records.extend(shuffled[n_train:])

    # Shuffle final lists for good measure
    rng.shuffle(train_records)
    rng.shuffle(val_records)

    return train_records, val_records


def write_tsv(records: list[dict], filepath: Path, columns: list[str]):
    """Write records to a TSV file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')

        # Write records
        for record in records:
            row = [str(record[col]) for col in columns]
            f.write('\t'.join(row) + '\n')


def print_class_counts(name: str, records: list[dict], label_key: str):
    """Print class distribution."""
    counts = defaultdict(int)
    for record in records:
        counts[record[label_key]] += 1

    print(f"\n{name}:")
    for label, count in sorted(counts.items(), key=lambda x: str(x[0])):
        print(f"  {label}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert PROCABLES raw files to TSV format for Gene Whisperer training.'
    )
    parser.add_argument(
        '--pnon',
        type=Path,
        required=True,
        help='Path to P-Non.txt (Stage 1: promoter/non-promoter)'
    )
    parser.add_argument(
        '--sw',
        type=Path,
        required=True,
        help='Path to SW.txt (Stage 2: strong/weak promoters)'
    )
    parser.add_argument(
        '--out_dir',
        type=Path,
        default=Path('data'),
        help='Output directory for TSV files (default: data/)'
    )
    parser.add_argument(
        '--train_ratio',
        type=float,
        default=0.8,
        help='Fraction of data for training (default: 0.8)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=1337,
        help='Random seed for reproducibility (default: 1337)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.pnon.exists():
        print(f"Error: P-Non file not found: {args.pnon}")
        sys.exit(1)
    if not args.sw.exists():
        print(f"Error: SW file not found: {args.sw}")
        sys.exit(1)
    if not 0 < args.train_ratio < 1:
        print(f"Error: train_ratio must be between 0 and 1, got {args.train_ratio}")
        sys.exit(1)

    print(f"Converting PROCABLES datasets...")
    print(f"  P-Non: {args.pnon}")
    print(f"  SW: {args.sw}")
    print(f"  Output: {args.out_dir}")
    print(f"  Train ratio: {args.train_ratio}")
    print(f"  Seed: {args.seed}")

    # Parse Stage 1 (P-Non.txt)
    print(f"\n--- Stage 1: Parsing P-Non.txt ---")
    stage1_records = parse_procables_file(args.pnon, 'stage1')
    print(f"Parsed {len(stage1_records)} Stage 1 records")

    # Parse Stage 2 (SW.txt)
    print(f"\n--- Stage 2: Parsing SW.txt ---")
    stage2_records = parse_procables_file(args.sw, 'stage2')
    print(f"Parsed {len(stage2_records)} Stage 2 records")

    # Stratified splits
    print(f"\n--- Performing stratified splits ---")

    stage1_train, stage1_val = stratified_split(
        stage1_records, 'is_promoter', args.train_ratio, args.seed
    )
    stage2_train, stage2_val = stratified_split(
        stage2_records, 'strength', args.train_ratio, args.seed
    )

    # Print class distributions
    print("\n=== Stage 1 Class Distribution ===")
    print_class_counts("Training", stage1_train, 'is_promoter')
    print_class_counts("Validation", stage1_val, 'is_promoter')

    print("\n=== Stage 2 Class Distribution ===")
    print_class_counts("Training", stage2_train, 'strength')
    print_class_counts("Validation", stage2_val, 'strength')

    # Write TSV files
    print(f"\n--- Writing TSV files ---")

    stage1_columns = ['sequence', 'is_promoter']
    stage2_columns = ['sequence', 'strength']

    out_dir = args.out_dir

    write_tsv(stage1_train, out_dir / 'stage1_train.tsv', stage1_columns)
    print(f"  Wrote {len(stage1_train)} records to {out_dir / 'stage1_train.tsv'}")

    write_tsv(stage1_val, out_dir / 'stage1_val.tsv', stage1_columns)
    print(f"  Wrote {len(stage1_val)} records to {out_dir / 'stage1_val.tsv'}")

    write_tsv(stage2_train, out_dir / 'stage2_train.tsv', stage2_columns)
    print(f"  Wrote {len(stage2_train)} records to {out_dir / 'stage2_train.tsv'}")

    write_tsv(stage2_val, out_dir / 'stage2_val.tsv', stage2_columns)
    print(f"  Wrote {len(stage2_val)} records to {out_dir / 'stage2_val.tsv'}")

    print("\nDone!")


if __name__ == '__main__':
    main()
