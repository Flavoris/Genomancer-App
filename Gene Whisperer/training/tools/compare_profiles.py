#!/usr/bin/env python3
"""Compare two profiling JSON reports and print before vs after summary.

Usage:
    python tools/compare_profiles.py <old_report.json> <new_report.json>

Prints:
    - dataloader ms
    - forward+backward ms
    - optimizer step ms
    - total step ms
    - % improvement
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load_report(path: Path) -> dict:
    """Load a profile report JSON file."""
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def format_change(old_val: float, new_val: float) -> str:
    """Format the change as improvement or regression."""
    if old_val == 0:
        return "N/A"
    pct_change = ((old_val - new_val) / old_val) * 100
    if pct_change > 0:
        return f"\033[32m-{pct_change:.1f}%\033[0m"  # Green for improvement
    elif pct_change < 0:
        return f"\033[31m+{abs(pct_change):.1f}%\033[0m"  # Red for regression
    else:
        return "0.0%"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare two profiling reports (before vs after)"
    )
    parser.add_argument("old_report", type=Path, help="Path to old/baseline report JSON")
    parser.add_argument("new_report", type=Path, help="Path to new/optimized report JSON")
    args = parser.parse_args()

    if not args.old_report.exists():
        print(f"Error: Old report not found: {args.old_report}", file=sys.stderr)
        return 1
    if not args.new_report.exists():
        print(f"Error: New report not found: {args.new_report}", file=sys.stderr)
        return 1

    old = load_report(args.old_report)
    new = load_report(args.new_report)

    # Extract timing values
    old_dataloader = old.get("avg_dataloader_time_ms", 0.0)
    new_dataloader = new.get("avg_dataloader_time_ms", 0.0)

    old_fwd_bwd = old.get("avg_forward_backward_time_ms", 0.0)
    new_fwd_bwd = new.get("avg_forward_backward_time_ms", 0.0)

    old_opt = old.get("avg_optimizer_step_time_ms", 0.0)
    new_opt = new.get("avg_optimizer_step_time_ms", 0.0)

    # Total step time (dataloader + forward+backward + optimizer_step)
    # Note: optimizer step happens less frequently (every grad_accum_steps),
    # so we include it directly for simplicity
    old_total = old_dataloader + old_fwd_bwd + old_opt
    new_total = new_dataloader + new_fwd_bwd + new_opt

    # Print comparison
    print("=" * 70)
    print("PROFILING COMPARISON: Before vs After")
    print("=" * 70)
    print(f"Old report: {args.old_report}")
    print(f"New report: {args.new_report}")
    print("-" * 70)
    print(f"{'Metric':<30} {'Before':>12} {'After':>12} {'Change':>12}")
    print("-" * 70)

    print(f"{'Dataloader (ms)':<30} {old_dataloader:>12.2f} {new_dataloader:>12.2f} {format_change(old_dataloader, new_dataloader):>12}")
    print(f"{'Forward+Backward (ms)':<30} {old_fwd_bwd:>12.2f} {new_fwd_bwd:>12.2f} {format_change(old_fwd_bwd, new_fwd_bwd):>12}")
    print(f"{'Optimizer Step (ms)':<30} {old_opt:>12.2f} {new_opt:>12.2f} {format_change(old_opt, new_opt):>12}")
    print("-" * 70)
    print(f"{'Total Step Time (ms)':<30} {old_total:>12.2f} {new_total:>12.2f} {format_change(old_total, new_total):>12}")
    print("=" * 70)

    # Overall improvement summary
    if old_total > 0:
        overall_improvement = ((old_total - new_total) / old_total) * 100
        if overall_improvement > 0:
            print(f"\n\033[32mOverall improvement: {overall_improvement:.1f}% faster\033[0m")
        elif overall_improvement < 0:
            print(f"\n\033[31mOverall regression: {abs(overall_improvement):.1f}% slower\033[0m")
        else:
            print("\nNo change in total step time")

    return 0


if __name__ == "__main__":
    sys.exit(main())
