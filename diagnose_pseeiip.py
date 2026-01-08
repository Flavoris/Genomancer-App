from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


def _load_compute_pseeiip():
    project_root = Path(__file__).resolve().parent
    training_dir = project_root / "Gene Whisperer" / "training"
    if not training_dir.exists():
        raise SystemExit(f"Training directory not found: {training_dir}")
    sys.path.insert(0, str(training_dir))
    from dataset import compute_pseeiip

    return compute_pseeiip


def _format_nonzero_preview(features: np.ndarray, limit: int = 10) -> str:
    nonzero_indices = np.flatnonzero(features)
    if nonzero_indices.size == 0:
        return "  (none)"
    lines = []
    for idx in nonzero_indices[:limit]:
        lines.append(f"  idx {idx:02d}: {features[idx]:.6f}")
    return "\n".join(lines)


def main() -> None:
    compute_pseeiip = _load_compute_pseeiip()
    test_seq = (
        "TTGACAATTTTTCTTGATAATGTAACTCACTTAATCTTGATAAATGCTATAATGTGTCG"
        "AAAAAAAAAAAAAAAAAAAAA"
    )
    features = compute_pseeiip(test_seq).detach().cpu().numpy()

    nonzero_indices = np.flatnonzero(features)
    print("PseEIIP Diagnostic Results:\n")
    print(f"Feature dimension: {features.shape[0]}")
    print(f"Non-zero elements: {nonzero_indices.size}")
    print(f"Sum of values: {features.sum():.6f}")
    print(f"Min value: {features.min():.6f}")
    print(f"Max value: {features.max():.6f}")
    print("First 10 non-zero indices and values:")
    print(_format_nonzero_preview(features))


if __name__ == "__main__":
    main()
