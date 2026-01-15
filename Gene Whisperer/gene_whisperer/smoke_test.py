#!/usr/bin/env python3
"""Smoke test for the Stage 1 ensemble predictor.

Runs a single sequence through the k=4 + k=6 ensemble and prints JSON output.

Usage:
    python -m gene_whisperer.smoke_test
    # or
    python gene_whisperer/smoke_test.py
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)
LOGGER = logging.getLogger("gene_whisperer.smoke_test")

# Project paths
_PACKAGE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _PACKAGE_DIR.parent


def run_smoke_test() -> int:
    """
    Run smoke test: predict one sequence using the ensemble.

    Returns:
        0 on success, 1 on failure
    """
    from gene_whisperer.predict import load_stage1_ensemble

    # Test sequence (E. coli promoter consensus region)
    # This is a known promoter-like sequence with TATA box and -10/-35 elements
    test_sequence = (
        "TTGACA"  # -35 element
        "NNNNNNNNNNNNNNNNNN"  # spacer
        "TATAAT"  # -10 element (TATA box)
        "NNNNNN"  # spacer to TSS
        "ATGCGTAACTGCGTAACTGCGTAACTGCGTAACTGCGTAACTGCGTAACT"  # downstream
    ).replace("N", "A")  # Replace N with A for simplicity

    LOGGER.info("Test sequence length: %d bp", len(test_sequence))
    LOGGER.info("Test sequence: %s...%s", test_sequence[:20], test_sequence[-20:])

    # Paths to artifacts
    config_path = _PROJECT_ROOT / "training" / "config.yaml"
    artifacts_dir = _PROJECT_ROOT / "artifacts"

    # Check for available checkpoints
    kmer_config = {}

    # Try k=4
    k4_checkpoint = artifacts_dir / "checkpoints" / "stage1_k4.pt"
    k4_vocab = artifacts_dir / "vocabs" / "k4_vocab.json"
    if k4_checkpoint.exists() and k4_vocab.exists():
        kmer_config[4] = {
            "checkpoint": str(k4_checkpoint),
            "vocab": str(k4_vocab),
        }
        LOGGER.info("Found k=4 model")
    else:
        LOGGER.warning("k=4 model not found (checkpoint=%s, vocab=%s)", k4_checkpoint.exists(), k4_vocab.exists())

    # Try k=6
    k6_checkpoint = artifacts_dir / "checkpoints" / "stage1_k6.pt"
    k6_vocab = artifacts_dir / "vocabs" / "k6_vocab.json"
    if k6_checkpoint.exists() and k6_vocab.exists():
        kmer_config[6] = {
            "checkpoint": str(k6_checkpoint),
            "vocab": str(k6_vocab),
        }
        LOGGER.info("Found k=6 model")
    else:
        LOGGER.warning("k=6 model not found (checkpoint=%s, vocab=%s)", k6_checkpoint.exists(), k6_vocab.exists())

    # Fallback to k=3 if neither k=4 nor k=6 available
    if not kmer_config:
        k3_checkpoint = artifacts_dir / "checkpoints" / "stage1_k3.pt"
        k3_vocab = artifacts_dir / "vocabs" / "k3_vocab.json"
        if k3_checkpoint.exists() and k3_vocab.exists():
            kmer_config[3] = {
                "checkpoint": str(k3_checkpoint),
                "vocab": str(k3_vocab),
            }
            LOGGER.info("Falling back to k=3 model")
        else:
            LOGGER.error("No models available for ensemble!")
            return 1

    if not config_path.exists():
        LOGGER.error("Config not found: %s", config_path)
        return 1

    try:
        LOGGER.info("Loading ensemble with k-mers: %s", list(kmer_config.keys()))
        predictor = load_stage1_ensemble(
            config_path=config_path,
            checkpoints_and_vocabs=kmer_config,
        )

        LOGGER.info("Running prediction...")
        result = predictor.predict(test_sequence)

        # Output as JSON (matches ensemble_infer.py output semantics)
        print("\n" + "=" * 60)
        print("SMOKE TEST RESULT:")
        print("=" * 60)
        print(json.dumps(result, indent=2))
        print("=" * 60 + "\n")

        # Validate result structure
        assert "probability" in result, "Missing 'probability' key"
        assert "label" in result, "Missing 'label' key"
        assert isinstance(result["probability"], float), "probability must be float"
        assert result["label"] in ("promoter", "non-promoter"), "Invalid label"
        assert 0.0 <= result["probability"] <= 1.0, "probability out of range"

        LOGGER.info("Smoke test PASSED!")
        LOGGER.info("Probability: %.4f", result["probability"])
        LOGGER.info("Label: %s", result["label"])
        return 0

    except Exception as e:
        LOGGER.exception("Smoke test FAILED: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(run_smoke_test())
