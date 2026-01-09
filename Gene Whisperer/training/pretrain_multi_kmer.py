"""
Multi-K-mer MLM Pretraining

Trains separate DNABERT-style models for each k-mer size (3, 4, 5, 6).
Each model learns representations specific to its vocabulary size.

Usage:
    python pretrain_multi_kmer.py --config config.yaml
    python pretrain_multi_kmer.py --config config.yaml --kmers 3 4  # Only k=3,4
    python pretrain_multi_kmer.py --config config.yaml --resume  # Resume from checkpoints

This script leverages the existing run_mlm_pretrain() function with per-kmer overrides.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from pretrain_mlm import run_mlm_pretrain
from early_stopping import PretrainingEarlyStopping

LOGGER = logging.getLogger("pretrain_multi_kmer")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Environment setup
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_file) as f:
        return yaml.safe_load(f)


def get_kmer_config(cfg: dict, kmer: int) -> dict:
    """
    Get configuration overrides for specific k-mer.

    Args:
        cfg: Base configuration dictionary
        kmer: K-mer size (3, 4, 5, or 6)

    Returns:
        Dictionary of config overrides for this k-mer
    """
    overrides = {}

    # Set k-mer specific parameters
    overrides["mlm_kmer"] = kmer

    # Calculate vocab size: 4^k + 3 special tokens ([MASK], <UNK>, <PAD>)
    vocab_size = (4 ** kmer) + 3
    overrides["vocab_size"] = vocab_size

    # Set vocab path
    vocab_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
    overrides["mlm_vocab_path"] = str(vocab_dir / f"mlm_k{kmer}_vocab.json")

    # Set encoder checkpoint paths
    checkpoints = cfg.get("checkpoints", {})
    encoder_by_k = checkpoints.get("mlm_encoder_ckpt_by_k", {})
    metadata_by_k = checkpoints.get("mlm_metadata_by_k", {})

    if kmer in encoder_by_k:
        overrides["mlm_encoder_path"] = encoder_by_k[kmer]
    else:
        overrides["mlm_encoder_path"] = f"../artifacts/mlm_encoder_k{kmer}.pt"

    # Apply k-mer specific overrides from config
    kmer_overrides = cfg.get("kmer_pretrain_overrides", {}).get(kmer, {})
    for key, value in kmer_overrides.items():
        overrides[key] = value

    return overrides


def check_existing_checkpoint(cfg: dict, kmer: int) -> Optional[Dict[str, Any]]:
    """
    Check if an existing checkpoint exists for this k-mer.

    Args:
        cfg: Configuration dictionary
        kmer: K-mer size

    Returns:
        Checkpoint metadata if exists, None otherwise
    """
    script_dir = Path(__file__).parent

    # Get encoder path from config
    checkpoints = cfg.get("checkpoints", {})
    encoder_by_k = checkpoints.get("mlm_encoder_ckpt_by_k", {})
    metadata_by_k = checkpoints.get("mlm_metadata_by_k", {})

    if kmer in encoder_by_k:
        encoder_path = Path(encoder_by_k[kmer])
    else:
        encoder_path = Path(f"../artifacts/mlm_encoder_k{kmer}.pt")

    if not encoder_path.is_absolute():
        encoder_path = (script_dir / encoder_path).resolve()

    metadata_path = encoder_path.with_suffix(".metadata.json")

    if encoder_path.exists() and metadata_path.exists():
        try:
            with open(metadata_path) as f:
                metadata = json.load(f)
            LOGGER.info(f"Found existing checkpoint for k={kmer}: {encoder_path}")
            return metadata
        except Exception as e:
            LOGGER.warning(f"Could not load metadata for k={kmer}: {e}")

    return None


def pretrain_single_kmer(
    cfg: dict,
    kmer: int,
    resume: bool = False,
    skip_if_exists: bool = True,
) -> Dict[str, Any]:
    """
    Pretrain MLM model for a single k-mer size.

    Args:
        cfg: Base configuration
        kmer: K-mer size (3, 4, 5, or 6)
        resume: Whether to resume from existing checkpoint
        skip_if_exists: If True, skip training if checkpoint already exists

    Returns:
        Dict with training results
    """
    LOGGER.info("=" * 70)
    LOGGER.info(f"PRETRAINING K={kmer} MLM MODEL")
    LOGGER.info("=" * 70)

    # Check for existing checkpoint
    existing_checkpoint = check_existing_checkpoint(cfg, kmer)

    if existing_checkpoint and skip_if_exists and not resume:
        LOGGER.info(f"Checkpoint exists for k={kmer}, skipping (use --resume to retrain)")
        return {
            "kmer": kmer,
            "status": "skipped",
            "reason": "checkpoint_exists",
            "best_val_loss": existing_checkpoint.get("best_val_loss"),
            "checkpoint_path": existing_checkpoint.get("encoder_path"),
        }

    # Get k-mer specific overrides
    kmer_overrides = get_kmer_config(cfg, kmer)

    # Log configuration
    LOGGER.info(f"K-mer: {kmer}")
    LOGGER.info(f"Vocab size: {kmer_overrides.get('vocab_size', 'default')}")
    LOGGER.info(f"Epochs: {kmer_overrides.get('mlm_epochs', cfg.get('mlm_epochs', 200))}")
    LOGGER.info(f"Learning rate: {kmer_overrides.get('mlm_lr', cfg.get('mlm_lr', 0.0002))}")
    LOGGER.info(f"Vocab path: {kmer_overrides.get('mlm_vocab_path')}")
    LOGGER.info(f"Encoder path: {kmer_overrides.get('mlm_encoder_path')}")

    # Verify vocabulary exists
    vocab_path = Path(kmer_overrides["mlm_vocab_path"])
    if not vocab_path.is_absolute():
        vocab_path = (Path(__file__).parent / vocab_path).resolve()

    if not vocab_path.exists():
        raise FileNotFoundError(
            f"Vocabulary not found: {vocab_path}\n"
            f"Run build_multi_kmer_vocabs.py first to generate vocabularies."
        )

    # Run pretraining
    start_time = time.time()
    try:
        result = run_mlm_pretrain(cfg, overrides=kmer_overrides)
        training_time = time.time() - start_time

        return {
            "kmer": kmer,
            "status": "completed",
            "mlm_loss": result.get("mlm_loss"),
            "mlm_steps": result.get("mlm_steps"),
            "encoder_ckpt_path": result.get("encoder_ckpt_path"),
            "pretrain_seconds": result.get("pretrain_seconds", training_time),
            "params": result.get("params"),
        }
    except Exception as e:
        LOGGER.error(f"Training failed for k={kmer}: {e}")
        return {
            "kmer": kmer,
            "status": "failed",
            "error": str(e),
        }


def pretrain_all_kmers(
    config_path: str,
    kmers: Optional[List[int]] = None,
    resume: bool = False,
    skip_existing: bool = True,
) -> Dict[int, Dict[str, Any]]:
    """
    Pretrain MLM models for all k-mer sizes sequentially.

    Args:
        config_path: Path to config.yaml
        kmers: List of k-mer sizes to pretrain (default: [3, 4, 5, 6])
        resume: Whether to resume from existing checkpoints
        skip_existing: If True, skip k-mers that already have checkpoints

    Returns:
        Dict mapping k-mer to training results
    """
    cfg = load_config(config_path)

    # Default k-mers from config or use default [3, 4, 5, 6]
    if kmers is None:
        kmers = cfg.get("pretrain_kmers", [3, 4, 5, 6])

    # Get multi-kmer config section
    multi_kmer_cfg = cfg.get("multi_kmer_pretrain", {})
    if not multi_kmer_cfg.get("multi_kmer_pretrain_enabled", True):
        LOGGER.warning("Multi-kmer pretraining is disabled in config")

    LOGGER.info("=" * 70)
    LOGGER.info("MULTI-KMER MLM PRETRAINING")
    LOGGER.info("=" * 70)
    LOGGER.info(f"K-mers to pretrain: {kmers}")
    LOGGER.info(f"Resume from checkpoints: {resume}")
    LOGGER.info(f"Skip existing checkpoints: {skip_existing}")
    LOGGER.info("=" * 70)

    results: Dict[int, Dict[str, Any]] = {}
    total_start_time = time.time()

    for i, kmer in enumerate(kmers):
        LOGGER.info(f"\n[{i+1}/{len(kmers)}] Starting k={kmer} pretraining...\n")

        try:
            result = pretrain_single_kmer(
                cfg,
                kmer,
                resume=resume,
                skip_if_exists=skip_existing and not resume,
            )
            results[kmer] = result

            # Log result
            status = result.get("status", "unknown")
            if status == "completed":
                LOGGER.info(
                    f"\n[OK] K={kmer} pretraining complete: "
                    f"loss={result.get('mlm_loss', 'N/A'):.4f}, "
                    f"steps={result.get('mlm_steps', 'N/A')}"
                )
            elif status == "skipped":
                LOGGER.info(f"\n[SKIPPED] K={kmer}: {result.get('reason', 'unknown')}")
            else:
                LOGGER.error(f"\n[FAILED] K={kmer}: {result.get('error', 'unknown error')}")

        except Exception as e:
            LOGGER.error(f"\n[FAILED] K={kmer} pretraining failed: {e}")
            results[kmer] = {"kmer": kmer, "status": "failed", "error": str(e)}

    total_time = time.time() - total_start_time

    # Print summary
    _print_summary(results, total_time)

    # Save summary to file
    _save_summary(results, total_time, config_path)

    return results


def _print_summary(results: Dict[int, Dict[str, Any]], total_time: float) -> None:
    """Print training summary."""
    LOGGER.info("\n" + "=" * 70)
    LOGGER.info("PRETRAINING SUMMARY")
    LOGGER.info("=" * 70)

    completed = sum(1 for r in results.values() if r.get("status") == "completed")
    skipped = sum(1 for r in results.values() if r.get("status") == "skipped")
    failed = sum(1 for r in results.values() if r.get("status") == "failed")

    LOGGER.info(f"Total time: {total_time/60:.1f} minutes")
    LOGGER.info(f"Completed: {completed}, Skipped: {skipped}, Failed: {failed}")
    LOGGER.info("-" * 70)

    for kmer in sorted(results.keys()):
        result = results[kmer]
        status = result.get("status", "unknown")

        if status == "completed":
            LOGGER.info(
                f"k={kmer}: loss={result.get('mlm_loss', 'N/A'):.4f}, "
                f"steps={result.get('mlm_steps', 'N/A')}, "
                f"time={result.get('pretrain_seconds', 0)/60:.1f}min"
            )
        elif status == "skipped":
            LOGGER.info(f"k={kmer}: SKIPPED ({result.get('reason', 'N/A')})")
        else:
            LOGGER.info(f"k={kmer}: FAILED ({result.get('error', 'unknown')[:50]})")

    LOGGER.info("=" * 70)


def _save_summary(
    results: Dict[int, Dict[str, Any]],
    total_time: float,
    config_path: str,
) -> None:
    """Save training summary to JSON file."""
    script_dir = Path(__file__).parent
    summary_path = (script_dir / "../artifacts/multi_kmer_pretrain_summary.json").resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "timestamp": datetime.now().isoformat(),
        "config_path": config_path,
        "total_time_seconds": total_time,
        "results": {str(k): v for k, v in results.items()},  # Convert keys to strings
    }

    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    LOGGER.info(f"Saved summary to {summary_path}")


def validate_vocabs_exist(cfg: dict, kmers: List[int]) -> List[int]:
    """
    Check which vocabularies exist.

    Args:
        cfg: Configuration dictionary
        kmers: List of k-mer sizes to check

    Returns:
        List of k-mers with missing vocabularies
    """
    script_dir = Path(__file__).parent
    vocab_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
    if not vocab_dir.is_absolute():
        vocab_dir = (script_dir / vocab_dir).resolve()

    missing = []
    for kmer in kmers:
        vocab_path = vocab_dir / f"mlm_k{kmer}_vocab.json"
        if not vocab_path.exists():
            missing.append(kmer)

    return missing


def main():
    parser = argparse.ArgumentParser(
        description="Multi-K-mer MLM Pretraining",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Train all k-mers (3, 4, 5, 6)
    python pretrain_multi_kmer.py --config config.yaml

    # Train only k=3 and k=4
    python pretrain_multi_kmer.py --config config.yaml --kmers 3 4

    # Resume training (don't skip existing)
    python pretrain_multi_kmer.py --config config.yaml --resume

    # Force retrain all
    python pretrain_multi_kmer.py --config config.yaml --force
        """,
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to config file (default: config.yaml)",
    )
    parser.add_argument(
        "--kmers",
        type=int,
        nargs="+",
        default=None,
        help="K-mer sizes to pretrain (default: from config or [3, 4, 5, 6])",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume training from existing checkpoints",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force retrain even if checkpoints exist",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate that vocabularies exist, don't train",
    )

    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    # Determine k-mers to train
    kmers = args.kmers or cfg.get("pretrain_kmers", [3, 4, 5, 6])

    # Validate vocabularies exist
    missing_vocabs = validate_vocabs_exist(cfg, kmers)
    if missing_vocabs:
        LOGGER.error(
            f"Missing vocabularies for k={missing_vocabs}. "
            f"Run build_multi_kmer_vocabs.py first."
        )
        if not args.validate_only:
            sys.exit(1)

    if args.validate_only:
        if missing_vocabs:
            LOGGER.info(f"Validation failed: missing vocabs for k={missing_vocabs}")
            sys.exit(1)
        else:
            LOGGER.info(f"Validation passed: all vocabularies exist for k={kmers}")
            sys.exit(0)

    # Run pretraining
    pretrain_all_kmers(
        config_path=args.config,
        kmers=kmers,
        resume=args.resume,
        skip_existing=not args.force and not args.resume,
    )


if __name__ == "__main__":
    main()
