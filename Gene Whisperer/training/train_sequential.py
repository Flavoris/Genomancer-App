"""Sequential Fine-Tuning for Gene Whisperer.

Implements msBERT's two-stage approach:
1. Fine-tune on promoter identification (Stage 1) for each k-mer
2. Initialize Stage 2 from Stage 1 and fine-tune on strength
3. Evaluate soft voting ensemble

This orchestrator trains models for each k-mer size independently,
then combines them using soft voting for maximum accuracy.

Usage:
    python train_sequential.py --config config.yaml
    python train_sequential.py --config config.yaml --kmers 4 6
"""
from __future__ import annotations

import argparse
import copy
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import torch
import yaml

LOGGER = logging.getLogger("gene_whisperer.sequential")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def load_config(config_path: Path) -> Dict:
    """Load YAML configuration file."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg or {}


def sequential_train(
    config_path: str,
    kmers: Optional[List[int]] = None,
    skip_stage2: bool = False,
) -> Dict:
    """
    Execute sequential fine-tuning.

    Args:
        config_path: Path to config.yaml
        kmers: List of k-mer sizes to train. Defaults to config's multi_scale_kmers.
        skip_stage2: Skip Stage 2 training (useful if no Stage 2 data available)

    Returns:
        Dict containing ensemble results and individual model checkpoints
    """
    from train_stage1 import run_stage1_training, load_config as load_stage1_config

    config_path = Path(config_path).resolve()
    cfg = load_stage1_config(config_path)

    # Get k-mers to train
    if kmers is None:
        kmers = cfg.get("multi_scale_kmers", [3, 4, 5, 6])

    print("=" * 70)
    print("SEQUENTIAL FINE-TUNING")
    print(f"K-mers: {kmers}")
    print("=" * 70)

    # Phase 1: Train Stage 1 models for each k-mer
    stage1_checkpoints: Dict[int, str] = {}
    stage1_results: Dict[int, Dict] = {}

    for kmer in kmers:
        print(f"\n[Phase 1] Training Stage 1 for k={kmer}")
        print("-" * 60)

        # Create overrides for this k-mer
        # Disable multi-scale since sequential training handles k-mers independently
        overrides = {
            "_config_path": str(config_path),
            "kmer": kmer,
            "stage1_checkpoint_name": f"stage1_k{kmer}_sequential.pt",
            "use_multi_scale": False,
        }

        # Check for k-mer specific MLM vocab
        vocab_cache_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
        mlm_vocab_path = vocab_cache_dir / f"mlm_k{kmer}_vocab.json"
        if mlm_vocab_path.exists():
            overrides["mlm_vocab_path"] = str(mlm_vocab_path)

        # Check for k-mer specific MLM checkpoint
        mlm_encoder_path = Path(cfg.get("mlm_encoder_path", "")).parent / f"mlm_encoder_k{kmer}.pt"
        if mlm_encoder_path.exists():
            overrides["mlm_encoder_checkpoint"] = str(mlm_encoder_path)

        # Run training
        results = run_stage1_training(cfg, overrides=overrides)

        stage1_checkpoints[kmer] = results["best_ckpt_path"]
        stage1_results[kmer] = results

        print(f"\nStage 1 k={kmer} Complete:")
        print(f"  Accuracy: {results['best_val_acc']:.4f} ({results['best_val_acc']*100:.2f}%)")
        print(f"  MCC: {results['best_val_mcc']:.4f}")
        print(f"  AUC: {results['best_val_auc']:.4f}")
        print(f"  Checkpoint: {results['best_ckpt_path']}")

    # Phase 2: Sequential fine-tuning for Stage 2 (if Stage 2 data exists)
    stage2_results: Dict[int, Dict] = {}

    if not skip_stage2 and cfg.get("stage2_train"):
        print(f"\n[Phase 2] Sequential Stage 2 Training")
        print("-" * 60)

        for kmer in kmers:
            print(f"\nInitializing Stage 2 k={kmer} from Stage 1 checkpoint")

            stage1_ckpt = stage1_checkpoints[kmer]

            try:
                from train_stage2 import run_stage2_training

                # Create Stage 2 overrides
                stage2_overrides = {
                    "_config_path": str(config_path),
                    "stage2_kmer": kmer,
                    "stage1_checkpoint": stage1_ckpt,
                    "stage2_checkpoint_name": f"stage2_k{kmer}_sequential.pt",
                }

                results = run_stage2_training(cfg, overrides=stage2_overrides)
                stage2_results[kmer] = results

                print(f"Stage 2 k={kmer}: Acc={results.get('best_val_acc', 0):.4f}")

            except ImportError:
                LOGGER.warning("train_stage2.run_stage2_training not available, skipping Stage 2")
                break
            except Exception as e:
                LOGGER.warning(f"Stage 2 training failed for k={kmer}: {e}")
                continue
    else:
        if skip_stage2:
            print("\n[Phase 2] Skipped (--skip_stage2 flag)")
        else:
            print("\n[Phase 2] Skipped (no stage2_train data configured)")

    # Phase 3: Ensemble evaluation
    print(f"\n[Phase 3] Ensemble Evaluation")
    print("-" * 60)

    from evaluate_ensemble import evaluate_soft_voting_ensemble

    ensemble_results = evaluate_soft_voting_ensemble(
        config_path=str(config_path),
        checkpoint_paths=stage1_checkpoints,
    )

    # Print final summary
    print("\n" + "=" * 70)
    print("SEQUENTIAL TRAINING COMPLETE")
    print("=" * 70)

    print("\nIndividual Model Results:")
    for kmer in kmers:
        r = stage1_results[kmer]
        print(f"  k={kmer}: Acc={r['best_val_acc']*100:.2f}%, MCC={r['best_val_mcc']:.4f}, AUC={r['best_val_auc']:.4f}")

    print(f"\nEnsemble Results (Soft Voting):")
    print(f"  Accuracy: {ensemble_results['accuracy']:.4f} ({ensemble_results['accuracy']*100:.2f}%)")
    print(f"  MCC: {ensemble_results['mcc']:.4f}")
    print(f"  AUC: {ensemble_results['auc']:.4f}")

    # Calculate improvement from ensemble
    best_individual_acc = max(r['best_val_acc'] for r in stage1_results.values())
    ensemble_acc = ensemble_results['accuracy']
    improvement = ensemble_acc - best_individual_acc

    print(f"\nEnsemble Improvement over Best Individual:")
    print(f"  +{improvement*100:.2f}% accuracy gain")

    print("=" * 70)

    return {
        "stage1_checkpoints": stage1_checkpoints,
        "stage1_results": stage1_results,
        "stage2_results": stage2_results,
        "ensemble_results": ensemble_results,
        "config_path": str(config_path),
    }


def main() -> None:
    """CLI entry point for sequential training."""
    parser = argparse.ArgumentParser(
        description="Sequential Fine-Tuning for Gene Whisperer"
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file",
    )
    parser.add_argument(
        "--kmers",
        type=int,
        nargs="+",
        default=None,
        help="K-mer sizes to train (e.g., --kmers 3 4 5 6)",
    )
    parser.add_argument(
        "--skip_stage2",
        action="store_true",
        help="Skip Stage 2 training (only train Stage 1 models)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save results JSON",
    )
    args = parser.parse_args()

    results = sequential_train(
        config_path=args.config,
        kmers=args.kmers,
        skip_stage2=args.skip_stage2,
    )

    # Save results if output path specified
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert paths to strings for JSON serialization
        serializable = {
            "stage1_checkpoints": results["stage1_checkpoints"],
            "stage1_results": results["stage1_results"],
            "stage2_results": results["stage2_results"],
            "ensemble_results": results["ensemble_results"],
            "config_path": results["config_path"],
            "timestamp": datetime.now().isoformat(),
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(serializable, f, indent=2)

        print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
