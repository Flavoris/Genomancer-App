#!/usr/bin/env python3
"""Quick experiment to test regularization impact on overfitting.

This script tests 3 different regularization levels and compares
their effects on the train/val accuracy gap. The goal is to find
settings where the gap is minimized while validation accuracy remains high.

Usage:
    python run_regularization_experiment.py
"""

import json
import sys
from pathlib import Path
from typing import Any

# Add parent to path for imports
script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))

from train_stage1 import run_stage1_training, load_config

# Configurations to test
EXPERIMENTS = [
    {
        "name": "baseline",
        "transformer_dropout": 0.12,
        "classifier_dropout": 0.15,
        "weight_decay": 0.02,
        "label_smoothing": 0.0,
    },
    {
        "name": "medium_reg",
        "transformer_dropout": 0.20,
        "classifier_dropout": 0.20,
        "weight_decay": 0.05,
        "label_smoothing": 0.1,
    },
    {
        "name": "strong_reg",
        "transformer_dropout": 0.25,
        "classifier_dropout": 0.30,
        "weight_decay": 0.08,
        "label_smoothing": 0.15,
    },
]


def run_experiment(
    base_config: dict, exp_config: dict, epochs: int = 30
) -> dict[str, Any] | None:
    """Run training with specified config overrides.

    Args:
        base_config: Base configuration loaded from config.yaml
        exp_config: Experiment-specific config overrides
        epochs: Number of training epochs

    Returns:
        Dictionary with experiment results, or None if failed
    """
    name = exp_config["name"]
    output_dir = script_dir / "results" / "reg_experiment" / name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build overrides
    overrides = {
        "epochs": epochs,
        "transformer_dropout": exp_config["transformer_dropout"],
        "label_smoothing": exp_config["label_smoothing"],
        "weight_decay": exp_config["weight_decay"],
        # Classifier dropout goes into simplified_model section
        "_classifier_dropout_override": exp_config["classifier_dropout"],
    }

    print(f"\n{'='*60}")
    print(f"Running experiment: {name}")
    print(f"{'='*60}")
    print(f"  transformer_dropout: {exp_config['transformer_dropout']}")
    print(f"  classifier_dropout:  {exp_config['classifier_dropout']}")
    print(f"  weight_decay:        {exp_config['weight_decay']}")
    print(f"  label_smoothing:     {exp_config['label_smoothing']}")
    print(f"  epochs:              {epochs}")
    print(f"  output_dir:          {output_dir}")
    print()

    try:
        # Update simplified_model classifier_dropout in base config
        if "simplified_model" not in base_config:
            base_config["simplified_model"] = {}
        base_config["simplified_model"]["classifier_dropout"] = exp_config["classifier_dropout"]

        # Run training
        result = run_stage1_training(base_config, overrides=overrides)

        # Save results
        metrics_file = output_dir / "final_metrics.json"
        with open(metrics_file, "w") as f:
            json.dump(result, f, indent=2)

        return result

    except Exception as e:
        print(f"ERROR: Experiment '{name}' failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Run all regularization experiments and compare results."""
    config_path = script_dir / "config.yaml"
    base_config = load_config(config_path)

    results = []

    for exp_config in EXPERIMENTS:
        # Reload config for each experiment to avoid contamination
        fresh_config = load_config(config_path)
        metrics = run_experiment(fresh_config, exp_config, epochs=30)

        if metrics:
            # Extract relevant metrics
            eval_metrics = metrics.get("evaluation_metrics", {})
            training_result = metrics.get("training_result", metrics)

            results.append({
                "name": exp_config["name"],
                "val_acc": eval_metrics.get("accuracy", training_result.get("best_val_accuracy", 0)),
                "val_mcc": eval_metrics.get("mcc", training_result.get("best_val_mcc", 0)),
                "best_epoch": training_result.get("best_epoch", 0),
                **exp_config
            })

    # Print comparison
    print("\n" + "="*80)
    print("REGULARIZATION EXPERIMENT RESULTS")
    print("="*80)
    print(f"{'Name':<15} {'Val Acc':<10} {'Val MCC':<10} {'Best Ep':<10} {'Dropout':<10}")
    print("-"*80)

    for r in results:
        print(
            f"{r['name']:<15} "
            f"{r['val_acc']:.4f}     "
            f"{r['val_mcc']:.4f}     "
            f"{r['best_epoch']:<10} "
            f"{r['transformer_dropout']}"
        )

    if results:
        # Find best by validation accuracy
        best = max(results, key=lambda x: x['val_acc'])
        print("\n" + "="*80)
        print(f"BEST CONFIG: {best['name']} with {best['val_acc']:.4f} accuracy")
        print("="*80)

        # Save summary
        summary_path = script_dir / "results" / "reg_experiment" / "summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump({
                "experiments": results,
                "best_config": best['name'],
                "best_accuracy": best['val_acc'],
            }, f, indent=2)
        print(f"\nSummary saved to: {summary_path}")
    else:
        print("\nNo experiments completed successfully.")


if __name__ == "__main__":
    main()
