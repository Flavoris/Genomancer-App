#!/usr/bin/env python3
"""
Quick Variant Sanity Check Script

Fast iteration tool for verifying model variants work correctly before
running full experiments. Runs 5 epochs with 500 samples each for rapid
feedback during development.

Target runtime: <5 minutes on GPU, <15 minutes on CPU
"""

import sys
import time
from contextlib import nullcontext
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from amp_utils import create_grad_scaler, get_amp_context
from dataset import build_dataloaders
from model import GeneWhispererStage1Legacy, create_model_variant
from seed_utils import set_global_seed
from train_stage1 import calculate_metrics, load_config

# Sanity check configuration
SANITY_CONFIG = {
    "epochs": 5,
    "max_samples": 500,
    "batch_size": 32,
    "learning_rate": 1e-4,
    "seed": 42,
}

VARIANTS_TO_TEST = ["transformer_only", "features_only", "combined"]


def create_sanity_dataloader(full_loader: DataLoader, max_samples: int, batch_size: int) -> DataLoader:
    """Create a smaller dataloader for quick sanity checks."""
    dataset = full_loader.dataset
    num_samples = min(max_samples, len(dataset))
    indices = list(range(num_samples))
    subset = Subset(dataset, indices)

    return DataLoader(
        subset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=0,  # Faster for small datasets
        pin_memory=torch.cuda.is_available(),
    )


def check_for_nan_inf(tensor: torch.Tensor, name: str) -> list[str]:
    """Check tensor for NaN or Inf values."""
    issues = []
    if torch.isnan(tensor).any():
        issues.append(f"NaN detected in {name}")
    if torch.isinf(tensor).any():
        issues.append(f"Inf detected in {name}")
    return issues


def check_gradients(model: nn.Module) -> tuple[bool, list[str]]:
    """Verify gradients flow correctly through the model."""
    issues = []
    has_grad = False

    for name, param in model.named_parameters():
        if param.requires_grad:
            if param.grad is not None:
                has_grad = True
                issues.extend(check_for_nan_inf(param.grad, f"gradient of {name}"))

    if not has_grad:
        issues.append("No gradients computed - check backward pass")

    return len(issues) == 0, issues


def run_sanity_epoch(
    loader: DataLoader,
    model: nn.Module,
    criterion: nn.Module,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    variant: str,
    is_training: bool = True,
    amp_ctx=None,
    scaler=None,
) -> dict:
    """Run a single epoch for sanity checking."""
    model.train() if is_training else model.eval()

    all_probs = []
    all_labels = []
    total_loss = 0.0
    num_batches = 0
    nan_inf_issues = []

    context = torch.enable_grad() if is_training else torch.no_grad()

    with context:
        for tokens, engineered, labels in loader:
            tokens = tokens.to(device)
            engineered = engineered.to(device)
            labels = labels.to(device).float()

            # Check inputs for NaN/Inf
            nan_inf_issues.extend(check_for_nan_inf(tokens.float(), "input tokens"))
            nan_inf_issues.extend(check_for_nan_inf(engineered, "engineered features"))

            if is_training:
                optimizer.zero_grad()

            # Forward pass with variant-specific inputs
            with amp_ctx if amp_ctx else nullcontext():
                if variant == "transformer_only":
                    logits = model(tokens, engineered_features=None)
                elif variant == "features_only":
                    logits = model(tokens=None, engineered_features=engineered)
                else:  # combined
                    logits = model(tokens, engineered_features=engineered)

                logits = logits.squeeze(-1)
                loss = criterion(logits, labels)

            # Check outputs for NaN/Inf
            nan_inf_issues.extend(check_for_nan_inf(logits, "model logits"))
            nan_inf_issues.extend(check_for_nan_inf(loss, "loss"))

            if is_training:
                if scaler is not None:
                    scaler.scale(loss).backward()
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    loss.backward()
                    optimizer.step()

            probs = torch.sigmoid(logits.detach()).cpu()
            all_probs.append(probs)
            all_labels.append(labels.detach().cpu())
            total_loss += loss.item()
            num_batches += 1

    # Concatenate all predictions
    all_probs = torch.cat(all_probs)
    all_labels = torch.cat(all_labels)

    # Calculate metrics (expects tensors, not numpy)
    metrics = calculate_metrics(all_probs, all_labels)
    metrics["loss"] = total_loss / max(num_batches, 1)
    metrics["nan_inf_issues"] = nan_inf_issues

    return metrics


def train_variant_sanity(
    variant: str,
    cfg: dict,
    train_loader: DataLoader,
    val_loader: DataLoader,
    device: torch.device,
    epochs: int,
) -> dict:
    """Train a single variant for sanity checking."""
    print(f"\nVariant: {variant}")

    # Create model
    if variant == "original":
        model = GeneWhispererStage1Legacy(
            vocab_size=cfg["vocab_size"],
            embedding_dim=cfg["embedding_dim"],
            max_seq_len=cfg["max_bp_len"] - cfg["kmer"] + 2,
            num_layers=cfg.get("transformer_layers", 6),
            num_heads=cfg.get("transformer_heads", 6),
            ff_dim=cfg.get("transformer_ff_dim", 1024),
            dropout=cfg.get("transformer_dropout", 0.1),
            pad_token_id=cfg.get("pad_token_id", 4098),
            engineered_dim=cfg.get("engineered_dim", 288),
        ).to(device)
    else:
        model = create_model_variant(variant, cfg).to(device)

    # Setup training
    criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=SANITY_CONFIG["learning_rate"],
        weight_decay=0.01,
    )

    # AMP setup
    use_amp = device.type == "cuda"
    amp_ctx = get_amp_context(enabled=use_amp, device=device, dtype="float16") if use_amp else None
    scaler = create_grad_scaler(device) if use_amp else None

    results = {
        "variant": variant,
        "epochs": [],
        "gradient_ok": True,
        "nan_inf_issues": [],
        "final_val_acc": 0.0,
        "final_val_mcc": 0.0,
    }

    for epoch in range(1, epochs + 1):
        # Training
        train_metrics = run_sanity_epoch(
            train_loader, model, criterion, optimizer, device, variant,
            is_training=True, amp_ctx=amp_ctx, scaler=scaler,
        )

        # Check gradients after first training batch
        if epoch == 1:
            grad_ok, grad_issues = check_gradients(model)
            results["gradient_ok"] = grad_ok
            if not grad_ok:
                results["nan_inf_issues"].extend(grad_issues)

        # Validation
        val_metrics = run_sanity_epoch(
            val_loader, model, criterion, optimizer, device, variant,
            is_training=False, amp_ctx=amp_ctx, scaler=scaler,
        )

        # Collect NaN/Inf issues
        results["nan_inf_issues"].extend(train_metrics.get("nan_inf_issues", []))
        results["nan_inf_issues"].extend(val_metrics.get("nan_inf_issues", []))

        # Print epoch results
        print(f"  Epoch {epoch}: train_acc={train_metrics['accuracy']:.2f}, "
              f"val_acc={val_metrics['accuracy']:.2f}")

        results["epochs"].append({
            "epoch": epoch,
            "train_acc": train_metrics["accuracy"],
            "train_mcc": train_metrics["mcc"],
            "val_acc": val_metrics["accuracy"],
            "val_mcc": val_metrics["mcc"],
        })

    # Store final metrics
    results["final_val_acc"] = val_metrics["accuracy"]
    results["final_val_mcc"] = val_metrics["mcc"]

    print(f"  Final: val_acc={results['final_val_acc']:.2f}, "
          f"val_mcc={results['final_val_mcc']:.2f}")

    return results


def run_sanity_check(config_path: Path | None = None) -> bool:
    """
    Run the full sanity check suite.

    Returns:
        True if all checks pass, False otherwise.
    """
    print("=" * 50)
    print("=== Quick Sanity Check ===")
    print(f"Running {SANITY_CONFIG['epochs']} epochs with "
          f"{SANITY_CONFIG['max_samples']} samples each...")
    print("=" * 50)

    start_time = time.time()

    # Set seed for reproducibility
    set_global_seed(SANITY_CONFIG["seed"])

    # Device setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load config
    if config_path is None:
        config_path = Path(__file__).parent / "config.yaml"

    if not config_path.exists():
        print(f"ERROR: Config file not found at {config_path}")
        return False

    cfg = load_config(config_path)

    # Build full dataloaders
    try:
        loaders = build_dataloaders(cfg)
        full_train_loader = loaders["stage1"]["train"]
        full_val_loader = loaders["stage1"]["val"]
    except Exception as e:
        print(f"ERROR: Failed to build dataloaders: {e}")
        return False

    # Create smaller dataloaders for sanity check
    train_loader = create_sanity_dataloader(
        full_train_loader,
        SANITY_CONFIG["max_samples"],
        SANITY_CONFIG["batch_size"],
    )
    val_loader = create_sanity_dataloader(
        full_val_loader,
        SANITY_CONFIG["max_samples"] // 5,  # Smaller val set
        SANITY_CONFIG["batch_size"],
    )

    print(f"Train samples: {len(train_loader.dataset)}")
    print(f"Val samples: {len(val_loader.dataset)}")

    # Train each variant
    all_results = {}
    errors = []

    for variant in VARIANTS_TO_TEST:
        try:
            results = train_variant_sanity(
                variant=variant,
                cfg=cfg,
                train_loader=train_loader,
                val_loader=val_loader,
                device=device,
                epochs=SANITY_CONFIG["epochs"],
            )
            all_results[variant] = results
        except Exception as e:
            error_msg = f"Variant {variant} failed: {e}"
            print(f"  ERROR: {error_msg}")
            errors.append(error_msg)
            all_results[variant] = {"error": str(e)}

    # Print summary
    elapsed = time.time() - start_time
    print("\n" + "=" * 50)
    print("=== Sanity Check Results ===")
    print("=" * 50)

    checks_passed = []
    checks_failed = []

    # Check 1: All variants train without error
    variants_ok = all(
        "error" not in all_results.get(v, {"error": True})
        for v in VARIANTS_TO_TEST
    )
    if variants_ok:
        checks_passed.append("All variants train without error")
    else:
        checks_failed.append("Some variants failed to train")

    # Check 2: Combined >= max(transformer_only, features_only)
    if variants_ok:
        transformer_acc = all_results["transformer_only"]["final_val_acc"]
        features_acc = all_results["features_only"]["final_val_acc"]
        combined_acc = all_results["combined"]["final_val_acc"]
        baseline_max = max(transformer_acc, features_acc)

        # Allow small tolerance for randomness
        if combined_acc >= baseline_max - 0.05:
            checks_passed.append(f"Combined ({combined_acc:.2f}) >= max(baselines) ({baseline_max:.2f})")
        else:
            checks_failed.append(
                f"Combined ({combined_acc:.2f}) < max(baselines) ({baseline_max:.2f}) "
                f"- fusion may not be working"
            )

    # Check 3: No NaN/Inf in any outputs
    all_nan_inf_issues = []
    for v, results in all_results.items():
        if isinstance(results, dict) and "nan_inf_issues" in results:
            all_nan_inf_issues.extend(results["nan_inf_issues"])

    if not all_nan_inf_issues:
        checks_passed.append("No NaN/Inf in any outputs")
    else:
        checks_failed.append(f"NaN/Inf detected: {all_nan_inf_issues[:3]}")

    # Check 4: Gradients flow correctly
    gradients_ok = all(
        all_results.get(v, {}).get("gradient_ok", False)
        for v in VARIANTS_TO_TEST
        if "error" not in all_results.get(v, {"error": True})
    )
    if gradients_ok:
        checks_passed.append("Gradients flow correctly")
    else:
        checks_failed.append("Gradient issues detected")

    # Print results
    for check in checks_passed:
        print(f"✓ {check}")

    for check in checks_failed:
        print(f"✗ {check}")

    print(f"\nTotal time: {elapsed:.1f}s")

    all_passed = len(checks_failed) == 0
    if all_passed:
        print("\n🎉 Ready for full experiment!")
    else:
        print(f"\n⚠️  {len(checks_failed)} check(s) failed - investigate before full run")

    return all_passed


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Quick variant sanity check")
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to config.yaml (default: training/config.yaml)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=SANITY_CONFIG["epochs"],
        help=f"Number of epochs (default: {SANITY_CONFIG['epochs']})",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=SANITY_CONFIG["max_samples"],
        help=f"Max training samples (default: {SANITY_CONFIG['max_samples']})",
    )

    args = parser.parse_args()

    # Update config from args
    SANITY_CONFIG["epochs"] = args.epochs
    SANITY_CONFIG["max_samples"] = args.samples

    success = run_sanity_check(args.config)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
