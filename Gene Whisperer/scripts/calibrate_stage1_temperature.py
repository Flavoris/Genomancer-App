#!/usr/bin/env python3
"""
Temperature Scaling Calibration for Stage 1 Model.

Fits a single temperature parameter T > 0 that calibrates model outputs
so that predicted probabilities match empirical frequencies.

Usage:
    python scripts/calibrate_stage1_temperature.py --config config.yaml --ckpt path/to/checkpoint.pt

Output:
    artifacts/calibration_stage1.json with {temperature, calibrated_threshold}
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

# Add training directory to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = SCRIPT_DIR.parent / "training"
ARTIFACTS_DIR = SCRIPT_DIR.parent / "artifacts"
sys.path.insert(0, str(TRAINING_DIR))

import yaml

from dataset import build_dataloaders
from ensemble_infer import build_model, load_checkpoint


def select_device() -> torch.device:
    """Select best available device: CUDA > MPS > CPU."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def infer_model_architecture(checkpoint_path: Path) -> dict:
    """
    Infer model architecture parameters from checkpoint state_dict.

    Returns a dict of config overrides to match the checkpoint's architecture.
    """
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = checkpoint.get("model_state", checkpoint)

    overrides = {}

    # Count encoder layers
    encoder_layers = [k for k in state.keys() if "_full_encoder.layers." in k]
    if encoder_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in encoder_layers)
        overrides["transformer_layers"] = len(layer_nums)

    # Count TCN levels
    tcn_layers = [k for k in state.keys() if "tcn.network." in k]
    if tcn_layers:
        tcn_nums = set(int(k.split("network.")[1].split(".")[0]) for k in tcn_layers)
        overrides["tcn_levels"] = len(tcn_nums)

    # Count multiscale convs (infer kernel sizes from conv weights)
    mscale_convs = [k for k in state.keys() if "multiscale.convs." in k and ".weight" in k]
    if mscale_convs:
        kernel_sizes = []
        for k in sorted(mscale_convs):
            weight = state[k]
            kernel_sizes.append(weight.shape[2])  # Conv1d weight shape: (out_ch, in_ch, kernel)
        overrides["multiscale_kernels"] = kernel_sizes

    # Infer post-CNN transformer layers
    post_cnn_layers = [k for k in state.keys() if "post_cnn_transformer.layers." in k]
    if post_cnn_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in post_cnn_layers)
        overrides["post_cnn_transformer_layers"] = len(layer_nums)

    # Infer embedding dim from embedding weight
    if "embedding.weight" in state:
        overrides["embedding_dim"] = state["embedding.weight"].shape[1]

    # Infer transformer FF dim from FFN weights
    ffn_keys = [k for k in state.keys() if "_full_encoder.layers.0.ffn.0.weight" in k]
    if ffn_keys:
        overrides["transformer_ff_dim"] = state[ffn_keys[0]].shape[0]

    # Infer TCN hidden from TCN output
    tcn_out_keys = [k for k in state.keys() if "tcn.network.0.conv1.conv.weight" in k]
    if tcn_out_keys:
        overrides["tcn_hidden"] = state[tcn_out_keys[0]].shape[0]

    # Infer multiscale channels
    if mscale_convs:
        overrides["multiscale_channels"] = state[sorted(mscale_convs)[0]].shape[0]

    # Infer engineered MLP dimensions
    eng_mlp_keys = [k for k in state.keys() if "engineered_mlp.mlp." in k and "weight" in k]
    if eng_mlp_keys:
        # First linear layer: input -> hidden
        first_linear = [k for k in eng_mlp_keys if "mlp.0.weight" in k]
        if first_linear:
            overrides["engineered_mlp_hidden"] = state[first_linear[0]].shape[0]
        # Last linear layer: hidden -> output
        last_linear = [k for k in eng_mlp_keys if "mlp.8.weight" in k]
        if last_linear:
            overrides["engineered_mlp_output"] = state[last_linear[0]].shape[0]

    # Infer fusion hidden
    fusion_keys = [k for k in state.keys() if "fusion.seq_proj.weight" in k]
    if fusion_keys:
        overrides["fusion_hidden"] = state[fusion_keys[0]].shape[0]

    return overrides


class TemperatureScaler(nn.Module):
    """
    Temperature scaling module that learns a single temperature parameter.

    The temperature T is parameterized as exp(log_temperature) to ensure T > 0.
    Calibrated logits = logits / T
    """

    def __init__(self):
        super().__init__()
        # Initialize log(T) = 0, so T = 1.0 (no scaling initially)
        self.log_temperature = nn.Parameter(torch.zeros(1))

    @property
    def temperature(self) -> float:
        return self.log_temperature.exp().item()

    def forward(self, logits: torch.Tensor) -> torch.Tensor:
        """Apply temperature scaling to logits."""
        temperature = self.log_temperature.exp()
        return logits / temperature


def collect_logits_and_labels(
    model: nn.Module,
    val_loader: torch.utils.data.DataLoader,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Collect all logits and labels from the validation set.

    Returns:
        logits: (N,) tensor of pre-sigmoid logits
        labels: (N,) tensor of binary labels
    """
    model.eval()
    all_logits = []
    all_labels = []

    with torch.no_grad():
        for batch in val_loader:
            tokens, engineered, labels = batch
            tokens = tokens.to(device)
            engineered = engineered.to(device)

            # Get logits via return_logits=True
            _, logits = model(tokens, engineered_features=engineered, return_logits=True)

            all_logits.append(logits.squeeze(-1).cpu())
            all_labels.append(labels)

    logits = torch.cat(all_logits)
    labels = torch.cat(all_labels).float()

    return logits, labels


def fit_temperature(
    logits: torch.Tensor,
    labels: torch.Tensor,
    max_iter: int = 200,
    lr: float = 0.01,
    verbose: bool = True,
) -> float:
    """
    Fit temperature parameter to minimize NLL (BCE with logits).

    We optimize log(T) using LBFGS for fast convergence.

    Args:
        logits: (N,) pre-sigmoid logits
        labels: (N,) binary labels
        max_iter: Maximum optimization iterations
        lr: Learning rate (for LBFGS, this is step size)
        verbose: Print progress

    Returns:
        Optimal temperature T
    """
    scaler = TemperatureScaler()
    criterion = nn.BCEWithLogitsLoss()

    # Use LBFGS for efficient optimization of single parameter
    optimizer = torch.optim.LBFGS(
        scaler.parameters(),
        lr=lr,
        max_iter=max_iter,
        line_search_fn="strong_wolfe",
    )

    def closure():
        optimizer.zero_grad()
        scaled_logits = scaler(logits)
        loss = criterion(scaled_logits, labels)
        loss.backward()
        return loss

    # Initial loss
    with torch.no_grad():
        initial_loss = criterion(logits, labels).item()

    if verbose:
        print(f"Initial NLL (T=1.0): {initial_loss:.6f}")

    # Optimize
    optimizer.step(closure)

    # Final loss
    with torch.no_grad():
        final_loss = criterion(scaler(logits), labels).item()

    temperature = scaler.temperature

    if verbose:
        print(f"Final NLL (T={temperature:.4f}): {final_loss:.6f}")
        print(f"NLL improvement: {(initial_loss - final_loss) * 100:.2f}%")

    return temperature


def compute_metrics_at_threshold(
    probs: np.ndarray, labels: np.ndarray, threshold: float
) -> dict:
    """Compute accuracy at a given threshold."""
    preds = (probs >= threshold).astype(int)
    correct = (preds == labels).sum()
    accuracy = correct / len(labels)
    return {"accuracy": accuracy, "threshold": threshold}


def find_best_threshold(probs: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """
    Find the threshold that maximizes accuracy on the validation set.

    Returns:
        (best_threshold, best_accuracy)
    """
    thresholds = np.linspace(0.01, 0.99, 99)
    best_thr = 0.5
    best_acc = 0.0

    for thr in thresholds:
        metrics = compute_metrics_at_threshold(probs, labels, thr)
        if metrics["accuracy"] > best_acc:
            best_acc = metrics["accuracy"]
            best_thr = thr

    return best_thr, best_acc


def main():
    parser = argparse.ArgumentParser(
        description="Fit temperature scaling for Stage 1 model calibration."
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config (default: config.yaml)",
    )
    parser.add_argument(
        "--ckpt",
        default=None,
        help="Path to Stage 1 checkpoint. If omitted, auto-picks from stage1_ckpt_by_k[kmer].",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output JSON path. Default: artifacts/calibration_stage1.json",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=200,
        help="Maximum LBFGS iterations (default: 200)",
    )
    args = parser.parse_args()

    # Resolve config path
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (TRAINING_DIR / config_path).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    # Load config
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    # Change to training directory for relative paths in config
    os.chdir(TRAINING_DIR)

    # Determine checkpoint path
    ckpt_path = args.ckpt
    if ckpt_path is None:
        kmer = int(cfg.get("kmer", 6))
        ckpt_by_k = cfg.get("stage1_ckpt_by_k", {})
        if kmer in ckpt_by_k:
            ckpt_path = ckpt_by_k[kmer]
        elif str(kmer) in ckpt_by_k:
            ckpt_path = ckpt_by_k[str(kmer)]
        else:
            print(f"ERROR: No checkpoint found for kmer={kmer} in stage1_ckpt_by_k.")
            print(f"Available keys: {list(ckpt_by_k.keys())}")
            print("Use --ckpt to specify a checkpoint path manually.")
            sys.exit(1)

    ckpt_path = Path(ckpt_path)
    if not ckpt_path.is_absolute():
        ckpt_path = (TRAINING_DIR / ckpt_path).resolve()
    if not ckpt_path.exists():
        print(f"ERROR: Checkpoint not found: {ckpt_path}")
        sys.exit(1)

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = ARTIFACTS_DIR / "calibration_stage1.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("TEMPERATURE SCALING CALIBRATION")
    print("=" * 70)
    print(f"Config: {config_path}")
    print(f"Checkpoint: {ckpt_path}")
    print(f"Output: {output_path}")

    # Infer model architecture from checkpoint and apply overrides
    print("\nInferring model architecture from checkpoint...")
    arch_overrides = infer_model_architecture(ckpt_path)
    if arch_overrides:
        print(f"Detected architecture overrides: {arch_overrides}")
        cfg.update(arch_overrides)

    # Select device
    device = select_device()
    print(f"Device: {device}")

    # Build dataloaders and get val loader
    print("\nLoading validation data...")
    loaders = build_dataloaders(cfg)
    val_loader = loaders["stage1"]["val"]
    vocab = val_loader.dataset.vocab

    print(f"Validation samples: {len(val_loader.dataset)}")
    print(f"Vocab size: {len(vocab.itos)}, k={vocab.k}")

    # Build and load model with inferred architecture
    print("\nBuilding model...")
    model = build_model(cfg, vocab, device)
    load_checkpoint(model, ckpt_path, device)
    model.eval()

    # Collect logits and labels
    print("\nCollecting logits from validation set...")
    logits, labels = collect_logits_and_labels(model, val_loader, device)
    print(f"Collected {len(logits)} samples")

    # Statistics before calibration
    probs_before = torch.sigmoid(logits).numpy()
    labels_np = labels.numpy().astype(int)

    print("\n" + "=" * 70)
    print("BEFORE CALIBRATION (T=1.0)")
    print("=" * 70)
    print(f"Probs range: [{probs_before.min():.4f}, {probs_before.max():.4f}]")
    print(f"Probs mean: {probs_before.mean():.4f}")
    print(f"Probs < 0.5: {(probs_before < 0.5).sum()} ({100*(probs_before < 0.5).mean():.1f}%)")
    print(f"Probs >= 0.5: {(probs_before >= 0.5).sum()} ({100*(probs_before >= 0.5).mean():.1f}%)")

    # Find best threshold before calibration
    best_thr_before, best_acc_before = find_best_threshold(probs_before, labels_np)
    acc_at_05_before = compute_metrics_at_threshold(probs_before, labels_np, 0.5)["accuracy"]
    print(f"\nAccuracy @ 0.5: {acc_at_05_before:.4f}")
    print(f"Best threshold: {best_thr_before:.4f} (accuracy: {best_acc_before:.4f})")

    # Fit temperature
    print("\n" + "=" * 70)
    print("FITTING TEMPERATURE")
    print("=" * 70)
    temperature = fit_temperature(logits, labels, max_iter=args.max_iter, verbose=True)

    # Statistics after calibration
    calibrated_logits = logits / temperature
    probs_after = torch.sigmoid(calibrated_logits).numpy()

    print("\n" + "=" * 70)
    print(f"AFTER CALIBRATION (T={temperature:.4f})")
    print("=" * 70)
    print(f"Probs range: [{probs_after.min():.4f}, {probs_after.max():.4f}]")
    print(f"Probs mean: {probs_after.mean():.4f}")
    print(f"Probs < 0.5: {(probs_after < 0.5).sum()} ({100*(probs_after < 0.5).mean():.1f}%)")
    print(f"Probs >= 0.5: {(probs_after >= 0.5).sum()} ({100*(probs_after >= 0.5).mean():.1f}%)")

    # Find best threshold after calibration
    best_thr_after, best_acc_after = find_best_threshold(probs_after, labels_np)
    acc_at_05_after = compute_metrics_at_threshold(probs_after, labels_np, 0.5)["accuracy"]
    print(f"\nAccuracy @ 0.5: {acc_at_05_after:.4f}")
    print(f"Best threshold: {best_thr_after:.4f} (accuracy: {best_acc_after:.4f})")

    # Summary
    print("\n" + "=" * 70)
    print("CALIBRATION SUMMARY")
    print("=" * 70)
    print(f"Temperature: {temperature:.4f}")
    print(f"Accuracy @ 0.5 improvement: {(acc_at_05_after - acc_at_05_before)*100:+.2f}%")

    if abs(best_thr_after - 0.5) < abs(best_thr_before - 0.5):
        print(f"Best threshold moved closer to 0.5: {best_thr_before:.4f} -> {best_thr_after:.4f}")

    # Interpretation
    if temperature > 1.0:
        print(f"\nInterpretation: Model was overconfident (T={temperature:.2f} > 1)")
        print("  Calibration spreads probabilities more evenly around 0.5")
    else:
        print(f"\nInterpretation: Model was underconfident (T={temperature:.2f} < 1)")
        print("  Calibration pushes probabilities toward 0 and 1")

    # Save calibration
    calibration = {
        "temperature": temperature,
        "calibrated_threshold": float(best_thr_after),
        "accuracy_at_05_before": float(acc_at_05_before),
        "accuracy_at_05_after": float(acc_at_05_after),
        "best_threshold_before": float(best_thr_before),
        "best_threshold_after": float(best_thr_after),
        "probs_range_before": [float(probs_before.min()), float(probs_before.max())],
        "probs_range_after": [float(probs_after.min()), float(probs_after.max())],
        "checkpoint": str(ckpt_path),
    }

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(calibration, f, indent=2)

    print(f"\nCalibration saved to: {output_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
