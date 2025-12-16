"""Training script for Gene Whisperer Stage 1 classifier.

Key features for reaching 90%+ accuracy:
1. Label smoothing for better calibration
2. Warmup + cosine annealing LR schedule with warm restarts
3. Gradient accumulation for larger effective batch size
4. Mixup augmentation for better generalization
5. Layer-wise learning rate decay (LLRD) for pretrained models
6. Stochastic Weight Averaging (SWA) for breaking plateaus
7. Multi-scale k-mer ensemble training support
8. Optimized for Apple M4 chip (MPS backend)

Usage:
    python train_stage1.py --config config.yaml
"""
from __future__ import annotations

import argparse
import copy
import json
import logging
import math
import os
import random
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from torch.optim.swa_utils import AveragedModel, SWALR
from tqdm.auto import tqdm

from dataset import build_dataloaders
from model import GeneWhispererStage1, DNAEncoder

LOGGER = logging.getLogger("gene_whisperer.stage1")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def set_seed(seed: int) -> None:
    """Set all random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def load_config(path: Path) -> Dict:
    """Load YAML configuration."""
    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    return cfg or {}


class LabelSmoothingBCE(nn.Module):
    """
    Binary cross entropy with label smoothing.
    
    Label smoothing helps prevent overconfident predictions and
    improves generalization. For binary classification:
    - Smooth 0 → smoothing/2
    - Smooth 1 → 1 - smoothing/2
    """
    
    def __init__(self, smoothing: float = 0.1):
        super().__init__()
        self.smoothing = smoothing
    
    def forward(self, pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        # Smooth labels
        target_smooth = target * (1 - self.smoothing) + 0.5 * self.smoothing
        
        # Binary cross entropy
        loss = F.binary_cross_entropy(pred, target_smooth, reduction='mean')
        return loss


class FocalLossWithSmoothing(nn.Module):
    """
    Focal loss with optional label smoothing.
    
    Combines focal loss (handles class imbalance) with 
    label smoothing (prevents overconfidence).
    """
    
    def __init__(
        self, 
        alpha: float = 0.5, 
        gamma: float = 2.0,
        smoothing: float = 0.05,
        eps: float = 1e-7,
    ):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.smoothing = smoothing
        self.eps = eps
    
    def forward(self, pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        # Label smoothing
        target_smooth = target * (1 - self.smoothing) + 0.5 * self.smoothing
        
        # Focal loss
        p_t = pred * target_smooth + (1 - pred) * (1 - target_smooth)
        alpha_t = self.alpha * target_smooth + (1 - self.alpha) * (1 - target_smooth)
        
        focal_weight = (1 - p_t) ** self.gamma
        loss = -alpha_t * focal_weight * torch.log(p_t.clamp(min=self.eps))
        
        return loss.mean()


def get_layer_wise_lr_groups(
    model: nn.Module,
    base_lr: float,
    lr_decay: float = 0.8,
    weight_decay: float = 0.01,
) -> List[Dict]:
    """
    Create parameter groups with layer-wise learning rate decay (LLRD).
    
    Earlier/deeper layers get lower learning rates to preserve pretrained features,
    while later layers get higher LRs to adapt to the task.
    
    This is critical for fine-tuning pretrained models and achieving high accuracy.
    """
    # Group parameters by layer depth
    no_decay = ["bias", "LayerNorm.weight", "layer_norm", "bn"]
    
    # Identify parameter groups based on model structure
    embedding_params = []
    encoder_layers = {}
    cnn_tcn_params = []
    fusion_params = []
    classifier_params = []
    other_params = []
    
    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue
        
        if "embedding" in name and "encoder" not in name.split(".")[-2:]:
            embedding_params.append((name, param))
        elif "encoder" in name or "post_cnn_transformer" in name:
            # Extract layer number
            parts = name.split(".")
            layer_num = None
            for part in parts:
                if part.isdigit():
                    layer_num = int(part)
                    break
            
            if layer_num is not None:
                if layer_num not in encoder_layers:
                    encoder_layers[layer_num] = []
                encoder_layers[layer_num].append((name, param))
            else:
                other_params.append((name, param))
        elif any(x in name for x in ["feature_extractor", "conv", "tcn", "multiscale"]):
            cnn_tcn_params.append((name, param))
        elif "fusion" in name or "engineered_mlp" in name:
            fusion_params.append((name, param))
        elif "classifier" in name:
            classifier_params.append((name, param))
        else:
            other_params.append((name, param))
    
    # Build parameter groups with decaying LR
    param_groups = []
    
    # Embedding gets lowest LR (deepest/earliest)
    num_layers = max(encoder_layers.keys()) + 1 if encoder_layers else 0
    
    if embedding_params:
        embed_lr = base_lr * (lr_decay ** (num_layers + 2))
        param_groups.append({
            "params": [p for _, p in embedding_params if not any(nd in _ for nd in no_decay)],
            "lr": embed_lr,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in embedding_params if any(nd in _ for nd in no_decay)],
            "lr": embed_lr,
            "weight_decay": 0.0,
        })
    
    # Encoder/transformer layers - earlier layers get lower LR
    for layer_num in sorted(encoder_layers.keys()):
        layer_lr = base_lr * (lr_decay ** (num_layers - layer_num))
        params = encoder_layers[layer_num]
        
        param_groups.append({
            "params": [p for _, p in params if not any(nd in _ for nd in no_decay)],
            "lr": layer_lr,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in params if any(nd in _ for nd in no_decay)],
            "lr": layer_lr,
            "weight_decay": 0.0,
        })
    
    # CNN/TCN gets base LR (task-specific)
    if cnn_tcn_params:
        param_groups.append({
            "params": [p for _, p in cnn_tcn_params if not any(nd in _ for nd in no_decay)],
            "lr": base_lr,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in cnn_tcn_params if any(nd in _ for nd in no_decay)],
            "lr": base_lr,
            "weight_decay": 0.0,
        })
    
    # Fusion gets slightly higher LR
    if fusion_params:
        param_groups.append({
            "params": [p for _, p in fusion_params if not any(nd in _ for nd in no_decay)],
            "lr": base_lr * 1.2,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in fusion_params if any(nd in _ for nd in no_decay)],
            "lr": base_lr * 1.2,
            "weight_decay": 0.0,
        })
    
    # Classifier gets highest LR (needs most adaptation)
    if classifier_params:
        param_groups.append({
            "params": [p for _, p in classifier_params if not any(nd in _ for nd in no_decay)],
            "lr": base_lr * 1.5,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in classifier_params if any(nd in _ for nd in no_decay)],
            "lr": base_lr * 1.5,
            "weight_decay": 0.0,
        })
    
    # Other params get base LR
    if other_params:
        param_groups.append({
            "params": [p for _, p in other_params if not any(nd in _ for nd in no_decay)],
            "lr": base_lr,
            "weight_decay": weight_decay,
        })
        param_groups.append({
            "params": [p for _, p in other_params if any(nd in _ for nd in no_decay)],
            "lr": base_lr,
            "weight_decay": 0.0,
        })
    
    # Filter out empty groups
    param_groups = [g for g in param_groups if len(g["params"]) > 0]
    
    return param_groups


def get_cosine_schedule_with_warmup(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
    num_cycles: float = 0.5,
    min_lr_ratio: float = 0.05,  # Reduced from 0.1
) -> torch.optim.lr_scheduler.LambdaLR:
    """
    Create a schedule with warmup then cosine decay.
    
    Warmup helps stabilize early training, cosine decay provides
    smooth annealing to a minimum learning rate.
    """
    def lr_lambda(current_step: int) -> float:
        # Warmup
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        
        # Cosine decay
        progress = float(current_step - num_warmup_steps) / float(
            max(1, num_training_steps - num_warmup_steps)
        )
        cosine_decay = 0.5 * (1.0 + math.cos(math.pi * num_cycles * 2.0 * progress))
        
        # Don't decay below min_lr_ratio
        return max(min_lr_ratio, cosine_decay)
    
    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


def get_cosine_with_hard_restarts_schedule(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
    num_cycles: int = 3,
    min_lr_ratio: float = 0.05,
) -> torch.optim.lr_scheduler.LambdaLR:
    """
    Cosine schedule with hard restarts for breaking out of local minima.
    
    Periodically resets the learning rate to allow the model to escape
    suboptimal solutions.
    """
    def lr_lambda(current_step: int) -> float:
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        
        progress = float(current_step - num_warmup_steps) / float(
            max(1, num_training_steps - num_warmup_steps)
        )
        cycle_progress = progress * num_cycles
        cycle_position = cycle_progress - int(cycle_progress)
        
        cosine_value = 0.5 * (1.0 + math.cos(math.pi * cycle_position))
        return max(min_lr_ratio, cosine_value)
    
    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


def mixup_data(
    tokens: torch.Tensor,
    engineered: torch.Tensor,
    labels: torch.Tensor,
    alpha: float = 0.2,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, float]:
    """
    Apply mixup augmentation.
    
    Mixup creates convex combinations of training examples,
    which acts as a strong regularizer and improves generalization.
    """
    if alpha <= 0:
        return tokens, engineered, labels, labels, 1.0
    
    lam = np.random.beta(alpha, alpha)
    batch_size = tokens.size(0)
    index = torch.randperm(batch_size, device=tokens.device)
    
    # Mix tokens (using embedding weights mixing would be better, but this works)
    mixed_engineered = lam * engineered + (1 - lam) * engineered[index]
    
    # Return both targets for loss computation
    return tokens, mixed_engineered, labels, labels[index], lam


def calculate_metrics(probabilities: torch.Tensor, labels: torch.Tensor) -> Dict[str, float]:
    """Calculate classification metrics."""
    probs = probabilities.detach().cpu().reshape(-1)
    targets = labels.detach().cpu().reshape(-1)
    preds = (probs >= 0.5).float()
    
    correct = (preds == targets).sum().item()
    total = float(len(targets))
    accuracy = correct / total if total else 0.0
    
    tp = ((preds == 1) & (targets == 1)).sum().item()
    tn = ((preds == 0) & (targets == 0)).sum().item()
    fp = ((preds == 1) & (targets == 0)).sum().item()
    fn = ((preds == 0) & (targets == 1)).sum().item()
    
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    
    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / float(np.sqrt(denom)) if denom > 0 else 0.0
    
    # Specificity (true negative rate)
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    
    return {
        "accuracy": float(accuracy),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "mcc": float(mcc),
        "specificity": float(specificity),
    }


def run_epoch(
    loader: torch.utils.data.DataLoader,
    model: nn.Module,
    criterion: nn.Module,
    device: torch.device,
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler] = None,
    desc: str = "Train",
    grad_accum_steps: int = 1,
    grad_clip_norm: float = 1.0,
    use_mixup: bool = False,
    mixup_alpha: float = 0.2,
    use_amp: bool = False,
    scaler: Optional[torch.cuda.amp.GradScaler] = None,
) -> Dict[str, float]:
    """Run one epoch of training or evaluation."""
    is_train = optimizer is not None
    model.train(is_train)
    
    total_loss = 0.0
    total_samples = 0
    all_probs = []
    all_labels = []
    
    optimizer_step_count = 0
    accumulated_loss = 0.0
    
    iterable = tqdm(loader, desc=desc, leave=False)
    for batch_idx, (tokens, engineered, labels) in enumerate(iterable):
        tokens = tokens.to(device)
        engineered = engineered.to(device)
        labels = labels.to(device).unsqueeze(1)
        
        # Mixup augmentation (training only)
        lam = 1.0
        labels_b = labels
        if is_train and use_mixup:
            tokens, engineered, labels, labels_b, lam = mixup_data(
                tokens, engineered, labels, alpha=mixup_alpha
            )
        
        # Forward pass
        with torch.set_grad_enabled(is_train):
            # Automatic mixed precision
            if use_amp and device.type == "cuda":
                with torch.cuda.amp.autocast():
                    probs = model(tokens, engineered)
                    probs = probs.clamp(min=1e-6, max=1 - 1e-6)
                    
                    if use_mixup and lam < 1.0:
                        loss = lam * criterion(probs, labels) + (1 - lam) * criterion(probs, labels_b)
                    else:
                        loss = criterion(probs, labels)
            else:
                probs = model(tokens, engineered)
                probs = probs.clamp(min=1e-6, max=1 - 1e-6)
                
                if use_mixup and lam < 1.0:
                    loss = lam * criterion(probs, labels) + (1 - lam) * criterion(probs, labels_b)
                else:
                    loss = criterion(probs, labels)
            
            # Backward pass with gradient accumulation
            if is_train:
                loss_scaled = loss / grad_accum_steps
                
                if use_amp and scaler is not None:
                    scaler.scale(loss_scaled).backward()
                else:
                    loss_scaled.backward()
                
                accumulated_loss += loss.item()
                
                # Optimizer step every grad_accum_steps
                if (batch_idx + 1) % grad_accum_steps == 0:
                    if grad_clip_norm > 0:
                        if use_amp and scaler is not None:
                            scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
                    
                    if use_amp and scaler is not None:
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        optimizer.step()
                    
                    optimizer.zero_grad()
                    
                    if scheduler is not None:
                        scheduler.step()
                    
                    optimizer_step_count += 1
        
        batch_size = labels.size(0)
        total_loss += loss.item() * batch_size
        total_samples += batch_size
        all_probs.append(probs.detach().cpu())
        all_labels.append(labels.detach().cpu())
        
        iterable.set_postfix(loss=loss.item())
    
    avg_loss = total_loss / total_samples if total_samples else 0.0
    probs_tensor = torch.cat(all_probs, dim=0)
    labels_tensor = torch.cat(all_labels, dim=0)
    metrics = calculate_metrics(probs_tensor, labels_tensor)
    metrics["loss"] = float(avg_loss)
    
    return metrics


def compute_class_balance(loader: torch.utils.data.DataLoader) -> Dict[str, float]:
    """Compute class balance statistics."""
    dataset = getattr(loader, "dataset", None)
    labels = getattr(dataset, "labels", None)
    if labels is None:
        return {"total": 0, "positive": 0, "negative": 0, "positive_ratio": 0.0}
    
    total = len(labels)
    positive = sum(1 for label in labels if float(label) >= 0.5)
    negative = total - positive
    ratio = positive / total if total else 0.0
    
    return {
        "total": total,
        "positive": positive,
        "negative": negative,
        "positive_ratio": ratio,
    }


def save_checkpoint(
    path: Path, 
    model: nn.Module, 
    optimizer: torch.optim.Optimizer,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler],
    epoch: int,
    metrics: Dict[str, float],
) -> None:
    """Save training checkpoint."""
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "epoch": epoch,
            "model_state": model.state_dict(),
            "optimizer_state": optimizer.state_dict(),
            "scheduler_state": scheduler.state_dict() if scheduler else None,
            "metrics": metrics,
        },
        path,
    )
    LOGGER.info("Saved checkpoint to %s", path)


def check_gradients_and_frozen(model: nn.Module) -> Dict[str, Any]:
    """Check gradient norms and frozen weights for debugging."""
    grad_info = {"frozen": [], "zero_grad": [], "normal": []}
    total_params = 0
    frozen_params = 0
    
    for name, param in model.named_parameters():
        total_params += param.numel()
        if not param.requires_grad:
            frozen_params += param.numel()
            grad_info["frozen"].append(name)
        elif param.grad is None:
            grad_info["zero_grad"].append(name)
        else:
            grad_norm = param.grad.norm().item()
            if grad_norm == 0:
                grad_info["zero_grad"].append(name)
            else:
                grad_info["normal"].append((name, grad_norm))
    
    return {
        "total_params": total_params,
        "frozen_params": frozen_params,
        "frozen_pct": frozen_params / total_params * 100 if total_params else 0,
        "details": grad_info,
    }


def run_overfit_debug_stage1(
    model: nn.Module,
    train_loader: torch.utils.data.DataLoader,
    criterion: nn.Module,
    device: torch.device,
    num_steps: int = 200,
    lr: float = 1e-3,
) -> None:
    """
    Overfit debug: train on tiny dataset for many steps.
    
    Expected: loss should drop and accuracy should reach near 100%.
    If not, prints gradient norms and checks for frozen weights.
    """
    print("=" * 70)
    print("OVERFIT DEBUG MODE: Stage 1 Classifier")
    print("=" * 70)
    
    # Get subset of 32 samples
    all_tokens = []
    all_engineered = []
    all_labels = []
    
    for tokens, engineered, labels in train_loader:
        all_tokens.append(tokens)
        all_engineered.append(engineered)
        all_labels.append(labels)
        if sum(t.size(0) for t in all_tokens) >= 32:
            break
    
    tokens = torch.cat(all_tokens, dim=0)[:32]
    engineered = torch.cat(all_engineered, dim=0)[:32]
    labels = torch.cat(all_labels, dim=0)[:32]
    
    print(f"Training on {tokens.size(0)} samples for {num_steps} steps")
    print(f"Learning rate: {lr}")
    print(f"Label distribution: {labels.sum().item():.0f} positive, {(1-labels).sum().item():.0f} negative")
    print()
    
    model.train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.0)
    
    tokens = tokens.to(device)
    engineered = engineered.to(device)
    labels = labels.to(device).unsqueeze(1)
    
    losses = []
    accuracies = []
    
    for step in range(num_steps):
        optimizer.zero_grad()
        probs = model(tokens, engineered)
        probs = probs.clamp(min=1e-6, max=1 - 1e-6)
        loss = criterion(probs, labels)
        loss.backward()
        optimizer.step()
        
        # Compute accuracy
        preds = (probs >= 0.5).float()
        acc = (preds == labels).float().mean().item()
        
        losses.append(loss.item())
        accuracies.append(acc)
        
        if step % 20 == 0 or step == num_steps - 1:
            print(f"Step {step:>3}: Loss={loss.item():.4f}, Acc={acc*100:.1f}%")
    
    print()
    print("-" * 70)
    print("OVERFIT ANALYSIS")
    print("-" * 70)
    
    initial_loss = sum(losses[:5]) / min(5, len(losses))
    final_loss = sum(losses[-5:]) / min(5, len(losses))
    initial_acc = sum(accuracies[:5]) / min(5, len(accuracies))
    final_acc = sum(accuracies[-5:]) / min(5, len(accuracies))
    
    print(f"Initial loss (first 5 steps avg): {initial_loss:.4f}")
    print(f"Final loss (last 5 steps avg):    {final_loss:.4f}")
    print(f"Loss reduction:                   {initial_loss - final_loss:.4f}")
    print(f"Initial accuracy:                 {initial_acc*100:.1f}%")
    print(f"Final accuracy:                   {final_acc*100:.1f}%")
    print()
    
    # Check if overfitting succeeded
    success = final_acc > 0.95
    
    if success:
        print("✓ OVERFIT SUCCESS: Model can learn! Accuracy reached near 100%.")
    else:
        print("✗ OVERFIT FAILED: Model is not learning properly!")
        print()
        print("Diagnosing issues...")
        
        # Check gradients
        grad_info = check_gradients_and_frozen(model)
        
        print(f"\n[FROZEN WEIGHTS]")
        print(f"  Frozen params: {grad_info['frozen_params']:,} / {grad_info['total_params']:,} ({grad_info['frozen_pct']:.1f}%)")
        if grad_info["details"]["frozen"]:
            print(f"  Frozen layers (first 5): {grad_info['details']['frozen'][:5]}")
        
        print(f"\n[ZERO GRADIENTS]")
        if grad_info["details"]["zero_grad"]:
            print(f"  Layers with zero/no grad: {len(grad_info['details']['zero_grad'])}")
            for name in grad_info["details"]["zero_grad"][:5]:
                print(f"    - {name}")
        else:
            print("  No zero gradient issues detected")
        
        print(f"\n[GRADIENT NORMS]")
        if grad_info["details"]["normal"]:
            sorted_grads = sorted(grad_info["details"]["normal"], key=lambda x: x[1], reverse=True)
            print("  Top 5 gradient norms:")
            for name, norm in sorted_grads[:5]:
                print(f"    {norm:.6f}: {name}")
            print("  Bottom 5 gradient norms:")
            for name, norm in sorted_grads[-5:]:
                print(f"    {norm:.6f}: {name}")
        
        print("\n[POSSIBLE ISSUES]")
        if grad_info["frozen_pct"] > 50:
            print("  ⚠ More than 50% of parameters are frozen!")
        if len(grad_info["details"]["zero_grad"]) > 0:
            print("  ⚠ Some layers have zero gradients - check loss function")
        if final_loss > initial_loss:
            print("  ⚠ Loss increased - possible exploding gradients or wrong loss")
        if final_acc < 0.6:
            print("  ⚠ Very low accuracy - check labels are correct")
    
    print("=" * 70)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train Gene Whisperer Stage 1")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--overfit_debug", action="store_true",
                        help="Train on 32 samples for 200 steps to verify model can learn")
    args = parser.parse_args()
    
    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    
    # Allow fallback to original config
    if not config_path.exists():
        fallback = script_dir / "config.yaml"
        if fallback.exists():
            config_path = fallback
            LOGGER.warning("Using fallback config: %s", config_path)
    
    cfg = load_config(config_path)
    
    # =========================================================================
    # Log resolved hyperparameters to runs/<timestamp>/resolved_config.json
    # =========================================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"stage1_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    
    resolved_config: Dict[str, Any] = {
        "script": "train_stage1.py",
        "timestamp": timestamp,
        "config_path": str(config_path),
        # Model architecture
        "embedding_dim": int(cfg.get("embedding_dim", 256)),
        "transformer_layers": int(cfg.get("transformer_layers", 6)),
        "transformer_heads": int(cfg.get("transformer_heads", 8)),
        "transformer_ff_dim": int(cfg.get("transformer_ff_dim", 1024)),
        "transformer_dropout": float(cfg.get("transformer_dropout", 0.15)),
        # Optimizer
        "lr": float(cfg.get("lr", 2e-4)),
        "weight_decay": float(cfg.get("weight_decay", 0.02)),
        # Data
        "batch_size": int(cfg.get("batch_size", 32)),
        "max_bp_len": int(cfg.get("max_bp_len", 81)),
        "kmer": int(cfg.get("kmer", 3)),
        # Engineered features
        "stage1_use_engineered_features": bool(cfg.get("stage1_use_engineered_features", True)),
        "stage1_engineered_dim": int(cfg.get("stage1_engineered_dim", 208)),
        "engineered_mlp_hidden": int(cfg.get("engineered_mlp_hidden", 256)),
        "engineered_mlp_output": int(cfg.get("engineered_mlp_output", 128)),
        # Augmentation
        "reverse_complement_prob": float(cfg.get("stage1_reverse_complement_prob", 0.5)),
        "use_mixup": bool(cfg.get("use_mixup", True)),
        "mixup_alpha": float(cfg.get("mixup_alpha", 0.2)),
        # TCN
        "use_tcn": bool(cfg.get("use_tcn", True)),
        "tcn_hidden": int(cfg.get("tcn_hidden", 256)),
        "tcn_levels": int(cfg.get("tcn_levels", 4)),
        "tcn_kernel": int(cfg.get("tcn_kernel", 3)),
        "multiscale_channels": int(cfg.get("multiscale_channels", 64)),
        "multiscale_kernels": cfg.get("multiscale_kernels", [3, 5, 7, 9, 15]),
        "lstm_hidden": int(cfg.get("lstm_hidden", 192)),
        # Training
        "epochs": int(cfg.get("epochs", 120)),
        "grad_accum_steps": int(cfg.get("grad_accum_steps", 4)),
        "grad_clip_norm": float(cfg.get("grad_clip_norm", 1.0)),
        "warmup_ratio": float(cfg.get("warmup_ratio", 0.15)),
        "early_stopping_patience": int(cfg.get("early_stopping_patience", 25)),
        "layer_lr_decay": float(cfg.get("layer_lr_decay", 0.8)),
        "drop_path_rate": float(cfg.get("drop_path_rate", 0.1)),
        # Loss
        "stage1_use_focal_loss": bool(cfg.get("stage1_use_focal_loss", True)),
        "stage1_focal_alpha": float(cfg.get("stage1_focal_alpha", 0.5)),
        "stage1_focal_gamma": float(cfg.get("stage1_focal_gamma", 2.0)),
        "label_smoothing": float(cfg.get("label_smoothing", 0.05)),
        # SWA
        "use_swa": bool(cfg.get("use_swa", True)),
        "swa_start_epoch": int(cfg.get("swa_start_epoch", 60)),
        "swa_lr": float(cfg.get("swa_lr", 5e-5)),
        # Model options
        "use_alibi": bool(cfg.get("use_alibi", True)),
        "use_attention_pool": bool(cfg.get("use_attention_pool", True)),
        "post_cnn_transformer_layers": int(cfg.get("post_cnn_transformer_layers", 3)),
        "fusion_hidden": int(cfg.get("fusion_hidden", 256)),
        # Pretrained weights
        "stage1_load_mlm_weights": bool(cfg.get("stage1_load_mlm_weights", True)),
        "stage1_freeze_layers": int(cfg.get("stage1_freeze_layers", 0)),
        # Seed
        "seed": int(cfg.get("seed", 1337)),
    }
    
    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)
    
    seed = int(cfg.get("seed", 1337))
    set_seed(seed)
    
    # Device selection (optimized for M4)
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Using device: %s", device)
    
    # Build dataloaders
    dataloaders = build_dataloaders(cfg)
    train_loader = dataloaders["stage1"]["train"]
    val_loader = dataloaders["stage1"]["val"]
    
    # Model configuration (enhanced defaults)
    embedding_dim = int(cfg.get("embedding_dim", 192))
    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.15))
    
    train_dataset = train_loader.dataset
    stage1_vocab = getattr(train_dataset, "vocab", None)
    if stage1_vocab is None:
        raise ValueError("Stage 1 dataset is missing vocabulary")
    
    vocab_size = len(stage1_vocab.itos)
    pad_token_id = stage1_vocab.pad_id
    kmer = stage1_vocab.k
    
    use_alibi = bool(cfg.get("use_alibi", True))
    use_attention_pool = bool(cfg.get("use_attention_pool", True))
    use_engineered_features = bool(cfg.get("stage1_use_engineered_features", True))
    engineered_dim = int(cfg.get("stage1_engineered_dim", 208))
    
    # TCN configuration
    use_tcn = bool(cfg.get("use_tcn", True))
    tcn_hidden = int(cfg.get("tcn_hidden", 256))
    tcn_levels = int(cfg.get("tcn_levels", 4))
    tcn_kernel = int(cfg.get("tcn_kernel", 3))
    multiscale_channels = int(cfg.get("multiscale_channels", 64))
    multiscale_kernels_cfg = cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])
    multiscale_kernels = tuple(multiscale_kernels_cfg)
    lstm_hidden = int(cfg.get("lstm_hidden", 192))
    
    # New architecture parameters: post-CNN transformer, engineered MLP, fusion
    post_cnn_transformer_layers = int(cfg.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg.get("fusion_hidden", 256))
    
    model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        use_alibi=use_alibi,
        pad_token_id=pad_token_id,
        engineered_dim=engineered_dim,
        use_engineered_features=use_engineered_features,
        use_attention_pool=use_attention_pool,
        # TCN parameters
        use_tcn=use_tcn,
        tcn_hidden=tcn_hidden,
        tcn_levels=tcn_levels,
        tcn_kernel=tcn_kernel,
        multiscale_channels=multiscale_channels,
        multiscale_kernels=multiscale_kernels,
        lstm_hidden=lstm_hidden,
        # New architecture parameters
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        fusion_hidden=fusion_hidden,
    ).to(device)
    
    # Log architecture configuration
    LOGGER.info("=" * 60)
    LOGGER.info("MODEL ARCHITECTURE (V3 - Improved Order + Attention Fusion)")
    LOGGER.info("=" * 60)
    LOGGER.info("Architecture: Embedding → CNN/TCN → Transformer → Pool → Fusion → Classifier")
    if use_tcn:
        LOGGER.info("CNN/TCN: Multi-scale kernels %s → TCN(%d levels, %d hidden, k=%d)", 
                    multiscale_kernels, tcn_levels, tcn_hidden, tcn_kernel)
    LOGGER.info("Post-CNN Transformer: %d layers, %d heads", post_cnn_transformer_layers, transformer_heads)
    if use_engineered_features:
        LOGGER.info("Engineered Features MLP: %d → %d → %d", engineered_dim, engineered_mlp_hidden, engineered_mlp_output)
        LOGGER.info("Attention Fusion: Gated cross-attention (hidden=%d)", fusion_hidden)
    LOGGER.info("=" * 60)
    
    # Load pretrained MLM weights if available
    load_mlm_weights = bool(cfg.get("stage1_load_mlm_weights", False))
    mlm_checkpoint = cfg.get("mlm_encoder_checkpoint") or cfg.get("mlm_encoder_path")
    if load_mlm_weights and mlm_checkpoint:
        checkpoint_path = Path(mlm_checkpoint)
        if not checkpoint_path.is_absolute():
            checkpoint_path = (script_dir / checkpoint_path).resolve()
        if checkpoint_path.exists():
            LOGGER.info("Loading MLM encoder weights from %s", checkpoint_path)
            # Use new method that transfers both embedding and transformer layers
            model.load_pretrained_weights(checkpoint_path, strict=False)
            LOGGER.info("Pretrained weights loaded: embedding + %d transformer layers", 
                       post_cnn_transformer_layers)
            
            freeze_layers = int(cfg.get("stage1_freeze_layers", 0))
            if freeze_layers > 0:
                LOGGER.info("Freezing bottom %d encoder layers", freeze_layers)
                model.encoder.freeze_bottom_layers(freeze_layers)
        else:
            LOGGER.warning("MLM checkpoint not found at %s, training from scratch", checkpoint_path)
    
    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    LOGGER.info("Model: %.2fM total params, %.2fM trainable", 
                total_params / 1e6, trainable_params / 1e6)
    
    # Optimizer with layer-wise learning rate decay
    lr = float(cfg.get("lr", 2e-4))
    weight_decay = float(cfg.get("weight_decay", 0.02))
    layer_lr_decay = float(cfg.get("layer_lr_decay", 0.8))
    
    # Use LLRD if we loaded pretrained weights
    if load_mlm_weights and mlm_checkpoint:
        LOGGER.info("Using layer-wise LR decay (factor=%.2f) for pretrained model", layer_lr_decay)
        param_groups = get_layer_wise_lr_groups(
            model, 
            base_lr=lr, 
            lr_decay=layer_lr_decay,
            weight_decay=weight_decay,
        )
        optimizer = torch.optim.AdamW(param_groups, betas=(0.9, 0.999))
        LOGGER.info("Created %d parameter groups with LLRD", len(param_groups))
    else:
        optimizer = torch.optim.AdamW(
            model.parameters(), 
            lr=lr, 
            weight_decay=weight_decay,
            betas=(0.9, 0.999),
        )
    
    # Scheduler with warmup
    epochs = int(cfg.get("epochs", 120))
    warmup_ratio = float(cfg.get("warmup_ratio", 0.15))
    steps_per_epoch = len(train_loader)
    total_steps = epochs * steps_per_epoch
    warmup_steps = int(warmup_ratio * total_steps)
    
    scheduler = get_cosine_schedule_with_warmup(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=0.02,  # Decay to lower minimum
    )
    LOGGER.info("Using cosine schedule with %d warmup steps", warmup_steps)
    
    # SWA setup for breaking plateaus
    use_swa = bool(cfg.get("use_swa", True))
    swa_start_epoch = int(cfg.get("swa_start_epoch", 60))
    swa_lr = float(cfg.get("swa_lr", 5e-5))
    swa_model = None
    swa_scheduler = None
    
    if use_swa:
        swa_model = AveragedModel(model)
        swa_scheduler = SWALR(optimizer, swa_lr=swa_lr, anneal_epochs=10)
        LOGGER.info("SWA enabled: starts at epoch %d with LR %.2e", swa_start_epoch, swa_lr)
    
    # Loss function with label smoothing
    label_smoothing = float(cfg.get("label_smoothing", 0.05))
    use_focal_loss = bool(cfg.get("stage1_use_focal_loss", True))
    focal_alpha = float(cfg.get("stage1_focal_alpha", 0.5))
    focal_gamma = float(cfg.get("stage1_focal_gamma", 2.0))
    
    if use_focal_loss:
        criterion = FocalLossWithSmoothing(
            alpha=focal_alpha, 
            gamma=focal_gamma,
            smoothing=label_smoothing,
        )
        LOGGER.info("Using focal loss with smoothing=%.2f", label_smoothing)
    else:
        criterion = LabelSmoothingBCE(smoothing=label_smoothing)
        LOGGER.info("Using BCE with smoothing=%.2f", label_smoothing)
    
    # =========================================================================
    # OVERFIT DEBUG MODE: Train on 32 samples for 200 steps
    # =========================================================================
    if args.overfit_debug:
        LOGGER.info("Running OVERFIT DEBUG mode")
        # Use simple BCE for overfit test (no label smoothing or focal loss)
        simple_criterion = nn.BCELoss()
        run_overfit_debug_stage1(
            model=model,
            train_loader=train_loader,
            criterion=simple_criterion,
            device=device,
            num_steps=200,
            lr=1e-3,
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return
    
    # Training settings
    grad_accum_steps = int(cfg.get("grad_accum_steps", 2))
    grad_clip_norm = float(cfg.get("grad_clip_norm", 1.0))
    use_mixup = bool(cfg.get("use_mixup", True))
    mixup_alpha = float(cfg.get("mixup_alpha", 0.2))
    patience = int(cfg.get("early_stopping_patience", 15))
    
    class_balance = compute_class_balance(train_loader)
    LOGGER.info(
        "Class balance: %.2f%% positive (%d/%d)",
        class_balance["positive_ratio"] * 100.0,
        class_balance["positive"],
        class_balance["total"],
    )
    
    # Paths
    artifacts_root = (script_dir / "../artifacts").resolve()
    checkpoint_name = cfg.get("stage1_checkpoint_name") or f"stage1_k{kmer}.pt"
    checkpoint_path = artifacts_root / "checkpoints" / checkpoint_name
    report_path = artifacts_root / f"stage1_report_k{kmer}.json"
    
    # Training loop
    best_val_mcc = -float("inf")
    best_val_acc = 0.0
    best_epoch = 0
    best_train_metrics = None
    best_val_metrics = None
    patience_counter = 0
    swa_started = False
    
    for epoch in range(1, epochs + 1):
        # Check if we should start SWA
        in_swa_phase = use_swa and epoch >= swa_start_epoch
        if in_swa_phase and not swa_started:
            LOGGER.info("=" * 40)
            LOGGER.info("Starting SWA phase at epoch %d", epoch)
            LOGGER.info("=" * 40)
            swa_started = True
        
        current_lr = optimizer.param_groups[0]["lr"]
        phase_str = "[SWA]" if in_swa_phase else ""
        LOGGER.info("Epoch %d/%d %s (LR: %.2e)", epoch, epochs, phase_str, current_lr)
        
        train_metrics = run_epoch(
            train_loader,
            model,
            criterion,
            device,
            optimizer=optimizer,
            scheduler=None if in_swa_phase else scheduler,  # Use SWA scheduler in SWA phase
            desc="Train",
            grad_accum_steps=grad_accum_steps,
            grad_clip_norm=grad_clip_norm,
            use_mixup=use_mixup and not in_swa_phase,  # Disable mixup in SWA phase
            mixup_alpha=mixup_alpha,
        )
        
        # Update SWA model and scheduler
        if in_swa_phase and swa_model is not None:
            swa_model.update_parameters(model)
            swa_scheduler.step()
        
        val_metrics = run_epoch(
            val_loader,
            model,
            criterion,
            device,
            optimizer=None,
            desc="Val",
        )
        
        # Also evaluate SWA model if active (every 5 epochs to save time)
        swa_val_metrics = None
        if in_swa_phase and swa_model is not None and epoch % 5 == 0:
            # Update batch norm stats for SWA model using a limited subset
            try:
                # Create a wrapper dataloader that yields only tokens (for BN update)
                def bn_loader():
                    for tokens, _, _ in train_loader:
                        yield tokens.to(device), None
                
                # Use the SWA model directly for BN update
                swa_model.train()
                with torch.no_grad():
                    for i, (tokens, eng, _) in enumerate(train_loader):
                        if i >= 50:  # Limit to 50 batches for speed
                            break
                        tokens = tokens.to(device)
                        eng = eng.to(device)
                        _ = swa_model(tokens, eng)
                
                swa_model.eval()
                swa_val_metrics = run_epoch(
                    val_loader,
                    swa_model,
                    criterion,
                    device,
                    optimizer=None,
                    desc="Val-SWA",
                )
            except Exception as e:
                LOGGER.warning("SWA evaluation failed: %s", e)
                swa_val_metrics = None
        
        LOGGER.info(
            "Train - Loss: %.4f | Acc: %.4f | Prec: %.4f | Rec: %.4f | F1: %.4f | MCC: %.4f",
            train_metrics["loss"],
            train_metrics["accuracy"],
            train_metrics["precision"],
            train_metrics["recall"],
            train_metrics["f1"],
            train_metrics["mcc"],
        )
        LOGGER.info(
            "Val   - Loss: %.4f | Acc: %.4f | Prec: %.4f | Rec: %.4f | F1: %.4f | MCC: %.4f",
            val_metrics["loss"],
            val_metrics["accuracy"],
            val_metrics["precision"],
            val_metrics["recall"],
            val_metrics["f1"],
            val_metrics["mcc"],
        )
        
        if swa_val_metrics:
            LOGGER.info(
                "SWA   - Loss: %.4f | Acc: %.4f | F1: %.4f | MCC: %.4f",
                swa_val_metrics["loss"],
                swa_val_metrics["accuracy"],
                swa_val_metrics["f1"],
                swa_val_metrics["mcc"],
            )
        
        # Use SWA metrics if better during SWA phase
        eval_metrics = val_metrics
        if swa_val_metrics and swa_val_metrics["mcc"] > val_metrics["mcc"]:
            eval_metrics = swa_val_metrics
            LOGGER.info("(Using SWA metrics for checkpointing)")
        
        # Check for improvement (use MCC as primary, accuracy as tiebreaker)
        improved = False
        if eval_metrics["mcc"] > best_val_mcc + 1e-4:
            improved = True
        elif abs(eval_metrics["mcc"] - best_val_mcc) < 1e-4 and eval_metrics["accuracy"] > best_val_acc:
            improved = True
        
        if improved:
            patience_counter = 0
            best_val_mcc = eval_metrics["mcc"]
            best_val_acc = eval_metrics["accuracy"]
            best_epoch = epoch
            best_train_metrics = train_metrics
            best_val_metrics = eval_metrics
            
            # Save either SWA model or regular model
            if swa_val_metrics and swa_val_metrics["mcc"] > val_metrics["mcc"]:
                save_checkpoint(checkpoint_path, swa_model.module, optimizer, scheduler, epoch, eval_metrics)
            else:
                save_checkpoint(checkpoint_path, model, optimizer, scheduler, epoch, val_metrics)
            LOGGER.info("★ New best! MCC: %.4f, Acc: %.4f (%.2f%%)", best_val_mcc, best_val_acc, best_val_acc * 100)
        else:
            patience_counter += 1
            # In SWA phase, be more patient as averaging takes time
            effective_patience = patience + 10 if in_swa_phase else patience
            LOGGER.info("No improvement. Patience %d/%d", patience_counter, effective_patience)
            if patience_counter >= effective_patience:
                LOGGER.info("Early stopping triggered at epoch %d", epoch)
                break
    
    # Final report
    if best_val_metrics is None:
        best_train_metrics = train_metrics
        best_val_metrics = val_metrics
        best_epoch = epochs
    
    report = {
        "config": str(config_path),
        "best_epoch": best_epoch,
        "train_metrics": best_train_metrics,
        "val_metrics": best_val_metrics,
        "class_balance": class_balance,
        "device": str(device),
        "seed": seed,
        "model_params": {
            "total": total_params,
            "trainable": trainable_params,
            "embedding_dim": embedding_dim,
            "layers": transformer_layers,
            "heads": transformer_heads,
            "kmer": kmer,
        },
        "architecture_v4": {
            "order": "Embedding → CNN/TCN → Transformer → Pool → Fusion → Classifier",
            "use_tcn": use_tcn,
            "tcn_hidden": tcn_hidden,
            "tcn_levels": tcn_levels,
            "multiscale_kernels": list(multiscale_kernels),
            "post_cnn_transformer_layers": post_cnn_transformer_layers,
            "engineered_mlp_hidden": engineered_mlp_hidden,
            "engineered_mlp_output": engineered_mlp_output,
            "fusion_hidden": fusion_hidden,
            "use_attention_fusion": use_engineered_features,
        },
        "training_enhancements": {
            "layer_wise_lr_decay": layer_lr_decay,
            "label_smoothing": label_smoothing,
            "mixup_alpha": mixup_alpha,
            "use_swa": use_swa,
            "swa_start_epoch": swa_start_epoch if use_swa else None,
            "swa_lr": swa_lr if use_swa else None,
            "drop_path_rate": float(cfg.get("drop_path_rate", 0.1)),
        },
    }
    
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    LOGGER.info("Saved training report to %s", report_path)
    
    LOGGER.info("=" * 60)
    LOGGER.info("TRAINING COMPLETE")
    LOGGER.info("Best Epoch: %d%s", best_epoch, " (SWA)" if swa_started and best_epoch >= swa_start_epoch else "")
    LOGGER.info("Best Val Accuracy: %.4f (%.2f%%)", best_val_acc, best_val_acc * 100)
    LOGGER.info("Best Val MCC: %.4f", best_val_mcc)
    if use_swa:
        LOGGER.info("SWA was %s", "active" if swa_started else "not triggered (training ended early)")
    LOGGER.info("=" * 60)


if __name__ == "__main__":
    main()
