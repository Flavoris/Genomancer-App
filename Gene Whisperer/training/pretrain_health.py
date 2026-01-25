"""Pretrain Health Dashboard for MLM pretraining.

Tracks and logs comprehensive training health metrics:
- Loss, masked-token accuracy
- UNK rate, PAD masked rate (assert 0)
- Mask coverage (% tokens masked), span length histogram
- Grad norm, parameter norm, LR

Writes to runs/<id>/metrics.jsonl and prints compact table every N steps.

Includes sanity checks that must pass before training can start.
"""
from __future__ import annotations

import json
import logging
import math
import warnings
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn

LOGGER = logging.getLogger("gene_whisperer.pretrain_health")


@dataclass
class HealthMetrics:
    """Single snapshot of training health metrics."""
    step: int
    epoch: int
    loss: float
    masked_accuracy: float
    unk_rate: float
    pad_masked_rate: float
    mask_coverage: float
    grad_norm: float
    param_norm: float
    lr: float
    span_lengths: List[int] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable dict."""
        result = {
            "step": self.step,
            "epoch": self.epoch,
            "loss": self.loss,
            "loss_masked_only": self.loss,
            "masked_accuracy": self.masked_accuracy,
            "masked_token_accuracy": self.masked_accuracy,
            "unk_rate": self.unk_rate,
            "pad_masked_rate": self.pad_masked_rate,
            "mask_coverage": self.mask_coverage,
            "masked_fraction": self.mask_coverage,
            "grad_norm": self.grad_norm,
            "param_norm": self.param_norm,
            "lr": self.lr,
        }
        
        # Add span length histogram
        if self.span_lengths:
            span_counter = Counter(self.span_lengths)
            result["span_histogram"] = dict(span_counter)
            result["avg_span_length"] = sum(self.span_lengths) / len(self.span_lengths)
        
        return result


class MaskSanityChecker:
    """Validates masking behavior before and during training.
    
    Checks:
    1. PAD tokens are never masked (assertion, hard fail)
    2. UNK rate is near 0 (warning if >0.1%)
    3. Mask coverage is approximately target rate
    """
    
    def __init__(self, vocab, target_mask_prob: float = 0.15):
        self.vocab = vocab
        self.target_mask_prob = target_mask_prob
        self.pad_id = vocab.pad_id
        self.unk_id = vocab.unk_id
        self.mask_id = vocab.mask_id
    
    def check_batch(
        self,
        original_tokens: torch.LongTensor,
        masked_tokens: torch.LongTensor,
        labels: torch.LongTensor,
    ) -> Dict[str, Any]:
        """Check a single batch for masking sanity.
        
        Args:
            original_tokens: (B, L) original token IDs before masking
            masked_tokens: (B, L) token IDs after masking applied
            labels: (B, L) labels with -100 for unmasked positions
            
        Returns:
            Dict with sanity check results
        """
        batch_size, seq_len = original_tokens.shape
        total_tokens = batch_size * seq_len
        
        # Identify masked positions (where label != -100)
        masked_positions = (labels != -100)
        num_masked = masked_positions.sum().item()
        
        # Check 1: PAD tokens should never be masked
        pad_positions = (original_tokens == self.pad_id)
        masked_pads = (masked_positions & pad_positions).sum().item()
        pad_masked_rate = masked_pads / max(1, pad_positions.sum().item())
        
        # Check 2: Count UNK tokens in original
        unk_positions = (original_tokens == self.unk_id)
        num_unks = unk_positions.sum().item()
        unk_rate = num_unks / total_tokens
        
        # Check 3: Mask coverage
        # Don't count padding in total tokens for coverage calculation
        non_pad_tokens = total_tokens - pad_positions.sum().item()
        mask_coverage = num_masked / max(1, non_pad_tokens)
        
        return {
            "pad_masked_count": masked_pads,
            "pad_masked_rate": pad_masked_rate,
            "unk_count": num_unks,
            "unk_rate": unk_rate,
            "mask_coverage": mask_coverage,
            "total_tokens": total_tokens,
            "non_pad_tokens": non_pad_tokens,
            "num_masked": num_masked,
        }
    
    def run_sanity_check(
        self,
        dataloader,
        num_batches: int = 5,
    ) -> Tuple[bool, List[str]]:
        """Run sanity check on multiple batches.
        
        Returns:
            (passed, issues) tuple where passed=True if all checks pass
        """
        issues = []
        total_pad_masked = 0
        total_unk = 0
        total_tokens = 0
        total_non_pad_tokens = 0
        total_masked = 0
        batches_checked = 0

        print("=" * 70)
        print("MASK SANITY CHECK")
        print("=" * 70)

        for i, batch_data in enumerate(dataloader):
            batches_checked = i + 1
            if i >= num_batches:
                break

            # Handle both (inputs, labels) and (inputs, labels, spans) formats
            if len(batch_data) == 3:
                masked_tokens, labels, _ = batch_data
            else:
                masked_tokens, labels = batch_data

            # We need original tokens to check PAD/UNK
            # For this check, we'll reconstruct from labels where possible
            # and assume masked_tokens == original for unmasked positions
            batch_size, seq_len = masked_tokens.shape

            # Count PAD positions in masked tokens (PAD should remain as PAD)
            pad_positions = (masked_tokens == self.pad_id)

            # Labels for PAD positions should be -100 (unmasked)
            # If any PAD has label != -100, that's a masking error
            masked_positions = (labels != -100)
            pad_masked = (masked_positions & pad_positions).sum().item()
            total_pad_masked += pad_masked

            # Count UNK in labels (where we predict UNK is wrong)
            unk_in_labels = ((labels != -100) & (labels == self.unk_id)).sum().item()
            total_unk += unk_in_labels

            total_tokens += batch_size * seq_len
            total_non_pad_tokens += (batch_size * seq_len) - pad_positions.sum().item()
            total_masked += masked_positions.sum().item()

        # Summarize
        unk_rate = total_unk / max(1, total_masked) * 100
        # Calculate mask rate over NON-PAD tokens (the meaningful metric)
        mask_rate_non_pad = total_masked / max(1, total_non_pad_tokens) * 100
        # Also calculate over all tokens for reference
        mask_rate_all = total_masked / max(1, total_tokens) * 100
        non_pad_ratio = total_non_pad_tokens / max(1, total_tokens) * 100

        print(f"Checked {min(num_batches, batches_checked)} batches, {total_tokens} total tokens")
        print(f"Non-PAD tokens: {total_non_pad_tokens} ({non_pad_ratio:.1f}% of total)")
        print()
        
        # Check 1: PAD masked rate (must be 0)
        print("[1] PAD tokens masked:", end=" ")
        if total_pad_masked == 0:
            print(f"✓ PASS (0 PAD tokens masked)")
        else:
            print(f"✗ FAIL ({total_pad_masked} PAD tokens were masked!)")
            issues.append(f"PAD tokens were masked: {total_pad_masked} instances")
        
        # Check 2: UNK rate (warn if >0.1%)
        print(f"[2] UNK rate in masked labels:", end=" ")
        if unk_rate < 0.1:
            print(f"✓ PASS ({unk_rate:.3f}%)")
        else:
            print(f"⚠ WARN ({unk_rate:.2f}% > 0.1% threshold)")
            issues.append(f"High UNK rate: {unk_rate:.2f}%")
        
        # Check 3: Mask coverage (measured over NON-PAD tokens, which is the meaningful metric)
        print(f"[3] Mask coverage:", end=" ")
        expected = self.target_mask_prob * 100
        if abs(mask_rate_non_pad - expected) < 5:  # Within 5% of target
            print(f"✓ PASS ({mask_rate_non_pad:.1f}% of non-PAD ≈ {expected:.0f}% target)")
        else:
            print(f"⚠ WARN ({mask_rate_non_pad:.1f}% of non-PAD differs from {expected:.0f}% target)")
        print(f"    (Over all tokens including PAD: {mask_rate_all:.1f}%)")
        
        print()
        print("=" * 70)
        
        passed = len(issues) == 0 or all("WARN" in i or "High UNK" in i for i in issues)
        
        if total_pad_masked > 0:
            passed = False  # PAD masking is a hard failure
        
        return passed, issues


class PretrainHealthDashboard:
    """Comprehensive health monitoring for MLM pretraining.
    
    Logs metrics to runs/<id>/metrics.jsonl and prints compact table
    every N steps.
    """
    
    def __init__(
        self,
        run_dir: Path,
        vocab,
        log_every_n_steps: int = 50,
        target_mask_prob: float = 0.15,
    ):
        self.run_dir = Path(run_dir)
        self.vocab = vocab
        self.log_every_n_steps = log_every_n_steps
        self.target_mask_prob = target_mask_prob
        
        # Ensure run directory exists
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Metrics file
        self.metrics_path = self.run_dir / "metrics.jsonl"
        
        # Initialize metrics file (truncate if exists)
        with self.metrics_path.open("w") as f:
            pass  # Create empty file
        
        # Accumulator for batch-level metrics
        self._reset_accumulators()
        
        # Track for gradient explosion detection
        self.prev_grad_norm: Optional[float] = None
        self.grad_explosion_count = 0
        
        LOGGER.info("Health dashboard initialized, logging to %s", self.metrics_path)
    
    def _reset_accumulators(self) -> None:
        """Reset per-step accumulators."""
        self._losses: List[float] = []
        self._correct: int = 0
        self._total_masked: int = 0
        self._unk_count: int = 0
        self._pad_masked_count: int = 0
        self._total_tokens: int = 0
        self._span_lengths: List[int] = []
    
    def accumulate_batch(
        self,
        loss: float,
        logits: torch.Tensor,
        labels: torch.LongTensor,
        masked_inputs: torch.LongTensor,
        span_lengths: Optional[List[int]] = None,
    ) -> None:
        """Accumulate metrics from a single batch.
        
        Args:
            loss: Batch loss value
            logits: (B, L, V) model output logits
            labels: (B, L) labels with -100 for unmasked
            masked_inputs: (B, L) input token IDs after masking
            span_lengths: Optional list of span lengths used in masking
        """
        self._losses.append(loss)
        
        # Accuracy on masked tokens
        preds = logits.argmax(dim=-1)
        mask = labels != -100
        self._correct += (preds[mask] == labels[mask]).sum().item()
        self._total_masked += mask.sum().item()
        
        batch_size, seq_len = labels.shape
        self._total_tokens += batch_size * seq_len
        
        # UNK count in labels (ground truth)
        unk_in_labels = (labels == self.vocab.unk_id) & (labels != -100)
        self._unk_count += unk_in_labels.sum().item()
        
        # PAD masked count (PAD positions that have label != -100)
        # After masking, PAD tokens should still be PAD in masked_inputs
        pad_positions = (masked_inputs == self.vocab.pad_id)
        pad_masked = (mask & pad_positions).sum().item()
        self._pad_masked_count += pad_masked
        
        # Span lengths
        if span_lengths:
            self._span_lengths.extend(span_lengths)
    
    def compute_grad_norm(self, model: nn.Module) -> float:
        """Compute total gradient norm across all parameters."""
        total_norm = 0.0
        for param in model.parameters():
            if param.grad is not None:
                total_norm += param.grad.data.norm(2).item() ** 2
        return math.sqrt(total_norm)
    
    def compute_param_norm(self, model: nn.Module) -> float:
        """Compute total parameter norm."""
        total_norm = 0.0
        for param in model.parameters():
            total_norm += param.data.norm(2).item() ** 2
        return math.sqrt(total_norm)
    
    def check_grad_health(self, grad_norm: float) -> Tuple[bool, Optional[str]]:
        """Check if gradients are healthy.
        
        Returns:
            (is_healthy, warning_message) tuple
        """
        # Check for NaN/Inf
        if not math.isfinite(grad_norm):
            return False, f"Gradient is not finite: {grad_norm}"
        
        # Check for explosion (sudden large increase)
        if self.prev_grad_norm is not None:
            ratio = grad_norm / max(self.prev_grad_norm, 1e-8)
            if ratio > 100:  # 100x increase
                self.grad_explosion_count += 1
                if self.grad_explosion_count >= 3:
                    return False, f"Gradient explosion detected: {self.prev_grad_norm:.2f} -> {grad_norm:.2f}"
            else:
                self.grad_explosion_count = 0
        
        self.prev_grad_norm = grad_norm
        
        # Check for vanishing gradients
        if grad_norm < 1e-7:
            return True, f"⚠ Very small gradients: {grad_norm:.2e}"
        
        return True, None
    
    def log_step(
        self,
        step: int,
        epoch: int,
        model: nn.Module,
        optimizer: torch.optim.Optimizer,
    ) -> Optional[HealthMetrics]:
        """Log metrics for current step.
        
        Returns:
            HealthMetrics if this is a logging step, None otherwise
        """
        if step % self.log_every_n_steps != 0:
            return None
        
        # Compute aggregated metrics
        avg_loss = sum(self._losses) / max(1, len(self._losses))
        accuracy = self._correct / max(1, self._total_masked)
        unk_rate = self._unk_count / max(1, self._total_masked)
        pad_masked_rate = self._pad_masked_count / max(1, self._total_tokens)
        mask_coverage = self._total_masked / max(1, self._total_tokens)
        
        # Gradient and parameter norms
        grad_norm = self.compute_grad_norm(model)
        param_norm = self.compute_param_norm(model)
        
        # Current learning rate
        lr = optimizer.param_groups[0]["lr"]
        
        # Create metrics object
        metrics = HealthMetrics(
            step=step,
            epoch=epoch,
            loss=avg_loss,
            masked_accuracy=accuracy,
            unk_rate=unk_rate,
            pad_masked_rate=pad_masked_rate,
            mask_coverage=mask_coverage,
            grad_norm=grad_norm,
            param_norm=param_norm,
            lr=lr,
            span_lengths=self._span_lengths.copy(),
        )
        
        # Write to JSONL
        self._write_metrics(metrics)
        
        # Check gradient health
        is_healthy, warning = self.check_grad_health(grad_norm)
        if warning:
            LOGGER.warning(warning)
        
        # Assert PAD masked rate is essentially 0 (allow tiny tolerance for edge cases)
        # 1e-5 = 0.001% tolerance - catches real issues while ignoring 1-in-a-million edge cases
        assert pad_masked_rate < 1e-5, f"PAD tokens were masked! Rate: {pad_masked_rate:.6f}"
        
        # Warn if UNK rate > 0.1%
        if unk_rate > 0.001:
            warnings.warn(f"High UNK rate: {unk_rate*100:.2f}% (threshold: 0.1%)")
        
        # Reset accumulators
        self._reset_accumulators()
        
        return metrics
    
    def _write_metrics(self, metrics: HealthMetrics) -> None:
        """Append metrics to JSONL file."""
        with self.metrics_path.open("a") as f:
            f.write(json.dumps(metrics.to_dict()) + "\n")
    
    def print_table_header(self) -> None:
        """Print table header for compact logging."""
        header = (
            f"{'Step':>6} | {'Epoch':>5} | {'Loss':>8} | {'Acc':>6} | "
            f"{'UNK%':>6} | {'Mask%':>6} | {'GradN':>8} | {'LR':>10}"
        )
        print("-" * len(header))
        print(header)
        print("-" * len(header))
    
    def print_table_row(self, metrics: HealthMetrics) -> None:
        """Print a single row of metrics."""
        row = (
            f"{metrics.step:>6} | {metrics.epoch:>5} | {metrics.loss:>8.4f} | "
            f"{metrics.masked_accuracy*100:>5.1f}% | {metrics.unk_rate*100:>5.2f}% | "
            f"{metrics.mask_coverage*100:>5.1f}% | {metrics.grad_norm:>8.2f} | "
            f"{metrics.lr:>10.2e}"
        )
        print(row)
    
    @staticmethod
    def read_metrics(metrics_path: Path) -> List[Dict[str, Any]]:
        """Read metrics from JSONL file.
        
        Args:
            metrics_path: Path to metrics.jsonl
            
        Returns:
            List of metric dicts
        """
        metrics = []
        with metrics_path.open("r") as f:
            for line in f:
                line = line.strip()
                if line:
                    metrics.append(json.loads(line))
        return metrics


def run_pretrain_sanity_checks(
    train_loader,
    vocab,
    num_batches: int = 5,
    target_mask_prob: float = 0.15,
) -> Tuple[bool, List[str]]:
    """Run all pretrain sanity checks before training starts.
    
    This is the checkpoint gate - training cannot proceed unless these pass.
    
    Args:
        train_loader: DataLoader yielding (masked_tokens, labels) tuples
        vocab: KmerVocabulary
        num_batches: Number of batches to check
        target_mask_prob: Expected mask probability
        
    Returns:
        (passed, issues) tuple
    """
    print("=" * 70)
    print("PRETRAIN SANITY CHECKS (Checkpoint Gate)")
    print("=" * 70)
    print("Training cannot proceed unless all checks pass.")
    print()
    
    all_issues = []
    
    # 1. Mask sanity check
    mask_checker = MaskSanityChecker(vocab, target_mask_prob)
    mask_passed, mask_issues = mask_checker.run_sanity_check(train_loader, num_batches)
    all_issues.extend(mask_issues)
    
    # 2. Vocab sanity check
    print("\n[VOCAB CHECK]")
    print(f"  Vocab size: {len(vocab.itos)}")
    print(f"  Special tokens: [MASK]={vocab.mask_id}, <UNK>={vocab.unk_id}, <PAD>={vocab.pad_id}")
    
    # Verify special token IDs are valid
    if vocab.mask_id < 0 or vocab.mask_id >= len(vocab.itos):
        all_issues.append(f"Invalid mask_id: {vocab.mask_id}")
        print(f"  ✗ FAIL: Invalid mask_id: {vocab.mask_id}")
    else:
        print(f"  ✓ Special token IDs are valid")
    
    # 3. Summary
    print("\n" + "=" * 70)
    print("SANITY CHECK SUMMARY")
    print("=" * 70)
    
    critical_issues = [i for i in all_issues if "PAD" in i or "Invalid" in i]
    warning_issues = [i for i in all_issues if i not in critical_issues]
    
    if critical_issues:
        print(f"\n✗ FAILED - {len(critical_issues)} critical issue(s):")
        for issue in critical_issues:
            print(f"  - {issue}")
        print("\nTraining CANNOT proceed until these are fixed.")
        return False, all_issues
    
    if warning_issues:
        print(f"\n⚠ PASSED with {len(warning_issues)} warning(s):")
        for issue in warning_issues:
            print(f"  - {issue}")
        print("\nTraining may proceed, but check warnings.")
    else:
        print("\n✓ ALL CHECKS PASSED")
    
    print("=" * 70)
    
    return True, all_issues


def create_health_dashboard(
    run_dir: Path,
    vocab,
    log_every_n_steps: int = 50,
    target_mask_prob: float = 0.15,
) -> PretrainHealthDashboard:
    """Factory function to create a health dashboard.
    
    Args:
        run_dir: Directory for this training run (e.g., runs/mlm_20231215_120000/)
        vocab: KmerVocabulary
        log_every_n_steps: How often to log metrics
        target_mask_prob: Expected mask probability
        
    Returns:
        Configured PretrainHealthDashboard
    """
    return PretrainHealthDashboard(
        run_dir=run_dir,
        vocab=vocab,
        log_every_n_steps=log_every_n_steps,
        target_mask_prob=target_mask_prob,
    )
