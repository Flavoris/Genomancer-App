"""
Early Stopping for MLM Pretraining

Monitors validation loss and accuracy with intelligent convergence detection:
- Rate-of-change detection (stops when improvement rate drops)
- Accuracy stagnation detection
- "Good enough" thresholds for early exit
- K-mer aware configuration (smaller k-mers need less training)
"""

from __future__ import annotations

import logging
from collections import deque
from pathlib import Path
from typing import Optional, Literal
import torch

LOGGER = logging.getLogger(__name__)


class EarlyStopping:
    """
    Early stopping to halt training when validation metric stops improving.

    Supports:
    - Patience-based stopping (wait N epochs for improvement)
    - Minimum delta (ignore tiny improvements)
    - Minimum epochs (don't stop before M epochs)
    - Best model checkpoint saving/restoring

    Example:
        early_stopper = EarlyStopping(patience=20, min_delta=0.001)

        for epoch in range(max_epochs):
            train_loss = train_one_epoch(...)
            val_loss = validate(...)

            early_stopper(val_loss, model, epoch, checkpoint_path)

            if early_stopper.should_stop:
                print(f"Early stopping at epoch {epoch}")
                break

        # Restore best model
        early_stopper.restore_best_weights(model)
    """

    def __init__(
        self,
        patience: int = 20,
        min_delta: float = 0.001,
        min_epochs: int = 30,
        mode: Literal["min", "max"] = "min",
        restore_best: bool = True,
        verbose: bool = True,
    ):
        """
        Args:
            patience: Number of epochs to wait for improvement before stopping
            min_delta: Minimum change in monitored metric to qualify as improvement
            min_epochs: Minimum epochs before early stopping can trigger
            mode: "min" for loss (lower is better), "max" for accuracy (higher is better)
            restore_best: Whether to restore best weights when training ends
            verbose: Whether to log early stopping events
        """
        self.patience = patience
        self.min_delta = min_delta
        self.min_epochs = min_epochs
        self.mode = mode
        self.restore_best = restore_best
        self.verbose = verbose

        # State
        self.best_score: Optional[float] = None
        self.best_epoch: int = 0
        self.epochs_without_improvement: int = 0
        self.should_stop: bool = False
        self.best_checkpoint_path: Optional[Path] = None

        # For tracking
        self.history: list[dict] = []

    def __call__(
        self,
        current_score: float,
        model: torch.nn.Module,
        epoch: int,
        checkpoint_dir: Optional[Path] = None,
    ) -> bool:
        """
        Check if training should stop and optionally save checkpoint.

        Args:
            current_score: Current validation metric value
            model: Model to potentially checkpoint
            epoch: Current epoch number
            checkpoint_dir: Directory to save best checkpoint

        Returns:
            True if training should stop, False otherwise
        """
        # Track history
        self.history.append({
            "epoch": epoch,
            "score": current_score,
            "best_score": self.best_score,
        })

        # Check if this is an improvement
        is_improvement = self._is_improvement(current_score)

        if is_improvement:
            if self.verbose and self.best_score is not None:
                improvement = abs(current_score - self.best_score)
                LOGGER.info(
                    f"Epoch {epoch}: Validation improved by {improvement:.6f} "
                    f"({self.best_score:.6f} -> {current_score:.6f})"
                )

            self.best_score = current_score
            self.best_epoch = epoch
            self.epochs_without_improvement = 0

            # Save checkpoint
            if checkpoint_dir is not None:
                self._save_checkpoint(model, checkpoint_dir, epoch, current_score)
        else:
            self.epochs_without_improvement += 1

            if self.verbose:
                LOGGER.info(
                    f"Epoch {epoch}: No improvement for {self.epochs_without_improvement} epochs "
                    f"(best: {self.best_score:.6f} at epoch {self.best_epoch})"
                )

        # Check stopping conditions (only stop AFTER min_epochs have passed)
        if epoch > self.min_epochs and self.epochs_without_improvement >= self.patience:
            self.should_stop = True
            if self.verbose:
                LOGGER.info(
                    f"Early stopping triggered at epoch {epoch}. "
                    f"Best score: {self.best_score:.6f} at epoch {self.best_epoch}"
                )

        return self.should_stop

    def _is_improvement(self, current_score: float) -> bool:
        """Check if current score is an improvement over best."""
        if self.best_score is None:
            return True

        if self.mode == "min":
            # For loss: lower is better, improvement = best - current > delta
            return (self.best_score - current_score) > self.min_delta
        else:
            # For accuracy: higher is better, improvement = current - best > delta
            return (current_score - self.best_score) > self.min_delta

    def _save_checkpoint(
        self,
        model: torch.nn.Module,
        checkpoint_dir: Path,
        epoch: int,
        score: float,
    ) -> None:
        """Save model checkpoint."""
        checkpoint_dir = Path(checkpoint_dir)
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        checkpoint_path = checkpoint_dir / "best_model.pt"

        torch.save({
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "best_score": score,
        }, checkpoint_path)

        self.best_checkpoint_path = checkpoint_path

        if self.verbose:
            LOGGER.debug(f"Saved best checkpoint to {checkpoint_path}")

    def restore_best_weights(self, model: torch.nn.Module) -> bool:
        """
        Restore model weights from best checkpoint.

        Returns:
            True if weights were restored, False otherwise
        """
        if not self.restore_best:
            return False

        if self.best_checkpoint_path is None or not self.best_checkpoint_path.exists():
            LOGGER.warning("No best checkpoint found to restore")
            return False

        checkpoint = torch.load(self.best_checkpoint_path, map_location="cpu")
        model.load_state_dict(checkpoint["model_state_dict"])

        if self.verbose:
            LOGGER.info(
                f"Restored best weights from epoch {checkpoint['epoch']} "
                f"(score: {checkpoint['best_score']:.6f})"
            )

        return True

    def get_summary(self) -> dict:
        """Get summary of early stopping state."""
        return {
            "best_score": self.best_score,
            "best_epoch": self.best_epoch,
            "epochs_without_improvement": self.epochs_without_improvement,
            "stopped_early": self.should_stop,
            "total_epochs_tracked": len(self.history),
        }


class PretrainingEarlyStopping(EarlyStopping):
    """
    Early stopping specifically for MLM pretraining with intelligent convergence detection.

    Features beyond basic EarlyStopping:
    - Monitors both validation loss AND accuracy
    - Rate-of-change detection (stops when improvement rate drops below threshold)
    - Accuracy stagnation detection over a sliding window
    - "Good enough" thresholds for early exit when performance is satisfactory
    - Tokenizer-aware configuration (BPE and k-mer have different convergence patterns)
    - Saves MLM-specific metadata with checkpoints
    """

    # Default tokenizer-specific settings (BPE and smaller k-mers converge faster)
    TOKENIZER_DEFAULTS = {
        "bpe": {"min_epochs": 10, "patience": 15, "target_accuracy": 0.75},
        3: {"min_epochs": 5, "patience": 8, "target_accuracy": 0.85},
        4: {"min_epochs": 8, "patience": 10, "target_accuracy": 0.82},
        5: {"min_epochs": 10, "patience": 12, "target_accuracy": 0.80},
        6: {"min_epochs": 12, "patience": 15, "target_accuracy": 0.78},
    }
    # Backwards compatibility alias
    KMER_DEFAULTS = TOKENIZER_DEFAULTS

    def __init__(
        self,
        patience: int = 20,
        min_delta: float = 0.001,
        min_epochs: int = 30,
        kmer: int = 6,
        tokenizer_type: Optional[str] = None,
        # New convergence detection parameters
        accuracy_stagnation_window: int = 5,
        accuracy_stagnation_threshold: float = 0.01,
        rate_of_change_window: int = 10,
        rate_of_change_min_improvement: float = 0.005,
        target_accuracy: Optional[float] = None,
        target_loss: float = 0.5,
        use_kmer_defaults: bool = True,
        **kwargs,
    ):
        """
        Args:
            patience: Epochs to wait for improvement before stopping
            min_delta: Minimum change in loss to qualify as improvement
            min_epochs: Minimum epochs before stopping can trigger
            kmer: K-mer size (used for k-mer aware defaults, ignored if tokenizer_type="bpe")
            tokenizer_type: "bpe" or "kmer" - determines defaults and checkpoint naming
            accuracy_stagnation_window: Epochs to check for accuracy stagnation
            accuracy_stagnation_threshold: Min accuracy improvement over window
            rate_of_change_window: Epochs to calculate improvement rate
            rate_of_change_min_improvement: Min improvement rate per epoch
            target_accuracy: Stop if accuracy exceeds this (good enough)
            target_loss: Stop if loss drops below this (good enough)
            use_kmer_defaults: Whether to use tokenizer-specific defaults
        """
        # Determine tokenizer type: explicit param > kmer=1 means BPE > default to kmer
        if tokenizer_type is not None:
            self.tokenizer_type = tokenizer_type
        elif kmer == 1:
            self.tokenizer_type = "bpe"
        else:
            self.tokenizer_type = "kmer"

        self.kmer = kmer if self.tokenizer_type == "kmer" else None

        # Apply tokenizer-specific defaults if enabled
        defaults_key = "bpe" if self.tokenizer_type == "bpe" else kmer
        if use_kmer_defaults and defaults_key in self.TOKENIZER_DEFAULTS:
            defaults = self.TOKENIZER_DEFAULTS[defaults_key]
            min_epochs = defaults["min_epochs"]
            patience = defaults["patience"]
            if target_accuracy is None:
                target_accuracy = defaults["target_accuracy"]
            if self.tokenizer_type == "bpe":
                LOGGER.info(
                    f"Using BPE defaults: min_epochs={min_epochs}, "
                    f"patience={patience}, target_accuracy={target_accuracy}"
                )
            else:
                LOGGER.info(
                    f"Using k={kmer} defaults: min_epochs={min_epochs}, "
                    f"patience={patience}, target_accuracy={target_accuracy}"
                )

        super().__init__(
            patience=patience,
            min_delta=min_delta,
            min_epochs=min_epochs,
            mode="min",  # Always minimize loss for pretraining
            **kwargs,
        )

        # Convergence detection settings
        self.accuracy_stagnation_window = accuracy_stagnation_window
        self.accuracy_stagnation_threshold = accuracy_stagnation_threshold
        self.rate_of_change_window = rate_of_change_window
        self.rate_of_change_min_improvement = rate_of_change_min_improvement
        self.target_accuracy = target_accuracy
        self.target_loss = target_loss

        # History tracking for convergence detection
        self.accuracy_history: deque = deque(maxlen=max(
            accuracy_stagnation_window, rate_of_change_window
        ) + 1)
        self.loss_history: deque = deque(maxlen=rate_of_change_window + 1)
        self.perplexity_history: list[float] = []

        # Stop reasons for logging
        self.stop_reason: Optional[str] = None

    def __call__(
        self,
        val_loss: float,
        model: torch.nn.Module,
        epoch: int,
        checkpoint_dir: Optional[Path] = None,
        val_perplexity: Optional[float] = None,
        val_accuracy: Optional[float] = None,
    ) -> bool:
        """
        Check if pretraining should stop with intelligent convergence detection.

        Args:
            val_loss: Validation loss
            model: MLM model
            epoch: Current epoch
            checkpoint_dir: Directory for checkpoints
            val_perplexity: Optional validation perplexity
            val_accuracy: Validation accuracy (required for smart stopping)
        """
        # If already stopped in a previous epoch, keep the original stop reason.
        if self.should_stop:
            return True

        # Track metrics history
        if val_perplexity is not None:
            self.perplexity_history.append(val_perplexity)

        if val_accuracy is not None:
            self.accuracy_history.append(val_accuracy)

        self.loss_history.append(val_loss)

        # Run base early stopping check (loss-based)
        base_should_stop = super().__call__(val_loss, model, epoch, checkpoint_dir)

        if base_should_stop:
            self.stop_reason = "patience_exceeded"
            return True

        # Skip additional checks until min_epochs passed
        if epoch < self.min_epochs:
            return False

        # Check "good enough" thresholds
        if self._check_good_enough(val_loss, val_accuracy, epoch):
            self.should_stop = True
            return True

        # Check accuracy stagnation
        if val_accuracy is not None and self._check_accuracy_stagnation(epoch):
            self.should_stop = True
            return True

        # Check rate of change (diminishing returns)
        if self._check_diminishing_returns(epoch):
            self.should_stop = True
            return True

        return False

    def _check_good_enough(
        self, val_loss: float, val_accuracy: Optional[float], epoch: int
    ) -> bool:
        """Check if we've reached 'good enough' performance to stop early."""
        # Check accuracy threshold
        if val_accuracy is not None and self.target_accuracy is not None:
            if val_accuracy >= self.target_accuracy:
                self.stop_reason = f"target_accuracy_reached ({val_accuracy:.1%} >= {self.target_accuracy:.1%})"
                LOGGER.info(
                    f"Epoch {epoch}: Target accuracy reached! "
                    f"{val_accuracy:.1%} >= {self.target_accuracy:.1%} - stopping early"
                )
                return True

        # Check loss threshold
        if val_loss <= self.target_loss:
            self.stop_reason = f"target_loss_reached ({val_loss:.4f} <= {self.target_loss:.4f})"
            LOGGER.info(
                f"Epoch {epoch}: Target loss reached! "
                f"{val_loss:.4f} <= {self.target_loss:.4f} - stopping early"
            )
            return True

        return False

    def _check_accuracy_stagnation(self, epoch: int) -> bool:
        """Check if accuracy has stagnated over the window."""
        if len(self.accuracy_history) < self.accuracy_stagnation_window:
            return False

        # Get accuracy values from window
        recent_accuracies = list(self.accuracy_history)[-self.accuracy_stagnation_window:]
        oldest = recent_accuracies[0]
        newest = recent_accuracies[-1]
        improvement = newest - oldest

        # Also check max improvement in window
        min_in_window = min(recent_accuracies)
        max_in_window = max(recent_accuracies)
        range_in_window = max_in_window - min_in_window

        if not (
            improvement < self.accuracy_stagnation_threshold
            and range_in_window < self.accuracy_stagnation_threshold * 2
        ):
            return False

        # Guardrail: don't stop on accuracy stagnation if loss is still improving meaningfully
        if len(self.loss_history) >= self.accuracy_stagnation_window:
            recent_losses = list(self.loss_history)[-self.accuracy_stagnation_window:]
            loss_improvement = recent_losses[0] - recent_losses[-1]
            min_loss_improvement = self.min_delta * self.accuracy_stagnation_window
            if loss_improvement >= min_loss_improvement:
                LOGGER.info(
                    "Epoch %d: Accuracy stagnation detected but loss improved by %.6f "
                    "(>= %.6f); continuing training.",
                    epoch,
                    loss_improvement,
                    min_loss_improvement,
                )
                return False

        self.stop_reason = (
            f"accuracy_stagnation (improvement={improvement:.4f} < "
            f"{self.accuracy_stagnation_threshold:.4f} over {self.accuracy_stagnation_window} epochs)"
        )
        LOGGER.info(
            f"Epoch {epoch}: Accuracy stagnation detected! "
            f"Improvement over last {self.accuracy_stagnation_window} epochs: "
            f"{improvement:.4f} (threshold: {self.accuracy_stagnation_threshold:.4f})"
        )
        return True

    def _check_diminishing_returns(self, epoch: int) -> bool:
        """Check if improvement rate has dropped below threshold (diminishing returns)."""
        if len(self.loss_history) < self.rate_of_change_window:
            return False

        # Calculate average improvement rate over window
        recent_losses = list(self.loss_history)[-self.rate_of_change_window:]
        first_half = recent_losses[:len(recent_losses)//2]
        second_half = recent_losses[len(recent_losses)//2:]

        avg_first = sum(first_half) / len(first_half)
        avg_second = sum(second_half) / len(second_half)

        # For loss, improvement is positive when avg_first > avg_second
        improvement_rate = (avg_first - avg_second) / (self.rate_of_change_window / 2)

        # Also check if we're oscillating (no net progress)
        net_improvement = recent_losses[0] - recent_losses[-1]

        if improvement_rate < self.rate_of_change_min_improvement and \
           net_improvement < self.rate_of_change_min_improvement * self.rate_of_change_window:
            self.stop_reason = (
                f"diminishing_returns (rate={improvement_rate:.6f}/epoch < "
                f"{self.rate_of_change_min_improvement:.6f})"
            )
            LOGGER.info(
                f"Epoch {epoch}: Diminishing returns detected! "
                f"Improvement rate: {improvement_rate:.6f}/epoch "
                f"(threshold: {self.rate_of_change_min_improvement:.6f})"
            )
            return True

        return False

    def _save_checkpoint(
        self,
        model: torch.nn.Module,
        checkpoint_dir: Path,
        epoch: int,
        score: float,
    ) -> None:
        """Save MLM checkpoint with metadata."""
        checkpoint_dir = Path(checkpoint_dir)
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        # Use tokenizer type for checkpoint naming
        if self.tokenizer_type == "bpe":
            checkpoint_path = checkpoint_dir / "mlm_encoder_bpe_best.pt"
            tokenizer_label = "BPE"
        else:
            checkpoint_path = checkpoint_dir / f"mlm_encoder_k{self.kmer}_best.pt"
            tokenizer_label = f"k={self.kmer}"

        # Get latest perplexity and accuracy if available
        perplexity = self.perplexity_history[-1] if self.perplexity_history else None
        accuracy = self.accuracy_history[-1] if self.accuracy_history else None

        torch.save({
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "best_score": score,  # Required for base class restore_best_weights
            "best_val_loss": score,
            "best_perplexity": perplexity,
            "best_accuracy": accuracy,
            "kmer": self.kmer,
            "tokenizer_type": self.tokenizer_type,
        }, checkpoint_path)

        self.best_checkpoint_path = checkpoint_path

        if self.verbose:
            perplexity_str = f"{perplexity:.2f}" if perplexity is not None else "N/A"
            accuracy_str = f"{accuracy:.1%}" if accuracy is not None else "N/A"
            LOGGER.info(
                f"Saved {tokenizer_label} MLM checkpoint (epoch {epoch}, "
                f"val_loss={score:.6f}, perplexity={perplexity_str}, accuracy={accuracy_str})"
            )

    def get_summary(self) -> dict:
        """Get summary of early stopping state including convergence info."""
        base_summary = super().get_summary()
        base_summary.update({
            "kmer": self.kmer,
            "tokenizer_type": self.tokenizer_type,
            "stop_reason": self.stop_reason,
            "final_accuracy": self.accuracy_history[-1] if self.accuracy_history else None,
            "accuracy_history_len": len(self.accuracy_history),
            "target_accuracy": self.target_accuracy,
            "target_loss": self.target_loss,
        })
        return base_summary
