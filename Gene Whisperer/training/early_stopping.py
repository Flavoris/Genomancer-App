"""
Early Stopping for MLM Pretraining

Monitors validation loss and stops training when no improvement is seen
for a specified number of epochs (patience).
"""

from __future__ import annotations

import logging
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
    Early stopping specifically for MLM pretraining.

    Adds additional checks relevant to pretraining:
    - Monitors both validation loss and perplexity
    - Saves MLM-specific metadata with checkpoints
    - Supports k-mer specific configuration
    """

    def __init__(
        self,
        patience: int = 20,
        min_delta: float = 0.001,
        min_epochs: int = 30,
        kmer: int = 6,
        **kwargs,
    ):
        super().__init__(
            patience=patience,
            min_delta=min_delta,
            min_epochs=min_epochs,
            mode="min",  # Always minimize loss for pretraining
            **kwargs,
        )
        self.kmer = kmer
        self.perplexity_history: list[float] = []

    def __call__(
        self,
        val_loss: float,
        model: torch.nn.Module,
        epoch: int,
        checkpoint_dir: Optional[Path] = None,
        val_perplexity: Optional[float] = None,
    ) -> bool:
        """
        Check if pretraining should stop.

        Args:
            val_loss: Validation loss
            model: MLM model
            epoch: Current epoch
            checkpoint_dir: Directory for checkpoints
            val_perplexity: Optional validation perplexity
        """
        # Track perplexity if provided
        if val_perplexity is not None:
            self.perplexity_history.append(val_perplexity)

        return super().__call__(val_loss, model, epoch, checkpoint_dir)

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

        checkpoint_path = checkpoint_dir / f"mlm_encoder_k{self.kmer}_best.pt"

        # Get perplexity if available
        perplexity = self.perplexity_history[-1] if self.perplexity_history else None

        torch.save({
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "best_val_loss": score,
            "best_perplexity": perplexity,
            "kmer": self.kmer,
        }, checkpoint_path)

        self.best_checkpoint_path = checkpoint_path

        if self.verbose:
            LOGGER.info(
                f"Saved k={self.kmer} MLM checkpoint (epoch {epoch}, "
                f"val_loss={score:.6f}, perplexity={perplexity:.2f if perplexity else 'N/A'})"
            )
