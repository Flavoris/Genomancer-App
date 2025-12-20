"""Progressive length curriculum for MLM pretraining.

This module implements curriculum learning for sequence length:
- Start training with shorter sequences (e.g., 32 tokens)
- Gradually ramp up to full length (e.g., 79 tokens) by a configurable % of training
- Helps the model learn local patterns first before global context

References:
- Curriculum Learning (Bengio et al., 2009)
- BERT pretraining benefits from progressive complexity
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple

import torch
import torch.nn.functional as F

LOGGER = logging.getLogger(__name__)


@dataclass
class LengthCurriculumConfig:
    """Configuration for length curriculum schedule.
    
    Attributes:
        enabled: Whether length curriculum is active
        start_tokens: Starting sequence length in tokens
        end_tokens: Final sequence length in tokens (full length)
        start_pct: % of total steps where ramp begins (0.0 = start immediately)
        end_pct: % of total steps where full length is reached
    """
    enabled: bool = True
    start_tokens: int = 32
    end_tokens: int = 79
    start_pct: float = 0.0
    end_pct: float = 0.3
    
    def __post_init__(self):
        """Validate configuration values."""
        if self.start_tokens < 8:
            raise ValueError(f"start_tokens must be >= 8, got {self.start_tokens}")
        if self.end_tokens <= self.start_tokens:
            raise ValueError(
                f"end_tokens ({self.end_tokens}) must be > start_tokens ({self.start_tokens})"
            )
        if not 0.0 <= self.start_pct < 1.0:
            raise ValueError(f"start_pct must be in [0, 1), got {self.start_pct}")
        if not 0.0 < self.end_pct <= 1.0:
            raise ValueError(f"end_pct must be in (0, 1], got {self.end_pct}")
        if self.end_pct <= self.start_pct:
            raise ValueError(
                f"end_pct ({self.end_pct}) must be > start_pct ({self.start_pct})"
            )
    
    @classmethod
    def from_config(cls, cfg: Dict[str, Any]) -> "LengthCurriculumConfig":
        """Create config from YAML/dict configuration."""
        return cls(
            enabled=bool(cfg.get("mlm_length_curriculum_enabled", False)),
            start_tokens=int(cfg.get("mlm_length_curriculum_start_tokens", 32)),
            end_tokens=int(cfg.get("mlm_length_curriculum_end_tokens", 79)),
            start_pct=float(cfg.get("mlm_length_curriculum_start_pct", 0.0)),
            end_pct=float(cfg.get("mlm_length_curriculum_end_pct", 0.3)),
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for logging/serialization."""
        return {
            "enabled": self.enabled,
            "start_tokens": self.start_tokens,
            "end_tokens": self.end_tokens,
            "start_pct": self.start_pct,
            "end_pct": self.end_pct,
        }


class LengthCurriculumScheduler:
    """Scheduler for progressive sequence length during training.
    
    Computes the current sequence length based on training progress.
    Uses linear interpolation between start and end lengths.
    
    Schedule:
    - Steps [0, start_pct * total): Use start_tokens length
    - Steps [start_pct * total, end_pct * total): Linear ramp from start to end
    - Steps [end_pct * total, total): Use end_tokens length (full)
    
    Example:
        >>> scheduler = LengthCurriculumScheduler(
        ...     total_steps=1000,
        ...     config=LengthCurriculumConfig(
        ...         start_tokens=32, end_tokens=79,
        ...         start_pct=0.0, end_pct=0.3
        ...     )
        ... )
        >>> scheduler.get_current_length(0)     # Returns 32
        >>> scheduler.get_current_length(150)   # Returns ~55 (midpoint of ramp)
        >>> scheduler.get_current_length(300)   # Returns 79 (full length)
        >>> scheduler.get_current_length(500)   # Returns 79 (stays at full)
    """
    
    def __init__(
        self,
        total_steps: int,
        config: LengthCurriculumConfig,
    ):
        self.total_steps = max(1, total_steps)
        self.config = config
        
        # Precompute step boundaries
        self.ramp_start_step = int(config.start_pct * total_steps)
        self.ramp_end_step = int(config.end_pct * total_steps)
        
        # Current state
        self._current_step = 0
        self._current_length = config.start_tokens
        
        if config.enabled:
            LOGGER.info(
                "Length curriculum: %d -> %d tokens over steps %d-%d (%.0f%%-%.0f%% of training)",
                config.start_tokens,
                config.end_tokens,
                self.ramp_start_step,
                self.ramp_end_step,
                config.start_pct * 100,
                config.end_pct * 100,
            )
    
    def get_current_length(self, step: Optional[int] = None) -> int:
        """Get the sequence length for a given training step.
        
        Args:
            step: Training step (uses internal state if None)
            
        Returns:
            Current sequence length in tokens
        """
        if not self.config.enabled:
            return self.config.end_tokens
        
        if step is None:
            step = self._current_step
        
        # Before ramp starts
        if step < self.ramp_start_step:
            return self.config.start_tokens
        
        # After ramp ends
        if step >= self.ramp_end_step:
            return self.config.end_tokens
        
        # During ramp: linear interpolation
        ramp_duration = self.ramp_end_step - self.ramp_start_step
        if ramp_duration <= 0:
            return self.config.end_tokens
        
        progress = (step - self.ramp_start_step) / ramp_duration
        length_range = self.config.end_tokens - self.config.start_tokens
        current_length = self.config.start_tokens + int(progress * length_range)
        
        return min(self.config.end_tokens, max(self.config.start_tokens, current_length))
    
    def step(self, step: Optional[int] = None) -> int:
        """Advance the scheduler and return the new length.
        
        Args:
            step: Explicit step to set (increments by 1 if None)
            
        Returns:
            Current sequence length after update
        """
        if step is not None:
            self._current_step = step
        else:
            self._current_step += 1
        
        self._current_length = self.get_current_length(self._current_step)
        return self._current_length
    
    @property
    def current_step(self) -> int:
        """Current training step."""
        return self._current_step
    
    @property
    def current_length(self) -> int:
        """Current sequence length."""
        return self._current_length
    
    def is_ramping(self) -> bool:
        """Returns True if currently in the ramping phase."""
        return (
            self.config.enabled
            and self.ramp_start_step <= self._current_step < self.ramp_end_step
        )
    
    def get_progress(self) -> Dict[str, Any]:
        """Get curriculum progress info for logging."""
        if not self.config.enabled:
            return {
                "enabled": False,
                "current_length": self.config.end_tokens,
            }
        
        phase = "warmup"
        if self._current_step >= self.ramp_end_step:
            phase = "full"
        elif self._current_step >= self.ramp_start_step:
            phase = "ramping"
        
        ramp_progress = 0.0
        if self.is_ramping():
            ramp_duration = max(1, self.ramp_end_step - self.ramp_start_step)
            ramp_progress = (self._current_step - self.ramp_start_step) / ramp_duration
        elif phase == "full":
            ramp_progress = 1.0
        
        return {
            "enabled": True,
            "phase": phase,
            "current_step": self._current_step,
            "current_length": self._current_length,
            "target_length": self.config.end_tokens,
            "ramp_progress": ramp_progress,
        }


class CurriculumMLMCollator:
    """MLM Collator with curriculum-based dynamic sequence length.
    
    Wraps the base masking logic and applies dynamic truncation/padding
    based on the current curriculum length.
    
    Key features:
    - Truncates sequences to current curriculum length
    - Pads shorter sequences with pad_token_id
    - Updates length when scheduler.step() is called
    - Logs length changes at configurable intervals
    
    Note: The scheduler should be stepped externally (in training loop)
    before each batch to get the updated length.
    """
    
    def __init__(
        self,
        vocab: Any,  # KmerVocabulary
        scheduler: LengthCurriculumScheduler,
        mask_prob: float = 0.15,
        track_spans: bool = False,
        debug: bool = False,
    ):
        self.vocab = vocab
        self.scheduler = scheduler
        self.mask_prob = mask_prob
        self.track_spans = track_spans
        self.debug = debug
        self._debug_ran = False
        self._last_logged_length: Optional[int] = None
    
    def __call__(
        self, batch: Iterable[torch.LongTensor]
    ) -> Tuple[torch.LongTensor, torch.LongTensor] | Tuple[torch.LongTensor, torch.LongTensor, List[int]]:
        """Collate batch with curriculum-aware length adjustment.
        
        Steps:
        1. Stack batch into tensor
        2. Truncate/pad to current curriculum length
        3. Apply span masking
        4. Return masked inputs and labels
        """
        from pretrain_mlm import mask_tokens_span
        
        # Stack original batch
        original_inputs = torch.stack(list(batch))
        batch_size, original_length = original_inputs.shape
        
        # Get current curriculum length
        current_length = self.scheduler.current_length
        
        # Log length changes
        if current_length != self._last_logged_length:
            if self._last_logged_length is not None:
                LOGGER.info(
                    "Curriculum length changed: %d -> %d tokens",
                    self._last_logged_length,
                    current_length,
                )
            self._last_logged_length = current_length
        
        # Adjust sequence length
        if current_length < original_length:
            # Truncate to curriculum length
            adjusted_inputs = original_inputs[:, :current_length]
        elif current_length > original_length:
            # Pad to curriculum length (shouldn't happen normally)
            pad_size = current_length - original_length
            adjusted_inputs = F.pad(
                original_inputs,
                (0, pad_size),
                value=self.vocab.pad_id,
            )
        else:
            adjusted_inputs = original_inputs
        
        # Apply span masking
        if self.track_spans:
            masked_inputs, labels, span_lengths = mask_tokens_span(
                adjusted_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=True
            )
        else:
            masked_inputs, labels = mask_tokens_span(
                adjusted_inputs, self.vocab, mask_prob=self.mask_prob, track_spans=False
            )
            span_lengths = []
        
        # Debug verification (runs once if enabled)
        if self.debug and not self._debug_ran:
            self._run_debug_verification(adjusted_inputs, masked_inputs, labels, current_length)
            self._debug_ran = True
        
        if self.track_spans:
            return masked_inputs, labels, span_lengths
        return masked_inputs, labels
    
    def _run_debug_verification(
        self,
        original: torch.LongTensor,
        masked_input: torch.LongTensor,
        labels: torch.LongTensor,
        current_length: int,
    ) -> None:
        """Verify curriculum collation is working correctly."""
        batch_size, seq_len = original.shape
        
        print("=" * 80)
        print("DEBUG MODE: Curriculum MLM Collator Verification")
        print("=" * 80)
        print(f"Batch size: {batch_size}")
        print(f"Current curriculum length: {current_length}")
        print(f"Actual sequence length: {seq_len}")
        print(f"Match: {'YES' if seq_len == current_length else 'NO - MISMATCH!'}")
        
        # Verify shapes match
        assert masked_input.shape == (batch_size, current_length), (
            f"Shape mismatch: expected ({batch_size}, {current_length}), "
            f"got {masked_input.shape}"
        )
        assert labels.shape == (batch_size, current_length), (
            f"Label shape mismatch: expected ({batch_size}, {current_length}), "
            f"got {labels.shape}"
        )
        
        print("Shape verification: PASSED")
        print("=" * 80)


def create_curriculum_collator(
    vocab: Any,
    total_steps: int,
    cfg: Dict[str, Any],
    mask_prob: float = 0.15,
    track_spans: bool = False,
    debug: bool = False,
) -> Tuple[CurriculumMLMCollator, LengthCurriculumScheduler]:
    """Factory function to create curriculum collator and scheduler.
    
    Args:
        vocab: KmerVocabulary for tokenization
        total_steps: Total training steps (for scheduler)
        cfg: Configuration dictionary
        mask_prob: Masking probability for MLM
        track_spans: Whether to track span lengths
        debug: Enable debug verification
        
    Returns:
        Tuple of (collator, scheduler)
    """
    config = LengthCurriculumConfig.from_config(cfg)
    scheduler = LengthCurriculumScheduler(total_steps, config)
    collator = CurriculumMLMCollator(
        vocab=vocab,
        scheduler=scheduler,
        mask_prob=mask_prob,
        track_spans=track_spans,
        debug=debug,
    )
    return collator, scheduler


def test_length_curriculum():
    """Basic unit test for length curriculum functionality."""
    print("\n" + "=" * 60)
    print("Testing LengthCurriculumScheduler")
    print("=" * 60)
    
    config = LengthCurriculumConfig(
        enabled=True,
        start_tokens=32,
        end_tokens=79,
        start_pct=0.0,
        end_pct=0.3,
    )
    
    total_steps = 1000
    scheduler = LengthCurriculumScheduler(total_steps, config)
    
    # Test key points
    test_cases = [
        (0, 32, "Start of training"),
        (150, None, "Midpoint of ramp"),  # Should be ~55
        (299, None, "Near end of ramp"),  # Should be ~78
        (300, 79, "End of ramp"),
        (500, 79, "After ramp"),
        (1000, 79, "End of training"),
    ]
    
    all_passed = True
    for step, expected, description in test_cases:
        actual = scheduler.get_current_length(step)
        if expected is not None:
            passed = actual == expected
            status = "PASS" if passed else "FAIL"
            print(f"  Step {step:4d}: {actual:3d} tokens (expected {expected:3d}) - {status} ({description})")
            if not passed:
                all_passed = False
        else:
            # Just verify it's in range
            in_range = config.start_tokens <= actual <= config.end_tokens
            status = "PASS" if in_range else "FAIL"
            print(f"  Step {step:4d}: {actual:3d} tokens (in range) - {status} ({description})")
            if not in_range:
                all_passed = False
    
    # Test schedule is monotonically increasing
    print("\nMonotonicity check (length should only increase):")
    prev_length = 0
    for step in range(0, total_steps + 1, 50):
        length = scheduler.get_current_length(step)
        if length < prev_length:
            print(f"  FAIL: Step {step} has length {length} < previous {prev_length}")
            all_passed = False
        prev_length = length
    
    if all_passed:
        print("  All monotonicity checks PASSED")
    
    # Test disabled curriculum
    print("\nDisabled curriculum check:")
    disabled_config = LengthCurriculumConfig(enabled=False, start_tokens=32, end_tokens=79)
    disabled_scheduler = LengthCurriculumScheduler(1000, disabled_config)
    for step in [0, 500, 1000]:
        length = disabled_scheduler.get_current_length(step)
        passed = length == 79
        print(f"  Step {step}: {length} tokens (expected 79) - {'PASS' if passed else 'FAIL'}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 60)
    print(f"Overall: {'ALL TESTS PASSED' if all_passed else 'SOME TESTS FAILED'}")
    print("=" * 60 + "\n")
    
    return all_passed


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    test_length_curriculum()

