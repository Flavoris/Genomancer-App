"""Positional Motif Attention Bias for DNA Promoter Detection.

Promoters have known regulatory motifs at specific positions relative to the TSS:
- TATA box: typically at -30 to -25 bp (positions ~45-55 in 81bp window centered on TSS)
- -10 box: typically at -12 to -7 bp (positions ~64-74)
- TSS: position 0 (positions ~75-81)

This module adds learnable position-aware attention bias that:
1. Initializes with slight bias toward known regulatory regions
2. Lets the model learn optimal attention patterns during training
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn


# Default promoter priors for 81bp sequences centered on TSS
# Format: {(start, end): bias_value}
# These values are small enough to be just a hint, not a strong prior
DEFAULT_PROMOTER_PRIORS: Dict[Tuple[int, int], float] = {
    (45, 55): 0.1,   # TATA box region (-35 to -25 bp)
    (64, 74): 0.1,   # -10 box region (-15 to -5 bp)
    (75, 81): 0.15,  # TSS region (most important)
}


class PositionalMotifBias(nn.Module):
    """Learnable attention bias for known promoter motif positions.

    Initializes with slight bias toward known regulatory regions,
    then lets the model learn optimal attention patterns.

    The bias is added to attention scores before softmax, allowing
    the model to naturally attend more to biologically important
    positions while still being fully learnable.

    Args:
        num_heads: Number of attention heads (each learns separate bias)
        max_seq_len: Maximum sequence length
        init_bias_positions: Optional dict mapping (start, end) -> bias_value
            for initial bias toward known motif positions
        learnable: If True, bias is learnable. If False, bias is fixed.
    """

    def __init__(
        self,
        num_heads: int,
        max_seq_len: int,
        init_bias_positions: Optional[Dict[Tuple[int, int], float]] = None,
        learnable: bool = True,
    ):
        super().__init__()
        self.num_heads = num_heads
        self.max_seq_len = max_seq_len

        # Learnable position-position bias (like ALiBi but learnable)
        # Shape: (H, L, L) - each head has its own bias matrix
        if learnable:
            self.bias = nn.Parameter(
                torch.zeros(num_heads, max_seq_len, max_seq_len)
            )
        else:
            self.register_buffer(
                "bias",
                torch.zeros(num_heads, max_seq_len, max_seq_len)
            )

        # Initialize with biological priors if provided
        if init_bias_positions:
            self._init_with_priors(init_bias_positions)

    def _init_with_priors(self, priors: Dict[Tuple[int, int], float]) -> None:
        """Initialize with known motif positions.

        For 81bp sequences centered on TSS:
        - TATA box: typically at -30 to -25 (positions ~45-50 in 81bp window)
        - -10 box: typically at -12 to -7 (positions ~64-69)
        - TSS: position 0 (position ~75-81)

        The bias is symmetric in the sense that we boost attention both
        TO these important positions (from all other positions) and
        FROM these positions (to all other positions), though the latter
        is weighted less heavily.
        """
        with torch.no_grad():
            for (start, end), bias_val in priors.items():
                # Clamp to valid range
                start = max(0, min(start, self.max_seq_len))
                end = max(0, min(end, self.max_seq_len))

                if start >= end:
                    continue

                # Boost attention TO these positions (from all positions)
                # This is the primary effect: all queries should attend to motifs
                self.bias[:, :, start:end] += bias_val

                # Also boost attention FROM these positions (to all positions)
                # This is secondary: motif positions should have broader context
                self.bias[:, start:end, :] += bias_val * 0.5

    def forward(self, attn_scores: torch.Tensor) -> torch.Tensor:
        """Add positional motif bias to attention scores before softmax.

        Args:
            attn_scores: Attention scores of shape (B, H, L, L)
                where B=batch, H=heads, L=sequence length

        Returns:
            Modified attention scores with positional bias added
        """
        B, H, L, _ = attn_scores.shape

        # Handle variable sequence lengths by slicing the bias
        if L > self.max_seq_len:
            raise ValueError(
                f"Sequence length {L} exceeds max_seq_len {self.max_seq_len}. "
                f"Increase max_seq_len in PositionalMotifBias."
            )

        # Slice bias to match actual sequence length
        bias_slice = self.bias[:, :L, :L]

        # Add bias (broadcast over batch dimension)
        return attn_scores + bias_slice.unsqueeze(0)

    def extra_repr(self) -> str:
        return f"num_heads={self.num_heads}, max_seq_len={self.max_seq_len}"


def create_promoter_motif_bias(
    num_heads: int,
    max_seq_len: int = 81,
    use_default_priors: bool = True,
    custom_priors: Optional[Dict[Tuple[int, int], float]] = None,
    learnable: bool = True,
) -> PositionalMotifBias:
    """Factory function to create PositionalMotifBias with promoter-specific priors.

    Args:
        num_heads: Number of attention heads
        max_seq_len: Maximum sequence length (default 81 for promoter detection)
        use_default_priors: If True, use DEFAULT_PROMOTER_PRIORS for initialization
        custom_priors: Optional custom priors to use instead of or in addition to defaults
        learnable: If True, bias parameters are learnable

    Returns:
        Configured PositionalMotifBias module
    """
    priors = {}

    if use_default_priors:
        priors.update(DEFAULT_PROMOTER_PRIORS)

    if custom_priors:
        priors.update(custom_priors)

    return PositionalMotifBias(
        num_heads=num_heads,
        max_seq_len=max_seq_len,
        init_bias_positions=priors if priors else None,
        learnable=learnable,
    )
