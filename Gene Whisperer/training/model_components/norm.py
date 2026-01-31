"""Normalization layers for Gene Whisperer models."""
import torch
import torch.nn as nn


class RMSNorm(nn.Module):
    """Root Mean Square Layer Normalization.

    Simpler and faster than LayerNorm:
    - No mean centering
    - No bias term
    Used in LLaMA, Mistral, and other modern LLMs.
    """

    def __init__(self, dim: int, eps: float = 1e-6):
        super().__init__()
        self.eps = eps
        self.weight = nn.Parameter(torch.ones(dim))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # RMS normalization
        rms = torch.sqrt(x.pow(2).mean(dim=-1, keepdim=True) + self.eps)
        x_normed = x / rms
        return x_normed * self.weight
