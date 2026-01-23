"""Forward hooks for ablation experiments.

Provides hooks to replace transformer output with random vectors,
isolating the contribution of engineered features.
"""
from __future__ import annotations

import torch
import torch.nn as nn


class TransformerOutputRandomizer:
    """Replaces pooled transformer output with random vectors.

    Supports both V5 (GeneWhispererStage1 with attention pooling)
    and V6 (GeneWhispererV6 with mean pooling via final_norm).
    """

    def __init__(self, model: nn.Module):
        self.model = model
        self._hook_target = self._find_hook_target()

    def attach(self) -> torch.utils.hooks.RemovableHook:
        """Register the forward hook and return a removable handle."""
        return self._hook_target.register_forward_hook(self._hook_fn)

    def _find_hook_target(self) -> nn.Module:
        """Identify the correct module to hook based on model architecture."""
        # V5: GeneWhispererStage1 has self.pool (AttentionPooling)
        if hasattr(self.model, "pool") and hasattr(self.model, "use_attention_pool"):
            if self.model.use_attention_pool:
                return self.model.pool

        # V6 and fallback: hook on final_norm (output feeds into mean pooling)
        if hasattr(self.model, "final_norm"):
            return self.model.final_norm

        raise ValueError(
            f"Cannot determine hook target for model type: {type(self.model).__name__}"
        )

    def _hook_fn(
        self,
        module: nn.Module,
        input_tensor: tuple,
        output: torch.Tensor,
    ) -> torch.Tensor:
        """Replace module output with random vectors of the same shape."""
        # Determine embedding dimension from the last dimension
        embed_dim = output.shape[-1]

        if output.dim() == 3:
            # Shape (B, L, D): final_norm output for V6
            # Use a single random vector per sample, repeated across L
            # so mean pooling yields that same random vector
            batch_size = output.shape[0]
            seq_len = output.shape[1]
            random_vec = torch.randn(
                batch_size, 1, embed_dim,
                device=output.device,
                dtype=output.dtype,
            )
            return random_vec.expand(batch_size, seq_len, embed_dim)

        # Shape (B, D): attention pooling output for V5
        return torch.randn_like(output)
