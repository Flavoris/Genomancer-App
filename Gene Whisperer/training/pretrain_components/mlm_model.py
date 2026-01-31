"""DNAMLM: Masked Language Model wrapper for DNA sequences."""
from __future__ import annotations

import logging
import math
from typing import List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

# Import encoder from parent training directory
import sys
from pathlib import Path
_training_dir = Path(__file__).resolve().parent.parent
if str(_training_dir) not in sys.path:
    sys.path.insert(0, str(_training_dir))

from model import DNAEncoder

LOGGER = logging.getLogger(__name__)


class DNAMLM(nn.Module):
    """Masked Language Model wrapper around DNAEncoder.

    Key improvements for lower loss:
    1. Weight tying between embedding and LM head (BERT-style)
    2. Bias in LM head for better output distribution modeling
    3. Optional output layer norm before projection
    4. Special token exclusion: [MASK], <UNK>, <PAD> are excluded from
       prediction space by setting their logits to -inf.
    """

    def __init__(
        self,
        encoder: DNAEncoder,
        vocab_size: int,
        special_token_ids: List[int] | None = None,
        tie_weights: bool = True,
        use_output_norm: bool = True,
        use_clip_projection: bool = False,
    ):
        super().__init__()
        self.encoder = encoder
        self.vocab_size = vocab_size
        self.tie_weights = tie_weights
        self.use_clip_projection = use_clip_projection

        self.special_token_ids = special_token_ids or []
        if self.special_token_ids:
            LOGGER.info(
                "Special token exclusion enabled: IDs %s will be masked from predictions",
                self.special_token_ids
            )

        self.use_output_norm = use_output_norm
        if use_output_norm:
            self.output_norm = nn.LayerNorm(encoder.embedding_dim)

        self.lm_head = nn.Linear(encoder.embedding_dim, vocab_size, bias=True)

        if use_clip_projection:
            self.logit_scale = nn.Parameter(torch.tensor(10.0))
            LOGGER.warning(
                "CLIP-style projection enabled. This is NOT recommended for MLM - "
                "consider setting mlm_use_clip_projection: false in config."
            )
        else:
            self.logit_scale = None
            LOGGER.info("Using standard MLM projection (recommended)")

        if tie_weights:
            self.lm_head.weight = encoder.embedding.weight
            LOGGER.info("Weight tying enabled: LM head shares embedding weights")
        else:
            nn.init.trunc_normal_(self.lm_head.weight, std=0.02, a=-0.04, b=0.04)

        nn.init.zeros_(self.lm_head.bias)

    def _mask_special_tokens(self, logits: torch.Tensor) -> torch.Tensor:
        """Set logits for special tokens to -inf so they can never be predicted."""
        if not self.special_token_ids:
            return logits

        logits = logits.clone()

        for token_id in self.special_token_ids:
            logits[:, :, token_id] = float('-inf')

        return logits

    def forward(
        self,
        token_ids: torch.LongTensor,
        labels: torch.LongTensor | None = None,
        key_padding_mask: torch.Tensor | None = None,
    ):
        """Forward pass for MLM.

        Args:
            token_ids: (B, L) input token indices
            labels: (B, L) target labels, -100 where we ignore
            key_padding_mask: Optional (B, L) boolean mask for padding

        Returns:
            (logits, loss) tuple where loss is None if labels not provided
        """
        x = self.encoder(token_ids, key_padding_mask=key_padding_mask)

        if self.use_output_norm:
            x = self.output_norm(x)

        if self.use_clip_projection:
            h = F.normalize(x, dim=-1)
            w = F.normalize(self.lm_head.weight, dim=-1)
            logits = self.logit_scale * (h @ w.T)
            if self.lm_head.bias is not None:
                logits = logits + self.lm_head.bias
        else:
            logits = self.lm_head(x)

        logits = self._mask_special_tokens(logits)

        if labels is None:
            return logits, None

        logits = logits.float()
        assert logits.dtype == torch.float32
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
        )
        return logits, loss

    def compute_accuracy(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute masked token prediction accuracy."""
        preds = logits.argmax(dim=-1)
        mask = labels != -100
        if mask.sum() == 0:
            return 0.0
        correct = (preds[mask] == labels[mask]).sum().item()
        total = mask.sum().item()
        return correct / total

    def compute_perplexity(
        self,
        logits: torch.Tensor,
        labels: torch.LongTensor,
    ) -> float:
        """Compute perplexity on masked tokens."""
        logits = logits.float()
        assert logits.dtype == torch.float32
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
            reduction='mean',
        )
        return math.exp(loss.item())
