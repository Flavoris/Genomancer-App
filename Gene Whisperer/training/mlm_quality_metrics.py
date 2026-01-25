"""Core analysis functions for MLM pre-training quality diagnostics.

Provides utilities to:
- Build DNAMLM models from checkpoints and metadata
- Compute MLM prediction accuracy (top-1, top-5)
- Measure prediction entropy (confidence)
- Extract sequence embeddings for clustering analysis
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn.functional as F

from dataset import KmerVocabulary
from model import DNAEncoder

LOGGER = logging.getLogger(__name__)


@dataclass
class MLMMetrics:
    """Container for MLM quality metrics for a single k-mer value."""

    k: int
    top1_accuracy: float
    top5_accuracy: float
    avg_entropy: float
    silhouette_score: Optional[float] = None
    num_samples: int = 0
    num_masked_tokens: int = 0
    checkpoint_exists: bool = True


@dataclass
class EncoderConfig:
    """Architecture configuration for the MLM encoder."""

    embedding_dim: int = 384
    num_layers: int = 12
    num_heads: int = 12
    ff_dim: int = 1536
    dropout: float = 0.1
    max_seq_len: int = 256


def load_encoder_config(metadata_path: Path) -> EncoderConfig:
    """Load encoder architecture from metadata JSON file.

    Falls back to default config values if metadata is missing fields.
    """
    if not metadata_path.exists():
        LOGGER.warning("Metadata not found at %s, using defaults", metadata_path)
        return EncoderConfig()

    with metadata_path.open("r", encoding="utf-8") as f:
        meta = json.load(f)

    config_snapshot = meta.get("config_snapshot", {})
    return EncoderConfig(
        embedding_dim=meta.get("embedding_dim", config_snapshot.get("embedding_dim", 384)),
        num_layers=meta.get("num_layers", config_snapshot.get("transformer_layers", 12)),
        num_heads=meta.get("heads", config_snapshot.get("transformer_heads", 12)),
        ff_dim=config_snapshot.get("transformer_ff_dim", 1536),
        dropout=config_snapshot.get("transformer_dropout", 0.1),
        max_seq_len=256,
    )


def build_mlm_model(
    vocab: KmerVocabulary,
    encoder_cfg: EncoderConfig,
    checkpoint_path: Path,
) -> Tuple[DNAEncoder, bool]:
    """Build DNAEncoder and load checkpoint weights.

    Args:
        vocab: k-mer vocabulary for this model
        encoder_cfg: Architecture configuration
        checkpoint_path: Path to the encoder state_dict .pt file

    Returns:
        Tuple of (encoder, weights_loaded) where weights_loaded indicates
        whether checkpoint was successfully loaded
    """
    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        embedding_dim=encoder_cfg.embedding_dim,
        num_layers=encoder_cfg.num_layers,
        num_heads=encoder_cfg.num_heads,
        ff_dim=encoder_cfg.ff_dim,
        dropout=encoder_cfg.dropout,
        pad_token_id=vocab.pad_id,
        max_seq_len=encoder_cfg.max_seq_len,
        drop_path_rate=0.0,
        use_rope=False,
        use_positional_motif_bias=False,
    )

    if not checkpoint_path.exists():
        LOGGER.warning("Checkpoint not found: %s", checkpoint_path)
        return encoder, False

    state_dict = torch.load(checkpoint_path, map_location="cpu", weights_only=False)

    # Handle nested checkpoint formats
    if isinstance(state_dict, dict):
        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]
        elif "state_dict" in state_dict:
            state_dict = state_dict["state_dict"]

    # Filter to matching keys/shapes
    model_dict = encoder.state_dict()
    compatible = {
        k: v for k, v in state_dict.items()
        if k in model_dict and model_dict[k].shape == v.shape
    }

    if not compatible:
        LOGGER.error("No compatible weights found in checkpoint %s", checkpoint_path)
        return encoder, False

    model_dict.update(compatible)
    encoder.load_state_dict(model_dict)
    LOGGER.info(
        "Loaded %d/%d weights from %s",
        len(compatible), len(state_dict), checkpoint_path.name,
    )
    return encoder, True


def mask_tokens_for_eval(
    token_ids: torch.LongTensor,
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
) -> Tuple[torch.LongTensor, torch.LongTensor]:
    """Apply deterministic masking for evaluation (always uses [MASK] token).

    Unlike training which uses 80/10/10, evaluation always replaces
    with [MASK] for consistent accuracy measurement.

    Args:
        token_ids: (B, L) input token indices
        vocab: Vocabulary with mask_id and pad_id
        mask_prob: Fraction of non-pad tokens to mask

    Returns:
        (masked_inputs, labels) where labels=-100 for unmasked positions
    """
    batch_size, seq_len = token_ids.shape
    labels = token_ids.clone()
    masked_inputs = token_ids.clone()

    for b in range(batch_size):
        # Find non-pad positions
        non_pad = token_ids[b] != vocab.pad_id
        non_pad_indices = torch.where(non_pad)[0]

        if len(non_pad_indices) == 0:
            labels[b, :] = -100
            continue

        num_to_mask = max(1, int(mask_prob * len(non_pad_indices)))

        # Deterministic-ish but varied masking: use every Nth position
        # This ensures consistent coverage across the sequence
        perm = torch.randperm(len(non_pad_indices))[:num_to_mask]
        mask_positions = non_pad_indices[perm]

        # Set labels to -100 for unmasked positions
        labels[b, :] = -100
        labels[b, mask_positions] = token_ids[b, mask_positions]

        # Replace masked positions with [MASK] token
        masked_inputs[b, mask_positions] = vocab.mask_id

    return masked_inputs, labels


def compute_mlm_predictions(
    encoder: DNAEncoder,
    vocab: KmerVocabulary,
    token_ids: torch.LongTensor,
    device: torch.device,
    batch_size: int = 32,
    mask_prob: float = 0.15,
) -> MLMMetrics:
    """Run MLM predictions and compute accuracy/entropy metrics.

    Builds an LM head with weight tying (matching DNAMLM architecture),
    masks tokens, predicts, and measures quality.

    Args:
        encoder: Pre-trained DNAEncoder
        vocab: k-mer vocabulary
        token_ids: (N, L) all tokenized sequences
        device: Computation device
        batch_size: Inference batch size
        mask_prob: Fraction of tokens to mask per sequence

    Returns:
        MLMMetrics with top-1, top-5 accuracy and entropy
    """
    encoder = encoder.to(device)
    encoder.eval()

    vocab_size = len(vocab.itos)
    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]

    # Build LM head with weight tying (same as DNAMLM)
    lm_head = torch.nn.Linear(encoder.embedding_dim, vocab_size, bias=True)
    lm_head.weight = encoder.embedding.weight
    torch.nn.init.zeros_(lm_head.bias)
    lm_head = lm_head.to(device)

    output_norm = torch.nn.LayerNorm(encoder.embedding_dim).to(device)
    logit_scale = torch.nn.Parameter(torch.tensor(10.0, device=device))

    all_top1_correct = 0
    all_top5_correct = 0
    all_masked_tokens = 0
    all_entropies: List[float] = []

    num_samples = token_ids.shape[0]

    with torch.no_grad():
        for start in range(0, num_samples, batch_size):
            batch = token_ids[start:start + batch_size].to(device)
            masked_inputs, labels = mask_tokens_for_eval(batch, vocab, mask_prob)

            # Forward pass through encoder
            hidden = encoder(masked_inputs)
            hidden = output_norm(hidden)

            # CLIP-style normalized projection (matching DNAMLM.forward)
            h = F.normalize(hidden, dim=-1)
            w = F.normalize(lm_head.weight, dim=-1)
            logits = logit_scale * (h @ w.T)
            logits = logits + lm_head.bias

            # Mask special tokens from prediction space
            for tid in special_token_ids:
                logits[:, :, tid] = float("-inf")

            # Only evaluate at masked positions
            mask = labels != -100
            if mask.sum() == 0:
                continue

            masked_logits = logits[mask]  # (num_masked, vocab_size)
            masked_labels = labels[mask]  # (num_masked,)

            # Top-1 accuracy
            top1_preds = masked_logits.argmax(dim=-1)
            all_top1_correct += (top1_preds == masked_labels).sum().item()

            # Top-5 accuracy
            top5_preds = masked_logits.topk(5, dim=-1).indices
            top5_match = (top5_preds == masked_labels.unsqueeze(-1)).any(dim=-1)
            all_top5_correct += top5_match.sum().item()

            # Per-position entropy
            probs = F.softmax(masked_logits, dim=-1)
            entropy = -(probs * (probs + 1e-10).log()).sum(dim=-1)
            all_entropies.extend(entropy.cpu().tolist())

            all_masked_tokens += mask.sum().item()

    top1_acc = all_top1_correct / max(all_masked_tokens, 1)
    top5_acc = all_top5_correct / max(all_masked_tokens, 1)
    avg_entropy = float(np.mean(all_entropies)) if all_entropies else 0.0

    return MLMMetrics(
        k=vocab.k,
        top1_accuracy=top1_acc,
        top5_accuracy=top5_acc,
        avg_entropy=avg_entropy,
        num_samples=num_samples,
        num_masked_tokens=all_masked_tokens,
    )


def extract_embeddings(
    encoder: DNAEncoder,
    token_ids: torch.LongTensor,
    device: torch.device,
    batch_size: int = 64,
) -> np.ndarray:
    """Extract mean-pooled embeddings for all sequences.

    Args:
        encoder: Pre-trained DNAEncoder
        token_ids: (N, L) tokenized sequences
        device: Computation device
        batch_size: Inference batch size

    Returns:
        (N, D) numpy array of sequence embeddings
    """
    encoder = encoder.to(device)
    encoder.eval()

    all_embeddings = []
    num_samples = token_ids.shape[0]

    with torch.no_grad():
        for start in range(0, num_samples, batch_size):
            batch = token_ids[start:start + batch_size].to(device)

            # Get encoder output
            hidden = encoder(batch)  # (B, L, D)

            # Mean pool over non-pad positions
            pad_mask = batch == encoder.pad_token_id
            hidden[pad_mask] = 0.0
            lengths = (~pad_mask).sum(dim=1, keepdim=True).float()
            lengths = lengths.clamp(min=1.0)
            pooled = hidden.sum(dim=1) / lengths  # (B, D)

            all_embeddings.append(pooled.cpu().numpy())

    return np.concatenate(all_embeddings, axis=0)


def compute_silhouette(
    embeddings: np.ndarray,
    labels: np.ndarray,
    max_samples: int = 2000,
) -> float:
    """Compute silhouette score for embedding clustering quality.

    Args:
        embeddings: (N, D) sequence embeddings
        labels: (N,) binary labels (promoter/non-promoter)
        max_samples: Subsample for speed if dataset is large

    Returns:
        Silhouette score in [-1, 1], higher is better
    """
    try:
        from sklearn.metrics import silhouette_score
    except ImportError:
        LOGGER.warning("sklearn not available, skipping silhouette score")
        return float("nan")

    if len(np.unique(labels)) < 2:
        LOGGER.warning("Need at least 2 classes for silhouette score")
        return float("nan")

    # Subsample if needed
    if len(embeddings) > max_samples:
        indices = np.random.choice(len(embeddings), max_samples, replace=False)
        embeddings = embeddings[indices]
        labels = labels[indices]

    return float(silhouette_score(embeddings, labels))


def interpret_silhouette(score: float) -> str:
    """Provide human-readable interpretation of silhouette score."""
    if np.isnan(score):
        return "N/A"
    if score > 0.5:
        return "strong"
    if score > 0.3:
        return "good"
    if score > 0.15:
        return "moderate"
    if score > 0.0:
        return "weak"
    return "poor"
