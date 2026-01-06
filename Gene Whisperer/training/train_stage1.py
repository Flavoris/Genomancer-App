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
import time
from collections import defaultdict
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
from torch.profiler import profile, ProfilerActivity, schedule, tensorboard_trace_handler
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler
from tqdm.auto import tqdm

from amp_utils import get_amp_context, log_amp_status, AMPContext
from dataset import (
    KmerVocabulary,
    build_dataloaders,
    compute_cksnap,
    compute_pseeiip,
    compute_pstnp,
    compute_tnc,
    random_base_deletion,
    random_base_insertion,
    random_base_substitution,
    reverse_complement,
)
from model import GeneWhispererStage1, DNAEncoder, MultiScaleEnsemble
from numerics import cap_model_param_norm_
from seed_utils import set_global_seed
from stage1_settings import log_stage1_training_configuration, resolve_stage1_lr

LOGGER = logging.getLogger("gene_whisperer.stage1")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


# =============================================================================
# KNOWLEDGE DISTILLATION UTILITIES
# =============================================================================

class DistillationLoss(nn.Module):
    """
    Knowledge distillation loss for binary classification.

    Combines standard cross-entropy with KL divergence between teacher/student
    logits, using temperature scaling to soften the distributions.

    For binary classification with logits, we convert logits to 2-class
    distributions: [1-sigmoid(logit), sigmoid(logit)] for proper KL.

    Loss = (1-alpha)*CE(student, labels) + alpha*T^2*KL(soft_student, soft_teacher)

    The T^2 factor compensates for the gradient magnitude reduction from temperature.
    """

    def __init__(
        self,
        alpha: float = 0.5,
        temperature: float = 2.0,
        base_criterion: Optional[nn.Module] = None,
    ):
        super().__init__()
        self.alpha = alpha
        self.temperature = temperature
        self.base_criterion = base_criterion or nn.BCEWithLogitsLoss()

    def forward(
        self,
        student_logits: torch.Tensor,
        teacher_logits: torch.Tensor,
        labels: torch.Tensor,
    ) -> Tuple[torch.Tensor, Dict[str, float]]:
        """
        Compute combined CE + KL distillation loss.

        Args:
            student_logits: (B, 1) student pre-sigmoid logits
            teacher_logits: (B, 1) teacher pre-sigmoid logits (detached)
            labels: (B, 1) ground truth labels
        Returns:
            Combined loss and dict with loss components
        """
        # Standard cross-entropy loss on logits
        ce_loss = self.base_criterion(student_logits, labels)

        # Convert binary logits to 2-class log-softmax with temperature
        # For binary: [logit_class0, logit_class1] = [-logit, logit]
        # This makes softmax([−l, l]) = [1−σ(l), σ(l)]
        student_logits_2class = torch.cat([-student_logits, student_logits], dim=-1)
        teacher_logits_2class = torch.cat([-teacher_logits, teacher_logits], dim=-1)

        # Apply temperature scaling
        student_soft = F.log_softmax(student_logits_2class / self.temperature, dim=-1)
        teacher_soft = F.softmax(teacher_logits_2class / self.temperature, dim=-1)

        # KL divergence (reduction='batchmean' is standard for KL)
        kl_loss = F.kl_div(student_soft, teacher_soft, reduction='batchmean')

        # Scale by T^2 to match gradient magnitudes (Hinton et al.)
        kl_loss = kl_loss * (self.temperature ** 2)

        # Combined loss
        total_loss = (1 - self.alpha) * ce_loss + self.alpha * kl_loss

        return total_loss, {
            "ce_loss": ce_loss.item(),
            "kl_loss": kl_loss.item(),
            "total_loss": total_loss.item(),
        }


def load_teacher_model(
    cfg: dict,
    teacher_ckpt_path: Path,
    device: torch.device,
) -> nn.Module:
    """
    Load a frozen teacher model for distillation.

    The teacher is loaded from a Stage 1 checkpoint and set to eval mode
    with all gradients disabled.

    Args:
        cfg: Configuration dict (used for model architecture params)
        teacher_ckpt_path: Path to teacher checkpoint
        device: Device to load model to
    Returns:
        Frozen teacher model in eval mode
    """
    from dataset import KmerVocabulary

    # Load checkpoint to get model configuration
    checkpoint = torch.load(teacher_ckpt_path, map_location=device, weights_only=False)

    # Get model state
    if isinstance(checkpoint, dict) and "model_state" in checkpoint:
        state_dict = checkpoint["model_state"]
    else:
        state_dict = checkpoint

    # Try to infer k-mer from checkpoint path or use config
    teacher_kmer = int(cfg.get("distill", {}).get("teacher_kmer", 6))

    # Build vocab for teacher
    vocab_cache_dir = Path(cfg.get("vocab_cache_dir", "../artifacts/vocabs"))
    mlm_vocab_path = cfg.get("mlm_vocab_path")
    mlm_kmer = int(cfg.get("mlm_kmer", teacher_kmer))
    if mlm_vocab_path and Path(mlm_vocab_path).exists() and mlm_kmer == teacher_kmer:
        teacher_vocab = KmerVocabulary.load(Path(mlm_vocab_path))
    else:
        teacher_vocab_path = vocab_cache_dir / f"k{teacher_kmer}_vocab.json"
        teacher_vocab = KmerVocabulary.load(teacher_vocab_path)

    # Create teacher model with same architecture
    embedding_dim = int(cfg.get("embedding_dim", 256))
    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 1024))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.15))

    teacher = GeneWhispererStage1(
        vocab_size=len(teacher_vocab.itos),
        kmer=teacher_kmer,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=teacher_vocab.pad_id,
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 256)),
        tcn_levels=int(cfg.get("tcn_levels", 4)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        fusion_hidden=int(cfg.get("fusion_hidden", 256)),
        max_seq_len=int(cfg.get("max_bp_len", 81)) - teacher_kmer + 1,
    ).to(device)

    # Load state dict
    # Handle potential key mismatches
    model_dict = teacher.state_dict()
    pretrained_dict = {}

    for k, v in state_dict.items():
        # Remove "module." prefix if present (from DataParallel)
        key = k.replace("module.", "") if k.startswith("module.") else k
        if key in model_dict and model_dict[key].shape == v.shape:
            pretrained_dict[key] = v

    model_dict.update(pretrained_dict)
    teacher.load_state_dict(model_dict)

    LOGGER.info(
        "Loaded teacher model from %s (%d/%d weights)",
        teacher_ckpt_path,
        len(pretrained_dict),
        len(state_dict),
    )

    # Freeze teacher and set to eval
    teacher.eval()
    for param in teacher.parameters():
        param.requires_grad = False

    return teacher


def verify_distillation_dry_run(
    student: nn.Module,
    teacher: nn.Module,
    train_loader: torch.utils.data.DataLoader,
    device: torch.device,
    teacher_vocab: Any,
    student_vocab: Any,
) -> bool:
    """
    Verify that distillation setup is correct with a dry run.

    Checks:
    1. Teacher forward pass works
    2. Logit shapes match between teacher and student
    3. Both produce valid outputs

    Args:
        student: Student model
        teacher: Frozen teacher model
        train_loader: Training dataloader (for student k-mer)
        device: Device
        teacher_vocab: Teacher's vocabulary
        student_vocab: Student's vocabulary
    Returns:
        True if verification passed
    """
    from dataset import KmerVocabulary, PromoterDatasetStage1

    LOGGER.info("=" * 60)
    LOGGER.info("DISTILLATION DRY-RUN VERIFICATION")
    LOGGER.info("=" * 60)

    student.eval()
    teacher.eval()

    # Get a batch from train loader
    batch_iter = iter(train_loader)
    tokens_student, engineered, labels = next(batch_iter)
    tokens_student = tokens_student.to(device)
    engineered = engineered.to(device)
    labels = labels.to(device)

    batch_size = tokens_student.size(0)

    # Get raw sequences from the dataset for re-tokenization
    # We need to tokenize the same sequences with teacher's k-mer
    dataset = train_loader.dataset

    # Get sequences for re-tokenization
    # The dataset should have the raw sequences accessible
    sequences = []
    for i in range(min(batch_size, len(dataset))):
        if hasattr(dataset, 'sequences'):
            sequences.append(dataset.sequences[i])
        elif hasattr(dataset, 'data') and hasattr(dataset.data[i], '__getitem__'):
            sequences.append(dataset.data[i][0])  # Assume (seq, label) format
        else:
            # Fallback: decode from tokens
            seq = student_vocab.decode(tokens_student[i].cpu().tolist())
            sequences.append(seq)

    # Tokenize sequences for teacher
    tokens_teacher = []
    for seq in sequences:
        tok = teacher_vocab.tokenize(seq, max_bp_len=len(seq))
        tokens_teacher.append(torch.tensor(tok, dtype=torch.long))

    # Pad to same length
    max_len = max(t.size(0) for t in tokens_teacher)
    tokens_teacher_padded = torch.zeros(len(tokens_teacher), max_len, dtype=torch.long)
    for i, t in enumerate(tokens_teacher):
        tokens_teacher_padded[i, :t.size(0)] = t

    tokens_teacher = tokens_teacher_padded.to(device)

    LOGGER.info("Student tokens shape: %s", tokens_student.shape)
    LOGGER.info("Teacher tokens shape: %s", tokens_teacher.shape)

    try:
        with torch.no_grad():
            # Student forward with logits
            student_logits = student(tokens_student, engineered)
            student_probs = torch.sigmoid(student_logits)

            # Teacher forward with logits
            teacher_logits = teacher(tokens_teacher, engineered)
            teacher_probs = torch.sigmoid(teacher_logits)

        LOGGER.info("Student logits shape: %s", student_logits.shape)
        LOGGER.info("Teacher logits shape: %s", teacher_logits.shape)
        LOGGER.info("Student probs range: [%.4f, %.4f]", student_probs.min().item(), student_probs.max().item())
        LOGGER.info("Teacher probs range: [%.4f, %.4f]", teacher_probs.min().item(), teacher_probs.max().item())
        LOGGER.info("Student logits range: [%.4f, %.4f]", student_logits.min().item(), student_logits.max().item())
        LOGGER.info("Teacher logits range: [%.4f, %.4f]", teacher_logits.min().item(), teacher_logits.max().item())

        # Check shapes match
        if student_logits.shape != teacher_logits.shape:
            LOGGER.error(
                "Logit shape mismatch: student=%s, teacher=%s",
                student_logits.shape,
                teacher_logits.shape,
            )
            return False

        # Quick distillation loss test
        distill_loss = DistillationLoss(alpha=0.5, temperature=2.0)
        loss, loss_dict = distill_loss(
            student_logits,
            teacher_logits,
            labels.unsqueeze(1),
        )

        LOGGER.info("Test distillation loss: %.4f (CE=%.4f, KL=%.4f)",
                    loss_dict["total_loss"], loss_dict["ce_loss"], loss_dict["kl_loss"])

        if not torch.isfinite(loss):
            LOGGER.error("Distillation loss is not finite: %s", loss.item())
            return False

        LOGGER.info("=" * 60)
        LOGGER.info("✓ DISTILLATION DRY-RUN PASSED")
        LOGGER.info("=" * 60)
        return True

    except Exception as e:
        LOGGER.error("Distillation dry-run failed with error: %s", e)
        import traceback
        traceback.print_exc()
        return False
    finally:
        student.train()
        # Teacher stays in eval mode


def apply_distill_guardrails(
    *,
    distill_enabled: bool,
    teacher_kmer: int,
    student_kmer: int,
    teacher_ckpt_path: Path,
    student_ckpt_path: Path,
) -> tuple[bool, str]:
    """
    Check distillation safety guardrails.

    Prevents distillation in unsafe scenarios:
    1. If distillation is already disabled
    2. If teacher and student have the same k-mer (no knowledge transfer benefit)
    3. If teacher and student checkpoint paths resolve to the same file

    Args:
        distill_enabled: Whether distillation was requested
        teacher_kmer: Teacher model's k-mer size
        student_kmer: Student model's k-mer size
        teacher_ckpt_path: Path to teacher checkpoint
        student_ckpt_path: Path to student checkpoint (output)

    Returns:
        Tuple of (enabled, reason):
        - (False, reason) if distillation should be disabled
        - (True, "") if distillation can proceed
    """
    if not distill_enabled:
        return (False, "distillation not enabled")

    if teacher_kmer == student_kmer:
        return (False, f"teacher k-mer ({teacher_kmer}) equals student k-mer ({student_kmer})")

    if teacher_ckpt_path.resolve() == student_ckpt_path.resolve():
        return (
            False,
            f"teacher checkpoint resolves to same path as student output: {teacher_ckpt_path.resolve()}",
        )

    return (True, "")


# Note: set_seed is now a wrapper for set_global_seed from seed_utils
def set_seed(seed: int) -> None:
    """Set all random seeds for reproducibility (wrapper for set_global_seed)."""
    set_global_seed(seed)


def load_config(path: Path) -> Dict:
    """Load YAML configuration."""
    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    return cfg or {}


def get_with_fallback(cfg: dict, stage_key: str, global_key: str, default: Any) -> Any:
    """Get stage-specific value, falling back to global if missing or null."""
    value = cfg.get(stage_key)
    if value is not None:
        return value
    return cfg.get(global_key, default)


class LabelSmoothingBCE(nn.Module):
    """
    Binary cross entropy with label smoothing.
    
    Label smoothing helps prevent overconfident predictions and
    improves generalization. For binary classification:
    - Smooth 0 → smoothing/2
    - Smooth 1 → 1 - smoothing/2
    
    Expects logits as input.
    """
    
    def __init__(self, smoothing: float = 0.1):
        super().__init__()
        self.smoothing = smoothing
    
    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        # Smooth labels
        target_smooth = target * (1 - self.smoothing) + 0.5 * self.smoothing
        
        # Binary cross entropy on logits
        loss = F.binary_cross_entropy_with_logits(logits, target_smooth, reduction='mean')
        return loss


class FocalLossWithSmoothing(nn.Module):
    """
    Focal loss with optional label smoothing.
    
    Combines focal loss (handles class imbalance) with 
    label smoothing (prevents overconfidence).
    
    Expects logits as input.
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
    
    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        # Label smoothing
        target_smooth = target * (1 - self.smoothing) + 0.5 * self.smoothing
        
        # Focal loss
        probs = torch.sigmoid(logits)
        p_t = probs * target_smooth + (1 - probs) * (1 - target_smooth)
        alpha_t = self.alpha * target_smooth + (1 - self.alpha) * (1 - target_smooth)
        
        focal_weight = (1 - p_t) ** self.gamma
        loss = -alpha_t * focal_weight * torch.log(p_t.clamp(min=self.eps))
        
        return loss.mean()


# =============================================================================
# MULTI-SCALE ENSEMBLE UTILITIES
# =============================================================================

class MultiKmerStage1Dataset(Dataset):
    """Stage 1 dataset that returns tokens for multiple k-mer vocabularies."""

    def __init__(
        self,
        sequences: List[str],
        labels: List[float],
        max_bp_len: int,
        vocabs: Dict[int, KmerVocabulary],
        use_engineered_features: bool = True,
        feature_enable_tnc: bool = True,
        feature_enable_pseeiip: bool = True,
        feature_enable_cksnap: bool = True,
        feature_enable_pstnp: bool = True,
        engineered_dim: int = 288,
        reverse_complement_prob: float = 0.5,
        base_substitution_prob: float = 0.0,
        base_indel_prob: float = 0.0,
        cache_engineered_features: bool = False,
        cache_tokens: bool = False,
    ):
        if len(sequences) != len(labels):
            raise ValueError(
                f"Sequences/labels length mismatch: {len(sequences)} vs {len(labels)}"
            )
        if not vocabs:
            raise ValueError("MultiKmerStage1Dataset requires at least one vocabulary")

        self.sequences = list(sequences)
        self.labels = list(labels)
        self.max_bp_len = int(max_bp_len)
        self.vocabs = {int(k): vocab for k, vocab in vocabs.items()}
        self.use_engineered_features = bool(use_engineered_features)
        self.feature_enable_tnc = bool(feature_enable_tnc)
        self.feature_enable_pseeiip = bool(feature_enable_pseeiip)
        self.feature_enable_cksnap = bool(feature_enable_cksnap)
        self.feature_enable_pstnp = bool(feature_enable_pstnp)
        self.engineered_dim = int(engineered_dim)
        self.reverse_complement_prob = max(0.0, min(1.0, float(reverse_complement_prob)))
        self.base_substitution_prob = max(0.0, min(0.2, float(base_substitution_prob)))
        self.base_indel_prob = max(0.0, min(0.1, float(base_indel_prob)))
        self._zero_engineered = torch.zeros(self.engineered_dim, dtype=torch.float32)

        self.cache_engineered_features = bool(cache_engineered_features)
        self.cache_tokens = bool(cache_tokens)
        self._cached_engineered: Optional[List[torch.FloatTensor]] = None
        self._cached_tokens: Optional[List[Dict[int, torch.LongTensor]]] = None

        has_augmentation = (
            self.reverse_complement_prob > 0
            or self.base_substitution_prob > 0
            or self.base_indel_prob > 0
        )

        if self.cache_engineered_features:
            if has_augmentation:
                LOGGER.warning(
                    "cache_engineered_features=True but augmentation is enabled; "
                    "caching is less effective with augmentation"
                )
            self._precompute_engineered_features()

        if self.cache_tokens:
            if has_augmentation:
                LOGGER.warning(
                    "cache_tokens=True but augmentation is enabled; "
                    "caching is less effective with augmentation"
                )
            self._precompute_tokens()

    def __len__(self) -> int:
        return len(self.sequences)

    def _compute_engineered_features(self, sequence: str) -> torch.FloatTensor:
        tnc = compute_tnc(sequence)
        pseeiip = compute_pseeiip(sequence)
        cksnap = compute_cksnap(sequence)
        pstnp = compute_pstnp(sequence)
        if not self.feature_enable_tnc:
            tnc = torch.zeros_like(tnc)
        if not self.feature_enable_pseeiip:
            pseeiip = torch.zeros_like(pseeiip)
        if not self.feature_enable_cksnap:
            cksnap = torch.zeros_like(cksnap)
        if not self.feature_enable_pstnp:
            pstnp = torch.zeros_like(pstnp)
        engineered = torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)
        if engineered.shape[-1] != self.engineered_dim:
            raise ValueError(
                f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
                f"expected {self.engineered_dim}"
            )
        return engineered

    def _build_engineered_features(self, sequence: str) -> torch.FloatTensor:
        if not self.use_engineered_features:
            engineered = self._zero_engineered.clone()
            if engineered.shape[-1] != self.engineered_dim:
                raise ValueError(
                    f"Engineered feature dim mismatch: got {engineered.shape[-1]}, "
                    f"expected {self.engineered_dim}"
                )
            return engineered
        return self._compute_engineered_features(sequence)

    def _augment_sequence(self, sequence: str) -> str:
        if self.reverse_complement_prob > 0 and random.random() < self.reverse_complement_prob:
            sequence = reverse_complement(sequence)
        if self.base_substitution_prob > 0:
            sequence = random_base_substitution(sequence, self.base_substitution_prob)
        if self.base_indel_prob > 0:
            if random.random() < 0.5:
                sequence = random_base_deletion(sequence, self.base_indel_prob)
            else:
                sequence = random_base_insertion(sequence, self.base_indel_prob)
        return sequence

    def _tokenize_all(self, sequence: str) -> Dict[int, torch.LongTensor]:
        return {
            k: vocab.tokenize(sequence, self.max_bp_len)
            for k, vocab in self.vocabs.items()
        }

    def _precompute_engineered_features(self) -> None:
        LOGGER.info("Precomputing engineered features for %d sequences...", len(self.sequences))
        self._cached_engineered = []
        for seq in self.sequences:
            engineered = self._build_engineered_features(seq)
            self._cached_engineered.append(engineered)
        LOGGER.info("Engineered features cached.")

    def _precompute_tokens(self) -> None:
        LOGGER.info("Pre-tokenizing %d sequences for %d k-mers...", len(self.sequences), len(self.vocabs))
        self._cached_tokens = []
        for seq in self.sequences:
            tokens_by_k = self._tokenize_all(seq)
            self._cached_tokens.append(tokens_by_k)
        LOGGER.info("Tokens cached.")

    def __getitem__(self, idx: int):
        sequence = self.sequences[idx]
        has_augmentation = (
            self.reverse_complement_prob > 0
            or self.base_substitution_prob > 0
            or self.base_indel_prob > 0
        )

        if has_augmentation:
            sequence = self._augment_sequence(sequence)
            tokens_by_k = self._tokenize_all(sequence)
            engineered = self._build_engineered_features(sequence)
        else:
            if self._cached_tokens is not None:
                tokens_by_k = self._cached_tokens[idx]
            else:
                tokens_by_k = self._tokenize_all(sequence)

            if self._cached_engineered is not None:
                engineered = self._cached_engineered[idx]
            else:
                engineered = self._build_engineered_features(sequence)

        label = torch.tensor(float(self.labels[idx]), dtype=torch.float32)
        return tokens_by_k, engineered, label


def make_weighted_sampler(
    labels: List[float],
    imbalance_threshold: float = 0.2,
) -> Optional[WeightedRandomSampler]:
    """Create a weighted sampler when class imbalance is large."""
    total = len(labels)
    if total == 0:
        return None
    pos = sum(1 for label in labels if label >= 0.5)
    neg = total - pos
    if pos == 0 or neg == 0:
        return None
    pos_fraction = pos / total
    if abs(pos_fraction - 0.5) <= imbalance_threshold:
        return None
    weights = torch.zeros(total, dtype=torch.double)
    neg_weight = total / (2.0 * neg)
    pos_weight = total / (2.0 * pos)
    for idx, label in enumerate(labels):
        weights[idx] = pos_weight if label >= 0.5 else neg_weight
    LOGGER.info(
        "Using WeightedRandomSampler (pos_frac=%.3f, threshold=%.2f)",
        pos_fraction,
        imbalance_threshold,
    )
    return WeightedRandomSampler(weights, num_samples=total, replacement=True)


def combine_multiscale_losses(
    ensemble_loss: torch.Tensor,
    individual_losses: List[torch.Tensor],
    individual_loss_weight: float,
) -> torch.Tensor:
    """Combine ensemble and per-model losses for co-training."""
    if individual_loss_weight < 0:
        raise ValueError("individual_loss_weight must be >= 0")
    if not individual_losses:
        return ensemble_loss
    return ensemble_loss + individual_loss_weight * torch.stack(individual_losses).sum()


def _normalize_kmers(kmers: List[int]) -> List[int]:
    normalized: List[int] = []
    seen = set()
    for kmer in kmers:
        kmer_int = int(kmer)
        if kmer_int not in seen:
            normalized.append(kmer_int)
            seen.add(kmer_int)
    if not normalized:
        raise ValueError("multi_scale_kmers must contain at least one k-mer size")
    return normalized


def _resolve_ensemble_weights(cfg_run: dict, kmers: List[int]) -> Optional[List[float]]:
    raw = cfg_run.get("ensemble_weights")
    if raw is None:
        return None
    if isinstance(raw, dict):
        weights: List[float] = []
        for kmer in kmers:
            value = raw.get(kmer)
            if value is None:
                value = raw.get(str(kmer))
            if value is None:
                raise ValueError(f"ensemble_weights missing entry for k={kmer}")
            weights.append(float(value))
        return weights
    if isinstance(raw, (list, tuple)):
        if len(raw) != len(kmers):
            raise ValueError("ensemble_weights length must match multi_scale_kmers")
        return [float(value) for value in raw]
    raise ValueError("ensemble_weights must be a dict or list when provided")


def _build_loader_kwargs(cfg_run: dict, device: torch.device) -> Dict[str, Any]:
    num_workers = int(cfg_run.get("num_workers", 0))
    persistent_workers = bool(cfg_run.get("persistent_workers", False)) and num_workers > 0
    prefetch_factor_raw = cfg_run.get("prefetch_factor", None)
    prefetch_factor = int(prefetch_factor_raw) if prefetch_factor_raw is not None and num_workers > 0 else None
    pin_memory = bool(cfg_run.get("pin_memory", False)) and device.type == "cuda"

    if num_workers > 0:
        LOGGER.info(
            "DataLoader config: num_workers=%d, persistent_workers=%s, prefetch_factor=%s, pin_memory=%s",
            num_workers,
            persistent_workers,
            prefetch_factor,
            pin_memory,
        )

    loader_kwargs: Dict[str, Any] = {
        "num_workers": num_workers,
        "pin_memory": pin_memory,
    }
    if persistent_workers:
        loader_kwargs["persistent_workers"] = True
    if prefetch_factor is not None:
        loader_kwargs["prefetch_factor"] = prefetch_factor
    return loader_kwargs


def _build_multiscale_vocabs(
    cfg_run: dict,
    kmers: List[int],
    sequences: List[str],
    script_dir: Path,
) -> Dict[int, KmerVocabulary]:
    vocab_cache_dir = Path(cfg_run.get("vocab_cache_dir", "../artifacts/vocabs"))
    if not vocab_cache_dir.is_absolute():
        vocab_cache_dir = (script_dir / vocab_cache_dir).resolve()
    vocab_cache_dir.mkdir(parents=True, exist_ok=True)

    mlm_vocab_path = cfg_run.get("mlm_vocab_path")
    mlm_kmer = int(cfg_run.get("mlm_kmer", 0))
    mlm_vocab_resolved: Optional[Path] = None
    if mlm_vocab_path:
        candidate = Path(mlm_vocab_path)
        if not candidate.is_absolute():
            candidate = (script_dir / candidate).resolve()
        if candidate.exists():
            mlm_vocab_resolved = candidate
        else:
            LOGGER.warning("mlm_vocab_path not found: %s", candidate)

    vocabs: Dict[int, KmerVocabulary] = {}
    for kmer in kmers:
        if mlm_vocab_resolved is not None and mlm_kmer == kmer:
            vocab = KmerVocabulary.load(mlm_vocab_resolved)
            if vocab.k != kmer:
                raise ValueError(
                    f"MLM vocab k={vocab.k} does not match requested k={kmer}"
                )
            cache_path = vocab_cache_dir / f"k{kmer}_vocab.json"
            try:
                vocab.save(cache_path)
            except Exception:
                pass
            LOGGER.info(
                "Loaded MLM vocab for k=%d from %s (vocab_size=%d)",
                kmer,
                mlm_vocab_resolved,
                len(vocab.itos),
            )
        else:
            cache_path = vocab_cache_dir / f"k{kmer}_vocab.json"
            if cache_path.exists():
                vocab = KmerVocabulary.load(cache_path)
                if vocab.k != kmer:
                    raise ValueError(
                        f"Vocabulary at {cache_path} was built for k={vocab.k}, expected {kmer}"
                    )
            else:
                vocab = KmerVocabulary.build_from_sequences(sequences, kmer)
                vocab.save(cache_path)
                LOGGER.info("Built k=%d vocabulary with %d entries", kmer, len(vocab.itos))
        vocabs[kmer] = vocab
    return vocabs


def _build_multiscale_dataloaders(
    cfg_run: dict,
    kmers: List[int],
    device: torch.device,
    script_dir: Path,
) -> Tuple[DataLoader, DataLoader, Dict[int, KmerVocabulary]]:
    base_cfg = dict(cfg_run)
    base_cfg["kmer"] = kmers[0]
    # Avoid MLM vocab guardrails here; multi-scale handles MLM loading per k-mer.
    base_cfg["stage1_load_mlm_weights"] = False
    base_dataloaders = build_dataloaders(base_cfg)
    train_base = base_dataloaders["stage1"]["train"].dataset
    val_base = base_dataloaders["stage1"]["val"].dataset

    sequences_train = getattr(train_base, "sequences", None)
    labels_train = getattr(train_base, "labels", None)
    sequences_val = getattr(val_base, "sequences", None)
    labels_val = getattr(val_base, "labels", None)

    if sequences_train is None or labels_train is None:
        raise ValueError("Stage1 train dataset missing sequences/labels")
    if sequences_val is None or labels_val is None:
        raise ValueError("Stage1 val dataset missing sequences/labels")

    vocabs = _build_multiscale_vocabs(cfg_run, kmers, sequences_train, script_dir)

    max_bp_len = int(cfg_run.get("max_bp_len", 81))
    stage1_use_engineered = bool(cfg_run.get("stage1_use_engineered_features", True))
    stage1_feature_enable_tnc = bool(cfg_run.get("stage1_feature_enable_tnc", True))
    stage1_feature_enable_pseeiip = bool(cfg_run.get("stage1_feature_enable_pseeiip", True))
    stage1_feature_enable_cksnap = bool(cfg_run.get("stage1_feature_enable_cksnap", True))
    stage1_feature_enable_pstnp = bool(cfg_run.get("stage1_feature_enable_pstnp", True))
    engineered_dim = int(cfg_run.get("engineered_dim", 288))
    stage1_rc_prob = float(cfg_run.get("stage1_reverse_complement_prob", 0.5))
    base_sub_prob = float(cfg_run.get("stage1_base_substitution_prob", 0.0))
    base_indel_prob = float(cfg_run.get("stage1_base_indel_prob", 0.0))
    cache_engineered_features = bool(cfg_run.get("cache_engineered_features", False))
    cache_tokens = bool(cfg_run.get("cache_tokens", False))

    train_dataset = MultiKmerStage1Dataset(
        sequences_train,
        labels_train,
        max_bp_len,
        vocabs=vocabs,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        feature_enable_cksnap=stage1_feature_enable_cksnap,
        feature_enable_pstnp=stage1_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        reverse_complement_prob=stage1_rc_prob,
        base_substitution_prob=base_sub_prob,
        base_indel_prob=base_indel_prob,
        cache_engineered_features=False,
        cache_tokens=False,
    )

    val_dataset = MultiKmerStage1Dataset(
        sequences_val,
        labels_val,
        max_bp_len,
        vocabs=vocabs,
        use_engineered_features=stage1_use_engineered,
        feature_enable_tnc=stage1_feature_enable_tnc,
        feature_enable_pseeiip=stage1_feature_enable_pseeiip,
        feature_enable_cksnap=stage1_feature_enable_cksnap,
        feature_enable_pstnp=stage1_feature_enable_pstnp,
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.0,
        base_substitution_prob=0.0,
        base_indel_prob=0.0,
        cache_engineered_features=cache_engineered_features,
        cache_tokens=cache_tokens,
    )

    batch_size = int(cfg_run.get("batch_size", 64))
    loader_kwargs = _build_loader_kwargs(cfg_run, device)
    sampler = make_weighted_sampler(train_dataset.labels)

    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=sampler is None,
        sampler=sampler,
        **loader_kwargs,
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        **loader_kwargs,
    )

    return train_loader, val_loader, vocabs


def _build_stage1_model_for_kmer(
    cfg_run: dict,
    vocab: KmerVocabulary,
    device: torch.device,
) -> GeneWhispererStage1:
    embedding_dim = int(cfg_run.get("embedding_dim", 192))
    transformer_layers = int(cfg_run.get("transformer_layers", 6))
    transformer_heads = int(cfg_run.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg_run.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg_run.get("transformer_dropout", 0.15))

    use_attention_pool = bool(cfg_run.get("use_attention_pool", True))
    use_engineered_features = bool(cfg_run.get("stage1_use_engineered_features", True))
    engineered_dim = int(cfg_run.get("engineered_dim", 288))

    use_tcn = bool(cfg_run.get("use_tcn", True))
    tcn_hidden = int(cfg_run.get("tcn_hidden", 256))
    tcn_levels = int(cfg_run.get("tcn_levels", 4))
    tcn_kernel = int(cfg_run.get("tcn_kernel", 3))
    multiscale_channels = int(cfg_run.get("multiscale_channels", 64))
    multiscale_kernels_cfg = cfg_run.get("multiscale_kernels", [3, 5, 7, 9, 15])
    multiscale_kernels = tuple(multiscale_kernels_cfg)
    lstm_hidden = int(cfg_run.get("lstm_hidden", 192))

    post_cnn_transformer_layers = int(cfg_run.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg_run.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg_run.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg_run.get("fusion_hidden", 256))
    stage1_drop_path_rate = get_with_fallback(cfg_run, "stage1_drop_path_rate", "drop_path_rate", 0.1)

    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        engineered_dim=engineered_dim,
        use_engineered_features=use_engineered_features,
        use_attention_pool=use_attention_pool,
        use_tcn=use_tcn,
        tcn_hidden=tcn_hidden,
        tcn_levels=tcn_levels,
        tcn_kernel=tcn_kernel,
        multiscale_channels=multiscale_channels,
        multiscale_kernels=multiscale_kernels,
        lstm_hidden=lstm_hidden,
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        fusion_hidden=fusion_hidden,
        drop_path_rate=stage1_drop_path_rate,
        max_seq_len=int(cfg_run.get("max_bp_len", 81)) - vocab.k + 1,
    ).to(device)

    return model


def _maybe_load_mlm_weights(
    model: GeneWhispererStage1,
    cfg_run: dict,
    *,
    kmer: int,
    script_dir: Path,
) -> bool:
    load_mlm_weights = bool(cfg_run.get("stage1_load_mlm_weights", False))
    transfer_mode = str(cfg_run.get("stage1_mlm_transfer_mode", "embed_only"))
    if transfer_mode not in {"none", "embed_only", "embed_plus_adapter"}:
        raise ValueError(f"Invalid stage1_mlm_transfer_mode: {transfer_mode}")

    if not load_mlm_weights or transfer_mode == "none":
        return False

    mlm_kmer = int(cfg_run.get("mlm_kmer", kmer))
    if kmer != mlm_kmer:
        LOGGER.info(
            "Skipping MLM load for k=%d (mlm_kmer=%d)",
            kmer,
            mlm_kmer,
        )
        return False

    mlm_checkpoint = cfg_run.get("mlm_encoder_path") or cfg_run.get("mlm_encoder_checkpoint")
    if not mlm_checkpoint:
        LOGGER.warning("stage1_load_mlm_weights=true but no MLM checkpoint configured")
        return False

    checkpoint_path = Path(mlm_checkpoint)
    if not checkpoint_path.is_absolute():
        checkpoint_path = (script_dir / checkpoint_path).resolve()
    if not checkpoint_path.exists():
        LOGGER.warning("MLM checkpoint not found at %s", checkpoint_path)
        return False

    LOGGER.info(
        "Loading MLM encoder weights for k=%d from %s (transfer_mode=%s)",
        kmer,
        checkpoint_path,
        transfer_mode,
    )
    model.load_pretrained_weights(checkpoint_path, strict=False, transfer_mode=transfer_mode)

    freeze_layers = int(cfg_run.get("stage1_freeze_layers", 0))
    if freeze_layers > 0:
        LOGGER.info("Freezing bottom %d encoder layers for k=%d", freeze_layers, kmer)
        model.encoder.freeze_bottom_layers(freeze_layers)

    return True


def _resolve_stage1_ckpt_paths(
    cfg_run: dict,
    kmers: List[int],
    script_dir: Path,
) -> Dict[int, Path]:
    ckpt_by_k = cfg_run.get("stage1_ckpt_by_k")
    resolved: Dict[int, Path] = {}
    for kmer in kmers:
        candidate = None
        if isinstance(ckpt_by_k, dict):
            candidate = ckpt_by_k.get(kmer) or ckpt_by_k.get(str(kmer))
        if not candidate:
            candidate = f"../artifacts/checkpoints/stage1_k{kmer}.pt"
        path = Path(candidate)
        if not path.is_absolute():
            path = (script_dir / path).resolve()
        resolved[kmer] = path
    return resolved


def _current_ensemble_weights(
    ensemble: MultiScaleEnsemble,
    kmers: List[int],
) -> List[float]:
    if ensemble.learnable_weights:
        weights = torch.softmax(ensemble.weights.detach().cpu(), dim=0)
        return [float(value) for value in weights.tolist()]
    if ensemble.weights is None:
        return [1.0 / len(kmers)] * len(kmers)
    return [float(value) for value in ensemble.weights.detach().cpu().tolist()]


def save_multiscale_checkpoints(
    checkpoint_paths: Dict[int, Path],
    models: Dict[int, nn.Module],
    optimizer: torch.optim.Optimizer,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler],
    epoch: int,
    metrics: Dict[str, float],
) -> None:
    for kmer, model in models.items():
        path = checkpoint_paths[kmer]
        save_checkpoint(path, model, optimizer, scheduler, epoch, metrics)


def save_ensemble_weights(
    path: Path,
    ensemble: MultiScaleEnsemble,
    kmers: List[int],
) -> None:
    payload = {
        "kmers": kmers,
        "weights": _current_ensemble_weights(ensemble, kmers),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def run_multiscale_epoch(
    loader: DataLoader,
    models: Dict[int, nn.Module],
    ensemble: MultiScaleEnsemble,
    ordered_kmers: List[int],
    criterion: nn.Module,
    device: torch.device,
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler] = None,
    desc: str = "Train",
    grad_accum_steps: int = 1,
    max_optimizer_steps: Optional[int] = None,
    grad_clip_norm: float = 1.0,
    amp_ctx: Optional[AMPContext] = None,
    scaler: Optional[torch.cuda.amp.GradScaler] = None,
    param_norm_warn: float = 150.0,
    param_norm_cap: float = 200.0,
    individual_loss_weight: float = 0.2,
) -> Dict[str, float]:
    """Run one epoch of multi-scale co-training."""
    is_train = optimizer is not None
    ensemble.train(is_train)

    use_autocast = amp_ctx is not None and amp_ctx.enabled
    use_cuda_scaler = scaler is not None and device.type == "cuda"
    force_fp32_loss = amp_ctx.force_fp32_loss if amp_ctx else True

    if is_train and max_optimizer_steps is not None and max_optimizer_steps <= 0:
        return {
            "accuracy": 0.0,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
            "mcc": 0.0,
            "roc_auc": 0.5,
            "specificity": 0.0,
            "loss": 0.0,
            "ensemble_loss": 0.0,
            "individual_loss": 0.0,
            "optimizer_steps": 0,
        }

    total_loss = 0.0
    total_ensemble_loss = 0.0
    total_individual_loss = 0.0
    total_samples = 0
    all_probs: List[torch.Tensor] = []
    all_labels: List[torch.Tensor] = []

    optimizer_step_count = 0
    prev_param_norm: Optional[float] = None
    hit_max_steps = False

    iterable = tqdm(loader, desc=desc, leave=False)
    for batch_idx, (tokens_by_k, engineered, labels) in enumerate(iterable):
        engineered = engineered.to(device)
        labels = labels.to(device).unsqueeze(1)

        tokens_for_k = {}
        for kmer in ordered_kmers:
            tokens = tokens_by_k.get(kmer)
            if tokens is None:
                raise KeyError(f"Missing tokens for k={kmer} in multiscale batch")
            tokens_for_k[kmer] = tokens.to(device)

        with torch.set_grad_enabled(is_train):
            if use_autocast and amp_ctx is not None:
                with amp_ctx.autocast():
                    logits_list = [
                        models[kmer](tokens_for_k[kmer], engineered)
                        for kmer in ordered_kmers
                    ]
            else:
                logits_list = [
                    models[kmer](tokens_for_k[kmer], engineered)
                    for kmer in ordered_kmers
                ]

            if force_fp32_loss:
                logits_for_loss = [logits.float() for logits in logits_list]
            else:
                logits_for_loss = logits_list

            individual_losses = [
                criterion(logits, labels) for logits in logits_for_loss
            ]

            # Manual soft voting to avoid a second forward pass through the ensemble.
            stacked = torch.stack(logits_for_loss, dim=0)
            probs = torch.sigmoid(stacked)
            if ensemble.learnable_weights:
                weights = torch.softmax(ensemble.weights, dim=0)
                weights = weights.to(device=probs.device, dtype=probs.dtype)
                probs = (weights.view(-1, 1, 1) * probs).sum(dim=0)
            elif ensemble.weights is None:
                probs = probs.mean(dim=0)
            else:
                weights = ensemble.weights.to(device=probs.device, dtype=probs.dtype)
                probs = (weights.view(-1, 1, 1) * probs).sum(dim=0)

            eps = torch.finfo(probs.dtype).eps
            ensemble_logits = torch.logit(probs.clamp(min=eps, max=1 - eps))
            ensemble_loss = criterion(ensemble_logits, labels)
            combined_loss = combine_multiscale_losses(
                ensemble_loss, individual_losses, individual_loss_weight
            )

            if is_train:
                loss_scaled = combined_loss / grad_accum_steps
                if use_cuda_scaler:
                    scaler.scale(loss_scaled).backward()
                else:
                    loss_scaled.backward()

                if (batch_idx + 1) % grad_accum_steps == 0:
                    if grad_clip_norm > 0:
                        if use_cuda_scaler:
                            scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(ensemble.parameters(), grad_clip_norm)

                    if use_cuda_scaler:
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        optimizer.step()

                    total_norm, scale = cap_model_param_norm_(ensemble, max_norm=param_norm_cap)
                    if scale < 1.0:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                            optimizer_step_count,
                            total_norm,
                            param_norm_cap,
                            scale,
                        )
                        prev_param_norm = param_norm_cap
                    else:
                        if (
                            param_norm_warn > 0
                            and (prev_param_norm is None or prev_param_norm <= param_norm_warn)
                            and total_norm > param_norm_warn
                        ):
                            LOGGER.warning(
                                "Step %d: param_norm=%.2f exceeded warn=%.2f",
                                optimizer_step_count,
                                total_norm,
                                param_norm_warn,
                            )
                        prev_param_norm = total_norm

                    optimizer.zero_grad()

                    if scheduler is not None:
                        scheduler.step()

                    optimizer_step_count += 1
                    if max_optimizer_steps is not None and optimizer_step_count >= max_optimizer_steps:
                        hit_max_steps = True

        batch_size = labels.size(0)
        total_loss += combined_loss.item() * batch_size
        total_ensemble_loss += ensemble_loss.item() * batch_size
        total_individual_loss += sum(loss.item() for loss in individual_losses) * batch_size
        total_samples += batch_size

        all_probs.append(probs.detach().float().cpu())
        all_labels.append(labels.detach().cpu())

        iterable.set_postfix(loss=combined_loss.item())

        if hit_max_steps:
            LOGGER.info("%s: reached max optimizer steps (%d)", desc, max_optimizer_steps)
            break

    avg_loss = total_loss / total_samples if total_samples else 0.0
    avg_ensemble_loss = total_ensemble_loss / total_samples if total_samples else 0.0
    avg_individual_loss = total_individual_loss / total_samples if total_samples else 0.0

    if total_samples == 0 or not all_probs:
        return {
            "accuracy": 0.0,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
            "mcc": 0.0,
            "roc_auc": 0.5,
            "specificity": 0.0,
            "loss": float(avg_loss),
            "ensemble_loss": float(avg_ensemble_loss),
            "individual_loss": float(avg_individual_loss),
            "optimizer_steps": int(optimizer_step_count),
        }

    probs_tensor = torch.cat(all_probs, dim=0)
    labels_tensor = torch.cat(all_labels, dim=0)
    metrics = calculate_metrics(probs_tensor, labels_tensor)
    metrics["loss"] = float(avg_loss)
    metrics["ensemble_loss"] = float(avg_ensemble_loss)
    metrics["individual_loss"] = float(avg_individual_loss)
    metrics["optimizer_steps"] = int(optimizer_step_count)

    return metrics


def train_multi_scale_ensemble(
    cfg_run: dict,
    device: torch.device,
    amp_ctx: AMPContext,
    scaler: Optional[torch.cuda.amp.GradScaler],
    run_dir: Path,
    resolved_config_path: Path,
    resolved_config: Dict[str, Any],
    *,
    epochs: int,
    ablation_max_steps: Optional[int],
    seed: int,
    config_id: str,
    start_time: float,
) -> Dict[str, Any]:
    """Co-train all k-mer models with a soft voting ensemble loss."""
    kmers = _normalize_kmers(cfg_run.get("multi_scale_kmers", [3, 4, 5, 6]))
    ensemble_method = str(cfg_run.get("ensemble_method", "soft_voting")).lower()
    if ensemble_method not in {"soft_voting", "soft-voting"}:
        raise ValueError(f"Unsupported ensemble_method: {ensemble_method}")

    if bool(cfg_run.get("overfit_debug", False)):
        LOGGER.warning("overfit_debug ignored in multi-scale mode")
    if bool(cfg_run.get("profile", False)):
        LOGGER.warning("profile ignored in multi-scale mode")
    if bool(cfg_run.get("distill", {}).get("enabled", False)):
        LOGGER.warning("distillation ignored in multi-scale mode")

    stage1_use_mixup = get_with_fallback(cfg_run, "stage1_use_mixup", "use_mixup", True)
    if stage1_use_mixup:
        LOGGER.warning("mixup disabled in multi-scale mode")
    stage1_use_swa = get_with_fallback(cfg_run, "stage1_use_swa", "use_swa", True)
    if stage1_use_swa:
        LOGGER.warning("SWA disabled in multi-scale mode")

    max_bp_len = int(cfg_run.get("max_bp_len", 81))
    lr = resolve_stage1_lr(cfg_run)

    script_dir = Path(__file__).resolve().parent
    train_loader, val_loader, vocabs = _build_multiscale_dataloaders(
        cfg_run, kmers, device, script_dir
    )

    batch_size = train_loader.batch_size or int(cfg_run.get("batch_size", 64))
    log_stage1_training_configuration(
        max_bp_len=max_bp_len,
        lr=lr,
        batch_size=batch_size,
        logger=LOGGER,
    )

    models: Dict[int, GeneWhispererStage1] = {}
    for kmer in kmers:
        model = _build_stage1_model_for_kmer(cfg_run, vocabs[kmer], device)
        _maybe_load_mlm_weights(model, cfg_run, kmer=kmer, script_dir=script_dir)
        models[kmer] = model

    total_params = sum(p.numel() for model in models.values() for p in model.parameters())
    trainable_params = sum(
        p.numel() for model in models.values() for p in model.parameters() if p.requires_grad
    )
    LOGGER.info(
        "Multi-scale models: %.2fM total params, %.2fM trainable",
        total_params / 1e6,
        trainable_params / 1e6,
    )

    resolved_config["system_info"]["device"] = str(device)
    resolved_config["system_info"]["stage1_multiscale_total_params"] = total_params
    resolved_config["system_info"]["stage1_multiscale_trainable_params"] = trainable_params
    resolved_config["multi_scale_kmers"] = kmers
    resolved_config["ensemble_method"] = ensemble_method
    resolved_config["ensemble_learn_weights"] = bool(cfg_run.get("ensemble_learn_weights", False))
    resolved_config["ensemble_weights"] = cfg_run.get("ensemble_weights")
    with resolved_config_path.open("w", encoding="utf-8") as handle:
        json.dump(resolved_config, handle, indent=2)

    ensemble_weights = _resolve_ensemble_weights(cfg_run, kmers)
    ensemble = MultiScaleEnsemble(
        [models[kmer] for kmer in kmers],
        weights=ensemble_weights,
        learnable_weights=bool(cfg_run.get("ensemble_learn_weights", False)),
    ).to(device)

    weight_decay = float(cfg_run.get("weight_decay", 0.02))
    optimizer_foreach = bool(cfg_run.get("optimizer_foreach", True))
    optimizer = torch.optim.AdamW(
        ensemble.parameters(),
        lr=lr,
        weight_decay=weight_decay,
        betas=(0.9, 0.999),
        foreach=optimizer_foreach,
    )

    steps_per_epoch = len(train_loader)
    total_steps = epochs * steps_per_epoch
    warmup_ratio = float(cfg_run.get("warmup_ratio", 0.15))
    warmup_steps = int(warmup_ratio * total_steps)
    scheduler = get_cosine_schedule_with_warmup(
        optimizer,
        num_warmup_steps=warmup_steps,
        num_training_steps=total_steps,
        min_lr_ratio=0.02,
    )

    stage1_label_smoothing = get_with_fallback(cfg_run, "stage1_label_smoothing", "label_smoothing", 0.05)
    stage1_use_focal_loss = bool(cfg_run.get("stage1_use_focal_loss", True))
    focal_alpha = float(cfg_run.get("stage1_focal_alpha", 0.5))
    focal_gamma = float(cfg_run.get("stage1_focal_gamma", 2.0))

    if stage1_use_focal_loss:
        criterion = FocalLossWithSmoothing(
            alpha=focal_alpha,
            gamma=focal_gamma,
            smoothing=stage1_label_smoothing,
        )
    else:
        criterion = LabelSmoothingBCE(smoothing=stage1_label_smoothing)

    grad_accum_steps = int(cfg_run.get("grad_accum_steps", 2))
    grad_clip_norm = float(cfg_run.get("grad_clip_norm", 1.0))
    param_norm_warn = float(cfg_run.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg_run.get("param_norm_cap", 200.0))
    patience = int(cfg_run.get("early_stopping_patience", 15))
    individual_loss_weight = float(cfg_run.get("ensemble_individual_loss_weight", 0.2))

    class_balance = compute_class_balance(train_loader)
    LOGGER.info(
        "Class balance: %.2f%% positive (%d/%d)",
        class_balance["positive_ratio"] * 100.0,
        class_balance["positive"],
        class_balance["total"],
    )

    checkpoint_paths = _resolve_stage1_ckpt_paths(cfg_run, kmers, script_dir)
    ensemble_weights_path = run_dir / "stage1_ensemble_weights.json"

    best_val_mcc = -float("inf")
    best_val_acc = 0.0
    best_epoch = 0
    best_train_metrics = None
    best_val_metrics = None
    patience_counter = 0
    global_optimizer_steps = 0

    for epoch in range(1, epochs + 1):
        current_lr = optimizer.param_groups[0]["lr"]
        LOGGER.info("Epoch %d/%d (LR: %.2e)", epoch, epochs, current_lr)

        remaining_steps = None if ablation_max_steps is None else max(ablation_max_steps - global_optimizer_steps, 0)

        train_metrics = run_multiscale_epoch(
            loader=train_loader,
            models=models,
            ensemble=ensemble,
            ordered_kmers=kmers,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            scheduler=scheduler,
            desc="Train-Multi",
            grad_accum_steps=grad_accum_steps,
            max_optimizer_steps=remaining_steps,
            grad_clip_norm=grad_clip_norm,
            amp_ctx=amp_ctx,
            scaler=scaler,
            param_norm_warn=param_norm_warn,
            param_norm_cap=param_norm_cap,
            individual_loss_weight=individual_loss_weight,
        )
        global_optimizer_steps += int(train_metrics.get("optimizer_steps", 0))
        reached_max_steps = ablation_max_steps is not None and global_optimizer_steps >= ablation_max_steps

        val_metrics = run_multiscale_epoch(
            loader=val_loader,
            models=models,
            ensemble=ensemble,
            ordered_kmers=kmers,
            criterion=criterion,
            device=device,
            optimizer=None,
            desc="Val-Multi",
            amp_ctx=amp_ctx,
            individual_loss_weight=individual_loss_weight,
        )

        LOGGER.info(
            "Train - Loss: %.4f (ensemble=%.4f, indiv=%.4f) | Acc: %.4f | F1: %.4f | MCC: %.4f",
            train_metrics["loss"],
            train_metrics.get("ensemble_loss", 0.0),
            train_metrics.get("individual_loss", 0.0),
            train_metrics["accuracy"],
            train_metrics["f1"],
            train_metrics["mcc"],
        )
        LOGGER.info(
            "Val   - Loss: %.4f | Acc@0.5: %.4f | Acc@best: %.4f (thr=%.2f) | MCC@best: %.4f",
            val_metrics["loss"],
            val_metrics["accuracy"],
            val_metrics.get("acc_best", val_metrics["accuracy"]),
            val_metrics.get("best_threshold", 0.5),
            val_metrics.get("mcc_best", val_metrics["mcc"]),
        )

        eval_mcc_best = val_metrics.get("mcc_best", val_metrics["mcc"])
        eval_acc_best = val_metrics.get("acc_best", val_metrics["accuracy"])
        improved = False
        if eval_mcc_best > best_val_mcc + 1e-4:
            improved = True
        elif abs(eval_mcc_best - best_val_mcc) < 1e-4 and eval_acc_best > best_val_acc:
            improved = True

        if improved:
            patience_counter = 0
            best_val_mcc = eval_mcc_best
            best_val_acc = eval_acc_best
            best_epoch = epoch
            best_train_metrics = train_metrics
            best_val_metrics = val_metrics

            save_multiscale_checkpoints(
                checkpoint_paths=checkpoint_paths,
                models=models,
                optimizer=optimizer,
                scheduler=scheduler,
                epoch=epoch,
                metrics=val_metrics,
            )
            save_ensemble_weights(ensemble_weights_path, ensemble, kmers)
            LOGGER.info(
                "New best! MCC@best: %.4f, Acc@best: %.4f (%.2f%%)",
                best_val_mcc,
                best_val_acc,
                best_val_acc * 100,
            )
        else:
            patience_counter += 1
            LOGGER.info("No improvement. Patience %d/%d", patience_counter, patience)
            if patience_counter >= patience:
                LOGGER.info("Early stopping triggered at epoch %d", epoch)
                break

        if reached_max_steps:
            LOGGER.info(
                "Reached ablation_max_steps_stage1=%d (optimizer steps); stopping after validation",
                ablation_max_steps,
            )
            break

    if best_val_metrics is None:
        best_train_metrics = train_metrics
        best_val_metrics = val_metrics
        best_epoch = epochs

    report = {
        "config": config_id,
        "best_epoch": best_epoch,
        "train_metrics": best_train_metrics,
        "val_metrics": best_val_metrics,
        "class_balance": class_balance,
        "device": str(device),
        "seed": seed,
        "multi_scale": {
            "kmers": kmers,
            "ensemble_method": ensemble_method,
            "ensemble_learn_weights": bool(cfg_run.get("ensemble_learn_weights", False)),
            "ensemble_weights": _current_ensemble_weights(ensemble, kmers),
            "individual_loss_weight": individual_loss_weight,
        },
        "model_params": {
            "total": total_params,
            "trainable": trainable_params,
        },
    }

    report_path = run_dir / "stage1_ensemble_report.json"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    LOGGER.info("Saved multi-scale training report to %s", report_path)

    stage1_seconds = float(time.perf_counter() - start_time)
    best_val_auc = float(best_val_metrics.get("roc_auc", 0.5)) if best_val_metrics else 0.5

    LOGGER.info("=" * 60)
    LOGGER.info("MULTI-SCALE TRAINING COMPLETE")
    LOGGER.info("Best Epoch: %d", best_epoch)
    LOGGER.info("Best Val Acc@best: %.4f (%.2f%%)", best_val_acc, best_val_acc * 100)
    LOGGER.info("Best Val MCC@best: %.4f", best_val_mcc)
    LOGGER.info("Ensemble weights saved to %s", ensemble_weights_path)
    LOGGER.info("=" * 60)

    return {
        "stage1_acc": float(best_val_acc),
        "stage1_auc": best_val_auc,
        "stage1_mcc": float(best_val_mcc),
        "best_val_acc": float(best_val_acc),
        "best_val_auc": best_val_auc,
        "best_val_mcc": float(best_val_mcc),
        "best_epoch": int(best_epoch),
        "best_ckpt_path": str(ensemble_weights_path),
        "stage1_seconds": stage1_seconds,
        "params": int(total_params),
    }

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


# =============================================================================
# EMBEDDING-SPACE MIXUP UTILITIES
# =============================================================================

# Flag to track if we've logged embedding mixup message this session
_embedding_mixup_logged = False


def unwrap_model(m: nn.Module) -> nn.Module:
    """Unwrap a model from DataParallel/DistributedDataParallel wrapper."""
    return m.module if hasattr(m, "module") else m


def mixup_batch_embeddings(
    model: nn.Module,
    tokens: torch.Tensor,
    engineered: torch.Tensor,
    labels: torch.Tensor,
    alpha: float,
):
    """
    Proper embedding-space mixup for GeneWhispererStage1.

    Returns:
      mixed_embeds: (B, L, D)
      mixed_engineered: (B, engineered_dim)
      y_a: (B, 1)
      y_b: (B, 1)
      lam: float
      padding_mask: Optional[(B, L)] True=pad
    """
    core = unwrap_model(model)

    # Default: no mixup
    if alpha is None or alpha <= 0 or tokens.size(0) < 2:
        with torch.no_grad():
            emb = core.embedding(tokens)
        pad_id = getattr(core, "pad_token_id", None)
        padding_mask = tokens.eq(pad_id) if pad_id is not None else None
        return emb, engineered, labels, labels, 1.0, padding_mask

    lam = float(np.random.beta(alpha, alpha))
    B = tokens.size(0)
    index = torch.randperm(B, device=tokens.device)

    # Embeddings (keep grad so embedding learns under mixup)
    emb_a = core.embedding(tokens)              # (B, L, D)
    emb_b = core.embedding(tokens[index])       # (B, L, D)
    mixed_embeds = lam * emb_a + (1.0 - lam) * emb_b

    # Engineered features
    if engineered is not None:
        engineered_b = engineered[index]
        mixed_engineered = lam * engineered + (1.0 - lam) * engineered_b
    else:
        mixed_engineered = None

    y_a = labels
    y_b = labels[index]

    # Padding mask (conservative OR)
    pad_id = getattr(core, "pad_token_id", None)
    if pad_id is not None:
        mask_a = tokens.eq(pad_id)
        mask_b = tokens[index].eq(pad_id)
        padding_mask = mask_a | mask_b
        if padding_mask.any():
            # Keep fully-padded positions zeroed after mixing.
            mixed_embeds = mixed_embeds.masked_fill(padding_mask.unsqueeze(-1), 0.0)
    else:
        padding_mask = None

    return mixed_embeds, mixed_engineered, y_a, y_b, lam, padding_mask


def mixup_data(
    tokens: torch.Tensor,
    engineered: torch.Tensor,
    labels: torch.Tensor,
    alpha: float = 0.2,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, float]:
    """
    DEPRECATED: Use mixup_batch_embeddings() instead.

    This function only mixed engineered features, not the token embeddings.
    It has been superseded by mixup_batch_embeddings() which properly mixes
    in embedding space.
    """
    import warnings
    warnings.warn(
        "mixup_data() is deprecated and does not properly mix sequences. "
        "Use mixup_batch_embeddings() for embedding-space mixup instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if alpha <= 0:
        return tokens, engineered, labels, labels, 1.0

    lam = np.random.beta(alpha, alpha)
    batch_size = tokens.size(0)
    index = torch.randperm(batch_size, device=tokens.device)

    # Mix tokens (using embedding weights mixing would be better, but this works)
    mixed_engineered = lam * engineered + (1 - lam) * engineered[index]

    # Return both targets for loss computation
    return tokens, mixed_engineered, labels, labels[index], lam


def _binary_roc_auc(probabilities: np.ndarray, targets: np.ndarray) -> float:
    """Compute ROC AUC for binary targets without external dependencies."""
    probs = np.asarray(probabilities, dtype=np.float64).reshape(-1)
    y = np.asarray(targets, dtype=np.float64).reshape(-1)
    y = (y >= 0.5).astype(np.int32)

    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return 0.5

    order = np.argsort(probs, kind="mergesort")  # stable for deterministic tie handling
    sorted_probs = probs[order]
    sorted_y = y[order]

    n = len(sorted_probs)
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_probs[j + 1] == sorted_probs[i]:
            j += 1
        avg_rank = (i + j + 2) / 2.0  # 1-based ranks
        ranks[i : j + 1] = avg_rank
        i = j + 1

    rank_sum_pos = float(ranks[sorted_y == 1].sum())
    auc = (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def _compute_metrics_at_threshold(
    probs: torch.Tensor, targets: torch.Tensor, threshold: float
) -> Dict[str, float]:
    """Compute classification metrics at a specific threshold."""
    preds = (probs >= threshold).float()

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

    specificity = tn / (tn + fp) if (tn + fp) else 0.0

    return {
        "accuracy": float(accuracy),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "mcc": float(mcc),
        "specificity": float(specificity),
    }


def compute_best_threshold(
    probs: torch.Tensor, labels: torch.Tensor, metric: str = "accuracy", num_thresholds: int = 99
) -> Tuple[float, float]:
    """
    Find the optimal decision threshold by sweeping thresholds.

    Args:
        probs: Predicted probabilities (1D tensor)
        labels: Ground truth labels 0/1 (1D tensor)
        metric: Metric to optimize ("accuracy" or "mcc")
        num_thresholds: Number of thresholds to try (default 99)

    Returns:
        (best_threshold, best_metric_value)
    """
    probs = probs.detach().cpu().reshape(-1)
    labels = labels.detach().cpu().reshape(-1)

    thresholds = np.linspace(0.01, 0.99, num_thresholds)
    best_thr = 0.5
    best_value = -float("inf")

    for thr in thresholds:
        metrics = _compute_metrics_at_threshold(probs, labels, thr)
        value = metrics[metric]
        if value > best_value:
            best_value = value
            best_thr = thr

    return float(best_thr), float(best_value)


def calculate_metrics(probabilities: torch.Tensor, labels: torch.Tensor) -> Dict[str, float]:
    """
    Calculate classification metrics at both fixed (0.5) and best thresholds.

    Returns metrics dict with:
    - Standard metrics at threshold=0.5 (accuracy, precision, recall, f1, mcc, etc.)
    - Best threshold metrics: acc_best, mcc_best, best_threshold
    """
    probs = probabilities.detach().cpu().reshape(-1)
    targets = labels.detach().cpu().reshape(-1)

    # Metrics at fixed threshold 0.5
    metrics_05 = _compute_metrics_at_threshold(probs, targets, 0.5)

    # Find best threshold and compute metrics there
    best_thr, _ = compute_best_threshold(probs, targets, metric="accuracy")
    metrics_best = _compute_metrics_at_threshold(probs, targets, best_thr)

    # ROC AUC (threshold-independent)
    roc_auc = _binary_roc_auc(probs.numpy(), targets.numpy())

    return {
        # Metrics at 0.5 threshold (for backward compatibility, these are the "primary" keys)
        "accuracy": metrics_05["accuracy"],
        "precision": metrics_05["precision"],
        "recall": metrics_05["recall"],
        "f1": metrics_05["f1"],
        "mcc": metrics_05["mcc"],
        "specificity": metrics_05["specificity"],
        "roc_auc": float(roc_auc),
        # Best threshold metrics (for honest evaluation)
        "acc_best": metrics_best["accuracy"],
        "mcc_best": metrics_best["mcc"],
        "f1_best": metrics_best["f1"],
        "best_threshold": best_thr,
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
    max_optimizer_steps: Optional[int] = None,
    grad_clip_norm: float = 1.0,
    use_mixup: bool = False,
    mixup_alpha: float = 0.2,
    amp_ctx: Optional[AMPContext] = None,
    scaler: Optional[torch.cuda.amp.GradScaler] = None,
    param_norm_warn: float = 150.0,
    param_norm_cap: float = 200.0,
    eval_checkpoint_path: Optional[Path] = None,
) -> Dict[str, float]:
    """Run one epoch of training or evaluation."""
    is_train = optimizer is not None
    model.train(is_train)

    # AMP setup: use autocast for forward pass, GradScaler for CUDA fp16
    use_autocast = amp_ctx is not None and amp_ctx.enabled
    use_cuda_scaler = scaler is not None and device.type == "cuda"
    force_fp32_loss = amp_ctx.force_fp32_loss if amp_ctx else True

    if is_train and use_cuda_scaler:
        LOGGER.info("%s: GradScaler scale = %.0f", desc, scaler.get_scale())

    if is_train and max_optimizer_steps is not None and max_optimizer_steps <= 0:
        return {
            "accuracy": 0.0,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
            "mcc": 0.0,
            "roc_auc": 0.5,
            "specificity": 0.0,
            "loss": 0.0,
            "optimizer_steps": 0,
        }
    
    total_loss = 0.0
    total_samples = 0
    all_probs = []
    all_labels = []
    
    optimizer_step_count = 0
    accumulated_loss = 0.0
    prev_param_norm: Optional[float] = None
    hit_max_steps = False

    logits_min: Optional[float] = None
    logits_max: Optional[float] = None
    logits_sum = 0.0
    logits_count = 0
    
    iterable = tqdm(loader, desc=desc, leave=False)
    for batch_idx, (tokens, engineered, labels) in enumerate(iterable):
        tokens = tokens.to(device)
        engineered = engineered.to(device)
        labels = labels.to(device).unsqueeze(1)
        
        # Embedding-space mixup augmentation (training only)
        lam = 1.0
        labels_b = labels
        mixed_emb = None
        mix_padding_mask = None
        if is_train and use_mixup:
            mixed_emb, engineered, labels, labels_b, lam, mix_padding_mask = mixup_batch_embeddings(
                model, tokens, engineered, labels, alpha=mixup_alpha
            )

        # Forward pass
        with torch.set_grad_enabled(is_train):
            # Automatic mixed precision: forward under autocast, loss in fp32 for stability
            if use_autocast and amp_ctx is not None:
                with amp_ctx.autocast():
                    if is_train and use_mixup and lam < 1.0:
                        # Use embedding route for mixup
                        core = unwrap_model(model)
                        logits = core.forward_from_embeds(mixed_emb, mix_padding_mask, engineered)
                    else:
                        logits = model(tokens, engineered)
            else:
                if is_train and use_mixup and lam < 1.0:
                    # Use embedding route for mixup
                    core = unwrap_model(model)
                    logits = core.forward_from_embeds(mixed_emb, mix_padding_mask, engineered)
                else:
                    logits = model(tokens, engineered)

            # Ensure logits are in fp32 for loss computation if force_fp32_loss is True
            if force_fp32_loss:
                logits = logits.float()

            if use_mixup and lam < 1.0:
                loss = lam * criterion(logits, labels) + (1 - lam) * criterion(logits, labels_b)
            else:
                loss = criterion(logits, labels)
            
            # Backward pass with gradient accumulation
            if is_train:
                loss_scaled = loss / grad_accum_steps

                if use_cuda_scaler:
                    scaler.scale(loss_scaled).backward()
                else:
                    loss_scaled.backward()

                accumulated_loss += loss.item()

                # Optimizer step every grad_accum_steps
                if (batch_idx + 1) % grad_accum_steps == 0:
                    if grad_clip_norm > 0:
                        if use_cuda_scaler:
                            scaler.unscale_(optimizer)
                        # Gradient clipping - returns the total norm BEFORE clipping
                        unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
                        # After clipping, effective norm is min(unclipped, max_norm)
                        clipped_grad_norm = min(unclipped_grad_norm.item(), grad_clip_norm)
                        
                        # Log gradient norms periodically (every 50 optimizer steps)
                        if optimizer_step_count % 50 == 0:
                            LOGGER.debug(
                                "Step %d: grad_norm (unclipped=%.4f, clipped=%.4f)",
                                optimizer_step_count, unclipped_grad_norm.item(), clipped_grad_norm
                            )
                    
                    if use_cuda_scaler:
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        optimizer.step()

                    total_norm, scale = cap_model_param_norm_(model, max_norm=param_norm_cap)
                    if scale < 1.0:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f exceeded cap=%.2f; scaled params by %.6f",
                            optimizer_step_count,
                            total_norm,
                            param_norm_cap,
                            scale,
                        )
                        prev_param_norm = param_norm_cap
                    else:
                        if (
                            param_norm_warn > 0
                            and (prev_param_norm is None or prev_param_norm <= param_norm_warn)
                            and total_norm > param_norm_warn
                        ):
                            LOGGER.warning(
                                "Step %d: param_norm=%.2f exceeded warn=%.2f",
                                optimizer_step_count,
                                total_norm,
                                param_norm_warn,
                            )
                        prev_param_norm = total_norm
                    
                    optimizer.zero_grad()
                    
                    if scheduler is not None:
                        scheduler.step()
                    
                    optimizer_step_count += 1
                    if max_optimizer_steps is not None and optimizer_step_count >= max_optimizer_steps:
                        hit_max_steps = True
        
        batch_size = labels.size(0)
        total_loss += loss.item() * batch_size
        total_samples += batch_size
        probs = torch.sigmoid(logits.detach().float())
        all_probs.append(probs.cpu())
        all_labels.append(labels.detach().cpu())

        logits_float = logits.detach().float()
        batch_min = float(logits_float.min().item())
        batch_max = float(logits_float.max().item())
        logits_min = batch_min if logits_min is None else min(logits_min, batch_min)
        logits_max = batch_max if logits_max is None else max(logits_max, batch_max)
        logits_sum += float(logits_float.sum().item())
        logits_count += logits_float.numel()
        
        iterable.set_postfix(loss=loss.item())

        if hit_max_steps:
            LOGGER.info("%s: reached max optimizer steps (%d); stopping early", desc, max_optimizer_steps)
            break
    
    avg_loss = total_loss / total_samples if total_samples else 0.0
    if total_samples == 0 or not all_probs:
        if not is_train:
            checkpoint_display = (
                str(eval_checkpoint_path)
                if eval_checkpoint_path is not None
                else "<in-memory (no checkpoint load)>"
            )
            LOGGER.info("what am I evaluating?")
            LOGGER.info("checkpoint path actually loaded: %s", checkpoint_display)
            LOGGER.info("number of val samples: %d", int(total_samples))
            LOGGER.info("label balance (%% positives): n/a")
            LOGGER.info("prob min/mean/max: n/a")
            LOGGER.info("logit min/mean/max: n/a")
        return {
            "accuracy": 0.0,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
            "mcc": 0.0,
            "roc_auc": 0.5,
            "specificity": 0.0,
            "loss": float(avg_loss),
            "optimizer_steps": int(optimizer_step_count),
        }
    probs_tensor = torch.cat(all_probs, dim=0)
    labels_tensor = torch.cat(all_labels, dim=0)
    metrics = calculate_metrics(probs_tensor, labels_tensor)
    metrics["loss"] = float(avg_loss)
    metrics["optimizer_steps"] = int(optimizer_step_count)

    if not is_train:
        checkpoint_display = (
            str(eval_checkpoint_path)
            if eval_checkpoint_path is not None
            else "<in-memory (no checkpoint load)>"
        )
        labels_flat = labels_tensor.view(-1)
        total_eval = int(labels_flat.numel())
        positive_count = int((labels_flat >= 0.5).sum().item())
        positive_pct = (positive_count / total_eval) * 100 if total_eval else 0.0
        probs_flat = probs_tensor.view(-1)
        prob_min = float(probs_flat.min().item())
        prob_mean = float(probs_flat.mean().item())
        prob_max = float(probs_flat.max().item())

        LOGGER.info("what am I evaluating?")
        LOGGER.info("checkpoint path actually loaded: %s", checkpoint_display)
        LOGGER.info("number of val samples: %d", total_eval)
        LOGGER.info(
            "label balance (%% positives): %.2f%% (%d/%d)",
            positive_pct,
            positive_count,
            total_eval,
        )
        LOGGER.info(
            "prob min/mean/max: %.6f/%.6f/%.6f",
            prob_min,
            prob_mean,
            prob_max,
        )
        if logits_count > 0 and logits_min is not None and logits_max is not None:
            logit_mean = logits_sum / logits_count
            LOGGER.info(
                "logit min/mean/max: %.6f/%.6f/%.6f",
                logits_min,
                logit_mean,
                logits_max,
            )
        else:
            LOGGER.info("logit min/mean/max: n/a")
    
    return metrics


def run_distillation_epoch(
    student_loader: torch.utils.data.DataLoader,
    student_model: nn.Module,
    teacher_model: nn.Module,
    distill_criterion: DistillationLoss,
    device: torch.device,
    student_vocab: Any,
    teacher_vocab: Any,
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler: Optional[torch.optim.lr_scheduler._LRScheduler] = None,
    desc: str = "Train-Distill",
    grad_accum_steps: int = 1,
    max_optimizer_steps: Optional[int] = None,
    grad_clip_norm: float = 1.0,
    amp_ctx: Optional[AMPContext] = None,
    scaler: Optional[torch.cuda.amp.GradScaler] = None,
    param_norm_warn: float = 150.0,
    param_norm_cap: float = 200.0,
) -> Dict[str, float]:
    """
    Run one epoch of distillation training.

    The key difference from run_epoch is that we:
    1. Get student tokens from the dataloader
    2. Re-tokenize the same sequences for the teacher (different k-mer)
    3. Run both models and compute distillation loss

    Args:
        student_loader: Dataloader yielding (tokens, engineered, labels) for student
        student_model: Student model being trained
        teacher_model: Frozen teacher model (in eval mode)
        distill_criterion: DistillationLoss instance
        device: Device to train on
        student_vocab: Student's KmerVocabulary for decoding
        teacher_vocab: Teacher's KmerVocabulary for re-tokenizing
        optimizer: Optimizer (None for validation)
        scheduler: LR scheduler
        desc: Description for progress bar
        grad_accum_steps: Gradient accumulation steps
        max_optimizer_steps: Max optimizer steps before stopping
        grad_clip_norm: Gradient clipping norm
        amp_ctx: AMP context for mixed precision
        scaler: CUDA GradScaler
        param_norm_warn: Warning threshold for parameter norm
        param_norm_cap: Cap for parameter norm
    Returns:
        Metrics dict with loss components
    """
    is_train = optimizer is not None
    student_model.train(is_train)
    teacher_model.eval()  # Teacher always in eval mode

    # AMP setup
    use_autocast = amp_ctx is not None and amp_ctx.enabled
    use_cuda_scaler = scaler is not None and device.type == "cuda"
    force_fp32_loss = amp_ctx.force_fp32_loss if amp_ctx else True

    if is_train and max_optimizer_steps is not None and max_optimizer_steps <= 0:
        return {
            "accuracy": 0.0, "precision": 0.0, "recall": 0.0, "f1": 0.0,
            "mcc": 0.0, "roc_auc": 0.5, "specificity": 0.0,
            "loss": 0.0, "ce_loss": 0.0, "kl_loss": 0.0, "optimizer_steps": 0,
        }

    total_loss = 0.0
    total_ce_loss = 0.0
    total_kl_loss = 0.0
    total_samples = 0
    all_probs = []
    all_labels = []

    optimizer_step_count = 0
    prev_param_norm: Optional[float] = None
    hit_max_steps = False

    # Get the dataset for sequence access
    dataset = student_loader.dataset

    iterable = tqdm(student_loader, desc=desc, leave=False)
    for batch_idx, (tokens_student, engineered, labels) in enumerate(iterable):
        batch_size = tokens_student.size(0)
        tokens_student = tokens_student.to(device)
        engineered = engineered.to(device)
        labels = labels.to(device).unsqueeze(1)

        # Re-tokenize sequences for teacher
        # Get raw sequences by decoding student tokens
        sequences = []
        for i in range(batch_size):
            # Decode student tokens back to sequence
            seq = student_vocab.decode(tokens_student[i].cpu().tolist())
            sequences.append(seq)

        # Tokenize for teacher
        tokens_teacher_list = []
        for seq in sequences:
            tok = teacher_vocab.tokenize(seq, max_bp_len=len(seq))
            tokens_teacher_list.append(torch.tensor(tok, dtype=torch.long))

        # Pad teacher tokens to same length
        max_len = max(t.size(0) for t in tokens_teacher_list)
        tokens_teacher = torch.full(
            (batch_size, max_len),
            teacher_vocab.pad_id,
            dtype=torch.long,
            device=device,
        )
        for i, t in enumerate(tokens_teacher_list):
            tokens_teacher[i, :t.size(0)] = t.to(device)

        # Forward pass
        with torch.set_grad_enabled(is_train):
            if use_autocast and amp_ctx is not None:
                with amp_ctx.autocast():
                    student_logits = student_model(tokens_student, engineered)
            else:
                student_logits = student_model(tokens_student, engineered)

            # Teacher forward with logits (no grad needed)
            with torch.no_grad():
                if use_autocast and amp_ctx is not None:
                    with amp_ctx.autocast():
                        teacher_logits = teacher_model(tokens_teacher, engineered)
                else:
                    teacher_logits = teacher_model(tokens_teacher, engineered)

            # Ensure fp32 for loss computation
            if force_fp32_loss:
                student_logits = student_logits.float()
                teacher_logits = teacher_logits.float()

            student_probs = torch.sigmoid(student_logits.detach().float())

            # Compute distillation loss
            loss, loss_dict = distill_criterion(
                student_logits, teacher_logits, labels
            )

            # Backward pass with gradient accumulation
            if is_train:
                loss_scaled = loss / grad_accum_steps

                if use_cuda_scaler:
                    scaler.scale(loss_scaled).backward()
                else:
                    loss_scaled.backward()

                # Optimizer step every grad_accum_steps
                if (batch_idx + 1) % grad_accum_steps == 0:
                    if grad_clip_norm > 0:
                        if use_cuda_scaler:
                            scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(student_model.parameters(), grad_clip_norm)

                    if use_cuda_scaler:
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        optimizer.step()

                    total_norm, scale = cap_model_param_norm_(student_model, max_norm=param_norm_cap)
                    if scale < 1.0:
                        LOGGER.warning(
                            "Step %d: param_norm=%.2f exceeded cap; scaled by %.6f",
                            optimizer_step_count, total_norm, scale,
                        )
                    elif param_norm_warn > 0 and total_norm > param_norm_warn:
                        if prev_param_norm is None or prev_param_norm <= param_norm_warn:
                            LOGGER.warning(
                                "Step %d: param_norm=%.2f exceeded warn=%.2f",
                                optimizer_step_count, total_norm, param_norm_warn,
                            )
                    prev_param_norm = total_norm

                    optimizer.zero_grad()

                    if scheduler is not None:
                        scheduler.step()

                    optimizer_step_count += 1
                    if max_optimizer_steps is not None and optimizer_step_count >= max_optimizer_steps:
                        hit_max_steps = True

        # Accumulate metrics
        total_loss += loss_dict["total_loss"] * batch_size
        total_ce_loss += loss_dict["ce_loss"] * batch_size
        total_kl_loss += loss_dict["kl_loss"] * batch_size
        total_samples += batch_size
        all_probs.append(student_probs.cpu())
        all_labels.append(labels.detach().cpu())

        iterable.set_postfix(
            loss=loss_dict["total_loss"],
            ce=loss_dict["ce_loss"],
            kl=loss_dict["kl_loss"],
        )

        if hit_max_steps:
            LOGGER.info("%s: reached max optimizer steps (%d)", desc, max_optimizer_steps)
            break

    # Compute final metrics
    avg_loss = total_loss / total_samples if total_samples else 0.0
    avg_ce_loss = total_ce_loss / total_samples if total_samples else 0.0
    avg_kl_loss = total_kl_loss / total_samples if total_samples else 0.0

    if total_samples == 0 or not all_probs:
        return {
            "accuracy": 0.0, "precision": 0.0, "recall": 0.0, "f1": 0.0,
            "mcc": 0.0, "roc_auc": 0.5, "specificity": 0.0,
            "loss": float(avg_loss), "ce_loss": float(avg_ce_loss),
            "kl_loss": float(avg_kl_loss), "optimizer_steps": int(optimizer_step_count),
        }

    probs_tensor = torch.cat(all_probs, dim=0)
    labels_tensor = torch.cat(all_labels, dim=0)
    metrics = calculate_metrics(probs_tensor, labels_tensor)
    metrics["loss"] = float(avg_loss)
    metrics["ce_loss"] = float(avg_ce_loss)
    metrics["kl_loss"] = float(avg_kl_loss)
    metrics["optimizer_steps"] = int(optimizer_step_count)

    return metrics


def run_profiling_epoch(
    loader: torch.utils.data.DataLoader,
    model: nn.Module,
    criterion: nn.Module,
    device: torch.device,
    optimizer: torch.optim.Optimizer,
    run_dir: Path,
    num_profile_steps: int = 50,
    grad_accum_steps: int = 1,
    grad_clip_norm: float = 1.0,
) -> Dict[str, float]:
    """
    Run profiling for a limited number of steps.

    Returns timing metrics:
    - avg_dataloader_time_ms: Average time to get a batch
    - avg_forward_backward_time_ms: Average time for forward+backward pass
    - avg_optimizer_step_time_ms: Average time for optimizer step
    """
    model.train()

    timing_stats: Dict[str, List[float]] = defaultdict(list)

    # Determine profiler activities based on device
    activities = [ProfilerActivity.CPU]
    if device.type == "cuda":
        activities.append(ProfilerActivity.CUDA)

    # Profile schedule: wait=5, warmup=5, active=20, repeat=1
    # Total steps needed: (5+5+20) * 1 = 30 steps minimum
    prof_schedule = schedule(
        wait=5,
        warmup=5,
        active=20,
        repeat=1,
    )

    profile_trace_path = run_dir / "profile_trace.json"

    LOGGER.info("=" * 60)
    LOGGER.info("PROFILING MODE")
    LOGGER.info("=" * 60)
    LOGGER.info("Running profiler for %d steps", num_profile_steps)
    LOGGER.info("Profile trace will be saved to: %s", profile_trace_path)

    step_count = 0
    accumulated_loss = 0.0

    with profile(
        activities=activities,
        schedule=prof_schedule,
        on_trace_ready=lambda p: p.export_chrome_trace(str(profile_trace_path)),
        record_shapes=True,
        profile_memory=True,
        with_stack=False,
    ) as prof:
        data_iter = iter(loader)

        while step_count < num_profile_steps:
            # Measure dataloader time
            t_data_start = time.perf_counter()
            try:
                tokens, engineered, labels = next(data_iter)
            except StopIteration:
                data_iter = iter(loader)
                tokens, engineered, labels = next(data_iter)
            t_data_end = time.perf_counter()
            timing_stats["dataloader"].append((t_data_end - t_data_start) * 1000)

            tokens = tokens.to(device)
            engineered = engineered.to(device)
            labels = labels.to(device).unsqueeze(1)

            # Measure forward + backward time
            t_fwd_start = time.perf_counter()
            logits = model(tokens, engineered)
            logits = logits.float()
            loss = criterion(logits, labels)
            loss_scaled = loss / grad_accum_steps
            loss_scaled.backward()
            accumulated_loss += loss.item()
            t_fwd_end = time.perf_counter()
            timing_stats["forward_backward"].append((t_fwd_end - t_fwd_start) * 1000)

            # Measure optimizer step time (every grad_accum_steps)
            if (step_count + 1) % grad_accum_steps == 0:
                t_opt_start = time.perf_counter()
                if grad_clip_norm > 0:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
                optimizer.step()
                optimizer.zero_grad()
                t_opt_end = time.perf_counter()
                timing_stats["optimizer_step"].append((t_opt_end - t_opt_start) * 1000)

            prof.step()
            step_count += 1

            if step_count % 10 == 0:
                LOGGER.info("Profiling step %d/%d", step_count, num_profile_steps)

    # Compute averages
    avg_dataloader = np.mean(timing_stats["dataloader"]) if timing_stats["dataloader"] else 0.0
    avg_fwd_bwd = np.mean(timing_stats["forward_backward"]) if timing_stats["forward_backward"] else 0.0
    avg_opt = np.mean(timing_stats["optimizer_step"]) if timing_stats["optimizer_step"] else 0.0

    LOGGER.info("-" * 60)
    LOGGER.info("PROFILING RESULTS")
    LOGGER.info("-" * 60)
    LOGGER.info("Average dataloader batch time:    %.2f ms", avg_dataloader)
    LOGGER.info("Average forward+backward time:    %.2f ms", avg_fwd_bwd)
    LOGGER.info("Average optimizer step time:      %.2f ms", avg_opt)
    LOGGER.info("Total time per batch (approx):    %.2f ms", avg_dataloader + avg_fwd_bwd)
    LOGGER.info("-" * 60)
    LOGGER.info("Profile trace exported to: %s", profile_trace_path)
    LOGGER.info("Open in Chrome at chrome://tracing/ to visualize")
    LOGGER.info("=" * 60)

    return {
        "avg_dataloader_time_ms": float(avg_dataloader),
        "avg_forward_backward_time_ms": float(avg_fwd_bwd),
        "avg_optimizer_step_time_ms": float(avg_opt),
        "num_profile_steps": num_profile_steps,
        "profile_trace_path": str(profile_trace_path),
    }


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

    # Extract best threshold info from metrics for deployment
    best_threshold = float(metrics.get("best_threshold", 0.5))
    val_acc_best_threshold = float(metrics.get("acc_best", metrics.get("accuracy", 0.0)))

    checkpoint = {
        "epoch": epoch,
        "model_state": model.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "scheduler_state": scheduler.state_dict() if scheduler else None,
        "metrics": metrics,
        # Top-level keys for deployment-safe inference
        "best_threshold": best_threshold,
        "val_acc_best_threshold": val_acc_best_threshold,
    }
    torch.save(checkpoint, path)
    LOGGER.info("Saved checkpoint to %s", path)
    LOGGER.info("  best_threshold=%.4f, val_acc_best_threshold=%.4f", best_threshold, val_acc_best_threshold)


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
        logits = model(tokens, engineered)
        loss = criterion(logits, labels)
        loss.backward()
        
        # Gradient clipping with norm logging
        unclipped_grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        clipped_grad_norm = min(unclipped_grad_norm.item(), 1.0)
        
        optimizer.step()
        
        # Compute accuracy
        probs = torch.sigmoid(logits.detach())
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
    parser.add_argument(
        "--overfit_debug",
        action="store_true",
        help="Train on 32 samples for 200 steps to verify model can learn",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Enable torch.profiler for ~50 steps and export chrome trace",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility. Overrides config value.",
    )
    parser.add_argument(
        "--kmer",
        type=int,
        default=None,
        help="K-mer size (e.g. 3, 4, 6). Overrides config value.",
    )
    parser.add_argument(
        "--vocab_cache_dir",
        type=str,
        default=None,
        help="Path to vocab cache directory. Overrides config value.",
    )
    parser.add_argument(
        "--checkpoint_name",
        type=str,
        default=None,
        help="Checkpoint filename (e.g. stage1_k4.pt). Overrides stage1_checkpoint_name in config.",
    )
    # Architecture overrides for lightweight models
    parser.add_argument(
        "--embedding_dim",
        type=int,
        default=None,
        help="Embedding dimension. Overrides config value.",
    )
    parser.add_argument(
        "--transformer_layers",
        type=int,
        default=None,
        help="Number of transformer encoder layers. Overrides config value.",
    )
    parser.add_argument(
        "--post_cnn_transformer_layers",
        type=int,
        default=None,
        help="Number of post-CNN transformer layers. Overrides config value.",
    )
    parser.add_argument(
        "--tcn_hidden",
        type=int,
        default=None,
        help="TCN hidden dimension. Overrides config value.",
    )
    parser.add_argument(
        "--multiscale_channels",
        type=int,
        default=None,
        help="Multiscale CNN channels. Overrides config value.",
    )
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

    overrides: Dict[str, Any] = {"_config_path": str(config_path)}
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.overfit_debug:
        overrides["overfit_debug"] = True
    if args.profile:
        overrides["profile"] = True
    if args.kmer is not None:
        overrides["kmer"] = args.kmer
    if args.vocab_cache_dir is not None:
        overrides["vocab_cache_dir"] = args.vocab_cache_dir
    if args.checkpoint_name is not None:
        overrides["stage1_checkpoint_name"] = args.checkpoint_name
    # Architecture overrides for lightweight models
    if args.embedding_dim is not None:
        overrides["embedding_dim"] = args.embedding_dim
    if args.transformer_layers is not None:
        overrides["transformer_layers"] = args.transformer_layers
    if args.post_cnn_transformer_layers is not None:
        overrides["post_cnn_transformer_layers"] = args.post_cnn_transformer_layers
    if args.tcn_hidden is not None:
        overrides["tcn_hidden"] = args.tcn_hidden
    if args.multiscale_channels is not None:
        overrides["multiscale_channels"] = args.multiscale_channels

    _ = run_stage1_training(cfg, overrides=overrides)
    return


def run_stage1_training(cfg: dict, *, overrides: dict | None = None) -> dict:
    """Trains Stage1 and returns best val metrics + timing + params."""
    start_time = time.perf_counter()

    cfg_run = copy.copy(cfg) if isinstance(cfg, dict) else dict(cfg or {})
    if overrides:
        cfg_run.update(overrides)

    if bool(cfg_run.get("ablation_disable_wandb", False)):
        os.environ["WANDB_MODE"] = "disabled"
        os.environ["WANDB_DISABLED"] = "true"

    script_dir = Path(__file__).resolve().parent
    config_path_value = cfg_run.get("_config_path") or cfg_run.get("config_path")
    config_id = str(config_path_value) if config_path_value else "<in-memory>"

    ablation_cfg = cfg_run.get("ablation")
    if not isinstance(ablation_cfg, dict):
        ablation_cfg = {}
    ablation_enabled = bool(ablation_cfg.get("enabled", False))

    epochs = int(cfg_run.get("epochs", 120))
    if ablation_enabled and ablation_cfg.get("stage1_epochs") is not None:
        epochs_raw = ablation_cfg.get("stage1_epochs")
        try:
            epochs = int(epochs_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.stage1_epochs must be an int, got {epochs_raw!r}") from exc
        if epochs <= 0:
            raise ValueError(f"ablation.stage1_epochs must be > 0, got {epochs}")

    ablation_max_steps: Optional[int] = None
    if ablation_enabled and ablation_cfg.get("stage1_max_steps") is not None:
        max_steps_raw = ablation_cfg.get("stage1_max_steps")
        try:
            max_steps_int = int(max_steps_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"ablation.stage1_max_steps must be an int, got {max_steps_raw!r}") from exc
        if max_steps_int < 0:
            raise ValueError(f"ablation.stage1_max_steps must be >= 0, got {max_steps_int}")
        if max_steps_int > 0:
            ablation_max_steps = max_steps_int
    else:
        ablation_max_steps_raw = cfg_run.get("ablation_max_steps_stage1")
        if ablation_max_steps_raw is not None:
            try:
                max_steps_int = int(ablation_max_steps_raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"ablation_max_steps_stage1 must be an int, got {ablation_max_steps_raw!r}"
                ) from exc
            if max_steps_int < 0:
                raise ValueError(f"ablation_max_steps_stage1 must be >= 0, got {max_steps_int}")
            if max_steps_int > 0:
                ablation_max_steps = max_steps_int

    stage1_max_train_samples = ablation_cfg.get("stage1_max_train_samples") if ablation_enabled else None
    stage1_max_val_samples = ablation_cfg.get("stage1_max_val_samples") if ablation_enabled else None

    max_bp_len = int(cfg_run.get("max_bp_len", 81))
    lr = resolve_stage1_lr(cfg_run)

    # =========================================================================
    # Log resolved hyperparameters to runs/<timestamp>/resolved_config.json
    # =========================================================================
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = (script_dir / "../runs" / f"stage1_{timestamp}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    
    resolved_config: Dict[str, Any] = {
        "script": "train_stage1.py",
        "timestamp": timestamp,
        "config_path": config_id,
        # System info for reproducibility
        "system_info": {
            "python_version": f"{os.sys.version_info.major}.{os.sys.version_info.minor}.{os.sys.version_info.micro}",
            "torch_version": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
            "mps_available": hasattr(torch.backends, "mps") and torch.backends.mps.is_available(),
        },
        # Ablations / debugging
        "ablation_enabled": ablation_enabled,
        "ablation_max_steps_stage1": ablation_max_steps,
        "ablation_stage1_epochs": epochs if ablation_enabled else None,
        "ablation_stage1_max_train_samples": stage1_max_train_samples if ablation_enabled else None,
        "ablation_stage1_max_val_samples": stage1_max_val_samples if ablation_enabled else None,
        "ablation_disable_wandb": bool(cfg_run.get("ablation_disable_wandb", False)),
        # Model architecture
        "embedding_dim": int(cfg_run.get("embedding_dim", 256)),
        "transformer_layers": int(cfg_run.get("transformer_layers", 6)),
        "transformer_heads": int(cfg_run.get("transformer_heads", 8)),
        "transformer_ff_dim": int(cfg_run.get("transformer_ff_dim", 1024)),
        "transformer_dropout": float(cfg_run.get("transformer_dropout", 0.15)),
        # Optimizer
        "lr": lr,
        "weight_decay": float(cfg_run.get("weight_decay", 0.02)),
        # Data
        "batch_size": int(cfg_run.get("batch_size", 32)),
        "max_bp_len": max_bp_len,
        "kmer": int(cfg_run.get("kmer", 3)),
        # Engineered features
        "stage1_use_engineered_features": bool(cfg_run.get("stage1_use_engineered_features", True)),
        "engineered_dim": int(cfg_run.get("engineered_dim", 288)),
        "engineered_mlp_hidden": int(cfg_run.get("engineered_mlp_hidden", 256)),
        "engineered_mlp_output": int(cfg_run.get("engineered_mlp_output", 128)),
        # Augmentation (effective Stage 1 values with fallbacks)
        "reverse_complement_prob": float(cfg_run.get("stage1_reverse_complement_prob", 0.5)),
        "stage1_use_mixup": bool(get_with_fallback(cfg_run, "stage1_use_mixup", "use_mixup", True)),
        "stage1_mixup_alpha": float(get_with_fallback(cfg_run, "stage1_mixup_alpha", "mixup_alpha", 0.2)),
        # TCN
        "use_tcn": bool(cfg_run.get("use_tcn", True)),
        "tcn_hidden": int(cfg_run.get("tcn_hidden", 256)),
        "tcn_levels": int(cfg_run.get("tcn_levels", 4)),
        "tcn_kernel": int(cfg_run.get("tcn_kernel", 3)),
        "multiscale_channels": int(cfg_run.get("multiscale_channels", 64)),
        "multiscale_kernels": cfg_run.get("multiscale_kernels", [3, 5, 7, 9, 15]),
        "lstm_hidden": int(cfg_run.get("lstm_hidden", 192)),
        # Training
        "epochs": int(epochs),
        "grad_accum_steps": int(cfg_run.get("grad_accum_steps", 4)),
        "grad_clip_norm": float(cfg_run.get("grad_clip_norm", 1.0)),
        "warmup_ratio": float(cfg_run.get("warmup_ratio", 0.15)),
        "early_stopping_patience": int(cfg_run.get("early_stopping_patience", 35)),
        "layer_lr_decay": float(cfg_run.get("layer_lr_decay", 0.8)),
        "stage1_drop_path_rate": float(get_with_fallback(cfg_run, "stage1_drop_path_rate", "drop_path_rate", 0.1)),
        # Loss (effective Stage 1 values with fallbacks)
        "stage1_use_focal_loss": bool(cfg_run.get("stage1_use_focal_loss", True)),
        "stage1_focal_alpha": float(cfg_run.get("stage1_focal_alpha", 0.5)),
        "stage1_focal_gamma": float(cfg_run.get("stage1_focal_gamma", 2.0)),
        "stage1_label_smoothing": float(get_with_fallback(cfg_run, "stage1_label_smoothing", "label_smoothing", 0.05)),
        # SWA (effective Stage 1 values with fallbacks)
        "stage1_use_swa": bool(get_with_fallback(cfg_run, "stage1_use_swa", "use_swa", True)),
        "stage1_swa_start_epoch": int(get_with_fallback(cfg_run, "stage1_swa_start_epoch", "swa_start_epoch", 60)),
        "stage1_swa_lr": float(get_with_fallback(cfg_run, "stage1_swa_lr", "swa_lr", 5e-5)),
        # Model options
        "use_attention_pool": bool(cfg_run.get("use_attention_pool", True)),
        "post_cnn_transformer_layers": int(cfg_run.get("post_cnn_transformer_layers", 3)),
        "fusion_hidden": int(cfg_run.get("fusion_hidden", 256)),
        # Pretrained weights
        "stage1_load_mlm_weights": bool(cfg_run.get("stage1_load_mlm_weights", True)),
        "stage1_mlm_transfer_mode": str(cfg_run.get("stage1_mlm_transfer_mode", "embed_only")),
        "stage1_freeze_layers": int(cfg_run.get("stage1_freeze_layers", 0)),
        # Multi-scale ensemble
        "use_multi_scale": bool(cfg_run.get("use_multi_scale", False)),
        "multi_scale_kmers": cfg_run.get("multi_scale_kmers", [3, 4, 5, 6]),
        "ensemble_method": cfg_run.get("ensemble_method", "soft_voting"),
        "ensemble_learn_weights": bool(cfg_run.get("ensemble_learn_weights", False)),
        "ensemble_weights": cfg_run.get("ensemble_weights"),
        "ensemble_individual_loss_weight": float(cfg_run.get("ensemble_individual_loss_weight", 0.2)),
        # Seed
        "seed": int(cfg_run.get("seed", 1337) or 1337),
        # AMP (Mixed Precision)
        "amp_enabled": bool(cfg_run.get("amp_enabled", True)),
        "amp_dtype": str(cfg_run.get("amp_dtype", "bfloat16")),
        "amp_force_fp32_loss": bool(cfg_run.get("amp_force_fp32_loss", True)),
    }
    
    resolved_config_path = run_dir / "resolved_config.json"
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)
    LOGGER.info("Saved resolved config to %s", resolved_config_path)
    
    seed = int(cfg_run.get("seed", 1337) or 1337)
    set_seed(seed)
    LOGGER.info("Using seed: %d", seed)
    
    # Device selection (optimized for M4)
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    LOGGER.info("Using device: %s", device)

    # AMP (Mixed Precision) setup
    amp_ctx = get_amp_context(device, cfg_run, logger=LOGGER)
    log_amp_status(amp_ctx, logger=LOGGER)

    # GradScaler for CUDA fp16 (not needed for MPS bfloat16)
    scaler = None
    if amp_ctx.enabled and device.type == "cuda" and amp_ctx.dtype == torch.float16:
        scaler = torch.cuda.amp.GradScaler()
        LOGGER.info("Using GradScaler for CUDA fp16")

    if bool(cfg_run.get("use_multi_scale", False)):
        LOGGER.info("Multi-scale co-training enabled")
        return train_multi_scale_ensemble(
            cfg_run=cfg_run,
            device=device,
            amp_ctx=amp_ctx,
            scaler=scaler,
            run_dir=run_dir,
            resolved_config_path=resolved_config_path,
            resolved_config=resolved_config,
            epochs=epochs,
            ablation_max_steps=ablation_max_steps,
            seed=seed,
            config_id=config_id,
            start_time=start_time,
        )

    # Build dataloaders
    dataloaders = build_dataloaders(cfg_run)
    train_loader = dataloaders["stage1"]["train"]
    val_loader = dataloaders["stage1"]["val"]

    batch_size = train_loader.batch_size or int(cfg_run.get("batch_size", 64))
    log_stage1_training_configuration(
        max_bp_len=max_bp_len,
        lr=lr,
        batch_size=batch_size,
        logger=LOGGER,
    )
    
    # Model configuration (enhanced defaults)
    embedding_dim = int(cfg_run.get("embedding_dim", 192))
    transformer_layers = int(cfg_run.get("transformer_layers", 6))
    transformer_heads = int(cfg_run.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg_run.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg_run.get("transformer_dropout", 0.15))
    
    train_dataset = train_loader.dataset
    stage1_vocab = getattr(train_dataset, "vocab", None)
    if stage1_vocab is None:
        raise ValueError("Stage 1 dataset is missing vocabulary")
    
    vocab_size = len(stage1_vocab.itos)
    pad_token_id = stage1_vocab.pad_id
    kmer = stage1_vocab.k
    
    use_attention_pool = bool(cfg_run.get("use_attention_pool", True))
    use_engineered_features = bool(cfg_run.get("stage1_use_engineered_features", True))
    engineered_dim = int(cfg_run.get("engineered_dim", 288))
    
    # TCN configuration
    use_tcn = bool(cfg_run.get("use_tcn", True))
    tcn_hidden = int(cfg_run.get("tcn_hidden", 256))
    tcn_levels = int(cfg_run.get("tcn_levels", 4))
    tcn_kernel = int(cfg_run.get("tcn_kernel", 3))
    multiscale_channels = int(cfg_run.get("multiscale_channels", 64))
    multiscale_kernels_cfg = cfg_run.get("multiscale_kernels", [3, 5, 7, 9, 15])
    multiscale_kernels = tuple(multiscale_kernels_cfg)
    lstm_hidden = int(cfg_run.get("lstm_hidden", 192))
    
    # New architecture parameters: post-CNN transformer, engineered MLP, fusion
    post_cnn_transformer_layers = int(cfg_run.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg_run.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg_run.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg_run.get("fusion_hidden", 256))

    # Compute stage1_drop_path_rate with fallback (needed before model creation)
    stage1_drop_path_rate = get_with_fallback(cfg_run, "stage1_drop_path_rate", "drop_path_rate", 0.1)

    model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
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
        # Stochastic depth rate (effective Stage 1 value)
        drop_path_rate=stage1_drop_path_rate,
        # Positional embedding max length
        max_seq_len=max_bp_len - kmer + 1,
    ).to(device)
    
    # Log architecture configuration
    LOGGER.info("=" * 60)
    LOGGER.info("MODEL ARCHITECTURE V4")
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
    load_mlm_weights = bool(cfg_run.get("stage1_load_mlm_weights", False))
    transfer_mode = str(cfg_run.get("stage1_mlm_transfer_mode", "embed_only"))
    if transfer_mode not in {"none", "embed_only", "embed_plus_adapter"}:
        raise ValueError(f"Invalid stage1_mlm_transfer_mode: {transfer_mode}")

    # Safety guard: ensure Stage 1 kmer matches mlm_kmer when loading MLM weights
    if load_mlm_weights and transfer_mode != "none":
        mlm_kmer = int(cfg_run.get("mlm_kmer", 3))
        if kmer != mlm_kmer:
            raise ValueError(
                f"K-mer mismatch: Stage 1 kmer={kmer} but mlm_kmer={mlm_kmer}. "
                f"When stage1_load_mlm_weights=true, these must match to ensure "
                f"proper weight transfer. Either set kmer={mlm_kmer} or disable "
                f"stage1_load_mlm_weights."
            )
    mlm_encoder_path = cfg_run.get("mlm_encoder_path")
    mlm_encoder_checkpoint = cfg_run.get("mlm_encoder_checkpoint")
    mlm_checkpoint = mlm_encoder_path or mlm_encoder_checkpoint
    should_load_mlm = load_mlm_weights and transfer_mode != "none"
    with torch.no_grad():
        emb_before = model.embedding.weight.detach().float().clone()
        emb_before_norm = emb_before.norm().item()
    load_attempted = False
    if should_load_mlm and mlm_checkpoint:
        checkpoint_path = Path(mlm_checkpoint)
        if not checkpoint_path.is_absolute():
            checkpoint_path = (script_dir / checkpoint_path).resolve()
        if checkpoint_path.exists():
            LOGGER.info(
                "Loading MLM encoder weights from %s (transfer_mode=%s)",
                checkpoint_path,
                transfer_mode,
            )
            load_attempted = True
            model.load_pretrained_weights(
                checkpoint_path, strict=False, transfer_mode=transfer_mode
            )
            if transfer_mode == "embed_plus_adapter":
                LOGGER.info(
                    "Pretrained weights loaded: embedding + %d transformer layers",
                    post_cnn_transformer_layers,
                )
            else:
                LOGGER.info("Pretrained weights loaded: embedding only")
            
            freeze_layers = int(cfg_run.get("stage1_freeze_layers", 0))
            if freeze_layers > 0:
                LOGGER.info("Freezing bottom %d encoder layers", freeze_layers)
                model.encoder.freeze_bottom_layers(freeze_layers)
        else:
            LOGGER.warning("MLM checkpoint not found at %s, training from scratch", checkpoint_path)
    elif load_mlm_weights and transfer_mode == "none":
        LOGGER.info("stage1_mlm_transfer_mode=none; skipping MLM weight load")
    with torch.no_grad():
        emb_after = model.embedding.weight.detach().float()
        emb_after_norm = emb_after.norm().item()
        diff_norm = (emb_after - emb_before).norm().item()
    LOGGER.info(
        "Load sanity: embedding norm before=%.6f after=%.6f diff_norm=%.6f (load_attempted=%s)",
        emb_before_norm,
        emb_after_norm,
        diff_norm,
        load_attempted,
    )
    if load_attempted and diff_norm < 1e-6:
        LOGGER.warning("Pretrain load had no effect (likely key mismatch or already loaded).")
    
    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    LOGGER.info("Model: %.2fM total params, %.2fM trainable",
                total_params / 1e6, trainable_params / 1e6)

    # Update resolved_config with device and model info, then re-save
    resolved_config["system_info"]["device"] = str(device)
    resolved_config["system_info"]["stage1_total_params"] = total_params
    resolved_config["system_info"]["stage1_trainable_params"] = trainable_params
    with resolved_config_path.open("w", encoding="utf-8") as f:
        json.dump(resolved_config, f, indent=2)

    # Optimizer with layer-wise learning rate decay
    weight_decay = float(cfg_run.get("weight_decay", 0.02))
    layer_lr_decay = float(cfg_run.get("layer_lr_decay", 0.8))
    optimizer_foreach = bool(cfg_run.get("optimizer_foreach", True))

    LOGGER.info("Optimizer: AdamW foreach=%s", optimizer_foreach)

    # Use LLRD if we loaded pretrained weights
    if should_load_mlm and mlm_checkpoint:
        LOGGER.info("Using layer-wise LR decay (factor=%.2f) for pretrained model", layer_lr_decay)
        param_groups = get_layer_wise_lr_groups(
            model,
            base_lr=lr,
            lr_decay=layer_lr_decay,
            weight_decay=weight_decay,
        )
        optimizer = torch.optim.AdamW(param_groups, betas=(0.9, 0.999), foreach=optimizer_foreach)
        LOGGER.info("Created %d parameter groups with LLRD", len(param_groups))
    else:
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay,
            betas=(0.9, 0.999),
            foreach=optimizer_foreach,
        )
    
    # Scheduler with warmup
    warmup_ratio = float(cfg_run.get("warmup_ratio", 0.15))
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

    # =========================================================================
    # Stage 1 effective hyperparameters (with stage1_ overrides and fallbacks)
    # =========================================================================
    stage1_label_smoothing = get_with_fallback(cfg_run, "stage1_label_smoothing", "label_smoothing", 0.05)
    stage1_use_mixup = get_with_fallback(cfg_run, "stage1_use_mixup", "use_mixup", True)
    stage1_mixup_alpha = get_with_fallback(cfg_run, "stage1_mixup_alpha", "mixup_alpha", 0.2)
    stage1_use_swa = get_with_fallback(cfg_run, "stage1_use_swa", "use_swa", True)
    stage1_swa_start_epoch = get_with_fallback(cfg_run, "stage1_swa_start_epoch", "swa_start_epoch", 60)
    stage1_swa_lr = get_with_fallback(cfg_run, "stage1_swa_lr", "swa_lr", 5e-5)
    stage1_drop_path_rate = get_with_fallback(cfg_run, "stage1_drop_path_rate", "drop_path_rate", 0.1)
    stage1_use_focal_loss = bool(cfg_run.get("stage1_use_focal_loss", True))

    # Log effective Stage 1 hyperparameters
    LOGGER.info("=" * 60)
    LOGGER.info("EFFECTIVE STAGE 1 HYPERPARAMETERS")
    LOGGER.info("=" * 60)
    LOGGER.info("  stage1_label_smoothing: %.3f", stage1_label_smoothing)
    LOGGER.info("  stage1_use_mixup: %s", stage1_use_mixup)
    LOGGER.info("  stage1_mixup_alpha: %.3f", stage1_mixup_alpha)
    LOGGER.info("  stage1_use_swa: %s", stage1_use_swa)
    LOGGER.info("  stage1_swa_start_epoch: %d", stage1_swa_start_epoch)
    LOGGER.info("  stage1_swa_lr: %.2e", stage1_swa_lr)
    LOGGER.info("  stage1_drop_path_rate: %.3f", stage1_drop_path_rate)
    LOGGER.info("  stage1_use_focal_loss: %s", stage1_use_focal_loss)
    LOGGER.info("=" * 60)

    # SWA setup for breaking plateaus (using effective Stage 1 values)
    swa_model = None
    swa_scheduler = None

    if stage1_use_swa:
        swa_model = AveragedModel(model)
        swa_scheduler = SWALR(optimizer, swa_lr=stage1_swa_lr, anneal_epochs=10)
        LOGGER.info("SWA enabled: starts at epoch %d with LR %.2e", stage1_swa_start_epoch, stage1_swa_lr)

    # Loss function with label smoothing (using effective Stage 1 values)
    focal_alpha = float(cfg_run.get("stage1_focal_alpha", 0.5))
    focal_gamma = float(cfg_run.get("stage1_focal_gamma", 2.0))

    if stage1_use_focal_loss:
        criterion = FocalLossWithSmoothing(
            alpha=focal_alpha,
            gamma=focal_gamma,
            smoothing=stage1_label_smoothing,
        )
        LOGGER.info("Using focal loss with smoothing=%.2f", stage1_label_smoothing)
    else:
        criterion = LabelSmoothingBCE(smoothing=stage1_label_smoothing)
        LOGGER.info("Using BCE with smoothing=%.2f", stage1_label_smoothing)
    
    # =========================================================================
    # OVERFIT DEBUG MODE: Train on 32 samples for 200 steps
    # =========================================================================
    if bool(cfg_run.get("overfit_debug", False)):
        LOGGER.info("Running OVERFIT DEBUG mode")
        # Use simple BCEWithLogitsLoss for overfit test (no label smoothing or focal loss)
        simple_criterion = nn.BCEWithLogitsLoss()
        run_overfit_debug_stage1(
            model=model,
            train_loader=train_loader,
            criterion=simple_criterion,
            device=device,
            num_steps=200,
            lr=1e-3,
        )
        LOGGER.info("Overfit debug complete. Exiting.")
        return {
            "overfit_debug": True,
            "stage1_acc": float("nan"),
            "stage1_auc": float("nan"),
            "stage1_mcc": float("nan"),
            "best_ckpt_path": None,
            "stage1_seconds": float(time.perf_counter() - start_time),
            "params": int(total_params),
            "best_val_acc": float("nan"),
            "best_val_auc": float("nan"),
            "best_val_mcc": float("nan"),
            "best_epoch": 0,
        }

    # =========================================================================
    # PROFILING MODE: Run profiler for ~50 steps
    # =========================================================================
    if bool(cfg_run.get("profile", False)):
        LOGGER.info("Running PROFILING mode")
        grad_accum_steps_prof = int(cfg_run.get("grad_accum_steps", 4))
        grad_clip_norm_prof = float(cfg_run.get("grad_clip_norm", 1.0))

        profile_metrics = run_profiling_epoch(
            loader=train_loader,
            model=model,
            criterion=criterion,
            device=device,
            optimizer=optimizer,
            run_dir=run_dir,
            num_profile_steps=50,
            grad_accum_steps=grad_accum_steps_prof,
            grad_clip_norm=grad_clip_norm_prof,
        )

        # Save profile metrics to run_dir
        profile_report_path = run_dir / "profile_report.json"
        with profile_report_path.open("w", encoding="utf-8") as f:
            json.dump(profile_metrics, f, indent=2)
        LOGGER.info("Saved profile report to %s", profile_report_path)
        LOGGER.info("Profiling complete. Exiting.")

        return {
            "profile": True,
            "profile_metrics": profile_metrics,
            "stage1_acc": float("nan"),
            "stage1_auc": float("nan"),
            "stage1_mcc": float("nan"),
            "best_ckpt_path": None,
            "stage1_seconds": float(time.perf_counter() - start_time),
            "params": int(total_params),
            "best_val_acc": float("nan"),
            "best_val_auc": float("nan"),
            "best_val_mcc": float("nan"),
            "best_epoch": 0,
        }

    # =========================================================================
    # KNOWLEDGE DISTILLATION SETUP
    # =========================================================================
    from dataset import KmerVocabulary

    distill_cfg = cfg_run.get("distill", {})
    distill_enabled = bool(distill_cfg.get("enabled", False))
    teacher_model = None
    teacher_vocab = None
    distill_criterion = None

    if distill_enabled:
        distill_alpha = float(distill_cfg.get("alpha", 0.5))
        distill_temperature = float(distill_cfg.get("temperature", 2.0))
        teacher_ckpt = distill_cfg.get("teacher_ckpt", "../artifacts/checkpoints/stage1_k6.pt")
        teacher_kmer = int(distill_cfg.get("teacher_kmer", 6))

        # Resolve teacher checkpoint path
        teacher_ckpt_path = Path(teacher_ckpt)
        if not teacher_ckpt_path.is_absolute():
            teacher_ckpt_path = (script_dir / teacher_ckpt_path).resolve()

        # Compute student checkpoint path for guardrail check
        artifacts_root = (script_dir / "../artifacts").resolve()
        checkpoint_name = cfg_run.get("stage1_checkpoint_name") or f"stage1_k{kmer}.pt"
        student_ckpt_path = artifacts_root / "checkpoints" / checkpoint_name

        # Apply distillation safety guardrails
        guardrail_ok, guardrail_reason = apply_distill_guardrails(
            distill_enabled=distill_enabled,
            teacher_kmer=teacher_kmer,
            student_kmer=kmer,
            teacher_ckpt_path=teacher_ckpt_path,
            student_ckpt_path=student_ckpt_path,
        )
        if not guardrail_ok:
            LOGGER.warning("Distillation guardrail triggered: %s", guardrail_reason)
            distill_enabled = False

        if not distill_enabled:
            pass  # Skip distillation setup, teacher_model/teacher_vocab/distill_criterion stay None
        elif not teacher_ckpt_path.exists():
            LOGGER.error("Teacher checkpoint not found: %s", teacher_ckpt_path)
            LOGGER.error("Distillation requires a trained teacher model. Disabling distillation.")
            distill_enabled = False
        else:
            LOGGER.info("=" * 60)
            LOGGER.info("KNOWLEDGE DISTILLATION ENABLED")
            LOGGER.info("=" * 60)
            LOGGER.info("Teacher checkpoint: %s", teacher_ckpt_path)
            LOGGER.info("Teacher k-mer: %d, Student k-mer: %d", teacher_kmer, kmer)
            LOGGER.info("Alpha: %.2f, Temperature: %.2f", distill_alpha, distill_temperature)
            LOGGER.info("=" * 60)

            # Update distill config with teacher_kmer for load_teacher_model
            distill_cfg["teacher_kmer"] = teacher_kmer
            cfg_run["distill"] = distill_cfg

            # Load teacher model
            teacher_model = load_teacher_model(cfg_run, teacher_ckpt_path, device)

            # Get teacher vocab
            vocab_cache_dir = Path(cfg_run.get("vocab_cache_dir", "../artifacts/vocabs"))
            if not vocab_cache_dir.is_absolute():
                vocab_cache_dir = (script_dir / vocab_cache_dir).resolve()
            teacher_vocab_path = vocab_cache_dir / f"k{teacher_kmer}_vocab.json"
            teacher_vocab = KmerVocabulary.load(teacher_vocab_path)

            # Create distillation criterion with base criterion (focal loss or BCE)
            distill_criterion = DistillationLoss(
                alpha=distill_alpha,
                temperature=distill_temperature,
                base_criterion=criterion,  # Use the same loss function for CE part
            )

            # Dry run verification
            dry_run_ok = verify_distillation_dry_run(
                student=model,
                teacher=teacher_model,
                train_loader=train_loader,
                device=device,
                teacher_vocab=teacher_vocab,
                student_vocab=stage1_vocab,
            )

            if not dry_run_ok:
                LOGGER.error("Distillation dry-run failed! Falling back to regular training.")
                distill_enabled = False
                teacher_model = None
                teacher_vocab = None
                distill_criterion = None
            else:
                LOGGER.info("Distillation setup complete. Training will use CE + KL loss.")

    # Training settings
    grad_accum_steps = int(cfg_run.get("grad_accum_steps", 2))
    grad_clip_norm = float(cfg_run.get("grad_clip_norm", 1.0))
    param_norm_warn = float(cfg_run.get("param_norm_warn", 150.0))
    param_norm_cap = float(cfg_run.get("param_norm_cap", 200.0))
    # Note: stage1_use_mixup and stage1_mixup_alpha are computed earlier with fallbacks
    patience = int(cfg_run.get("early_stopping_patience", 15))

    global_optimizer_steps = 0

    class_balance = compute_class_balance(train_loader)
    LOGGER.info(
        "Class balance: %.2f%% positive (%d/%d)",
        class_balance["positive_ratio"] * 100.0,
        class_balance["positive"],
        class_balance["total"],
    )
    
    # Paths
    artifacts_root = (script_dir / "../artifacts").resolve()
    checkpoint_name = cfg_run.get("stage1_checkpoint_name") or f"stage1_k{kmer}.pt"
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
        # Check if we should start SWA (using effective Stage 1 values)
        in_swa_phase = stage1_use_swa and epoch >= stage1_swa_start_epoch
        if in_swa_phase and not swa_started:
            LOGGER.info("=" * 40)
            LOGGER.info("Starting SWA phase at epoch %d", epoch)
            LOGGER.info("=" * 40)
            swa_started = True
        
        current_lr = optimizer.param_groups[0]["lr"]
        phase_str = "[SWA]" if in_swa_phase else ""
        LOGGER.info("Epoch %d/%d %s (LR: %.2e)", epoch, epochs, phase_str, current_lr)
        
        remaining_steps = None if ablation_max_steps is None else max(ablation_max_steps - global_optimizer_steps, 0)

        # Use distillation training if enabled
        if distill_enabled and teacher_model is not None:
            train_metrics = run_distillation_epoch(
                student_loader=train_loader,
                student_model=model,
                teacher_model=teacher_model,
                distill_criterion=distill_criterion,
                device=device,
                student_vocab=stage1_vocab,
                teacher_vocab=teacher_vocab,
                optimizer=optimizer,
                scheduler=None if in_swa_phase else scheduler,
                desc="Train-Distill",
                grad_accum_steps=grad_accum_steps,
                max_optimizer_steps=remaining_steps,
                grad_clip_norm=grad_clip_norm,
                amp_ctx=amp_ctx,
                scaler=scaler,
                param_norm_warn=param_norm_warn,
                param_norm_cap=param_norm_cap,
            )
        else:
            train_metrics = run_epoch(
                train_loader,
                model,
                criterion,
                device,
                optimizer=optimizer,
                scheduler=None if in_swa_phase else scheduler,  # Use SWA scheduler in SWA phase
                desc="Train",
                grad_accum_steps=grad_accum_steps,
                max_optimizer_steps=remaining_steps,
                grad_clip_norm=grad_clip_norm,
                param_norm_warn=param_norm_warn,
                param_norm_cap=param_norm_cap,
                use_mixup=stage1_use_mixup and not in_swa_phase,  # Disable mixup in SWA phase
                mixup_alpha=stage1_mixup_alpha,
                amp_ctx=amp_ctx,
                scaler=scaler,
            )
        global_optimizer_steps += int(train_metrics.get("optimizer_steps", 0))
        reached_max_steps = ablation_max_steps is not None and global_optimizer_steps >= ablation_max_steps
        
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
            amp_ctx=amp_ctx,
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
                    amp_ctx=amp_ctx,
                )
            except Exception as e:
                LOGGER.warning("SWA evaluation failed: %s", e)
                swa_val_metrics = None
        
        # Log training metrics
        if distill_enabled and "ce_loss" in train_metrics:
            LOGGER.info(
                "Train - Loss: %.4f (CE: %.4f, KL: %.4f) | Acc: %.4f | F1: %.4f | MCC: %.4f",
                train_metrics["loss"],
                train_metrics.get("ce_loss", 0.0),
                train_metrics.get("kl_loss", 0.0),
                train_metrics["accuracy"],
                train_metrics["f1"],
                train_metrics["mcc"],
            )
        else:
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
            "Val   - Loss: %.4f | Acc@0.5: %.4f | Acc@best: %.4f (thr=%.2f) | MCC@best: %.4f",
            val_metrics["loss"],
            val_metrics["accuracy"],
            val_metrics.get("acc_best", val_metrics["accuracy"]),
            val_metrics.get("best_threshold", 0.5),
            val_metrics.get("mcc_best", val_metrics["mcc"]),
        )
        
        if swa_val_metrics:
            LOGGER.info(
                "SWA   - Loss: %.4f | Acc@0.5: %.4f | Acc@best: %.4f (thr=%.2f) | MCC@best: %.4f",
                swa_val_metrics["loss"],
                swa_val_metrics["accuracy"],
                swa_val_metrics.get("acc_best", swa_val_metrics["accuracy"]),
                swa_val_metrics.get("best_threshold", 0.5),
                swa_val_metrics.get("mcc_best", swa_val_metrics["mcc"]),
            )
        
        # Use SWA metrics if better during SWA phase (compare using best-threshold MCC)
        eval_metrics = val_metrics
        swa_mcc_best = swa_val_metrics.get("mcc_best", swa_val_metrics["mcc"]) if swa_val_metrics else -1
        val_mcc_best = val_metrics.get("mcc_best", val_metrics["mcc"])
        if swa_val_metrics and swa_mcc_best > val_mcc_best:
            eval_metrics = swa_val_metrics
            LOGGER.info("(Using SWA metrics for checkpointing)")

        # Check for improvement using best-threshold metrics (MCC@best primary, Acc@best tiebreaker)
        eval_mcc_best = eval_metrics.get("mcc_best", eval_metrics["mcc"])
        eval_acc_best = eval_metrics.get("acc_best", eval_metrics["accuracy"])
        improved = False
        if eval_mcc_best > best_val_mcc + 1e-4:
            improved = True
        elif abs(eval_mcc_best - best_val_mcc) < 1e-4 and eval_acc_best > best_val_acc:
            improved = True

        if improved:
            patience_counter = 0
            best_val_mcc = eval_mcc_best
            best_val_acc = eval_acc_best
            best_epoch = epoch
            best_train_metrics = train_metrics
            best_val_metrics = eval_metrics

            # Save either SWA model or regular model
            if swa_val_metrics and swa_mcc_best > val_mcc_best:
                save_checkpoint(checkpoint_path, swa_model.module, optimizer, scheduler, epoch, eval_metrics)
            else:
                save_checkpoint(checkpoint_path, model, optimizer, scheduler, epoch, val_metrics)
            LOGGER.info(
                "★ New best! MCC@best: %.4f, Acc@best: %.4f (%.2f%%), thr=%.2f",
                best_val_mcc, best_val_acc, best_val_acc * 100,
                eval_metrics.get("best_threshold", 0.5),
            )
        else:
            patience_counter += 1
            # In SWA phase, be more patient as averaging takes time
            effective_patience = patience + 10 if in_swa_phase else patience
            LOGGER.info("No improvement. Patience %d/%d", patience_counter, effective_patience)
            if patience_counter >= effective_patience:
                LOGGER.info("Early stopping triggered at epoch %d", epoch)
                break

        if reached_max_steps:
            LOGGER.info(
                "Reached ablation_max_steps_stage1=%d (optimizer steps); stopping after validation",
                ablation_max_steps,
            )
            break
    
    # Final report
    if best_val_metrics is None:
        best_train_metrics = train_metrics
        best_val_metrics = val_metrics
        best_epoch = epochs
    
    report = {
        "config": config_id,
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
            "stage1_label_smoothing": stage1_label_smoothing,
            "stage1_mixup_alpha": stage1_mixup_alpha,
            "stage1_use_swa": stage1_use_swa,
            "stage1_swa_start_epoch": stage1_swa_start_epoch if stage1_use_swa else None,
            "stage1_swa_lr": stage1_swa_lr if stage1_use_swa else None,
            "stage1_drop_path_rate": stage1_drop_path_rate,
            "ablation_max_steps_stage1": ablation_max_steps,
            "ablation_disable_wandb": bool(cfg_run.get("ablation_disable_wandb", False)),
        },
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    LOGGER.info("Saved training report to %s", report_path)

    LOGGER.info("=" * 60)
    LOGGER.info("TRAINING COMPLETE")
    LOGGER.info("Best Epoch: %d%s", best_epoch, " (SWA)" if swa_started and best_epoch >= stage1_swa_start_epoch else "")
    LOGGER.info("Best Val Acc@best: %.4f (%.2f%%)", best_val_acc, best_val_acc * 100)
    LOGGER.info("Best Val MCC@best: %.4f", best_val_mcc)
    if best_val_metrics:
        LOGGER.info("Best threshold: %.2f", best_val_metrics.get("best_threshold", 0.5))
        LOGGER.info("Acc@0.5: %.4f (%.2f%%)", best_val_metrics.get("accuracy", 0), best_val_metrics.get("accuracy", 0) * 100)
    if stage1_use_swa:
        LOGGER.info("SWA was %s", "active" if swa_started else "not triggered (training ended early)")
    LOGGER.info("=" * 60)

    stage1_seconds = float(time.perf_counter() - start_time)
    best_val_auc = float(best_val_metrics.get("roc_auc", 0.5)) if best_val_metrics else 0.5
    return {
        "stage1_acc": float(best_val_acc),
        "stage1_auc": best_val_auc,
        "stage1_mcc": float(best_val_mcc),
        "best_val_acc": float(best_val_acc),
        "best_val_auc": best_val_auc,
        "best_val_mcc": float(best_val_mcc),
        "best_epoch": int(best_epoch),
        "best_ckpt_path": str(checkpoint_path),
        "stage1_seconds": stage1_seconds,
        "params": int(total_params),
    }


if __name__ == "__main__":
    main()
