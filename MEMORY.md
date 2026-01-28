# Gene Whisperer Project Memory (Condensed)

_Last summarized: 2026-01-28_

## Overview
- DNA promoter classification with a transformer-only architecture.
- Stage 1: promoter vs non-promoter. Stage 2: strong vs weak promoters.
- Primary model variants: v5 (simplified transformer-only) and v6 (modern transformer with RoPE/SwiGLU/RMSNorm/positional bias/GQA).

## Architecture Highlights
- Embedding -> N transformer blocks -> attention pooling -> classifier.
- Optional engineered features concatenated pre-classifier in combined variants.
- V6: ModernTransformerBlock with RoPE attention, optional GQA, SwiGLU FFN, RMSNorm, positional motif bias (first layer).
- RoPE and positional motif bias are compatible; relative position bias is disabled when RoPE is on.

## Key Files
- `Gene Whisperer/training/config.yaml`
- `Gene Whisperer/training/model.py` (v5 + common blocks)
- `Gene Whisperer/training/model_v6.py` (v6)
- `Gene Whisperer/training/train_stage1.py` (Stage 1 training)
- `Gene Whisperer/training/pretrain_mlm.py` (MLM pretraining)
- `Gene Whisperer/training/tta.py` (TTA/RC)
- `Gene Whisperer/training/ensemble_infer.py`, `evaluate_ensemble.py`, `ensemble_weights.py`

## Current Defaults / Recommended Config (as of 2026-01-28)
- Regularization: transformer_dropout=0.20, classifier_dropout=0.20, weight_decay=0.05, label_smoothing=0.1, mixup_alpha=0.2.
- Early stopping: patience=10, min_delta=0.001, restore_best=true (SWA extends patience +10).
- Model arch: use_rope=true, rope_base=10000.0, ffn_type="swiglu", ffn_mult=2.67, norm_type="rmsnorm".
- Positional motif bias enabled with TATA, -10, TSS regions.

## Inference/Robustness
- TTA: average forward + reverse complement probabilities; configurable aggregation (mean or geometric mean).
- Bidirectional consistency loss for training (alpha=0.3) to enforce strand agreement.

## Ensemble
- Soft-voting weights optimized and persisted in `artifacts/ensemble_weights.json`.
- `ensemble_infer.py` can load optimized weights and apply TTA.

## Data Pipeline & Performance
- BPE tokenizer uses greedy longest-match (fast) rather than merge-iteration.
- Dataset pre-caches tokens and engineered features for forward and reverse complement.
- DataLoader prefetch_factor increased to 4 to improve GPU utilization.

## MLM Pretraining (BPE) — Critical Fixes Summary
- BPE vocab: 4096 tokens; default `mlm_max_token_len=64` for 234bp windows.
- Fixed padding bug: tokenization now pads to `max_token_len`, not `window_bp`.
- Fixed 80/10/10 masking: random replacement excludes special tokens (IDs 0-4).
- Masking config now passes `max_span_len`, `span_distribution`, `mean_span`, `exclude_special_from_labels`.
- Removed label_smoothing in MLM loss due to `-inf` logits for special tokens.
- Disabled CLIP-style projection for MLM (`mlm_use_clip_projection=false`).
- Modern arch params (RoPE/SwiGLU/RMSNorm/ffn_mult) are now passed into MLM encoders.

## Known Gotchas
- RoPE and relative position bias are mutually exclusive (bias is ignored when RoPE is on).
- For BPE MLM, `mlm_max_span_len=1` is important to avoid masking too much context.

## Recent Decisions / Notes
- 2026-01-28: Biopython is helpful for parsing/annotation pipelines, but not required for core training/inference. Prefer keeping core loops dependency-light and performance-focused.
- 2026-01-28: Fixed BPE MLM padding in non-streaming paths to use `mlm_max_token_len` (token count) instead of `mlm_window_size` (bp), and corrected dry-run sample tokenization arg; added unit test for padding length.
