# Gene Whisperer Project Memory

## Project Overview
Gene Whisperer is a DNA promoter classification model using a transformer-only architecture. The model classifies DNA sequences as promoters or non-promoters (Stage 1) and strong vs weak promoters (Stage 2).

## Architecture
- **Simplified Architecture v5.0**: Embedding → 12-Layer Transformer → Attention Pool → Classifier
- Removed CNN/TCN components after ablation studies showed they hurt performance
- Achieves 84.1% accuracy with 23.9M parameters (vs 83.3% with 38.4M in legacy architecture)

## Recent Changes

### 2026-01-19: Early Stopping with Best Model Restoration
Implemented a robust `EarlyStopping` class in train_stage1.py that saves and restores best model weights:

**New Features:**
- `EarlyStopping` class that tracks best model state_dict in memory
- Automatically restores best weights when early stopping triggers
- Configurable via config.yaml parameters

**Config Changes:**
| Parameter | Old Value | New Value | Purpose |
|-----------|-----------|-----------|---------|
| `early_stopping_patience` | 15 | 10 | More aggressive stopping |
| `early_stopping_min_delta` | N/A | 0.001 | Minimum improvement threshold |
| `early_stopping_restore_best` | N/A | true | Always restore best weights |

**How it works:**
1. After each validation epoch, `early_stopping(mcc_score, epoch, model)` is called
2. If MCC improves by `min_delta`, the model's state_dict is saved in memory
3. If no improvement for `patience` epochs, training stops and best weights are restored
4. SWA phase extends patience by +10 epochs for averaging to stabilize

### 2026-01-19: Stronger Regularization to Combat Overfitting
Training curves showed severe overfitting (train acc 97% vs val acc 80%). Updated config.yaml with stronger regularization:

| Parameter | Old Value | New Value | Purpose |
|-----------|-----------|-----------|---------|
| `transformer_dropout` | 0.12 | 0.20 | More dropout in transformer layers |
| `classifier_dropout` | 0.15 | 0.25 | More dropout before classifier |
| `weight_decay` | 0.01 | 0.05 | Stronger L2 regularization |
| `label_smoothing` | 0.0 | 0.1 | Prevents overconfident predictions |
| `mixup_alpha` | 0.1 | 0.2 | Stronger data augmentation |

**Expected Outcome**: Train/val accuracy gap should reduce from ~17% to <10%. Target: train acc 85%, val acc 82%.

## Key Files
- `Gene Whisperer/training/config.yaml` - Main configuration file
- `Gene Whisperer/training/train_stage1.py` - Stage 1 training script
- `Gene Whisperer/training/model.py` - Model definitions (simplified architecture)
- `Gene Whisperer/training/model_legacy.py` - Legacy CNN/TCN components (deprecated)

## Training Tips
1. Label smoothing is already implemented in `LabelSmoothingBCE` class in train_stage1.py
2. Config values are read via `get_with_fallback()` function which checks stage-specific keys first
3. The model uses mixup augmentation when `use_mixup: true`
4. SWA (Stochastic Weight Averaging) is available for breaking plateaus

### 2026-01-19: Gradient Monitoring for Training Stability
Added comprehensive gradient monitoring to detect training instabilities early.

**New Function: `compute_gradient_stats(model)`**
Returns dictionary with:
- `grad_norm`: L2 norm of all gradients
- `max_grad`: Maximum absolute gradient value
- `min_grad`: Minimum absolute gradient value
- `nan_grads`: Count of parameters with NaN/Inf gradients
- `param_count`: Number of parameters with valid gradients

**Integration Points:**
1. Main training loop - logs every 50 optimizer steps
2. Knowledge distillation loop - logs every 50 optimizer steps
3. Overfit test loop - logs every 20 steps with loss/accuracy

**Logging Behavior:**
- INFO level: `Step X: grad_norm=Y, max=Z, nan_count=N`
- WARNING: Triggered when `grad_norm > 100` (exploding gradients)
- ERROR: Triggered when NaN/Inf gradients are detected

### 2026-01-19: Regularization Experiment Results
Ran `run_regularization_experiment.py` to find optimal regularization settings.

**Results:**
| Name | Val Acc | Val MCC | Best Epoch | Trans. Dropout | Class. Dropout |
|------|---------|---------|------------|----------------|----------------|
| baseline | 85.88% | 0.718 | 15 | 0.12 | 0.15 |
| **medium_reg** | **86.55%** | **0.731** | 20 | 0.20 | 0.20 |
| strong_reg | 85.29% | 0.706 | 28 | 0.25 | 0.30 |

**Winner: medium_reg** - Best balance of regularization. Strong_reg was too aggressive.

**Updated config.yaml with optimal settings:**
| Parameter | Value |
|-----------|-------|
| transformer_dropout | 0.20 |
| classifier_dropout | 0.20 |
| weight_decay | 0.05 |
| label_smoothing | 0.1 |

## Future Considerations
- If regularization is still insufficient, consider:
  - Adding feature dropout before classifier (FeatureDropout layer)
  - Increasing drop_path_rate (currently 0.1)
  - Reducing model capacity (fewer transformer layers)

### 2026-01-20: Optimized Ensemble Weights
Implemented optimized soft-voting weights for multi-kmer Stage 1 ensembles.

**What changed:**
- Added `ensemble_weights.py` with grid/random search optimization, weight persistence, and load helpers.
- `evaluate_ensemble.py` now applies sigmoid to logits, compares baseline vs optimized weights, and writes weights to `artifacts/ensemble_weights.json`.
- `ensemble_infer.py` can load optimized weights via config or CLI (`--optimized-weights`).
- Added `ensemble` config block (`optimize_weights`, `optimization_metric`, `min_weight`, `weights_path`).
- Updated config-related tests to reflect current regularization settings and added tests for weight optimization utilities.

### 2026-01-20: Test-Time Augmentation (TTA) with Reverse Complement
Implemented TTA to average predictions over forward and reverse complement DNA strands.

**Why TTA?**
DNA is double-stranded, so a sequence and its reverse complement encode the same biological information. By averaging predictions from both orientations, we improve model robustness. Expected improvement: +0.3-0.5% accuracy.

**New Files:**
- `Gene Whisperer/training/tta.py` - Core TTA implementation
  - `get_reverse_complement_tokens()` - Converts tokenized sequence to its reverse complement
  - `predict_with_tta()` - Single prediction with TTA
  - `predict_batch_with_tta()` - Batch prediction returning forward, RC, and averaged probs
  - `TTAWrapper` - Convenient wrapper class for TTA-enabled models

- `Gene Whisperer/training/tests/test_tta.py` - Comprehensive test suite (18 tests)

**Modified Files:**
- `Gene Whisperer/training/ensemble_infer.py`
  - Added `--tta` and `--tta-aggregation` CLI flags
  - Reads TTA settings from config `inference.use_tta`
  - Stores vocab for each model to enable TTA tokenization
  - TTA details included in `--per_model` output (forward/RC probabilities)

- `Gene Whisperer/training/config.yaml`
  - Added `inference` section with TTA settings:
    ```yaml
    inference:
      use_tta: true
      tta_aggregation: "mean"  # or "geometric_mean"
    ```

**Usage:**
```bash
# Via CLI flag
python ensemble_infer.py --sequence "ATGCATGC..." --tta

# Via config (inference.use_tta: true)
python ensemble_infer.py --sequence "ATGCATGC..."

# With geometric mean aggregation
python ensemble_infer.py --sequence "ATGCATGC..." --tta --tta-aggregation geometric_mean

# See per-model TTA details
python ensemble_infer.py --sequence "ATGCATGC..." --tta --per_model
```

**Output with TTA:**
```json
{
  "probability": 0.85,
  "label": "promoter",
  "threshold": 0.5,
  "tta_enabled": true,
  "tta_aggregation": "mean",
  "per_model": [
    {"k": 6, "probability": 0.85, "forward": 0.87, "reverse_complement": 0.83}
  ]
}
```

**Key Properties:**
1. Double reverse complement returns original sequence (verified by tests)
2. Engineered features are recomputed for RC sequence for accuracy
3. TTA can be disabled at runtime with config or by omitting `--tta` flag

### 2026-01-20: SwiGLU Activation and RMSNorm
Implemented SwiGLU feed-forward network and RMSNorm as modern transformer components (used in LLaMA, PaLM, Mistral).

**Why SwiGLU?**
SwiGLU provides better gradient flow and representation capacity than standard GELU FFN. The formula is:
```
SwiGLU(x) = (Swish(W1*x) * W3*x) * W2
```
Uses SiLU (Swish) activation with a gating mechanism for improved expressiveness.

**Why RMSNorm?**
RMSNorm is simpler and faster than LayerNorm:
- No mean centering (only RMS normalization)
- No bias term
- Provides similar training stability with lower computational cost

**New Files:**
- `Gene Whisperer/training/tests/test_swiglu_rmsnorm.py` - Comprehensive test suite (24 tests)

**Modified Files:**
- `Gene Whisperer/training/model.py`
  - Added `RMSNorm` class for RMS layer normalization
  - Added `SwiGLU` class with the LLaMA-style implementation
  - Updated `PreNormTransformerLayer` with `ffn_type` and `norm_type` parameters
  - Updated `DNAEncoder` and all model classes to pass through new parameters
  - Updated `create_model_variant()` factory function

- `Gene Whisperer/training/config.yaml`
  - Added `ffn_type: "swiglu"` (options: "gelu", "glu", "swiglu")
  - Added `ffn_mult: 2.67` (SwiGLU expansion ratio)
  - Added `norm_type: "rmsnorm"` (options: "layernorm", "rmsnorm")

**Config Options:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `ffn_type` | "swiglu" | FFN type: "gelu" (standard), "glu" (GLUFFN), "swiglu" |
| `ffn_mult` | 2.67 | FFN expansion multiplier (for SwiGLU, lower due to gating) |
| `norm_type` | "rmsnorm" | Normalization: "layernorm" or "rmsnorm" |

**Backward Compatibility:**
- Legacy `use_glu_ffn` parameter still works (maps to `ffn_type="glu"`)
- Default `norm_type="layernorm"` preserves existing behavior if not specified
- All existing checkpoints remain compatible

**Expected Improvement:** +0.5-1% accuracy from better gradient flow

### 2026-01-19: Rotary Position Embedding (RoPE)
Implemented Rotary Position Embedding (RoPE) as an alternative to absolute position embeddings.

**Why RoPE?**
RoPE from RoFormer/LLaMA rotates query and key vectors based on their positions instead of adding position embeddings. This captures relative positions better than absolute embeddings, which is important for DNA motif detection where the relative distance between bases matters more than their absolute position.

**New Files:**
- `Gene Whisperer/training/rope.py` - Core RoPE implementation
  - `RotaryEmbedding` class - Precomputes and caches cos/sin for position rotation
  - `rotate_half()` - Helper function to rotate half the hidden dimensions
  - `apply_rotary_emb()` - Applies rotation to query and key tensors
  - `RoPEAttention` - Standalone attention module with RoPE (for reference)

**Modified Files:**
- `Gene Whisperer/training/model.py`
  - Added `use_rope` and `rope_base` parameters to `DNAEncoder`
  - Added RoPE support to `PreNormTransformerLayer`
  - Updated all model classes (`GeneWhispererStage1`, `GeneWhispererTransformerOnly`, `GeneWhispererCombined`, `GeneWhispererStage2`) to pass RoPE parameters
  - Updated `create_model_variant()` factory function

- `Gene Whisperer/training/config.yaml`
  - Added `use_rope: true` (enables RoPE)
  - Added `rope_base: 10000.0` (frequency base)
  - Set `use_relative_position_bias: false` (disabled when using RoPE)

**Config Options:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_rope` | true | Enable Rotary Position Embedding |
| `rope_base` | 10000.0 | Base for frequency computation |

**How it Works:**
1. When `use_rope=True`, position embeddings are NOT added to input embeddings
2. Instead, RoPE rotates Q and K vectors in each attention layer based on position
3. The rotation encodes relative position information directly in the attention computation
4. RoPE is applied before computing attention scores: `attn = softmax((RoPE(Q) @ RoPE(K)^T) / sqrt(d))`

**Expected Improvement:** +1-2% accuracy due to better relative position understanding

**Compatibility:**
- RoPE is mutually exclusive with `use_relative_position_bias` - when RoPE is enabled, relative position bias is ignored
- Checkpoints are backward compatible (position embeddings are still created but not used when RoPE is enabled)
