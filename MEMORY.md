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

### 2026-01-20: Bidirectional Consistency Training
Implemented bidirectional consistency training to enforce agreement between forward and reverse complement DNA strand predictions during training.

**Why Consistency Training?**
DNA is double-stranded, so a sequence and its reverse complement encode the same biological information. By enforcing that model predictions agree for both orientations during training, we improve model robustness and strand-handling capabilities.

**Formula:**
```
total_loss = classification_loss + alpha * consistency_loss

Where:
- classification_loss = (BCE(fwd_logits, labels) + BCE(rc_logits, labels)) / 2
- consistency_loss = MSE(sigmoid(fwd_logits), sigmoid(rc_logits))
```

**New Components:**

1. **`BidirectionalConsistencyLoss` class** (`train_stage1.py`)
   - Combined classification + consistency loss
   - Returns total loss and detailed metrics (cls_loss, consistency_loss)
   - Supports label smoothing

2. **`get_reverse_complement_tokens_for_training()`** (`train_stage1.py`)
   - Computes reverse complement tokens during training
   - Decodes tokens → computes RC → re-tokenizes

3. **`run_consistency_epoch()`** (`train_stage1.py`)
   - Training loop variant for consistency training
   - Computes both forward and RC predictions per batch
   - Optionally recomputes engineered features for RC sequences

**Config Options:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `stage1_use_consistency_loss` | true | Enable consistency training |
| `stage1_consistency_alpha` | 0.3 | Weight for consistency term (0.0-1.0) |

**Usage:**
```yaml
loss:
  stage1_use_consistency_loss: true
  stage1_consistency_alpha: 0.3
```

**Log Output:**
```
Train - Loss: 0.4521 (CLS: 0.4123, Cons: 0.0398) | Acc: 0.8542 | F1: 0.8521 | MCC: 0.7084
```

**Expected Improvement:** +0.3-0.5% accuracy from better strand handling

**Interaction with Other Features:**
- Compatible with SWA (consistency training used in regular phase, not SWA phase)
- Mutually exclusive with distillation training (distillation takes priority)
- Works with all model variants (original, transformer_only, combined)

### 2026-01-20: Positional Motif Attention Bias
Implemented learnable position-aware attention bias for known promoter regulatory regions.

**Why Positional Motif Bias?**
Promoters have known regulatory motifs at specific positions relative to the TSS:
- TATA box: -30 to -25 bp (positions ~45-55 in 81bp window)
- -10 box: -12 to -7 bp (positions ~64-74)
- TSS: position 0 (positions ~75-81)

By initializing attention bias toward these biologically important positions, the model can learn to focus on regulatory regions more effectively.

**New Files:**
- `Gene Whisperer/training/positional_motif_bias.py` - Core implementation
  - `PositionalMotifBias` class with learnable (H, L, L) bias matrix
  - `DEFAULT_PROMOTER_PRIORS` - Default bias values for promoter regions
  - `create_promoter_motif_bias()` - Factory function

- `Gene Whisperer/training/tests/test_positional_motif_bias.py` - Test suite (21 tests)

**Modified Files:**
- `Gene Whisperer/training/model.py`
  - Added `motif_bias` parameter to `PreNormTransformerLayer.forward()`
  - Added `use_positional_motif_bias` and `motif_regions` to `DNAEncoder`
  - Updated all model classes (`GeneWhispererStage1`, `GeneWhispererTransformerOnly`, `GeneWhispererCombined`, `GeneWhispererStage2`)
  - Updated `create_model_variant()` factory function

- `Gene Whisperer/training/config.yaml`
  - Added `use_positional_motif_bias: true`
  - Added `motif_regions` with TATA box, -10 box, and TSS positions

**Config Options:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_positional_motif_bias` | true | Enable positional motif attention bias |
| `motif_regions.tata_box` | [45, 55] | TATA box region |
| `motif_regions.minus_10_box` | [64, 74] | -10 box region |
| `motif_regions.tss` | [75, 81] | TSS region (highest bias) |

**How it Works:**
1. A learnable bias matrix (H, L, L) is created for each attention head
2. Bias is initialized with small values toward known motif positions:
   - TO columns: All positions attend more to motif regions (bias=0.1-0.15)
   - FROM rows: Motif positions have broader attention context (bias=0.05-0.075)
3. The bias is added to attention scores before softmax in every transformer layer
4. Compatible with RoPE (motif bias handles absolute positions, RoPE handles relative)

**Usage:**
```yaml
model:
  use_positional_motif_bias: true
  motif_regions:
    tata_box: [45, 55]
    minus_10_box: [64, 74]
    tss: [75, 81]
```

**Expected Improvement:** +0.5-1% accuracy from position-aware attention

### 2026-01-20: GeneWhispererV6 Modern Architecture
Implemented GeneWhispererV6, a modern LLM-style architecture that integrates all improvements into a clean, unified model.

**Why V6?**
V6 combines all the modern transformer improvements (RoPE, SwiGLU, RMSNorm, positional motif bias) with optional Grouped Query Attention (GQA) in a simplified, efficient architecture following the patterns used in LLaMA, Mistral, and PaLM.

**New Files:**
- `Gene Whisperer/training/model_v6.py` - Core V6 implementation
  - `ModernTransformerBlock` - Unified block with all modern features
  - `GeneWhispererV6` - Main model class
  - `create_v6_model()` - Factory function

- `Gene Whisperer/training/tests/test_model_v6.py` - Comprehensive test suite (30 tests)

**Modified Files:**
- `Gene Whisperer/training/model.py`
  - Added import for V6 model components
  - Updated `create_model_variant()` to support "v6" variant

- `Gene Whisperer/training/config.yaml`
  - Added `model_version: "v5"` option (v5 or v6)
  - Added `use_gqa: false` for Grouped Query Attention
  - Added `num_kv_heads: 4` for GQA configuration
  - Added convenience flags: `use_swiglu`, `use_rmsnorm`, `use_positional_bias`

**Architecture:**
```
Input -> Embedding -> [ModernTransformerBlock x N] -> Pool -> [+Features] -> Classifier
```

Each ModernTransformerBlock contains:
- Pre-norm with RMSNorm (or LayerNorm)
- RoPE-enhanced multi-head attention (or GQA)
- Optional positional motif bias (first layer only)
- SwiGLU FFN (or standard GELU)

**Config Options:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `model_version` | "v5" | Model version: "v5" (simplified) or "v6" (modern) |
| `model_variant` | "stage1" | Model variant: transformer_only, features_only, combined, v6 |
| `use_gqa` | false | Enable Grouped Query Attention |
| `num_kv_heads` | 4 | Number of KV heads for GQA (must divide transformer_heads) |
| `use_swiglu` | true | Use SwiGLU FFN in V6 |
| `use_rmsnorm` | true | Use RMSNorm in V6 |
| `use_positional_bias` | true | Use positional motif bias in V6 |

**Usage:**
```yaml
model:
  model_variant: "v6"  # Use V6 modern architecture
  use_rope: true
  use_swiglu: true
  use_rmsnorm: true
  use_positional_bias: true
  use_gqa: false  # Optional: enable for parameter reduction
  num_kv_heads: 4  # Only used if use_gqa: true
```

**GQA Parameter Reduction:**
With GQA (num_kv_heads=4, transformer_heads=12), KV projection parameters are reduced by ~67%, resulting in ~10-15% total parameter reduction with minimal accuracy impact.

**Expected Total Improvement from V6:** +2-4% accuracy over baseline transformer

### 2026-01-21: V6 load_pretrained_weights Bug Fix
Fixed TypeError when running Stage 1 training with V6 architecture.

**Issue:**
```
TypeError: GeneWhispererV6.load_pretrained_weights() got an unexpected keyword argument 'transfer_mode'
```

**Root Cause:**
The `GeneWhispererV6.load_pretrained_weights()` method in `model_v6.py` was missing the `transfer_mode` parameter that other model classes support. The training code in `train_stage1.py:3176` passes this parameter when loading MLM encoder weights.

**Fix:**
Updated `GeneWhispererV6.load_pretrained_weights()` to accept `transfer_mode` parameter with the same signature as other models:
- Added `transfer_mode: str = "embed_only"` parameter
- Added validation for supported modes: "embed_only", "embed_plus_adapter", "none"
- Added early return when `transfer_mode="none"`

**Files Modified:**
- `Gene Whisperer/training/model_v6.py` - Added transfer_mode parameter to load_pretrained_weights()

### 2026-01-22: Ensemble Inference Architecture Alignment
Updated `ensemble_infer.build_model()` to mirror the simplified Stage 1 architecture used in training.

**Changes:**
- Added pooling/classifier defaults from `simplified_model` config.
- Pulled `stage1_drop_path_rate` with fallback to global `drop_path_rate`.
- Wired RoPE + SwiGLU/RMSNorm settings (`use_rope`, `rope_base`, `ffn_type`, `norm_type`, `ffn_mult`).
- Added focused tests for build_model parameter wiring.

**Files Modified:**
- `Gene Whisperer/training/ensemble_infer.py`
- `Gene Whisperer/training/tests/test_ensemble_infer_build_model.py`

### 2026-01-23: Transformer Contribution Diagnostic Tool
Added a diagnostic script to measure whether the transformer is learning useful patterns or if predictive power comes entirely from engineered features.

**New Files:**
- `Gene Whisperer/training/diagnose_transformer_contribution.py` - CLI entry point
- `Gene Whisperer/training/experiment_runner.py` - Core experiment logic (model building, checkpoint loading, evaluation)
- `Gene Whisperer/training/hooks.py` - Forward hooks for ablation (randomizes transformer output)
- `Gene Whisperer/training/tests/test_diagnose_transformer.py` - 27 tests

**Three Experiments:**
| Experiment | What it tests | Method |
|-----------|---------------|--------|
| A (Full model) | Baseline accuracy | Normal forward pass |
| B (Transformer only) | Transformer contribution | Zero out engineered_features tensor |
| C (Features only) | Engineered features contribution | Replace pooled transformer output with random vectors via forward hook |

**Interpretation:**
- Exp B ~50-55%: Transformer NOT learning → improve pre-training
- Exp B >> 50%: Transformer IS learning useful patterns
- Exp C ≈ Exp A: Engineered features do all the work

**Usage:**
```bash
cd "Gene Whisperer/training"
python diagnose_transformer_contribution.py --config config.yaml
python diagnose_transformer_contribution.py --checkpoint ../artifacts/checkpoints/stage1_k6.pt
python diagnose_transformer_contribution.py --max-samples 500  # Quick test
```

**Architecture Support:**
- V5 (GeneWhispererStage1): Hooks on attention pooling module
- V6 (GeneWhispererV6): Hooks on final_norm, produces consistent per-sample random vectors so mean pooling yields random output

### 2026-01-23: Checkpoint Compatibility Verifier
Added a CLI verification helper to compare checkpoint keys/shapes against the
architecture produced by `ensemble_infer.build_model`.

**Changes:**
- `verify_checkpoint_compatibility()` prints missing/unexpected keys and shape mismatches.
- Added `--verify` CLI entrypoint for quick checks.
- Added tests covering the verification workflow.

**Files Modified:**
- `Gene Whisperer/training/ensemble_infer.py`
- `Gene Whisperer/training/tests/test_verify_checkpoint_compatibility.py`

### 2026-01-23: MLM Pre-training Quality Diagnostic Tool
Added a diagnostic script to analyze whether MLM pre-training is producing useful representations across different k-mer values.

**New Files:**
- `Gene Whisperer/training/diagnose_mlm_quality.py` - CLI entry point
- `Gene Whisperer/training/mlm_quality_metrics.py` - Core analysis logic (model building, predictions, embeddings, silhouette)
- `Gene Whisperer/training/tests/test_mlm_quality.py` - 30 tests

**What it measures:**
| Metric | Purpose |
|--------|---------|
| Top-1 Accuracy | How often the correct k-mer is the top prediction |
| Top-5 Accuracy | How often the correct k-mer is in top 5 predictions |
| Avg Entropy | Whether predictions are confident (low) or diffuse (high) |
| Silhouette Score | Whether embeddings cluster by promoter label |

**Architecture Details:**
- Loads encoder from `../artifacts/mlm_encoder_k{3,4,5,6}.pt` checkpoints
- Reads architecture params from `mlm_encoder_k{k}.metadata.json` (dim=256, layers=8, heads=8)
- Rebuilds DNAMLM-style LM head with weight tying and CLIP-style normalization
- Uses deterministic [MASK]-only masking (no 80/10/10) for consistent evaluation
- Mean-pools encoder output for embedding extraction

**Usage:**
```bash
cd "Gene Whisperer/training"
python diagnose_mlm_quality.py
python diagnose_mlm_quality.py --kmers 3 4
python diagnose_mlm_quality.py --max-samples 200 --mask-prob 0.20
python diagnose_mlm_quality.py --save-embeddings  # Saves t-SNE plots
```

**Note:** Requires MLM encoder checkpoints to exist (generated by `pretrain_mlm.py`). Currently only metadata files exist; the actual `.pt` checkpoints need to be generated by running pre-training.

### 2026-01-24: BPE Tokenization (DNABERT-2 Style)
Implemented Byte Pair Encoding tokenizer as an alternative to k-mer tokenization, following the DNABERT-2 paper (Zhou et al., ICLR 2024).

**Why BPE?**
DNABERT-2 demonstrated that BPE tokenization is superior to k-mer for genome foundation models:
- **No information leakage** (unlike overlapping k-mer where adjacent tokens reveal masked content)
- **No sample inefficiency** (unlike non-overlapping k-mer where a 1bp shift completely changes tokenization)
- **~5x sequence length reduction** (reduces quadratic attention cost)
- **Variable-length tokens** capture biologically meaningful motifs (TATA, CAAT, etc.)
- **Naturally challenging MLM** (model must predict both token length and content)

**New Files:**
- `Gene Whisperer/training/bpe_tokenizer.py` - Core BPE implementation
  - `DNABPETokenizer` class with train, tokenize, decode, save, load
  - Default vocab_size=4096 (optimal per DNABERT-2 empirical analysis)
  - Special tokens: [PAD], [UNK], [CLS], [SEP], [MASK]
  - `encode_with_special()` for [CLS]+tokens+[SEP] wrapping
  - `compression_ratio()` and `get_token_lengths()` for analysis

- `Gene Whisperer/training/tests/test_bpe_tokenizer.py` - 35 tests covering:
  - Initialization and special token IDs
  - BPE training algorithm (merge ordering, frequency thresholds, dirty input)
  - Tokenization (merge application, case normalization, empty input)
  - Decode roundtrip correctness
  - Save/load persistence
  - Compression statistics
  - Biological motif capture (TATA, repeat regions)

**Algorithm:**
1. Initialize vocab with special tokens + {A, T, C, G}
2. Split all sequences into character-level tokens
3. Count adjacent pair frequencies across corpus
4. Merge most frequent pair into new token
5. Repeat until vocab_size reached or min_frequency not met
6. Tokenization applies learned merges in priority order (greedy)

**Usage:**
```python
from bpe_tokenizer import DNABPETokenizer

# Train on genomic sequences
tok = DNABPETokenizer(vocab_size=4096)
tok.train(sequences, min_frequency=2)

# Tokenize
ids = tok.tokenize("ATCGATCGATCG")
ids_special = tok.encode_with_special("ATCGATCG")

# Check compression
ratio = tok.compression_ratio("ATCGATCGATCG")  # e.g., 3.5x

# Save/load
tok.save("bpe_vocab.json")
loaded = DNABPETokenizer.load("bpe_vocab.json")
```

**Key Design Decisions:**
- Vocab size 4096 chosen per DNABERT-2 Figure 3c (best performance/efficiency tradeoff)
- Greedy merge-order tokenization (applies merges in training priority order)
- Non-ATCG characters stripped during sanitization (handles N, lowercase, etc.)
- JSON serialization for vocabulary portability

### 2026-01-24: BPE Pipeline Integration (Complete)
Fully integrated BPE tokenization across the entire training pipeline, replacing k-mer tokenization.

**Files Modified:**
| File | Changes |
|------|---------|
| `bpe_tokenizer.py` | Added `tokenize_and_pad()`, `random_token_id()`, `itos` property |
| `train_bpe_vocab.py` | **NEW** - Standalone script to train BPE vocabulary from FASTA |
| `dataset.py` | Replaced `KmerVocabulary` with `DNABPETokenizer`, updated datasets |
| `model.py` | Removed `kmer` parameter from all model classes |
| `model_v6.py` | Updated `create_v6_model()` to use `max_token_len` instead of kmer |
| `pretrain_mlm.py` | Replaced k-mer MLM with BPE-based MLM |
| `train_stage1.py` | Updated vocab handling, model creation, checkpoint paths |
| `ensemble_infer.py` | Simplified to single BPE model (removed multi-kmer ensemble) |
| `tta.py` | Updated reverse complement to use BPE decode/encode |
| `config.yaml` | Added BPE settings, removed kmer settings |

**Key Configuration Changes:**
| Parameter | Old (k-mer) | New (BPE) |
|-----------|-------------|-----------|
| `vocab_size` | 4099 | 4096 |
| `pad_token_id` | 4098 | 0 |
| `kmer` | 6 | *removed* |
| `max_bp_len` | 81 | *replaced by max_token_len* |
| `max_token_len` | N/A | 24 (for 81bp sequences) |
| `mlm_max_token_len` | N/A | 64 (for 234bp MLM sequences) |
| `mlm_max_span_len` | 3 | 1 (independent token masking) |
| `use_positional_bias` | true | false (positions don't map to bp) |

**New API Methods:**
```python
# Tokenize with padding to fixed length
tokens = vocab.tokenize_and_pad(sequence, max_token_len=24)  # -> torch.LongTensor

# Get random token ID for MLM replacement
rand_id = vocab.random_token_id()  # excludes special tokens

# Index-to-string list for compatibility
tokens_list = vocab.itos  # List[str]
```

**Training Workflow:**
```bash
# 1. Train BPE vocabulary (run once)
python train_bpe_vocab.py --config config.yaml

# 2. Run MLM pre-training with BPE
python pretrain_mlm.py --config config.yaml

# 3. Fine-tune Stage 1
python train_stage1.py --config config.yaml

# 4. Inference
python ensemble_infer.py --sequence "ATGC..." --config config.yaml
```

**Checkpoint Naming:**
- MLM encoder: `mlm_encoder_bpe.pt` (was `mlm_encoder_k6.pt`)
- Stage 1: `stage1_bpe.pt` (was `stage1_k6.pt`)
- BPE vocab: `bpe_vocab.json` (new)

**Test Status:** All 401 tests passing.

### 2026-01-25: BPE Vocabulary Training Completed
Successfully trained BPE vocabulary with 4096 tokens on DNA sequences.

**Training Results:**
- Vocabulary size: 4096 tokens
- Compression ratio: 4.97x average (1000 sample sequences)
- Average tokens per 234bp MLM window: 47.2
- Max tokens per MLM window: 53
- Token length distribution: 1-mer to 25-mer

**Recommended vs Configured max_token_len:**
| Sequence Type | Recommended | Configured | Status |
|---------------|-------------|------------|--------|
| 81bp fine-tuning | 22 | 24 | ✓ OK |
| 234bp MLM pre-training | 58 | 64 | ✓ OK |

**File Location:**
- BPE vocab saved to: `Gene Whisperer/artifacts/vocabs/bpe_vocab.json`
- Config path: `bpe_vocab_path: ../artifacts/vocabs/bpe_vocab.json`

**Next Steps:**
1. Run MLM pre-training: `python pretrain_mlm.py --config config.yaml`
2. Fine-tune Stage 1: `python train_stage1.py --config config.yaml`

### 2026-01-25: Data Loading Performance Optimization
Fixed critical data loading bottleneck where GPU was idle 63% of the time waiting for data.

**Problem Identified:**
Telemetry showed `data=6215ms [63%], compute=3651ms [37%]` per step, meaning the GPU was only doing useful work 1/3 of the time. Root causes:
1. BPE tokenization applied 4000 merge rules sequentially per sequence
2. Engineered features (TNC, PseEIIP, CKSNAP, PSTNP) computed on every `__getitem__` call
3. Reverse complement augmentation (50%) invalidated caching
4. Low prefetch_factor (2) in DataLoader

**Optimizations Implemented:**

| Component | Before | After | Improvement |
|-----------|--------|-------|-------------|
| BPE tokenization | O(merges × seq_len) | O(seq_len × max_token_len) | ~100x faster |
| Token caching | None (on-the-fly) | Pre-cached at init | No computation at runtime |
| Feature caching | None (on-the-fly) | Pre-cached at init | No computation at runtime |
| RC augmentation | Tokenize on-demand | Pre-cache both fwd+rev | Just tensor indexing |
| prefetch_factor | 2 | 4 | 2x more batches prefetched |

**Files Modified:**
| File | Changes |
|------|---------|
| `bpe_tokenizer.py` | Replaced O(merges) tokenization with greedy longest-match |
| `dataset.py` | Pre-cache tokens+features for fwd and RC at init |
| `config.yaml` | Increased `prefetch_factor` from 2 to 4 |

**New BPE Tokenization Algorithm:**
```python
# Old: Apply 4000 merges sequentially
for left, right in self.merges:  # 4000 iterations
    tokens = self._apply_merge_to_tokens(tokens, left, right)

# New: Greedy longest-match (single pass)
while i < seq_len:
    for length in range(max_len, 0, -1):  # Max ~25 iterations
        if sequence[i:i+length] in self.vocab:
            tokens.append(vocab[substring])
            i += length
            break
```

**New Dataset Caching Strategy:**
```python
# Pre-compute at __init__:
# - tokens_fwd[N, max_token_len] for all sequences
# - tokens_rev[N, max_token_len] for reverse complements
# - engineered_fwd[N, 288] for all sequences
# - engineered_rev[N, 288] for reverse complements

# At __getitem__ (fast path):
use_reverse = random.random() < 0.5
tokens = tokens_rev[idx] if use_reverse else tokens_fwd[idx]
engineered = engineered_rev[idx] if use_reverse else engineered_fwd[idx]
```

**Expected Results:**
- GPU utilization should increase from 37% to 80-95%
- Training speed should be 2-3x faster
- Memory usage increases slightly due to caching (~2x dataset size)

**Trade-offs:**
- Increased initialization time (pre-computation happens once at startup)
- Higher memory usage (caching both fwd and rev versions)
- Base substitution augmentation (rare, usually disabled) still requires on-the-fly computation

### 2026-01-25: MLM Pretraining Fixes (BPE Naming & CLIP Projection)
Fixed multiple issues discovered after first BPE pretraining run that achieved only 11% accuracy.

**Issues Fixed:**

1. **Checkpoint Naming Bug (k=6 / k=1 confusion)**
   - `PretrainingEarlyStopping` was logging "Saved k=6 MLM checkpoint" even when using BPE
   - `get_kmer_pretrain_config()` was generating `mlm_encoder_k1.pt` paths for BPE
   - **Fix**: Added `tokenizer_type` parameter to `PretrainingEarlyStopping`, updated checkpoint naming to use "BPE" label
   - **Fix**: Modified `get_kmer_pretrain_config()` to preserve BPE paths from config when `kmer=1`

2. **Low Accuracy Root Cause: CLIP-Style Projection**
   - DNAMLM was using CLIP-style normalized cosine similarity projection:
     ```python
     h = F.normalize(x, dim=-1)
     w = F.normalize(self.lm_head.weight, dim=-1)
     logits = logit_scale * (h @ w.T)  # Bounded to [-scale, +scale]
     ```
   - This is designed for contrastive learning (image-text), NOT masked language modeling
   - The normalization constrains logits to a bounded range, limiting learning capacity
   - **Fix**: Added `use_clip_projection` parameter (default: `False`)
   - Standard MLM projection now used by default: `logits = lm_head(x)`

3. **Length Curriculum Warning**
   - Warning "Length curriculum disabled for streaming datasets" always appeared
   - This is expected behavior (length curriculum incompatible with streaming)
   - **Fix**: Set `mlm_length_curriculum_enabled: false` in config since streaming is used

**Files Modified:**
| File | Changes |
|------|---------|
| `pretrain_mlm.py` | Added `use_clip_projection` param to DNAMLM, updated `get_kmer_pretrain_config()` for BPE paths |
| `early_stopping.py` | Added `tokenizer_type` param to `PretrainingEarlyStopping`, updated checkpoint naming/logging |
| `config.yaml` | Added `mlm_use_clip_projection: false`, disabled length curriculum |

**Config Changes:**
| Parameter | Old Value | New Value | Purpose |
|-----------|-----------|-----------|---------|
| `mlm_use_clip_projection` | N/A (always on) | `false` | Use standard MLM projection |
| `mlm_length_curriculum_enabled` | `true` | `false` | Avoid warning with streaming |

**Expected Improvement:**
- Accuracy should improve significantly (target: 70-85%) with standard projection
- No more confusing k=6/k=1 naming in logs
- No more length curriculum warning

**Next Steps:**
1. Re-run MLM pretraining with fixed configuration
2. Monitor accuracy - should see steady improvement above 50% within first few epochs
3. If accuracy still low, investigate masking strategy or learning rate

### 2026-01-26: CRITICAL FIX - BPE Token Padding Bug (Root Cause of 11% Accuracy)

**The Real Problem:**
After first fix attempt, model still stuck at 11.4% accuracy. Investigation revealed:
- `masked_fraction` was only **2.67%** instead of expected **15%**
- This indicated massive over-padding

**Root Cause Found:**
In `GenomeMLMIterableDataset._sample_tokens()` (line 632):
```python
# BUG: Using window_bp (234) instead of max_token_len (64)
tokens = self.vocab.tokenize_and_pad(window, self.window_bp)  # WRONG!
```

This caused:
- 234bp DNA sequence → ~60 BPE tokens
- But padded to **234 tokens** (using `window_bp` instead of `max_token_len`)
- Result: **174 PAD tokens** out of 234 = **74% padding!**
- Only 26% non-PAD tokens × 15% mask rate = **3.9% masked** (matches observed 2.67%)

**The Fix:**
1. Added `max_token_len` parameter to `GenomeMLMIterableDataset.__init__`
2. Changed tokenization call to use `self.max_token_len` instead of `self.window_bp`
3. Updated all call sites to pass `mlm_max_token_len` from config (default: 64)
4. Added logging to show both `window_bp` and `max_token_len` during training

**Files Modified:**
| File | Change |
|------|--------|
| `pretrain_mlm.py:499` | Added `max_token_len` parameter to `GenomeMLMIterableDataset.__init__` |
| `pretrain_mlm.py:527-530` | Store `max_token_len` with fallback to `window_bp` for legacy compatibility |
| `pretrain_mlm.py:632` | Use `self.max_token_len` for tokenization padding |
| `pretrain_mlm.py:3712` | Read `mlm_max_token_len` from config for streaming training |
| `pretrain_mlm.py:3793,3807` | Pass `max_token_len` to train/val dataset constructors |
| `pretrain_mlm.py:2802,2861,2875` | Same for dry run function |

**Expected Results After Fix:**
- `masked_fraction` should be ~15% (was 2.67%)
- With proper masking, accuracy should improve significantly (target: 70-85%)
- Perplexity should drop from ~770 to under 50

**Key Insight:**
The bug was subtle: `window_bp` (DNA base pairs) vs `max_token_len` (BPE tokens) have different meanings:
- `window_bp=234`: DNA sequence length in base pairs
- `max_token_len=64`: Target token sequence length after BPE compression (~4x)

Using DNA length for token padding created massive over-padding, severely limiting what could be masked during MLM training.

### 2026-01-26: Mode Collapse Fix (10-11% Accuracy Ceiling)

**Symptom:**
After padding fix, model still stuck at 10-11% accuracy despite improved masked_fraction (~10%).

**Diagnosis:**
- Loss = 7.77, but Accuracy = 10.8%
- This indicates **mode collapse**: model predicting same dominant token for all positions
- That token appears ~10% in labels, giving 10% accuracy by chance
- The model wasn't learning meaningful representations

**Root Cause Analysis:**
The combination of:
1. Low learning rate (0.0002) - model settling into local minimum
2. Weight tying (mlm_tie_weights: true) - constraining input/output representations
3. No label smoothing - allowing overconfident predictions

**Fixes Applied:**
| Parameter | Old Value | New Value | Rationale |
|-----------|-----------|-----------|-----------|
| `mlm_lr` | 0.0002 | 0.0005 | Higher LR to escape local minima |
| `mlm_tie_weights` | true | false | Allow different input/output representations |
| Label smoothing | 0.0 | 0.1 | Prevent overconfident predictions |

**Files Modified:**
| File | Change |
|------|--------|
| `config.yaml:211` | `mlm_lr: 0.0005` |
| `config.yaml:220` | `mlm_tie_weights: false` |
| `pretrain_mlm.py:2180` | Added `label_smoothing=0.1` to cross_entropy |

**Expected Results:**
- Model should start learning and accuracy should steadily increase
- Watch for accuracy >30% within first 5-10 epochs to confirm fix is working
- Target accuracy: 70-85% at convergence
