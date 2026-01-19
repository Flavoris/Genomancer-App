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
