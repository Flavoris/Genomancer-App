# Genomancer Long-Term Memory

## Project Overview
Genomancer (Gene Whisperer) is a DNA sequence analysis tool that uses transformer-based models for promoter prediction. The current (2026-02-04) rebuild uses:
- Custom BPE tokenizer (DNABERT2-inspired) with 4096 vocab
- Lightweight transformer encoder (default: 6 layers, 256 embedding dim) for mobile-friendly deployment
- Two-stage training: MLM pretraining → supervised fine-tuning (stage1 promoter vs non, stage2 strong vs weak)
- PROCABLES-style engineered features (CKSNAP, TNC, PseEIIP, PSTNP, GC/AT) fused with transformer embeddings

## Key Architecture Decisions

### Current Config Locations (2026-02-04)
- `gene_whisperer/configs/pretrain.yaml` - MLM pretraining config
- `gene_whisperer/configs/finetune.yaml` - stage1/stage2 fine-tuning config
- `gene_whisperer/artifacts/bpe_tokenizer.json` - saved tokenizer
- `gene_whisperer/artifacts/finetune/pstnp_stage*.json` - saved PSTNP matrices

### Colab Pretraining Notes (2026-02-08)
- MLM pretraining now prints startup/device and per-batch progress logs.
- Pretraining config now includes `training.samples_per_epoch`, `training.log_interval`, and `training.num_workers`.
- FASTA loading supports `mlm.max_bases_per_file` for faster smoke runs on very large genomes.
- Config paths are resolved relative to the YAML file location to avoid Colab cwd issues.

### Tokenizer Reliability Notes (2026-02-10)
- BPE tokenizer training now logs merge progress to avoid "silent" long-running Colab cells.
- Tokenizer fitting now supports bounded sampled corpora via:
  - `mlm.tokenizer_max_bases`
  - `mlm.tokenizer_max_sequences`
  - `mlm.tokenizer_window_size`
- These tokenizer limits improve runtime stability and do not cap MLM pretraining data unless `mlm.max_bases_per_file` is set.

### MLM Stop/Save Policy (2026-02-10)
- MLM pretraining no longer stops at a fixed small epoch count by default.
- Training now uses early stopping on loss with configurable `min_epochs`, `early_stopping_patience`, and `early_stopping_min_delta`.
- Checkpoint policy defaults to best-only (`save_best_only: true`), writing `mlm_best.pt` instead of saving every epoch checkpoint.

### MLM Convergence Upgrade (2026-02-12)
- Pretraining now uses warmup + cosine LR decay, gradient clipping, and AMP controls from config.
- Optimizer now excludes bias/normalization vectors from weight decay (AdamW best practice) to improve convergence.
- MLM dataset now:
  - samples genomes proportional to sequence length (`mlm.sample_by_length`)
  - seeds RNG per worker to avoid duplicated random streams with `num_workers > 0`
  - guarantees at least one masked token per sample
  - avoids reserved tokens for 10% random replacement in MLM masking.
- MLM head upgraded to BERT-style transform (`Linear -> GELU -> LayerNorm -> Dropout -> vocab projection + bias`) with tied token embeddings.
- Default pretrain config moved to a stronger recipe (larger encoder + longer schedule + lower LR + larger effective batch) to push loss lower before early stopping.

### MLM Plateau Mitigation (2026-02-13)
- Added ambiguity-aware masking controls in pretraining config:
  - `mlm.mask_ambiguous_tokens`
  - `mlm.min_masked_tokens`
- MLM dataset now supports excluding tokens containing `N` from mask targets (default off for noisy ambiguous regions) and enforcing a minimum masked-target count per sample.
- Default pretrain schedule tuned for deeper convergence:
  - lower LR (`2e-4`)
  - longer horizon (`epochs: 320`, `min_epochs: 80`, `patience: 50`, `min_delta: 2e-4`)
  - larger sample budget (`samples_per_epoch: 90000`)
  - stronger accumulation (`grad_accum_steps: 3`)
  - reduced dropout (`0.05`) for better fit in MLM stage.

### Transformer Warning Handling (2026-02-21)
- PyTorch warning observed in Colab:
  - `enable_nested_tensor is True, but self.use_nested_tensor is False because encoder_layer.norm_first was True`
- Not a correctness bug; it indicates nested-tensor optimization is incompatible with pre-norm encoder layers.
- Kept pre-norm (`norm_first=True`) for stability and explicitly disabled nested tensor in transformer construction to remove noisy warnings.

### MLM Plateau Around Loss ~6 (2026-02-21)
- Added low-signal window handling in MLM dataset:
  - `mlm.min_maskable_tokens`
  - `mlm.resample_attempts`
- Dataset now resamples windows when too few maskable targets are present (common in ambiguous `N`-heavy regions), reducing wasted updates.
- Masking logic now prevents masking all targetable tokens in a sample, keeping at least one context token when possible.
- Pretraining loop now:
  - skips no-target batches safely
  - tracks `trained_batches`, `skipped_batches`, and `supervised_tokens` per epoch
  - raises a clear error if an epoch has zero supervised targets.

### Tokenizer Reuse + Tail-Control Fix (2026-03-13)
- Existing tokenizer files were previously reused blindly, which meant new tokenizer settings in Colab had no effect if `bpe_tokenizer.json` already existed on Drive.
- Tokenizer files now carry metadata and are retrained automatically when config-critical settings do not match:
  - `mlm.vocab_size`
  - `mlm.tokenizer_min_freq`
  - `mlm.tokenizer_max_token_length`
- Added tokenizer controls to reduce ultra-rare long BPE tokens that can keep MLM loss artificially high:
  - `mlm.tokenizer_min_freq`
  - `mlm.tokenizer_max_token_length`
  - `mlm.tokenizer_retrain_if_mismatch`
- Default tokenizer sampling budget increased (`tokenizer_max_bases`, `tokenizer_max_sequences`) to build a more stable vocabulary across organisms.

### Pretraining Window Utilization Fix (2026-03-13)
- MLM pretraining had been using a short `234 bp` window inherited from promoter-task sequence length, which underused the BPE token budget during whole-genome pretraining.
- Pretraining config now uses a larger MLM window (`mlm.window_size: 1024`) and a minimum tokenized-length guard (`mlm.min_tokenized_tokens`).
- MLM dataset now resamples windows if BPE tokenization yields too few non-padding tokens, not just too few maskable targets.
- Epoch loss reporting now uses token-weighted averaging instead of averaging per-batch means, preventing variable-target batches from skewing the reported loss and early-stopping signal.

### Why 18% MLM Accuracy is Hard to Improve
1. **DNA is fundamentally harder than text for MLM**
   - Random baseline: 0.02% (1/4091 tokens)
   - Text has strong grammar patterns; DNA has weaker biological patterns
   - BPE tokens may not capture biological motifs well

2. **Masking strategy creates a hard problem**
   - 80% [MASK]: must predict from context (hardest)
   - 10% random: must ignore misleading input
   - 10% unchanged: theoretically easy but model doesn't know which is which

3. **Diagnostic findings (2026-02-01)**
   - Model CAN overfit to 100% on fixed batch (architecture works)
   - Model struggles on random DNA (no patterns to learn)
   - Real DNA reaches 18% (learning biological patterns)

## Legacy Notes (Before Reset)
The previous Gene Whisperer folder was moved out of the repo for a clean rebuild. Older notes about MLM accuracy plateaus (18%) and 12-layer/384-dim configs are retained only as historical context.

## Important Files
- `gene_whisperer/tokenization/bpe.py` - BPE tokenizer implementation
- `gene_whisperer/features/engineered.py` - engineered features + PSTNP
- `gene_whisperer/models/` - transformer, MLM head, promoter classifier
- `gene_whisperer/training/pretrain_mlm.py` - MLM pretraining script
- `gene_whisperer/training/finetune_promoter.py` - fine-tuning script

## Test Commands
\`\`\`bash
# Run all tests
python -m pytest -q
\`\`\`
