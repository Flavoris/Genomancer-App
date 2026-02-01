# Genomancer Long-Term Memory

## Project Overview
Genomancer (Gene Whisperer) is a DNA sequence analysis tool that uses transformer-based models for promoter prediction. The project uses:
- BPE tokenization (DNABERT-2 style) with 4096 vocabulary
- 12-layer transformer encoder with 384 embedding dim
- Two-stage training: MLM pretraining → supervised fine-tuning

## Key Architecture Decisions

### MLM Pretraining Settings (Updated 2026-02-01)
- **LR**: 0.0005 (increased from 0.0001 - too conservative for DNA)
- **Weight Decay**: 0.01 (reduced from 0.05 - was slowing learning)
- **Effective Batch**: 512 (reduced from 2048 for faster learning)
- **Warmup**: 2% of training (reduced from 6%)
- **Architecture**: BERT-style (LayerNorm, GELU, no RoPE) for bidirectional MLM

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

## Current Issues

### MLM Accuracy Plateau (18%)
- **Status**: Being addressed
- **Root Cause**: Combination of low LR, large effective batch, and inherent difficulty of DNA MLM
- **Fix Applied**: Increased LR 5x, reduced batch 4x, faster warmup

## Important Files
- config.yaml - Main training configuration
- pretrain_mlm.py - MLM pretraining script
- pretrain_components/mlm_model.py - DNAMLM model
- bpe_tokenizer.py - BPE tokenizer
- diagnose_mlm_issue.py - Diagnostic tool

## Test Commands
\`\`\`bash
# Run all tests
python -m pytest tests/ -x -q

# Run MLM-specific tests
python -m pytest training/tests/test_mlm_quality.py training/tests/test_colab_run_mlm.py -x -q

# Run diagnostic
python Gene\ Whisperer/training/diagnose_mlm_issue.py
\`\`\`
