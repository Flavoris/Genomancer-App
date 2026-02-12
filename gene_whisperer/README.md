# Gene Whisperer (Genomancer)

This folder houses the new Gene Whisperer model pipeline:

- **BPE tokenization** (DNABERT2-inspired) via `gene_whisperer/tokenization/bpe.py`
- **Engineered features** (PROCABLES-inspired) via `gene_whisperer/features/engineered.py`
- **Transformer encoder** + feature fusion via `gene_whisperer/models/`
- **Training scripts** for MLM pretraining and promoter fine-tuning via `gene_whisperer/training/`

## Quick start (local)

```bash
# Pretrain MLM
python gene_whisperer/training/pretrain_mlm.py --config gene_whisperer/configs/pretrain.yaml

# Fine-tune stage1 (promoter vs non)
python gene_whisperer/training/finetune_promoter.py --config gene_whisperer/configs/finetune.yaml --task stage1

# Fine-tune stage2 (strong vs weak)
python gene_whisperer/training/finetune_promoter.py --config gene_whisperer/configs/finetune.yaml --task stage2
```

## Notes

- Genome FASTA paths are configured in `gene_whisperer/configs/pretrain.yaml`.
- `training.samples_per_epoch` controls how many random MLM windows are trained per epoch.
- `mlm.max_bases_per_file` can cap FASTA read size for faster Colab bring-up runs.
- `mlm.sample_by_length` keeps genome sampling proportional to sequence length instead of equally weighting contigs.
- Tokenizer training uses sampled windows (`mlm.tokenizer_max_bases`, `mlm.tokenizer_max_sequences`, `mlm.tokenizer_window_size`) to avoid multi-hour startup on full genomes.
- MLM pretraining uses early stopping (`training.min_epochs`, `training.early_stopping_patience`, `training.early_stopping_min_delta`) and saves `mlm_best.pt` as the best checkpoint.
- MLM optimization uses warmup+cosine LR (`training.warmup_ratio`, `training.min_lr_ratio`), AMP (`training.use_amp`), and gradient clipping (`training.max_grad_norm`).
- The tokenizer is saved to `gene_whisperer/artifacts/bpe_tokenizer.json`.
- PSTNP matrices are saved to `gene_whisperer/artifacts/finetune/pstnp_stage*.json`.
- Keep configs small for iPhone deployment; the defaults target a lightweight encoder.
