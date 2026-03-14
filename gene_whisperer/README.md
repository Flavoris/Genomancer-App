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
- `mlm.mask_ambiguous_tokens` controls whether tokens containing `N` are mask targets (default `false` to reduce ambiguous-noise loss).
- `mlm.min_masked_tokens` enforces a minimum number of supervised mask positions per sample.
- `mlm.min_maskable_tokens` and `mlm.resample_attempts` resample windows to avoid low-signal batches with too few maskable targets.
- `mlm.min_tokenized_tokens` keeps BPE-pretrained windows dense enough to use the MLM token budget instead of training mostly on short tokenized spans.
- Pretraining `mlm.window_size` should be much larger than promoter sequence length; the MLM stage benefits from longer whole-genome context.
- Tokenizer training now samples multiple random windows across long genomes instead of taking only one window per FASTA record, which produces a more representative BPE vocabulary.
- Keep tokenizer corpus budgets modest (`mlm.tokenizer_max_bases`, `mlm.tokenizer_max_sequences`, `mlm.tokenizer_window_size`) because the current pure-Python BPE trainer is intentionally simple and scales with sampled corpus size.
- Tokenizer quality is controlled with `mlm.tokenizer_min_freq` and `mlm.tokenizer_max_token_length` to avoid a long tail of ultra-rare BPE tokens.
- Existing tokenizer files are automatically retrained when tokenizer metadata does not match the requested config (`mlm.tokenizer_retrain_if_mismatch`).
- MLM pretraining uses early stopping (`training.min_epochs`, `training.early_stopping_patience`, `training.early_stopping_min_delta`) and saves `mlm_best.pt` as the best checkpoint.
- MLM optimization uses warmup+cosine LR (`training.warmup_ratio`, `training.min_lr_ratio`), AMP (`training.use_amp`), and gradient clipping (`training.max_grad_norm`).
- The tokenizer is saved to `gene_whisperer/artifacts/bpe_tokenizer.json`.
- PSTNP matrices are saved to `gene_whisperer/artifacts/finetune/pstnp_stage*.json`.
- Keep configs small for iPhone deployment; the defaults target a lightweight encoder.
