# Cursor Prompt Script: Ship Gene Whisperer (Transformer, 2-stage)

Use these copy-paste prompts (in order) with Cursor to build and integrate the Gene Whisperer Transformer feature. Swap paths if your workspace differs.

---

## 1) Scaffold ML workspace
```
Act as a senior ML engineer. Create the following under Gene Whisperer/:
- data/ (empty placeholders: stage1_train.csv, stage1_val.csv, stage2_train.csv, stage2_val.csv, stage2_test.csv)
- training/ (Python code)
- artifacts/ (empty; will hold checkpoints/CoreML)
No code yet; just ensure dirs exist.
```

## 2) Add Python env + deps
```
Create a setup note in Gene Whisperer/training/README.md with the exact commands to:
- python -m venv ../.venv
- source ../.venv/bin/activate
- pip install torch torchvision torchaudio coremltools pandas numpy scikit-learn pyyaml tqdm
Keep it concise.
```

## 3) Write config.yaml
```
Add Gene Whisperer/training/config.yaml with defaults:
max_bp_len: 81
batch_size: 64
lr: 0.001
weight_decay: 0.0005
epochs: 30
early_stopping_patience: 5
train_val_split: 0.9
delimiter: "\t"
stage1_train: ../data/stage1_train.tsv
stage1_val: ../data/stage1_val.tsv
stage2_train: ../data/stage2_train.tsv
stage2_val: ../data/stage2_val.tsv
stage2_test: ../data/stage2_test.tsv
stage2_strength_positive: 1   # 1 = strong
stage2_strength_negative: 0   # 0 = weak
metrics_stage1: ["accuracy","precision","recall","mcc"]
metrics_stage2: ["accuracy","mcc","roc_auc"]
```

## 4) Implement dataset.py (tokenizer + features)
```
Create Gene Whisperer/training/dataset.py:
- Build 64-size 3-mer vocab (AAA..TTT) with deterministic ordering.
- tokenize_kmers(seq,max_bp_len): uppercase/strip, pad/truncate bp to max_bp_len, make overlapping 3-mers, map to ids, pad/truncate to max_tokens=max_bp_len-2, return LongTensor.
- compute_tnc(seq): 64-dim normalized tri-nucleotide frequency.
- compute_pseeiip(seq): 64-dim PseEIIP features using published EIIP constants; unknowns -> zero.
- compute_cksnap(seq, max_gap=5): 96-dim composition of k-spaced nucleotide pairs (16 pairs × 6 gaps).
- compute_pstnp(seq): 64-dim position-weighted trinucleotide frequencies (emphasize -10/-35 regions).
- TSV/CSV loading using cfg.delimiter.
- PromoterDatasetStage1 returns (tokens, engineered[288], label_float) for columns sequence,is_promoter (0/1).
- PromoterDatasetStage2 returns (tokens, engineered[288], label_float) for columns sequence,strength using cfg.stage2_strength_positive (1=strong) and cfg.stage2_strength_negative (0=weak).
- build_dataloaders(cfg): returns train/val/test loaders for both stages; random split if val missing; use WeightedRandomSampler when imbalance >20%.
Include minimal validation/logging and seeds.
```

## 5) Implement model.py (Transformer + heads)
```
Create Gene Whisperer/training/model.py:
- TransformerBackbone: embedding_dim=48, vocab=64, learned positional embeddings up to max_tokens, 4 encoder layers with LayerNorm -> MHA (heads=4, dropout=0.1) + residual -> LayerNorm -> FFN 48->64->48 + residual.
- PromoterClassifier head: Conv1d(48->128,k=5,pad=2)+ReLU -> BiLSTM(128 per dir) -> global max pool -> Dense128+ReLU -> Dropout0.5 -> Dense1 (logits).
- StrengthClassifier head: Conv1d(48->96,k=5,pad=2)+ReLU -> LSTM(128) -> global avg pool -> concat with TNC/PseEIIP (each 64-dim) -> Dense128+ReLU -> Dropout0.5 -> Dense1 (logits).
- GeneWhispererStage1(backbone+head) and GeneWhispererStage2(backbone+head) share backbone; Stage2 exposes freeze_lower_layers() to freeze bottom half (layers 0-1).
No custom ops.
```

## 6) Implement train_stage1.py
```
Create Gene Whisperer/training/train_stage1.py:
- CLI: python train_stage1.py --config config.yaml
- Loads cfg, dataloaders, backbone+Stage1 head.
- Loss: BCEWithLogitsLoss; metrics: accuracy, precision, recall, MCC.
- Early stopping on val loss/MCC; checkpoints best to ../artifacts/checkpoints/stage1_best.pt.
- Save ../artifacts/stage1_report.json (train/val metrics, class balance).
Use tqdm; set seeds for reproducibility.
```

## 7) Implement train_stage2.py
```
Create Gene Whisperer/training/train_stage2.py:
- CLI: python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt
- Load backbone weights; freeze bottom half of Transformer layers (layers 0-1); optionally freeze embeddings via a flag.
- Train on Stage2 dataset with aux features.
- Loss: BCEWithLogitsLoss; metrics: accuracy, MCC, ROC-AUC.
- Early stopping; save best to ../artifacts/checkpoints/stage2_best.pt and ../artifacts/stage2_report.json.
```

## 8) Implement inference.py
```
Create Gene Whisperer/training/inference.py:
- Functions predict_promoter(sequence) and predict_strength(sequence).
- CLI: python inference.py --sequence ACTG... --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt --stage2_ckpt ../artifacts/checkpoints/stage2_best.pt --config config.yaml
- If Stage1 predicts non-promoter, skip Stage2 and print that it was skipped.
```

## 9) Implement export_coreml.py (two models)
```
Create Gene Whisperer/training/export_coreml.py:
- Load stage1_best.pt and stage2_best.pt.
- Dummy inputs: tokens (1,max_tokens) plus engineered (1,128) for Stage1/Stage2.
- Convert with coremltools.convert to two models:
  - GeneWhispererPromoter.mlmodel/mlpackage
  - GeneWhispererStrength.mlmodel/mlpackage
- Set metadata: author, short_description, input_desc, class_labels.
- Save to ../artifacts/.
```

## 10) Add data placeholder note (you do this)
```
Add a note in Gene Whisperer/training/README.md reminding:
- You must drop CSVs at data/stage1_train.csv, stage1_val.csv (sequence,is_promoter=0/1).
- You must drop CSVs at data/stage2_train.csv, stage2_val.csv (sequence,strength in {strong,weak}).
```

## 11) Swift predictor
```
Add apps/ios/Genomancer/Genomancer/App/GeneWhispererPredictor.swift:
- Shared 3-mer tokenizer mirroring Python (pad/truncate to max_bp_len, overlapping 3-mers -> ids).
- Feature builders for TNC and PseEIIP (64-dim each) matching Python.
- Load two CoreML models (promoter + strength).
- predictPromoter(sequence) -> (isPromoter: Bool, score: Double); predictStrength(sequence) -> (isStrong: Bool, score: Double).
- Caller should skip predictStrength if predictPromoter is false.
```

## 12) Swift UI
```
Add apps/ios/Genomancer/Genomancer/App/GeneWhispererView.swift:
- @State sequence; validate A/C/G/T/N only, length <= max_bp_len.
- Button "Predict": runs Stage1; if promoter, runs Stage2; show scores.
- Display validation errors and a loading indicator.
```

## 13) Wire navigation
```
Update HomeView (apps/ios/Genomancer/Genomancer/App/HomeView.swift) to include a NavigationLink/card to GeneWhispererView labeled "Gene Whisperer (AI)" with a subtitle.
```

## 14) Add models to Xcode target
```
Create apps/ios/Genomancer/Genomancer/Resources/GeneWhisperer/ and add both CoreML models (promoter + strength) to the iOS target resources.
```

## 15) Quick runbook
```
- User step: place stage1 and stage2 TSVs into Gene Whisperer/data/ (delimiter is tab; strength: 1=strong, 0=weak).
- From Gene Whisperer/training:
  python train_stage1.py --config config.yaml
  python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt
  python export_coreml.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt --stage2_ckpt ../artifacts/checkpoints/stage2_best.pt
- Copy both exported models into Resources/GeneWhisperer and build/run the app; smoke test with sequences for promoter/non-promoter and strong/weak.
```
