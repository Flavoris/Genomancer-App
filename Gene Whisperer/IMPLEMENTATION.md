# Gene Whisperer Implementation Guide (Transformer, 2-Stage)

Step-by-step plan for adding the Gene Whisperer ML feature (PyTorch -> CoreML) to Genomancer with a shared Transformer backbone and two classification stages: (1) promoter detection; (2) strong vs weak promoters. All training code and artifacts live under `Gene Whisperer/`, with only the compiled `.mlmodel/.mlpackage` files and Swift glue added to the iOS app.

## 0) Prereqs
- Python 3.10+ with PyTorch (CPU is fine), `coremltools`, `pandas`, `numpy`, `scikit-learn`, `pyyaml`, `tqdm`.
- Xcode 15+.

## 1) Directory Scaffold (in repo root)
```
Gene Whisperer/
  data/
    stage1_train.csv   # <-- YOU provide
    stage1_val.csv     # <-- YOU provide
    stage2_train.csv   # <-- YOU provide
    stage2_val.csv     # <-- YOU provide
    stage2_test.csv    # optional
  training/
    config.yaml
    dataset.py
    model.py
    train_stage1.py
    train_stage2.py
    inference.py
    export_coreml.py
  artifacts/
    checkpoints/
      stage1_best.pt
      stage2_best.pt
    GeneWhispererPromoter.mlpackage / .mlmodel
    GeneWhispererStrength.mlpackage / .mlmodel
```

## 2) Provide Datasets (you do this before training)
- Stage 1 (promoter vs non-promoter): TSV/CSV with columns `sequence` and `is_promoter` where `is_promoter` is 0/1.
- Stage 2 (strong vs weak promoters): TSV/CSV with columns `sequence` and `strength`, where `strength` is 1 for strong and 0 for weak (only promoter sequences). If you prefer string labels, adjust the config mapping.
- Sequences: uppercase A/C/G/T/N only. Length ~81 bp (variable allowed). Any N/invalid tokens will be treated as unknown.
- Drop your files into `Gene Whisperer/data/` with the names above (or update `config.yaml` paths). Set the delimiter in `config.yaml` if you use TSV.

## 3) Python Env
- From repo root:
  - `python -m venv Gene\ Whisperer/.venv`
  - `source Gene\ Whisperer/.venv/bin/activate`
  - `pip install torch torchvision torchaudio coremltools pandas numpy scikit-learn pyyaml tqdm`

## 4) Config
- Create `Gene Whisperer/training/config.yaml`:
```yaml
max_bp_len: 81           # bp length cap; tokens = max_bp_len - 2 after 3-mer
batch_size: 64
lr: 0.001
weight_decay: 0.0005
epochs: 30
early_stopping_patience: 5
train_val_split: 0.9     # used only if val is missing
delimiter: "\t"          # set to "," for CSV
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

## 5) Tokenization & Features
- Add `Gene Whisperer/training/dataset.py`:
  - Build 64-size 3-mer vocab (AAA...TTT) with deterministic index order.
  - `tokenize_kmers(seq, max_bp_len)`: uppercase, strip, pad/truncate bp to `max_bp_len`, create overlapping 3-mers, map to vocab indices, pad/truncate tokens to `max_tokens = max_bp_len-2`.
  - `compute_tnc(seq)`: 64-dim normalized tri-nucleotide frequency vector.
  - `compute_pseeiip(seq)`: 64-dim PseEIIP vector using published EIIP values; unknowns -> zero.
  - `compute_cksnap(seq, max_gap=5)`: 96-dim k-spaced dinucleotide composition (16 pairs × 6 gaps).
  - `compute_pstnp(seq)`: 64-dim position-weighted trinucleotide propensity (emphasize -10/-35 regions).
  - TSV/CSV loading: respect `delimiter` in config.
  - `PromoterDatasetStage1`: yields (tokens_tensor[seq_len], engineered[288], label_float) for `is_promoter` 0/1.
  - `PromoterDatasetStage2`: yields (tokens_tensor[seq_len], engineered[288], label_float) for `strength` using `stage2_strength_positive`/`stage2_strength_negative` mapping (e.g., 1=strong, 0=weak).
  - `build_dataloaders(cfg)`: DataLoaders for both stages; supports random split if val missing; uses class-weighted sampling when imbalance >20%.

## 6) Model
- Add `Gene Whisperer/training/model.py`:
  - Shared Transformer encoder (`TransformerBackbone`) with:
    - 3-mer embedding dim 48.
    - 4 encoder layers; each: LayerNorm -> Multi-Head Self-Attention (num_heads=4, hidden_dim=48, dropout=0.1) -> residual; LayerNorm -> FFN (64 -> 48) -> residual.
    - Positional embeddings (learned) sized to `max_tokens`.
  - Stage 1 head (`PromoterClassifier`):
    - Input: contextual embeddings.
    - 1D Conv (in_channels=48, out_channels=128, kernel_size=5, padding=2) + ReLU.
    - BiLSTM with 128 units per direction.
    - Global max pool.
    - Dense 128 + ReLU -> Dropout 0.5 -> Dense 1 (logits).
  - Stage 2 head (`StrengthClassifier`):
    - Input: contextual embeddings.
    - 1D Conv (48 -> 96, kernel_size=5, padding=2) + ReLU.
    - Uni-directional LSTM 128 units.
    - Global average pool.
    - Concatenate with TNC + PseEIIP (64-dim each) to form a single feature vector.
    - Dense 128 + ReLU -> Dropout 0.5 -> Dense 1 (logits).
  - `GeneWhispererStage1` and `GeneWhispererStage2` wrappers share the backbone; Stage 2 exposes `freeze_lower_layers()` to freeze bottom half of encoder layers before training.

## 7) Training Stage 1 (Promoter detection)
- Add `Gene Whisperer/training/train_stage1.py`:
  - CLI: `python train_stage1.py --config config.yaml`.
  - Loads cfg, builds backbone+Stage1 head, dataloaders.
  - Loss: BCEWithLogitsLoss; metrics: accuracy, precision, recall, MCC.
  - Early stopping on val loss or MCC; save best to `../artifacts/checkpoints/stage1_best.pt`.
  - Save `../artifacts/stage1_report.json` with metrics, class distribution.

## 8) Training Stage 2 (Strong vs Weak)
- Add `Gene Whisperer/training/train_stage2.py`:
  - CLI: `python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt`.
  - Loads backbone with Stage1 weights, freezes bottom half of Transformer layers (e.g., layers 0-1), keeps embeddings trainable or frozen per a flag.
  - Uses Stage2 dataset with aux features.
  - Loss: BCEWithLogitsLoss; metrics: accuracy, MCC, ROC-AUC.
  - Early stopping; save best to `../artifacts/checkpoints/stage2_best.pt` and `../artifacts/stage2_report.json`.

## 9) Inference Helper (Python)
- Add `Gene Whisperer/training/inference.py`:
  - Functions: `predict_promoter(sequence) -> (score, is_promoter)` using Stage1 model; `predict_strength(sequence) -> (score, is_strong)` using Stage2 model (assumes promoter).
  - CLI: `python inference.py --sequence ACTG... --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt --stage2_ckpt ../artifacts/checkpoints/stage2_best.pt --config config.yaml`.
  - If Stage1 is negative, skip Stage2 and print that it was skipped.

## 10) Export to CoreML
- Add `Gene Whisperer/training/export_coreml.py`:
  - Loads `stage1_best.pt` and `stage2_best.pt`.
  - Uses dummy inputs: token tensor `(1, max_tokens)` plus engineered features `(1, 128)` for Stage1/Stage2.
  - Converts to two models via `coremltools.convert`:
    - `GeneWhispererPromoter.mlmodel`/`.mlpackage`
    - `GeneWhispererStrength.mlmodel`/`.mlpackage`
  - Sets metadata (author, short descriptions, task names, input descriptions: tokens and features) and class labels (`promoter`/`non-promoter`, `strong`/`weak`).
  - Saves into `Gene Whisperer/artifacts/`.

## 11) iOS Integration
- Bundle both exported models into `apps/ios/Genomancer/Genomancer/Resources/GeneWhisperer/` and add to the Xcode target.
- Add Swift glue `GeneWhispererPredictor.swift`:
  - Shared 3-mer tokenizer mirroring Python (pad/truncate to `max_bp_len`, overlapping 3-mers -> indices).
  - Feature builders for TNC and PseEIIP (64-dim each); ensure parity with Python.
  - `predictPromoter(sequence:String) -> (isPromoter: Bool, score: Double)`.
  - `predictStrength(sequence:String) -> (isStrong: Bool, score: Double)` that calls Stage2 model; caller should skip if Stage1 is negative.
- Add UI `GeneWhispererView.swift`:
  - TextEditor for sequence input with validation (A/C/G/T/N, length <= max_bp_len).
  - Button runs Stage1; if positive, runs Stage2; show scores and labels.
  - Inline validation/errors; lightweight loading state.
- Wire navigation from `HomeView` to the new screen (“Gene Whisperer (AI)”).

## 12) Testing/Verification
- Python: unit tests for tokenizer/feature parity (vs Swift reference vectors), quick overfit test on a tiny batch, inference skip-logic when Stage1 is negative.
- iOS: simulator test with sample sequences covering promoter/non-promoter and strong/weak; verify model load and feature computation parity.

## 13) Quick command sequence (after files and data are in place)
```
cd "Gene Whisperer"
source .venv/bin/activate
cd training
python train_stage1.py --config config.yaml
python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt
python export_coreml.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt --stage2_ckpt ../artifacts/checkpoints/stage2_best.pt
```
- Then copy both exported models into the iOS Resources path and run the app.

## 14) Deliverables to expect
- Training code under `Gene Whisperer/training/`.
- Artifacts under `Gene Whisperer/artifacts/` (checkpoints + two CoreML models).
- Swift files: predictor and UI for Gene Whisperer, plus navigation entry.
- Updated instructions for dataset drop-in locations.
