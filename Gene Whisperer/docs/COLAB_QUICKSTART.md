# Colab Quickstart

Run Gene Whisperer training in Google Colab with GPU acceleration.

---

## Cell 1: Setup

Mount Drive, clone the repo, and run the bootstrap script.

```python
# Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

# Clone the repository
!git clone https://github.com/YOUR_USERNAME/Genomancer.git /content/Genomancer

# Run bootstrap (installs requirements, checks GPU)
%cd /content/Genomancer
!bash "Gene Whisperer/scripts/colab_bootstrap.sh"

# If you store weights on Drive and want to skip Git LFS downloads:
# !SKIP_LFS_PULL=1 bash "Gene Whisperer/scripts/colab_bootstrap.sh"
```

---

## Training Pipeline Overview

Gene Whisperer now uses **BPE tokenization** (DNABERT-2 style) instead of k-mer tokenization. The training pipeline has 4 steps:

```
1. Train BPE Vocabulary → 2. MLM Pretraining → 3. Stage 1 Training → 4. Stage 2 Training
        (one-time)           (self-supervised)    (promoter vs non)    (strong vs weak)
```

---

## Cell 2: Train BPE Vocabulary (One-Time)

**Run this once** to create the BPE vocabulary from your genomic FASTA data.

```python
%cd /content/Genomancer/"Gene Whisperer"/training

# Train BPE vocabulary (takes ~5-10 minutes)
!python train_bpe_vocab.py --config config.yaml
```

This creates `artifacts/vocabs/bpe_vocab.json` with 4096 tokens optimized for DNA sequences.

**Options:**
```python
# Custom vocab size (default: 4096, recommended by DNABERT-2)
!python train_bpe_vocab.py --config config.yaml --vocab-size 8192

# Custom FASTA data
!python train_bpe_vocab.py --fasta /path/to/genomes.fna --output custom_vocab.json

# Limit training windows (faster, for testing)
!python train_bpe_vocab.py --config config.yaml --max-windows 100000
```

---

## Cell 3: MLM Pretraining

Self-supervised masked language model pretraining on genomic sequences.

### Option A: Using Colab Script (Recommended)

```python
!bash "Gene Whisperer/scripts/colab_run_mlm.sh" --run_name mlm_bpe_run
```

### Option B: Direct Python Command

```python
%cd /content/Genomancer/"Gene Whisperer"/training

# Run MLM pretraining with BPE
!python pretrain_mlm.py --config config.yaml
```

**Output:** Creates `artifacts/mlm_encoder_bpe.pt` — the pretrained encoder weights.

---

## Cell 4: Stage 1 Training (Promoter vs Non-Promoter)

Binary classification: is this sequence a promoter?

### Option A: Using Colab Script (Recommended)

```python
!bash "Gene Whisperer/scripts/colab_run_stage1.sh" --run_name stage1_bpe_run
```

### Option B: Direct Python Command

```python
%cd /content/Genomancer/"Gene Whisperer"/training

# Train Stage 1 classifier
!python train_stage1.py --config config.yaml
```

**With a pretrained MLM encoder:**

The config automatically looks for `artifacts/mlm_encoder_bpe.pt`. To use a custom checkpoint:

```python
# Set environment variable for custom MLM checkpoint
import os
os.environ['MLM_ENCODER_CHECKPOINT'] = '/path/to/mlm_encoder_bpe.pt'
!python train_stage1.py --config config.yaml
```

**Output:** Creates `artifacts/checkpoints/stage1_bpe.pt`

---

## Cell 5: Stage 2 Training (Strong vs Weak Promoters)

Binary classification: is this promoter strong or weak?

### Option A: Using Colab Script (Recommended)

```python
!bash "Gene Whisperer/scripts/colab_run_stage2.sh" --run_name stage2_bpe_run
```

### Option B: Direct Python Command

```python
%cd /content/Genomancer/"Gene Whisperer"/training

# Train Stage 2 classifier
!python train_stage2.py --config config.yaml
```

**With a Stage 1 checkpoint:**

```python
!python train_stage2.py --config config.yaml \
    --stage1_ckpt "Gene Whisperer/artifacts/checkpoints/stage1_bpe.pt"
```

**Output:** Creates `artifacts/checkpoints/stage2_bpe.pt`

---

## Running Inference

After training, run inference on new sequences:

```python
%cd /content/Genomancer/"Gene Whisperer"/training

# Single sequence inference
!python ensemble_infer.py --sequence "ATGCATGCATGC..." --config config.yaml

# With Test-Time Augmentation (TTA) for better accuracy
!python ensemble_infer.py --sequence "ATGCATGCATGC..." --config config.yaml --tta
```

---

## Where Artifacts Appear

All outputs are saved to Google Drive at:

```
/content/drive/MyDrive/GeneWhispererRuns/<run_name>/
├── artifacts/       # Vocabularies, reports, encoder weights
├── checkpoints/     # Model checkpoints (.pt files)
└── logs/            # Training logs
```

Default location: `/content/drive/MyDrive/GeneWhispererRuns/`

To use a custom location:
```python
!bash "Gene Whisperer/scripts/colab_run_mlm.sh" --run_name my_run --drive_dir /content/drive/MyDrive/CustomFolder
```

---

## Key Configuration (config.yaml)

The BPE pipeline uses these settings in `config.yaml`:

```yaml
# BPE Tokenization
bpe_vocab_path: ../artifacts/vocabs/bpe_vocab.json
vocab_size: 4096              # DNABERT-2 optimal
pad_token_id: 0               # BPE pad token
max_token_len: 24             # For 81bp sequences (~4x compression)
mlm_max_token_len: 64         # For 234bp MLM sequences

# MLM Pretraining
mlm_encoder_checkpoint: ../artifacts/mlm_encoder_bpe.pt
mlm_max_span_len: 1           # Independent token masking for BPE

# Checkpoints
stage1_checkpoint_name: stage1_bpe.pt
```

---

## Script Arguments

All Colab scripts accept:

| Argument | Description | Default |
|----------|-------------|---------|
| `--run_name` | Name for the run directory | `<stage>_YYYYMMDD_HHMMSS` |
| `--drive_dir` | Google Drive output directory | `/content/drive/MyDrive/GeneWhispererRuns` |
| `--config` | Path to config YAML | `Gene Whisperer/training/config.yaml` |

Stage-specific arguments:

| Script | Argument | Description |
|--------|----------|-------------|
| `colab_run_stage1.sh` | `--stage1_ckpt` | MLM encoder checkpoint for initialization |
| `colab_run_stage2.sh` | `--stage1_ckpt` | Stage 1 checkpoint for initialization |

---

## Full Training Example

Complete Colab notebook cells to train from scratch:

```python
# Cell 1: Setup
from google.colab import drive
drive.mount('/content/drive')
!git clone https://github.com/YOUR_USERNAME/Genomancer.git /content/Genomancer
%cd /content/Genomancer
!bash "Gene Whisperer/scripts/colab_bootstrap.sh"
```

```python
# Cell 2: Train BPE Vocabulary (run once)
%cd /content/Genomancer/"Gene Whisperer"/training
!python train_bpe_vocab.py --config config.yaml
```

```python
# Cell 3: MLM Pretraining
!python pretrain_mlm.py --config config.yaml
```

```python
# Cell 4: Stage 1 Training
!python train_stage1.py --config config.yaml
```

```python
# Cell 5: Stage 2 Training
!python train_stage2.py --config config.yaml
```

```python
# Cell 6: Test Inference
!python ensemble_infer.py \
    --sequence "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" \
    --config config.yaml --tta
```

---

## Troubleshooting

### Missing GPU

If you see "No NVIDIA GPU detected":

1. Go to **Runtime > Change runtime type**
2. Set **Hardware accelerator** to **GPU** (T4 or higher)
3. Click **Save** and re-run both cells

### BPE Vocabulary Not Found

If you see `FileNotFoundError: BPE vocabulary not found`:

```python
# Train the BPE vocabulary first
%cd /content/Genomancer/"Gene Whisperer"/training
!python train_bpe_vocab.py --config config.yaml
```

### Requirements Path Issues (Spaces in Folder Name)

The folder is named `Gene Whisperer` (with a space). If you get path errors:

- Ensure you're using quoted paths: `"Gene Whisperer/scripts/..."`
- Don't rename the folder; the scripts expect this exact name
- If cloning manually, keep the repo structure intact

### Git LFS Quota Errors

If you see LFS errors or checkout failures, rerun bootstrap with:

```python
!SKIP_LFS_PULL=1 bash "Gene Whisperer/scripts/colab_bootstrap.sh"
```

Then copy your Drive checkpoints into `Gene Whisperer/artifacts/` as needed.

### Small Reference Genome Files

If `B_subtilis_genome` or `S_cerevisiae_genome` show up as ~100-byte files, they were checked out as placeholders.
The bootstrap script now repairs these by re-downloading clean FASTA files. To skip that step:

```python
!SKIP_REFERENCE_GENOME_FIX=1 bash "Gene Whisperer/scripts/colab_bootstrap.sh"
```

---

## Migration from K-mer to BPE

If you have existing k-mer checkpoints (e.g., `mlm_encoder_k6.pt`, `stage1_k6.pt`), note that:

- **BPE and k-mer checkpoints are NOT compatible** — you cannot load k-mer weights into a BPE model
- BPE provides better performance with ~5x sequence compression
- Retraining with BPE is recommended for best results

Key changes:
- `--kmer` argument is removed (BPE uses variable-length tokens)
- Vocab size: 4099 → 4096
- Pad token ID: 4098 → 0
- Checkpoint names: `*_k6.pt` → `*_bpe.pt`
