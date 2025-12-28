# Colab Quickstart

Run Gene Whisperer training in Google Colab with two cells.

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
```

---

## Cell 2: Run Training

### Option A: MLM Pretraining

```python
!bash "Gene Whisperer/scripts/colab_run_mlm.sh" --run_name my_first_run
```

### Option B: Stage 1 Training (Promoter vs Non-Promoter)

```python
!bash "Gene Whisperer/scripts/colab_run_stage1.sh" --run_name stage1_run
```

With a pretrained MLM encoder:

```python
!bash "Gene Whisperer/scripts/colab_run_stage1.sh" \
    --run_name stage1_run \
    --kmer 4 \
    --stage1_ckpt "Gene Whisperer/artifacts/mlm_encoder_k4.pt"
```

### Option C: Stage 2 Training (Strong vs Weak Promoters)

```python
!bash "Gene Whisperer/scripts/colab_run_stage2.sh" --run_name stage2_run
```

With a pretrained Stage 1 checkpoint:

```python
!bash "Gene Whisperer/scripts/colab_run_stage2.sh" \
    --run_name stage2_run \
    --kmer 4 \
    --stage1_ckpt "Gene Whisperer/artifacts/checkpoints/stage1_k4.pt"
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

## Script Arguments

All three scripts (`colab_run_mlm.sh`, `colab_run_stage1.sh`, `colab_run_stage2.sh`) accept:

| Argument | Description | Default |
|----------|-------------|---------|
| `--run_name` | Name for the run directory | `<stage>_YYYYMMDD_HHMMSS` |
| `--drive_dir` | Google Drive output directory | `/content/drive/MyDrive/GeneWhispererRuns` |
| `--config` | Path to config YAML | `Gene Whisperer/training/config.yaml` |
| `--kmer` | K-mer size (3, 4, or 6) | From config |

Stage-specific arguments:

| Script | Argument | Description |
|--------|----------|-------------|
| `colab_run_stage1.sh` | `--stage1_ckpt` | MLM encoder checkpoint for initialization |
| `colab_run_stage2.sh` | `--stage1_ckpt` | Stage 1 checkpoint for initialization |

---

## Troubleshooting

### Missing GPU

If you see "No NVIDIA GPU detected":

1. Go to **Runtime > Change runtime type**
2. Set **Hardware accelerator** to **GPU** (T4 or higher)
3. Click **Save** and re-run both cells

### Requirements Path Issues (Spaces in Folder Name)

The folder is named `Gene Whisperer` (with a space). If you get path errors:

- Ensure you're using quoted paths: `"Gene Whisperer/scripts/..."`
- Don't rename the folder; the scripts expect this exact name
- If cloning manually, keep the repo structure intact

### Resume Training

All scripts look for existing checkpoints in your Drive run directory. Currently, the training scripts do not support `--resume` directly.

To continue from a checkpoint manually:

1. Note the checkpoint path from your previous run
2. For Stage 1/2: pass the checkpoint via `--stage1_ckpt`
3. Checkpoints from interrupted runs are automatically saved to Drive

When `--resume` support is added, scripts will automatically resume from the latest checkpoint found in `<drive_dir>/<run_name>/checkpoints/`.
