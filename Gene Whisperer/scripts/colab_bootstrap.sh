#!/usr/bin/env bash
# chmod +x scripts/colab_bootstrap.sh
#
# Colab Bootstrap Script for Gene Whisperer
# Sets up the environment for running Gene Whisperer in Google Colab.
#
# Usage (in Colab):
#   !bash "Gene Whisperer/scripts/colab_bootstrap.sh"
#
# This script is idempotent - safe to run multiple times.

set -euo pipefail

# -----------------------------------------------------------------------------
# Colab Detection
# -----------------------------------------------------------------------------
is_colab() {
    # Check for COLAB_RELEASE_TAG environment variable or /content directory
    if [[ -n "${COLAB_RELEASE_TAG:-}" ]] || [[ -d "/content" ]]; then
        return 0
    fi
    return 1
}

if ! is_colab; then
    echo "============================================================"
    echo "Not running in Google Colab."
    echo ""
    echo "This script is designed for Google Colab environments."
    echo "If you're running locally, please follow the standard setup:"
    echo "  1. pip install -r 'Gene Whisperer/requirements.txt'"
    echo "  2. cd 'Gene Whisperer/training'"
    echo "  3. python train_stage1.py --config config.yaml"
    echo "============================================================"
    exit 0
fi

echo "============================================================"
echo "Gene Whisperer - Colab Bootstrap"
echo "============================================================"
echo ""

# -----------------------------------------------------------------------------
# Navigate to Repository Root
# -----------------------------------------------------------------------------
echo "[1/8] Navigating to repository root..."

# Get the git root directory
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"

if [[ -z "$REPO_ROOT" ]]; then
    echo "ERROR: Not inside a git repository."
    echo "Please clone the repository first:"
    echo "  !git clone <your-repo-url>"
    exit 1
fi

cd "$REPO_ROOT"
echo "  Repository root: $REPO_ROOT"
echo ""

# -----------------------------------------------------------------------------
# Install Git LFS and Pull Large Files (Optional)
# -----------------------------------------------------------------------------
echo "[2/8] Setting up Git LFS..."

SKIP_LFS_PULL="${SKIP_LFS_PULL:-0}"

if [[ "$SKIP_LFS_PULL" == "1" ]]; then
    echo "  SKIP_LFS_PULL=1 -> skipping Git LFS setup."
else
    if grep -q "filter=lfs" "$REPO_ROOT/.gitattributes" 2>/dev/null; then
        # Check if git-lfs is installed
        if ! command -v git-lfs &> /dev/null; then
            echo "  Installing Git LFS..."
            apt-get update -qq && apt-get install -qq -y git-lfs
        fi

        # Initialize LFS and pull actual file content
        git lfs install --skip-smudge 2>/dev/null || true
        echo "  Pulling LFS files (genome data, checkpoints)..."
        if ! git lfs pull; then
            echo "  WARNING: Git LFS pull failed; continuing without LFS files."
            echo "  If you need weights, copy them from Drive into Gene Whisperer/artifacts."
        fi
        echo "  Git LFS setup complete."
    else
        echo "  No LFS filters configured; skipping Git LFS pull."
    fi
fi
echo ""

# -----------------------------------------------------------------------------
# Install Requirements
# -----------------------------------------------------------------------------
echo "[3/8] Installing Python requirements..."

REQUIREMENTS_PATH="$REPO_ROOT/Gene Whisperer/requirements.txt"

if [[ -f "$REQUIREMENTS_PATH" ]]; then
    pip install -q -r "$REQUIREMENTS_PATH"
    echo "  Requirements installed from: $REQUIREMENTS_PATH"
else
    echo "  WARNING: requirements.txt not found at: $REQUIREMENTS_PATH"
    echo "  Attempting to install core dependencies manually..."
    pip install -q torch numpy pandas scikit-learn pyyaml tqdm matplotlib
fi
echo ""

# -----------------------------------------------------------------------------
# Create Required Directories
# -----------------------------------------------------------------------------
echo "[4/8] Creating required directories..."

GW_DIR="$REPO_ROOT/Gene Whisperer"

declare -a DIRS=(
    "$GW_DIR/artifacts"
    "$GW_DIR/artifacts/checkpoints"
    "$GW_DIR/artifacts/vocabs"
    "$GW_DIR/data"
)

for dir in "${DIRS[@]}"; do
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        echo "  Created: $dir"
    else
        echo "  Exists:  $dir"
    fi
done
echo ""

# -----------------------------------------------------------------------------
# Download Large Genome Files from Drive
# -----------------------------------------------------------------------------
echo "[5/8] Downloading large genome files from Google Drive..."

GENOME_DOWNLOAD_SCRIPT="$REPO_ROOT/Gene Whisperer/scripts/colab_download_genomes.sh"

if [[ -f "$GENOME_DOWNLOAD_SCRIPT" ]]; then
    bash "$GENOME_DOWNLOAD_SCRIPT"
else
    echo "  WARNING: Genome download script not found: $GENOME_DOWNLOAD_SCRIPT"
    echo "  Large genome files (human_genome, japanese_rice_genome) may not be available."
fi
echo ""

# -----------------------------------------------------------------------------
# Verify Genome Files
# -----------------------------------------------------------------------------
echo "[6/8] Verifying genome files..."
echo ""

DATA_DIR="$GW_DIR/data"

# Expected genome files for MLM pretraining
declare -a EXPECTED_GENOMES=(
    "GCF_000005845.2_ASM584v2_genomic.fna"
    "B_subtilis_genome"
    "S_cerevisiae_genome"
    "human_genome"
    "pseudomonas_genome"
    "japanese_rice_genome"
)

echo "  Genome file status:"
echo "  -------------------"
MISSING_COUNT=0
for genome in "${EXPECTED_GENOMES[@]}"; do
    GENOME_PATH="$DATA_DIR/$genome"
    if [[ -f "$GENOME_PATH" ]]; then
        SIZE=$(ls -lh "$GENOME_PATH" 2>/dev/null | awk '{print $5}')
        echo "  [OK]      $genome ($SIZE)"
    else
        echo "  [MISSING] $genome"
        MISSING_COUNT=$((MISSING_COUNT + 1))
    fi
done
echo ""

if [[ $MISSING_COUNT -gt 0 ]]; then
    echo "  WARNING: $MISSING_COUNT genome file(s) missing!"
    echo "  MLM pretraining will only use available genomes."
    echo "  To add missing genomes, upload them to Google Drive:"
    echo "    /MyDrive/GeneWhispererData/"
else
    echo "  All expected genome files are present."
fi
echo ""

# -----------------------------------------------------------------------------
# GPU Status
# -----------------------------------------------------------------------------
echo "[7/8] GPU Status..."

if command -v nvidia-smi &> /dev/null; then
    nvidia-smi
else
    echo "  No NVIDIA GPU detected."
fi
echo ""

# -----------------------------------------------------------------------------
# Environment Info
# -----------------------------------------------------------------------------
echo "[8/8] Environment Information..."

echo "  Python version: $(python --version 2>&1)"

python -c "
import torch
print(f'  PyTorch version: {torch.__version__}')
print(f'  CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'  CUDA version: {torch.version.cuda}')
    print(f'  GPU device: {torch.cuda.get_device_name(0)}')
"
echo ""

# -----------------------------------------------------------------------------
# Training Commands
# -----------------------------------------------------------------------------
echo "============================================================"
echo "Setup Complete! Training Commands:"
echo "============================================================"
echo ""
echo "MLM Pretraining:"
echo "  cd \"$GW_DIR/training\" && python pretrain_mlm.py --config config.yaml"
echo ""
echo "Stage 1 Training:"
echo "  cd \"$GW_DIR/training\" && python train_stage1.py --config config.yaml"
echo ""
echo "Stage 2 Training:"
echo "  cd \"$GW_DIR/training\" && python train_stage2.py --config config.yaml"
echo ""
echo "============================================================"
echo "Tip: To use a specific checkpoint for Stage 2:"
echo "  python train_stage2.py --config config.yaml --stage1_ckpt ../artifacts/checkpoints/stage1_best.pt"
echo "============================================================"
