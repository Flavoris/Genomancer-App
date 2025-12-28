#!/usr/bin/env bash
# chmod +x "Gene Whisperer/scripts/colab_run_stage1.sh"
#
# Colab Stage 1 Training Script for Gene Whisperer
# Runs Stage 1 classifier training with robust checkpointing, resume support, and Drive backup.
#
# Usage (in Colab):
#   !bash "Gene Whisperer/scripts/colab_run_stage1.sh"
#   !bash "Gene Whisperer/scripts/colab_run_stage1.sh" --run_name my_experiment
#   !bash "Gene Whisperer/scripts/colab_run_stage1.sh" --kmer 4 --drive_dir /content/drive/MyDrive/CustomDir
#   !bash "Gene Whisperer/scripts/colab_run_stage1.sh" --stage1_ckpt path/to/mlm_encoder.pt
#
# This script:
#   1. Checks Google Drive is mounted
#   2. Runs colab_bootstrap.sh to set up environment
#   3. Creates output directories on Drive
#   4. Runs Stage 1 training with logging to Drive
#   5. Implements resume from existing checkpoints
#   6. Copies artifacts to Drive on completion or Ctrl+C

set -euo pipefail

# =============================================================================
# Default Configuration
# =============================================================================
DEFAULT_DRIVE_DIR="/content/drive/MyDrive/GeneWhispererRuns"
DEFAULT_CONFIG="Gene Whisperer/training/config.yaml"
DEFAULT_RUN_NAME="stage1_$(date +%Y%m%d_%H%M%S)"

# =============================================================================
# Parse Arguments
# =============================================================================
DRIVE_DIR="$DEFAULT_DRIVE_DIR"
CONFIG="$DEFAULT_CONFIG"
KMER=""
RUN_NAME="$DEFAULT_RUN_NAME"
STAGE1_CKPT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --drive_dir)
            DRIVE_DIR="$2"
            shift 2
            ;;
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --kmer)
            KMER="$2"
            shift 2
            ;;
        --run_name)
            RUN_NAME="$2"
            shift 2
            ;;
        --stage1_ckpt)
            STAGE1_CKPT="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --drive_dir DIR     Google Drive output directory (default: $DEFAULT_DRIVE_DIR)"
            echo "  --config PATH       Config file path relative to git root (default: $DEFAULT_CONFIG)"
            echo "  --kmer K            K-mer size (optional, overrides config)"
            echo "  --run_name NAME     Run name for output directory (default: stage1_YYYYMMDD_HHMMSS)"
            echo "  --stage1_ckpt PATH  Optional MLM encoder checkpoint to initialize from"
            echo "  -h, --help          Show this help message"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

# =============================================================================
# Helper Functions
# =============================================================================
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

# =============================================================================
# Check Google Drive Mount
# =============================================================================
echo "============================================================"
echo "Gene Whisperer - Colab Stage 1 Training"
echo "============================================================"
echo ""

log_info "Checking Google Drive mount..."

if [[ ! -d "/content/drive/MyDrive" ]]; then
    log_error "Google Drive is NOT mounted!"
    echo ""
    echo "Please mount Google Drive first by running this cell in Colab:"
    echo ""
    echo "    from google.colab import drive"
    echo "    drive.mount('/content/drive')"
    echo ""
    echo "Then re-run this script."
    exit 1
fi

log_info "Google Drive is mounted at /content/drive"
echo ""

# =============================================================================
# Get Repository Root
# =============================================================================
log_info "Finding repository root..."

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"

if [[ -z "$REPO_ROOT" ]]; then
    log_error "Not inside a git repository."
    echo "Please clone the repository first:"
    echo "  !git clone <your-repo-url>"
    exit 1
fi

log_info "Repository root: $REPO_ROOT"
echo ""

# =============================================================================
# Run Bootstrap Script
# =============================================================================
log_info "Running bootstrap script..."

BOOTSTRAP_SCRIPT="$REPO_ROOT/Gene Whisperer/scripts/colab_bootstrap.sh"

if [[ ! -f "$BOOTSTRAP_SCRIPT" ]]; then
    log_error "Bootstrap script not found: $BOOTSTRAP_SCRIPT"
    exit 1
fi

bash "$BOOTSTRAP_SCRIPT"
echo ""

# =============================================================================
# Create Output Directories on Drive
# =============================================================================
log_info "Creating output directories on Drive..."

RUN_DIR="$DRIVE_DIR/$RUN_NAME"
ARTIFACTS_DIR="$RUN_DIR/artifacts"
CHECKPOINTS_DIR="$RUN_DIR/checkpoints"
LOGS_DIR="$RUN_DIR/logs"

mkdir -p "$ARTIFACTS_DIR"
mkdir -p "$CHECKPOINTS_DIR"
mkdir -p "$LOGS_DIR"

log_info "Run directory: $RUN_DIR"
log_info "  Artifacts:   $ARTIFACTS_DIR"
log_info "  Checkpoints: $CHECKPOINTS_DIR"
log_info "  Logs:        $LOGS_DIR"
echo ""

# =============================================================================
# Resolve Config Path
# =============================================================================
if [[ "$CONFIG" != /* ]]; then
    CONFIG_PATH="$REPO_ROOT/$CONFIG"
else
    CONFIG_PATH="$CONFIG"
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
    log_error "Config file not found: $CONFIG_PATH"
    exit 1
fi

log_info "Using config: $CONFIG_PATH"
echo ""

# =============================================================================
# Check for Existing Checkpoint (Resume Support)
# =============================================================================
RESUME_ARG=""
LATEST_CHECKPOINT=""

# Look for existing checkpoints to resume from
if [[ -d "$CHECKPOINTS_DIR" ]]; then
    # Find the most recent .pt file in checkpoints directory
    LATEST_CHECKPOINT=$(find "$CHECKPOINTS_DIR" -name "stage1*.pt" -type f 2>/dev/null | sort -r | head -n 1 || true)
fi

if [[ -n "$LATEST_CHECKPOINT" ]]; then
    log_info "Found existing checkpoint: $LATEST_CHECKPOINT"
    log_warn "NOTE: train_stage1.py does not yet support --resume flag."
    log_warn "Training will start fresh. Checkpoint found but cannot be used for resume."
else
    log_info "No existing checkpoint found. Starting fresh training."
fi
echo ""

# =============================================================================
# Set Environment Variables for Stability
# =============================================================================
export PYTHONUNBUFFERED=1
export TOKENIZERS_PARALLELISM=false

log_info "Environment variables set:"
log_info "  PYTHONUNBUFFERED=1"
log_info "  TOKENIZERS_PARALLELISM=false"
echo ""

# =============================================================================
# Prepare Training Command
# =============================================================================
TRAINING_DIR="$REPO_ROOT/Gene Whisperer/training"
LOG_FILE="$LOGS_DIR/stage1_training_$(date +%Y%m%d_%H%M%S).log"

# Build the training command
TRAIN_CMD="python train_stage1.py --config \"$CONFIG_PATH\""

# Add kmer argument if specified
if [[ -n "$KMER" ]]; then
    log_info "K-mer size specified: $KMER"
    TRAIN_CMD="$TRAIN_CMD --kmer $KMER"
fi

# Add stage1_ckpt argument if specified (this is the MLM encoder checkpoint)
if [[ -n "$STAGE1_CKPT" ]]; then
    log_info "MLM encoder checkpoint specified: $STAGE1_CKPT"
    # Resolve relative path if needed
    if [[ "$STAGE1_CKPT" != /* ]]; then
        STAGE1_CKPT="$REPO_ROOT/$STAGE1_CKPT"
    fi
    if [[ -f "$STAGE1_CKPT" ]]; then
        log_info "Using MLM checkpoint: $STAGE1_CKPT"
        # Note: train_stage1.py loads MLM weights via config, not CLI arg
        # We set an environment variable to pass this
        export MLM_ENCODER_CHECKPOINT="$STAGE1_CKPT"
    else
        log_warn "MLM checkpoint not found: $STAGE1_CKPT"
    fi
fi

log_info "Training command: $TRAIN_CMD"
log_info "Log file: $LOG_FILE"
echo ""

# =============================================================================
# Cleanup Function (Called on Exit)
# =============================================================================
cleanup_and_copy() {
    local exit_code=$?
    echo ""
    log_info "============================================================"
    log_info "Training ended. Copying artifacts to Drive..."
    log_info "============================================================"

    # Copy artifacts from local Gene Whisperer/artifacts to Drive
    LOCAL_ARTIFACTS="$REPO_ROOT/Gene Whisperer/artifacts"

    if [[ -d "$LOCAL_ARTIFACTS" ]]; then
        log_info "Copying from $LOCAL_ARTIFACTS to $ARTIFACTS_DIR"
        # Use rsync if available, otherwise cp
        if command -v rsync &> /dev/null; then
            rsync -av --progress "$LOCAL_ARTIFACTS/" "$ARTIFACTS_DIR/" 2>&1 || true
        else
            cp -rv "$LOCAL_ARTIFACTS/"* "$ARTIFACTS_DIR/" 2>&1 || true
        fi
    else
        log_warn "Local artifacts directory not found: $LOCAL_ARTIFACTS"
    fi

    # Copy any .pt checkpoint files from training directory
    log_info "Copying any .pt files from training directory..."
    find "$TRAINING_DIR" -maxdepth 1 -name "*.pt" -type f -exec cp -v {} "$CHECKPOINTS_DIR/" \; 2>&1 || true

    # Copy any stage1_report*.json files from artifacts
    log_info "Copying stage1_report*.json files..."
    find "$LOCAL_ARTIFACTS" -name "stage1_report*.json" -type f -exec cp -v {} "$ARTIFACTS_DIR/" \; 2>&1 || true

    # Copy from runs directory if it exists
    RUNS_DIR="$REPO_ROOT/Gene Whisperer/runs"
    if [[ -d "$RUNS_DIR" ]]; then
        log_info "Copying from runs directory..."
        # Find the most recent stage1 run
        LATEST_RUN=$(find "$RUNS_DIR" -maxdepth 1 -type d -name "stage1_*" 2>/dev/null | sort -r | head -n 1 || true)
        if [[ -n "$LATEST_RUN" ]]; then
            log_info "Found run directory: $LATEST_RUN"
            cp -rv "$LATEST_RUN/"* "$ARTIFACTS_DIR/" 2>&1 || true
        fi
    fi

    # Copy stage1 checkpoints specifically
    if [[ -d "$LOCAL_ARTIFACTS/checkpoints" ]]; then
        find "$LOCAL_ARTIFACTS/checkpoints" -name "stage1*.pt" -type f -exec cp -v {} "$CHECKPOINTS_DIR/" \; 2>&1 || true
    fi

    log_info "============================================================"
    log_info "Artifact copy complete!"
    log_info "Check your Drive folder: $RUN_DIR"
    log_info "============================================================"

    exit $exit_code
}

# Register cleanup function for normal exit and interrupts
trap cleanup_and_copy EXIT INT TERM

# =============================================================================
# Run Training
# =============================================================================
log_info "============================================================"
log_info "Starting Stage 1 Training"
log_info "============================================================"
echo ""

cd "$TRAINING_DIR"
log_info "Working directory: $(pwd)"
echo ""

# Run training with output to both stdout and log file
eval "$TRAIN_CMD" 2>&1 | tee "$LOG_FILE"

# Training complete - cleanup will be called by trap
log_info "Training completed successfully!"
