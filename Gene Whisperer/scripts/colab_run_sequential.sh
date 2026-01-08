#!/usr/bin/env bash
# chmod +x "Gene Whisperer/scripts/colab_run_sequential.sh"
#
# Colab Sequential Training Script for Gene Whisperer
# Implements msBERT's two-stage approach for 90%+ accuracy:
#   1. Train Stage 1 models for each k-mer (3, 4, 5, 6)
#   2. Evaluate soft voting ensemble
#
# Usage (in Colab):
#   !bash "Gene Whisperer/scripts/colab_run_sequential.sh"
#   !bash "Gene Whisperer/scripts/colab_run_sequential.sh" --run_name my_experiment
#   !bash "Gene Whisperer/scripts/colab_run_sequential.sh" --kmers 4 6
#   !bash "Gene Whisperer/scripts/colab_run_sequential.sh" --skip_stage2
#
# This script:
#   1. Checks Google Drive is mounted
#   2. Runs colab_bootstrap.sh to set up environment
#   3. Creates output directories on Drive
#   4. Runs sequential training (all k-mers + ensemble)
#   5. Copies all artifacts and checkpoints to Drive

set -euo pipefail

# =============================================================================
# Default Configuration
# =============================================================================
DEFAULT_DRIVE_DIR="/content/drive/MyDrive/GeneWhispererRuns"
DEFAULT_CONFIG="Gene Whisperer/training/config.yaml"
DEFAULT_RUN_NAME="sequential_$(date +%Y%m%d_%H%M%S)"

# =============================================================================
# Parse Arguments
# =============================================================================
DRIVE_DIR="$DEFAULT_DRIVE_DIR"
CONFIG="$DEFAULT_CONFIG"
KMERS=""
RUN_NAME="$DEFAULT_RUN_NAME"
SKIP_STAGE2=""

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
        --kmers)
            # Collect all k-mer arguments
            shift
            KMERS=""
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                KMERS="$KMERS $1"
                shift
            done
            ;;
        --run_name)
            RUN_NAME="$2"
            shift 2
            ;;
        --skip_stage2)
            SKIP_STAGE2="--skip_stage2"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Sequential Fine-Tuning for Gene Whisperer (msBERT approach)"
            echo ""
            echo "Options:"
            echo "  --drive_dir DIR     Google Drive output directory (default: $DEFAULT_DRIVE_DIR)"
            echo "  --config PATH       Config file path relative to git root (default: $DEFAULT_CONFIG)"
            echo "  --kmers K1 K2 ...   K-mer sizes to train (default: 3 4 5 6 from config)"
            echo "  --run_name NAME     Run name for output directory (default: sequential_YYYYMMDD_HHMMSS)"
            echo "  --skip_stage2       Skip Stage 2 training (only train Stage 1 models)"
            echo "  -h, --help          Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Train all k-mers (3, 4, 5, 6) with ensemble evaluation"
            echo "  !bash \"Gene Whisperer/scripts/colab_run_sequential.sh\""
            echo ""
            echo "  # Train only k=4 and k=6"
            echo "  !bash \"Gene Whisperer/scripts/colab_run_sequential.sh\" --kmers 4 6"
            echo ""
            echo "  # Custom run name"
            echo "  !bash \"Gene Whisperer/scripts/colab_run_sequential.sh\" --run_name my_experiment"
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
echo "Gene Whisperer - Sequential Fine-Tuning (msBERT approach)"
echo "============================================================"
echo ""
echo "This script will train Stage 1 models for multiple k-mer sizes"
echo "and evaluate them as a soft voting ensemble for 90%+ accuracy."
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
RESULTS_DIR="$RUN_DIR/results"

mkdir -p "$ARTIFACTS_DIR"
mkdir -p "$CHECKPOINTS_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "$RESULTS_DIR"

log_info "Run directory: $RUN_DIR"
log_info "  Artifacts:   $ARTIFACTS_DIR"
log_info "  Checkpoints: $CHECKPOINTS_DIR"
log_info "  Logs:        $LOGS_DIR"
log_info "  Results:     $RESULTS_DIR"
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
LOG_FILE="$LOGS_DIR/sequential_training_$(date +%Y%m%d_%H%M%S).log"
RESULTS_FILE="$RESULTS_DIR/sequential_results.json"

# Build the training command
TRAIN_CMD="python train_sequential.py --config \"$CONFIG_PATH\" --output \"$RESULTS_FILE\""

# Add kmers argument if specified
if [[ -n "$KMERS" ]]; then
    log_info "K-mer sizes specified:$KMERS"
    TRAIN_CMD="$TRAIN_CMD --kmers$KMERS"
else
    log_info "Using default k-mer sizes from config (typically 3, 4, 5, 6)"
fi

# Add skip_stage2 flag if specified
if [[ -n "$SKIP_STAGE2" ]]; then
    log_info "Stage 2 training will be skipped"
    TRAIN_CMD="$TRAIN_CMD $SKIP_STAGE2"
fi

log_info "Training command: $TRAIN_CMD"
log_info "Log file: $LOG_FILE"
log_info "Results file: $RESULTS_FILE"
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

    # Copy from runs directory if it exists
    RUNS_DIR="$REPO_ROOT/Gene Whisperer/runs"
    if [[ -d "$RUNS_DIR" ]]; then
        log_info "Copying from runs directory..."
        # Find all stage1 runs from this session
        find "$RUNS_DIR" -maxdepth 1 -type d -name "stage1_*" 2>/dev/null | while read -r run; do
            log_info "Copying run: $run"
            RUN_BASENAME=$(basename "$run")
            mkdir -p "$ARTIFACTS_DIR/$RUN_BASENAME"
            cp -rv "$run/"* "$ARTIFACTS_DIR/$RUN_BASENAME/" 2>&1 || true
        done
    fi

    # Copy all stage1 checkpoints (for each k-mer)
    if [[ -d "$LOCAL_ARTIFACTS/checkpoints" ]]; then
        log_info "Copying all stage1 checkpoints..."
        find "$LOCAL_ARTIFACTS/checkpoints" -name "stage1*.pt" -type f -exec cp -v {} "$CHECKPOINTS_DIR/" \; 2>&1 || true
    fi

    # Copy vocabs for reproducibility
    if [[ -d "$LOCAL_ARTIFACTS/vocabs" ]]; then
        log_info "Copying vocabularies..."
        mkdir -p "$ARTIFACTS_DIR/vocabs"
        cp -rv "$LOCAL_ARTIFACTS/vocabs/"* "$ARTIFACTS_DIR/vocabs/" 2>&1 || true
    fi

    log_info "============================================================"
    log_info "Artifact copy complete!"
    log_info "Check your Drive folder: $RUN_DIR"
    log_info ""
    log_info "Contents:"
    log_info "  $CHECKPOINTS_DIR/ - All k-mer checkpoints"
    log_info "  $RESULTS_DIR/     - Ensemble results JSON"
    log_info "  $LOGS_DIR/        - Training logs"
    log_info "  $ARTIFACTS_DIR/   - Vocabs and run artifacts"
    log_info "============================================================"

    exit $exit_code
}

# Register cleanup function for normal exit and interrupts
trap cleanup_and_copy EXIT INT TERM

# =============================================================================
# Run Training
# =============================================================================
log_info "============================================================"
log_info "Starting Sequential Fine-Tuning"
log_info "============================================================"
echo ""
log_info "This will train Stage 1 models for each k-mer size,"
log_info "then evaluate them as a soft voting ensemble."
echo ""

cd "$TRAINING_DIR"
log_info "Working directory: $(pwd)"
echo ""

# Run training with output to both stdout and log file
eval "$TRAIN_CMD" 2>&1 | tee "$LOG_FILE"

# Training complete - cleanup will be called by trap
log_info "============================================================"
log_info "Sequential training completed successfully!"
log_info "============================================================"
log_info ""
log_info "Results saved to: $RESULTS_FILE"
log_info "View ensemble metrics in the results JSON file."
