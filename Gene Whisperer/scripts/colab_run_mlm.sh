#!/usr/bin/env bash
# chmod +x "Gene Whisperer/scripts/colab_run_mlm.sh"
#
# Colab MLM Pretraining Script for Gene Whisperer
# Runs MLM pretraining with robust checkpointing, resume support, and Drive backup.
#
# Usage (in Colab):
#   # Auto mode (uses config; multi-kmer sequential if enabled):
#   !bash "Gene Whisperer/scripts/colab_run_mlm.sh" --run_name my_experiment
#
#   # Force single-tokenizer (BPE) pretraining:
#   !bash "Gene Whisperer/scripts/colab_run_mlm.sh" --single_kmer --run_name my_single_run
#
#   # Force multi-kmer sequential pretraining (k=3,4,5,6):
#   !bash "Gene Whisperer/scripts/colab_run_mlm.sh" --multi_kmer --run_name my_multi_run
#
#   # Multi-kmer with specific k-mers:
#   !bash "Gene Whisperer/scripts/colab_run_mlm.sh" --multi_kmer --kmers 3 4 --run_name my_run
#
#   # Resume interrupted multi-kmer training:
#   !bash "Gene Whisperer/scripts/colab_run_mlm.sh" --multi_kmer --resume --run_name my_run
#
# This script:
#   1. Checks Google Drive is mounted
#   2. Runs colab_bootstrap.sh to set up environment
#   3. Creates output directories on Drive
#   4. Builds missing multi-kmer vocabularies (if needed)
#   5. Runs MLM pretraining (single or multi-kmer) with logging to Drive
#   5. Implements resume from existing checkpoints
#   6. Copies artifacts to Drive on completion or Ctrl+C

set -euo pipefail

# =============================================================================
# Default Configuration
# =============================================================================
DEFAULT_DRIVE_DIR="/content/drive/MyDrive/GeneWhispererRuns"
DEFAULT_CONFIG="Gene Whisperer/training/config.yaml"
DEFAULT_RUN_NAME="mlm_$(date +%Y%m%d_%H%M%S)"

# =============================================================================
# Parse Arguments
# =============================================================================
DRIVE_DIR="$DEFAULT_DRIVE_DIR"
CONFIG="$DEFAULT_CONFIG"
KMER=""
RUN_NAME="$DEFAULT_RUN_NAME"
PRETRAIN_MODE="auto"
KMERS=""
RESUME=false
FORCE=false
SKIP_VOCAB_BUILD=false

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
            if [[ "$PRETRAIN_MODE" == "auto" ]]; then
                PRETRAIN_MODE="single"
            fi
            shift 2
            ;;
        --run_name)
            RUN_NAME="$2"
            shift 2
            ;;
        --multi_kmer)
            PRETRAIN_MODE="multi"
            shift
            ;;
        --kmers)
            shift
            KMERS=""
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                KMERS="$KMERS $1"
                shift
            done
            KMERS="${KMERS# }"
            ;;
        --single_kmer)
            PRETRAIN_MODE="single"
            shift
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --skip_vocab_build)
            SKIP_VOCAB_BUILD=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --drive_dir DIR   Google Drive output directory (default: $DEFAULT_DRIVE_DIR)"
            echo "  --config PATH     Config file path relative to git root (default: $DEFAULT_CONFIG)"
            echo "  --kmer K          Legacy k-mer size flag (ignored for BPE single-tokenizer mode)"
            echo "  --run_name NAME   Run name for output directory (default: mlm_YYYYMMDD_HHMMSS)"
            echo "  --multi_kmer      Force multi-kmer sequential pretraining"
            echo "  --single_kmer     Force single-tokenizer (BPE) pretraining"
            echo "  --kmers K1 K2     Specific k-mers to train (e.g., 3 4). Multi-kmer mode only"
            echo "  --resume          Resume from existing checkpoints (single) or skip existing (multi)"
            echo "  --force           Force retrain even if checkpoints exist (for multi-kmer mode)"
            echo "  --skip_vocab_build  Skip building missing multi-kmer vocabularies"
            echo "  -h, --help        Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Auto mode (uses config; multi-kmer sequential if enabled):"
            echo "  $0 --run_name my_experiment"
            echo ""
            echo "  # Force multi-kmer sequential pretraining (k=3,4,5,6):"
            echo "  $0 --multi_kmer --run_name my_multi_run"
            echo ""
            echo "  # Multi-kmer with specific k-mers:"
            echo "  $0 --multi_kmer --kmers 3 4 --run_name my_run"
            echo ""
            echo "  # Force single-tokenizer (BPE) pretraining:"
            echo "  $0 --single_kmer --run_name my_single_run"
            echo ""
            echo "  # Resume interrupted multi-kmer training:"
            echo "  $0 --multi_kmer --resume --run_name my_run"
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

resolve_config_settings() {
    python - <<'PY'
import os
import shlex
from pathlib import Path

import yaml

config_path = Path(os.environ["CONFIG_PATH"])
cfg = yaml.safe_load(config_path.read_text()) or {}
multi_cfg = cfg.get("multi_kmer_pretrain", {}) or {}

enabled = cfg.get("multi_kmer_pretrain_enabled")
if enabled is None:
    enabled = multi_cfg.get("multi_kmer_pretrain_enabled")
if enabled is None:
    enabled = False

pretrain_mode = cfg.get("pretrain_mode") or multi_cfg.get("pretrain_mode") or ""
kmers = cfg.get("pretrain_kmers") or multi_cfg.get("pretrain_kmers") or [3, 4, 5, 6]
vocab_dir = cfg.get("vocab_cache_dir") or cfg.get("vocab_paths", {}).get("vocab_cache_dir") or "../artifacts/vocabs"

training_dir = Path(os.environ["TRAINING_DIR"])
vocab_dir = Path(vocab_dir)
if not vocab_dir.is_absolute():
    vocab_dir = (training_dir / vocab_dir).resolve()

print(f"CONFIG_MULTI_KMER_ENABLED={str(bool(enabled)).lower()}")
print(f"CONFIG_PRETRAIN_MODE={shlex.quote(str(pretrain_mode))}")
print(f"CONFIG_PRETRAIN_KMERS={shlex.quote(' '.join(str(k) for k in kmers))}")
print(f"CONFIG_VOCAB_DIR={shlex.quote(str(vocab_dir))}")
PY
}

ensure_multi_kmer_vocabs() {
    if [[ "$PRETRAIN_MODE" != "multi" ]]; then
        return
    fi
    if [[ "$SKIP_VOCAB_BUILD" == true ]]; then
        log_warn "Skipping vocabulary build (--skip_vocab_build)"
        return
    fi
    if [[ -z "$RESOLVED_KMERS" ]]; then
        log_warn "No k-mers resolved for vocabulary check"
        return
    fi

    local missing_kmers
    missing_kmers="$(CONFIG_VOCAB_DIR="$CONFIG_VOCAB_DIR" RESOLVED_KMERS="$RESOLVED_KMERS" python - <<'PY'
import os
from pathlib import Path

vocab_dir = Path(os.environ["CONFIG_VOCAB_DIR"])
kmers = os.environ.get("RESOLVED_KMERS", "").split()
missing = [k for k in kmers if not (vocab_dir / f"mlm_k{k}_vocab.json").exists()]
print(" ".join(missing))
PY
)"

    if [[ -n "$missing_kmers" ]]; then
        log_info "Missing vocabularies for k-mers: $missing_kmers"
        log_info "Building vocabularies in: $CONFIG_VOCAB_DIR"
        python "$TRAINING_DIR/build_multi_kmer_vocabs.py" --output-dir "$CONFIG_VOCAB_DIR" --kmers $missing_kmers
    else
        log_info "All vocabularies present in: $CONFIG_VOCAB_DIR"
    fi
}

# =============================================================================
# Check Google Drive Mount
# =============================================================================
echo "============================================================"
echo "Gene Whisperer - Colab MLM Pretraining"
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
# Resolve Training Paths and Config-Driven Settings
# =============================================================================
TRAINING_DIR="$REPO_ROOT/Gene Whisperer/training"
eval "$(CONFIG_PATH="$CONFIG_PATH" TRAINING_DIR="$TRAINING_DIR" resolve_config_settings)"

if [[ "$PRETRAIN_MODE" == "auto" ]]; then
    if [[ "$CONFIG_MULTI_KMER_ENABLED" == "true" ]]; then
        PRETRAIN_MODE="multi"
    else
        PRETRAIN_MODE="single"
    fi
fi

if [[ "$PRETRAIN_MODE" == "multi" && -n "$KMER" ]]; then
    log_error "--kmer is only valid in single-tokenizer mode (BPE); flag is legacy"
    exit 1
fi

if [[ "$PRETRAIN_MODE" == "single" && -n "$KMERS" ]]; then
    log_error "--kmers is only valid in multi-kmer mode"
    exit 1
fi

RESOLVED_KMERS="$KMERS"
if [[ -z "$RESOLVED_KMERS" ]]; then
    RESOLVED_KMERS="$CONFIG_PRETRAIN_KMERS"
fi

if [[ "$PRETRAIN_MODE" == "multi" ]]; then
    log_info "Mode: Multi-kmer sequential pretraining"
    log_info "K-mers: $RESOLVED_KMERS"
    if [[ -n "$CONFIG_PRETRAIN_MODE" && "$CONFIG_PRETRAIN_MODE" != "sequential" ]]; then
        log_warn "Config pretrain_mode=$CONFIG_PRETRAIN_MODE (script runs sequential mode)"
    fi
else
    log_info "Mode: Single-tokenizer pretraining (BPE)"
fi
echo ""

# =============================================================================
# Check for Existing Checkpoint (Resume Support)
# =============================================================================
RESUME_ARG=""
LATEST_CHECKPOINT=""

if [[ "$PRETRAIN_MODE" == "single" ]]; then
    if [[ -d "$ARTIFACTS_DIR/checkpoints" ]]; then
        LATEST_CHECKPOINT="$(ls -t "$ARTIFACTS_DIR"/checkpoints/checkpoint_step_*.pt 2>/dev/null | head -n 1 || true)"
    fi

    if [[ -n "$LATEST_CHECKPOINT" ]]; then
        log_info "Found training checkpoint: $LATEST_CHECKPOINT"
        if [[ "$RESUME" == true ]]; then
            RESUME_ARG="--resume_path \"$LATEST_CHECKPOINT\""
            log_info "Resume enabled for single-tokenizer pretraining"
        else
            log_info "Resume disabled; starting fresh"
        fi
    else
        log_info "No training checkpoint found. Starting fresh training."
    fi
else
    if [[ "$RESUME" == true ]]; then
        log_info "Resume mode: ON (multi-kmer skips existing checkpoints)"
    fi
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
LOG_FILE="$LOGS_DIR/mlm_training_$(date +%Y%m%d_%H%M%S).log"

# Build the training command based on mode
ensure_multi_kmer_vocabs

if [[ "$PRETRAIN_MODE" == "multi" ]]; then
    # Multi-kmer mode: use pretrain_multi_kmer.py
    log_info "Mode: Multi-kmer pretraining"
    TRAIN_CMD="python pretrain_multi_kmer.py --config \"$CONFIG_PATH\""

    # Add resolved k-mers
    if [[ -n "$RESOLVED_KMERS" ]]; then
        log_info "K-mers to train: $RESOLVED_KMERS"
        TRAIN_CMD="$TRAIN_CMD --kmers $RESOLVED_KMERS"
    fi

    # Add resume flag if specified
    if [[ "$RESUME" == true ]]; then
        log_info "Resume mode: ON (will continue from existing checkpoints)"
        TRAIN_CMD="$TRAIN_CMD --resume"
    fi

    # Add force flag if specified
    if [[ "$FORCE" == true ]]; then
        log_info "Force mode: ON (will retrain all k-mers)"
        TRAIN_CMD="$TRAIN_CMD --force"
    fi
else
    # Single-tokenizer mode (BPE): use pretrain_mlm.py
    log_info "Mode: Single-tokenizer pretraining (BPE)"
    TRAIN_CMD="python pretrain_mlm.py --config \"$CONFIG_PATH\""

    # Add kmer argument if specified
    # Note: pretrain_mlm.py may not support --kmer directly; it reads from config
    if [[ -n "$KMER" ]]; then
        log_info "K-mer size specified: $KMER"
        log_warn "Note: --kmer may not be supported by pretrain_mlm.py; check config instead."
        # TRAIN_CMD="$TRAIN_CMD --kmer $KMER"
    fi

    # Add resume argument if available
    if [[ -n "$RESUME_ARG" ]]; then
        TRAIN_CMD="$TRAIN_CMD $RESUME_ARG"
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

    # Copy any .json report files from training directory
    log_info "Copying any .json report files..."
    find "$TRAINING_DIR" -maxdepth 1 -name "*.json" -type f -exec cp -v {} "$ARTIFACTS_DIR/" \; 2>&1 || true

    # Copy mlm encoder checkpoints for all k-mers (3, 4, 5, 6)
    for k in 3 4 5 6; do
        if [[ -f "$LOCAL_ARTIFACTS/mlm_encoder_k${k}.pt" ]]; then
            log_info "Copying mlm_encoder_k${k}.pt..."
            cp -v "$LOCAL_ARTIFACTS/mlm_encoder_k${k}.pt" "$CHECKPOINTS_DIR/" 2>&1 || true
        fi
        # Also copy metadata files
        if [[ -f "$LOCAL_ARTIFACTS/mlm_encoder_k${k}.metadata.json" ]]; then
            cp -v "$LOCAL_ARTIFACTS/mlm_encoder_k${k}.metadata.json" "$ARTIFACTS_DIR/" 2>&1 || true
        fi
    done

    # Copy multi-kmer training summary if it exists
    if [[ -f "$LOCAL_ARTIFACTS/multi_kmer_pretrain_summary.json" ]]; then
        log_info "Copying multi-kmer training summary..."
        cp -v "$LOCAL_ARTIFACTS/multi_kmer_pretrain_summary.json" "$ARTIFACTS_DIR/" 2>&1 || true
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
log_info "Starting MLM Pretraining"
log_info "============================================================"
echo ""

cd "$TRAINING_DIR"
log_info "Working directory: $(pwd)"
echo ""

# Run training with output to both stdout and log file
eval "$TRAIN_CMD" 2>&1 | tee "$LOG_FILE"

# Training complete - cleanup will be called by trap
log_info "Training completed successfully!"
