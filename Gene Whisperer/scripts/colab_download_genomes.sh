#!/usr/bin/env bash
# chmod +x "Gene Whisperer/scripts/colab_download_genomes.sh"
#
# Colab Genome Download Script for Gene Whisperer
# Downloads large genome files from Google Drive to the data directory.
#
# Usage (in Colab):
#   !bash "Gene Whisperer/scripts/colab_download_genomes.sh"
#
# Prerequisites:
#   1. Mount Google Drive in Colab:
#        from google.colab import drive
#        drive.mount('/content/drive')
#   2. Upload genome files to: /content/drive/MyDrive/GeneWhispererData/
#      Required files:
#        - human_genome (3.1 GB)
#        - japanese_rice_genome (373 MB)
#
# This script is idempotent - safe to run multiple times.

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================
# Default source directory on Google Drive
DEFAULT_DRIVE_DATA_DIR="/content/drive/MyDrive/GeneWhispererData"

# Parse optional --drive_data_dir argument
DRIVE_DATA_DIR="$DEFAULT_DRIVE_DATA_DIR"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --drive_data_dir)
            DRIVE_DATA_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --drive_data_dir DIR   Google Drive data directory (default: $DEFAULT_DRIVE_DATA_DIR)"
            echo "  -h, --help             Show this help message"
            exit 0
            ;;
        *)
            shift
            ;;
    esac
done

# Large genome files that need to be downloaded from Drive
declare -a LARGE_GENOMES=(
    "human_genome"
    "japanese_rice_genome"
)

# =============================================================================
# Helper Functions
# =============================================================================
log_info() {
    echo "[GENOME] $(date '+%H:%M:%S') - $*"
}

log_warn() {
    echo "[GENOME] $(date '+%H:%M:%S') - WARNING: $*"
}

log_error() {
    echo "[GENOME] $(date '+%H:%M:%S') - ERROR: $*" >&2
}

# =============================================================================
# Colab Detection
# =============================================================================
is_colab() {
    if [[ -n "${COLAB_RELEASE_TAG:-}" ]] || [[ -d "/content" ]]; then
        return 0
    fi
    return 1
}

if ! is_colab; then
    log_info "Not running in Colab - skipping genome download."
    log_info "Genome files should already be present locally."
    exit 0
fi

# =============================================================================
# Check Google Drive Mount
# =============================================================================
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

log_info "Google Drive is mounted."

# =============================================================================
# Get Repository Root and Data Directory
# =============================================================================
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"

if [[ -z "$REPO_ROOT" ]]; then
    log_error "Not inside a git repository."
    exit 1
fi

DATA_DIR="$REPO_ROOT/Gene Whisperer/data"

if [[ ! -d "$DATA_DIR" ]]; then
    log_info "Creating data directory: $DATA_DIR"
    mkdir -p "$DATA_DIR"
fi

# =============================================================================
# Check Drive Data Directory
# =============================================================================
if [[ ! -d "$DRIVE_DATA_DIR" ]]; then
    log_warn "Drive data directory not found: $DRIVE_DATA_DIR"
    echo ""
    echo "Please create this directory on Google Drive and upload genome files:"
    echo "  1. Go to Google Drive"
    echo "  2. Create folder: GeneWhispererData"
    echo "  3. Upload these files to the folder:"
    for genome in "${LARGE_GENOMES[@]}"; do
        echo "     - $genome"
    done
    echo ""
    echo "Or specify a custom directory with --drive_data_dir"
    exit 1
fi

log_info "Drive data directory: $DRIVE_DATA_DIR"

# =============================================================================
# Download Genome Files
# =============================================================================
echo ""
log_info "============================================================"
log_info "Downloading Large Genome Files from Google Drive"
log_info "============================================================"
echo ""

DOWNLOAD_COUNT=0
SKIP_COUNT=0
MISSING_COUNT=0

for genome in "${LARGE_GENOMES[@]}"; do
    SOURCE_FILE="$DRIVE_DATA_DIR/$genome"
    DEST_FILE="$DATA_DIR/$genome"

    # Check if already exists locally
    if [[ -f "$DEST_FILE" ]]; then
        LOCAL_SIZE=$(stat -f%z "$DEST_FILE" 2>/dev/null || stat -c%s "$DEST_FILE" 2>/dev/null || echo "0")
        log_info "[$genome] Already exists locally ($(numfmt --to=iec $LOCAL_SIZE 2>/dev/null || echo "${LOCAL_SIZE} bytes"))"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # Check if exists on Drive
    if [[ ! -f "$SOURCE_FILE" ]]; then
        log_warn "[$genome] Not found on Drive: $SOURCE_FILE"
        MISSING_COUNT=$((MISSING_COUNT + 1))
        continue
    fi

    # Get file size for progress info
    SOURCE_SIZE=$(stat -f%z "$SOURCE_FILE" 2>/dev/null || stat -c%s "$SOURCE_FILE" 2>/dev/null || echo "unknown")
    log_info "[$genome] Copying from Drive ($(numfmt --to=iec $SOURCE_SIZE 2>/dev/null || echo "$SOURCE_SIZE bytes"))..."

    # Copy with progress
    if command -v pv &> /dev/null; then
        pv "$SOURCE_FILE" > "$DEST_FILE"
    else
        cp "$SOURCE_FILE" "$DEST_FILE"
    fi

    log_info "[$genome] Download complete!"
    DOWNLOAD_COUNT=$((DOWNLOAD_COUNT + 1))
done

# =============================================================================
# Summary
# =============================================================================
echo ""
log_info "============================================================"
log_info "Genome Download Summary"
log_info "============================================================"
log_info "  Downloaded: $DOWNLOAD_COUNT"
log_info "  Skipped (already present): $SKIP_COUNT"
log_info "  Missing on Drive: $MISSING_COUNT"

if [[ $MISSING_COUNT -gt 0 ]]; then
    echo ""
    log_warn "Some genome files are missing from Google Drive."
    log_warn "Please upload them to: $DRIVE_DATA_DIR"
    log_warn "MLM pretraining may fail or skip these genomes."
fi

echo ""
log_info "Data directory contents:"
ls -lh "$DATA_DIR"/ | head -20

log_info "============================================================"
