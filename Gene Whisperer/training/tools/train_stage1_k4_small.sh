#!/usr/bin/env bash
#
# train_stage1_k4_small.sh - Lightweight Stage 1 training for k=4
#
# Uses a smaller model configuration optimized for MacBook:
#   - embedding_dim=192
#   - transformer_layers=0 (rely on TCN/multiscale + post adapter)
#   - post_cnn_transformer_layers=1
#   - tcn_hidden=192
#   - multiscale_channels=48
#
# Usage:
#   ./tools/train_stage1_k4_small.sh [--overfit_debug] [--config CONFIG_PATH]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRAINING_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_PATH="${TRAINING_DIR}/config.yaml"

# Collect extra args to pass through
EXTRA_ARGS=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --overfit_debug)
            EXTRA_ARGS+=("--overfit_debug")
            shift
            ;;
        *)
            # Pass through any other arguments
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

# Verify config exists
if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "Error: Config file not found: $CONFIG_PATH"
    exit 1
fi

echo "============================================================"
echo "STAGE 1 TRAINING - k=4 SMALL (Lightweight)"
echo "============================================================"
echo "Config:       $CONFIG_PATH"
echo "Training dir: $TRAINING_DIR"
echo ""
echo "Architecture overrides:"
echo "  embedding_dim=192"
echo "  transformer_layers=0"
echo "  post_cnn_transformer_layers=1"
echo "  tcn_hidden=192"
echo "  multiscale_channels=48"
echo "  checkpoint_name=stage1_k4.pt"
echo "============================================================"

cd "$TRAINING_DIR"

python train_stage1.py \
    --config "$CONFIG_PATH" \
    --kmer 4 \
    --embedding_dim 192 \
    --transformer_layers 0 \
    --post_cnn_transformer_layers 1 \
    --tcn_hidden 192 \
    --multiscale_channels 48 \
    --checkpoint_name stage1_k4.pt \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"

echo ""
echo "============================================================"
echo "STAGE 1 k=4 SMALL TRAINING COMPLETE"
echo "============================================================"
