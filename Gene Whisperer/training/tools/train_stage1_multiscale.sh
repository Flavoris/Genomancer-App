#!/usr/bin/env bash
#
# train_stage1_multiscale.sh - Train both k=4 and k=6 models, then run ensemble evaluation
#
# This script:
#   1. Trains k=4 small model -> saves checkpoint to stage1_k4.pt
#   2. Trains k=6 small model -> saves checkpoint to stage1_k6.pt
#   3. Runs eval_stage1_ensemble.py to evaluate:
#      - k=4 alone
#      - k=6 alone
#      - Ensemble average of probabilities (0.5*(p4+p6))
#
# Usage:
#   ./tools/train_stage1_multiscale.sh [--config CONFIG_PATH]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRAINING_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_PATH="${TRAINING_DIR}/config.yaml"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        *)
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
echo "MULTI-SCALE STAGE 1 TRAINING PIPELINE"
echo "============================================================"
echo "Config:       $CONFIG_PATH"
echo "Training dir: $TRAINING_DIR"
echo ""
echo "This script will:"
echo "  1. Train k=4 small model"
echo "  2. Train k=6 small model"
echo "  3. Evaluate ensemble (k4, k6, and combined)"
echo "============================================================"
echo ""

cd "$TRAINING_DIR"

# Step 1: Train k=4 model
echo "============================================================"
echo "STEP 1/3: Training k=4 model"
echo "============================================================"
./tools/train_stage1_k4_small.sh --config "$CONFIG_PATH"
echo ""

# Step 2: Train k=6 model
echo "============================================================"
echo "STEP 2/3: Training k=6 model"
echo "============================================================"
./tools/train_stage1_k6_small.sh --config "$CONFIG_PATH"
echo ""

# Step 3: Run ensemble evaluation
echo "============================================================"
echo "STEP 3/3: Running ensemble evaluation"
echo "============================================================"
python eval_stage1_ensemble.py --config "$CONFIG_PATH"
echo ""

echo "============================================================"
echo "MULTI-SCALE TRAINING COMPLETE"
echo "============================================================"
echo "Checkpoints saved to:"
echo "  - artifacts/checkpoints/stage1_k4.pt"
echo "  - artifacts/checkpoints/stage1_k6.pt"
echo ""
echo "Ensemble report saved to runs/<timestamp>/stage1_ensemble_report.json"
echo "============================================================"
