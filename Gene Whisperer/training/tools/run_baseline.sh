#!/usr/bin/env bash
#
# run_baseline.sh - Run full Gene Whisperer baseline training pipeline
#
# This script:
# 1. Runs train_stage1.py with config.yaml (uses seed from config)
# 2. Runs train_stage2.py with config.yaml
# 3. Copies final JSON reports to runs/baseline/
#
# Usage:
#   ./tools/run_baseline.sh [--config CONFIG_PATH]
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
            echo "Unknown option: $1"
            echo "Usage: $0 [--config CONFIG_PATH]"
            exit 1
            ;;
    esac
done

# Verify config exists
if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "Error: Config file not found: $CONFIG_PATH"
    exit 1
fi

echo "============================================================"
echo "GENE WHISPERER BASELINE TRAINING"
echo "============================================================"
echo "Config:       $CONFIG_PATH"
echo "Training dir: $TRAINING_DIR"
echo "============================================================"

# Print system info first
echo ""
echo "[1/4] Printing system info..."
python "${SCRIPT_DIR}/print_system_info.py" --config "$CONFIG_PATH"

# Create baseline output directory
BASELINE_DIR="${TRAINING_DIR}/../runs/baseline"
mkdir -p "$BASELINE_DIR"
echo ""
echo "Baseline output directory: $BASELINE_DIR"

# Run Stage 1 training
echo ""
echo "[2/4] Running Stage 1 training..."
echo "============================================================"
cd "$TRAINING_DIR"
python train_stage1.py --config "$CONFIG_PATH"
STAGE1_EXIT=$?

if [[ $STAGE1_EXIT -ne 0 ]]; then
    echo "Error: Stage 1 training failed with exit code $STAGE1_EXIT"
    exit $STAGE1_EXIT
fi

# Run Stage 2 training
echo ""
echo "[3/4] Running Stage 2 training..."
echo "============================================================"
python train_stage2.py --config "$CONFIG_PATH"
STAGE2_EXIT=$?

if [[ $STAGE2_EXIT -ne 0 ]]; then
    echo "Error: Stage 2 training failed with exit code $STAGE2_EXIT"
    exit $STAGE2_EXIT
fi

# Copy reports to baseline directory
echo ""
echo "[4/4] Copying reports to baseline directory..."
echo "============================================================"

ARTIFACTS_DIR="${TRAINING_DIR}/../artifacts"

# Find and copy Stage 1 report (format: stage1_report_k{kmer}.json)
for report in "${ARTIFACTS_DIR}"/stage1_report_k*.json; do
    if [[ -f "$report" ]]; then
        cp "$report" "$BASELINE_DIR/"
        echo "Copied: $(basename "$report")"
    fi
done

# Find and copy Stage 2 report (format: stage2_report_k{kmer}.json)
for report in "${ARTIFACTS_DIR}"/stage2_report_k*.json; do
    if [[ -f "$report" ]]; then
        cp "$report" "$BASELINE_DIR/"
        echo "Copied: $(basename "$report")"
    fi
done

# Also copy the config used for this baseline run
cp "$CONFIG_PATH" "$BASELINE_DIR/config_baseline.yaml"
echo "Copied: config_baseline.yaml"

# Generate a summary timestamp file
TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
cat > "$BASELINE_DIR/baseline_info.txt" << EOF
Gene Whisperer Baseline Run
============================
Timestamp: $TIMESTAMP
Config:    $CONFIG_PATH

Reports:
$(ls -la "$BASELINE_DIR"/*.json 2>/dev/null || echo "  (no JSON reports found)")
EOF

echo ""
echo "============================================================"
echo "BASELINE TRAINING COMPLETE"
echo "============================================================"
echo "Reports saved to: $BASELINE_DIR"
echo ""
ls -la "$BASELINE_DIR"
