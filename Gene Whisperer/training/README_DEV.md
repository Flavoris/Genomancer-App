# Gene Whisperer Development Guide

## How to Reproduce Baseline

This section documents how to run a reproducible baseline training run.

### Prerequisites

1. Ensure you have the required Python packages installed:
   ```bash
   pip install torch pyyaml tqdm numpy
   ```

2. Ensure the data files are in place:
   - `../data/stage1_train.tsv`
   - `../data/stage1_val.tsv`
   - `../data/stage2_train.tsv`
   - `../data/stage2_val.tsv`

### Quick Start: Full Baseline Run

Run the complete baseline training pipeline (Stage 1 + Stage 2):

```bash
cd Gene\ Whisperer/training
./tools/run_baseline.sh
```

This script will:
1. Print system info (Python version, PyTorch version, device, model params)
2. Run Stage 1 training with `config.yaml`
3. Run Stage 2 training with `config.yaml`
4. Copy final reports to `runs/baseline/`

### Manual Steps

#### 1. Print System Info

```bash
cd Gene\ Whisperer/training
python tools/print_system_info.py --config config.yaml
```

This prints:
- Python version
- PyTorch version
- Device selected (cuda/mps/cpu)
- Stage 1 model parameter count

#### 2. Run Stage 1 Training

```bash
python train_stage1.py --config config.yaml
```

The seed is read from `config.yaml` (default: 1337).

#### 3. Run Stage 2 Training

```bash
python train_stage2.py --config config.yaml
```

### Overfit Debug Mode

To verify the model can learn (sanity check):

```bash
python train_stage1.py --config config.yaml --overfit_debug
```

This trains on 32 samples for 200 steps. Expected behavior:
- Loss should decrease
- Accuracy should reach near 100%

If overfit fails, the script will diagnose potential issues (frozen weights, zero gradients, etc.)

### Configuration

The main configuration file is `config.yaml`. Key settings:

- `seed`: Random seed for reproducibility (default: 1337)
- `epochs`: Number of training epochs
- `batch_size`: Batch size for training
- `lr`: Learning rate

### Output Files

Training runs create:
- `runs/stage1_<timestamp>/resolved_config.json` - Full resolved config with system info
- `runs/stage2_<timestamp>/resolved_config.json` - Full resolved config with system info
- `../artifacts/stage1_report_k*.json` - Stage 1 training report
- `../artifacts/stage2_report_k*.json` - Stage 2 training report
- `../artifacts/checkpoints/stage1_k*.pt` - Stage 1 model checkpoint
- `../artifacts/checkpoints/stage2_k*.pt` - Stage 2 model checkpoint

The baseline script copies reports to `runs/baseline/` for easy comparison.
