# Gene Whisperer Training

## Model Architecture

Gene Whisperer v5.0 uses a simplified transformer-based architecture:

```
Input (DNA sequence) → k-mer Tokenization →
Embedding (384-dim) → 12-Layer Transformer →
Attention Pooling → [Concat Engineered Features] →
Classifier → Promoter Prediction
```

### Why Simplified?

Our ablation studies showed that CNN/TCN layers actually hurt performance:

- **With CNN/TCN**: 83.3% accuracy
- **Without CNN/TCN**: 84.1% accuracy (+0.8%)

The transformer alone, combined with engineered features, outperforms the more complex architecture while using 38% fewer parameters.

## Performance

| Metric | Value |
|--------|-------|
| Accuracy | 84.1% |
| MCC | 0.683 |
| AUC | 0.899 |
| Parameters | 23.9M |

## Setup

Setup virtual environment and install dependencies:

```bash
python -m venv ../.venv
source ../.venv/bin/activate
pip install torch torchvision torchaudio coremltools pandas numpy scikit-learn pyyaml tqdm
```

## Quick Start

```bash
# Train Stage 1 (promoter detection)
python train_stage1.py --config config.yaml

# Train Stage 2 (strong vs weak)
python train_stage2.py --config config.yaml

# Export to CoreML
python export_coreml.py --config config.yaml
```

## Migrating from v4.0

If you have checkpoints from the old CNN/TCN architecture:

1. The transformer weights will transfer automatically
2. CNN/TCN weights will be ignored (they hurt performance anyway)
3. You may want to fine-tune for 10-20 epochs to adapt

```python
from model import GeneWhispererStage1

model = GeneWhispererStage1(...)
model.load_legacy_checkpoint("path/to/v4_checkpoint.pt")
```
