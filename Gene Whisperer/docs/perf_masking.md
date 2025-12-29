# Masking Performance Benchmark - Baseline

## Test Information

- **Date/Time**: 2025-12-28
- **Machine**: Apple M4 (arm64), macOS Darwin 25.1.0
- **Python**: Local (CPU)

## Command Lines Used

```bash
cd "Gene Whisperer"
python scripts/bench_masking.py --seq_len=81 --iters=200
python scripts/bench_masking.py --seq_len=234 --iters=200
```

## Configuration

- Batch size: 64
- Iterations: 200
- K-mer size: 6
- Vocabulary size: 4099 tokens
- Warmup iterations: 10
- Device: CPU

## Results

### seq_len=81

| Metric | Value |
|--------|-------|
| Avg ms/iter | 181.246 ms |
| Tokens/sec | 28,602 |
| Min ms/iter | 170.314 ms |
| Max ms/iter | 198.881 ms |
| Total time | 36.249 s |

### seq_len=234

| Metric | Value |
|--------|-------|
| Avg ms/iter | 1489.612 ms |
| Tokens/sec | 10,054 |
| Min ms/iter | 1418.059 ms |
| Max ms/iter | 1594.864 ms |
| Total time | 297.922 s |

## Ratio Analysis

| Metric | seq_len=81 | seq_len=234 | Ratio (234/81) |
|--------|------------|-------------|----------------|
| Sequence length | 81 | 234 | 2.89x |
| Avg ms/iter | 181.246 | 1489.612 | 8.22x |
| Tokens/sec | 28,602 | 10,054 | 0.35x (2.84x slower) |

**Observation**: Time scales ~8.2x for a ~2.9x increase in sequence length, indicating worse than linear scaling (approximately O(n^2.8)).
