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

---

## After Prompt 2: Rejection Sampling Optimization

### Change Summary

Replaced the O(L * placements * scan) span placement with fast rejection sampling:
- Instead of scanning every possible start position to build `valid_starts`, now randomly samples up to 16 start positions per span length
- Added fallback fill: if rejection sampling falls short, randomly selects remaining non-masked non-pad positions
- Preserved all existing rules: no touching spans, PAD safety, span length fallback

### Results

#### seq_len=81

| Metric | Before | After | Speedup |
|--------|--------|-------|---------|
| Avg ms/iter | 181.246 ms | 2.622 ms | **69x** |
| Tokens/sec | 28,602 | 1,976,743 | **69x** |
| Min ms/iter | 170.314 ms | 2.435 ms | 70x |
| Max ms/iter | 198.881 ms | 2.955 ms | 67x |
| Total time | 36.249 s | 0.524 s | 69x |

#### seq_len=234

| Metric | Before | After | Speedup |
|--------|--------|-------|---------|
| Avg ms/iter | 1489.612 ms | 7.088 ms | **210x** |
| Tokens/sec | 10,054 | 2,112,876 | **210x** |
| Min ms/iter | 1418.059 ms | 6.693 ms | 212x |
| Max ms/iter | 1594.864 ms | 7.488 ms | 213x |
| Total time | 297.922 s | 1.418 s | 210x |

### Scaling Analysis

| Metric | seq_len=81 | seq_len=234 | Ratio (234/81) |
|--------|------------|-------------|----------------|
| Sequence length | 81 | 234 | 2.89x |
| Avg ms/iter | 2.622 ms | 7.088 ms | 2.70x |
| Tokens/sec | 1,976,743 | 2,112,876 | 1.07x (faster!) |

**Observation**: Time now scales ~2.7x for a ~2.9x increase in sequence length, indicating approximately linear scaling O(n). The rejection sampling approach eliminates the quadratic scan overhead.

---

## After Prompt 3: Vectorized Special-Token Exclusion & Random-Token Replacement

### Change Summary

Vectorized two remaining Python loops in `mask_tokens_span()`:

**A) Special-token label exclusion** - Replaced per-token loop with single vectorized mask:
```python
# Before:
special_token_ids = {mask_id, unk_id, pad_id}
for special_id in special_token_ids:
    labels[inputs == special_id] = -100

# After:
special_mask = (inputs == mask_id) | (inputs == unk_id) | (inputs == pad_id)
labels[special_mask] = -100
```

**B) Random token replacement** - Replaced Python list comprehension with torch sampling:
```python
# Before:
random_ids = torch.tensor([random.choice(vocab._base_token_ids) for _ in range(num_rand)], ...)

# After:
base_ids = torch.as_tensor(vocab._base_token_ids, device=device, dtype=inputs.dtype)
rand_idx = torch.randint(0, base_ids.numel(), (num_rand,), device=device)
random_ids = base_ids[rand_idx]
```

### Results

#### seq_len=81

| Metric | Before (Prompt 2) | After (Prompt 3) | Change |
|--------|-------------------|------------------|--------|
| Avg ms/iter | 2.622 ms | 2.702 ms | +3.1% |
| Tokens/sec | 1,976,743 | 1,918,738 | -2.9% |
| Min ms/iter | 2.435 ms | 2.471 ms | +1.5% |
| Max ms/iter | 2.955 ms | 3.065 ms | +3.7% |
| Total time | 0.524 s | 0.540 s | +3.1% |

#### seq_len=234

| Metric | Before (Prompt 2) | After (Prompt 3) | Change |
|--------|-------------------|------------------|--------|
| Avg ms/iter | 7.088 ms | 7.218 ms | +1.8% |
| Tokens/sec | 2,112,876 | 2,074,745 | -1.8% |
| Min ms/iter | 6.693 ms | 6.813 ms | +1.8% |
| Max ms/iter | 7.488 ms | 15.215 ms | +103% (outlier) |
| Total time | 1.418 s | 1.444 s | +1.8% |

### Analysis

The vectorization changes show negligible performance difference (within ~3% noise margin). This is expected because:

1. **Special-token loop was already O(3)** - Only 3 iterations, so vectorizing provides minimal benefit
2. **Random token replacement is ~10% of masked tokens** - With span masking at 15% and 10% random replacement, only ~1.5% of tokens trigger this path
3. **Torch tensor creation overhead** - `torch.as_tensor()` has fixed overhead that offsets gains for small arrays
4. **Already at ~2M tokens/sec** - The previous rejection sampling optimization achieved the major speedup; these are micro-optimizations

**Conclusion**: The code is now cleaner and more idiomatic PyTorch, but performance is equivalent. The 69-210x speedup from Prompt 2 remains the dominant optimization.
