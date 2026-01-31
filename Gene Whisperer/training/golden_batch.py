"""Golden batch utilities for deterministic reproducibility testing.

A "golden batch" is a saved snapshot of:
- Input tensors (masked tokens)
- Label tensors  
- Expected loss value from a single forward/backward pass

This enables verification that:
1. Same seed produces identical batch data
2. Same seed + same model = same loss (within tolerance)
3. Different seeds produce different batches

Usage:
    # Generate golden batch
    python golden_batch.py --config config.yaml --seed 42 --generate
    
    # Test against golden batch  
    python golden_batch.py --config config.yaml --seed 42 --test
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import torch
import torch.nn as nn
import yaml

LOGGER = logging.getLogger("gene_whisperer.golden_batch")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Default path for golden batch file
DEFAULT_GOLDEN_BATCH_PATH = Path(__file__).parent / "../artifacts/golden_batch.pt"


def save_golden_batch(
    inputs: torch.Tensor,
    labels: torch.Tensor,
    loss: float,
    seed: int,
    path: Path,
    metadata: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Save a golden batch checkpoint.
    
    Args:
        inputs: Input tensor (masked tokens) [B, L]
        labels: Label tensor [B, L] 
        loss: Expected loss value from forward pass
        seed: Seed used to generate this batch
        path: Path to save the checkpoint
        metadata: Optional additional metadata (model config, etc.)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    
    checkpoint = {
        "inputs": inputs.cpu().clone(),
        "labels": labels.cpu().clone(),
        "loss": float(loss),
        "seed": int(seed),
        "metadata": metadata or {},
    }
    
    torch.save(checkpoint, path)
    LOGGER.info("Saved golden batch to %s (seed=%d, loss=%.6f)", path, seed, loss)


def load_golden_batch(path: Path) -> Dict[str, Any]:
    """
    Load a golden batch checkpoint.
    
    Args:
        path: Path to the golden batch file
        
    Returns:
        Dict with keys: inputs, labels, loss, seed, metadata
    """
    if not path.exists():
        raise FileNotFoundError(f"Golden batch not found at {path}")
    
    checkpoint = torch.load(path, map_location="cpu", weights_only=True)
    LOGGER.info(
        "Loaded golden batch from %s (seed=%d, loss=%.6f)",
        path, checkpoint["seed"], checkpoint["loss"]
    )
    return checkpoint


def generate_golden_batch(
    dataloader: torch.utils.data.DataLoader,
    model: nn.Module,
    device: torch.device,
    seed: int,
    path: Path,
    metadata: Optional[Dict[str, Any]] = None,
) -> Tuple[torch.Tensor, torch.Tensor, float]:
    """
    Generate and save a golden batch from the first batch of a dataloader.
    
    Performs one forward/backward pass and saves the batch + loss.
    
    Args:
        dataloader: DataLoader to get batch from
        model: Model to run forward/backward pass
        device: Device to run on
        seed: Seed used for reproducibility
        path: Path to save golden batch
        metadata: Optional metadata to include
        
    Returns:
        Tuple of (inputs, labels, loss)
    """
    from seed_utils import set_global_seed
    
    # Reset seed before getting batch
    set_global_seed(seed)
    
    model.train()
    model.to(device)
    
    # Get first batch
    batch_iter = iter(dataloader)
    batch_data = next(batch_iter)
    
    # Handle different batch formats (with or without span_lengths)
    if len(batch_data) == 3:
        inputs, labels, _ = batch_data
    else:
        inputs, labels = batch_data
    
    inputs = inputs.to(device)
    labels = labels.to(device)
    
    # Forward pass
    logits, loss = model(inputs, labels)
    
    # Backward pass (to verify gradients are deterministic too)
    loss.backward()
    
    loss_value = loss.item()
    
    # Save
    save_golden_batch(
        inputs=inputs,
        labels=labels,
        loss=loss_value,
        seed=seed,
        path=path,
        metadata=metadata,
    )
    
    return inputs.cpu(), labels.cpu(), loss_value


def test_golden_batch(
    model: nn.Module,
    device: torch.device,
    path: Path,
    tolerance: float = 1e-5,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Test that a model produces the same loss on the golden batch.
    
    Args:
        model: Model to test
        device: Device to run on
        path: Path to golden batch file
        tolerance: Absolute tolerance for loss comparison
        
    Returns:
        Tuple of (passed, results_dict)
    """
    checkpoint = load_golden_batch(path)
    
    inputs = checkpoint["inputs"].to(device)
    labels = checkpoint["labels"].to(device)
    expected_loss = checkpoint["loss"]
    seed = checkpoint["seed"]
    
    # Reset seed for deterministic behavior
    from seed_utils import set_global_seed
    set_global_seed(seed)
    
    model.train()
    model.to(device)
    model.zero_grad()
    
    # Forward pass
    logits, loss = model(inputs, labels)
    
    # Backward pass
    loss.backward()
    
    actual_loss = loss.item()
    diff = abs(actual_loss - expected_loss)
    passed = diff <= tolerance
    
    results = {
        "passed": passed,
        "expected_loss": expected_loss,
        "actual_loss": actual_loss,
        "difference": diff,
        "tolerance": tolerance,
        "seed": seed,
    }
    
    if passed:
        LOGGER.info(
            "✓ Golden batch test PASSED (loss=%.6f, expected=%.6f, diff=%.2e)",
            actual_loss, expected_loss, diff
        )
    else:
        LOGGER.error(
            "✗ Golden batch test FAILED (loss=%.6f, expected=%.6f, diff=%.2e > tol=%.2e)",
            actual_loss, expected_loss, diff, tolerance
        )
    
    return passed, results


def compare_batch_compositions(
    dataloader: torch.utils.data.DataLoader,
    seed1: int,
    seed2: int,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Verify that different seeds produce different batch compositions.
    
    Args:
        dataloader: DataLoader to test
        seed1: First seed
        seed2: Second seed (should be different from seed1)
        
    Returns:
        Tuple of (batches_differ, comparison_results)
    """
    from seed_utils import set_global_seed
    
    # Get batch with seed1
    set_global_seed(seed1)
    batch_iter1 = iter(dataloader)
    batch1 = next(batch_iter1)
    if len(batch1) == 3:
        inputs1, labels1, _ = batch1
    else:
        inputs1, labels1 = batch1
    
    # Get batch with seed2
    set_global_seed(seed2)
    batch_iter2 = iter(dataloader)
    batch2 = next(batch_iter2)
    if len(batch2) == 3:
        inputs2, labels2, _ = batch2
    else:
        inputs2, labels2 = batch2
    
    # Compare
    inputs_differ = not torch.equal(inputs1, inputs2)
    labels_differ = not torch.equal(labels1, labels2)
    
    # At least one should differ if seeds are different
    batches_differ = inputs_differ or labels_differ
    
    results = {
        "seed1": seed1,
        "seed2": seed2,
        "inputs_differ": inputs_differ,
        "labels_differ": labels_differ,
        "batches_differ": batches_differ,
    }
    
    if batches_differ:
        LOGGER.info(
            "✓ Batch composition test PASSED (seeds %d vs %d produce different batches)",
            seed1, seed2
        )
    else:
        LOGGER.error(
            "✗ Batch composition test FAILED (seeds %d vs %d produce identical batches!)",
            seed1, seed2
        )
    
    return batches_differ, results


def verify_same_seed_same_batch(
    dataloader_factory,
    seed: int,
    num_trials: int = 3,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Verify that the same seed produces the same batch across multiple trials.
    
    Args:
        dataloader_factory: Callable that creates a fresh DataLoader
        seed: Seed to test
        num_trials: Number of times to regenerate and compare
        
    Returns:
        Tuple of (all_match, comparison_results)
    """
    from seed_utils import set_global_seed
    
    batches = []
    
    for trial in range(num_trials):
        set_global_seed(seed)
        loader = dataloader_factory()
        batch_iter = iter(loader)
        batch = next(batch_iter)
        
        if len(batch) == 3:
            inputs, labels, _ = batch
        else:
            inputs, labels = batch
        
        batches.append((inputs.clone(), labels.clone()))
    
    # Compare all batches to first
    reference_inputs, reference_labels = batches[0]
    all_match = True
    mismatches = []
    
    for trial, (inputs, labels) in enumerate(batches[1:], start=2):
        inputs_match = torch.equal(reference_inputs, inputs)
        labels_match = torch.equal(reference_labels, labels)
        
        if not (inputs_match and labels_match):
            all_match = False
            mismatches.append({
                "trial": trial,
                "inputs_match": inputs_match,
                "labels_match": labels_match,
            })
    
    results = {
        "seed": seed,
        "num_trials": num_trials,
        "all_match": all_match,
        "mismatches": mismatches,
    }
    
    if all_match:
        LOGGER.info(
            "✓ Same-seed consistency test PASSED (%d trials with seed %d)",
            num_trials, seed
        )
    else:
        LOGGER.error(
            "✗ Same-seed consistency test FAILED (%d mismatches in %d trials)",
            len(mismatches), num_trials
        )
    
    return all_match, results


def main() -> None:
    """CLI for golden batch generation and testing."""
    parser = argparse.ArgumentParser(description="Golden batch utilities")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--generate", 
        action="store_true", 
        help="Generate and save golden batch"
    )
    parser.add_argument(
        "--test", 
        action="store_true", 
        help="Test model against golden batch"
    )
    parser.add_argument(
        "--path",
        type=str,
        default=None,
        help="Path to golden batch file"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-5,
        help="Loss tolerance for testing"
    )
    args = parser.parse_args()
    
    if not args.generate and not args.test:
        parser.error("Must specify --generate or --test")
    
    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    
    # Determine golden batch path
    if args.path:
        golden_path = Path(args.path)
        if not golden_path.is_absolute():
            golden_path = script_dir / golden_path
    else:
        golden_path = DEFAULT_GOLDEN_BATCH_PATH.resolve()
    
    # Import here to avoid circular dependencies
    from seed_utils import set_global_seed
    from pretrain_mlm import prepare_dataset, MLMDataset, MLMCollator, DNAMLM
    from model import DNAEncoder
    from torch.utils.data import DataLoader, random_split
    
    set_global_seed(args.seed)
    
    # Prepare data
    token_tensors, vocab = prepare_dataset(cfg)
    full_dataset = MLMDataset(token_tensors)
    
    # Split dataset deterministically
    val_ratio = float(cfg.get("mlm_val_ratio", 0.1))
    val_size = int(len(full_dataset) * val_ratio)
    train_size = len(full_dataset) - val_size
    train_dataset, _ = random_split(
        full_dataset,
        [train_size, val_size],
        generator=torch.Generator().manual_seed(args.seed)
    )
    
    batch_size = int(cfg.get("mlm_batch_size", 128))
    collator = MLMCollator(vocab, mask_prob=0.15, track_spans=True)
    
    from seed_utils import get_deterministic_dataloader_kwargs
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=0,
        drop_last=True,
        **get_deterministic_dataloader_kwargs(args.seed),
    )
    
    # Build model
    device = torch.device("cpu")  # Use CPU for determinism
    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
    
    embedding_dim = int(cfg.get("mlm_embedding_dim", cfg.get("embedding_dim", 256)))
    transformer_layers = int(cfg.get("mlm_transformer_layers", cfg.get("transformer_layers", 8)))
    transformer_heads = int(cfg.get("mlm_transformer_heads", cfg.get("transformer_heads", 8)))
    transformer_ff_dim = int(cfg.get("mlm_transformer_ff_dim", cfg.get("transformer_ff_dim", embedding_dim * 4)))
    transformer_dropout = float(cfg.get("mlm_transformer_dropout", cfg.get("transformer_dropout", 0.1)))
    
    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        drop_path_rate=0.0,
    )
    
    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=len(vocab.itos),
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,  # BERT-style Dense+GELU+LN before projection
    )
    
    if args.generate:
        LOGGER.info("Generating golden batch with seed %d...", args.seed)
        metadata = {
            "config": str(config_path),
            "embedding_dim": embedding_dim,
            "transformer_layers": transformer_layers,
            "transformer_heads": transformer_heads,
            "batch_size": batch_size,
        }
        generate_golden_batch(
            dataloader=train_loader,
            model=model,
            device=device,
            seed=args.seed,
            path=golden_path,
            metadata=metadata,
        )
    
    if args.test:
        LOGGER.info("Testing against golden batch...")
        passed, results = test_golden_batch(
            model=model,
            device=device,
            path=golden_path,
            tolerance=args.tolerance,
        )
        
        if not passed:
            raise AssertionError(
                f"Golden batch test failed! Expected loss {results['expected_loss']:.6f}, "
                f"got {results['actual_loss']:.6f} (diff={results['difference']:.2e})"
            )


if __name__ == "__main__":
    main()
