"""Deterministic reproducibility utilities for Gene Whisperer.

This module provides comprehensive seed setting for all random sources:
- Python's random module
- NumPy's random generators  
- PyTorch (CPU, CUDA, MPS)
- DataLoader worker seeds

Usage:
    from seed_utils import set_global_seed, get_deterministic_dataloader_kwargs
    
    set_global_seed(42)
    loader = DataLoader(dataset, **get_deterministic_dataloader_kwargs(42))
"""
from __future__ import annotations

import logging
import os
import random
from typing import Any, Callable, Dict, Optional

import numpy as np
import torch

LOGGER = logging.getLogger("gene_whisperer.seed_utils")


def set_global_seed(seed: int) -> None:
    """
    Set seeds for all random number generators to ensure reproducibility.
    
    Covers:
    - Python's random module
    - NumPy's random generators
    - PyTorch CPU random
    - PyTorch CUDA random (all GPUs)
    - PyTorch MPS random (Apple Silicon)
    - CUDNN deterministic settings
    
    Args:
        seed: Integer seed value (recommended range: 0 to 2^32-1)
    """
    # Python random
    random.seed(seed)
    
    # NumPy random (both legacy and new-style)
    np.random.seed(seed)
    
    # PyTorch CPU
    torch.manual_seed(seed)
    
    # PyTorch CUDA (all GPUs)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    
    # PyTorch MPS (Apple Silicon)
    if hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
        try:
            torch.mps.manual_seed(seed)
        except Exception:
            pass  # MPS seed may not be available in all PyTorch versions
    
    # CUDNN deterministic settings for reproducibility
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    
    # Set environment variable for any child processes
    os.environ["PYTHONHASHSEED"] = str(seed)
    os.environ["GENE_WHISPERER_SEED"] = str(seed)
    
    LOGGER.debug("Global seed set to %d", seed)


def get_worker_init_fn(base_seed: int) -> Callable[[int], None]:
    """
    Create a worker_init_fn for DataLoader that ensures each worker
    has a deterministic but unique seed.
    
    This is critical for reproducibility when using num_workers > 0.
    Each worker gets: base_seed + worker_id
    
    Args:
        base_seed: Base seed value
        
    Returns:
        A function suitable for DataLoader's worker_init_fn parameter
    """
    def worker_init_fn(worker_id: int) -> None:
        worker_seed = base_seed + worker_id
        random.seed(worker_seed)
        np.random.seed(worker_seed)
        torch.manual_seed(worker_seed)
    
    return worker_init_fn


def get_deterministic_generator(seed: int) -> torch.Generator:
    """
    Create a deterministic PyTorch Generator for use with DataLoader.
    
    Args:
        seed: Seed value for the generator
        
    Returns:
        A torch.Generator initialized with the given seed
    """
    generator = torch.Generator()
    generator.manual_seed(seed)
    return generator


def get_deterministic_dataloader_kwargs(
    seed: int,
    num_workers: int = 0,
) -> Dict[str, Any]:
    """
    Get DataLoader kwargs for deterministic behavior.
    
    Returns a dict with:
    - generator: Seeded Generator for shuffling
    - worker_init_fn: Per-worker seed initialization
    
    Usage:
        loader = DataLoader(
            dataset,
            shuffle=True,
            num_workers=4,
            **get_deterministic_dataloader_kwargs(42, num_workers=4)
        )
    
    Args:
        seed: Base seed for reproducibility
        num_workers: Number of DataLoader workers
        
    Returns:
        Dict of kwargs to pass to DataLoader
    """
    kwargs: Dict[str, Any] = {
        "generator": get_deterministic_generator(seed),
    }
    
    if num_workers > 0:
        kwargs["worker_init_fn"] = get_worker_init_fn(seed)
    
    return kwargs


def verify_reproducibility(seed: int, num_checks: int = 5) -> bool:
    """
    Quick verification that seeding is working correctly.
    
    Sets seed, generates random values, resets seed, generates again,
    and verifies they match.
    
    Args:
        seed: Seed to test
        num_checks: Number of random values to compare
        
    Returns:
        True if reproducibility is verified, False otherwise
    """
    # First pass
    set_global_seed(seed)
    python_vals1 = [random.random() for _ in range(num_checks)]
    numpy_vals1 = np.random.rand(num_checks).tolist()
    torch_vals1 = torch.rand(num_checks).tolist()
    
    # Second pass (reset and regenerate)
    set_global_seed(seed)
    python_vals2 = [random.random() for _ in range(num_checks)]
    numpy_vals2 = np.random.rand(num_checks).tolist()
    torch_vals2 = torch.rand(num_checks).tolist()
    
    # Verify match
    python_match = python_vals1 == python_vals2
    numpy_match = numpy_vals1 == numpy_vals2
    torch_match = torch_vals1 == torch_vals2
    
    all_match = python_match and numpy_match and torch_match
    
    if not all_match:
        LOGGER.warning(
            "Reproducibility check failed: python=%s, numpy=%s, torch=%s",
            python_match, numpy_match, torch_match
        )
    
    return all_match

