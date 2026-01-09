from __future__ import annotations

import warnings
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from amp_utils import create_grad_scaler, GradScalerType


def test_create_grad_scaler_returns_none_for_non_cuda() -> None:
    assert create_grad_scaler("cpu") is None
    assert create_grad_scaler("mps") is None


def test_create_grad_scaler_for_cuda_returns_instance() -> None:
    # Suppress warnings when CUDA is unavailable in CPU-only test environments.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="torch.cuda.amp.GradScaler is enabled, but CUDA is not available.*",
        )
        scaler = create_grad_scaler("cuda")

    assert scaler is not None
    assert isinstance(scaler, GradScalerType)
