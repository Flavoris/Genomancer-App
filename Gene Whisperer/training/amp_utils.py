"""
amp_utils.py - Mixed Precision (AMP) utilities for MPS and CPU training.

Provides helper functions to create autocast contexts for Apple MPS devices
without requiring GradScaler (which is not supported on MPS).

Usage:
    from amp_utils import get_amp_context, log_amp_status

    # At startup
    amp_ctx = get_amp_context(device, cfg)
    log_amp_status(amp_ctx, logger)

    # In training loop
    with amp_ctx.autocast():
        outputs = model(inputs)
    # Loss computed in fp32 if amp_force_fp32_loss is True
    loss = criterion(outputs.float(), targets)
"""

import inspect
import logging
from contextlib import nullcontext
from dataclasses import dataclass
from typing import Any, ContextManager, Dict, Optional

import torch

LOGGER = logging.getLogger("gene_whisperer.amp_utils")

try:
    from torch.amp import GradScaler as AmpGradScaler
    _GRAD_SCALER_SOURCE = "torch.amp"
except Exception:  # pragma: no cover - fallback for older torch
    from torch.cuda.amp import GradScaler as AmpGradScaler
    _GRAD_SCALER_SOURCE = "torch.cuda.amp"

GradScalerType = AmpGradScaler


def _resolve_grad_scaler_device_kwarg() -> Optional[str]:
    """Detect the supported device kwarg for torch.amp.GradScaler across versions."""
    if _GRAD_SCALER_SOURCE != "torch.amp":
        return None
    try:
        params = inspect.signature(AmpGradScaler).parameters
    except (TypeError, ValueError):
        return None
    if "device" in params:
        return "device"
    if "device_type" in params:
        return "device_type"
    return None


_GRAD_SCALER_DEVICE_KWARG = _resolve_grad_scaler_device_kwarg()


@dataclass
class AMPContext:
    """Container for AMP settings and autocast context."""

    enabled: bool
    device_type: str
    dtype: Optional[torch.dtype]
    force_fp32_loss: bool

    def autocast(self) -> ContextManager:
        """Return the appropriate autocast context manager."""
        if not self.enabled or self.dtype is None:
            return nullcontext()
        return torch.autocast(self.device_type, dtype=self.dtype)

    def to_dict(self) -> Dict[str, Any]:
        """Return dict representation for logging/config."""
        dtype_str = None
        if self.dtype is not None:
            dtype_str = str(self.dtype).replace("torch.", "")
        return {
            "enabled": self.enabled,
            "device_type": self.device_type,
            "dtype": dtype_str,
            "force_fp32_loss": self.force_fp32_loss,
        }


def _parse_amp_dtype(dtype_str: str) -> torch.dtype:
    """Parse dtype string to torch.dtype."""
    dtype_map = {
        "float16": torch.float16,
        "fp16": torch.float16,
        "bfloat16": torch.bfloat16,
        "bf16": torch.bfloat16,
    }
    dtype_lower = dtype_str.lower().strip()
    if dtype_lower not in dtype_map:
        raise ValueError(
            f"Invalid amp_dtype: {dtype_str!r}. Must be one of: float16, bfloat16"
        )
    return dtype_map[dtype_lower]


def _test_mps_autocast(dtype: torch.dtype) -> bool:
    """Test if MPS autocast works with the given dtype."""
    try:
        with torch.autocast("mps", dtype=dtype):
            x = torch.zeros(1, device="mps", dtype=torch.float32)
            _ = x + x
        return True
    except Exception:
        return False


def get_amp_context(
    device: torch.device,
    cfg: Dict[str, Any],
    *,
    logger: Optional[logging.Logger] = None,
) -> AMPContext:
    """
    Create an AMP context based on device and config.

    Args:
        device: The torch device being used
        cfg: Configuration dict with amp_enabled, amp_dtype, amp_force_fp32_loss
        logger: Optional logger for warnings

    Returns:
        AMPContext with appropriate settings for the device
    """
    log = logger or LOGGER

    amp_enabled = bool(cfg.get("amp_enabled", True))
    amp_dtype_str = str(cfg.get("amp_dtype", "bfloat16"))
    force_fp32_loss = bool(cfg.get("amp_force_fp32_loss", True))

    # Default: disabled
    if not amp_enabled:
        return AMPContext(
            enabled=False,
            device_type=device.type,
            dtype=None,
            force_fp32_loss=force_fp32_loss,
        )

    # CPU: autocast not beneficial, disable
    if device.type == "cpu":
        log.info("AMP: disabled on CPU (no speedup benefit)")
        return AMPContext(
            enabled=False,
            device_type="cpu",
            dtype=None,
            force_fp32_loss=force_fp32_loss,
        )

    # Parse requested dtype
    try:
        requested_dtype = _parse_amp_dtype(amp_dtype_str)
    except ValueError as e:
        log.warning("AMP: %s, falling back to bfloat16", e)
        requested_dtype = torch.bfloat16

    # MPS device
    if device.type == "mps":
        # bfloat16 is preferred on MPS (wider dynamic range, no GradScaler needed)
        # float16 works but may have more numerical issues
        if _test_mps_autocast(requested_dtype):
            return AMPContext(
                enabled=True,
                device_type="mps",
                dtype=requested_dtype,
                force_fp32_loss=force_fp32_loss,
            )
        else:
            # Try fallback to bfloat16 if float16 was requested
            if requested_dtype == torch.float16:
                log.warning(
                    "AMP: float16 autocast not available on MPS, trying bfloat16"
                )
                if _test_mps_autocast(torch.bfloat16):
                    return AMPContext(
                        enabled=True,
                        device_type="mps",
                        dtype=torch.bfloat16,
                        force_fp32_loss=force_fp32_loss,
                    )
            log.warning("AMP: autocast not available on this MPS device, disabling")
            return AMPContext(
                enabled=False,
                device_type="mps",
                dtype=None,
                force_fp32_loss=force_fp32_loss,
            )

    # CUDA device
    if device.type == "cuda":
        # Note: For CUDA, GradScaler is typically used with float16
        # This helper provides autocast context; GradScaler should be handled separately
        return AMPContext(
            enabled=True,
            device_type="cuda",
            dtype=requested_dtype,
            force_fp32_loss=force_fp32_loss,
        )

    # Unknown device type
    log.warning("AMP: unknown device type %s, disabling", device.type)
    return AMPContext(
        enabled=False,
        device_type=device.type,
        dtype=None,
        force_fp32_loss=force_fp32_loss,
    )


def log_amp_status(
    amp_ctx: AMPContext,
    logger: Optional[logging.Logger] = None,
) -> None:
    """
    Log AMP status at startup.

    Args:
        amp_ctx: The AMP context to log
        logger: Logger to use (defaults to module logger)
    """
    log = logger or LOGGER
    dtype_str = "N/A"
    if amp_ctx.dtype is not None:
        dtype_str = str(amp_ctx.dtype).replace("torch.", "")

    log.info(
        "AMP: enabled=%s dtype=%s device=%s force_fp32_loss=%s",
        amp_ctx.enabled,
        dtype_str,
        amp_ctx.device_type,
        amp_ctx.force_fp32_loss,
    )


def create_grad_scaler(device_type: str) -> Optional[AmpGradScaler]:
    """
    Create a CUDA GradScaler using the non-deprecated torch.amp API when available.

    Args:
        device_type: Torch device type ("cuda", "mps", "cpu", ...)

    Returns:
        A GradScaler for CUDA training, or None when scaling is not applicable.
    """
    if device_type != "cuda":
        return None
    if _GRAD_SCALER_DEVICE_KWARG:
        return AmpGradScaler(**{_GRAD_SCALER_DEVICE_KWARG: device_type})
    return AmpGradScaler()
