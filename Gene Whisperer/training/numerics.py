"""
numerics.py - NaN/Inf detection and numerical stability helpers.

Provides "tripwire" functions that fail fast when non-finite values
appear in tensors, gradients, or model parameters during training.
This prevents silent NaN/Inf propagation that can corrupt training.

Usage:
    from numerics import assert_finite, detect_nan_in_model, NumericsError

    # In training loop:
    logits, loss = model(inputs, labels)
    assert_finite(logits, "logits")
    assert_finite(loss, "loss")
    
    loss.backward()
    detect_nan_in_model(model, context={"epoch": epoch, "step": step, ...})
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn

LOGGER = logging.getLogger("gene_whisperer.numerics")


class NumericsError(RuntimeError):
    """
    Raised when non-finite values (NaN or Inf) are detected in tensors.
    
    Contains detailed diagnostics including:
    - Where the issue was detected
    - Training context (epoch, step, lr, etc.)
    - Model parameter norms
    - Top gradient norms
    """
    pass


@dataclass
class GradientInfo:
    """Information about a parameter's gradient."""
    name: str
    grad_norm: float
    param_norm: float
    has_nan: bool
    has_inf: bool
    shape: Tuple[int, ...]


def assert_finite(
    tensor: torch.Tensor,
    name: str,
    context: Optional[Dict[str, Any]] = None,
    allow_neg_inf: bool = True,
) -> None:
    """
    Assert that a tensor contains only finite values (no NaN or problematic Inf).
    
    Args:
        tensor: The tensor to check
        name: Name of the tensor (for error messages)
        context: Optional dict with training context (epoch, step, lr, etc.)
        allow_neg_inf: If True, allows -inf (used for masking logits). Default True.
        
    Raises:
        NumericsError: If tensor contains NaN or positive Inf values
        
    Note:
        Negative infinity (-inf) is allowed by default because it's commonly used
        to mask out logits for tokens that should never be predicted (e.g., special
        tokens like [MASK], [PAD], [UNK]). After softmax, -inf becomes 0 probability.
    """
    if tensor is None:
        return
    
    has_nan = torch.isnan(tensor).any().item()
    
    # Check for positive infinity (always problematic)
    pos_inf_mask = torch.isinf(tensor) & (tensor > 0)
    has_pos_inf = pos_inf_mask.any().item()
    
    # Check for negative infinity (only if not allowed)
    neg_inf_mask = torch.isinf(tensor) & (tensor < 0)
    has_neg_inf = neg_inf_mask.any().item() if not allow_neg_inf else False
    
    has_problematic_inf = has_pos_inf or has_neg_inf
    
    if has_nan or has_problematic_inf:
        # Count non-finite values
        nan_count = torch.isnan(tensor).sum().item()
        pos_inf_count = pos_inf_mask.sum().item()
        neg_inf_count = neg_inf_mask.sum().item() if not allow_neg_inf else 0
        inf_count = pos_inf_count + neg_inf_count
        total_elements = tensor.numel()
        
        # Build diagnostic message
        issues = []
        if has_nan:
            issues.append(f"NaN ({nan_count}/{total_elements} elements)")
        if has_pos_inf:
            issues.append(f"+Inf ({pos_inf_count}/{total_elements} elements)")
        if has_neg_inf:
            issues.append(f"-Inf ({neg_inf_count}/{total_elements} elements)")
        
        msg_parts = [
            f"\n{'='*70}",
            f"NUMERICAL INSTABILITY DETECTED in '{name}'",
            f"{'='*70}",
            f"Issues: {', '.join(issues)}",
            f"Tensor shape: {tuple(tensor.shape)}",
            f"Tensor dtype: {tensor.dtype}",
        ]
        
        # Add finite value statistics if any exist
        finite_mask = torch.isfinite(tensor)
        if finite_mask.any():
            finite_vals = tensor[finite_mask]
            msg_parts.extend([
                f"Finite value stats:",
                f"  Min: {finite_vals.min().item():.6e}",
                f"  Max: {finite_vals.max().item():.6e}",
                f"  Mean: {finite_vals.mean().item():.6e}",
                f"  Std: {finite_vals.std().item():.6e}",
            ])
        else:
            msg_parts.append("No finite values in tensor!")
        
        # Add context if provided
        if context:
            msg_parts.append("\nTraining context:")
            for key, value in context.items():
                if isinstance(value, float):
                    msg_parts.append(f"  {key}: {value:.6e}")
                else:
                    msg_parts.append(f"  {key}: {value}")
        
        msg_parts.append("="*70)
        
        raise NumericsError("\n".join(msg_parts))


def compute_gradient_info(model: nn.Module) -> List[GradientInfo]:
    """
    Collect gradient information for all parameters with gradients.
    
    Args:
        model: The model to inspect
        
    Returns:
        List of GradientInfo for each parameter with a gradient
    """
    grad_infos = []
    
    for name, param in model.named_parameters():
        if param.grad is not None:
            grad = param.grad
            grad_infos.append(GradientInfo(
                name=name,
                grad_norm=grad.norm().item(),
                param_norm=param.norm().item(),
                has_nan=torch.isnan(grad).any().item(),
                has_inf=torch.isinf(grad).any().item(),
                shape=tuple(grad.shape),
            ))
    
    return grad_infos


def compute_total_grad_norm(model: nn.Module, norm_type: float = 2.0) -> float:
    """
    Compute the total gradient norm across all parameters.
    
    Args:
        model: The model
        norm_type: Type of norm (default L2)
        
    Returns:
        Total gradient norm, or 0.0 if no gradients
    """
    total_norm = 0.0
    for param in model.parameters():
        if param.grad is not None:
            total_norm += param.grad.norm(norm_type).item() ** norm_type
    
    return total_norm ** (1.0 / norm_type) if total_norm > 0 else 0.0


def compute_total_param_norm(model: nn.Module, norm_type: float = 2.0) -> float:
    """
    Compute the total parameter norm across all parameters.
    
    Args:
        model: The model
        norm_type: Type of norm (default L2)
        
    Returns:
        Total parameter norm
    """
    total_norm = 0.0
    for param in model.parameters():
        total_norm += param.norm(norm_type).item() ** norm_type
    
    return total_norm ** (1.0 / norm_type) if total_norm > 0 else 0.0


def cap_model_param_norm_(
    model: nn.Module,
    max_norm: float,
    *,
    norm_type: float = 2.0,
    eps: float = 1e-12,
) -> Tuple[float, float]:
    """
    Cap model parameter norm by scaling all parameters in-place.

    This is a numerical safety valve (not a replacement for proper optimization
    settings). It computes the total parameter norm and, if it exceeds
    `max_norm`, rescales all parameters by `max_norm / total_norm`.

    Returns:
        (total_norm_before, scale_applied). `scale_applied` is 1.0 if no scaling
        occurred, otherwise < 1.0.
    """
    if max_norm <= 0:
        return compute_total_param_norm(model, norm_type=norm_type), 1.0

    with torch.no_grad():
        total_norm = compute_total_param_norm(model, norm_type=norm_type)
        if total_norm > max_norm:
            scale = max_norm / max(total_norm, eps)
            for param in model.parameters():
                param.mul_(scale)
            return total_norm, scale
        return total_norm, 1.0


def detect_nan_in_model(
    model: nn.Module,
    context: Optional[Dict[str, Any]] = None,
    check_params: bool = True,
    check_grads: bool = True,
) -> None:
    """
    Scan model parameters and gradients for NaN/Inf values.
    
    If any non-finite values are found, raises NumericsError with a
    detailed diagnostic report including:
    - Which parameters/gradients are affected
    - Training context (epoch, step, lr, etc.)
    - Total grad norm and param norm
    - Top-5 largest gradient norms (to identify exploding gradients)
    
    Args:
        model: The model to scan
        context: Optional dict with training context (epoch, step, lr, etc.)
        check_params: Whether to check parameter values
        check_grads: Whether to check gradient values
        
    Raises:
        NumericsError: If any NaN/Inf values are found
    """
    issues = []
    
    # Check parameters
    if check_params:
        for name, param in model.named_parameters():
            if torch.isnan(param).any():
                issues.append(f"PARAM NaN: {name}")
            if torch.isinf(param).any():
                issues.append(f"PARAM Inf: {name}")
    
    # Check gradients
    grad_infos = []
    if check_grads:
        for name, param in model.named_parameters():
            if param.grad is not None:
                grad = param.grad
                grad_info = GradientInfo(
                    name=name,
                    grad_norm=grad.norm().item(),
                    param_norm=param.norm().item(),
                    has_nan=torch.isnan(grad).any().item(),
                    has_inf=torch.isinf(grad).any().item(),
                    shape=tuple(grad.shape),
                )
                grad_infos.append(grad_info)
                
                if grad_info.has_nan:
                    issues.append(f"GRAD NaN: {name}")
                if grad_info.has_inf:
                    issues.append(f"GRAD Inf: {name}")
    
    # If no issues, return early
    if not issues:
        return
    
    # Build detailed error report
    total_grad_norm = compute_total_grad_norm(model)
    total_param_norm = compute_total_param_norm(model)
    
    # Sort by gradient norm descending to find top-5
    grad_infos_sorted = sorted(grad_infos, key=lambda x: x.grad_norm, reverse=True)
    top_5_grads = grad_infos_sorted[:5]
    
    msg_parts = [
        f"\n{'='*70}",
        "NUMERICAL INSTABILITY DETECTED IN MODEL",
        f"{'='*70}",
        "",
        f"Issues found: {len(issues)}",
    ]
    
    for issue in issues[:10]:  # Show first 10 issues
        msg_parts.append(f"  • {issue}")
    if len(issues) > 10:
        msg_parts.append(f"  ... and {len(issues) - 10} more")
    
    msg_parts.extend([
        "",
        "Gradient Statistics:",
        f"  Total grad norm: {total_grad_norm:.6e}",
        f"  Total param norm: {total_param_norm:.6e}",
        "",
        "Top-5 Largest Gradient Norms:",
    ])
    
    for i, gi in enumerate(top_5_grads, 1):
        status = ""
        if gi.has_nan:
            status += " [NaN!]"
        if gi.has_inf:
            status += " [Inf!]"
        msg_parts.append(
            f"  {i}. {gi.name}: grad_norm={gi.grad_norm:.6e}, "
            f"param_norm={gi.param_norm:.6e}{status}"
        )
    
    # Add context if provided
    if context:
        msg_parts.extend(["", "Training Context:"])
        for key, value in context.items():
            if isinstance(value, float):
                msg_parts.append(f"  {key}: {value:.6e}")
            else:
                msg_parts.append(f"  {key}: {value}")
    
    msg_parts.extend([
        "",
        "Suggested fixes:",
        "  1. Reduce learning rate",
        "  2. Add/increase gradient clipping",
        "  3. Check for loss explosion (label smoothing, weight decay)",
        "  4. Verify input data is properly normalized",
        "  5. Check for division by zero in loss computation",
        f"{'='*70}",
    ])
    
    raise NumericsError("\n".join(msg_parts))


def create_training_context(
    epoch: int,
    step: int,
    lr: float,
    model: nn.Module,
    batch_idx: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Create a context dictionary for error reporting.
    
    Args:
        epoch: Current epoch number
        step: Global training step
        lr: Current learning rate
        model: The model (for computing norms)
        batch_idx: Optional batch index within epoch
        
    Returns:
        Dict with training context for error messages
    """
    context = {
        "epoch": epoch,
        "global_step": step,
        "learning_rate": lr,
        "grad_norm": compute_total_grad_norm(model),
        "param_norm": compute_total_param_norm(model),
    }
    
    if batch_idx is not None:
        context["batch_idx"] = batch_idx
    
    return context


class NumericsChecker:
    """
    Stateful numerics checker for use in training loops.
    
    Provides convenient methods to check tensors and model state
    with automatic context tracking.
    
    Usage:
        checker = NumericsChecker()
        
        for epoch in range(epochs):
            for batch_idx, (inputs, labels) in enumerate(loader):
                checker.set_context(epoch=epoch, step=step, lr=lr, batch_idx=batch_idx)
                
                logits, loss = model(inputs, labels)
                checker.check_forward(logits, loss)
                
                loss.backward()
                checker.check_backward(model)
                
                optimizer.step()
    """
    
    def __init__(self, enabled: bool = True):
        """
        Initialize the numerics checker.
        
        Args:
            enabled: If False, all checks become no-ops (for performance)
        """
        self.enabled = enabled
        self._context: Dict[str, Any] = {}
        self._model: Optional[nn.Module] = None
    
    def set_context(
        self,
        epoch: int,
        step: int,
        lr: float,
        batch_idx: Optional[int] = None,
        model: Optional[nn.Module] = None,
        **extra: Any,
    ) -> None:
        """
        Update the current training context.
        
        Args:
            epoch: Current epoch
            step: Global step
            lr: Current learning rate
            batch_idx: Batch index within epoch
            model: Model reference (for norm computations)
            **extra: Additional context values
        """
        self._context = {
            "epoch": epoch,
            "global_step": step,
            "learning_rate": lr,
        }
        if batch_idx is not None:
            self._context["batch_idx"] = batch_idx
        self._context.update(extra)
        
        if model is not None:
            self._model = model
    
    def check_forward(
        self,
        logits: torch.Tensor,
        loss: torch.Tensor,
    ) -> None:
        """
        Check forward pass outputs for NaN/Inf.
        
        Args:
            logits: Model output logits
            loss: Computed loss value
        """
        if not self.enabled:
            return
        
        context = self._get_context_with_norms()
        assert_finite(logits, "logits", context=context)
        assert_finite(loss, "loss", context=context)
    
    def check_backward(
        self,
        model: nn.Module,
        *,
        check_params: bool = True,
        check_grads: bool = True,
    ) -> None:
        """
        Check gradients after backward pass for NaN/Inf.
        
        Args:
            model: The model with computed gradients
            check_params: Whether to check parameter values
            check_grads: Whether to check gradient values
        """
        if not self.enabled:
            return
        
        context = self._get_context_with_norms()
        detect_nan_in_model(model, context=context, check_params=check_params, check_grads=check_grads)
    
    def check_tensor(self, tensor: torch.Tensor, name: str) -> None:
        """
        Check an arbitrary tensor for NaN/Inf.
        
        Args:
            tensor: Tensor to check
            name: Name for error messages
        """
        if not self.enabled:
            return
        
        context = self._get_context_with_norms()
        assert_finite(tensor, name, context=context)
    
    def _get_context_with_norms(self) -> Dict[str, Any]:
        """Get context dict with model norms if model is set."""
        context = dict(self._context)
        if self._model is not None:
            context["grad_norm"] = compute_total_grad_norm(self._model)
            context["param_norm"] = compute_total_param_norm(self._model)
        return context


# Convenience function for quick checks without a NumericsChecker instance
def check_training_step(
    model: nn.Module,
    logits: torch.Tensor,
    loss: torch.Tensor,
    epoch: int,
    step: int,
    lr: float,
    batch_idx: Optional[int] = None,
    after_backward: bool = False,
) -> None:
    """
    One-shot check for a training step.
    
    Call after forward to check logits and loss.
    Call after backward (with after_backward=True) to also check gradients.
    
    Args:
        model: The model
        logits: Output logits from forward pass
        loss: Computed loss
        epoch: Current epoch
        step: Global step
        lr: Current learning rate
        batch_idx: Optional batch index
        after_backward: If True, also check gradients
    """
    context = {
        "epoch": epoch,
        "global_step": step,
        "learning_rate": lr,
        "grad_norm": compute_total_grad_norm(model),
        "param_norm": compute_total_param_norm(model),
    }
    if batch_idx is not None:
        context["batch_idx"] = batch_idx
    
    assert_finite(logits, "logits", context=context)
    assert_finite(loss, "loss", context=context)
    
    if after_backward:
        detect_nan_in_model(model, context=context)
