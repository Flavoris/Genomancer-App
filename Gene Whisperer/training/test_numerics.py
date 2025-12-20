"""
Test suite for numerics.py NaN/Inf tripwire functionality.

Verifies that:
1. assert_finite passes for finite tensors
2. assert_finite raises NumericsError for NaN tensors
3. assert_finite raises NumericsError for Inf tensors
4. detect_nan_in_model catches NaN/Inf in gradients
5. NumericsChecker works in a simulated training loop
"""

import sys
from pathlib import Path

import torch
import torch.nn as nn

# Add training directory to path
sys.path.insert(0, str(Path(__file__).parent))

from numerics import (
    assert_finite,
    detect_nan_in_model,
    NumericsChecker,
    NumericsError,
    compute_total_grad_norm,
    compute_total_param_norm,
    cap_model_param_norm_,
)


def test_assert_finite_passes_for_finite_tensors():
    """assert_finite should not raise for finite tensors."""
    print("Test 1: assert_finite with finite tensors... ", end="")
    
    # Various finite tensors
    tensors = [
        torch.zeros(10),
        torch.ones(5, 5),
        torch.randn(3, 4, 5),
        torch.tensor([1.0, -1.0, 0.0, 1e-10, 1e10]),
        torch.tensor([[1.0]]),
    ]
    
    for i, t in enumerate(tensors):
        try:
            assert_finite(t, f"test_tensor_{i}")
        except NumericsError as e:
            print(f"FAILED - raised error for finite tensor: {e}")
            raise AssertionError(f"raised error for finite tensor: {e}") from e
    
    print("PASSED")


def test_assert_finite_raises_for_nan():
    """assert_finite should raise NumericsError for NaN tensors."""
    print("Test 2: assert_finite raises for NaN... ", end="")
    
    nan_tensor = torch.tensor([1.0, float('nan'), 3.0])
    
    try:
        assert_finite(nan_tensor, "nan_test")
    except NumericsError as e:
        # Verify error message contains expected info
        msg = str(e)
        assert "NaN" in msg, f"error message missing 'NaN': {msg}"
        assert "nan_test" in msg, f"error message missing tensor name: {msg}"
    else:
        print("FAILED - did not raise error")
        raise AssertionError("did not raise error")
    print("PASSED")


def test_assert_finite_raises_for_inf():
    """assert_finite should raise NumericsError for Inf tensors."""
    print("Test 3: assert_finite raises for Inf... ", end="")
    
    inf_tensor = torch.tensor([1.0, float('inf'), 3.0])
    
    try:
        assert_finite(inf_tensor, "inf_test")
    except NumericsError as e:
        msg = str(e)
        assert "Inf" in msg, f"error message missing 'Inf': {msg}"
    else:
        print("FAILED - did not raise error")
        raise AssertionError("did not raise error")
    print("PASSED")


def test_assert_finite_raises_for_neg_inf():
    """assert_finite should raise NumericsError for -Inf tensors."""
    print("Test 4: assert_finite raises for -Inf... ", end="")
    
    neg_inf_tensor = torch.tensor([1.0, float('-inf'), 3.0])
    
    try:
        assert_finite(neg_inf_tensor, "neg_inf_test", allow_neg_inf=False)
    except NumericsError as e:
        msg = str(e)
        assert "Inf" in msg, f"error message missing 'Inf': {msg}"
    else:
        print("FAILED - did not raise error")
        raise AssertionError("did not raise error")
    print("PASSED")


def test_assert_finite_with_context():
    """assert_finite should include context in error message."""
    print("Test 5: assert_finite includes context... ", end="")
    
    nan_tensor = torch.tensor([float('nan')])
    context = {
        "epoch": 5,
        "global_step": 1234,
        "learning_rate": 1e-4,
        "batch_idx": 42,
    }
    
    try:
        assert_finite(nan_tensor, "logits", context=context)
    except NumericsError as e:
        msg = str(e)
        # Check context is in message
        assert "epoch" in msg.lower() and "5" in msg, "missing epoch in message"
        assert "1234" in msg, "missing step in message"
    else:
        print("FAILED - did not raise error")
        raise AssertionError("did not raise error")
    print("PASSED")


def test_detect_nan_in_model_passes_for_clean_model():
    """detect_nan_in_model should not raise for clean gradients."""
    print("Test 6: detect_nan_in_model with clean model... ", end="")
    
    model = nn.Linear(10, 5)
    x = torch.randn(3, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()
    
    try:
        detect_nan_in_model(model)
    except NumericsError as e:
        print(f"FAILED - raised error for clean model: {e}")
        raise AssertionError(f"raised error for clean model: {e}") from e
    print("PASSED")


def test_detect_nan_in_model_catches_nan_grad():
    """detect_nan_in_model should catch NaN in gradients."""
    print("Test 7: detect_nan_in_model catches NaN gradient... ", end="")
    
    model = nn.Linear(10, 5)
    x = torch.randn(3, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()
    
    # Inject NaN into gradient
    model.weight.grad[0, 0] = float('nan')
    
    try:
        detect_nan_in_model(model)
    except NumericsError as e:
        msg = str(e)
        assert "GRAD NaN" in msg, f"error message missing 'GRAD NaN': {msg}"
        assert "weight" in msg, "error message missing parameter name"
        # Check for top-5 gradient norms
        assert "Top-5" in msg, "error message missing 'Top-5' gradient norms"
    else:
        print("FAILED - did not catch NaN gradient")
        raise AssertionError("did not catch NaN gradient")
    print("PASSED")


def test_detect_nan_in_model_catches_inf_param():
    """detect_nan_in_model should catch Inf in parameters."""
    print("Test 8: detect_nan_in_model catches Inf parameter... ", end="")
    
    model = nn.Linear(10, 5)
    
    # Inject Inf into parameter
    with torch.no_grad():
        model.weight[0, 0] = float('inf')
    
    try:
        detect_nan_in_model(model, check_params=True, check_grads=False)
    except NumericsError as e:
        msg = str(e)
        assert "PARAM Inf" in msg, f"error message missing 'PARAM Inf': {msg}"
    else:
        print("FAILED - did not catch Inf parameter")
        raise AssertionError("did not catch Inf parameter")
    print("PASSED")


def test_detect_nan_in_model_includes_context():
    """detect_nan_in_model should include training context in error."""
    print("Test 9: detect_nan_in_model includes context... ", end="")
    
    model = nn.Linear(10, 5)
    x = torch.randn(3, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()
    
    # Inject NaN
    model.weight.grad[0, 0] = float('nan')
    
    context = {
        "epoch": 10,
        "global_step": 5000,
        "learning_rate": 2e-5,
    }
    
    try:
        detect_nan_in_model(model, context=context)
    except NumericsError as e:
        msg = str(e)
        # Check context fields
        assert "epoch" in msg.lower(), "missing 'epoch' in context"
        assert "5000" in msg, "missing step value"
        assert "learning_rate" in msg.lower(), "missing 'learning_rate'"
    else:
        print("FAILED - did not raise error")
        raise AssertionError("did not raise error")
    print("PASSED")


def test_numerics_checker_class():
    """NumericsChecker should work in a simulated training loop."""
    print("Test 10: NumericsChecker in training loop... ", end="")
    
    model = nn.Linear(10, 5)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
    checker = NumericsChecker(enabled=True)
    
    # Simulate a few clean training steps
    for step in range(3):
        x = torch.randn(4, 10)
        
        checker.set_context(
            epoch=1,
            step=step,
            lr=0.01,
            batch_idx=step,
            model=model,
        )
        
        logits = model(x)
        loss = logits.sum()
        
        try:
            checker.check_forward(logits, loss)
        except NumericsError:
            print(f"FAILED - raised error on clean forward at step {step}")
            raise AssertionError(f"raised error on clean forward at step {step}")
        
        optimizer.zero_grad()
        loss.backward()
        
        try:
            checker.check_backward(model)
        except NumericsError:
            print(f"FAILED - raised error on clean backward at step {step}")
            raise AssertionError(f"raised error on clean backward at step {step}")
        
        optimizer.step()
    
    print("PASSED")


def test_numerics_checker_catches_nan_loss():
    """NumericsChecker should catch NaN loss."""
    print("Test 11: NumericsChecker catches NaN loss... ", end="")
    
    model = nn.Linear(10, 5)
    checker = NumericsChecker(enabled=True)
    
    checker.set_context(epoch=1, step=0, lr=0.01, model=model)
    
    logits = torch.randn(4, 5)
    loss = torch.tensor(float('nan'))
    
    try:
        checker.check_forward(logits, loss)
    except NumericsError as e:
        assert "loss" in str(e).lower(), "error not about loss"
    else:
        print("FAILED - did not catch NaN loss")
        raise AssertionError("did not catch NaN loss")
    print("PASSED")


def test_numerics_checker_disabled():
    """NumericsChecker with enabled=False should not check."""
    print("Test 12: NumericsChecker disabled mode... ", end="")
    
    model = nn.Linear(10, 5)
    checker = NumericsChecker(enabled=False)
    
    checker.set_context(epoch=1, step=0, lr=0.01, model=model)
    
    # These should NOT raise even with NaN
    nan_logits = torch.tensor([float('nan')])
    nan_loss = torch.tensor(float('nan'))
    
    try:
        checker.check_forward(nan_logits, nan_loss)
    except NumericsError:
        print("FAILED - disabled checker still raised error")
        raise AssertionError("disabled checker still raised error")
    print("PASSED")


def test_compute_norms():
    """Test gradient and parameter norm computation."""
    print("Test 13: compute_total_grad_norm and compute_total_param_norm... ", end="")
    
    model = nn.Linear(10, 5)
    x = torch.randn(3, 10)
    y = model(x)
    loss = y.sum()
    loss.backward()
    
    grad_norm = compute_total_grad_norm(model)
    param_norm = compute_total_param_norm(model)
    
    assert grad_norm > 0, f"grad_norm should be positive: {grad_norm}"
    assert param_norm > 0, f"param_norm should be positive: {param_norm}"
    
    print(f"PASSED (grad_norm={grad_norm:.4f}, param_norm={param_norm:.4f})")


def test_cap_model_param_norm():
    """cap_model_param_norm_ should rescale parameters when norm is too large."""
    print("Test 14: cap_model_param_norm_ rescales parameters... ", end="")
    
    class TinyModel(nn.Module):
        def __init__(self):
            super().__init__()
            self.w = nn.Parameter(torch.ones(256))
    
        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return x * self.w.sum()
    
    model = TinyModel()
    
    with torch.no_grad():
        model.w.fill_(1000.0)
    
    before = compute_total_param_norm(model)
    cap = 200.0
    total_norm, scale = cap_model_param_norm_(model, max_norm=cap)
    after = compute_total_param_norm(model)
    
    assert before > cap, f"setup norm should exceed cap: before={before:.4f}, cap={cap}"
    assert total_norm > cap, f"reported total_norm should exceed cap: total_norm={total_norm:.4f}, cap={cap}"
    assert scale < 1.0, f"expected scaling (<1), got scale={scale}"
    assert after <= cap * 1.001, f"norm not capped: after={after:.4f}, cap={cap}"
    
    print(f"PASSED (before={before:.2f}, after={after:.2f}, scale={scale:.6f})")


def run_all_tests():
    """Run all tests and report results."""
    print("=" * 60)
    print("NUMERICS.PY TEST SUITE")
    print("=" * 60)
    print()
    
    tests = [
        test_assert_finite_passes_for_finite_tensors,
        test_assert_finite_raises_for_nan,
        test_assert_finite_raises_for_inf,
        test_assert_finite_raises_for_neg_inf,
        test_assert_finite_with_context,
        test_detect_nan_in_model_passes_for_clean_model,
        test_detect_nan_in_model_catches_nan_grad,
        test_detect_nan_in_model_catches_inf_param,
        test_detect_nan_in_model_includes_context,
        test_numerics_checker_class,
        test_numerics_checker_catches_nan_loss,
        test_numerics_checker_disabled,
        test_compute_norms,
        test_cap_model_param_norm,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"EXCEPTION: {e}")
            failed += 1
    
    print()
    print("=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
