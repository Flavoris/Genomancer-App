"""Unit tests for weight tying verification.

Tests ensure:
1. LM head uses the exact same Parameter object as token embedding
2. Optimizer sees tied weights only once (no double-optimization)
3. Weight norms remain stable during training
"""
from __future__ import annotations

import pytest
import torch
import torch.nn as nn
from torch.optim import AdamW

# Import the modules under test
from model import DNAEncoder
from pretrain_mlm import DNAMLM, get_layer_wise_lr_groups


def create_test_model(
    vocab_size: int = 67,
    embedding_dim: int = 64,
    num_layers: int = 2,
    num_heads: int = 4,
    ff_dim: int = 128,
    tie_weights: bool = True,
) -> DNAMLM:
    """Create a small test model for unit tests."""
    encoder = DNAEncoder(
        vocab_size=vocab_size,
        kmer=3,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.1,
        use_alibi=True,
        pad_token_id=vocab_size - 1,  # Assume last token is PAD
    )
    
    model = DNAMLM(
        encoder=encoder,
        vocab_size=vocab_size,
        special_token_ids=[vocab_size - 1, vocab_size - 2, vocab_size - 3],
        tie_weights=tie_weights,
        use_output_norm=True,
    )
    
    return model


class TestWeightTying:
    """Tests for proper weight tying between embedding and LM head."""
    
    def test_weight_tying_same_parameter_object(self):
        """
        CRITICAL TEST: Verify LM head uses the exact same Parameter object
        as token embedding.
        
        This is the core requirement for weight tying - they must be the
        SAME object, not just equal values.
        """
        model = create_test_model(tie_weights=True)
        
        # Get the Parameter objects
        embedding_weight = model.encoder.embedding.weight
        lm_head_weight = model.lm_head.weight
        
        # CRITICAL: They must be the EXACT SAME OBJECT
        assert id(embedding_weight) == id(lm_head_weight), (
            f"Weight tying FAILED: embedding.weight id={id(embedding_weight)} "
            f"!= lm_head.weight id={id(lm_head_weight)}. "
            "These must be the exact same Parameter object."
        )
        
        # Also verify they have the same data_ptr (underlying storage)
        assert embedding_weight.data_ptr() == lm_head_weight.data_ptr(), (
            "Weight tying FAILED: underlying storage is different!"
        )
        
        print("✓ Weight tying verified: embedding.weight is lm_head.weight")
    
    def test_weight_tying_disabled(self):
        """Verify weights are NOT tied when tie_weights=False."""
        model = create_test_model(tie_weights=False)
        
        embedding_weight = model.encoder.embedding.weight
        lm_head_weight = model.lm_head.weight
        
        # They should be DIFFERENT objects
        assert id(embedding_weight) != id(lm_head_weight), (
            "tie_weights=False but weights are still the same object!"
        )
        
        print("✓ Verified: weights are NOT tied when tie_weights=False")
    
    def test_weight_modification_propagates(self):
        """
        Verify that modifying one weight modifies the other (they are tied).
        """
        model = create_test_model(tie_weights=True)
        
        # Store original value
        original_val = model.encoder.embedding.weight[0, 0].item()
        
        # Modify via lm_head
        with torch.no_grad():
            model.lm_head.weight[0, 0] = 12345.0
        
        # Verify change propagated to embedding
        new_val = model.encoder.embedding.weight[0, 0].item()
        assert new_val == 12345.0, (
            f"Modification did not propagate! embedding.weight[0,0] = {new_val}"
        )
        
        # Modify via embedding
        with torch.no_grad():
            model.encoder.embedding.weight[0, 0] = -9999.0
        
        # Verify change propagated to lm_head
        lm_val = model.lm_head.weight[0, 0].item()
        assert lm_val == -9999.0, (
            f"Modification did not propagate! lm_head.weight[0,0] = {lm_val}"
        )
        
        print("✓ Weight modifications propagate correctly (bidirectionally)")


class TestOptimizerParameterGroups:
    """Tests for optimizer parameter deduplication with tied weights."""
    
    def test_optimizer_sees_tied_weight_once(self):
        """
        CRITICAL TEST: Verify optimizer sees the tied weight parameter
        only ONCE, not twice.
        
        Double-optimization would cause the tied weight to receive 2x
        the gradient updates, leading to instability.
        """
        model = create_test_model(tie_weights=True)
        
        # Get parameter groups using the production function
        param_groups = get_layer_wise_lr_groups(
            model=model,
            base_lr=1e-4,
            lr_decay=0.9,
            weight_decay=0.01,
        )
        
        # Collect all parameter IDs in the optimizer
        param_ids = []
        for group in param_groups:
            for param in group['params']:
                param_ids.append(id(param))
        
        # Get the tied weight's ID
        tied_weight_id = id(model.encoder.embedding.weight)
        
        # Count how many times the tied weight appears
        tied_weight_count = param_ids.count(tied_weight_id)
        
        assert tied_weight_count == 1, (
            f"Tied weight appears {tied_weight_count} times in optimizer param groups! "
            "Must appear exactly once. Double-optimization causes instability."
        )
        
        print(f"✓ Optimizer sees tied weight exactly once (id={tied_weight_id})")
    
    def test_no_duplicate_parameters_in_optimizer(self):
        """Verify NO parameter appears more than once in optimizer groups."""
        model = create_test_model(tie_weights=True)
        
        param_groups = get_layer_wise_lr_groups(
            model=model,
            base_lr=1e-4,
            lr_decay=0.9,
            weight_decay=0.01,
        )
        
        # Collect all parameter IDs
        param_ids = []
        for group in param_groups:
            for param in group['params']:
                param_ids.append(id(param))
        
        # Check for duplicates
        unique_ids = set(param_ids)
        
        if len(param_ids) != len(unique_ids):
            # Find the duplicates
            from collections import Counter
            counts = Counter(param_ids)
            duplicates = [(pid, cnt) for pid, cnt in counts.items() if cnt > 1]
            
            # Try to identify which parameters are duplicated
            id_to_name = {}
            for name, param in model.named_parameters():
                id_to_name[id(param)] = name
            
            dup_names = [f"{id_to_name.get(pid, 'unknown')} ({cnt}x)" 
                        for pid, cnt in duplicates]
            
            pytest.fail(
                f"Found duplicate parameters in optimizer! "
                f"Total params: {len(param_ids)}, Unique: {len(unique_ids)}. "
                f"Duplicates: {dup_names}"
            )
        
        print(f"✓ No duplicate parameters ({len(unique_ids)} unique params)")
    
    def test_all_model_parameters_covered(self):
        """
        Verify all trainable model parameters are in the optimizer,
        accounting for weight tying (tied weight should appear once).
        """
        model = create_test_model(tie_weights=True)
        
        param_groups = get_layer_wise_lr_groups(
            model=model,
            base_lr=1e-4,
            lr_decay=0.9,
            weight_decay=0.01,
        )
        
        # Collect parameter IDs from optimizer
        optimizer_param_ids = set()
        for group in param_groups:
            for param in group['params']:
                optimizer_param_ids.add(id(param))
        
        # Collect unique parameter IDs from model
        # With tie_weights=True, embedding.weight and lm_head.weight are same
        model_param_ids = set(id(p) for p in model.parameters())
        
        # Every unique parameter should be in optimizer
        missing = model_param_ids - optimizer_param_ids
        
        if missing:
            # Identify missing parameters
            id_to_name = {id(p): n for n, p in model.named_parameters()}
            missing_names = [id_to_name.get(pid, f"unknown-{pid}") for pid in missing]
            raise AssertionError(f"Parameters missing from optimizer: {missing_names}")
        
        print(f"✓ All {len(model_param_ids)} unique model parameters are in optimizer")
    
    def test_optimizer_param_count_with_and_without_tying(self):
        """
        Verify that with weight tying, optimizer has one fewer parameter
        than without (since embedding.weight and lm_head.weight are same).
        """
        model_tied = create_test_model(tie_weights=True)
        model_untied = create_test_model(tie_weights=False)
        
        # Count unique parameters with tying
        tied_params = set(model_tied.parameters())
        tied_count = len(tied_params)
        
        # Count unique parameters without tying  
        untied_params = set(model_untied.parameters())
        untied_count = len(untied_params)
        
        # With tying, should have exactly 1 fewer param (the shared weight)
        expected_difference = 1
        actual_difference = untied_count - tied_count
        
        assert actual_difference == expected_difference, (
            f"Expected {expected_difference} fewer params with weight tying, "
            f"but difference is {actual_difference}. "
            f"Tied: {tied_count}, Untied: {untied_count}"
        )
        
        print(f"✓ Weight tying reduces param count by 1 ({untied_count} → {tied_count})")


class TestGradientFlow:
    """Tests for gradient flow through tied weights."""
    
    def test_gradient_accumulates_correctly(self):
        """
        Verify gradients from both embedding and LM head loss contribute
        to the tied weight's gradient.
        """
        model = create_test_model(tie_weights=True)
        model.train()
        
        # Create dummy input
        batch_size, seq_len = 2, 10
        vocab_size = 67
        
        inputs = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
        labels = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
        
        # Forward and backward
        logits, loss = model(inputs, labels)
        loss.backward()
        
        # The tied weight should have gradients
        tied_weight = model.encoder.embedding.weight
        
        assert tied_weight.grad is not None, (
            "Tied weight has no gradient after backward pass!"
        )
        
        # Gradient should be non-zero (model is learning)
        grad_norm = tied_weight.grad.norm().item()
        assert grad_norm > 0, (
            f"Tied weight gradient is zero! norm={grad_norm}"
        )
        
        # Both references should show the same gradient (they're the same object)
        assert model.lm_head.weight.grad is tied_weight.grad, (
            "lm_head.weight.grad is not the same as embedding.weight.grad!"
        )
        
        print(f"✓ Gradient flows correctly through tied weights (norm={grad_norm:.6f})")
    
    def test_single_optimizer_step_updates_both(self):
        """
        Verify a single optimizer step updates both embedding and lm_head
        (since they're the same parameter).
        """
        model = create_test_model(tie_weights=True)
        model.train()
        
        # Get initial weight values
        initial_embed = model.encoder.embedding.weight[0, 0].item()
        initial_lm = model.lm_head.weight[0, 0].item()
        
        assert initial_embed == initial_lm, "Initial values should be equal (same param)"
        
        # Create optimizer
        optimizer = AdamW(model.parameters(), lr=0.1)  # High LR for visible change
        
        # Forward + backward
        batch_size, seq_len = 2, 10
        vocab_size = 67
        inputs = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
        labels = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
        
        logits, loss = model(inputs, labels)
        loss.backward()
        optimizer.step()
        
        # Get updated values
        updated_embed = model.encoder.embedding.weight[0, 0].item()
        updated_lm = model.lm_head.weight[0, 0].item()
        
        # They should still be equal (same param)
        assert updated_embed == updated_lm, (
            f"After optimizer step, values diverged! "
            f"embedding: {updated_embed}, lm_head: {updated_lm}"
        )
        
        # Value should have changed
        assert updated_embed != initial_embed, (
            "Weight did not change after optimizer step!"
        )
        
        print(f"✓ Single optimizer step updates both references ({initial_embed:.4f} → {updated_embed:.4f})")


class TestWeightNormStability:
    """Tests for weight norm stability during training."""
    
    def test_weight_norm_does_not_explode(self):
        """
        Verify that tied embedding/LM weight norm doesn't explode
        over multiple training steps.
        
        This was a symptom of double-optimization.
        """
        model = create_test_model(tie_weights=True)
        model.train()
        
        # Use layer-wise LR decay optimizer (production setup)
        param_groups = get_layer_wise_lr_groups(
            model=model,
            base_lr=3e-4,
            lr_decay=0.9,
            weight_decay=0.01,
        )
        optimizer = AdamW(param_groups, betas=(0.9, 0.98), eps=1e-6)
        
        # Track weight norms
        norms = []
        initial_norm = model.encoder.embedding.weight.norm().item()
        norms.append(initial_norm)
        
        # Run multiple training steps
        batch_size, seq_len = 4, 20
        vocab_size = 67
        num_steps = 100
        
        for step in range(num_steps):
            optimizer.zero_grad()
            
            inputs = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
            labels = torch.randint(0, vocab_size - 3, (batch_size, seq_len))
            
            logits, loss = model(inputs, labels)
            loss.backward()
            
            # Gradient clipping (as in production)
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            
            current_norm = model.encoder.embedding.weight.norm().item()
            norms.append(current_norm)
        
        final_norm = norms[-1]
        max_norm = max(norms)
        
        # Weight norm shouldn't explode (>10x initial is suspicious)
        explosion_threshold = initial_norm * 10
        assert max_norm < explosion_threshold, (
            f"Weight norm exploded! Initial: {initial_norm:.2f}, "
            f"Max: {max_norm:.2f}, Final: {final_norm:.2f}. "
            "This may indicate double-optimization of tied weights."
        )
        
        # Also check it doesn't collapse to near-zero
        collapse_threshold = initial_norm * 0.01
        assert final_norm > collapse_threshold, (
            f"Weight norm collapsed! Initial: {initial_norm:.2f}, "
            f"Final: {final_norm:.2f}"
        )
        
        print(f"✓ Weight norm stable over {num_steps} steps: "
              f"{initial_norm:.2f} → {final_norm:.2f} (max: {max_norm:.2f})")


def run_all_tests():
    """Run all weight tying tests and print summary."""
    print("=" * 70)
    print("WEIGHT TYING VERIFICATION TESTS")
    print("=" * 70)
    print()
    
    # Test classes to run
    test_classes = [
        TestWeightTying,
        TestOptimizerParameterGroups,
        TestGradientFlow,
        TestWeightNormStability,
    ]
    
    passed = 0
    failed = 0
    
    for test_class in test_classes:
        print(f"\n{test_class.__name__}")
        print("-" * 50)
        
        instance = test_class()
        for method_name in dir(instance):
            if method_name.startswith("test_"):
                try:
                    getattr(instance, method_name)()
                    passed += 1
                except AssertionError as e:
                    print(f"✗ {method_name}: FAILED")
                    print(f"  {e}")
                    failed += 1
                except Exception as e:
                    print(f"✗ {method_name}: ERROR")
                    print(f"  {type(e).__name__}: {e}")
                    failed += 1
    
    print("\n" + "=" * 70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)
    
    return failed == 0


if __name__ == "__main__":
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)

