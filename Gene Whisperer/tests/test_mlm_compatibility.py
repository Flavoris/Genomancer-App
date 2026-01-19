"""
Test MLM checkpoint compatibility with the simplified GeneWhispererStage1 model.

Verifies that existing MLM checkpoints can be loaded into the new simplified
architecture without retraining.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pytest
import torch
import torch.nn as nn

# Add training directory to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = SCRIPT_DIR.parent / "training"
sys.path.insert(0, str(TRAINING_DIR))

from model import DNAEncoder, GeneWhispererStage1


# Paths to artifacts
ARTIFACTS_DIR = SCRIPT_DIR.parent / "artifacts"
MLM_CHECKPOINT_K6 = ARTIFACTS_DIR / "mlm_encoder_k6.pt"
MLM_METADATA_K6 = ARTIFACTS_DIR / "mlm_encoder_k6.metadata.json"


def load_metadata(metadata_path: Path) -> dict[str, Any]:
    """Load checkpoint metadata JSON."""
    if not metadata_path.exists():
        return {}
    with metadata_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def extract_state_dict(checkpoint: dict) -> dict[str, torch.Tensor]:
    """Extract the encoder state dict from a checkpoint."""
    # Handle nested checkpoint formats
    if "model_state_dict" in checkpoint:
        state = checkpoint["model_state_dict"]
    elif "encoder_state_dict" in checkpoint:
        state = checkpoint["encoder_state_dict"]
    elif "state_dict" in checkpoint:
        state = checkpoint["state_dict"]
    else:
        # Assume checkpoint IS the state dict
        state = checkpoint

    # Filter to only tensor values (skip metadata)
    return {k: v for k, v in state.items() if isinstance(v, torch.Tensor)}


def detect_relative_position_bias(state_dict: dict[str, torch.Tensor]) -> bool:
    """Detect if checkpoint uses relative position bias from keys."""
    return any("relative_position_bias" in k for k in state_dict.keys())


def detect_num_layers(state_dict: dict[str, torch.Tensor]) -> int:
    """Detect number of transformer layers from checkpoint keys."""
    layer_indices = set()
    for key in state_dict.keys():
        if key.startswith("layers."):
            try:
                idx = int(key.split(".")[1])
                layer_indices.add(idx)
            except (ValueError, IndexError):
                pass
    return len(layer_indices) if layer_indices else 8  # Default to 8


def detect_embedding_dim(state_dict: dict[str, torch.Tensor]) -> int:
    """Detect embedding dimension from checkpoint tensors."""
    if "embedding.weight" in state_dict:
        return state_dict["embedding.weight"].shape[1]
    return 256  # Default


def detect_num_heads(state_dict: dict[str, torch.Tensor]) -> int:
    """Detect number of attention heads from checkpoint tensors."""
    # Check relative position bias (has shape [num_buckets, num_heads])
    rel_pos_key = "relative_position_bias.relative_attention_bias.weight"
    if rel_pos_key in state_dict:
        return state_dict[rel_pos_key].shape[1]
    # Check attention projection shapes
    for key in state_dict.keys():
        if "q_proj.weight" in key:
            # q_proj has shape [embed_dim, embed_dim], heads = embed_dim / head_dim
            # but we can't easily determine head_dim, so fallback
            pass
    return 8  # Default


def create_model_matching_checkpoint(
    metadata: dict[str, Any],
    state_dict: dict[str, torch.Tensor],
) -> GeneWhispererStage1:
    """Create a model with dimensions matching the checkpoint."""
    # Prefer to detect from checkpoint, fallback to metadata
    embedding_dim = detect_embedding_dim(state_dict) or metadata.get("embedding_dim", 256)
    num_layers = detect_num_layers(state_dict) or metadata.get("num_layers", 8)
    num_heads = detect_num_heads(state_dict) or metadata.get("heads", 8)
    kmer = metadata.get("kmer", 6)
    vocab_size = 4**kmer + 3

    # Detect relative position bias from checkpoint keys
    use_relative_position_bias = detect_relative_position_bias(state_dict)

    return GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=embedding_dim * 4,
        dropout=0.1,
        pad_token_id=0,
        use_engineered_features=False,
        use_attention_pool=True,
        use_relative_position_bias=use_relative_position_bias,
    )


@pytest.fixture
def mlm_checkpoint() -> dict[str, torch.Tensor] | None:
    """Load the MLM checkpoint if it exists."""
    if not MLM_CHECKPOINT_K6.exists():
        pytest.skip(f"MLM checkpoint not found: {MLM_CHECKPOINT_K6}")
    raw = torch.load(MLM_CHECKPOINT_K6, map_location="cpu", weights_only=False)
    return extract_state_dict(raw)


@pytest.fixture
def mlm_metadata() -> dict[str, Any]:
    """Load MLM checkpoint metadata."""
    return load_metadata(MLM_METADATA_K6)


@pytest.fixture
def model_from_metadata(
    mlm_metadata: dict[str, Any],
    mlm_checkpoint: dict[str, torch.Tensor],
) -> GeneWhispererStage1:
    """Create a model with dimensions matching the MLM checkpoint."""
    return create_model_matching_checkpoint(mlm_metadata, mlm_checkpoint)


class TestMLMCheckpointStructure:
    """Test 1: Verify MLM checkpoint structure and keys."""

    def test_checkpoint_is_dict(self, mlm_checkpoint: dict) -> None:
        """MLM checkpoint should be a state dict."""
        assert isinstance(mlm_checkpoint, dict)
        assert len(mlm_checkpoint) > 0

    def test_checkpoint_contains_expected_keys(self, mlm_checkpoint: dict) -> None:
        """MLM checkpoint should contain encoder component keys."""
        expected_prefixes = {
            "embedding.",
            "pos_embedding.",
            "embed_norm.",
            "layers.",
            "final_norm.",
        }

        found_prefixes = set()
        for key in mlm_checkpoint.keys():
            for prefix in expected_prefixes:
                if key.startswith(prefix):
                    found_prefixes.add(prefix)

        missing = expected_prefixes - found_prefixes
        assert not missing, f"Missing expected key prefixes: {missing}"

    def test_checkpoint_has_no_cnn_tcn_keys(self, mlm_checkpoint: dict) -> None:
        """MLM checkpoint should not contain CNN/TCN components."""
        forbidden_patterns = ["cnn", "tcn", "conv", "causal", "multiscale"]

        for key in mlm_checkpoint.keys():
            key_lower = key.lower()
            for pattern in forbidden_patterns:
                assert pattern not in key_lower, (
                    f"Found forbidden pattern '{pattern}' in key: {key}"
                )

    def test_list_all_checkpoint_keys(self, mlm_checkpoint: dict) -> None:
        """Print all keys in the checkpoint for inspection."""
        print("\n" + "=" * 60)
        print("MLM CHECKPOINT KEYS:")
        print("=" * 60)
        for key in sorted(mlm_checkpoint.keys()):
            shape = tuple(mlm_checkpoint[key].shape)
            print(f"  {key}: {shape}")
        print("=" * 60)
        print(f"Total keys: {len(mlm_checkpoint)}")


class TestModelWeightCompatibility:
    """Test 2-3: Verify weight shapes match between checkpoint and model."""

    def test_model_encoder_keys_match_checkpoint(
        self,
        mlm_checkpoint: dict,
        model_from_metadata: GeneWhispererStage1,
    ) -> None:
        """Model encoder keys should match checkpoint keys."""
        model_keys = set(model_from_metadata.encoder.state_dict().keys())
        checkpoint_keys = set(mlm_checkpoint.keys())

        # Keys that are in both
        matching = model_keys & checkpoint_keys
        # Keys only in model
        model_only = model_keys - checkpoint_keys
        # Keys only in checkpoint
        checkpoint_only = checkpoint_keys - model_keys

        print("\n" + "=" * 60)
        print("KEY COMPATIBILITY REPORT:")
        print("=" * 60)
        print(f"Matching keys: {len(matching)}")
        print(f"Model-only keys: {len(model_only)}")
        print(f"Checkpoint-only keys: {len(checkpoint_only)}")

        if model_only:
            print("\nKeys in model but not checkpoint:")
            for key in sorted(model_only):
                print(f"  - {key}")

        if checkpoint_only:
            print("\nKeys in checkpoint but not model:")
            for key in sorted(checkpoint_only):
                print(f"  - {key}")

        # At least 90% of checkpoint keys should match
        match_ratio = len(matching) / len(checkpoint_keys) if checkpoint_keys else 0
        assert match_ratio >= 0.9, f"Only {match_ratio*100:.1f}% keys match"

    def test_weight_shapes_match(
        self,
        mlm_checkpoint: dict,
        model_from_metadata: GeneWhispererStage1,
    ) -> None:
        """Weight shapes should match between checkpoint and model."""
        model_state = model_from_metadata.encoder.state_dict()

        shape_mismatches = []
        shape_matches = 0

        for key, ckpt_tensor in mlm_checkpoint.items():
            if key in model_state:
                model_shape = tuple(model_state[key].shape)
                ckpt_shape = tuple(ckpt_tensor.shape)

                if model_shape == ckpt_shape:
                    shape_matches += 1
                else:
                    shape_mismatches.append((key, model_shape, ckpt_shape))

        print("\n" + "=" * 60)
        print("SHAPE COMPATIBILITY REPORT:")
        print("=" * 60)
        print(f"Matching shapes: {shape_matches}")
        print(f"Mismatched shapes: {len(shape_mismatches)}")

        if shape_mismatches:
            print("\nShape mismatches:")
            for key, model_shape, ckpt_shape in shape_mismatches:
                print(f"  {key}: model={model_shape}, checkpoint={ckpt_shape}")

        assert not shape_mismatches, f"Found {len(shape_mismatches)} shape mismatches"


class TestWeightLoading:
    """Test 4: Load weights and verify successful transfer."""

    def test_load_pretrained_weights(
        self,
        model_from_metadata: GeneWhispererStage1,
    ) -> None:
        """Model should successfully load pretrained MLM weights."""
        if not MLM_CHECKPOINT_K6.exists():
            pytest.skip(f"MLM checkpoint not found: {MLM_CHECKPOINT_K6}")

        # Record original embedding values
        orig_embedding = model_from_metadata.encoder.embedding.weight.clone()

        # Load weights
        model_from_metadata.load_pretrained_weights(
            MLM_CHECKPOINT_K6,
            strict=False,
            transfer_mode="embed_only",
        )

        # Verify embedding changed
        new_embedding = model_from_metadata.encoder.embedding.weight
        assert not torch.allclose(orig_embedding, new_embedding), (
            "Embedding weights did not change after loading"
        )

    def test_count_loaded_weights(
        self,
        mlm_checkpoint: dict,
        model_from_metadata: GeneWhispererStage1,
    ) -> None:
        """Count how many weights loaded successfully."""
        model_state = model_from_metadata.encoder.state_dict()

        loaded = 0
        skipped = 0
        total_params_loaded = 0
        total_params_skipped = 0

        for key, ckpt_tensor in mlm_checkpoint.items():
            if key in model_state:
                model_shape = tuple(model_state[key].shape)
                ckpt_shape = tuple(ckpt_tensor.shape)

                if model_shape == ckpt_shape:
                    loaded += 1
                    total_params_loaded += ckpt_tensor.numel()
                else:
                    skipped += 1
                    total_params_skipped += ckpt_tensor.numel()
            else:
                skipped += 1
                total_params_skipped += ckpt_tensor.numel()

        print("\n" + "=" * 60)
        print("WEIGHT LOADING REPORT:")
        print("=" * 60)
        print(f"Weights loaded: {loaded}")
        print(f"Weights skipped: {skipped}")
        print(f"Parameters loaded: {total_params_loaded:,}")
        print(f"Parameters skipped: {total_params_skipped:,}")

        load_ratio = loaded / len(mlm_checkpoint) if mlm_checkpoint else 0
        print(f"Load ratio: {load_ratio*100:.1f}%")

        # Should load at least 95% of weights
        assert load_ratio >= 0.95, (
            f"Only {load_ratio*100:.1f}% weights loaded, expected >= 95%"
        )


class TestLoadedModelFunctionality:
    """Test 5: Verify loaded model works correctly."""

    def test_forward_pass_after_loading(
        self,
        model_from_metadata: GeneWhispererStage1,
        mlm_metadata: dict,
    ) -> None:
        """Model should produce valid outputs after loading weights."""
        if not MLM_CHECKPOINT_K6.exists():
            pytest.skip(f"MLM checkpoint not found: {MLM_CHECKPOINT_K6}")

        model_from_metadata.load_pretrained_weights(
            MLM_CHECKPOINT_K6,
            strict=False,
        )
        model_from_metadata.eval()

        # Create random input
        kmer = mlm_metadata.get("kmer", 6)
        vocab_size = 4**kmer + 3
        batch_size = 4
        seq_len = 64

        token_ids = torch.randint(3, vocab_size, (batch_size, seq_len))

        with torch.no_grad():
            output = model_from_metadata(token_ids)

        # Check output is valid
        assert output.shape == (batch_size, 1), f"Unexpected output shape: {output.shape}"
        assert not torch.isnan(output).any(), "Output contains NaN values"
        assert not torch.isinf(output).any(), "Output contains Inf values"

    def test_output_distribution_changes(
        self,
        mlm_metadata: dict,
        mlm_checkpoint: dict,
    ) -> None:
        """Output distribution should change after loading pretrained weights."""
        if not MLM_CHECKPOINT_K6.exists():
            pytest.skip(f"MLM checkpoint not found: {MLM_CHECKPOINT_K6}")

        kmer = mlm_metadata.get("kmer", 6)
        vocab_size = 4**kmer + 3

        # Create fresh model matching checkpoint architecture
        model = create_model_matching_checkpoint(mlm_metadata, mlm_checkpoint)
        model.eval()

        # Create fixed input for comparison
        torch.manual_seed(42)
        batch_size = 8
        seq_len = 64
        token_ids = torch.randint(3, vocab_size, (batch_size, seq_len))

        # Get output before loading
        with torch.no_grad():
            output_before = model(token_ids).clone()

        # Load pretrained weights
        model.load_pretrained_weights(MLM_CHECKPOINT_K6, strict=False)

        # Get output after loading
        with torch.no_grad():
            output_after = model(token_ids)

        # Outputs should differ
        assert not torch.allclose(output_before, output_after, atol=1e-3), (
            "Output distribution did not change after loading weights"
        )

        # Compute statistics
        print("\n" + "=" * 60)
        print("OUTPUT DISTRIBUTION COMPARISON:")
        print("=" * 60)
        print(f"Before loading - mean: {output_before.mean():.4f}, std: {output_before.std():.4f}")
        print(f"After loading  - mean: {output_after.mean():.4f}, std: {output_after.std():.4f}")
        print(f"Max absolute difference: {(output_after - output_before).abs().max():.4f}")


class TestCompatibilitySummary:
    """Test 6: Print final compatibility summary."""

    def test_print_compatibility_summary(
        self,
        mlm_checkpoint: dict,
        mlm_metadata: dict,
    ) -> None:
        """Print comprehensive compatibility summary."""
        if not MLM_CHECKPOINT_K6.exists():
            pytest.skip(f"MLM checkpoint not found: {MLM_CHECKPOINT_K6}")

        # Detect actual dimensions from checkpoint
        embedding_dim = detect_embedding_dim(mlm_checkpoint)
        num_layers = detect_num_layers(mlm_checkpoint)
        num_heads = detect_num_heads(mlm_checkpoint)
        kmer = mlm_metadata.get("kmer", 6)
        vocab_size = 4**kmer + 3

        # Create model matching checkpoint architecture
        model = create_model_matching_checkpoint(mlm_metadata, mlm_checkpoint)

        model_state = model.encoder.state_dict()

        # Count compatible weights
        compatible = 0
        incompatible = 0
        compatible_params = 0
        total_params = 0

        for key, ckpt_tensor in mlm_checkpoint.items():
            total_params += ckpt_tensor.numel()
            if key in model_state:
                if tuple(model_state[key].shape) == tuple(ckpt_tensor.shape):
                    compatible += 1
                    compatible_params += ckpt_tensor.numel()
                else:
                    incompatible += 1
            else:
                incompatible += 1

        # Print summary
        print("\n")
        print("=" * 70)
        print("MLM CHECKPOINT COMPATIBILITY SUMMARY")
        print("=" * 70)
        print(f"Checkpoint: {MLM_CHECKPOINT_K6.name}")
        print(f"Model: GeneWhispererStage1 (simplified)")
        print("-" * 70)
        print(f"Checkpoint dimensions:")
        print(f"  - embedding_dim: {embedding_dim}")
        print(f"  - num_layers: {num_layers}")
        print(f"  - num_heads: {num_heads}")
        print(f"  - kmer: {kmer}")
        print(f"  - vocab_size: {vocab_size}")
        print("-" * 70)
        print(f"Weight compatibility:")
        print(f"  - Compatible weights: {compatible}/{len(mlm_checkpoint)}")
        print(f"  - Incompatible weights: {incompatible}")
        print(f"  - Parameters transferable: {compatible_params:,}/{total_params:,}")
        print("-" * 70)

        transfer_ratio = compatible_params / total_params if total_params > 0 else 0
        weight_ratio = compatible / len(mlm_checkpoint) if mlm_checkpoint else 0

        if weight_ratio >= 0.95:
            print("RESULT: MLM checkpoint compatible - no retraining needed")
            print("=" * 70)
        else:
            print("RESULT: Consider retraining MLM for optimal performance")
            print(f"        Only {weight_ratio*100:.1f}% of weights transfer successfully")
            print("=" * 70)

        # Assert high compatibility
        assert weight_ratio >= 0.95, (
            f"Compatibility too low: {weight_ratio*100:.1f}% < 95%"
        )


def run_compatibility_check() -> None:
    """Run the compatibility check as a standalone script."""
    import sys

    print("\n" + "=" * 70)
    print("RUNNING MLM CHECKPOINT COMPATIBILITY CHECK")
    print("=" * 70)

    if not MLM_CHECKPOINT_K6.exists():
        print(f"ERROR: MLM checkpoint not found: {MLM_CHECKPOINT_K6}")
        sys.exit(1)

    # Load checkpoint and metadata
    raw_checkpoint = torch.load(MLM_CHECKPOINT_K6, map_location="cpu", weights_only=False)
    checkpoint = extract_state_dict(raw_checkpoint)
    metadata = load_metadata(MLM_METADATA_K6)

    # Extract dimensions
    embedding_dim = metadata.get("embedding_dim", 256)
    num_layers = metadata.get("num_layers", 8)
    num_heads = metadata.get("heads", 8)
    kmer = metadata.get("kmer", 6)
    vocab_size = 4**kmer + 3

    print(f"\nCheckpoint: {MLM_CHECKPOINT_K6.name}")
    print(f"Dimensions: embed={embedding_dim}, layers={num_layers}, heads={num_heads}, k={kmer}")

    # Create model matching checkpoint architecture
    model = create_model_matching_checkpoint(metadata, checkpoint)

    # Analyze compatibility
    model_state = model.encoder.state_dict()

    print("\n--- Checkpoint Keys ---")
    for key in sorted(checkpoint.keys()):
        shape = tuple(checkpoint[key].shape)
        in_model = key in model_state
        shape_match = False
        if in_model:
            shape_match = tuple(model_state[key].shape) == shape
        status = "OK" if in_model and shape_match else "SKIP" if in_model else "MISSING"
        print(f"  [{status}] {key}: {shape}")

    # Count stats
    compatible = sum(
        1 for k, v in checkpoint.items()
        if k in model_state and tuple(model_state[k].shape) == tuple(v.shape)
    )
    compatible_params = sum(
        v.numel() for k, v in checkpoint.items()
        if k in model_state and tuple(model_state[k].shape) == tuple(v.shape)
    )
    total_params = sum(v.numel() for v in checkpoint.values())

    weight_ratio = compatible / len(checkpoint)
    param_ratio = compatible_params / total_params

    print(f"\n--- Summary ---")
    print(f"Compatible weights: {compatible}/{len(checkpoint)} ({weight_ratio*100:.1f}%)")
    print(f"Transferable params: {compatible_params:,}/{total_params:,} ({param_ratio*100:.1f}%)")

    # Test loading
    print("\n--- Testing Weight Loading ---")
    orig_embed = model.encoder.embedding.weight.clone()
    model.load_pretrained_weights(MLM_CHECKPOINT_K6, strict=False)
    new_embed = model.encoder.embedding.weight
    embed_changed = not torch.allclose(orig_embed, new_embed)
    print(f"Embedding weights changed: {embed_changed}")

    # Test forward pass
    print("\n--- Testing Forward Pass ---")
    model.eval()
    test_input = torch.randint(3, vocab_size, (4, 64))
    with torch.no_grad():
        output = model(test_input)
    print(f"Output shape: {output.shape}")
    print(f"Output valid: {not torch.isnan(output).any() and not torch.isinf(output).any()}")

    # Final verdict
    print("\n" + "=" * 70)
    if weight_ratio >= 0.95:
        print("MLM checkpoint compatible - no retraining needed")
    else:
        print("Consider retraining MLM for optimal performance")
    print("=" * 70)


if __name__ == "__main__":
    run_compatibility_check()
