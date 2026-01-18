"""Unit tests for embed_plus_adapter transfer mode.

Tests ensure:
1. A fake MLM checkpoint can be created and loaded
2. embed_plus_adapter properly transfers embeddings AND transformer layers
3. Transferred weights match exactly between checkpoint and model
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
import torch

from model import DNAEncoder, GeneWhispererStage1Legacy, PostCNNTransformerAdapter


def create_fake_mlm_checkpoint(
    vocab_size: int = 67,
    embedding_dim: int = 256,
    num_layers: int = 8,
    num_heads: int = 8,
    ff_dim: int = 1024,
    pad_token_id: int = 66,
) -> dict:
    """
    Create a fake MLM encoder state dict with random weights.

    Returns the encoder and its state dict for verification.
    """
    encoder = DNAEncoder(
        vocab_size=vocab_size,
        kmer=6,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.1,
        pad_token_id=pad_token_id,
    )

    # Initialize with specific random values for reproducibility
    torch.manual_seed(42)
    for param in encoder.parameters():
        param.data.uniform_(-1.0, 1.0)

    return encoder.state_dict(), encoder


def create_stage1_model(
    vocab_size: int = 67,
    embedding_dim: int = 256,
    num_layers: int = 6,
    num_heads: int = 8,
    ff_dim: int = 1024,
    pad_token_id: int = 66,
    post_cnn_transformer_layers: int = 3,
) -> GeneWhispererStage1Legacy:
    """Create a Stage 1 model with matching dimensions."""
    return GeneWhispererStage1Legacy(
        vocab_size=vocab_size,
        kmer=6,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.15,
        pad_token_id=pad_token_id,
        engineered_dim=128,
        use_engineered_features=True,
        use_attention_pool=True,
        use_tcn=True,
        tcn_hidden=256,
        tcn_levels=4,
        tcn_kernel=3,
        multiscale_channels=64,
        multiscale_kernels=(3, 5, 7, 9, 15),
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        engineered_mlp_hidden=256,
        engineered_mlp_output=128,
        fusion_hidden=256,
    )


class TestEmbedPlusAdapterTransfer:
    """Tests for embed_plus_adapter weight transfer mode."""

    def test_embed_only_does_not_change_post_cnn_transformer(self):
        """
        Verify embed_only mode does NOT modify post_cnn_transformer weights.
        """
        # Create checkpoint and model
        checkpoint_state, _ = create_fake_mlm_checkpoint()
        model = create_stage1_model()

        # Save post_cnn_transformer state before loading
        before_state = {
            name: param.clone()
            for name, param in model.post_cnn_transformer.named_parameters()
        }

        # Save checkpoint to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            # Load with embed_only mode
            model.load_pretrained_weights(
                checkpoint_path,
                strict=False,
                transfer_mode="embed_only"
            )

            # Verify post_cnn_transformer weights are unchanged
            for name, param in model.post_cnn_transformer.named_parameters():
                before_param = before_state[name]
                assert torch.allclose(param, before_param, atol=1e-7), (
                    f"embed_only changed post_cnn_transformer weight: {name}"
                )

            print("embed_only: post_cnn_transformer weights unchanged")

        finally:
            checkpoint_path.unlink()

    def test_embed_plus_adapter_changes_post_cnn_transformer(self):
        """
        Verify embed_plus_adapter mode DOES modify post_cnn_transformer weights.
        """
        # Create checkpoint and model
        checkpoint_state, _ = create_fake_mlm_checkpoint()
        model = create_stage1_model()

        # Save post_cnn_transformer state before loading
        before_state = {
            name: param.clone()
            for name, param in model.post_cnn_transformer.named_parameters()
        }

        # Save checkpoint to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            # Load with embed_plus_adapter mode
            model.load_pretrained_weights(
                checkpoint_path,
                strict=False,
                transfer_mode="embed_plus_adapter"
            )

            # Verify at least some post_cnn_transformer weights changed
            changed_count = 0
            for name, param in model.post_cnn_transformer.named_parameters():
                before_param = before_state[name]
                if not torch.allclose(param, before_param, atol=1e-7):
                    changed_count += 1

            assert changed_count > 0, (
                "embed_plus_adapter did NOT change any post_cnn_transformer weights!"
            )

            print(f"embed_plus_adapter: {changed_count} post_cnn_transformer params changed")

        finally:
            checkpoint_path.unlink()

    def test_embedding_weights_match_after_transfer(self):
        """
        Verify embedding weights match exactly after embed_plus_adapter transfer.
        """
        # Create checkpoint and model
        checkpoint_state, mlm_encoder = create_fake_mlm_checkpoint()
        model = create_stage1_model()

        # Save checkpoint to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            # Load with embed_plus_adapter mode
            model.load_pretrained_weights(
                checkpoint_path,
                strict=False,
                transfer_mode="embed_plus_adapter"
            )

            # Verify embedding weights match
            checkpoint_embedding = checkpoint_state["embedding.weight"]
            model_embedding = model.embedding.weight

            assert torch.allclose(model_embedding, checkpoint_embedding, atol=1e-7), (
                "Embedding weights do not match after transfer!"
            )

            print("Embedding weights match exactly")

        finally:
            checkpoint_path.unlink()

    def test_transformer_layer_weights_match_after_transfer(self):
        """
        CRITICAL TEST: Verify transformer layer weights from the checkpoint
        are correctly transferred to post_cnn_transformer layers.

        The transfer works as follows:
        1. Checkpoint is loaded into _full_encoder (which may have fewer layers)
        2. Top N layers from _full_encoder are copied to post_cnn_transformer

        For a checkpoint with 8 layers and Stage 1 with 6 encoder layers:
        - Checkpoint layers 0-5 are loaded into _full_encoder layers 0-5
        - _full_encoder layers 3-5 (top 3) are copied to adapter layers 0-2
        - So adapter gets checkpoint layers 3, 4, 5 (not 5, 6, 7)
        """
        # Use matching layer count for simpler verification
        num_mlm_layers = 6
        num_stage1_encoder_layers = 6
        num_post_cnn_layers = 3

        checkpoint_state, mlm_encoder = create_fake_mlm_checkpoint(num_layers=num_mlm_layers)
        model = create_stage1_model(
            num_layers=num_stage1_encoder_layers,
            post_cnn_transformer_layers=num_post_cnn_layers
        )

        # Save checkpoint to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            # Load with embed_plus_adapter mode
            model.load_pretrained_weights(
                checkpoint_path,
                strict=False,
                transfer_mode="embed_plus_adapter"
            )

            # Verify transformer layer weights match
            # _full_encoder has 6 layers, adapter wants 3
            # load_pretrained_layers copies from TOP of _full_encoder:
            # - src_idx = len(encoder_layers) - num_layers + i
            # - For i=0: src_idx = 6 - 3 + 0 = 3 -> adapter[0] gets encoder[3]
            # - For i=1: src_idx = 6 - 3 + 1 = 4 -> adapter[1] gets encoder[4]
            # - For i=2: src_idx = 6 - 3 + 2 = 5 -> adapter[2] gets encoder[5]
            #
            # Since checkpoint has 6 layers and _full_encoder has 6 layers,
            # all checkpoint layers 0-5 are loaded into _full_encoder.
            # So adapter[i] should match mlm_encoder[3+i]

            matched_layers = 0
            for adapter_idx in range(num_post_cnn_layers):
                # Adapter layer i gets _full_encoder layer (6 - 3 + i) = 3 + i
                encoder_src_idx = num_stage1_encoder_layers - num_post_cnn_layers + adapter_idx

                # Get adapter layer
                adapter_layer = model.post_cnn_transformer.layers[adapter_idx]

                # Get corresponding encoder layer from original checkpoint encoder
                encoder_layer = mlm_encoder.layers[encoder_src_idx]

                # Compare key weights (q_proj, k_proj, v_proj, out_proj)
                for weight_name in ["q_proj.weight", "k_proj.weight", "v_proj.weight", "out_proj.weight"]:
                    adapter_weight = dict(adapter_layer.named_parameters())[weight_name]
                    encoder_weight = dict(encoder_layer.named_parameters())[weight_name]

                    if torch.allclose(adapter_weight, encoder_weight, atol=1e-7):
                        pass  # Match!
                    else:
                        pytest.fail(
                            f"Layer weight mismatch! adapter[{adapter_idx}].{weight_name} "
                            f"!= encoder[{encoder_src_idx}].{weight_name}"
                        )

                matched_layers += 1

            assert matched_layers == num_post_cnn_layers, (
                f"Only {matched_layers}/{num_post_cnn_layers} layers matched!"
            )

            print(f"All {matched_layers} transformer layers match exactly")

        finally:
            checkpoint_path.unlink()

    def test_transfer_mode_none_skips_loading(self):
        """Verify transfer_mode='none' skips all weight loading."""
        checkpoint_state, _ = create_fake_mlm_checkpoint()
        model = create_stage1_model()

        # Save all weights before
        before_state = {
            name: param.clone()
            for name, param in model.named_parameters()
        }

        # Save checkpoint to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            # Load with none mode
            model.load_pretrained_weights(
                checkpoint_path,
                strict=False,
                transfer_mode="none"
            )

            # Verify NO weights changed
            for name, param in model.named_parameters():
                before_param = before_state[name]
                assert torch.allclose(param, before_param, atol=1e-7), (
                    f"transfer_mode='none' changed weight: {name}"
                )

            print("transfer_mode='none' correctly skips all weight loading")

        finally:
            checkpoint_path.unlink()

    def test_invalid_transfer_mode_raises(self):
        """Verify invalid transfer_mode raises ValueError."""
        checkpoint_state, _ = create_fake_mlm_checkpoint()
        model = create_stage1_model()

        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            torch.save(checkpoint_state, f.name)
            checkpoint_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Unsupported transfer_mode"):
                model.load_pretrained_weights(
                    checkpoint_path,
                    strict=False,
                    transfer_mode="invalid_mode"
                )

            print("Invalid transfer_mode correctly raises ValueError")

        finally:
            checkpoint_path.unlink()


def run_all_tests():
    """Run all embed_plus_adapter transfer tests and print summary."""
    print("=" * 70)
    print("EMBED_PLUS_ADAPTER TRANSFER VERIFICATION TESTS")
    print("=" * 70)
    print()

    test_class = TestEmbedPlusAdapterTransfer

    passed = 0
    failed = 0

    print(f"{test_class.__name__}")
    print("-" * 50)

    instance = test_class()
    for method_name in dir(instance):
        if method_name.startswith("test_"):
            try:
                getattr(instance, method_name)()
                passed += 1
            except AssertionError as e:
                print(f"  FAILED {method_name}")
                print(f"    {e}")
                failed += 1
            except Exception as e:
                print(f"  ERROR {method_name}")
                print(f"    {type(e).__name__}: {e}")
                failed += 1

    print("\n" + "=" * 70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)

    return failed == 0


if __name__ == "__main__":
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)
