"""Verification tests for the simplified GeneWhispererStage1."""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Dict, Optional
from unittest import SkipTest

import torch

# Add training directory to path for local imports.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "training"))

from model import GeneWhispererStage1, GeneWhispererStage1Legacy


def _count_parameters(model: torch.nn.Module) -> int:
    return sum(param.numel() for param in model.parameters())


def _small_model_kwargs() -> Dict[str, object]:
    return {
        "vocab_size": 256,
        "kmer": 6,
        "embedding_dim": 64,
        "num_layers": 2,
        "num_heads": 4,
        "ff_dim": 128,
        "dropout": 0.0,
        "pad_token_id": 0,
        "engineered_dim": 288,
        "use_engineered_features": True,
        "use_attention_pool": True,
        "engineered_mlp_hidden": 64,
        "engineered_mlp_output": 32,
        "engineered_mlp_input_dropout": 0.0,
        "engineered_mlp_pre_norm": True,
        "engineered_mlp_use_gated": True,
        "engineered_mlp_use_residual": True,
        "drop_path_rate": 0.0,
    }


def _artifacts_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "artifacts"


def _select_mlm_checkpoint() -> Optional[Path]:
    artifacts_dir = _artifacts_dir()
    for kmer in (6, 5, 4, 3):
        ckpt_path = artifacts_dir / f"mlm_encoder_k{kmer}.pt"
        if ckpt_path.exists():
            return ckpt_path
    return None


def _load_checkpoint_state(ckpt_path: Path) -> Dict[str, torch.Tensor]:
    checkpoint = torch.load(ckpt_path, map_location="cpu")
    if isinstance(checkpoint, dict):
        if "model_state_dict" in checkpoint:
            return checkpoint["model_state_dict"]
        if "state_dict" in checkpoint:
            return checkpoint["state_dict"]
    return checkpoint


def _infer_layer_count(state_dict: Dict[str, torch.Tensor]) -> int:
    layer_ids = set()
    for key in state_dict:
        if not key.startswith("layers."):
            continue
        parts = key.split(".")
        if len(parts) > 1 and parts[1].isdigit():
            layer_ids.add(int(parts[1]))
    if not layer_ids:
        return 0
    return max(layer_ids) + 1


def _infer_num_heads(state_dict: Dict[str, torch.Tensor], embedding_dim: int) -> int:
    rel_key = "relative_position_bias.relative_attention_bias.weight"
    if rel_key in state_dict:
        return int(state_dict[rel_key].shape[1])
    for candidate in (12, 8, 6, 4, 3, 2, 1):
        if embedding_dim % candidate == 0:
            return candidate
    return 1


def _infer_ff_dim(state_dict: Dict[str, torch.Tensor], embedding_dim: int) -> int:
    standard_key = "layers.0.ffn.0.weight"
    if standard_key in state_dict:
        return int(state_dict[standard_key].shape[0])
    glu_key = "layers.0.ffn.gate_proj.weight"
    if glu_key in state_dict:
        hidden_dim = int(state_dict[glu_key].shape[0])
        base = int(hidden_dim * 3 / 2)
        for ff_dim in range(max(1, base - 16), base + 17):
            projected = int(ff_dim * 2 / 3)
            projected = ((projected + 7) // 8) * 8
            if projected == hidden_dim:
                return ff_dim
        return base
    return embedding_dim * 4


def _infer_kmer_from_path(ckpt_path: Path) -> int:
    for token in ckpt_path.stem.split("_"):
        if token.startswith("k") and token[1:].isdigit():
            return int(token[1:])
    return 3


def test_model_instantiation() -> None:
    """Validate instantiation paths and parameter reduction."""
    default_model = GeneWhispererStage1()
    with_engineered = GeneWhispererStage1(engineered_dim=288, use_engineered_features=True)
    without_engineered = GeneWhispererStage1(engineered_dim=0, use_engineered_features=False)

    assert with_engineered.engineered_mlp is not None, "Engineered MLP should be initialized"
    assert without_engineered.engineered_mlp is None, "Engineered MLP should be disabled"

    small_model = GeneWhispererStage1(**_small_model_kwargs())
    simplified_params = _count_parameters(small_model)
    legacy_kwargs = dict(_small_model_kwargs())
    legacy_kwargs.update(
        {
            "use_tcn": True,
            "tcn_hidden": 64,
            "tcn_levels": 2,
            "tcn_kernel": 3,
            "multiscale_channels": 16,
            "multiscale_kernels": (3, 5),
            "post_cnn_transformer_layers": 1,
            "fusion_hidden": 64,
        }
    )
    stage1_params = _count_parameters(GeneWhispererStage1Legacy(**legacy_kwargs))
    assert simplified_params < stage1_params, "Simplified model should have fewer parameters"

    ratio = simplified_params / stage1_params
    assert ratio <= 0.85, f"Expected fewer params; ratio={ratio:.2f}"


def test_forward_pass() -> None:
    """Check forward outputs for shape and numerical stability."""
    torch.manual_seed(0)
    model = GeneWhispererStage1(**_small_model_kwargs())
    model.eval()

    tokens = torch.randint(1, model.embedding.num_embeddings, (4, 76))
    engineered = torch.randn(4, 288)

    with torch.no_grad():
        output = model(tokens, engineered)

    assert output.shape == (4, 1), f"Expected (4, 1), got {output.shape}"
    assert torch.isfinite(output).all(), "Output contains NaN/Inf values"


def test_gradient_flow() -> None:
    """Ensure gradients exist and are finite for all parameters."""
    torch.manual_seed(0)
    model = GeneWhispererStage1(**_small_model_kwargs())
    model.train()

    tokens = torch.randint(1, model.embedding.num_embeddings, (4, 76))
    engineered = torch.randn(4, 288)
    logits = model(tokens, engineered)
    labels = torch.randn(4, 1)
    loss = torch.nn.functional.mse_loss(logits, labels)
    loss.backward()

    for name, param in model.named_parameters():
        assert param.grad is not None, f"Missing gradient for {name}"
        assert torch.isfinite(param.grad).all(), f"Non-finite gradient for {name}"


def test_mlm_weight_loading() -> None:
    """Verify that MLM checkpoint loading updates encoder weights."""
    ckpt_path = _select_mlm_checkpoint()
    if ckpt_path is None:
        raise SkipTest("No MLM checkpoint found; skipping weight-load verification")

    state_dict = _load_checkpoint_state(ckpt_path)
    embedding_weight = state_dict.get("embedding.weight")
    if embedding_weight is None:
        raise SkipTest("Checkpoint missing embedding weights; skipping")

    vocab_size, embedding_dim = embedding_weight.shape
    num_layers = _infer_layer_count(state_dict)
    if num_layers < 1:
        raise SkipTest("Checkpoint missing transformer layers; skipping")
    num_heads = _infer_num_heads(state_dict, embedding_dim)
    ff_dim = _infer_ff_dim(state_dict, embedding_dim)
    use_glu_ffn = "layers.0.ffn.gate_proj.weight" in state_dict
    use_relative_position_bias = (
        "relative_position_bias.relative_attention_bias.weight" in state_dict
    )
    rel_buckets = None
    if use_relative_position_bias:
        rel_key = "relative_position_bias.relative_attention_bias.weight"
        rel_buckets = int(state_dict[rel_key].shape[0])

    model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=_infer_kmer_from_path(ckpt_path),
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.0,
        pad_token_id=0,
        engineered_dim=0,
        use_engineered_features=False,
        use_attention_pool=False,
        drop_path_rate=0.0,
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=rel_buckets or 32,
        relative_position_max_distance=128,
        use_glu_ffn=use_glu_ffn,
    )

    before_embed = model.embedding.weight.detach().clone()
    before_qproj = model.encoder.layers[0].q_proj.weight.detach().clone()

    model.load_pretrained_weights(state_dict, strict=False, transfer_mode="embed_only")

    after_embed = model.embedding.weight.detach()
    after_qproj = model.encoder.layers[0].q_proj.weight.detach()

    assert not torch.allclose(before_embed, after_embed), "Embedding weights did not change"
    assert not torch.allclose(before_qproj, after_qproj), "Transformer weights did not change"


def test_save_load_and_compile() -> None:
    """Confirm save/load round-trip and torch.compile compatibility."""
    model = GeneWhispererStage1(**_small_model_kwargs())
    model.eval()

    tokens = torch.randint(1, model.embedding.num_embeddings, (2, 16))
    engineered = torch.randn(2, 288)

    with torch.no_grad():
        baseline = model(tokens, engineered)

    with tempfile.TemporaryDirectory() as tmpdir:
        ckpt_path = Path(tmpdir) / "simplified_model.pt"
        torch.save(model.state_dict(), ckpt_path)

        reloaded = GeneWhispererStage1(**_small_model_kwargs())
        reloaded.load_state_dict(torch.load(ckpt_path, map_location="cpu"))
        reloaded.eval()

        with torch.no_grad():
            reloaded_out = reloaded(tokens, engineered)

    assert torch.allclose(baseline, reloaded_out, atol=1e-6), "Reloaded output mismatch"

    if hasattr(torch, "compile"):
        compiled = torch.compile(reloaded, backend="eager")
        with torch.no_grad():
            compiled_out = compiled(tokens, engineered)
        assert compiled_out.shape == (2, 1), "Compiled output has wrong shape"
        assert torch.isfinite(compiled_out).all(), "Compiled output contains NaN/Inf"
        assert torch.allclose(compiled_out, reloaded_out, atol=1e-5), (
            "Compiled output differs from eager output"
        )


def _run_test(name: str, fn) -> bool:
    try:
        fn()
    except SkipTest as exc:
        print(f"PASS: {name} (skipped: {exc})")
        return True
    except Exception as exc:  # noqa: BLE001 - clear output for script execution
        print(f"FAIL: {name} ({type(exc).__name__}: {exc})")
        return False
    print(f"PASS: {name}")
    return True


if __name__ == "__main__":
    tests = [
        ("model instantiation", test_model_instantiation),
        ("forward pass", test_forward_pass),
        ("gradient flow", test_gradient_flow),
        ("MLM weight loading", test_mlm_weight_loading),
        ("save/load and torch.compile", test_save_load_and_compile),
    ]
    all_passed = True
    for test_name, test_fn in tests:
        all_passed &= _run_test(test_name, test_fn)

    if not all_passed:
        raise SystemExit(1)
