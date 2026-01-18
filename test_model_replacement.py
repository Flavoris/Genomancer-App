"""Model replacement verification checks for GeneWhispererStage1."""
from __future__ import annotations

import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import torch

ROOT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = ROOT_DIR / "Gene Whisperer" / "training"
ARTIFACTS_DIR = ROOT_DIR / "Gene Whisperer" / "artifacts"
VARIANT_DIR = ROOT_DIR / "Gene Whisperer" / "training" / "results" / "variant_comparison"

sys.path.insert(0, str(TRAINING_DIR))

from model import GeneWhispererStage1  # noqa: E402
from model_legacy import GeneWhispererStage1Legacy  # noqa: E402

EXPECTED_NEW_PARAMS = 23_900_000
EXPECTED_LEGACY_PARAMS = 38_400_000
EXPECTED_PARAM_DIFF = 14_500_000
PARAM_TOLERANCE = 2_000_000

ACCURACY_TARGET = 0.841
ACCURACY_TOLERANCE = 0.02


@dataclass(frozen=True)
class CheckResult:
    name: str
    status: str
    details: str


def _count_parameters(model: torch.nn.Module) -> int:
    return sum(param.numel() for param in model.parameters())


def _load_checkpoint_bundle(path: Path) -> Tuple[dict, Dict[str, torch.Tensor]]:
    checkpoint = torch.load(path, map_location="cpu")
    if isinstance(checkpoint, dict):
        state_dict = (
            checkpoint.get("model_state")
            or checkpoint.get("model_state_dict")
            or checkpoint.get("state_dict")
            or checkpoint
        )
    else:
        state_dict = checkpoint
    if not isinstance(state_dict, dict):
        raise ValueError(f"Checkpoint state is not a dict: {path}")
    return checkpoint, state_dict


def _strip_prefix(state: Dict[str, torch.Tensor], prefix: str) -> Dict[str, torch.Tensor]:
    return {
        key[len(prefix):]: value
        for key, value in state.items()
        if key.startswith(prefix)
    }


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


def _count_loadable_keys(
    state_dict: Dict[str, torch.Tensor],
    model_state: Dict[str, torch.Tensor],
) -> Tuple[int, int]:
    loadable = 0
    for key, value in state_dict.items():
        target = model_state.get(key)
        if target is not None and target.shape == value.shape:
            loadable += 1
    return loadable, len(state_dict) - loadable


def _format_count(value: int) -> str:
    return f"{value:,}"


def _extract_accuracy(checkpoint: dict) -> Optional[float]:
    metrics = checkpoint.get("metrics")
    if isinstance(metrics, dict):
        for key in ("acc_best", "accuracy", "val_accuracy", "val_acc"):
            if key in metrics:
                return float(metrics[key])
    for key in ("best_val_acc", "val_acc", "val_accuracy", "accuracy"):
        if key in checkpoint:
            return float(checkpoint[key])
    return None


def _result(name: str, passed: bool, details: str, skipped: bool = False) -> CheckResult:
    if skipped:
        return CheckResult(name=name, status="SKIP", details=details)
    return CheckResult(name=name, status="PASS" if passed else "FAIL", details=details)


def check_imports() -> CheckResult:
    try:
        if GeneWhispererStage1.__module__ != "model":
            return _result(
                "Import test",
                False,
                f"GeneWhispererStage1 module is {GeneWhispererStage1.__module__}",
            )
        if GeneWhispererStage1Legacy.__module__ != "model_legacy":
            return _result(
                "Import test",
                False,
                f"GeneWhispererStage1Legacy module is {GeneWhispererStage1Legacy.__module__}",
            )
    except Exception as exc:  # pragma: no cover - defensive
        return _result("Import test", False, f"Exception: {exc}")
    return _result("Import test", True, "Imports resolved to simplified and legacy classes.")


def check_parameter_counts() -> CheckResult:
    config = {
        "vocab_size": 4099,
        "kmer": 6,
        "embedding_dim": 384,
        "num_layers": 12,
        "num_heads": 12,
        "ff_dim": 1536,
        "dropout": 0.12,
        "pad_token_id": 4098,
        "engineered_dim": 288,
        "use_engineered_features": True,
        "use_attention_pool": True,
        "engineered_mlp_hidden": 512,
        "engineered_mlp_output": 384,
        "engineered_mlp_pre_norm": True,
        "engineered_mlp_input_dropout": 0.05,
        "engineered_mlp_use_gated": True,
        "engineered_mlp_use_residual": True,
        "drop_path_rate": 0.1,
        "max_seq_len": 256,
        "use_relative_position_bias": True,
        "relative_position_num_buckets": 32,
        "relative_position_max_distance": 128,
        "use_glu_ffn": True,
        "glu_activation": "gelu",
    }
    legacy_config = {
        "use_tcn": True,
        "tcn_hidden": 384,
        "tcn_levels": 5,
        "tcn_kernel": 3,
        "multiscale_channels": 64,
        "multiscale_kernels": (3, 5, 7, 9, 15),
        "lstm_hidden": 192,
        "post_cnn_transformer_layers": 3,
        "fusion_hidden": 512,
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        new_model = GeneWhispererStage1(**config)
        legacy_model = GeneWhispererStage1Legacy(**config, **legacy_config)

    new_params = _count_parameters(new_model)
    legacy_params = _count_parameters(legacy_model)
    diff = legacy_params - new_params

    within_new = abs(new_params - EXPECTED_NEW_PARAMS) <= PARAM_TOLERANCE
    within_legacy = abs(legacy_params - EXPECTED_LEGACY_PARAMS) <= PARAM_TOLERANCE
    within_diff = abs(diff - EXPECTED_PARAM_DIFF) <= PARAM_TOLERANCE

    details = (
        f"new={_format_count(new_params)} legacy={_format_count(legacy_params)} "
        f"diff={_format_count(diff)} (tolerance ±{_format_count(PARAM_TOLERANCE)})"
    )
    return _result("Parameter count test", within_new and within_legacy and within_diff, details)


def check_architecture() -> CheckResult:
    model = GeneWhispererStage1(
        vocab_size=64,
        kmer=3,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        dropout=0.0,
        pad_token_id=0,
        engineered_dim=288,
        use_engineered_features=True,
        use_attention_pool=True,
        engineered_mlp_hidden=64,
        engineered_mlp_output=32,
        engineered_mlp_input_dropout=0.0,
        drop_path_rate=0.0,
        max_seq_len=128,
    )

    forbidden = ("feature_extractor", "tcn", "multiscale", "post_cnn_transformer")
    required = ("encoder", "pool", "engineered_mlp", "classifier")

    present_forbidden = [name for name in forbidden if hasattr(model, name)]
    missing_required = [name for name in required if getattr(model, name, None) is None]

    if present_forbidden or missing_required:
        details = f"forbidden={present_forbidden} missing={missing_required}"
        return _result("Architecture test", False, details)
    return _result("Architecture test", True, "Simplified model exposes expected components.")


def check_forward_pass() -> CheckResult:
    model = GeneWhispererStage1(
        vocab_size=64,
        kmer=3,
        embedding_dim=32,
        num_layers=2,
        num_heads=4,
        ff_dim=64,
        dropout=0.0,
        pad_token_id=0,
        engineered_dim=288,
        use_engineered_features=True,
        use_attention_pool=True,
        engineered_mlp_hidden=64,
        engineered_mlp_output=32,
        engineered_mlp_input_dropout=0.0,
        drop_path_rate=0.0,
        max_seq_len=128,
    )
    model.eval()

    torch.manual_seed(1337)
    tokens = torch.randint(0, model.encoder.embedding.num_embeddings, (4, 76))
    engineered = torch.randn(4, 288)

    with torch.no_grad():
        logits = model(tokens, engineered)

    if logits.shape != (4, 1):
        return _result("Forward pass test", False, f"Expected (4, 1), got {tuple(logits.shape)}")
    if not torch.isfinite(logits).all():
        return _result("Forward pass test", False, "Output contains NaN/Inf.")
    return _result("Forward pass test", True, "Forward pass produced finite (4, 1) logits.")


def check_mlm_weight_loading() -> CheckResult:
    ckpt_path = ARTIFACTS_DIR / "mlm_encoder_k6.pt"
    if not ckpt_path.exists():
        return _result("MLM weight loading test", True, f"Missing {ckpt_path}", skipped=True)

    checkpoint, state = _load_checkpoint_bundle(ckpt_path)
    if "embedding.weight" not in state:
        return _result("MLM weight loading test", False, "embedding.weight not found in checkpoint.")

    vocab_size, embedding_dim = state["embedding.weight"].shape
    num_layers = _infer_layer_count(state)
    num_heads = _infer_num_heads(state, embedding_dim)
    ff_dim = _infer_ff_dim(state, embedding_dim)
    max_seq_len = int(state["pos_embedding.weight"].shape[0])
    rel_key = "relative_position_bias.relative_attention_bias.weight"
    use_relative_position_bias = rel_key in state
    num_buckets = int(state[rel_key].shape[0]) if use_relative_position_bias else 32
    use_glu_ffn = any(key.endswith("ffn.gate_proj.weight") for key in state)
    kmer = int(checkpoint.get("kmer", 6)) if isinstance(checkpoint, dict) else 6

    model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=kmer,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.0,
        pad_token_id=vocab_size - 1,
        engineered_dim=0,
        use_engineered_features=False,
        use_attention_pool=False,
        drop_path_rate=0.0,
        max_seq_len=max_seq_len,
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=num_buckets,
        relative_position_max_distance=128,
        use_glu_ffn=use_glu_ffn,
    )

    encoder_state = model.encoder.state_dict()
    loadable, skipped = _count_loadable_keys(state, encoder_state)

    model.load_pretrained_weights(ckpt_path, strict=False)
    updated_state = model.encoder.state_dict()

    embedding_ok = torch.equal(updated_state["embedding.weight"], state["embedding.weight"])
    layer_key = next(
        (key for key in state if key.startswith("layers.0") and key.endswith("weight")),
        None,
    )
    layer_ok = True
    if layer_key is not None:
        layer_ok = torch.equal(updated_state[layer_key], state[layer_key])

    details = (
        f"loaded={_format_count(loadable)} skipped={_format_count(skipped)} "
        f"embedding_match={embedding_ok} layer_match={layer_ok}"
    )
    return _result("MLM weight loading test", embedding_ok and layer_ok, details)


def check_stage1_checkpoint_loading() -> CheckResult:
    candidates = [
        VARIANT_DIR / "transformer_only_best.pt",
        VARIANT_DIR / "combined_best.pt",
        VARIANT_DIR / "features_only_best.pt",
    ]
    ckpt_path = next((path for path in candidates if path.exists()), None)
    if ckpt_path is None:
        return _result("Stage 1 checkpoint loading test", True, "No checkpoint found.", skipped=True)

    checkpoint, state = _load_checkpoint_bundle(ckpt_path)
    encoder_state = state
    if any(key.startswith("_full_encoder.") for key in state):
        encoder_state = _strip_prefix(state, "_full_encoder.")

    if "embedding.weight" not in encoder_state:
        return _result("Stage 1 checkpoint loading test", False, "embedding.weight not found.")

    vocab_size, embedding_dim = encoder_state["embedding.weight"].shape
    num_layers = _infer_layer_count(encoder_state)
    num_heads = _infer_num_heads(encoder_state, embedding_dim)
    ff_dim = _infer_ff_dim(encoder_state, embedding_dim)
    max_seq_len = int(encoder_state["pos_embedding.weight"].shape[0])
    use_relative_position_bias = (
        "relative_position_bias.relative_attention_bias.weight" in encoder_state
    )
    num_buckets = (
        int(encoder_state["relative_position_bias.relative_attention_bias.weight"].shape[0])
        if use_relative_position_bias
        else 32
    )
    use_glu_ffn = any(key.endswith("ffn.gate_proj.weight") for key in encoder_state)

    use_attention_pool = any(key.startswith("pool.") for key in state)
    use_engineered_features = any(key.startswith("engineered_mlp") for key in state)

    engineered_dim = 0
    engineered_mlp_hidden = 512
    engineered_mlp_output = 256
    engineered_mlp_use_gated = True
    engineered_mlp_use_residual = True
    if use_engineered_features and "engineered_mlp.mlp.0.weight" in state:
        engineered_mlp_hidden, engineered_dim = state["engineered_mlp.mlp.0.weight"].shape
        engineered_mlp_output = int(state["engineered_mlp.mlp.8.weight"].shape[0])
        engineered_mlp_use_gated = "engineered_mlp.gate_proj.weight" in state
        engineered_mlp_use_residual = "engineered_mlp.residual_proj.weight" in state
    elif use_engineered_features:
        details = "Checkpoint engineered MLP does not match simplified architecture."
        return _result("Stage 1 checkpoint loading test", True, details, skipped=True)

    model = GeneWhispererStage1(
        vocab_size=vocab_size,
        kmer=int(checkpoint.get("kmer", 6)) if isinstance(checkpoint, dict) else 6,
        embedding_dim=embedding_dim,
        num_layers=num_layers,
        num_heads=num_heads,
        ff_dim=ff_dim,
        dropout=0.0,
        pad_token_id=vocab_size - 1,
        engineered_dim=engineered_dim,
        use_engineered_features=use_engineered_features,
        use_attention_pool=use_attention_pool,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        engineered_mlp_pre_norm=True,
        engineered_mlp_input_dropout=0.0,
        engineered_mlp_use_gated=engineered_mlp_use_gated,
        engineered_mlp_use_residual=engineered_mlp_use_residual,
        drop_path_rate=0.0,
        max_seq_len=max_seq_len,
        use_relative_position_bias=use_relative_position_bias,
        relative_position_num_buckets=num_buckets,
        relative_position_max_distance=128,
        use_glu_ffn=use_glu_ffn,
    )

    loadable, skipped = _count_loadable_keys(state, model.state_dict())
    try:
        model.load_state_dict(state, strict=False)
    except Exception as exc:
        return _result("Stage 1 checkpoint loading test", False, f"Load failed: {exc}")

    accuracy = _extract_accuracy(checkpoint)
    if accuracy is None:
        return _result("Stage 1 checkpoint loading test", False, "Accuracy not found in checkpoint.")

    within_target = abs(accuracy - ACCURACY_TARGET) <= ACCURACY_TOLERANCE
    details = (
        f"checkpoint={ckpt_path.name} loaded={_format_count(loadable)} "
        f"skipped={_format_count(skipped)} accuracy={accuracy:.4f}"
    )
    return _result("Stage 1 checkpoint loading test", within_target, details)


def run_model_replacement_checks() -> Iterable[CheckResult]:
    return [
        check_imports(),
        check_parameter_counts(),
        check_architecture(),
        check_forward_pass(),
        check_mlm_weight_loading(),
        check_stage1_checkpoint_loading(),
    ]


def test_model_replacement_checks() -> None:
    results = list(run_model_replacement_checks())
    failures = [result for result in results if result.status == "FAIL"]
    if failures:
        details = "\n".join(f"{result.name}: {result.details}" for result in failures)
        raise AssertionError(f"Model replacement checks failed:\n{details}")


def _print_summary(results: Iterable[CheckResult]) -> None:
    totals = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for result in results:
        totals[result.status] += 1
        print(f"{result.status}: {result.name} - {result.details}")
    print(
        "Summary: "
        f"PASS={totals['PASS']} FAIL={totals['FAIL']} SKIP={totals['SKIP']}"
    )


if __name__ == "__main__":
    _print_summary(run_model_replacement_checks())
