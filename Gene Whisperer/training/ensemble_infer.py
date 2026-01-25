"""Promoter prediction via Stage-1 model with BPE tokenization."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import torch
import yaml

from bpe_tokenizer import DNABPETokenizer
from dataset import compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc, reverse_complement
from tta import get_reverse_complement_tokens, predict_with_tta, TTAWrapper
from model import GeneWhispererStage1

LOGGER = logging.getLogger("gene_whisperer.ensemble")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

# Default calibration file path (relative to training directory)
DEFAULT_CALIBRATION_PATH = Path(__file__).resolve().parent.parent / "artifacts" / "calibration_stage1.json"


def load_calibration(calibration_path: Optional[Path] = None) -> Optional[Dict]:
    """
    Load temperature calibration from JSON file.

    Args:
        calibration_path: Path to calibration JSON. If None, uses default path.

    Returns:
        Dict with 'temperature' and 'calibrated_threshold', or None if not found.
    """
    path = calibration_path or DEFAULT_CALIBRATION_PATH
    if not path.exists():
        LOGGER.debug("No calibration file found at %s", path)
        return None

    try:
        with path.open("r", encoding="utf-8") as f:
            calibration = json.load(f)

        if "temperature" not in calibration:
            LOGGER.warning("Calibration file missing 'temperature' key: %s", path)
            return None

        LOGGER.info(
            "Loaded calibration: T=%.4f, threshold=%.4f from %s",
            calibration["temperature"],
            calibration.get("calibrated_threshold", 0.5),
            path,
        )
        return calibration
    except (json.JSONDecodeError, IOError) as e:
        LOGGER.warning("Failed to load calibration from %s: %s", path, e)
        return None



def apply_calibration(
    logits: torch.Tensor,
    calibration: Optional[Dict],
) -> Tuple[torch.Tensor, float]:
    """
    Apply temperature scaling to logits and return calibrated probabilities.

    Args:
        logits: Raw model logits (pre-sigmoid)
        calibration: Dict with 'temperature' key, or None for no calibration

    Returns:
        (probs, threshold): Calibrated probabilities and decision threshold
    """
    if calibration is None:
        # No calibration: use raw sigmoid and default threshold
        probs = torch.sigmoid(logits)
        return probs, 0.5

    temperature = calibration["temperature"]
    threshold = calibration.get("calibrated_threshold", 0.5)

    # Apply temperature scaling
    calibrated_logits = logits / temperature
    probs = torch.sigmoid(calibrated_logits)

    return probs, threshold


LEGACY_STATE_PREFIXES = (
    "feature_extractor.",
    "post_cnn_transformer.",
    "multiscale.",
    "tcn.",
    "fusion.",
)


def _extract_state_dict(checkpoint: dict) -> dict:
    if "model_state" in checkpoint:
        return checkpoint["model_state"]
    if "model_state_dict" in checkpoint:
        return checkpoint["model_state_dict"]
    if "state_dict" in checkpoint:
        return checkpoint["state_dict"]
    return checkpoint


def _normalize_state_dict(state: dict) -> dict:
    if any(key.startswith("module.") for key in state.keys()):
        return {key.removeprefix("module."): value for key, value in state.items()}
    return state


def _is_legacy_state(state: dict) -> bool:
    return any(key.startswith(LEGACY_STATE_PREFIXES) for key in state.keys())


def _filter_compatible_state(
    model: torch.nn.Module,
    state: dict,
) -> tuple[dict, list[str], list[str]]:
    model_state = model.state_dict()
    compatible: dict = {}
    skipped: list[str] = []
    for key, value in state.items():
        target = model_state.get(key)
        if target is None:
            skipped.append(key)
            continue
        if getattr(target, "shape", None) != getattr(value, "shape", None):
            skipped.append(key)
            continue
        compatible[key] = value
    missing = [key for key in model_state.keys() if key not in compatible]
    return compatible, skipped, missing


def infer_model_architecture(checkpoint_path: Path) -> dict:
    """
    Infer model architecture parameters from checkpoint state_dict.

    Returns a dict of config overrides to match the checkpoint's architecture.
    """
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = _normalize_state_dict(_extract_state_dict(checkpoint))

    overrides = {}

    # Count encoder layers
    encoder_layers = [k for k in state.keys() if "_full_encoder.layers." in k]
    if encoder_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in encoder_layers)
        overrides["transformer_layers"] = len(layer_nums)

    # Infer embedding dim from embedding weight
    if "embedding.weight" in state:
        overrides["embedding_dim"] = state["embedding.weight"].shape[1]
    elif "_full_encoder.embedding.weight" in state:
        overrides["embedding_dim"] = state["_full_encoder.embedding.weight"].shape[1]

    if "_full_encoder.pos_embedding.weight" in state:
        overrides["max_seq_len"] = state["_full_encoder.pos_embedding.weight"].shape[0]

    # Infer transformer FF dim from FFN weights
    ffn_keys = [k for k in state.keys() if "_full_encoder.layers.0.ffn.0.weight" in k]
    if ffn_keys:
        overrides["transformer_ff_dim"] = state[ffn_keys[0]].shape[0]

    # Infer engineered MLP dimensions
    eng_mlp_keys = [k for k in state.keys() if "engineered_mlp.mlp." in k and "weight" in k]
    if eng_mlp_keys:
        # First linear layer: input -> hidden
        first_linear = [k for k in eng_mlp_keys if "mlp.0.weight" in k]
        if first_linear:
            overrides["engineered_mlp_hidden"] = state[first_linear[0]].shape[0]
        # Last linear layer: hidden -> output
        last_linear = [k for k in eng_mlp_keys if "mlp.8.weight" in k]
        if last_linear:
            overrides["engineered_mlp_output"] = state[last_linear[0]].shape[0]

    return overrides


def load_checkpoint(model: GeneWhispererStage1, path: Path, device: torch.device) -> None:
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    state_dict = _normalize_state_dict(_extract_state_dict(checkpoint))

    if _is_legacy_state(state_dict):
        LOGGER.warning(
            "Legacy Stage 1 checkpoint detected; loading compatible encoder weights only."
        )
        legacy_state = {
            key: value
            for key, value in state_dict.items()
            if key.startswith(("_full_encoder.", "embedding."))
        }
        state_dict = legacy_state

    compatible, skipped, missing = _filter_compatible_state(model, state_dict)
    model.load_state_dict(compatible, strict=False)
    LOGGER.info(
        "Loaded %d tensors from %s (skipped=%d, missing=%d)",
        len(compatible),
        path,
        len(skipped),
        len(missing),
    )
    if skipped:
        preview = ", ".join(skipped[:5])
        suffix = "…" if len(skipped) > 5 else ""
        LOGGER.warning("Skipped keys (preview): %s%s", preview, suffix)

    # Attach best_threshold from checkpoint for deployment-safe inference
    if isinstance(checkpoint, dict) and "best_threshold" in checkpoint:
        model.best_threshold = float(checkpoint["best_threshold"])
        LOGGER.info("Loaded best_threshold=%.4f from checkpoint", model.best_threshold)
    else:
        model.best_threshold = None
        LOGGER.warning("No best_threshold found in checkpoint %s", path)

    model.to(device)
    model.eval()
    LOGGER.info("Loaded model from %s", path)


def load_vocab(vocab_path: Path) -> DNABPETokenizer:
    return DNABPETokenizer.load(str(vocab_path))


def _resolve_pooling_type(pooling_type: Optional[str]) -> str:
    """Normalize pooling type to supported values."""
    normalized = str(pooling_type or "attention").strip().lower()
    if normalized in {"attention", "mean"}:
        return normalized
    if normalized in {"cls", "max"}:
        LOGGER.warning(
            "pooling_type '%s' is not supported; using mean pooling",
            normalized,
        )
        return "mean"
    LOGGER.warning(
        "pooling_type '%s' is not recognized; using attention pooling",
        normalized,
    )
    return "attention"


def build_model(cfg: Dict, vocab: DNABPETokenizer, device: torch.device) -> GeneWhispererStage1:
    """Build simplified GeneWhispererStage1 model."""
    simplified_cfg = cfg.get("simplified_model") if isinstance(cfg.get("simplified_model"), dict) else {}
    pooling_type = _resolve_pooling_type(simplified_cfg.get("pooling_type", "attention"))
    classifier_hidden = simplified_cfg.get("classifier_hidden")
    if classifier_hidden is None:
        classifier_hidden = 256
    classifier_dropout = simplified_cfg.get("classifier_dropout")
    if classifier_dropout is None:
        classifier_dropout = 0.15

    max_token_len = int(cfg.get("max_token_len", 24))
    max_seq_len = int(cfg.get("max_seq_len", max_token_len))
    stage1_drop_path_rate = cfg.get("stage1_drop_path_rate")
    if stage1_drop_path_rate is None:
        stage1_drop_path_rate = cfg.get("drop_path_rate", 0.1)
    ffn_type = cfg.get("ffn_type")
    if ffn_type is not None:
        ffn_type = str(ffn_type)
    norm_type = cfg.get("norm_type")
    norm_type = "layernorm" if norm_type is None else str(norm_type)
    ffn_mult = cfg.get("ffn_mult")
    if ffn_mult is not None:
        ffn_mult = float(ffn_mult)

    model = GeneWhispererStage1(
        vocab_size=len(vocab),
        embedding_dim=int(cfg.get("embedding_dim", 384)),
        num_layers=int(cfg.get("transformer_layers", 12)),
        num_heads=int(cfg.get("transformer_heads", 12)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1536)),
        dropout=float(cfg.get("transformer_dropout", 0.12)),
        pad_token_id=vocab.pad_id,
        engineered_dim=int(cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        pooling_type=pooling_type,
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        classifier_hidden=int(classifier_hidden),
        classifier_dropout=float(classifier_dropout),
        drop_path_rate=float(stage1_drop_path_rate),
        max_seq_len=max_seq_len,
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", False)),
        relative_position_num_buckets=int(cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", False)),
        glu_activation=str(cfg.get("glu_activation", "gelu")),
        use_rope=bool(cfg.get("use_rope", False)),
        rope_base=float(cfg.get("rope_base", 10000.0)),
        ffn_type=ffn_type,
        norm_type=norm_type,
        ffn_mult=ffn_mult,
    ).to(device)
    return model


def compute_engineered_features(sequence: str, cfg: Optional[Dict] = None) -> torch.Tensor:
    """Compute engineered features: TNC(64) + PseEIIP(64) + CKSNAP(96) + PSTNP(64) = 288."""
    tnc = compute_tnc(sequence)          # 64-dim
    pseeiip = compute_pseeiip(sequence)  # 64-dim
    cksnap = compute_cksnap(sequence)    # 96-dim
    pstnp = compute_pstnp(sequence)      # 64-dim
    if cfg is not None:
        if not bool(cfg.get("stage1_feature_enable_tnc", True)):
            tnc = torch.zeros_like(tnc)
        if not bool(cfg.get("stage1_feature_enable_pseeiip", True)):
            pseeiip = torch.zeros_like(pseeiip)
        if not bool(cfg.get("stage1_feature_enable_cksnap", True)):
            cksnap = torch.zeros_like(cksnap)
        if not bool(cfg.get("stage1_feature_enable_pstnp", True)):
            pstnp = torch.zeros_like(pstnp)
    return torch.cat([tnc, pseeiip, cksnap, pstnp], dim=0)


def prepare_inputs(
    sequence: str, vocab: DNABPETokenizer, max_token_len: int, cfg: Optional[Dict] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    tokens = vocab.tokenize_and_pad(sequence, max_token_len)
    engineered = compute_engineered_features(sequence, cfg)
    return tokens.unsqueeze(0), engineered.unsqueeze(0)


def load_stage1_report_weight(report_path: Path) -> Optional[float]:
    if not report_path.exists():
        return None
    try:
        with report_path.open("r", encoding="utf-8") as handle:
            report = json.load(handle)
    except (json.JSONDecodeError, IOError) as exc:
        LOGGER.warning("Failed to load report %s: %s", report_path, exc)
        return None

    val_metrics = report.get("val_metrics") or {}
    for key in ("mcc", "mcc_best", "acc_best", "val_acc@best", "accuracy"):
        if key in val_metrics:
            return float(val_metrics[key])

    for key in ("val_mcc", "val_mcc_best", "val_acc@best", "val_accuracy", "val_acc_best"):
        if key in report:
            return float(report[key])

    return None


def resolve_ensemble_weights(checkpoints: List[str], artifacts_dir: Path) -> Optional[List[float]]:
    """Load weights from stage1 report for multi-checkpoint ensemble."""
    report_path = artifacts_dir / "stage1_report_bpe.json"
    weight = load_stage1_report_weight(report_path)
    if weight is None:
        return None
    # For single model, weight is 1.0
    return [1.0 / len(checkpoints)] * len(checkpoints)


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage-1 BPE model inference")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--sequence", required=True, help="Raw DNA sequence (A/C/G/T)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Promoter decision threshold")
    parser.add_argument("--checkpoint", type=str, default=None, help="Path to model checkpoint")
    parser.add_argument(
        "--calibration",
        type=str,
        default=None,
        help="Path to calibration JSON. If not specified, auto-loads from artifacts/calibration_stage1.json if it exists.",
    )
    parser.add_argument(
        "--no-calibration",
        action="store_true",
        help="Disable calibration even if calibration file exists.",
    )
    parser.add_argument(
        "--tta",
        action="store_true",
        help="Enable test-time augmentation using reverse complement.",
    )
    parser.add_argument(
        "--tta-aggregation",
        type=str,
        choices=["mean", "geometric_mean"],
        default="mean",
        help="TTA aggregation method (default: mean).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    sequence = args.sequence.strip().upper()
    max_token_len = int(cfg.get("max_token_len", 24))
    device = select_device()

    # Load BPE vocabulary
    bpe_vocab_path = cfg.get("bpe_vocab_path", "../artifacts/vocabs/bpe_vocab.json")
    vocab_path = Path(bpe_vocab_path)
    if not vocab_path.is_absolute():
        vocab_path = (script_dir / vocab_path).resolve()
    if not vocab_path.exists():
        raise FileNotFoundError(f"BPE vocab not found: {vocab_path}")
    vocab = load_vocab(vocab_path)

    # Load model checkpoint
    checkpoint_path = args.checkpoint
    if checkpoint_path is None:
        checkpoint_name = cfg.get("stage1_checkpoint_name", "stage1_bpe.pt")
        checkpoint_path = (script_dir / "../artifacts/checkpoints" / checkpoint_name).resolve()
    else:
        checkpoint_path = Path(checkpoint_path).resolve()
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    # Infer architecture from checkpoint to handle config mismatches
    arch_overrides = infer_model_architecture(checkpoint_path)
    model_cfg = {**cfg, **arch_overrides}
    model = build_model(model_cfg, vocab, device)
    load_checkpoint(model, checkpoint_path, device)

    # Prepare inputs
    tokens, engineered = prepare_inputs(sequence, vocab, max_token_len, cfg)
    tokens = tokens.to(device)
    engineered = engineered.to(device)

    # Load calibration if available and not disabled
    calibration = None
    if not args.no_calibration:
        calibration_path = Path(args.calibration) if args.calibration else None
        calibration = load_calibration(calibration_path)

    # Check TTA settings from config or CLI
    inference_cfg = cfg.get("inference") if isinstance(cfg.get("inference"), dict) else {}
    use_tta = args.tta or bool(inference_cfg.get("use_tta", False))
    tta_aggregation = args.tta_aggregation or inference_cfg.get("tta_aggregation", "mean")

    if use_tta:
        LOGGER.info("TTA enabled with %s aggregation", tta_aggregation)

    # Perform inference
    with torch.no_grad():
        # Determine decision threshold
        if calibration is None:
            model_threshold = getattr(model, "best_threshold", None)
            thr = model_threshold if model_threshold is not None else cfg.get("stage1_decision_threshold", args.threshold)
            LOGGER.info("Using decision threshold: %.4f (no calibration)", thr)

        if use_tta:
            probs = predict_with_tta(
                model=model,
                tokens=tokens,
                engineered_features=engineered,
                vocab=vocab,
                max_token_len=max_token_len,
                compute_features_fn=lambda seq: compute_engineered_features(seq, cfg),
                aggregation=tta_aggregation,
            )
        else:
            logits = model(tokens, engineered_features=engineered)
            probs = torch.sigmoid(logits)

        if calibration is not None:
            eps = 1e-7
            logits_from_probs = torch.log(probs / (1 - probs + eps) + eps)
            probs, _ = apply_calibration(logits_from_probs, calibration)
            thr = calibration.get("calibrated_threshold", 0.5)
            LOGGER.info(
                "Using calibrated inference: T=%.4f, threshold=%.4f",
                calibration["temperature"],
                thr,
            )

        prob = probs.squeeze().item()

    label = "promoter" if prob >= thr else "non-promoter"
    result = {"probability": prob, "label": label, "threshold": thr}
    if calibration is not None:
        result["calibrated"] = True
        result["temperature"] = calibration["temperature"]
    if use_tta:
        result["tta_enabled"] = True
        result["tta_aggregation"] = tta_aggregation
    print(json.dumps(result))


def check_best_threshold(checkpoint_path: str, config_path: str = "config.yaml") -> None:
    """
    Unit check: verify that best_threshold loads correctly from a checkpoint.

    Usage:
        python ensemble_infer.py --check-threshold path/to/checkpoint.pt
    """
    import sys

    script_dir = Path(__file__).resolve().parent
    ckpt_path = Path(checkpoint_path).resolve()
    cfg_path = Path(config_path)
    if not cfg_path.is_absolute():
        cfg_path = (script_dir / cfg_path).resolve()

    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    device = select_device()

    bpe_vocab_path = cfg.get("bpe_vocab_path", "../artifacts/vocabs/bpe_vocab.json")
    vocab_path = Path(bpe_vocab_path)
    if not vocab_path.is_absolute():
        vocab_path = (script_dir / vocab_path).resolve()
    if not vocab_path.exists():
        LOGGER.error("BPE vocab not found at %s", vocab_path)
        sys.exit(1)

    vocab = load_vocab(vocab_path)
    model = build_model(cfg, vocab, device)
    load_checkpoint(model, ckpt_path, device)

    # Assert best_threshold is loaded
    assert hasattr(model, "best_threshold"), "model.best_threshold attribute not set!"
    assert model.best_threshold is not None, "model.best_threshold is None!"

    LOGGER.info("Unit check passed: model.best_threshold = %.4f", model.best_threshold)


def verify_checkpoint_compatibility(
    checkpoint_path: str,
    config_path: str = "config.yaml",
) -> bool:
    """Verify that build_model creates architecture compatible with checkpoint."""
    import yaml

    script_dir = Path(__file__).resolve().parent
    cfg_path = Path(config_path)
    if not cfg_path.is_absolute():
        cfg_path = (script_dir / cfg_path).resolve()

    with cfg_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    ckpt_path = Path(checkpoint_path).resolve()
    if not ckpt_path.exists():
        print(f"Checkpoint not found: {ckpt_path}")
        return False

    checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    ckpt_state = _normalize_state_dict(_extract_state_dict(checkpoint))

    # Load BPE vocab
    bpe_vocab_path = cfg.get("bpe_vocab_path", "../artifacts/vocabs/bpe_vocab.json")
    vocab_path = Path(bpe_vocab_path)
    if not vocab_path.is_absolute():
        vocab_path = (script_dir / vocab_path).resolve()
    if not vocab_path.exists():
        print(f"BPE vocab not found: {vocab_path}")
        return False

    vocab = load_vocab(vocab_path)

    # Build model
    device = torch.device("cpu")
    model = build_model(cfg, vocab, device)
    model_state = model.state_dict()

    # Compare keys
    ckpt_keys = set(ckpt_state.keys())
    model_keys = set(model_state.keys())

    missing = model_keys - ckpt_keys
    unexpected = ckpt_keys - model_keys

    print(f"Checkpoint keys: {len(ckpt_keys)}")
    print(f"Model keys: {len(model_keys)}")
    print(f"Missing in checkpoint: {len(missing)}")
    print(f"Unexpected in checkpoint: {len(unexpected)}")

    if missing:
        print(f"Missing keys (first 10): {list(missing)[:10]}")
    if unexpected:
        print(f"Unexpected keys (first 10): {list(unexpected)[:10]}")

    # Check shape compatibility
    shape_mismatches = []
    for key in ckpt_keys & model_keys:
        if ckpt_state[key].shape != model_state[key].shape:
            shape_mismatches.append((key, ckpt_state[key].shape, model_state[key].shape))

    if shape_mismatches:
        print(f"Shape mismatches: {len(shape_mismatches)}")
        for key, ckpt_shape, model_shape in shape_mismatches[:5]:
            print(f"  {key}: ckpt={ckpt_shape}, model={model_shape}")

    if not missing and not unexpected and not shape_mismatches:
        print("✓ Checkpoint is FULLY COMPATIBLE with model architecture!")
        return True

    print("✗ Checkpoint has INCOMPATIBILITIES")
    return False


if __name__ == "__main__":
    import sys

    # Check for --check-threshold flag
    if len(sys.argv) >= 3 and sys.argv[1] == "--check-threshold":
        config_arg = "config.yaml"
        if len(sys.argv) >= 4:
            config_arg = sys.argv[3]
        check_best_threshold(sys.argv[2], config_arg)
    elif len(sys.argv) >= 3 and sys.argv[1] == "--verify":
        config_arg = "config.yaml"
        if len(sys.argv) >= 4:
            config_arg = sys.argv[3]
        ok = verify_checkpoint_compatibility(sys.argv[2], config_arg)
        if not ok:
            sys.exit(1)
    else:
        main()
