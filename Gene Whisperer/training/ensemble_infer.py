"""Promoter prediction via Stage-1 ensemble over multiple k-mer models."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import torch
import yaml

from dataset import KmerVocabulary, compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc
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


def load_vocab(vocab_path: Path) -> KmerVocabulary:
    return KmerVocabulary.load(vocab_path)


def build_model(cfg: Dict, vocab: KmerVocabulary, device: torch.device) -> GeneWhispererStage1:
    """Build simplified GeneWhispererStage1 model."""
    simplified_cfg = cfg.get("simplified_model") if isinstance(cfg.get("simplified_model"), dict) else {}
    pooling_type = simplified_cfg.get("pooling_type")
    classifier_hidden = simplified_cfg.get("classifier_hidden")
    classifier_dropout = simplified_cfg.get("classifier_dropout")

    max_bp_len = int(cfg.get("max_bp_len", 81))
    max_seq_len = int(cfg.get("max_seq_len", max_bp_len - vocab.k + 1))

    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
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
        classifier_hidden=classifier_hidden,
        classifier_dropout=classifier_dropout,
        drop_path_rate=float(cfg.get("drop_path_rate", 0.1)),
        max_seq_len=max_seq_len,
        use_relative_position_bias=bool(cfg.get("use_relative_position_bias", False)),
        relative_position_num_buckets=int(cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(cfg.get("use_glu_ffn", False)),
        glu_activation=str(cfg.get("glu_activation", "gelu")),
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
    sequence: str, vocab: KmerVocabulary, max_bp_len: int, cfg: Optional[Dict] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    tokens = vocab.tokenize(sequence, max_bp_len)
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


def resolve_ensemble_weights(kmers: List[int], artifacts_dir: Path) -> Optional[List[float]]:
    weights = []
    missing = []
    for k in kmers:
        report_path = artifacts_dir / f"stage1_report_k{k}.json"
        weight = load_stage1_report_weight(report_path)
        if weight is None:
            missing.append(k)
            weight = 1.0
        weights.append(float(weight))

    if missing:
        LOGGER.warning(
            "Missing validation metrics for k=%s; defaulting weights to 1.0 for those models",
            ",".join(str(k) for k in missing),
        )

    weight_sum = sum(weights)
    if weight_sum == 0.0:
        LOGGER.warning("Weighted ensemble requested but weights sum to 0; using uniform mean instead")
        return None

    return [weight / weight_sum for weight in weights]


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage-1 ensemble inference")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    parser.add_argument("--sequence", required=True, help="Raw DNA sequence (A/C/G/T)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Promoter decision threshold")
    parser.add_argument("--kmer_configs", type=str, help="JSON mapping of k to checkpoint/vocab paths")
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
        "--per_model",
        action="store_true",
        help="Print each model's probability alongside the ensemble probability.",
    )
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Weight ensemble outputs using validation metrics from artifacts/stage1_report_k{k}.json.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    sequence = args.sequence.strip().upper()
    max_bp_len = int(cfg.get("max_bp_len", 81))
    device = select_device()

    if args.kmer_configs:
        kmer_config = json.loads(args.kmer_configs)
    else:
        kmer_config = {}
        kmer_list = cfg.get("multi_scale_kmers", [6])
        stage1_ckpt_by_k = cfg.get("stage1_ckpt_by_k", {}) or {}
        vocab_cache_dir = cfg.get("vocab_cache_dir")
        if not vocab_cache_dir:
            LOGGER.warning("Missing vocab_cache_dir in config; cannot resolve vocab paths.")
        for raw_k in kmer_list:
            try:
                k = int(raw_k)
            except (TypeError, ValueError):
                LOGGER.warning("Skipping invalid k-mer value: %s", raw_k)
                continue
            ckpt_rel = stage1_ckpt_by_k.get(k, stage1_ckpt_by_k.get(str(k)))
            if not ckpt_rel or not vocab_cache_dir:
                LOGGER.warning("Skipping k=%s (missing checkpoint or vocab config)", k)
                continue
            ckpt_path = Path(ckpt_rel)
            if not ckpt_path.is_absolute():
                ckpt_path = (script_dir / ckpt_path).resolve()
            vocab_dir = Path(vocab_cache_dir)
            if not vocab_dir.is_absolute():
                vocab_dir = (script_dir / vocab_dir).resolve()
            vocab_path = (vocab_dir / f"k{k}_vocab.json").resolve()
            if not vocab_path.exists() or not ckpt_path.exists():
                LOGGER.warning("Skipping k=%s (missing vocab or checkpoint)", k)
                continue
            kmer_config[str(k)] = {"checkpoint": str(ckpt_path), "vocab": str(vocab_path)}

    models = []
    batch_inputs = {}
    for k_str, paths in kmer_config.items():
        try:
            k = int(k_str)
        except (TypeError, ValueError):
            LOGGER.warning("Skipping invalid k-mer entry: %s", k_str)
            continue
        vocab_path = Path(paths["vocab"]).resolve()
        checkpoint_path = Path(paths["checkpoint"]).resolve()
        if not vocab_path.exists() or not checkpoint_path.exists():
            LOGGER.warning("Skipping k=%s (missing vocab or checkpoint)", k_str)
            continue
        vocab = load_vocab(vocab_path)
        # Infer architecture from checkpoint to handle config mismatches
        arch_overrides = infer_model_architecture(checkpoint_path)
        model_cfg = {**cfg, **arch_overrides}
        model = build_model(model_cfg, vocab, device)
        load_checkpoint(model, checkpoint_path, device)
        tokens, engineered = prepare_inputs(sequence, vocab, max_bp_len, cfg)
        tokens = tokens.to(device)
        engineered = engineered.to(device)
        models.append((k, model))
        batch_inputs[k] = (tokens, engineered)

    if not models:
        raise RuntimeError("No models loaded for ensemble inference")

    ensemble_weights = None
    use_weighted = False
    if args.weighted:
        artifacts_dir = script_dir.parent / "artifacts"
        kmers = [k for k, _ in models]
        ensemble_weights = resolve_ensemble_weights(kmers, artifacts_dir)
        if ensemble_weights is None:
            LOGGER.warning("Weighted ensemble requested but weights unavailable; using unweighted averaging")
        else:
            use_weighted = True
            weight_summary = ", ".join(
                f"k{k}={weight:.4f}" for k, weight in zip(kmers, ensemble_weights)
            )
            LOGGER.info("Using weighted ensemble (%s)", weight_summary)

    # Load calibration if available and not disabled
    calibration = None
    if not args.no_calibration:
        calibration_path = Path(args.calibration) if args.calibration else None
        calibration = load_calibration(calibration_path)

    # Perform inference with calibration support
    with torch.no_grad():
        per_model_probs = []
        all_probs = []
        if calibration is None:
            # Original behavior: use ensemble without calibration
            # Determine decision threshold: use average of model thresholds or fallback
            model_thresholds = [
                getattr(m, "best_threshold", None) for m in models if getattr(m, "best_threshold", None) is not None
            ]
            if model_thresholds:
                thr = sum(model_thresholds) / len(model_thresholds)
            else:
                thr = cfg.get("stage1_decision_threshold", args.threshold)
            LOGGER.info("Using decision threshold: %.4f (no calibration)", thr)

        for k, model in models:
            tokens, engineered = batch_inputs.get(k, (None, None))
            if tokens is None:
                continue
            logits = model(tokens, engineered_features=engineered)
            if calibration is not None:
                probs, _ = apply_calibration(logits, calibration)
            else:
                probs = torch.sigmoid(logits)
            all_probs.append(probs)
            per_model_probs.append((k, probs.squeeze().item()))

        if not all_probs:
            raise RuntimeError("No model outputs available for ensemble inference")

        stacked = torch.stack(all_probs, dim=0)
        if use_weighted and ensemble_weights is not None:
            weights = torch.tensor(ensemble_weights, device=stacked.device, dtype=stacked.dtype)
            probs = (weights.view(-1, 1, 1) * stacked).sum(dim=0)
        else:
            probs = stacked.mean(dim=0)
        prob = probs.squeeze().item()

        if calibration is not None:
            thr = calibration.get("calibrated_threshold", 0.5)
            LOGGER.info(
                "Using calibrated inference: T=%.4f, threshold=%.4f",
                calibration["temperature"],
                thr,
            )

    label = "promoter" if prob >= thr else "non-promoter"
    result = {"probability": prob, "label": label, "threshold": thr}
    if calibration is not None:
        result["calibrated"] = True
        result["temperature"] = calibration["temperature"]
    if args.per_model:
        result["per_model"] = [{"k": k, "probability": p} for k, p in per_model_probs]
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

    # Infer k-mer from checkpoint name (e.g., stage1_k6.pt -> k=6)
    ckpt_name = ckpt_path.stem
    kmer = 6  # default
    if "_k" in ckpt_name:
        try:
            kmer = int(ckpt_name.split("_k")[-1].split(".")[0])
        except ValueError:
            pass

    vocab_path = script_dir / f"../artifacts/vocabs/k{kmer}_vocab.json"
    if not vocab_path.exists():
        LOGGER.error("Vocab not found at %s", vocab_path)
        sys.exit(1)

    vocab = load_vocab(vocab_path.resolve())
    model = build_model(cfg, vocab, device)
    load_checkpoint(model, ckpt_path, device)

    # Assert best_threshold is loaded
    assert hasattr(model, "best_threshold"), "model.best_threshold attribute not set!"
    assert model.best_threshold is not None, "model.best_threshold is None!"

    LOGGER.info("✓ Unit check passed: model.best_threshold = %.4f", model.best_threshold)


if __name__ == "__main__":
    import sys

    # Check for --check-threshold flag
    if len(sys.argv) >= 3 and sys.argv[1] == "--check-threshold":
        config_arg = "config.yaml"
        if len(sys.argv) >= 4:
            config_arg = sys.argv[3]
        check_best_threshold(sys.argv[2], config_arg)
    else:
        main()
