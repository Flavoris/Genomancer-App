"""Promoter prediction via Stage-1 ensemble over multiple k-mer models."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

import torch
import yaml

from dataset import KmerVocabulary, compute_tnc, compute_pseeiip
from model import GeneWhispererStage1, MultiScaleEnsemble

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


def infer_model_architecture(checkpoint_path: Path) -> dict:
    """
    Infer model architecture parameters from checkpoint state_dict.

    Returns a dict of config overrides to match the checkpoint's architecture.
    """
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = checkpoint.get("model_state", checkpoint)

    overrides = {}

    # Count encoder layers
    encoder_layers = [k for k in state.keys() if "_full_encoder.layers." in k]
    if encoder_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in encoder_layers)
        overrides["transformer_layers"] = len(layer_nums)

    # Count TCN levels
    tcn_layers = [k for k in state.keys() if "tcn.network." in k]
    if tcn_layers:
        tcn_nums = set(int(k.split("network.")[1].split(".")[0]) for k in tcn_layers)
        overrides["tcn_levels"] = len(tcn_nums)

    # Count multiscale convs (infer kernel sizes from conv weights)
    mscale_convs = [k for k in state.keys() if "multiscale.convs." in k and ".weight" in k]
    if mscale_convs:
        kernel_sizes = []
        for k in sorted(mscale_convs):
            weight = state[k]
            kernel_sizes.append(weight.shape[2])  # Conv1d weight shape: (out_ch, in_ch, kernel)
        overrides["multiscale_kernels"] = kernel_sizes

    # Infer post-CNN transformer layers
    post_cnn_layers = [k for k in state.keys() if "post_cnn_transformer.layers." in k]
    if post_cnn_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in post_cnn_layers)
        overrides["post_cnn_transformer_layers"] = len(layer_nums)

    # Infer embedding dim from embedding weight
    if "embedding.weight" in state:
        overrides["embedding_dim"] = state["embedding.weight"].shape[1]

    # Infer transformer FF dim from FFN weights
    ffn_keys = [k for k in state.keys() if "_full_encoder.layers.0.ffn.0.weight" in k]
    if ffn_keys:
        overrides["transformer_ff_dim"] = state[ffn_keys[0]].shape[0]

    # Infer TCN hidden from TCN output
    tcn_out_keys = [k for k in state.keys() if "tcn.network.0.conv1.conv.weight" in k]
    if tcn_out_keys:
        overrides["tcn_hidden"] = state[tcn_out_keys[0]].shape[0]

    # Infer multiscale channels
    if mscale_convs:
        overrides["multiscale_channels"] = state[sorted(mscale_convs)[0]].shape[0]

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

    # Infer fusion hidden
    fusion_keys = [k for k in state.keys() if "fusion.seq_proj.weight" in k]
    if fusion_keys:
        overrides["fusion_hidden"] = state[fusion_keys[0]].shape[0]

    return overrides


def load_checkpoint(model: GeneWhispererStage1, path: Path, device: torch.device) -> None:
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    state_dict = checkpoint.get("model_state", checkpoint)

    # Use strict=True to ensure all weights load correctly
    # Log any issues for debugging
    try:
        incompatible = model.load_state_dict(state_dict, strict=True)
        if incompatible.missing_keys:
            LOGGER.warning("Missing keys: %s", incompatible.missing_keys)
        if incompatible.unexpected_keys:
            LOGGER.warning("Unexpected keys: %s", incompatible.unexpected_keys)
    except RuntimeError as e:
        LOGGER.error("Failed to load checkpoint with strict=True: %s", e)
        LOGGER.info("Attempting load with strict=False...")
        incompatible = model.load_state_dict(state_dict, strict=False)
        if incompatible.missing_keys:
            LOGGER.warning("Missing keys (not loaded): %s", incompatible.missing_keys[:10])
        if incompatible.unexpected_keys:
            LOGGER.warning("Unexpected keys (ignored): %s", incompatible.unexpected_keys[:10])

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
    """Build GeneWhispererStage1 model with V3 architecture parameters."""
    embedding_dim = int(cfg.get("embedding_dim", 192))
    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.2))
    
    # TCN parameters
    use_tcn = bool(cfg.get("use_tcn", True))
    tcn_hidden = int(cfg.get("tcn_hidden", 256))
    tcn_levels = int(cfg.get("tcn_levels", 4))
    tcn_kernel = int(cfg.get("tcn_kernel", 3))
    multiscale_channels = int(cfg.get("multiscale_channels", 64))
    multiscale_kernels = tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15]))
    lstm_hidden = int(cfg.get("lstm_hidden", 192))
    
    # New V3 architecture parameters
    post_cnn_transformer_layers = int(cfg.get("post_cnn_transformer_layers", 3))
    engineered_mlp_hidden = int(cfg.get("engineered_mlp_hidden", 256))
    engineered_mlp_output = int(cfg.get("engineered_mlp_output", 128))
    fusion_hidden = int(cfg.get("fusion_hidden", 256))
    
    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        pad_token_id=vocab.pad_id,
        engineered_dim=int(cfg.get("engineered_dim", 128)),
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        # TCN parameters
        use_tcn=use_tcn,
        tcn_hidden=tcn_hidden,
        tcn_levels=tcn_levels,
        tcn_kernel=tcn_kernel,
        multiscale_channels=multiscale_channels,
        multiscale_kernels=multiscale_kernels,
        lstm_hidden=lstm_hidden,
        # V3 architecture parameters
        post_cnn_transformer_layers=post_cnn_transformer_layers,
        engineered_mlp_hidden=engineered_mlp_hidden,
        engineered_mlp_output=engineered_mlp_output,
        fusion_hidden=fusion_hidden,
        # Positional embedding max length
        max_seq_len=int(cfg.get("max_bp_len", 81)) - vocab.k + 1,
    ).to(device)
    return model


def compute_engineered_features(sequence: str, cfg: Optional[Dict] = None) -> torch.Tensor:
    """Compute engineered features: TNC(64) + PseEIIP(64) = 128."""
    tnc = compute_tnc(sequence)          # 64-dim
    pseeiip = compute_pseeiip(sequence)  # 64-dim
    if cfg is not None:
        if not bool(cfg.get("stage1_feature_enable_tnc", True)):
            tnc = torch.zeros_like(tnc)
        if not bool(cfg.get("stage1_feature_enable_pseeiip", True)):
            pseeiip = torch.zeros_like(pseeiip)
    return torch.cat([tnc, pseeiip], dim=0)


def prepare_inputs(
    sequence: str, vocab: KmerVocabulary, max_bp_len: int, cfg: Optional[Dict] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    tokens = vocab.tokenize(sequence, max_bp_len)
    engineered = compute_engineered_features(sequence, cfg)
    return tokens.unsqueeze(0), engineered.unsqueeze(0)


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
        # fallback to default naming convention
        kmer_config = {
            "3": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k3.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k3_vocab.json").resolve()),
            },
            "4": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k4.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k4_vocab.json").resolve()),
            },
            "6": {
                "checkpoint": str((script_dir / "../artifacts/checkpoints/stage1_k6.pt").resolve()),
                "vocab": str((script_dir / "../artifacts/vocabs/k6_vocab.json").resolve()),
            },
        }

    models = []
    batch_inputs = {}
    for k_str, paths in kmer_config.items():
        k = int(k_str)
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
        models.append(model)
        batch_inputs[k] = (tokens, engineered)

    if not models:
        raise RuntimeError("No models loaded for ensemble inference")

    # Load calibration if available and not disabled
    calibration = None
    if not args.no_calibration:
        calibration_path = Path(args.calibration) if args.calibration else None
        calibration = load_calibration(calibration_path)

    # Perform inference with calibration support
    with torch.no_grad():
        if calibration is not None:
            # Use logits and apply temperature scaling
            all_probs = []
            for model in models:
                encoder = getattr(model, "encoder", None)
                k = getattr(encoder, "k", None)
                if k is None or k not in batch_inputs:
                    continue
                tokens, engineered = batch_inputs[k]
                _, logits = model(tokens, engineered_features=engineered, return_logits=True)
                # Apply temperature scaling
                probs, _ = apply_calibration(logits, calibration)
                all_probs.append(probs)

            # Average calibrated probabilities
            stacked = torch.stack(all_probs, dim=0)
            prob = stacked.mean(dim=0).squeeze().item()
            thr = calibration.get("calibrated_threshold", 0.5)
            LOGGER.info("Using calibrated inference: T=%.4f, threshold=%.4f", calibration["temperature"], thr)
        else:
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

            ensemble = MultiScaleEnsemble(models).to(device)
            outputs = ensemble(batch_inputs)
            prob = outputs.squeeze().item()

    label = "promoter" if prob >= thr else "non-promoter"
    result = {"probability": prob, "label": label, "threshold": thr}
    if calibration is not None:
        result["calibrated"] = True
        result["temperature"] = calibration["temperature"]
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
