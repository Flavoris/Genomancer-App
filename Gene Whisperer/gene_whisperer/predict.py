"""Ensemble-ready Stage 1 inference API for Gene Whisperer.

Public API:
    load_stage1_ensemble(config_path, checkpoints_and_vocabs) -> Stage1Predictor
    Stage1Predictor.predict(sequence) -> dict {probability, label}

Example:
    >>> predictor = load_stage1_ensemble(
    ...     "training/config.yaml",
    ...     {
    ...         4: {"checkpoint": "artifacts/checkpoints/stage1_k4.pt",
    ...             "vocab": "artifacts/vocabs/k4_vocab.json"},
    ...         6: {"checkpoint": "artifacts/checkpoints/stage1_k6.pt",
    ...             "vocab": "artifacts/vocabs/k6_vocab.json"},
    ...     }
    ... )
    >>> result = predictor.predict("ACGTACGTACGT...")
    >>> print(result)  # {"probability": 0.87, "label": "promoter"}
"""
from __future__ import annotations

import json
import logging
import math
import statistics
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Union

import torch
import yaml

# Add training directory to path for imports
_PACKAGE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _PACKAGE_DIR.parent
_TRAINING_DIR = _PROJECT_ROOT / "training"
if str(_TRAINING_DIR) not in sys.path:
    sys.path.insert(0, str(_TRAINING_DIR))

from dataset import KmerVocabulary, compute_cksnap, compute_pseeiip, compute_pstnp, compute_tnc
from model import GeneWhispererStage1, MultiScaleEnsemble

LOGGER = logging.getLogger("gene_whisperer.predict")


def _select_device() -> torch.device:
    """Select best available device (CUDA > MPS > CPU)."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _compute_engineered_features(
    sequence: str,
    cfg: Optional[Dict] = None,
) -> torch.Tensor:
    """Compute engineered features: TNC(64) + PseEIIP(64) + CKSNAP(96) + PSTNP(64) = 288."""
    tnc = compute_tnc(sequence)
    pseeiip = compute_pseeiip(sequence)
    cksnap = compute_cksnap(sequence)
    pstnp = compute_pstnp(sequence)
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


def _load_vocab(vocab_path: Path) -> KmerVocabulary:
    """Load a KmerVocabulary from JSON file."""
    return KmerVocabulary.load(vocab_path)


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


def _infer_architecture_from_checkpoint(
    checkpoint_path: Path,
    device: torch.device,
) -> Dict:
    """
    Infer simplified model architecture parameters from checkpoint state dict.

    This allows loading checkpoints trained with different configs by extracting
    shapes for the transformer backbone and engineered MLP.
    """
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    state = _normalize_state_dict(_extract_state_dict(checkpoint))

    arch: Dict[str, object] = {}

    emb_key = "embedding.weight"
    if emb_key not in state:
        emb_key = "_full_encoder.embedding.weight"
    if emb_key in state:
        arch["embedding_dim"] = state[emb_key].shape[1]

    encoder_layers = [k for k in state.keys() if "_full_encoder.layers." in k]
    if encoder_layers:
        layer_nums = set(int(k.split("layers.")[1].split(".")[0]) for k in encoder_layers)
        arch["transformer_layers"] = len(layer_nums)

    if "_full_encoder.pos_embedding.weight" in state:
        arch["max_seq_len"] = state["_full_encoder.pos_embedding.weight"].shape[0]

    ffn_key = "_full_encoder.layers.0.ffn.0.weight"
    if ffn_key in state:
        arch["transformer_ff_dim"] = state[ffn_key].shape[0]

    eng_keys = [k for k in state.keys() if "engineered_mlp.mlp." in k and "weight" in k]
    if eng_keys:
        first_linear = [k for k in eng_keys if "mlp.0.weight" in k]
        if first_linear:
            arch["engineered_mlp_hidden"] = state[first_linear[0]].shape[0]
        last_linear = [k for k in eng_keys if "mlp.8.weight" in k]
        if last_linear:
            arch["engineered_mlp_output"] = state[last_linear[0]].shape[0]

    if _is_legacy_state(state):
        arch["legacy_checkpoint"] = True

    LOGGER.info("Inferred architecture from checkpoint: %s", arch)
    return arch


def _build_model(
    cfg: Dict,
    vocab: KmerVocabulary,
    device: torch.device,
    arch_overrides: Optional[Dict] = None,
) -> torch.nn.Module:
    """Build Stage 1 model with config parameters.

    Args:
        cfg: Base configuration dictionary
        vocab: Vocabulary for this k-mer size
        device: Torch device
        arch_overrides: Optional architecture overrides inferred from checkpoint
    """
    # Merge config with overrides (overrides take precedence)
    effective_cfg = dict(cfg)
    if arch_overrides:
        effective_cfg.update(arch_overrides)
        if "legacy_checkpoint" in effective_cfg:
            effective_cfg.pop("legacy_checkpoint", None)

    simplified_cfg = (
        effective_cfg.get("simplified_model")
        if isinstance(effective_cfg.get("simplified_model"), dict)
        else {}
    )
    pooling_type = simplified_cfg.get("pooling_type")
    classifier_hidden = simplified_cfg.get("classifier_hidden")
    classifier_dropout = simplified_cfg.get("classifier_dropout")

    max_bp_len = int(effective_cfg.get("max_bp_len", 81))
    max_seq_len = int(effective_cfg.get("max_seq_len", max_bp_len - vocab.k + 1))

    model = GeneWhispererStage1(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=int(effective_cfg.get("embedding_dim", 384)),
        num_layers=int(effective_cfg.get("transformer_layers", 12)),
        num_heads=int(effective_cfg.get("transformer_heads", 12)),
        ff_dim=int(effective_cfg.get("transformer_ff_dim", 1536)),
        dropout=float(effective_cfg.get("transformer_dropout", 0.12)),
        pad_token_id=vocab.pad_id,
        engineered_dim=int(effective_cfg.get("engineered_dim", 288)),
        use_engineered_features=bool(effective_cfg.get("stage1_use_engineered_features", True)),
        pooling_type=pooling_type,
        engineered_mlp_hidden=int(effective_cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(effective_cfg.get("engineered_mlp_output", 128)),
        classifier_hidden=classifier_hidden,
        classifier_dropout=classifier_dropout,
        drop_path_rate=float(effective_cfg.get("drop_path_rate", 0.1)),
        max_seq_len=max_seq_len,
        use_relative_position_bias=bool(effective_cfg.get("use_relative_position_bias", False)),
        relative_position_num_buckets=int(effective_cfg.get("relative_position_num_buckets", 32)),
        relative_position_max_distance=int(effective_cfg.get("relative_position_max_distance", 128)),
        use_glu_ffn=bool(effective_cfg.get("use_glu_ffn", False)),
        glu_activation=str(effective_cfg.get("glu_activation", "gelu")),
    ).to(device)
    return model


def _load_checkpoint(
    model: torch.nn.Module,
    path: Path,
    device: torch.device,
) -> Optional[float]:
    """Load model weights from checkpoint file and return best_threshold if present."""
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    state_dict = _normalize_state_dict(_extract_state_dict(checkpoint))

    if _is_legacy_state(state_dict):
        LOGGER.warning(
            "Legacy Stage 1 checkpoint detected; loading compatible encoder weights only."
        )
        state_dict = {
            key: value
            for key, value in state_dict.items()
            if key.startswith(("_full_encoder.", "embedding."))
        }

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

    best_threshold = None
    if isinstance(checkpoint, dict) and "best_threshold" in checkpoint:
        best_threshold = float(checkpoint["best_threshold"])
        model.best_threshold = best_threshold
        LOGGER.info("Loaded best_threshold=%.4f from %s", best_threshold, path)
    else:
        model.best_threshold = None
        LOGGER.debug("No best_threshold found in %s", path)

    model.to(device)
    model.eval()
    LOGGER.info("Loaded model from %s", path)
    return best_threshold


def _resolve_decision_threshold(
    requested_threshold: Optional[float],
    checkpoint_thresholds: List[float],
) -> float:
    """Select a decision threshold from user input or checkpoint metadata."""
    if requested_threshold is not None:
        return float(requested_threshold)

    candidates: List[float] = []
    for threshold in checkpoint_thresholds:
        if threshold is None:
            continue
        try:
            value = float(threshold)
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            candidates.append(value)

    if candidates:
        chosen = float(statistics.median(candidates))
        if len(candidates) > 1:
            LOGGER.info("Using median best_threshold=%.4f from %d checkpoints", chosen, len(candidates))
        else:
            LOGGER.info("Using checkpoint best_threshold=%.4f", chosen)
        return chosen

    fallback = 0.5
    LOGGER.info("No checkpoint threshold found; using default %.2f", fallback)
    return fallback


class Stage1Predictor:
    """
    Callable predictor for Stage 1 promoter classification.

    Uses an ensemble of models with different k-mer sizes for robust prediction.

    Attributes:
        threshold: Decision threshold for promoter classification
    """

    def __init__(
        self,
        ensemble: MultiScaleEnsemble,
        vocabs: Dict[int, KmerVocabulary],
        cfg: Dict,
        device: torch.device,
        threshold: float = 0.5,
    ):
        """
        Initialize the predictor.

        Args:
            ensemble: MultiScaleEnsemble containing all k-mer models
            vocabs: Dict mapping k-mer size to vocabulary
            cfg: Configuration dictionary
            device: Torch device for inference
            threshold: Decision threshold for promoter classification
        """
        self.ensemble = ensemble
        self.vocabs = vocabs
        self.cfg = cfg
        self.device = device
        self.threshold = threshold
        self.max_bp_len = int(cfg.get("max_bp_len", 81))

    def _prepare_inputs(
        self,
        sequence: str,
        vocab: KmerVocabulary,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Tokenize and compute features for a single sequence."""
        tokens = vocab.tokenize(sequence, self.max_bp_len)
        engineered = _compute_engineered_features(sequence, self.cfg)
        return tokens.unsqueeze(0), engineered.unsqueeze(0)

    def predict(self, sequence: str) -> Dict[str, Union[float, str]]:
        """
        Predict promoter probability for a DNA sequence.

        Args:
            sequence: Raw DNA sequence (A/C/G/T characters)

        Returns:
            Dict with keys:
                - probability: float between 0 and 1
                - label: "promoter" or "non-promoter"
        """
        sequence = sequence.strip().upper()

        batch_inputs: Dict[int, tuple[torch.Tensor, torch.Tensor]] = {}
        for k, vocab in self.vocabs.items():
            tokens, engineered = self._prepare_inputs(sequence, vocab)
            tokens = tokens.to(self.device)
            engineered = engineered.to(self.device)
            batch_inputs[k] = (tokens, engineered)

        with torch.no_grad():
            outputs = self.ensemble(batch_inputs)

        prob = outputs.squeeze().item()
        label = "promoter" if prob >= self.threshold else "non-promoter"

        return {"probability": prob, "label": label}

    def __call__(self, sequence: str) -> Dict[str, Union[float, str]]:
        """Alias for predict() to make the predictor callable."""
        return self.predict(sequence)


def load_stage1_ensemble(
    config_path: Union[str, Path],
    checkpoints_and_vocabs: Dict[int, Dict[str, str]],
    device: Optional[torch.device] = None,
    threshold: Optional[float] = None,
) -> Stage1Predictor:
    """
    Load a Stage 1 ensemble predictor from config and checkpoints.

    Args:
        config_path: Path to config.yaml file
        checkpoints_and_vocabs: Dict mapping k-mer size to paths, e.g.:
            {
                4: {"checkpoint": "path/to/stage1_k4.pt", "vocab": "path/to/k4_vocab.json"},
                6: {"checkpoint": "path/to/stage1_k6.pt", "vocab": "path/to/k6_vocab.json"},
            }
        device: Optional torch device (auto-selects if None)
        threshold: Decision threshold override. If None, uses checkpoint best_threshold.

    Returns:
        Stage1Predictor callable that takes a sequence and returns predictions

    Raises:
        FileNotFoundError: If config, checkpoint, or vocab files don't exist
        RuntimeError: If no models could be loaded

    Example:
        >>> predictor = load_stage1_ensemble(
        ...     "training/config.yaml",
        ...     {4: {"checkpoint": "stage1_k4.pt", "vocab": "k4_vocab.json"},
        ...      6: {"checkpoint": "stage1_k6.pt", "vocab": "k6_vocab.json"}},
        ... )
        >>> result = predictor("ATGCATGCATGC...")
        >>> print(result["probability"], result["label"])
    """
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    if device is None:
        device = _select_device()

    LOGGER.info("Using device: %s", device)

    models: List[torch.nn.Module] = []
    vocabs: Dict[int, KmerVocabulary] = {}
    checkpoint_thresholds: List[float] = []

    for k, paths in checkpoints_and_vocabs.items():
        k = int(k)
        checkpoint_path = Path(paths["checkpoint"]).resolve()
        vocab_path = Path(paths["vocab"]).resolve()

        if not checkpoint_path.exists():
            LOGGER.warning("Skipping k=%d: checkpoint not found at %s", k, checkpoint_path)
            continue
        if not vocab_path.exists():
            LOGGER.warning("Skipping k=%d: vocab not found at %s", k, vocab_path)
            continue

        vocab = _load_vocab(vocab_path)
        vocabs[k] = vocab

        # Infer architecture from checkpoint to handle different training configs
        arch_overrides = _infer_architecture_from_checkpoint(checkpoint_path, device)

        model = _build_model(cfg, vocab, device, arch_overrides=arch_overrides)
        loaded_threshold = _load_checkpoint(model, checkpoint_path, device)
        if loaded_threshold is not None:
            checkpoint_thresholds.append(loaded_threshold)
        models.append(model)

        LOGGER.info("Loaded k=%d model from %s", k, checkpoint_path)

    if not models:
        raise RuntimeError("No models could be loaded for ensemble inference")

    ensemble = MultiScaleEnsemble(models).to(device)
    ensemble.eval()

    LOGGER.info("Ensemble ready with %d models", len(models))

    resolved_threshold = _resolve_decision_threshold(threshold, checkpoint_thresholds)

    return Stage1Predictor(
        ensemble=ensemble,
        vocabs=vocabs,
        cfg=cfg,
        device=device,
        threshold=resolved_threshold,
    )


def main() -> None:
    """CLI entry point for Stage 1 ensemble prediction."""
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Predict promoter probability using Stage 1 ensemble"
    )
    parser.add_argument(
        "--config",
        default=str(_TRAINING_DIR / "config.yaml"),
        help="Path to config.yaml",
    )
    parser.add_argument(
        "--sequence",
        required=True,
        help="DNA sequence to classify (A/C/G/T)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Promoter decision threshold override (default: checkpoint best_threshold or 0.5)",
    )
    parser.add_argument(
        "--kmer-configs",
        type=str,
        help="JSON string mapping k to checkpoint/vocab paths",
    )
    args = parser.parse_args()

    # Default k-mer configurations (k=4 and k=6 ensemble)
    if args.kmer_configs:
        kmer_config = json.loads(args.kmer_configs)
    else:
        artifacts_dir = _PROJECT_ROOT / "artifacts"
        kmer_config = {
            4: {
                "checkpoint": str(artifacts_dir / "checkpoints" / "stage1_k4.pt"),
                "vocab": str(artifacts_dir / "vocabs" / "k4_vocab.json"),
            },
            6: {
                "checkpoint": str(artifacts_dir / "checkpoints" / "stage1_k6.pt"),
                "vocab": str(artifacts_dir / "vocabs" / "k6_vocab.json"),
            },
        }

    predictor = load_stage1_ensemble(
        config_path=args.config,
        checkpoints_and_vocabs=kmer_config,
        threshold=args.threshold,
    )

    result = predictor.predict(args.sequence)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
