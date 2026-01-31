"""Stage 1 classifier models for Gene Whisperer."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple, Union

import torch
import torch.nn as nn

from .encoder import DNAEncoder, AttentionPooling, _load_checkpoint_file
from .features import EngineeredFeatureMLP

LOGGER = logging.getLogger(__name__)


def _mean_pool_sequence(
    token_features: torch.Tensor,
    padding_mask: Optional[torch.Tensor],
) -> torch.Tensor:
    """Mean-pool sequence features with optional padding mask."""
    if padding_mask is None:
        return token_features.mean(dim=1)
    keep_mask = (~padding_mask).unsqueeze(-1).type_as(token_features)
    summed = (token_features * keep_mask).sum(dim=1)
    denom = keep_mask.sum(dim=1).clamp(min=1.0)
    return summed / denom


class GeneWhispererStage1(nn.Module):
    """Stage 1 promoter classifier with simplified transformer-only architecture.

    Architecture:
        Embedding -> 12-Layer Transformer -> Attention Pool -> Classifier

    Optionally concatenates engineered features (TNC, PseEIIP, CKSaap, PSTNPss)
    before the final classifier for improved accuracy.
    """

    def __init__(
        self,
        vocab_size: int = 4096,
        embedding_dim: int = 192,
        num_layers: int = 12,
        num_heads: int = 8,
        ff_dim: int = 384,
        dropout: float = 0.15,
        pad_token_id: Optional[int] = 0,
        encoder: Optional[DNAEncoder] = None,
        engineered_dim: int = 288,
        use_engineered_features: bool = True,
        use_attention_pool: bool = True,
        pooling_type: Optional[str] = None,
        engineered_mlp_hidden: int = 512,
        engineered_mlp_output: int = 256,
        engineered_mlp_pre_norm: bool = True,
        engineered_mlp_input_dropout: float = 0.05,
        engineered_mlp_use_gated: bool = True,
        engineered_mlp_use_residual: bool = True,
        classifier_hidden: Optional[int] = None,
        classifier_dropout: Optional[float] = None,
        drop_path_rate: float = 0.1,
        max_seq_len: int = 24,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
        use_rope: bool = False,
        rope_base: float = 10000.0,
        ffn_type: str = None,
        norm_type: str = "layernorm",
        ffn_mult: float = None,
        use_positional_motif_bias: bool = False,
        motif_regions: Optional[dict] = None,
    ):
        super().__init__()

        if pooling_type is not None:
            normalized_pooling = pooling_type.strip().lower()
            if normalized_pooling == "attention":
                use_attention_pool = True
            elif normalized_pooling == "mean":
                use_attention_pool = False
            else:
                raise ValueError(f"Unsupported pooling_type: {pooling_type}")

        self.pad_token_id = pad_token_id
        self.use_attention_pool = use_attention_pool
        self.engineered_dim = engineered_dim
        self.use_engineered_features = use_engineered_features and engineered_dim > 0

        self._full_encoder = encoder or DNAEncoder(
            vocab_size=vocab_size,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
            use_rope=use_rope,
            rope_base=rope_base,
            ffn_type=ffn_type,
            norm_type=norm_type,
            ffn_mult=ffn_mult,
            use_positional_motif_bias=use_positional_motif_bias,
            motif_regions=motif_regions,
        )
        self.embedding = self._full_encoder.embedding
        encoder_dim = self._full_encoder.embedding_dim

        if use_attention_pool:
            if encoder_dim % 4 != 0:
                raise ValueError(
                    f"embedding_dim ({encoder_dim}) must be divisible by 4 for attention pooling"
                )
            self.pool = AttentionPooling(encoder_dim, num_heads=4)
        else:
            self.pool = None

        if self.use_engineered_features:
            self.engineered_mlp = EngineeredFeatureMLP(
                input_dim=engineered_dim,
                hidden_dim=engineered_mlp_hidden,
                output_dim=engineered_mlp_output,
                dropout=dropout,
                use_pre_norm=engineered_mlp_pre_norm,
                input_dropout=engineered_mlp_input_dropout,
                use_gated=engineered_mlp_use_gated,
                use_residual=engineered_mlp_use_residual,
            )
            classifier_in = encoder_dim + self.engineered_mlp.output_dim
        else:
            self.engineered_mlp = None
            classifier_in = encoder_dim

        hidden_dim = classifier_hidden or encoder_dim
        classifier_drop = dropout if classifier_dropout is None else classifier_dropout
        self.classifier = nn.Sequential(
            nn.Linear(classifier_in, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(classifier_drop),
            nn.Linear(hidden_dim, 1),
        )

    @property
    def encoder(self) -> DNAEncoder:
        """Compatibility property for MLM weight loading."""
        return self._full_encoder

    def load_pretrained_weights(
        self,
        checkpoint_path: Union[str, Path],
        strict: bool = False,
        transfer_mode: str = "embed_only",
    ) -> None:
        """Load pretrained MLM weights into the transformer encoder."""
        if transfer_mode == "none":
            LOGGER.info("Skipping MLM weight loading (transfer_mode=none)")
            return
        if transfer_mode not in {"embed_only", "embed_plus_adapter"}:
            raise ValueError(f"Unsupported transfer_mode: {transfer_mode}")

        self._full_encoder.load_mlm_weights(checkpoint_path, strict=strict)
        if transfer_mode == "embed_plus_adapter":
            LOGGER.info("Simplified model ignores adapter transfer; loaded encoder weights only.")

    def load_legacy_checkpoint(
        self,
        checkpoint: Union[str, Path, dict],
        strict: bool = False,
    ) -> int:
        """Load legacy Stage 1 checkpoints while ignoring incompatible weights."""
        if isinstance(checkpoint, (str, Path)):
            state_dict = _load_checkpoint_file(Path(checkpoint))
        else:
            state_dict = checkpoint

        if "model_state_dict" in state_dict:
            state_dict = state_dict["model_state_dict"]
        elif "model_state" in state_dict:
            state_dict = state_dict["model_state"]
        elif "state_dict" in state_dict:
            state_dict = state_dict["state_dict"]

        LOGGER.warning(
            "Loading legacy Stage 1 checkpoint into simplified architecture; "
            "CNN/TCN and adapter weights will be ignored."
        )

        model_state = self.state_dict()
        compatible: dict = {}
        for key, value in state_dict.items():
            target = model_state.get(key)
            if target is None:
                continue
            if getattr(target, "shape", None) != getattr(value, "shape", None):
                continue
            compatible[key] = value

        missing = set(model_state.keys()) - set(compatible.keys())
        skipped = set(state_dict.keys()) - set(compatible.keys())

        if strict and missing:
            raise RuntimeError(
                f"strict=True but {len(missing)} model parameters not in checkpoint."
            )

        self.load_state_dict(compatible, strict=False)
        LOGGER.info(
            "Loaded %d compatible tensors (missing=%d, skipped=%d)",
            len(compatible),
            len(missing),
            len(skipped),
        )
        return len(compatible)

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        """Generate padding mask from tokens."""
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def _encode_tokens(
        self,
        tokens: torch.Tensor,
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        """Encode tokens with the full transformer and return token features + mask."""
        seq_len = tokens.size(1)
        max_len = self._full_encoder.pos_embedding.num_embeddings
        if seq_len > max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds max_seq_len={max_len}. "
                "Increase max_seq_len or shorten inputs."
            )
        padding_mask = self._get_padding_mask(tokens)
        token_features = self._full_encoder(tokens, key_padding_mask=padding_mask)
        return token_features, padding_mask

    def _mean_pool(
        self,
        token_features: torch.Tensor,
        padding_mask: Optional[torch.Tensor],
    ) -> torch.Tensor:
        if padding_mask is None:
            return token_features.mean(dim=1)
        keep_mask = (~padding_mask).unsqueeze(-1).type_as(token_features)
        summed = (token_features * keep_mask).sum(dim=1)
        denom = keep_mask.sum(dim=1).clamp(min=1.0)
        return summed / denom

    def _pool_sequence(
        self,
        token_features: torch.Tensor,
        padding_mask: Optional[torch.Tensor],
    ) -> torch.Tensor:
        if self.use_attention_pool:
            return self.pool(token_features, padding_mask)
        return self._mean_pool(token_features, padding_mask)

    def extract_features(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Extract pooled sequence features and optional fused features."""
        token_features, padding_mask = self._encode_tokens(tokens)
        seq_features = self._pool_sequence(token_features, padding_mask)

        if engineered_features is not None:
            assert engineered_features.size(-1) == self.engineered_dim, (
                f"Engineered feature dim mismatch: got {engineered_features.size(-1)}, "
                f"expected {self.engineered_dim}"
            )

        if self.use_engineered_features and engineered_features is not None:
            eng_features = self.engineered_mlp(engineered_features)
            fused = torch.cat([seq_features, eng_features], dim=-1)
            return seq_features, fused

        return seq_features, seq_features

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass for simplified Stage 1 classification."""
        _, fused = self.extract_features(tokens, engineered_features)
        logits = self.classifier(fused)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits


class GeneWhispererTransformerOnly(nn.Module):
    """Variant A: Transformer-only Stage 1 model with mean pooling."""

    def __init__(
        self,
        vocab_size: int = 4096,
        embedding_dim: int = 384,
        num_layers: int = 12,
        num_heads: int = 12,
        ff_dim: int = 1536,
        dropout: float = 0.1,
        pad_token_id: Optional[int] = 0,
        max_seq_len: int = 256,
        classifier_hidden: Optional[int] = None,
        classifier_dropout: Optional[float] = None,
        drop_path_rate: float = 0.1,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
        use_rope: bool = False,
        rope_base: float = 10000.0,
        ffn_type: str = None,
        norm_type: str = "layernorm",
        ffn_mult: float = None,
        use_positional_motif_bias: bool = False,
        motif_regions: Optional[dict] = None,
    ):
        super().__init__()
        self.pad_token_id = pad_token_id

        self.encoder = DNAEncoder(
            vocab_size=vocab_size,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
            use_rope=use_rope,
            rope_base=rope_base,
            ffn_type=ffn_type,
            norm_type=norm_type,
            ffn_mult=ffn_mult,
            use_positional_motif_bias=use_positional_motif_bias,
            motif_regions=motif_regions,
        )
        self.embedding = self.encoder.embedding

        hidden_dim = classifier_hidden or embedding_dim
        classifier_drop = dropout if classifier_dropout is None else classifier_dropout
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(classifier_drop),
            nn.Linear(hidden_dim, 1),
        )

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: Optional[torch.Tensor] = None,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        padding_mask = self._get_padding_mask(tokens)
        token_features = self.encoder(tokens, key_padding_mask=padding_mask)
        pooled = _mean_pool_sequence(token_features, padding_mask)
        logits = self.classifier(pooled)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits


class GeneWhispererFeaturesOnly(nn.Module):
    """Variant B: Features-only Stage 1 model without transformer."""

    def __init__(
        self,
        engineered_dim: int = 288,
        hidden_dim: int = 512,
        second_hidden_dim: int = 256,
        dropout: float = 0.3,
    ):
        super().__init__()
        self.engineered_dim = engineered_dim

        self.mlp = nn.Sequential(
            nn.Linear(engineered_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, second_hidden_dim),
            nn.LayerNorm(second_hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(second_hidden_dim, 1),
        )

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: torch.Tensor,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        if engineered_features is None:
            raise ValueError("FeaturesOnly model requires engineered_features")

        logits = self.mlp(engineered_features)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits


class GeneWhispererCombined(nn.Module):
    """Variant C: Combined transformer + engineered features model."""

    def __init__(
        self,
        vocab_size: int = 4096,
        embedding_dim: int = 384,
        num_layers: int = 12,
        num_heads: int = 12,
        ff_dim: int = 1536,
        dropout: float = 0.1,
        pad_token_id: Optional[int] = 0,
        max_seq_len: int = 256,
        engineered_dim: int = 288,
        engineered_mlp_hidden: int = 256,
        engineered_mlp_output: int = 128,
        classifier_hidden: Optional[int] = None,
        classifier_dropout: Optional[float] = None,
        drop_path_rate: float = 0.1,
        use_relative_position_bias: bool = False,
        relative_position_num_buckets: int = 32,
        relative_position_max_distance: int = 128,
        use_glu_ffn: bool = False,
        glu_activation: str = "gelu",
        use_rope: bool = False,
        rope_base: float = 10000.0,
        ffn_type: str = None,
        norm_type: str = "layernorm",
        ffn_mult: float = None,
        use_positional_motif_bias: bool = False,
        motif_regions: Optional[dict] = None,
    ):
        super().__init__()
        self.pad_token_id = pad_token_id
        self.engineered_dim = engineered_dim

        self.encoder = DNAEncoder(
            vocab_size=vocab_size,
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout=dropout,
            pad_token_id=pad_token_id,
            max_seq_len=max_seq_len,
            drop_path_rate=drop_path_rate,
            use_relative_position_bias=use_relative_position_bias,
            relative_position_num_buckets=relative_position_num_buckets,
            relative_position_max_distance=relative_position_max_distance,
            use_glu_ffn=use_glu_ffn,
            glu_activation=glu_activation,
            use_rope=use_rope,
            rope_base=rope_base,
            ffn_type=ffn_type,
            norm_type=norm_type,
            ffn_mult=ffn_mult,
            use_positional_motif_bias=use_positional_motif_bias,
            motif_regions=motif_regions,
        )
        self.embedding = self.encoder.embedding

        self.engineered_mlp = EngineeredFeatureMLP(
            input_dim=engineered_dim,
            hidden_dim=engineered_mlp_hidden,
            output_dim=engineered_mlp_output,
            dropout=dropout,
        )

        combined_dim = embedding_dim + engineered_mlp_output
        hidden_dim = classifier_hidden or combined_dim
        classifier_drop = dropout if classifier_dropout is None else classifier_dropout

        self.classifier = nn.Sequential(
            nn.Linear(combined_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(classifier_drop),
            nn.Linear(hidden_dim, 1),
        )

    def _get_padding_mask(self, tokens: torch.Tensor) -> Optional[torch.Tensor]:
        if self.pad_token_id is not None:
            return tokens.eq(self.pad_token_id)
        return None

    def forward(
        self,
        tokens: torch.Tensor,
        engineered_features: torch.Tensor,
        return_logits: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        if engineered_features is None:
            raise ValueError("Combined model requires engineered_features")

        padding_mask = self._get_padding_mask(tokens)
        token_features = self.encoder(tokens, key_padding_mask=padding_mask)
        seq_features = _mean_pool_sequence(token_features, padding_mask)

        eng_features = self.engineered_mlp(engineered_features)
        combined = torch.cat([seq_features, eng_features], dim=-1)

        logits = self.classifier(combined)

        if return_logits:
            probs = torch.sigmoid(logits)
            return probs, logits

        return logits
