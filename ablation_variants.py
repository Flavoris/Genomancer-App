from __future__ import annotations

from typing import Any

Variant = tuple[str, dict[str, Any]]


STAGE1_VARIANTS: list[Variant] = [
    ("no_attention_pool", {"use_attention_pool": False}),
    ("no_tcn", {"use_tcn": False}),
    ("no_postcnn_transformer", {"post_cnn_transformer_layers": 0}),
    ("no_engineered_features", {"stage1_use_engineered_features": False}),
    ("no_tnc", {"stage1_feature_enable_tnc": False}),
    ("no_pseeiip", {"stage1_feature_enable_pseeiip": False}),
    ("no_cksnap", {"stage1_feature_enable_cksnap": False}),
    ("no_pstnp", {"stage1_feature_enable_pstnp": False}),
]

STAGE2_VARIANTS: list[Variant] = [
    ("no_attention_pool", {"use_attention_pool": False}),
    ("no_tcn", {"use_tcn": False}),
    ("no_postcnn_transformer", {"post_cnn_transformer_layers": 0}),
    ("no_engineered_features", {"stage2_use_engineered_features": False}),
    ("no_tnc", {"stage2_feature_enable_tnc": False}),
    ("no_pseeiip", {"stage2_feature_enable_pseeiip": False}),
    ("no_cksnap", {"stage2_feature_enable_cksnap": False}),
    ("no_pstnp", {"stage2_feature_enable_pstnp": False}),
]

PIPELINE_VARIANTS: list[Variant] = [
    ("baseline", {}),
    ("no_attention_pool", {"use_attention_pool": False}),
    ("no_tcn", {"use_tcn": False}),
    ("no_postcnn_transformer", {"post_cnn_transformer_layers": 0}),
    (
        "no_engineered_features",
        {"stage1_use_engineered_features": False, "stage2_use_engineered_features": False},
    ),
    (
        "no_tnc",
        {"stage1_feature_enable_tnc": False, "stage2_feature_enable_tnc": False},
    ),
    (
        "no_pseeiip",
        {"stage1_feature_enable_pseeiip": False, "stage2_feature_enable_pseeiip": False},
    ),
    (
        "no_cksnap",
        {"stage1_feature_enable_cksnap": False, "stage2_feature_enable_cksnap": False},
    ),
    (
        "no_pstnp",
        {"stage1_feature_enable_pstnp": False, "stage2_feature_enable_pstnp": False},
    ),
]
