from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from dataset import (
    KmerVocabulary,
    PromoterDatasetStage1,
    build_dataloaders,
    compute_pseeiip,
    compute_tnc,
)
from ensemble_infer import compute_engineered_features


def _manual_engineered_features(sequence: str) -> torch.Tensor:
    tnc = compute_tnc(sequence)
    pseeiip = compute_pseeiip(sequence)
    return torch.cat([tnc, pseeiip], dim=0)


def test_stage1_feature_group_zeroing_keeps_dim_and_segments() -> None:
    sequence = "ACGT" * 30
    vocab = KmerVocabulary.build_from_sequences([sequence], k=3)
    df = pd.DataFrame({"sequence": [sequence], "is_promoter": [1.0]})

    base_ds = PromoterDatasetStage1(
        df,
        max_bp_len=81,
        vocab=vocab,
        use_engineered_features=True,
        reverse_complement_prob=0.0,
    )
    _, base_features, _ = base_ds[0]
    manual = _manual_engineered_features(sequence)

    assert base_features.shape == (128,)
    assert torch.allclose(base_features, manual)

    tnc_off = PromoterDatasetStage1(
        df,
        max_bp_len=81,
        vocab=vocab,
        use_engineered_features=True,
        feature_enable_tnc=False,
        reverse_complement_prob=0.0,
    )
    _, tnc_off_features, _ = tnc_off[0]
    assert tnc_off_features.shape == (128,)
    assert tnc_off_features[:64].abs().sum().item() == 0.0
    assert torch.allclose(tnc_off_features[64:], base_features[64:])

    pseeiip_off = PromoterDatasetStage1(
        df,
        max_bp_len=81,
        vocab=vocab,
        use_engineered_features=True,
        feature_enable_pseeiip=False,
        reverse_complement_prob=0.0,
    )
    _, pseeiip_off_features, _ = pseeiip_off[0]
    assert pseeiip_off_features.shape == (128,)
    assert pseeiip_off_features[64:].abs().sum().item() == 0.0
    assert torch.allclose(pseeiip_off_features[:64], base_features[:64])


def test_stage1_use_engineered_features_false_zeros_everything() -> None:
    sequence = "ACGT" * 30
    vocab = KmerVocabulary.build_from_sequences([sequence], k=3)
    df = pd.DataFrame({"sequence": [sequence], "is_promoter": [0.0]})

    ds = PromoterDatasetStage1(
        df,
        max_bp_len=81,
        vocab=vocab,
        use_engineered_features=False,
        feature_enable_tnc=True,
        feature_enable_pseeiip=True,
        reverse_complement_prob=0.0,
    )
    _, engineered, _ = ds[0]
    assert engineered.shape == (128,)
    assert engineered.abs().sum().item() == 0.0


def test_ensemble_infer_engineered_features_mirrors_cfg_toggles() -> None:
    sequence = "ACGT" * 30
    base = compute_engineered_features(sequence)
    assert base.shape == (128,)

    cfg = {"stage1_feature_enable_tnc": False}
    tnc_off = compute_engineered_features(sequence, cfg)
    assert tnc_off.shape == (128,)
    assert tnc_off[:64].abs().sum().item() == 0.0
    assert torch.allclose(tnc_off[64:], base[64:])

    cfg_pseeiip_off = {"stage1_feature_enable_pseeiip": False}
    pseeiip_off = compute_engineered_features(sequence, cfg_pseeiip_off)
    assert pseeiip_off.shape == (128,)
    assert pseeiip_off[64:].abs().sum().item() == 0.0
    assert torch.allclose(pseeiip_off[:64], base[:64])


def test_build_dataloaders_passes_feature_toggles(tmp_path: Path) -> None:
    sequence = "ACGT" * 30

    stage1_train = tmp_path / "stage1_train.tsv"
    stage1_val = tmp_path / "stage1_val.tsv"
    stage2_train = tmp_path / "stage2_train.tsv"
    stage2_val = tmp_path / "stage2_val.tsv"

    pd.DataFrame({"sequence": [sequence], "is_promoter": [1.0]}).to_csv(stage1_train, sep="\t", index=False)
    pd.DataFrame({"sequence": [sequence], "is_promoter": [0.0]}).to_csv(stage1_val, sep="\t", index=False)
    pd.DataFrame({"sequence": [sequence], "strength": [1.0]}).to_csv(stage2_train, sep="\t", index=False)
    pd.DataFrame({"sequence": [sequence], "strength": [0.0]}).to_csv(stage2_val, sep="\t", index=False)

    cfg = {
        "max_bp_len": 81,
        "batch_size": 1,
        "delimiter": "\t",
        "train_val_split": 0.9,
        "num_workers": 0,
        "pin_memory": False,
        "kmer": 3,
        "vocab_cache_dir": str(tmp_path / "vocabs"),
        "stage1_train": str(stage1_train),
        "stage1_val": str(stage1_val),
        "stage2_train": str(stage2_train),
        "stage2_val": str(stage2_val),
        "stage2_strength_positive": 1,
        "stage2_strength_negative": 0,
        "stage1_use_engineered_features": True,
        "stage1_feature_enable_tnc": False,
        "stage1_feature_enable_pseeiip": False,
        "stage1_reverse_complement_prob": 0.0,
    }

    loaders = build_dataloaders(cfg)
    train_ds = loaders["stage1"]["train"].dataset
    val_ds = loaders["stage1"]["val"].dataset

    assert train_ds.feature_enable_tnc is False
    assert train_ds.feature_enable_pseeiip is False

    assert val_ds.feature_enable_tnc is False
    assert val_ds.feature_enable_pseeiip is False
