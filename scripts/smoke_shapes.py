#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import torch
import yaml


def _resolve_config_path(training_dir: Path) -> Path:
    config_path = training_dir / "config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found at {config_path}")
    return config_path


def _load_cfg(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    training_dir = repo_root / "gene_whisperer" / "training"
    if not training_dir.exists():
        raise FileNotFoundError(f"Training dir not found: {training_dir}")

    sys.path.insert(0, str(training_dir))
    from dataset import KmerVocabulary, PromoterDatasetStage1  # noqa: E402
    from model import GeneWhispererStage1Legacy  # noqa: E402

    config_path = _resolve_config_path(training_dir)
    cfg = _load_cfg(config_path)

    stage1_train = cfg.get("stage1_train")
    if not stage1_train:
        raise ValueError("stage1_train is not configured in config.yaml")
    stage1_train_path = Path(stage1_train)
    if not stage1_train_path.is_absolute():
        stage1_train_path = (training_dir / stage1_train_path).resolve()
    if not stage1_train_path.exists():
        raise FileNotFoundError(f"Stage1 train file not found: {stage1_train_path}")

    delimiter = cfg.get("delimiter", "\t")
    df = pd.read_csv(stage1_train_path, delimiter=delimiter)
    if "sequence" not in df:
        raise ValueError("stage1_train.tsv missing 'sequence' column")
    if len(df) < 2:
        raise ValueError(f"Need at least 2 rows in {stage1_train_path}, found {len(df)}")
    df = df.head(2).copy()
    if "is_promoter" not in df:
        df["is_promoter"] = 0.0

    max_bp_len = int(cfg.get("max_bp_len", 81))
    kmer = int(cfg.get("kmer", 3))
    vocab = KmerVocabulary.build_from_sequences(df["sequence"].astype(str).tolist(), k=kmer)

    engineered_dim = int(cfg.get("engineered_dim", 288))
    dataset = PromoterDatasetStage1(
        df,
        max_bp_len=max_bp_len,
        vocab=vocab,
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        feature_enable_tnc=bool(cfg.get("stage1_feature_enable_tnc", True)),
        feature_enable_pseeiip=bool(cfg.get("stage1_feature_enable_pseeiip", True)),
        engineered_dim=engineered_dim,
        reverse_complement_prob=0.0,
    )

    tokens_list = []
    engineered_list = []
    for idx in range(2):
        tokens, engineered, _ = dataset[idx]
        tokens_list.append(tokens)
        engineered_list.append(engineered)

    tokens_batch = torch.stack(tokens_list, dim=0)
    engineered_batch = torch.stack(engineered_list, dim=0)

    model = GeneWhispererStage1Legacy(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=int(cfg.get("embedding_dim", 256)),
        num_layers=int(cfg.get("transformer_layers", 6)),
        num_heads=int(cfg.get("transformer_heads", 8)),
        ff_dim=int(cfg.get("transformer_ff_dim", 1024)),
        dropout=float(cfg.get("transformer_dropout", 0.15)),
        pad_token_id=vocab.pad_id,
        engineered_dim=engineered_dim,
        use_engineered_features=bool(cfg.get("stage1_use_engineered_features", True)),
        use_attention_pool=bool(cfg.get("use_attention_pool", True)),
        use_tcn=bool(cfg.get("use_tcn", True)),
        tcn_hidden=int(cfg.get("tcn_hidden", 256)),
        tcn_levels=int(cfg.get("tcn_levels", 4)),
        tcn_kernel=int(cfg.get("tcn_kernel", 3)),
        multiscale_channels=int(cfg.get("multiscale_channels", 64)),
        multiscale_kernels=tuple(cfg.get("multiscale_kernels", [3, 5, 7, 9, 15])),
        lstm_hidden=int(cfg.get("lstm_hidden", 192)),
        post_cnn_transformer_layers=int(cfg.get("post_cnn_transformer_layers", 3)),
        engineered_mlp_hidden=int(cfg.get("engineered_mlp_hidden", 256)),
        engineered_mlp_output=int(cfg.get("engineered_mlp_output", 128)),
        fusion_hidden=int(cfg.get("fusion_hidden", 256)),
    )

    model.eval()
    with torch.no_grad():
        logits = model(tokens_batch, engineered_batch)

    assert logits.shape[0] == 2, f"Expected 2 logits, got {logits.shape}"
    print(f"logits_shape={tuple(logits.shape)} engineered_dim={engineered_batch.size(-1)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
