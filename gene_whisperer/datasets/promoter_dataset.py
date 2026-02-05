"""Dataset for promoter classification and strength prediction."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import List, Sequence

import torch
from torch.utils.data import Dataset

from gene_whisperer.features.engineered import FeatureExtractor
from gene_whisperer.tokenization.bpe import BPETokenizer


def _load_tsv(path: Path) -> List[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


class PromoterDataset(Dataset[dict[str, torch.Tensor]]):
    """Promoter dataset for stage1 or stage2 classification."""

    def __init__(
        self,
        path: Path,
        tokenizer: BPETokenizer,
        feature_extractor: FeatureExtractor,
        max_length: int = 128,
        task: str = "stage1",
    ) -> None:
        self.rows = _load_tsv(path)
        if not self.rows:
            raise ValueError(f"No rows found in {path}")
        self.tokenizer = tokenizer
        self.feature_extractor = feature_extractor
        self.max_length = max_length
        self.task = task

        if task not in ("stage1", "stage2"):
            raise ValueError("task must be stage1 or stage2")

    def __len__(self) -> int:
        return len(self.rows)

    def _label_from_row(self, row: dict[str, str]) -> int:
        if self.task == "stage1":
            return int(row.get("is_promoter", 0))
        strength = row.get("strength", "weak").strip().lower()
        return 1 if strength == "strong" else 0

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        row = self.rows[idx]
        sequence = row.get("sequence", "")
        token_ids, attention_mask = self.tokenizer.encode(
            sequence, add_special_tokens=True, max_length=self.max_length, pad_to_max=True
        )
        engineered = self.feature_extractor.transform(sequence)
        label = self._label_from_row(row)

        return {
            "input_ids": torch.tensor(token_ids, dtype=torch.long),
            "attention_mask": torch.tensor(attention_mask, dtype=torch.long),
            "engineered": torch.tensor(engineered, dtype=torch.float32),
            "label": torch.tensor(label, dtype=torch.float32),
        }
