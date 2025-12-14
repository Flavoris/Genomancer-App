"""Masked language model pretraining on E. coli MG1655 genome."""
from __future__ import annotations

import argparse
import logging
import math
import os
import random
from pathlib import Path
from typing import Iterable, List, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from torch.optim import AdamW
from torch.utils.data import DataLoader, Dataset

from dataset import KmerVocabulary
from model import DNAEncoder

LOGGER = logging.getLogger("gene_whisperer.pretrain_mlm")
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

ALLOWED_BASES = {"A", "C", "G", "T"}


def set_seed(seed: int) -> None:
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    if hasattr(torch, "mps") and torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)


def read_fasta_sequence(path: Path) -> str:
    sequence_parts: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)


def extract_windows(sequence: str, window_size: int, stride: int) -> List[str]:
    windows: List[str] = []
    max_start = len(sequence) - window_size + 1
    for start in range(0, max_start, stride):
        window = sequence[start : start + window_size]
        if not set(window) <= ALLOWED_BASES:
            continue
        windows.append(window)
    return windows


class MLMDataset(Dataset):
    def __init__(self, token_tensors: List[torch.LongTensor]):
        self.samples = token_tensors

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> torch.LongTensor:
        return self.samples[idx]


def mask_tokens_span(
    inputs: torch.LongTensor,
    vocab: KmerVocabulary,
    mask_prob: float = 0.15,
    max_span_len: int = 3,
) -> Tuple[torch.LongTensor, torch.LongTensor]:
    """
    DNABERT-style span masking.
    
    Instead of independent per-token masking, this selects contiguous spans
    until approximately mask_prob of tokens are covered, then applies the
    standard 80/10/10 BERT masking within those spans.
    
    Args:
        inputs: (B, L) input token indices
        vocab: KmerVocabulary with mask_id and pad_id
        mask_prob: Target fraction of tokens to mask
        max_span_len: Maximum length of each masked span
        
    Returns:
        (masked_inputs, labels) tuple where labels has -100 for unmasked positions
    """
    device = inputs.device
    labels = inputs.clone()
    batch_size, seq_len = inputs.shape
    vocab_size = len(vocab.itos)

    # Initialize all positions as unmasked
    masked = torch.zeros_like(inputs, dtype=torch.bool, device=device)

    for b in range(batch_size):
        num_to_mask = max(1, int(mask_prob * seq_len))
        covered = 0
        attempts = 0
        max_attempts = num_to_mask * 10  # Safety limit to prevent infinite loops
        
        while covered < num_to_mask and attempts < max_attempts:
            span_len = random.randint(1, max_span_len)
            start = random.randint(0, max(0, seq_len - span_len))
            end = min(seq_len, start + span_len)
            
            # Avoid double-masking: only count newly masked positions
            span_mask = ~masked[b, start:end]
            if span_mask.any():
                masked[b, start:end] = True
                covered += span_mask.sum().item()
            
            attempts += 1
            
            # Safety break if we've covered enough
            if covered >= num_to_mask or covered >= seq_len:
                break

    # Don't mask padding positions
    if hasattr(vocab, "pad_id"):
        pad_mask = inputs.eq(vocab.pad_id)
        masked[pad_mask] = False

    # Set labels: -100 for unmasked positions (ignored in loss)
    labels[~masked] = -100

    # Apply 80/10/10 masking strategy within masked spans
    # 80% -> [MASK] token
    # 10% -> random token
    # 10% -> unchanged (keep original)
    rand = torch.rand_like(inputs, dtype=torch.float32, device=device)
    mask_token_indices = masked & (rand < 0.8)
    random_token_indices = masked & (rand >= 0.8) & (rand < 0.9)
    # Remaining 10% (rand >= 0.9) stay unchanged

    inputs = inputs.clone()
    inputs[mask_token_indices] = vocab.mask_id

    if random_token_indices.any():
        flat = random_token_indices.view(-1)
        num_rand = flat.sum().item()
        random_ids = torch.randint(
            low=0,
            high=vocab_size,
            size=(num_rand,),
            device=device,
            dtype=inputs.dtype,
        )
        inputs.view(-1)[flat] = random_ids

    return inputs, labels


def collate_mlm(
    batch: Iterable[torch.LongTensor],
    vocab: KmerVocabulary,
) -> Tuple[torch.LongTensor, torch.LongTensor]:
    """Collate batch and apply DNABERT-style span masking."""
    inputs = torch.stack(list(batch))
    masked, labels = mask_tokens_span(inputs, vocab)
    return masked, labels


class DNAMLM(nn.Module):
    """Masked Language Model wrapper around DNAEncoder."""
    
    def __init__(self, encoder: DNAEncoder, vocab_size: int):
        super().__init__()
        self.encoder = encoder
        self.lm_head = nn.Linear(encoder.embedding_dim, vocab_size, bias=False)

    def forward(
        self,
        token_ids: torch.LongTensor,
        labels: torch.LongTensor | None = None,
        key_padding_mask: torch.Tensor | None = None,
    ):
        """
        Forward pass for MLM.
        
        Args:
            token_ids: (B, L) input token indices
            labels: (B, L) target labels, -100 where we ignore
            key_padding_mask: Optional (B, L) boolean mask for padding
            
        Returns:
            (logits, loss) tuple where loss is None if labels not provided
        """
        x = self.encoder(token_ids, key_padding_mask=key_padding_mask)
        logits = self.lm_head(x)  # (B, L, vocab_size)
        
        if labels is None:
            return logits, None
        
        loss = F.cross_entropy(
            logits.view(-1, logits.size(-1)),
            labels.view(-1),
            ignore_index=-100,
        )
        return logits, loss


def load_or_build_vocab(windows: List[str], k: int, vocab_path: Path) -> KmerVocabulary:
    if vocab_path.exists():
        LOGGER.info("Loading existing vocabulary from %s", vocab_path)
        return KmerVocabulary.load(vocab_path)
    vocab = KmerVocabulary.build_from_sequences(windows, k=k)
    vocab.save(vocab_path)
    LOGGER.info("Saved vocabulary with %d entries to %s", len(vocab.itos), vocab_path)
    return vocab


def prepare_dataset(cfg) -> Tuple[List[torch.LongTensor], KmerVocabulary]:
    fasta_path = Path(cfg.get("mlm_fasta_path")).resolve()
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found at {fasta_path}")
    window_size = int(cfg.get("mlm_window_size", 81))
    stride = int(cfg.get("mlm_stride", 20))
    k = int(cfg.get("mlm_kmer", 3))
    vocab_path = Path(cfg.get("mlm_vocab_path", f"../artifacts/vocabs/k{k}_mlm_vocab.json")).resolve()

    LOGGER.info("Reading genome from %s", fasta_path)
    genome_sequence = read_fasta_sequence(fasta_path)
    windows = extract_windows(genome_sequence, window_size, stride)
    if not windows:
        raise ValueError("No valid windows generated from genome")
    LOGGER.info("Extracted %d windows", len(windows))
    vocab = load_or_build_vocab(windows, k=k, vocab_path=vocab_path)
    token_tensors = [vocab.tokenize(window, window_size) for window in windows]
    return token_tensors, vocab


def select_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main() -> None:
    parser = argparse.ArgumentParser(description="Pretrain DNA MLM on E. coli genome")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = (script_dir / config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    seed = int(cfg.get("seed", 1337))
    set_seed(seed)

    token_tensors, vocab = prepare_dataset(cfg)
    dataset = MLMDataset(token_tensors)
    batch_size = int(cfg.get("mlm_batch_size", 256))
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=lambda batch: collate_mlm(batch, vocab),
        num_workers=int(cfg.get("num_workers", 0)),
        pin_memory=True,
    )

    device = select_device()
    LOGGER.info("Using device: %s", device)
    embedding_dim = int(cfg.get("embedding_dim", 192))
    transformer_layers = int(cfg.get("transformer_layers", 6))
    transformer_heads = int(cfg.get("transformer_heads", 8))
    transformer_ff_dim = int(cfg.get("transformer_ff_dim", 384))
    transformer_dropout = float(cfg.get("transformer_dropout", 0.15))

    encoder = DNAEncoder(
        vocab_size=len(vocab.itos),
        kmer=vocab.k,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=transformer_ff_dim,
        dropout=transformer_dropout,
        use_alibi=bool(cfg.get("use_alibi", True)),
        pad_token_id=vocab.pad_id,
    )
    model = DNAMLM(encoder, vocab_size=len(vocab.itos)).to(device)

    optimizer = AdamW(model.parameters(), lr=float(cfg.get("mlm_lr", 1e-4)), weight_decay=float(cfg.get("mlm_weight_decay", 0.01)))
    epochs = int(cfg.get("mlm_epochs", 15))

    autocast_device = None
    if device.type == "cuda":
        autocast_device = "cuda"
    elif device.type == "mps":
        autocast_device = "mps"

    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0
        total_samples = 0
        for batch_idx, (inputs, labels) in enumerate(dataloader, start=1):
            inputs = inputs.to(device)
            labels = labels.to(device)
            optimizer.zero_grad()
            if autocast_device:
                with torch.autocast(autocast_device, dtype=torch.float16):
                    _, loss = model(inputs, labels)
            else:
                _, loss = model(inputs, labels)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += loss.item() * inputs.size(0)
            total_samples += inputs.size(0)
            if batch_idx % 50 == 0:
                LOGGER.info("Epoch %d Batch %d Loss %.4f", epoch, batch_idx, loss.item())
        avg_loss = total_loss / max(1, total_samples)
        LOGGER.info("Epoch %d average loss %.4f", epoch, avg_loss)

    encoder_path = Path(cfg.get("mlm_encoder_path", f"../artifacts/mlm_encoder_k{vocab.k}.pt")).resolve()
    encoder_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.encoder.state_dict(), encoder_path)
    LOGGER.info("Saved pretrained encoder to %s", encoder_path)


if __name__ == "__main__":
    main()
