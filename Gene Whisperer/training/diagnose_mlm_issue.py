#!/usr/bin/env python3
"""
Fast diagnostic to identify MLM training issues.

This script runs minimal experiments to identify the root cause of
the ~18% accuracy plateau:

1. Random baseline: What accuracy would random guessing achieve?
2. Overfit test: Can the model learn at all?
3. Token distribution: Is masking biased?
4. Prediction analysis: What is the model actually predicting?

Run: python diagnose_mlm_issue.py
"""

import json
import logging
import math
import sys
from pathlib import Path
from collections import Counter

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import AdamW

# Setup path
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from bpe_tokenizer import DNABPETokenizer
from model import DNAEncoder
from pretrain_components.mlm_model import DNAMLM
from pretrain_mlm import mask_tokens_span

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
LOGGER = logging.getLogger(__name__)


def load_vocab(vocab_path: Path) -> DNABPETokenizer:
    """Load BPE vocabulary."""
    vocab = DNABPETokenizer()
    with open(vocab_path, "r") as f:
        vocab_data = json.load(f)
    vocab.vocab = vocab_data["vocab"]
    # Build id_to_token from vocab (invert the mapping)
    vocab.id_to_token = {v: k for k, v in vocab_data["vocab"].items()}
    vocab.merges = [tuple(m) for m in vocab_data["merges"]]
    return vocab


def generate_random_dna(length: int) -> str:
    """Generate random DNA sequence."""
    import random
    return "".join(random.choices("ACGT", k=length))


def test_random_baseline(vocab: DNABPETokenizer, n_samples: int = 1000):
    """Test what accuracy random guessing achieves."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 1: RANDOM BASELINE")
    LOGGER.info("=" * 60)

    # Get vocabulary stats
    vocab_size = len(vocab.vocab)
    special_token_ids = {vocab.pad_id, vocab.unk_id, vocab.mask_id, 2, 3}  # PAD, UNK, MASK, CLS, SEP
    normal_tokens = [i for i in range(vocab_size) if i not in special_token_ids]
    n_normal_tokens = len(normal_tokens)

    LOGGER.info("Vocabulary size: %d", vocab_size)
    LOGGER.info("Normal (non-special) tokens: %d", n_normal_tokens)
    LOGGER.info("Special tokens: %s", sorted(special_token_ids))

    # Random baseline: 1/n_normal_tokens
    random_accuracy = 1.0 / n_normal_tokens
    LOGGER.info("Random baseline accuracy: %.4f (%.2f%%)", random_accuracy, random_accuracy * 100)

    # Generate random sequences and count token distribution
    all_tokens = []
    for _ in range(n_samples):
        dna = generate_random_dna(234)  # Window size
        tokens = vocab.tokenize(dna)
        all_tokens.extend(tokens)

    # Token frequency distribution
    token_counts = Counter(all_tokens)
    total_tokens = len(all_tokens)

    LOGGER.info("\nToken distribution (top 20 most common):")
    for token_id, count in token_counts.most_common(20):
        token_str = vocab.id_to_token.get(token_id, f"<{token_id}>")
        pct = count / total_tokens * 100
        LOGGER.info("  Token %4d (%8s): %6d (%.2f%%)", token_id, token_str, count, pct)

    # Check how skewed the distribution is
    top1_pct = token_counts.most_common(1)[0][1] / total_tokens
    top10_pct = sum(c for _, c in token_counts.most_common(10)) / total_tokens
    top100_pct = sum(c for _, c in token_counts.most_common(100)) / total_tokens

    LOGGER.info("\nDistribution concentration:")
    LOGGER.info("  Top 1 token: %.2f%%", top1_pct * 100)
    LOGGER.info("  Top 10 tokens: %.2f%%", top10_pct * 100)
    LOGGER.info("  Top 100 tokens: %.2f%%", top100_pct * 100)

    # If model always predicts the most common token
    most_common_accuracy = top1_pct
    LOGGER.info("\nIf model predicts most common token: %.2f%% accuracy", most_common_accuracy * 100)

    return {
        "vocab_size": vocab_size,
        "n_normal_tokens": n_normal_tokens,
        "random_accuracy": random_accuracy,
        "top1_concentration": top1_pct,
        "top10_concentration": top10_pct,
        "most_common_accuracy": most_common_accuracy,
    }


def test_masking(vocab: DNABPETokenizer, n_samples: int = 100):
    """Test masking behavior."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 2: MASKING ANALYSIS")
    LOGGER.info("=" * 60)

    # Generate sequences and apply masking
    masked_token_ids = []
    label_token_ids = []

    for _ in range(n_samples):
        dna = generate_random_dna(234)
        tokens = vocab.tokenize(dna)

        # Pad to max length
        max_len = 53  # mlm_max_token_len from config
        if len(tokens) < max_len:
            tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
        else:
            tokens = tokens[:max_len]

        input_ids = torch.tensor([tokens])

        # Apply masking
        masked_ids, labels = mask_tokens_span(
            input_ids,
            vocab,
            mask_prob=0.15,
            max_span_len=1,
            exclude_special_from_labels=True,
        )

        # Collect masked positions
        mask = labels[0] != -100
        if mask.sum() > 0:
            masked_token_ids.extend(masked_ids[0][mask].tolist())
            label_token_ids.extend(labels[0][mask].tolist())

    # Check what tokens are being masked
    label_counts = Counter(label_token_ids)
    masked_counts = Counter(masked_token_ids)

    LOGGER.info("Labels (ground truth) distribution (top 10):")
    for token_id, count in label_counts.most_common(10):
        token_str = vocab.id_to_token.get(token_id, f"<{token_id}>")
        LOGGER.info("  Token %4d (%8s): %6d", token_id, token_str, count)

    LOGGER.info("\nMasked input distribution (top 10):")
    for token_id, count in masked_counts.most_common(10):
        token_str = vocab.id_to_token.get(token_id, f"<{token_id}>")
        LOGGER.info("  Token %4d (%8s): %6d", token_id, token_str, count)

    # Check if [MASK] token is being used
    mask_token_count = masked_counts.get(vocab.mask_id, 0)
    total_masked = len(masked_token_ids)
    LOGGER.info("\n[MASK] token usage: %d / %d (%.2f%%)",
                mask_token_count, total_masked, mask_token_count / max(1, total_masked) * 100)

    # Check for special tokens in labels (should be none)
    special_in_labels = sum(label_counts.get(i, 0) for i in [0, 1, 2, 3, 4])
    LOGGER.info("Special tokens in labels: %d (should be 0)", special_in_labels)

    return {
        "total_masked": total_masked,
        "mask_token_usage": mask_token_count / max(1, total_masked),
        "special_in_labels": special_in_labels,
    }


def test_overfit(vocab: DNABPETokenizer, device: torch.device, n_steps: int = 200):
    """Test if the model can overfit on a tiny batch."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 3: OVERFIT TEST (Can the model learn at all?)")
    LOGGER.info("=" * 60)

    # Create a tiny model for fast testing
    embedding_dim = 128
    transformer_layers = 2
    transformer_heads = 4
    vocab_size = len(vocab.vocab)

    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=embedding_dim * 4,
        dropout=0.0,  # No dropout for overfit test
        pad_token_id=vocab.pad_id,
        use_rope=False,  # BERT-style
        ffn_type="gelu",
        norm_type="layernorm",
    )

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=vocab_size,
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,
    ).to(device)

    # Generate a small fixed batch
    batch_size = 8
    max_len = 53

    fixed_sequences = []
    for i in range(batch_size):
        # Use same seed for reproducibility
        dna = generate_random_dna(234)
        tokens = vocab.tokenize(dna)
        if len(tokens) < max_len:
            tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
        else:
            tokens = tokens[:max_len]
        fixed_sequences.append(tokens)

    input_ids = torch.tensor(fixed_sequences, device=device)

    # Apply masking once and keep it fixed
    masked_ids, labels = mask_tokens_span(
        input_ids,
        vocab,
        mask_prob=0.15,
        max_span_len=1,
        exclude_special_from_labels=True,
    )
    masked_ids = masked_ids.to(device)
    labels = labels.to(device)

    n_masked = (labels != -100).sum().item()
    LOGGER.info("Fixed batch: %d sequences, %d masked tokens", batch_size, n_masked)

    # Optimizer with high learning rate for fast overfitting
    optimizer = AdamW(model.parameters(), lr=5e-3, weight_decay=0.0)

    # Training loop
    LOGGER.info("Training for %d steps...", n_steps)

    losses = []
    accuracies = []

    for step in range(n_steps):
        optimizer.zero_grad()

        logits, loss = model(masked_ids, labels)

        loss.backward()
        optimizer.step()

        # Compute accuracy
        with torch.no_grad():
            preds = logits.argmax(dim=-1)
            mask = labels != -100
            correct = (preds[mask] == labels[mask]).sum().item()
            total = mask.sum().item()
            acc = correct / max(1, total)

        losses.append(loss.item())
        accuracies.append(acc)

        if step % 50 == 0 or step == n_steps - 1:
            LOGGER.info("  Step %3d: loss=%.4f, accuracy=%.4f (%.2f%%)",
                       step, loss.item(), acc, acc * 100)

    final_acc = accuracies[-1]
    peak_acc = max(accuracies)

    LOGGER.info("\nOverfit test results:")
    LOGGER.info("  Initial accuracy: %.4f (%.2f%%)", accuracies[0], accuracies[0] * 100)
    LOGGER.info("  Final accuracy: %.4f (%.2f%%)", final_acc, final_acc * 100)
    LOGGER.info("  Peak accuracy: %.4f (%.2f%%)", peak_acc, peak_acc * 100)

    if peak_acc > 0.90:
        LOGGER.info("  ✅ Model CAN learn (overfit successful)")
    elif peak_acc > 0.50:
        LOGGER.info("  ⚠️ Model learning but slow (check learning rate or architecture)")
    else:
        LOGGER.info("  ❌ Model CANNOT learn (architecture issue or gradient problem)")

    return {
        "initial_accuracy": accuracies[0],
        "final_accuracy": final_acc,
        "peak_accuracy": peak_acc,
        "losses": losses,
        "accuracies": accuracies,
    }


def test_prediction_analysis(vocab: DNABPETokenizer, device: torch.device):
    """Analyze what predictions the model makes."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 4: PREDICTION ANALYSIS (What is the model outputting?)")
    LOGGER.info("=" * 60)

    # Create model with same config as training
    embedding_dim = 384
    transformer_layers = 12
    transformer_heads = 12
    vocab_size = len(vocab.vocab)

    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=1536,
        dropout=0.1,
        pad_token_id=vocab.pad_id,
        use_rope=False,  # Config says mlm_use_rope: false
        ffn_type="gelu",  # Config says mlm_ffn_type: gelu
        norm_type="layernorm",  # Config says mlm_norm_type: layernorm
    )

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=vocab_size,
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,
    ).to(device)

    model.eval()

    # Generate test batch
    batch_size = 32
    max_len = 53

    sequences = []
    for _ in range(batch_size):
        dna = generate_random_dna(234)
        tokens = vocab.tokenize(dna)
        if len(tokens) < max_len:
            tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
        else:
            tokens = tokens[:max_len]
        sequences.append(tokens)

    input_ids = torch.tensor(sequences, device=device)
    masked_ids, labels = mask_tokens_span(
        input_ids,
        vocab,
        mask_prob=0.15,
        max_span_len=1,
        exclude_special_from_labels=True,
    )
    masked_ids = masked_ids.to(device)
    labels = labels.to(device)

    with torch.no_grad():
        logits, _ = model(masked_ids, labels)
        probs = F.softmax(logits, dim=-1)
        preds = logits.argmax(dim=-1)

    # Analyze predictions at masked positions
    mask = labels != -100
    pred_tokens = preds[mask].cpu().tolist()
    true_tokens = labels[mask].cpu().tolist()

    # Prediction distribution
    pred_counts = Counter(pred_tokens)
    LOGGER.info("Prediction distribution (untrained model, top 10):")
    for token_id, count in pred_counts.most_common(10):
        token_str = vocab.id_to_token.get(token_id, f"<{token_id}>")
        LOGGER.info("  Token %4d (%8s): %6d", token_id, token_str, count)

    # Check if predictions are concentrated on a few tokens
    total_preds = len(pred_tokens)
    top1_pct = pred_counts.most_common(1)[0][1] / total_preds if total_preds > 0 else 0
    LOGGER.info("\nPrediction concentration:")
    LOGGER.info("  Top 1 prediction: %.2f%%", top1_pct * 100)

    # Entropy of predictions (higher = more diverse)
    pred_probs = probs[mask]  # (N_masked, vocab_size)
    entropy = -(pred_probs * pred_probs.log()).sum(dim=-1).mean().item()
    max_entropy = math.log(vocab_size)
    LOGGER.info("  Average entropy: %.2f (max possible: %.2f)", entropy, max_entropy)

    # Check confidence distribution
    max_probs = pred_probs.max(dim=-1).values
    LOGGER.info("  Mean max probability: %.4f", max_probs.mean().item())
    LOGGER.info("  Min max probability: %.4f", max_probs.min().item())
    LOGGER.info("  Max max probability: %.4f", max_probs.max().item())

    return {
        "top1_concentration": top1_pct,
        "entropy": entropy,
        "max_entropy": max_entropy,
        "mean_max_prob": max_probs.mean().item(),
    }


def test_gradient_flow(vocab: DNABPETokenizer, device: torch.device):
    """Test if gradients are flowing correctly through the model."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 5: GRADIENT FLOW ANALYSIS")
    LOGGER.info("=" * 60)

    # Create model
    embedding_dim = 384
    transformer_layers = 12
    transformer_heads = 12
    vocab_size = len(vocab.vocab)

    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=1536,
        dropout=0.1,
        pad_token_id=vocab.pad_id,
        use_rope=False,
        ffn_type="gelu",
        norm_type="layernorm",
    )

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=vocab_size,
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,
    ).to(device)

    # Generate batch
    batch_size = 8
    max_len = 53

    sequences = []
    for _ in range(batch_size):
        dna = generate_random_dna(234)
        tokens = vocab.tokenize(dna)
        if len(tokens) < max_len:
            tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
        else:
            tokens = tokens[:max_len]
        sequences.append(tokens)

    input_ids = torch.tensor(sequences, device=device)
    masked_ids, labels = mask_tokens_span(
        input_ids,
        vocab,
        mask_prob=0.15,
        max_span_len=1,
        exclude_special_from_labels=True,
    )
    masked_ids = masked_ids.to(device)
    labels = labels.to(device)

    # Forward + backward
    model.zero_grad()
    logits, loss = model(masked_ids, labels)
    loss.backward()

    LOGGER.info("Loss: %.4f", loss.item())

    # Check gradients for key parameters
    LOGGER.info("\nGradient norms by parameter group:")

    param_groups = {
        "embedding": [(n, p) for n, p in model.named_parameters() if "embedding" in n],
        "transform": [(n, p) for n, p in model.named_parameters() if "transform" in n],
        "encoder_layers": [(n, p) for n, p in model.named_parameters() if "layers" in n],
        "lm_head": [(n, p) for n, p in model.named_parameters() if "lm_head" in n],
    }

    has_zero_grads = False
    has_nan_grads = False

    for group_name, params in param_groups.items():
        total_norm = 0.0
        n_params = 0

        for name, param in params:
            if param.grad is not None:
                grad_norm = param.grad.norm().item()
                if math.isnan(grad_norm):
                    has_nan_grads = True
                    LOGGER.warning("  %s: NaN gradient!", name)
                    continue
                total_norm += grad_norm ** 2
                n_params += 1

        if n_params > 0:
            total_norm = math.sqrt(total_norm)
            LOGGER.info("  %s: %.6f (%d params)", group_name, total_norm, n_params)
            if total_norm == 0:
                has_zero_grads = True

    if has_nan_grads:
        LOGGER.error("  ❌ NaN gradients detected! Training is unstable.")
    elif has_zero_grads:
        LOGGER.warning("  ⚠️ Some components have zero gradients.")
    else:
        LOGGER.info("  ✅ Gradients are flowing through all components.")

    return {
        "loss": loss.item(),
        "has_zero_grads": has_zero_grads,
        "has_nan_grads": has_nan_grads,
    }


def test_realistic_training(vocab: DNABPETokenizer, device: torch.device, n_epochs: int = 5):
    """Run a realistic mini-training to see learning dynamics."""
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("TEST 6: REALISTIC TRAINING SIMULATION")
    LOGGER.info("=" * 60)

    # Use settings similar to actual training
    embedding_dim = 384
    transformer_layers = 6  # Smaller for speed
    transformer_heads = 12
    vocab_size = len(vocab.vocab)
    batch_size = 32
    max_len = 53
    n_batches_per_epoch = 100

    LOGGER.info("Creating model with %d layers, %d dim...", transformer_layers, embedding_dim)

    encoder = DNAEncoder(
        vocab_size=vocab_size,
        embedding_dim=embedding_dim,
        num_layers=transformer_layers,
        num_heads=transformer_heads,
        ff_dim=embedding_dim * 4,
        dropout=0.1,
        pad_token_id=vocab.pad_id,
        use_rope=False,
        ffn_type="gelu",
        norm_type="layernorm",
    )

    special_token_ids = [vocab.mask_id, vocab.unk_id, vocab.pad_id]
    model = DNAMLM(
        encoder,
        vocab_size=vocab_size,
        special_token_ids=special_token_ids,
        tie_weights=True,
        use_output_norm=True,
        use_transform_layer=True,
    ).to(device)

    total_params = sum(p.numel() for p in model.parameters())
    LOGGER.info("Model parameters: %.2fM", total_params / 1e6)

    # Test different learning rates
    learning_rates = [1e-3, 5e-4, 1e-4, 5e-5]

    for lr in learning_rates:
        LOGGER.info("\n--- Testing LR = %.0e ---", lr)

        # Reset model
        model.apply(lambda m: m.reset_parameters() if hasattr(m, 'reset_parameters') else None)

        optimizer = AdamW(model.parameters(), lr=lr, weight_decay=0.01)

        epoch_accs = []

        for epoch in range(n_epochs):
            model.train()
            total_loss = 0.0
            total_correct = 0
            total_masked = 0

            for batch_idx in range(n_batches_per_epoch):
                # Generate fresh batch
                sequences = []
                for _ in range(batch_size):
                    dna = generate_random_dna(234)
                    tokens = vocab.tokenize(dna)
                    if len(tokens) < max_len:
                        tokens = tokens + [vocab.pad_id] * (max_len - len(tokens))
                    else:
                        tokens = tokens[:max_len]
                    sequences.append(tokens)

                input_ids = torch.tensor(sequences, device=device)
                masked_ids, labels = mask_tokens_span(
                    input_ids,
                    vocab,
                    mask_prob=0.15,
                    max_span_len=1,
                    exclude_special_from_labels=True,
                )
                masked_ids = masked_ids.to(device)
                labels = labels.to(device)

                optimizer.zero_grad()
                logits, loss = model(masked_ids, labels)
                loss.backward()

                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)

                optimizer.step()

                total_loss += loss.item()

                # Accuracy
                with torch.no_grad():
                    preds = logits.argmax(dim=-1)
                    mask = labels != -100
                    total_correct += (preds[mask] == labels[mask]).sum().item()
                    total_masked += mask.sum().item()

            avg_loss = total_loss / n_batches_per_epoch
            acc = total_correct / max(1, total_masked)
            epoch_accs.append(acc)

            LOGGER.info("  Epoch %d: loss=%.4f, acc=%.4f (%.2f%%)",
                       epoch + 1, avg_loss, acc, acc * 100)

        # Check if learning is happening
        improvement = epoch_accs[-1] - epoch_accs[0]
        LOGGER.info("  Improvement from epoch 1 to %d: %.4f", n_epochs, improvement)

        if improvement < 0.01:
            LOGGER.warning("  ⚠️ Minimal improvement with LR=%.0e", lr)
        else:
            LOGGER.info("  ✅ Learning with LR=%.0e", lr)

    return {}


def main():
    """Run all diagnostic tests."""
    LOGGER.info("=" * 60)
    LOGGER.info("MLM TRAINING DIAGNOSTIC")
    LOGGER.info("=" * 60)

    # Setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LOGGER.info("Device: %s", device)

    # Load vocabulary
    vocab_path = SCRIPT_DIR.parent / "artifacts" / "vocabs" / "bpe_vocab.json"
    if not vocab_path.exists():
        LOGGER.error("Vocabulary not found at %s", vocab_path)
        LOGGER.info("Please run BPE vocabulary training first.")
        return

    vocab = load_vocab(vocab_path)
    LOGGER.info("Loaded vocabulary with %d tokens", len(vocab.vocab))

    # Run tests
    results = {}

    # Test 1: Random baseline
    results["random_baseline"] = test_random_baseline(vocab)

    # Test 2: Masking analysis
    results["masking"] = test_masking(vocab)

    # Test 3: Overfit test
    results["overfit"] = test_overfit(vocab, device, n_steps=200)

    # Test 4: Prediction analysis
    results["predictions"] = test_prediction_analysis(vocab, device)

    # Test 5: Gradient flow
    results["gradients"] = test_gradient_flow(vocab, device)

    # Test 6: Realistic training
    results["realistic"] = test_realistic_training(vocab, device, n_epochs=5)

    # Summary
    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("DIAGNOSTIC SUMMARY")
    LOGGER.info("=" * 60)

    random_acc = results["random_baseline"]["random_accuracy"]
    most_common_acc = results["random_baseline"]["most_common_accuracy"]
    peak_overfit_acc = results["overfit"]["peak_accuracy"]

    LOGGER.info("Random baseline accuracy: %.2f%%", random_acc * 100)
    LOGGER.info("Most-common-token accuracy: %.2f%%", most_common_acc * 100)
    LOGGER.info("Peak overfit accuracy: %.2f%%", peak_overfit_acc * 100)

    # Diagnose the issue
    LOGGER.info("\n" + "-" * 60)
    LOGGER.info("DIAGNOSIS:")

    if peak_overfit_acc < 0.30:
        LOGGER.error("❌ Model cannot overfit - fundamental architecture or training issue")
        LOGGER.info("   Check: Learning rate, gradient flow, weight initialization")
    elif abs(results["random_baseline"]["most_common_accuracy"] - 0.18) < 0.05:
        LOGGER.warning("⚠️ 18% is close to 'always predict most common token' accuracy")
        LOGGER.info("   The model may have mode collapsed to predicting the same tokens")
        LOGGER.info("   Check: Token distribution in predictions, entropy of outputs")
    else:
        LOGGER.info("   Model can learn but something is preventing training progress")
        LOGGER.info("   Check: Learning rate schedule, batch size, data quality")

    LOGGER.info("\n" + "=" * 60)


if __name__ == "__main__":
    main()
