from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.tokenization.bpe import BPETokenizer


def test_bpe_train_encode_decode_roundtrip(tmp_path: Path) -> None:
    sequences = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
    tokenizer = BPETokenizer.train(sequences, vocab_size=20, min_freq=2)

    token_ids, mask = tokenizer.encode(
        "ACGTACGT", add_special_tokens=True, max_length=12, pad_to_max=True
    )
    assert len(token_ids) == 12
    assert len(mask) == 12
    assert token_ids[0] == tokenizer.cls_token_id
    assert token_ids[mask.index(0) - 1] == tokenizer.sep_token_id if 0 in mask else True

    decoded = tokenizer.decode(token_ids)
    assert "ACGTACGT" in decoded

    save_path = tmp_path / "bpe.json"
    tokenizer.save(save_path)
    loaded = BPETokenizer.load(save_path)
    assert loaded.vocab == tokenizer.vocab
    assert loaded.metadata == tokenizer.metadata


def test_bpe_train_respects_max_token_length() -> None:
    sequences = ["ACGTACGTACGT"] * 8
    tokenizer = BPETokenizer.train(
        sequences,
        vocab_size=32,
        min_freq=2,
        max_token_length=4,
    )

    merged_tokens = [
        token
        for token in tokenizer.vocab
        if token not in tokenizer.reserved_tokens
    ]
    assert max(len(token) for token in merged_tokens) <= 4
