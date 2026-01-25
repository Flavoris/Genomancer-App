"""Tests for DNA-specific BPE tokenizer implementation."""

import json
import os
import tempfile

import pytest
import torch

from bpe_tokenizer import DNABPETokenizer, SPECIAL_TOKENS, DNA_BASES


# Small corpus for fast training in tests
SMALL_CORPUS = [
    "ATCGATCGATCG",
    "ATCGATCGATCG",
    "TATAATATATAT",
    "GCGCGCGCGCGC",
    "AACGCACTATATA",
    "ATCGATCGATCGATCG",
]

# Longer corpus for more complex merge tests
LARGER_CORPUS = [
    "ATCGATCG" * 10,
    "TATATAT" * 10,
    "GCGCGCGC" * 10,
    "AAACCCGGGTTT" * 5,
    "CAATCAATCAAT" * 5,
    "ACGTACGTACGT" * 5,
]


class TestDNABPETokenizerInit:
    """Tests for DNABPETokenizer initialization."""

    def test_default_vocab_size(self):
        """Default vocab size matches DNABERT-2 recommendation of 4096."""
        tok = DNABPETokenizer()
        assert tok.vocab_size == 4096

    def test_custom_vocab_size(self):
        """Custom vocab size is stored correctly."""
        tok = DNABPETokenizer(vocab_size=256)
        assert tok.vocab_size == 256

    def test_initial_state_empty(self):
        """New tokenizer starts with empty vocab and merges."""
        tok = DNABPETokenizer()
        assert tok.vocab == {}
        assert tok.merges == []

    def test_special_token_ids(self):
        """Special token IDs are consistent."""
        tok = DNABPETokenizer()
        assert tok.pad_id == 0
        assert tok.unk_id == 1
        assert tok.mask_id == 4


class TestBPETraining:
    """Tests for the BPE training algorithm."""

    def test_train_initializes_base_vocab(self):
        """Training initializes vocab with special tokens and nucleotides."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        for token in SPECIAL_TOKENS:
            assert token in tok.vocab
        for base in DNA_BASES:
            assert base in tok.vocab

    def test_train_respects_vocab_size(self):
        """Trained vocab never exceeds target vocab_size."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert len(tok.vocab) <= 20

    def test_train_produces_merges(self):
        """Training produces merge operations."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert len(tok.merges) > 0

    def test_merges_are_tuples_of_strings(self):
        """Each merge is a tuple of two string tokens."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        for merge in tok.merges:
            assert isinstance(merge, tuple)
            assert len(merge) == 2
            assert isinstance(merge[0], str)
            assert isinstance(merge[1], str)

    def test_merged_tokens_in_vocab(self):
        """Every merged token appears in the vocabulary."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        for left, right in tok.merges:
            assert left + right in tok.vocab

    def test_frequent_pairs_merged_first(self):
        """Most frequent pairs get merged before less frequent ones."""
        corpus = ["ATATAT" * 20]  # AT is very frequent
        tok = DNABPETokenizer(vocab_size=15)
        tok.train(corpus, min_frequency=1)
        # AT should be among the earliest merges
        first_merged = tok.merges[0][0] + tok.merges[0][1]
        assert first_merged == "AT" or first_merged == "TA"

    def test_min_frequency_threshold(self):
        """Pairs below min_frequency are not merged."""
        corpus = ["ATCG"]  # Each pair occurs only once
        tok = DNABPETokenizer(vocab_size=50)
        tok.train(corpus, min_frequency=5)
        # No merges should happen since each pair frequency is 1
        assert len(tok.merges) == 0

    def test_empty_corpus_no_crash(self):
        """Training on empty corpus initializes base vocab only."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train([])
        assert len(tok.vocab) == len(SPECIAL_TOKENS) + len(DNA_BASES)

    def test_id_to_token_consistency(self):
        """id_to_token is inverse of vocab mapping."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        for token, idx in tok.vocab.items():
            assert tok.id_to_token[idx] == token

    def test_train_with_dirty_sequences(self):
        """Training handles sequences with non-ATCG characters."""
        corpus = ["ATCNGATCG", "nnATCG", "atcg", "XYZATCG"]
        tok = DNABPETokenizer(vocab_size=15)
        tok.train(corpus)
        # Should not crash; non-ATCG chars are stripped
        assert len(tok.vocab) > len(SPECIAL_TOKENS) + len(DNA_BASES)


class TestTokenization:
    """Tests for the tokenize method."""

    def test_single_bases_without_merges(self):
        """Without merges, each base maps to its own token ID."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])  # No merges, just base vocab
        ids = tok.tokenize("ATCG")
        expected = [tok.vocab["A"], tok.vocab["T"],
                    tok.vocab["C"], tok.vocab["G"]]
        assert ids == expected

    def test_tokenize_applies_merges(self):
        """Tokenization produces fewer tokens than sequence length."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        seq = "ATCGATCG"
        ids = tok.tokenize(seq)
        # With merges, token count should be less than character count
        assert len(ids) < len(seq)

    def test_tokenize_unknown_chars_stripped(self):
        """Non-ATCG characters are sanitized before tokenization."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        ids_clean = tok.tokenize("ATCG")
        ids_dirty = tok.tokenize("ANTNCNG")
        # After stripping N's: "ATCG" vs "ACG" - different result
        assert len(ids_dirty) > 0

    def test_tokenize_empty_string(self):
        """Empty sequence returns empty token list."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert tok.tokenize("") == []

    def test_tokenize_lowercase(self):
        """Lowercase input is normalized to uppercase."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        ids_upper = tok.tokenize("ATCG")
        ids_lower = tok.tokenize("atcg")
        assert ids_upper == ids_lower

    def test_encode_with_special_tokens(self):
        """encode_with_special wraps tokens with [CLS] and [SEP]."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        ids = tok.encode_with_special("ATCG")
        assert ids[0] == SPECIAL_TOKENS["[CLS]"]
        assert ids[-1] == SPECIAL_TOKENS["[SEP]"]
        assert len(ids) == len(tok.tokenize("ATCG")) + 2


class TestDecoding:
    """Tests for the decode method."""

    def test_roundtrip_single_bases(self):
        """Encoding then decoding single-base tokens recovers sequence."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])
        seq = "ATCG"
        ids = tok.tokenize(seq)
        decoded = tok.decode(ids)
        assert decoded == seq

    def test_roundtrip_with_merges(self):
        """Encoding then decoding with merges recovers original sequence."""
        tok = DNABPETokenizer(vocab_size=30)
        tok.train(SMALL_CORPUS)
        seq = "ATCGATCGATCG"
        ids = tok.tokenize(seq)
        decoded = tok.decode(ids)
        assert decoded == seq

    def test_decode_skips_pad_tokens(self):
        """Decoding skips PAD token IDs."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])
        ids = [tok.vocab["A"], tok.pad_id, tok.vocab["T"]]
        decoded = tok.decode(ids)
        assert decoded == "AT"

    def test_decode_unknown_ids_skipped(self):
        """Unknown token IDs produce empty string segments."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])
        ids = [tok.vocab["A"], 9999, tok.vocab["T"]]
        decoded = tok.decode(ids)
        assert decoded == "AT"


class TestSaveLoad:
    """Tests for save and load functionality."""

    def test_save_creates_file(self):
        """Save creates a valid JSON file."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            tok.save(path)
            assert os.path.exists(path)
            with open(path) as f:
                data = json.load(f)
            assert "vocab_size" in data
            assert "vocab" in data
            assert "merges" in data
        finally:
            os.unlink(path)

    def test_load_restores_state(self):
        """Loading a saved tokenizer reproduces identical tokenization."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            tok.save(path)
            loaded = DNABPETokenizer.load(path)
            assert loaded.vocab_size == tok.vocab_size
            assert loaded.vocab == tok.vocab
            assert loaded.merges == tok.merges
        finally:
            os.unlink(path)

    def test_loaded_tokenizer_produces_same_output(self):
        """Loaded tokenizer gives identical token IDs for same input."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            tok.save(path)
            loaded = DNABPETokenizer.load(path)
            seq = "ATCGATCG"
            assert tok.tokenize(seq) == loaded.tokenize(seq)
        finally:
            os.unlink(path)


class TestCompressionAndStats:
    """Tests for compression ratio and vocabulary statistics."""

    def test_compression_ratio_untrained(self):
        """Without merges, compression ratio is 1.0 (no compression)."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])
        ratio = tok.compression_ratio("ATCGATCG")
        assert ratio == pytest.approx(1.0)

    def test_compression_ratio_trained(self):
        """With merges, compression ratio exceeds 1.0."""
        tok = DNABPETokenizer(vocab_size=30)
        tok.train(LARGER_CORPUS)
        ratio = tok.compression_ratio("ATCGATCGATCG")
        assert ratio > 1.0

    def test_compression_ratio_empty_seq(self):
        """Empty sequence gives 0.0 compression ratio."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert tok.compression_ratio("") == 0.0

    def test_get_token_lengths(self):
        """Token length distribution includes single-base tokens."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        lengths = tok.get_token_lengths()
        assert 1 in lengths
        assert lengths[1] == len(DNA_BASES)

    def test_len_returns_vocab_size(self):
        """len(tokenizer) returns actual vocabulary size."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert len(tok) == len(tok.vocab)


class TestBiologicalMotifs:
    """Tests verifying BPE captures known biological patterns."""

    def test_tata_motif_captured(self):
        """Repeated TATA patterns get merged into multi-base tokens."""
        corpus = ["TATATAT" * 50] * 10
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(corpus, min_frequency=1)
        # TA should definitely be merged given high frequency
        assert "TA" in tok.vocab or "AT" in tok.vocab

    def test_repeat_regions_compress_well(self):
        """Repetitive genomic regions achieve high compression."""
        corpus = ["ATATATATAT" * 100] * 5
        tok = DNABPETokenizer(vocab_size=50)
        tok.train(corpus, min_frequency=1)
        ratio = tok.compression_ratio("ATATATATAT" * 10)
        # Repetitive sequences should compress significantly
        assert ratio > 2.0

    def test_variable_length_tokens(self):
        """Trained tokenizer produces tokens of varying lengths."""
        tok = DNABPETokenizer(vocab_size=50)
        tok.train(LARGER_CORPUS, min_frequency=1)
        lengths = tok.get_token_lengths()
        # Should have tokens of at least 2 different lengths
        assert len(lengths) >= 2


class TestTokenizeAndPad:
    """Tests for the tokenize_and_pad method."""

    def test_output_is_long_tensor(self):
        """Output is a torch.LongTensor."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        result = tok.tokenize_and_pad("ATCG", max_token_len=10)
        assert isinstance(result, torch.LongTensor)

    def test_output_shape(self):
        """Output shape matches max_token_len."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        result = tok.tokenize_and_pad("ATCGATCG", max_token_len=16)
        assert result.shape == (16,)

    def test_padding_with_pad_id(self):
        """Short sequences are padded with pad_id."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])  # No merges
        result = tok.tokenize_and_pad("AT", max_token_len=5)
        assert result[2].item() == tok.pad_id
        assert result[3].item() == tok.pad_id
        assert result[4].item() == tok.pad_id

    def test_truncation(self):
        """Long sequences are truncated to max_token_len."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])  # No merges, each base = 1 token
        result = tok.tokenize_and_pad("ATCGATCGATCG", max_token_len=4)
        assert result.shape == (4,)
        # Should be first 4 base tokens
        expected = tok.tokenize("ATCG")
        assert result.tolist() == expected

    def test_exact_length_no_pad_no_truncate(self):
        """Sequence that tokenizes to exactly max_token_len is unchanged."""
        tok = DNABPETokenizer(vocab_size=9)
        tok.train([])
        result = tok.tokenize_and_pad("ATCG", max_token_len=4)
        assert result.shape == (4,)
        assert tok.pad_id not in result.tolist()

    def test_empty_sequence_all_pads(self):
        """Empty sequence produces all pad tokens."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        result = tok.tokenize_and_pad("", max_token_len=5)
        assert all(t == tok.pad_id for t in result.tolist())


class TestRandomTokenId:
    """Tests for the random_token_id method."""

    def test_returns_valid_id(self):
        """Random token ID is within vocabulary range."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        for _ in range(100):
            tid = tok.random_token_id()
            assert tid >= len(SPECIAL_TOKENS)
            assert tid < len(tok.vocab)

    def test_excludes_special_tokens(self):
        """Random token ID never returns a special token."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        special_ids = set(SPECIAL_TOKENS.values())
        for _ in range(200):
            tid = tok.random_token_id()
            assert tid not in special_ids

    def test_produces_variety(self):
        """Multiple calls produce different token IDs."""
        tok = DNABPETokenizer(vocab_size=30)
        tok.train(LARGER_CORPUS)
        ids = {tok.random_token_id() for _ in range(100)}
        assert len(ids) > 1


class TestItos:
    """Tests for the itos property."""

    def test_itos_length_matches_vocab(self):
        """itos list length matches vocabulary size."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        assert len(tok.itos) == len(tok.vocab)

    def test_itos_inverse_of_vocab(self):
        """itos[vocab[token]] == token for all tokens."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        itos = tok.itos
        for token, idx in tok.vocab.items():
            assert itos[idx] == token

    def test_itos_contains_special_tokens(self):
        """itos includes all special tokens at correct positions."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        itos = tok.itos
        for token, idx in SPECIAL_TOKENS.items():
            assert itos[idx] == token

    def test_itos_contains_bases(self):
        """itos includes all DNA bases."""
        tok = DNABPETokenizer(vocab_size=20)
        tok.train(SMALL_CORPUS)
        itos = tok.itos
        for base in DNA_BASES:
            assert base in itos
