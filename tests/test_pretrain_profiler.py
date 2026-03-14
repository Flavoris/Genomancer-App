from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.tokenization.bpe import BPETokenizer
from gene_whisperer.training.pretrain_profiler import (
    build_preflight_profile,
    load_sequences_with_sources,
    write_json_report,
)


def _build_config() -> SimpleNamespace:
    return SimpleNamespace(
        window_size=32,
        max_length=32,
        mask_prob=0.15,
        samples_per_epoch=32,
        seed=7,
        sample_by_length=True,
        mask_ambiguous_tokens=False,
        min_masked_tokens=2,
        min_maskable_tokens=4,
        min_tokenized_tokens=8,
        resample_attempts=4,
        tokenizer_max_bases=64,
        tokenizer_max_sequences=4,
        tokenizer_window_size=16,
    )


def test_build_preflight_profile_reports_window_and_tokenizer_stats() -> None:
    config = _build_config()
    sequences = ["ACGT" * 32, "TGCA" * 32]
    sources = ["human", "yeast"]
    tokenizer = BPETokenizer.train(sequences=sequences, vocab_size=32)

    report = build_preflight_profile(
        config=config,
        sequences=sequences,
        tokenizer=tokenizer,
        sequence_sources=sources,
        source_summary={"human": {"sequence_count": 1}, "yeast": {"sequence_count": 1}},
        sample_count=8,
    )

    assert report["mlm_windows"]["sample_count"] == 8
    assert report["mlm_windows"]["tokenized_tokens_total"] > 0
    assert 0.0 <= report["mlm_windows"]["maskable_token_rate"] <= 1.0
    assert report["tokenizer"]["corpus_coverage"]["sampled_sequences"] == 4
    assert set(report["tokenizer"]["corpus_coverage"]["sources"]) == {"human", "yeast"}


def test_load_sequences_with_sources_groups_by_configured_path(tmp_path: Path) -> None:
    source_a = tmp_path / "human.fna"
    source_b = tmp_path / "yeast.fna"
    source_a.write_text(">a\nACGTACGT\n", encoding="utf-8")
    source_b.write_text(">b\nTTTTCCCC\n", encoding="utf-8")

    sequences, sources, summary = load_sequences_with_sources([source_a, source_b])

    assert sequences == ["ACGTACGT", "TTTTCCCC"]
    assert sources == ["human.fna", "yeast.fna"]
    assert summary["human.fna"]["bases"] == 8


def test_write_json_report_serializes_profile(tmp_path: Path) -> None:
    output_path = tmp_path / "profile.json"
    write_json_report(output_path, {"ok": True, "count": 3})

    assert json.loads(output_path.read_text()) == {"ok": True, "count": 3}
