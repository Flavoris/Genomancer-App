"""Generate a preflight profile for MLM pretraining inputs."""
from __future__ import annotations

import argparse
from pathlib import Path

from gene_whisperer.training.pretrain_config import load_pretrain_config
from gene_whisperer.training.pretrain_profiler import (
    build_preflight_profile,
    load_sequences_with_sources,
    write_json_report,
)
from gene_whisperer.training.pretrain_runtime import maybe_train_tokenizer


def main() -> None:
    parser = argparse.ArgumentParser(description="Profile MLM pretraining inputs")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs" / "pretrain.yaml",
        help="Path to pretraining config",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=512,
        help="Number of MLM windows to sample for the preflight report",
    )
    args = parser.parse_args()

    config = load_pretrain_config(args.config.resolve())
    sequences, sequence_sources, source_summary = load_sequences_with_sources(
        paths=config.fasta_paths,
        max_bases_per_file=config.max_bases_per_file,
        verbose=True,
    )
    tokenizer = maybe_train_tokenizer(config, sequences)
    report = build_preflight_profile(
        config=config,
        sequences=sequences,
        tokenizer=tokenizer,
        sequence_sources=sequence_sources,
        source_summary=source_summary,
        sample_count=args.samples,
    )
    output_path = args.output or (config.output_dir / "preflight_profile.json")
    write_json_report(output_path, report)
    print(f"Saved preflight profile: {output_path}", flush=True)


if __name__ == "__main__":
    main()
