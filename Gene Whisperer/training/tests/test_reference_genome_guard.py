"""Unit tests for reference genome validation helpers."""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path


def _load_module():
    script_path = (
        Path(__file__).resolve().parents[3]
        / "Gene Whisperer"
        / "scripts"
        / "ensure_reference_genomes.py"
    )
    module_name = "ensure_reference_genomes"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_is_lfs_pointer_detects_signature(tmp_path: Path) -> None:
    module = _load_module()
    pointer_path = tmp_path / "B_subtilis_genome"
    pointer_path.write_text(
        "version https://git-lfs.github.com/spec/v1\n"
        "oid sha256:deadbeef\n"
        "size 1234\n",
        encoding="utf-8",
    )

    assert module.is_lfs_pointer(pointer_path)


def test_check_genome_file_flags_small(tmp_path: Path) -> None:
    module = _load_module()
    genome_path = tmp_path / "S_cerevisiae_genome"
    genome_path.write_text(">chrI\nACGT\n", encoding="utf-8")

    status = module.check_genome_file(genome_path, min_bytes=1024, enforce_min_size=True)
    assert not status.valid
    assert status.reason == "too_small"


def test_fetch_url_to_path_supports_gzip(tmp_path: Path) -> None:
    module = _load_module()
    source_path = tmp_path / "source.fna.gz"
    dest_path = tmp_path / "dest.fna"
    payload = ">chrI\nACGTACGT\n"

    with gzip.open(source_path, "wt", encoding="utf-8") as handle:
        handle.write(payload)

    module.fetch_url_to_path(source_path.as_uri(), dest_path, compressed=True)

    assert dest_path.read_text(encoding="utf-8") == payload
