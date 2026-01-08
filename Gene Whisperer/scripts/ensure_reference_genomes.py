#!/usr/bin/env python3
"""Ensure bundled reference genomes are present and not truncated."""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

LFS_POINTER_SIGNATURE = b"version https://git-lfs.github.com/spec/v1"
USER_AGENT = "GeneWhisperer/1.0 (colab bootstrap)"


@dataclass(frozen=True)
class GenomeSource:
    """Metadata for a reference genome asset."""

    name: str
    filename: str
    min_bytes: int
    url: str
    compressed: bool


@dataclass(frozen=True)
class GenomeStatus:
    """Validation status for a genome file on disk."""

    path: Path
    size_bytes: int
    valid: bool
    reason: str


REFERENCE_GENOMES = (
    GenomeSource(
        name="B. subtilis",
        filename="B_subtilis_genome",
        min_bytes=1_000_000,
        url=(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=nuccore&id=AL009126.3&rettype=fasta&retmode=text"
        ),
        compressed=False,
    ),
    GenomeSource(
        name="S. cerevisiae",
        filename="S_cerevisiae_genome",
        min_bytes=5_000_000,
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/"
            "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
        ),
        compressed=True,
    ),
)


def log(message: str) -> None:
    """Print a standardized log line."""
    print(f"[GENOME] {message}")


def read_head(path: Path, limit: int = 200) -> bytes:
    """Read the first few bytes of a file safely."""
    try:
        with path.open("rb") as handle:
            return handle.read(limit)
    except OSError:
        return b""


def is_lfs_pointer(path: Path) -> bool:
    """Detect Git LFS pointer files by signature."""
    head = read_head(path)
    return LFS_POINTER_SIGNATURE in head


def check_genome_file(path: Path, min_bytes: int, enforce_min_size: bool) -> GenomeStatus:
    """Validate genome file existence and size integrity."""
    if not path.exists():
        return GenomeStatus(path=path, size_bytes=0, valid=False, reason="missing")

    size_bytes = path.stat().st_size
    if is_lfs_pointer(path):
        return GenomeStatus(path=path, size_bytes=size_bytes, valid=False, reason="lfs_pointer")

    if enforce_min_size and size_bytes < min_bytes:
        return GenomeStatus(path=path, size_bytes=size_bytes, valid=False, reason="too_small")

    return GenomeStatus(path=path, size_bytes=size_bytes, valid=True, reason="ok")


def format_bytes(size_bytes: int) -> str:
    """Human-readable byte formatting for logs."""
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(size_bytes)
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            return f"{size:.1f}{unit}" if unit != "B" else f"{int(size)}B"
        size /= 1024.0
    return f"{size_bytes}B"


def fetch_url_to_path(url: str, dest_path: Path, compressed: bool) -> None:
    """Download a URL to dest_path, optionally inflating gzip content."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    with urllib.request.urlopen(request, timeout=60) as response:
        with tempfile.NamedTemporaryFile(
            dir=dest_path.parent, delete=False, prefix=dest_path.name, suffix=".tmp"
        ) as tmp_handle:
            tmp_path = Path(tmp_handle.name)
            try:
                if compressed:
                    with gzip.GzipFile(fileobj=response) as gz_handle:
                        shutil.copyfileobj(gz_handle, tmp_handle)
                else:
                    shutil.copyfileobj(response, tmp_handle)
            except Exception:
                tmp_handle.close()
                tmp_path.unlink(missing_ok=True)
                raise

    os.replace(tmp_path, dest_path)


def attempt_download(
    source: GenomeSource,
    dest_dir: Path,
    enforce_min_size: bool,
) -> bool:
    """Download a genome file and validate its integrity."""
    dest_path = dest_dir / source.filename
    try:
        fetch_url_to_path(source.url, dest_path, source.compressed)
    except Exception as exc:
        log(f"{source.filename}: download failed ({exc})")
        return False

    status = check_genome_file(dest_path, source.min_bytes, enforce_min_size)
    if status.valid:
        log(f"{source.filename}: download ok ({format_bytes(status.size_bytes)})")
        return True

    log(f"{source.filename}: download invalid ({status.reason})")
    return False


def ensure_reference_genomes(
    data_dir: Path,
    allow_download: bool,
    enforce_min_size: bool,
) -> list[GenomeStatus]:
    """Validate and optionally repair reference genome files."""
    statuses: list[GenomeStatus] = []

    for source in REFERENCE_GENOMES:
        genome_path = data_dir / source.filename
        status = check_genome_file(genome_path, source.min_bytes, enforce_min_size)
        if status.valid:
            log(f"{source.filename}: ok ({format_bytes(status.size_bytes)})")
            statuses.append(status)
            continue

        log(f"{source.filename}: {status.reason} ({format_bytes(status.size_bytes)})")
        if allow_download:
            log(f"{source.filename}: attempting download from {source.url}")
            attempt_download(source, data_dir, enforce_min_size)
            status = check_genome_file(genome_path, source.min_bytes, enforce_min_size)
        statuses.append(status)

    return statuses


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse CLI arguments."""
    script_dir = Path(__file__).resolve().parent
    default_data_dir = script_dir.parent / "data"

    parser = argparse.ArgumentParser(
        description="Validate and repair bundled reference genome files."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=default_data_dir,
        help="Path to the Gene Whisperer data directory.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Only report issues without downloading.",
    )
    parser.add_argument(
        "--skip-min-size",
        action="store_true",
        help="Skip minimum file size validation.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if any file remains invalid.",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Entry point for CLI usage."""
    args = parse_args(argv)
    statuses = ensure_reference_genomes(
        data_dir=args.data_dir,
        allow_download=not args.skip_download,
        enforce_min_size=not args.skip_min_size,
    )

    invalid = [status for status in statuses if not status.valid]
    if invalid:
        log("Reference genomes still missing or invalid.")
        for status in invalid:
            log(
                f"  - {status.path.name}: {status.reason} "
                f"({format_bytes(status.size_bytes)})"
            )
        if args.strict:
            return 1
    else:
        log("All reference genomes look good.")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
