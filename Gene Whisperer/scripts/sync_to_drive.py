#!/usr/bin/env python3
"""
Sync artifacts to Google Drive (or any destination folder).

Copies files recursively from src to dst, preserving subfolder structure.
Only overwrites if source is newer or file sizes differ.
"""

import argparse
import os
import shutil
from pathlib import Path


def should_copy(src_path: Path, dst_path: Path) -> bool:
    """Return True if src should be copied to dst (newer or different size)."""
    if not dst_path.exists():
        return True

    src_stat = src_path.stat()
    dst_stat = dst_path.stat()

    # Copy if sizes differ
    if src_stat.st_size != dst_stat.st_size:
        return True

    # Copy if source is newer
    if src_stat.st_mtime > dst_stat.st_mtime:
        return True

    return False


def sync_directory(src: Path, dst: Path) -> tuple[int, int]:
    """
    Recursively sync src to dst.

    Returns:
        Tuple of (files_copied, bytes_copied)
    """
    files_copied = 0
    bytes_copied = 0

    if not src.exists():
        print(f"Error: Source directory does not exist: {src}")
        return files_copied, bytes_copied

    # Walk through source directory
    for src_path in src.rglob("*"):
        if src_path.is_file():
            # Compute relative path and destination
            rel_path = src_path.relative_to(src)
            dst_path = dst / rel_path

            if should_copy(src_path, dst_path):
                # Ensure parent directory exists
                dst_path.parent.mkdir(parents=True, exist_ok=True)

                # Copy the file
                shutil.copy2(src_path, dst_path)

                file_size = src_path.stat().st_size
                files_copied += 1
                bytes_copied += file_size

                print(f"  Copied: {rel_path} ({file_size:,} bytes)")

    return files_copied, bytes_copied


def format_bytes(num_bytes: int) -> str:
    """Format bytes as human-readable string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if abs(num_bytes) < 1024:
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1024
    return f"{num_bytes:.1f} TB"


def main():
    parser = argparse.ArgumentParser(
        description="Sync artifacts to Google Drive (or any destination folder).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sync_to_drive.py --dst /content/drive/MyDrive/artifacts
  python sync_to_drive.py --src ./checkpoints --dst /content/drive/MyDrive/checkpoints
        """
    )
    parser.add_argument(
        "--src",
        type=str,
        default="Gene Whisperer/artifacts",
        help="Source directory to sync (default: Gene Whisperer/artifacts)"
    )
    parser.add_argument(
        "--dst",
        type=str,
        required=True,
        help="Destination directory (required)"
    )

    args = parser.parse_args()

    src = Path(args.src)
    dst = Path(args.dst)

    print(f"Syncing: {src} -> {dst}")
    print("-" * 50)

    files_copied, bytes_copied = sync_directory(src, dst)

    print("-" * 50)
    print(f"Summary: {files_copied} file(s) copied, {format_bytes(bytes_copied)} total")


if __name__ == "__main__":
    main()
