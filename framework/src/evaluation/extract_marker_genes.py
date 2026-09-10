"""
Extract marker-gene identifiers from the MaxBin2 marker archive.

Usage:
    python framework/src/evaluation/extract_marker_genes.py
    python framework/src/evaluation/extract_marker_genes.py --archive framework/results/maxbin2/baseline.marker_of_each_bin.tar.gz
"""

from __future__ import annotations

import argparse
import tarfile
from pathlib import Path


def parse_fasta_markers(text: str):
    """Return marker identifiers from FASTA headers in a .marker.fasta file."""
    markers = []
    for line in text.splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            if not header:
                continue
            marker = header.split()[0]
            if marker:
                markers.append(marker)
    return markers


def extract_marker_genes(archive_path: Path):
    """Return a dict of bin_name -> list of marker genes."""
    archive_path = Path(archive_path)
    if not archive_path.exists():
        return {}

    marker_map = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            if not member.name.endswith(".marker.fasta"):
                continue

            extracted = archive.extractfile(member)
            if extracted is None:
                continue

            text = extracted.read().decode("utf-8", errors="replace")
            markers = parse_fasta_markers(text)
            bin_name = Path(member.name).name.replace(".marker.fasta", "")
            marker_map[bin_name] = markers

    return marker_map


def main():
    root = Path(__file__).resolve().parents[2]
    default_archive = root / "results" / "maxbin2" / "baseline.marker_of_each_bin.tar.gz"

    parser = argparse.ArgumentParser(description="Print marker genes from MaxBin2 marker archive.")
    parser.add_argument("--archive", type=Path, default=default_archive, help="Path to the MaxBin2 marker archive")
    parser.add_argument("--limit", type=int, default=20, help="Number of marker IDs to print per bin")
    args = parser.parse_args()

    marker_map = extract_marker_genes(args.archive)
    if not marker_map:
        print(f"Marker archive not found: {args.archive}")
        print()
        print("The MaxBin2 marker archive has not been generated yet.")
        print("Fix the pipeline first:")
        print("  1. Install MaxBin2 so 'run_MaxBin.pl' is available")
        print("  2. Run MaxBin2 on your assembly and abundance file")
        print("  3. Generate the archive baseline.marker_of_each_bin.tar.gz")
        print()
        print("Example:")
        print('  run_MaxBin.pl -contig "data/assemblies/metagem_1500/final.contigs.fa" -abund "data/processed/metagem_1500/depth.txt" -out "framework/results/maxbin2/bin" -thread 2')
        print("Then rerun this script.")
        return

    print(f"Found marker archive: {args.archive}")
    print(f"Bins with marker sets: {len(marker_map)}")
    print()

    for bin_name in sorted(marker_map):
        markers = marker_map[bin_name]
        print(f"{bin_name}: {len(markers)} marker genes")
        for marker in markers[: args.limit]:
            print(f"  - {marker}")
        if len(markers) > args.limit:
            print(f"  ... and {len(markers) - args.limit} more")
        print()


if __name__ == "__main__":
    main()
