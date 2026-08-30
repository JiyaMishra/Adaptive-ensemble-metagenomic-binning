from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import tarfile
import re


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"

REFINED = RESULTS / "refined_assignments.csv"
DECISIONS = RESULTS / "adaptive_decisions.csv"
MARKER_ARCHIVE = RESULTS / "maxbin2" / "baseline.marker_of_each_bin.tar.gz"

OUTPUT = RESULTS / "quality_evaluation.csv"


def parse_fasta_markers(text):
    """
    Return marker-gene names found in a marker FASTA file.

    We only need the FASTA headers because each marker sequence
    represents the presence of a marker gene.
    """
    markers = []

    for line in text.splitlines():
        if line.startswith(">"):
            header = line[1:].strip()

            # Keep the first token as the marker identifier.
            marker = header.split()[0]

            if marker:
                markers.append(marker)

    return markers


def load_marker_reference():
    """
    Read the MaxBin2 marker archive and determine the expected
    marker set represented by the reference bins.
    """
    if not MARKER_ARCHIVE.exists():
        raise FileNotFoundError(MARKER_ARCHIVE)

    marker_sets = {}

    with tarfile.open(MARKER_ARCHIVE, "r:gz") as archive:
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

            bin_name = Path(member.name).name.replace(
                ".marker.fasta", ""
            )

            marker_sets[bin_name] = set(markers)

    if not marker_sets:
        raise RuntimeError("No marker FASTA files found in archive.")

    return marker_sets


def marker_quality(assignments, marker_sets):
    """
    Calculate marker-based quality statistics for each bin.

    This is intentionally conservative:
    the marker archive is used as the biological reference,
    while the assignment file determines which contigs belong
    to each predicted bin.

    Because marker FASTAs describe marker sequences rather than
    directly mapping every contig to a marker, this function
    reports structural marker-reference availability rather than
    pretending to perform a full CheckM-style genome quality test.
    """

    rows = []

    for bin_name, group in assignments.groupby("bin"):
        if pd.isna(bin_name):
            continue

        bin_name = str(bin_name)

        # Map a predicted MaxBin bin name to the corresponding
        # marker reference where possible.
        normalized = bin_name.replace("maxbin_baseline.", "baseline.")

        marker_set = marker_sets.get(normalized, set())

        rows.append(
            {
                "bin": bin_name,
                "contigs": len(group),
                "marker_reference_genes": len(marker_set),
                "marker_reference_available": int(bool(marker_set)),
            }
        )

    return pd.DataFrame(rows)


def build_assignment_table(path, column):
    df = pd.read_csv(path)

    if "contig_id" not in df.columns:
        raise ValueError(f"{path} is missing contig_id")

    if column not in df.columns:
        raise ValueError(f"{path} is missing {column}")

    out = df[["contig_id", column]].copy()
    out = out.rename(columns={column: "bin"})

    out = out[out["bin"].notna()].copy()

    return out


def main():
    print("===== LOAD INPUTS =====")

    if not REFINED.exists():
        raise FileNotFoundError(REFINED)

    if not DECISIONS.exists():
        raise FileNotFoundError(DECISIONS)

    print("Refined:", REFINED)
    print("Decisions:", DECISIONS)
    print("Marker archive:", MARKER_ARCHIVE)

    print("\n===== LOAD MARKER REFERENCE =====")

    marker_sets = load_marker_reference()

    print("Marker reference bins:", len(marker_sets))
    print(
        "Marker genes represented:",
        sum(len(x) for x in marker_sets.values())
    )

    print("\n===== BUILD REFINED ASSIGNMENTS =====")

    refined_df = pd.read_csv(REFINED)

    refined = refined_df[
        ["contig_id", "refined_bin"]
    ].copy()

    refined = refined.rename(columns={"refined_bin": "bin"})
    refined = refined[refined["bin"].notna()].copy()

    print("Refined assigned contigs:", len(refined))
    print("Refined bins:", refined["bin"].nunique())

    print("\n===== BUILD ORIGINAL ASSIGNMENTS =====")

    decisions = pd.read_csv(DECISIONS)

    original = decisions[
        ["contig_id", "selected_bin"]
    ].copy()

    original = original.rename(columns={"selected_bin": "bin"})
    original = original[original["bin"].notna()].copy()

    print("Original assigned contigs:", len(original))
    print("Original bins:", original["bin"].nunique())

    print("\n===== MARKER REFERENCE COVERAGE =====")

    refined_quality = marker_quality(
        refined,
        marker_sets
    )

    original_quality = marker_quality(
        original,
        marker_sets
    )

    refined_quality["stage"] = "refined"
    original_quality["stage"] = "original"

    result = pd.concat(
        [original_quality, refined_quality],
        ignore_index=True
    )

    result = result[
        [
            "stage",
            "bin",
            "contigs",
            "marker_reference_genes",
            "marker_reference_available",
        ]
    ]

    result.to_csv(OUTPUT, index=False)

    print("\n===== QUALITY EVALUATION COMPLETE =====")
    print("Rows:", len(result))
    print("Original bins evaluated:",
          (result["stage"] == "original").sum())
    print("Refined bins evaluated:",
          (result["stage"] == "refined").sum())

    print("\n===== MARKER REFERENCE AVAILABILITY =====")
    print(
        result.groupby("stage")[
            "marker_reference_available"
        ].agg(["count", "sum", "mean"])
    )

    print("\nSaved:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
