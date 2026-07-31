"""
build_feature_matrix.py
-----------------------
Reads contigs, extracts biological features,
and builds the final feature matrix.
"""

from pathlib import Path

import pandas as pd

from readfasta import FASTAReader
from extract_features import extract_all_features


# =====================================================
# Configuration
# =====================================================

MIN_CONTIG_LENGTH = 1500


# =====================================================
# Main
# =====================================================

def build_feature_matrix():

    project_root = Path(__file__).resolve().parents[1]

    fasta_path = (
        project_root
        / "data"
        / "assemblies"
        / "ERR1018195.fasta"
    )

    output_folder = (
        project_root
        / "data"
        / "featurematrix"
    )

    output_folder.mkdir(parents=True, exist_ok=True)

    output_csv = output_folder / "featurematrix.csv"

    print("=" * 60)
    print("Reading FASTA file...")
    print("=" * 60)

    reader = FASTAReader(fasta_path)

    contigs = reader.read_contigs()

    print(f"Total contigs read : {len(contigs)}")

    # ------------------------------------------
    # Filter very short contigs
    # ------------------------------------------

    contigs = [
        c for c in contigs
        if c["length"] >= MIN_CONTIG_LENGTH
    ]

    print(f"After filtering (<{MIN_CONTIG_LENGTH} bp): {len(contigs)}")

    # ------------------------------------------
    # Extract Features
    # ------------------------------------------

    print("\nExtracting biological features...")

    features = extract_all_features(contigs)

    # ------------------------------------------
    # Create DataFrame
    # ------------------------------------------

    df = pd.DataFrame(features)

    # ------------------------------------------
    # Sort columns
    # ------------------------------------------

    fixed_columns = [
        "contig_id",
        "length",
        "coverage",
        "gc_content",
        "gc_skew",
        "at_skew",
        "entropy",
        "longest_homopolymer",
        "A_percent",
        "T_percent",
        "G_percent",
        "C_percent",
        "N_percent",
    ]

    tnf_columns = sorted(
        [c for c in df.columns if c.startswith("TNF_")]
    )

    df = df[
        fixed_columns +
        tnf_columns
    ]

    # ------------------------------------------
    # Save
    # ------------------------------------------

    df.to_csv(output_csv, index=False)

    print("\nFeature matrix saved successfully!")

    print(output_csv)

    # ------------------------------------------
    # Dataset Summary
    # ------------------------------------------

    print("\n========== Dataset Summary ==========")

    print(f"Total contigs : {len(df)}")

    print(f"Total features: {len(df.columns)}")

    print(f"Average length : {df['length'].mean():.2f}")

    print(f"Average GC     : {df['gc_content'].mean():.2f}")

    print(f"Average cov    : {df['coverage'].mean():.2f}")

    print("=====================================")

    return df


# =====================================================
# Run
# =====================================================

if __name__ == "__main__":

    build_feature_matrix()