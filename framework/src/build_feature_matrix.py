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
    / "metagem_1500"
    / "final.contigs.fa"
)

    output_folder = (
    project_root
    / "data"
    / "processed"
    / "metagem_1500"
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

     # ------------------------------------------
    # Extract Features
    # ------------------------------------------

    print("\nExtracting biological features...")

    features = extract_all_features(contigs)


    # ------------------------------------------
    # Merge Coverage
    # ------------------------------------------

    features_df = pd.DataFrame(features)

    from load_depth import load_depth

    depth_path = (
        project_root
        / "data"
        / "processed"
        / "metagem_1500"
        / "depth.txt"
    )

    depth_df = load_depth(depth_path)


    features_df = features_df.merge(
        depth_df,
        on="contig_id",
        how="left"
    )


    features_df["coverage"] = features_df["coverage_y"]


    features_df = features_df.drop(
        columns=["coverage_x", "coverage_y"]
    )


    # ------------------------------------------
    # Create DataFrame
    # ------------------------------------------

    df = features_df


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