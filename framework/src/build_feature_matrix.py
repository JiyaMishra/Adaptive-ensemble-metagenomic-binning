"""
Reads contigs, extracts biological features,
and builds the final feature matrix.
"""

from pathlib import Path
import pandas as pd

from readfasta import FASTAReader
from extract_features import extract_all_features
from loaddepth import load_depth

# =====================================================
# Configuration
# =====================================================

MIN_CONTIG_LENGTH = 1500

# =====================================================
# Main
# =====================================================

def build_feature_matrix():

    # CloudArchitechture/
    project_root = Path(__file__).resolve().parents[2]

    fasta_path = (
        project_root
        / "data"
        / "assemblies"
        / "metagem_1500"
        / "final.contigs.fa"
    )

    depth_path = (
        project_root
        / "data"
        / "processed"
        / "metagem_1500"
        / "depth.txt"
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
    # Filter short contigs
    # ------------------------------------------

    contigs = [
        c
        for c in contigs
        if c["length"] >= MIN_CONTIG_LENGTH
    ]

    print(f"After filtering (<{MIN_CONTIG_LENGTH} bp): {len(contigs)}")

    # ------------------------------------------
    # Feature extraction
    # ------------------------------------------

    print("\nExtracting biological features...")

    features = extract_all_features(contigs)

    features_df = pd.DataFrame(features)

    # ------------------------------------------
    # Merge coverage
    # ------------------------------------------

    depth_df = load_depth(depth_path)

    features_df = features_df.merge(
        depth_df,
        on="contig_id",
        how="left",
    )

    if "coverage_y" in features_df.columns:

        features_df["coverage"] = features_df["coverage_y"]

        drop_cols = []

        if "coverage_x" in features_df.columns:
            drop_cols.append("coverage_x")

        drop_cols.append("coverage_y")

        features_df = features_df.drop(columns=drop_cols)

    df = features_df

    # ------------------------------------------
    # Column ordering
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

    fixed_columns = [
        c for c in fixed_columns
        if c in df.columns
    ]

    tnf_columns = sorted(
        [
            c
            for c in df.columns
            if c.startswith("TNF_")
        ]
    )

    other_columns = [
        c
        for c in df.columns
        if c not in fixed_columns + tnf_columns
    ]

    df = df[
        fixed_columns
        + other_columns
        + tnf_columns
    ]

    # ------------------------------------------
    # Save
    # ------------------------------------------

    df.to_csv(
        output_csv,
        index=False,
    )

    print()
    print("=" * 60)
    print("Feature matrix created successfully!")
    print("=" * 60)
    print(f"Contigs : {len(df)}")
    print(f"Features: {len(df.columns)}")
    print(f"Saved to:")
    print(output_csv)


if __name__ == "__main__":
    build_feature_matrix()