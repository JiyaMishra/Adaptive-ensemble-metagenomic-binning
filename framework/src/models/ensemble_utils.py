"""
Utility functions for the adaptive ensemble.

Loads feature matrix and outputs from all binners
and merges them into one dataframe.
"""

from pathlib import Path
import pandas as pd

# CloudArchitechture/
PROJECT_ROOT = Path(__file__).resolve().parents[3]


def load_feature_matrix():

    path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "metagem_1500"
        / "featurematrix.csv"
    )

    df = pd.read_csv(path)

    if "Contig_ID" in df.columns:
        df = df.rename(columns={"Contig_ID": "contig_id"})

    return df


def load_metabat():

    path = (
        PROJECT_ROOT
        / "framework"
        / "results"
        / "metabat2"
        / "metabat_assignments.csv"
    )

    df = pd.read_csv(path)
    df.columns = ["contig_id", "metabat_bin"]

    return df


def load_maxbin():

    path = (
        PROJECT_ROOT
        / "framework"
        / "results"
        / "maxbin2"
        / "maxbin_assignments.csv"
    )

    df = pd.read_csv(path)
    df.columns = ["contig_id", "maxbin_bin"]

    return df


def load_vamb():

    path = (
        PROJECT_ROOT
        / "framework"
        / "results"
        / "vamb"
        / "vamb_assignments.csv"
    )

    df = pd.read_csv(path)
    df.columns = ["contig_id", "vamb_bin"]

    return df


def build_ensemble_dataframe():

    features = load_feature_matrix()
    metabat = load_metabat()
    maxbin = load_maxbin()
    vamb = load_vamb()

    df = features.merge(
        metabat,
        on="contig_id",
        how="left",
    )

    df = df.merge(
        maxbin,
        on="contig_id",
        how="left",
    )

    df = df.merge(
        vamb,
        on="contig_id",
        how="left",
    )

    return df


if __name__ == "__main__":

    df = build_ensemble_dataframe()

    print("\nRows:", len(df))
    print("Columns:", len(df.columns))
    print("\nFirst five rows:\n")
    print(df.head())