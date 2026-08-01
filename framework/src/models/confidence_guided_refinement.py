"""
Confidence-Guided Refinement

Reads:
- adaptive weights
- MetaBAT2 assignments
- MaxBin2 assignments
- VAMB assignments

Chooses the final bin assignment for every contig.

Output:
framework/results/ensemble/final_bins.csv
"""

from pathlib import Path
import os

import pandas as pd


# =====================================================
# Paths
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS = PROJECT_ROOT / "results"

WEIGHTS = (
    RESULTS
    / "ensemble"
    / "weights.csv"
)

METABAT = (
    RESULTS
    / "metabat2"
    / "metabat_assignments.csv"
)

MAXBIN = (
    RESULTS
    / "maxbin2"
    / "maxbin_assignments.csv"
)

VAMB = (
    RESULTS
    / "vamb"
    / "vamb_assignments.csv"
)

OUTPUT = (
    RESULTS
    / "ensemble"
    / "final_bins.csv"
)


# =====================================================
# Main
# =====================================================

def main():

    weights = pd.read_csv(WEIGHTS)

    metabat = pd.read_csv(METABAT)
    metabat.columns = ["contig_id", "metabat_bin"]

    maxbin = pd.read_csv(MAXBIN)
    maxbin.columns = ["contig_id", "maxbin_bin"]

    vamb = pd.read_csv(VAMB)
    vamb.columns = ["contig_id", "vamb_bin"]

    df = weights.merge(
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

    final = []

    for _, row in df.iterrows():

        scores = {
            "MetaBAT2": row["metabat_weight"],
            "MaxBin2": row["maxbin_weight"],
            "VAMB": row["vamb_weight"],
        }

        winner = max(
            scores,
            key=scores.get,
        )

        if winner == "MetaBAT2":
            chosen_bin = row["metabat_bin"]

        elif winner == "MaxBin2":
            chosen_bin = row["maxbin_bin"]

        else:
            chosen_bin = row["vamb_bin"]

        confidence = scores[winner]

        final.append(
            {
                "contig_id": row["contig_id"],
                "selected_tool": winner,
                "final_bin": chosen_bin,
                "confidence": round(confidence, 4),
            }
        )

    final_df = pd.DataFrame(final)

    os.makedirs(
        OUTPUT.parent,
        exist_ok=True,
    )

    final_df.to_csv(
        OUTPUT,
        index=False,
    )

    print()
    print("Confidence-guided refinement complete!")
    print(f"Contigs : {len(final_df)}")
    print(f"Saved to:\n{OUTPUT}")


if __name__ == "__main__":
    main()