"""
Adaptive weighted ensemble for metagenomic binning.

Combines MetaBAT2, MaxBin2 and VAMB
using dynamic feature-based weights.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from ensemble_utils import build_ensemble_dataframe


PROJECT_ROOT = Path(__file__).resolve().parents[3]

OUTPUT_DIR = (
    PROJECT_ROOT
    / "framework"
    / "results"
    / "ensemble"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =====================================================
# Dynamic Weights
# =====================================================

def compute_weights(row):
    """
    Compute adaptive weights using
    biological features.
    """

    length = row["length"]
    coverage = row["coverage"]
    gc = row["gc_content"]

    # default
    metabat = 0.34
    maxbin = 0.33
    vamb = 0.33

    # Long contigs
    if length >= 5000:
        metabat += 0.15

    # High coverage
    if coverage >= 20:
        maxbin += 0.10

    # Medium GC
    if 40 <= gc <= 60:
        vamb += 0.10

    total = metabat + maxbin + vamb

    return {
        "metabat": metabat / total,
        "maxbin": maxbin / total,
        "vamb": vamb / total,
    }


# =====================================================
# Weighted Voting
# =====================================================

def weighted_vote(row):

    weights = compute_weights(row)

    scores = {}

    mapping = [
        ("metabat_bin", weights["metabat"]),
        ("maxbin_bin", weights["maxbin"]),
        ("vamb_bin", weights["vamb"]),
    ]

    for column, weight in mapping:

        value = row[column]

        if pd.isna(value):
            continue

        scores[value] = scores.get(value, 0.0) + weight

    if len(scores) == 0:
        return pd.Series(
            {
                "ensemble_bin": np.nan,
                "confidence": 0.0,
            }
        )

    winner = max(scores, key=scores.get)

    confidence = scores[winner] / sum(scores.values())

    return pd.Series(
        {
            "ensemble_bin": winner,
            "confidence": round(confidence, 3),
        }
    )


# =====================================================
# Main
# =====================================================

def run_ensemble():

    df = build_ensemble_dataframe()

    predictions = df.apply(
        weighted_vote,
        axis=1,
    )

    final = pd.concat(
        [
            df,
            predictions,
        ],
        axis=1,
    )

    output = (
        OUTPUT_DIR
        / "ensemble_predictions.csv"
    )

    final.to_csv(
        output,
        index=False,
    )

    print()
    print("=" * 60)
    print("Adaptive Ensemble Complete")
    print("=" * 60)
    print("Contigs :", len(final))
    print("Output  :", output)
    print()
    print(final[
        [
            "contig_id",
            "metabat_bin",
            "maxbin_bin",
            "vamb_bin",
            "ensemble_bin",
            "confidence",
        ]
    ].head())


if __name__ == "__main__":
    run_ensemble()