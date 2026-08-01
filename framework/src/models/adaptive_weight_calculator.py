"""
Adaptive Weight Calculator

Reads the merged dataframe produced by ensemble_utils.py
and computes adaptive weights for:

- MetaBAT2
- MaxBin2
- VAMB

Output:
framework/results/ensemble/weights.csv
"""

from pathlib import Path
import os

import pandas as pd

from ensemble_utils import build_ensemble_dataframe


# =====================================================
# Normalize three scores
# =====================================================

def normalize_weights(metabat, maxbin, vamb):

    total = metabat + maxbin + vamb

    if total == 0:
        return (
            1 / 3,
            1 / 3,
            1 / 3,
        )

    return (
        metabat / total,
        maxbin / total,
        vamb / total,
    )


# =====================================================
# Weight calculation
# =====================================================

def calculate_weights(df):

    weights = []

    for _, row in df.iterrows():

        # -----------------------------
        # Read biological features
        # -----------------------------

        length = row["length"]
        coverage = row["coverage"]
        entropy = row["entropy"]

        # -----------------------------
        # Initial scores
        # -----------------------------

        metabat_score = 1.0
        maxbin_score = 1.0
        vamb_score = 1.0

        # -----------------------------
        # Long contigs favour MetaBAT2
        # -----------------------------

        if length >= 10000:
            metabat_score += 2

        elif length >= 5000:
            metabat_score += 1

        # -----------------------------
        # High coverage favours MaxBin2
        # -----------------------------

        if coverage >= 20:
            maxbin_score += 2

        elif coverage >= 10:
            maxbin_score += 1

        # -----------------------------
        # High entropy favours VAMB
        # -----------------------------

        if entropy >= 1.95:
            vamb_score += 2

        elif entropy >= 1.85:
            vamb_score += 1

        # -----------------------------
        # Normalize
        # -----------------------------

        m, x, v = normalize_weights(
            metabat_score,
            maxbin_score,
            vamb_score,
        )

        weights.append(
            {
                "contig_id": row["contig_id"],
                "metabat_weight": round(m, 4),
                "maxbin_weight": round(x, 4),
                "vamb_weight": round(v, 4),
            }
        )

    return pd.DataFrame(weights)


# =====================================================
# Main
# =====================================================

def main():

    PROJECT_ROOT = Path(__file__).resolve().parents[2]

    output_dir = (
        PROJECT_ROOT
        
        / "results"
        / "ensemble"
    )

    os.makedirs(output_dir, exist_ok=True)

    output_csv = output_dir / "weights.csv"

    print("Loading merged dataframe...")

    df = build_ensemble_dataframe()

    print("Calculating adaptive weights...")

    weights = calculate_weights(df)

    weights.to_csv(
        output_csv,
        index=False,
    )

    print()
    print("Adaptive weights created successfully!")
    print(f"Contigs : {len(weights)}")
    print(f"Saved to:\n{output_csv}")


if __name__ == "__main__":
    main()