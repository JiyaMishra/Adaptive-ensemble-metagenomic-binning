"""
normalize_features.py
---------------------
Normalizes numerical biological features while
keeping metadata unchanged.
"""

from pathlib import Path

import pandas as pd
from sklearn.preprocessing import StandardScaler


# =====================================================
# Main
# =====================================================

def normalize_features():

    project_root = Path(__file__).resolve().parents[1]

    input_csv = (
    project_root
    / "data"
    / "processed"
    / "metagem_1500"
    / "featurematrix.csv"
)

    output_folder = (
    project_root
    / "data"
    / "processed"
    / "metagem_1500"
)

    output_folder.mkdir(parents=True, exist_ok=True)

    output_csv = output_folder / "normalized_featurematrix.csv"

    print("=" * 60)
    print("Loading feature matrix...")
    print("=" * 60)

    df = pd.read_csv(input_csv)

    # -----------------------------------------
    # Columns NOT to normalize
    # -----------------------------------------

   metadata = [
    "contig_id",
    "length",
    "coverage"
]

    feature_columns = [
        c for c in df.columns
        if c not in metadata
    ]

    scaler = StandardScaler()

    normalized = scaler.fit_transform(df[feature_columns])

    normalized_df = pd.DataFrame(
        normalized,
        columns=feature_columns
    )

    final_df = pd.concat(
        [
            df[metadata],
            normalized_df
        ],
        axis=1
    )

    final_df.to_csv(
        output_csv,
        index=False
    )

    print("\nNormalization complete.")

    print(f"Saved to:\n{output_csv}")

    print("\nSummary")
    print("-----------------------------")
    print(f"Rows    : {len(final_df)}")
    print(f"Columns : {len(final_df.columns)}")

    return final_df


# =====================================================
# Run
# =====================================================

if __name__ == "__main__":

    normalize_features()