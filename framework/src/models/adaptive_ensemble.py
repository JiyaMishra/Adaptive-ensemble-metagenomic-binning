from pathlib import Path
import pandas as pd
from collections import Counter

PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
OUTPUT = PROJECT_ROOT / "framework/results/adaptive_ensemble_assignments.csv"


def choose_consensus(row):
    bins = [
        row.get("metabat_bin"),
        row.get("maxbin_bin"),
        row.get("vamb_bin"),
    ]

    bins = [b for b in bins if pd.notna(b)]

    if not bins:
        return None

    counts = Counter(bins)

    # Majority vote
    best_bin, best_count = counts.most_common(1)[0]

    if best_count >= 2:
        return best_bin

    # No agreement: retain the available MetaBAT assignment
    # as the fallback because it generally provides the
    # largest baseline assignment set.
    if pd.notna(row.get("metabat_bin")):
        return row["metabat_bin"]

    if pd.notna(row.get("maxbin_bin")):
        return row["maxbin_bin"]

    return row["vamb_bin"]


def main():
    df = pd.read_csv(INPUT)

    required = {
        "contig_id",
        "metabat_bin",
        "maxbin_bin",
        "vamb_bin",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["adaptive_bin"] = df.apply(choose_consensus, axis=1)

    df["agreement_count"] = df[
        ["metabat_bin", "maxbin_bin", "vamb_bin"]
    ].apply(
        lambda row: len(
            set(x for x in row if pd.notna(x))
        ),
        axis=1,
    )

    # More useful agreement score:
    def agreement_score(row):
        bins = [x for x in row if pd.notna(x)]

        if len(bins) <= 1:
            return 1.0

        counts = Counter(bins)
        return max(counts.values()) / len(bins)

    df["agreement_score"] = df[
        ["metabat_bin", "maxbin_bin", "vamb_bin"]
    ].apply(agreement_score, axis=1)

    df.to_csv(OUTPUT, index=False)

    print("Saved:", OUTPUT)
    print("Rows:", len(df))
    print("Unique contigs:", df["contig_id"].nunique())
    print(
        "Assigned adaptive bins:",
        df["adaptive_bin"].notna().sum()
    )

    print("\nAgreement distribution:")
    print(df["agreement_score"].value_counts().sort_index())

    print("\nAdaptive bin count:")
    print(df["adaptive_bin"].nunique())


if __name__ == "__main__":
    main()
