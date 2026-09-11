from pathlib import Path
import pandas as pd
from collections import Counter

from adaptive_weight_calculator import calculate_weights

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

    # No agreement: use the current adaptive per-tool weights.  The weights can
    # contain conservative feedback from the preceding XAI batch.
    candidates = [
        ("metabat_bin", "metabat_weight"),
        ("maxbin_bin", "maxbin_weight"),
        ("vamb_bin", "vamb_weight"),
    ]
    available = [
        (row[bin_column], float(row.get(weight_column, 0.0)))
        for bin_column, weight_column in candidates
        if pd.notna(row.get(bin_column))
    ]
    return max(available, key=lambda item: item[1])[0]


def main(feature_weights=None, input_path=INPUT, output_path=OUTPUT):
    """Run consensus using optional XAI-updated weights for this batch."""
    df = pd.read_csv(input_path)

    required = {
        "contig_id",
        "metabat_bin",
        "maxbin_bin",
        "vamb_bin",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    weights = calculate_weights(df, feature_weights=feature_weights)
    df = df.merge(weights, on="contig_id", how="left", validate="one_to_one")
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

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print("Saved:", output_path)
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
