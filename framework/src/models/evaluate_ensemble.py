from pathlib import Path
from collections import Counter
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/adaptive_ensemble_assignments.csv"
OUTPUT = PROJECT_ROOT / "framework/results/evaluation_summary.csv"


def evaluate(df, bin_column):
    assigned = df[df[bin_column].notna()].copy()

    if assigned.empty:
        return {
            "method": bin_column,
            "contigs_assigned": 0,
            "num_bins": 0,
            "assignment_rate": 0.0,
            "largest_bin_contigs": 0,
            "singleton_bins": 0,
        }

    counts = assigned[bin_column].value_counts()

    return {
        "method": bin_column,
        "contigs_assigned": len(assigned),
        "num_bins": assigned[bin_column].nunique(),
        "assignment_rate": len(assigned) / len(df),
        "largest_bin_contigs": int(counts.max()),
        "singleton_bins": int((counts == 1).sum()),
    }


def main():
    df = pd.read_csv(INPUT)

    methods = [
        "metabat_bin",
        "maxbin_bin",
        "vamb_bin",
        "adaptive_bin",
    ]

    results = [evaluate(df, method) for method in methods]

    result_df = pd.DataFrame(results)

    result_df.to_csv(OUTPUT, index=False)

    print("\n===== EVALUATION SUMMARY =====")
    print(result_df.to_string(index=False))

    print("\nSaved:", OUTPUT)


if __name__ == "__main__":
    main()
