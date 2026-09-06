from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
CONFIDENCE = PROJECT_ROOT / "framework/results/confidence_scores.csv"
OUTPUT = PROJECT_ROOT / "framework/results/adaptive_decisions.csv"


def choose_assignment(row):

    assignments = {
        "MetaBAT2": row["metabat_bin"],
        "MaxBin2": row["maxbin_bin"],
        "VAMB": row["vamb_bin"],
    }

    available = {
        tool: value
        for tool, value in assignments.items()
        if pd.notna(value)
    }

    if not available:
        return pd.Series({
            "selected_bin": np.nan,
            "selected_method": "UNASSIGNED",
            "decision_reason": "no_tool_assignment",
        })

    # Exact agreement is the strongest evidence.
    counts = {}

    for value in available.values():
        counts[value] = counts.get(value, 0) + 1

    max_support = max(counts.values())

    majority_bins = [
        b for b, count in counts.items()
        if count == max_support
    ]

    if max_support >= 2 and len(majority_bins) == 1:
        selected = majority_bins[0]

        supporters = [
            tool for tool, value in available.items()
            if value == selected
        ]

        if max_support == len(available):
            reason = "all_available_tools_agree"
            method = "CONSENSUS"
        else:
            reason = f"majority_support_{max_support}_tools"
            method = "+".join(supporters)

        return pd.Series({
            "selected_bin": selected,
            "selected_method": method,
            "decision_reason": reason,
        })

    # No bin-level consensus.
    # Use confidence as the adaptive selector.
    confidence = float(row["confidence_score"])
    agreement = float(row["tool_agreement_score"])

    # Strong agreement + high confidence:
    # prefer MetaBAT2 when available because it produced
    # substantially larger coherent baseline bins on this dataset.
    if confidence >= 0.75 and agreement >= 0.60:
        if pd.notna(row["metabat_bin"]):
            return pd.Series({
                "selected_bin": row["metabat_bin"],
                "selected_method": "MetaBAT2",
                "decision_reason": "high_confidence_adaptive_selection",
            })

    # Moderate confidence:
    # prefer the method with an available assignment according
    # to a deterministic tool-priority fallback.
    if confidence >= 0.50:
        priority = [
            ("MetaBAT2", row["metabat_bin"]),
            ("MaxBin2", row["maxbin_bin"]),
            ("VAMB", row["vamb_bin"]),
        ]

        for tool, value in priority:
            if pd.notna(value):
                return pd.Series({
                    "selected_bin": value,
                    "selected_method": tool,
                    "decision_reason": "medium_confidence_adaptive_selection",
                })

    # Low confidence:
    # preserve the assignment but explicitly mark it uncertain.
    priority = [
        ("MaxBin2", row["maxbin_bin"]),
        ("MetaBAT2", row["metabat_bin"]),
        ("VAMB", row["vamb_bin"]),
    ]

    for tool, value in priority:
        if pd.notna(value):
            return pd.Series({
                "selected_bin": value,
                "selected_method": tool,
                "decision_reason": "low_confidence_fallback",
            })

    return pd.Series({
        "selected_bin": np.nan,
        "selected_method": "UNASSIGNED",
        "decision_reason": "no_valid_assignment",
    })


def main():
    print("Loading ensemble:")
    print(INPUT)

    ensemble = pd.read_csv(INPUT)
    confidence = pd.read_csv(CONFIDENCE)

    print("Ensemble rows:", len(ensemble))
    print("Confidence rows:", len(confidence))

    required = {
        "contig_id",
        "metabat_bin",
        "maxbin_bin",
        "vamb_bin",
    }

    if not required.issubset(ensemble.columns):
        raise ValueError(
            f"Missing ensemble columns: "
            f"{required - set(ensemble.columns)}"
        )

    confidence_cols = [
        "contig_id",
        "assigned_tool_count",
        "tool_agreement_score",
        "biological_consistency_score",
        "confidence_score",
        "confidence_category",
    ]

    if not set(confidence_cols).issubset(confidence.columns):
        raise ValueError(
            f"Missing confidence columns: "
            f"{set(confidence_cols) - set(confidence.columns)}"
        )

    df = ensemble.merge(
        confidence[confidence_cols],
        on="contig_id",
        how="inner",
        validate="one_to_one",
    )

    if len(df) != len(ensemble):
        raise ValueError("Merge changed ensemble row count.")

    decisions = df.apply(
        choose_assignment,
        axis=1,
    )

    df = pd.concat([df, decisions], axis=1)

    print("")
    print("===== ADAPTIVE DECISION SUMMARY =====")

    print(
        df["selected_method"]
        .value_counts(dropna=False)
        .to_string()
    )

    print("")
    print("===== DECISION REASONS =====")

    print(
        df["decision_reason"]
        .value_counts(dropna=False)
        .to_string()
    )

    print("")
    print("===== CONFIDENCE BY DECISION =====")

    print(
        df.groupby("selected_method")["confidence_score"]
        .agg(["count", "mean", "median"])
        .round(4)
    )

    print("")
    print("===== ASSIGNMENT COVERAGE =====")

    print("Input contigs :", len(df))
    print("Assigned      :", df["selected_bin"].notna().sum())
    print("Unassigned    :", df["selected_bin"].isna().sum())

    df.to_csv(OUTPUT, index=False)

    print("")
    print("Saved:", OUTPUT)


if __name__ == "__main__":
    main()
