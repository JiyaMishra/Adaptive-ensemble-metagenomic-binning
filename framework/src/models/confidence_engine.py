from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
OUTPUT = PROJECT_ROOT / "framework/results/confidence_scores.csv"


TOOLS = [
    "metabat_bin",
    "maxbin_bin",
    "vamb_bin",
]

BIO_FEATURES = [
    "coverage",
    "gc_content",
    "entropy",
    "length",
]


def build_bin_memberships(df):
    """
    Build:
        tool -> bin -> set(contigs)

    This lets us compare bins across different binners
    using overlap rather than comparing their bin names.
    """

    memberships = {}

    for tool in TOOLS:

        memberships[tool] = defaultdict(set)

        assigned = df[df[tool].notna()][
            ["contig_id", tool]
        ]

        for _, row in assigned.iterrows():

            memberships[tool][
                str(row[tool])
            ].add(
                str(row["contig_id"])
            )

    return memberships


def jaccard(set_a, set_b):
    """
    Jaccard similarity between two contig sets.
    """

    if not set_a or not set_b:
        return 0.0

    intersection = len(set_a & set_b)
    union = len(set_a | set_b)

    if union == 0:
        return 0.0

    return intersection / union


def calculate_tool_agreement(row, memberships):
    """
    Compare the bin containing this contig from one tool
    against bins from the other tools.

    We use bin overlap/Jaccard similarity rather than
    comparing tool-specific bin IDs.
    """

    evidence = []

    for i, tool_a in enumerate(TOOLS):

        bin_a = row.get(tool_a)

        if pd.isna(bin_a):
            continue

        set_a = memberships[tool_a].get(
            str(bin_a),
            set()
        )

        for tool_b in TOOLS[i + 1:]:

            bin_b = row.get(tool_b)

            if pd.isna(bin_b):
                continue

            set_b = memberships[tool_b].get(
                str(bin_b),
                set()
            )

            evidence.append(
                jaccard(set_a, set_b)
            )

    if not evidence:
        return 0.0

    return float(np.mean(evidence))


def biological_consistency(row, df, memberships):
    """Calculate biological similarity to the assigned bin profile."""

    scores = []

    for tool in TOOLS:

        bin_value = row.get(tool)

        if pd.isna(bin_value):
            continue

        members = memberships[tool].get(str(bin_value), set())

        if len(members) < 3:
            continue

        subset = df[
            df["contig_id"].astype(str).isin(members)
        ]

        feature_scores = []

        for feature in BIO_FEATURES:

            if feature not in df.columns:
                continue

            values = pd.to_numeric(
                subset[feature], errors="coerce"
            ).dropna()

            value = pd.to_numeric(
                pd.Series([row.get(feature)]),
                errors="coerce"
            ).iloc[0]

            if values.empty or pd.isna(value):
                continue

            median = values.median()

            mad = np.median(np.abs(values - median))

            scale = 1.4826 * mad

            # If all values are identical, the contig is
            # perfectly consistent with the bin profile.
            if scale < 1e-9:
                score = 1.0 if abs(value - median) < 1e-9 else 0.0
            else:
                deviation = abs(value - median) / scale
                score = np.exp(-0.5 * deviation)

            feature_scores.append(float(score))

        if feature_scores:
            scores.append(float(np.mean(feature_scores)))

    if not scores:
        return 0.0

    return float(np.mean(scores))



def calculate_confidence(row, agreement, biological):

    assigned_tools = sum(
        pd.notna(row.get(tool))
        for tool in TOOLS
    )

    if assigned_tools == 0:
        return 0.0

    # One-tool assignments have limited evidence.
    if assigned_tools == 1:
        agreement_evidence = 0.35
    else:
        agreement_evidence = agreement

    assignment_coverage = assigned_tools / len(TOOLS)

    confidence = (
        0.50 * agreement_evidence
        + 0.35 * biological
        + 0.15 * assignment_coverage
    )

    return float(np.clip(confidence, 0.0, 1.0))



def confidence_category(score):

    if score >= 0.75:
        return "HIGH"

    if score >= 0.50:
        return "MEDIUM"

    return "LOW"


def main():

    print("Loading:", INPUT)

    df = pd.read_csv(INPUT)

    required = {
        "contig_id",
        "metabat_bin",
        "maxbin_bin",
        "vamb_bin",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    print("Input rows:", len(df))

    memberships = build_bin_memberships(df)

    results = []

    for _, row in df.iterrows():

        agreement = calculate_tool_agreement(
            row,
            memberships
        )

        biological = biological_consistency(
            row,
            df,
            memberships
        )

        confidence = calculate_confidence(
            row,
            agreement,
            biological
        )

        assigned_tools = [
            tool
            for tool in TOOLS
            if pd.notna(row.get(tool))
        ]

        results.append(
            {
                "contig_id":
                    row["contig_id"],

                "metabat_bin":
                    row["metabat_bin"],

                "maxbin_bin":
                    row["maxbin_bin"],

                "vamb_bin":
                    row["vamb_bin"],

                "assigned_tool_count":
                    len(assigned_tools),

                "tool_agreement_score":
                    round(agreement, 4),

                "biological_consistency_score":
                    round(biological, 4),

                "confidence_score":
                    round(confidence, 4),

                "confidence_category":
                    confidence_category(confidence),
            }
        )

    result = pd.DataFrame(results)

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
        OUTPUT,
        index=False
    )

    print("")
    print("===== CONFIDENCE ENGINE COMPLETE =====")
    print("Rows:", len(result))
    print(
        "High confidence:",
        (result["confidence_category"] == "HIGH").sum()
    )
    print(
        "Medium confidence:",
        (result["confidence_category"] == "MEDIUM").sum()
    )
    print(
        "Low confidence:",
        (result["confidence_category"] == "LOW").sum()
    )

    print("")
    print("Confidence statistics:")
    print(
        result["confidence_score"].describe()
    )

    print("")
    print("===== SAMPLE RESULTS =====")
    print(
        result[
            [
                "contig_id",
                "assigned_tool_count",
                "tool_agreement_score",
                "biological_consistency_score",
                "confidence_score",
                "confidence_category",
            ]
        ].head(10).to_string(index=False)
    )

    print("")
    print("Saved:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
