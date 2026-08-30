from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/adaptive_decisions.csv"
OUTPUT = PROJECT_ROOT / "framework/results/refined_assignments.csv"


FEATURES = [
    "Coverage",
    "GC_Content",
    "GC_Skew",
    "AT_Skew",
    "Entropy",
]


def robust_scale(series):
    series = pd.to_numeric(series, errors="coerce")

    median = series.median()
    mad = (series - median).abs().median()

    if pd.isna(mad) or mad == 0:
        std = series.std()

        if pd.isna(std) or std == 0:
            return pd.Series(0.0, index=series.index)

        return (series - median) / std

    return (series - median) / (1.4826 * mad)


def build_bin_profiles(df):
    profiles = {}

    assigned = df[df["selected_bin"].notna()].copy()

    for bin_name, group in assigned.groupby("selected_bin"):

        profile = {
            "size": len(group),
        }

        for feature in FEATURES:
            if feature not in group.columns:
                continue

            values = pd.to_numeric(
                group[feature],
                errors="coerce"
            ).dropna()

            if values.empty:
                continue

            profile[feature] = {
                "median": values.median(),
                "mad": (values - values.median()).abs().median(),
            }

        profiles[bin_name] = profile

    return profiles


def feature_distance(row, profile):
    distances = []

    for feature in FEATURES:

        if feature not in profile:
            continue

        value = pd.to_numeric(
            pd.Series([row[feature]]),
            errors="coerce"
        ).iloc[0]

        if pd.isna(value):
            continue

        median = profile[feature]["median"]
        mad = profile[feature]["mad"]

        scale = max(
            1.4826 * mad,
            1e-6
        )

        distances.append(
            abs(value - median) / scale
        )

    if not distances:
        return np.inf

    return float(np.mean(distances))


def refine(df):

    profiles = build_bin_profiles(df)

    result = df.copy()

    result["refined_bin"] = result["selected_bin"]

    result["refinement_action"] = np.where(
        result["selected_bin"].notna(),
        "retained",
        "unassigned",
    )

    result["refinement_distance"] = np.nan

    assigned = result["selected_bin"].notna()

    for idx in result.index[assigned]:

        row = result.loc[idx]

        selected_bin = row["selected_bin"]
        selected_profile = profiles.get(selected_bin)

        if selected_profile is None:
            continue

        current_distance = feature_distance(
            row,
            selected_profile
        )

        result.at[idx, "refinement_distance"] = current_distance

        # Protect high-confidence decisions.
        if row["confidence_category"] == "HIGH":
            result.at[idx, "refinement_action"] = (
                "protected_high_confidence"
            )
            continue

        # Only inspect bins with meaningful support.
        candidate_bins = [
            name
            for name, profile in profiles.items()
            if profile["size"] >= 5
            and name != selected_bin
        ]

        best_bin = selected_bin
        best_distance = current_distance

        for candidate in candidate_bins:

            distance = feature_distance(
                row,
                profiles[candidate]
            )

            if distance < best_distance:
                best_distance = distance
                best_bin = candidate

        # Require a meaningful improvement.
        if (
            best_bin != selected_bin
            and best_distance < current_distance * 0.60
        ):
            result.at[idx, "refined_bin"] = best_bin
            result.at[idx, "refinement_action"] = (
                "reassigned_feature_coherence"
            )
            result.at[idx, "refinement_distance"] = best_distance

    return result


def main():

    print("Loading:")
    print(INPUT)

    df = pd.read_csv(INPUT)

    print("Input rows:", len(df))

    required = [
        "contig_id",
        "selected_bin",
        "confidence_category",
    ]

    missing = [
        col for col in required
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    print("\n===== BUILD BIN PROFILES =====")

    profiles = build_bin_profiles(df)

    print("Profiles:", len(profiles))

    result = refine(df)

    print("\n===== REFINEMENT SUMMARY =====")

    print(
        result["refinement_action"]
        .value_counts(dropna=False)
    )

    print("\n===== ASSIGNMENT COUNTS =====")

    print(
        "Original assigned:",
        result["selected_bin"].notna().sum()
    )

    print(
        "Refined assigned:",
        result["refined_bin"].notna().sum()
    )

    print(
        "Reassigned:",
        (
            result["refined_bin"]
            != result["selected_bin"]
        ).sum()
    )

    print("\n===== REFINED BIN COUNT =====")

    print(
        "Original bins:",
        result["selected_bin"].dropna().nunique()
    )

    print(
        "Refined bins:",
        result["refined_bin"].dropna().nunique()
    )

    result.to_csv(
        OUTPUT,
        index=False
    )

    print("\nSaved:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
