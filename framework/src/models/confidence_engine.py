from pathlib import Path
from collections import defaultdict
import json
import warnings
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import find_peaks


PROJECT_ROOT = Path(__file__).resolve().parents[3]

INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
OUTPUT = PROJECT_ROOT / "framework/results/confidence_scores.csv"
CALIBRATION_OUTPUT = PROJECT_ROOT / "framework/results/confidence_calibration.json"


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


def build_bin_memberships(df: pd.DataFrame) -> dict[str, defaultdict[str, set[str]]]:
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


def jaccard(set_a: set[str], set_b: set[str]) -> float:
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


def calculate_tool_agreement(
    row: pd.Series,
    memberships: Mapping[str, Mapping[str, set[str]]],
) -> float:
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


def biological_consistency(
    row: pd.Series,
    df: pd.DataFrame,
    memberships: Mapping[str, Mapping[str, set[str]]],
) -> float:
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



def calculate_confidence(row: pd.Series, agreement: float, biological: float) -> float:

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



def _kde_bimodality(scores: np.ndarray) -> tuple[bool, float | None]:
    """Identify a KDE valley between the two most prominent score modes."""

    if scores.size < 3 or np.unique(scores).size < 2:
        return False, None

    try:
        density = stats.gaussian_kde(scores)
        # Grid density scales with observed data size, not a fixed score cutoff.
        grid_size = max(scores.size, np.unique(scores).size)
        grid = np.linspace(scores.min(), scores.max(), grid_size)
        values = density(grid)
        peaks, _ = find_peaks(values)
    except (np.linalg.LinAlgError, ValueError):
        return False, None

    if peaks.size < 2:
        return False, None

    strongest = peaks[np.argsort(values[peaks])[-2:]]
    left_peak, right_peak = np.sort(strongest)
    valley = left_peak + np.argmin(values[left_peak : right_peak + 1])
    return True, float(grid[valley])


def calculate_adaptive_thresholds(scores: Sequence[float] | np.ndarray) -> dict[str, Any]:
    """Derive confidence boundaries solely from the observed score distribution.

    Quartiles provide robust tail boundaries.  Their multipliers are derived from
    the empirical skewness so a longer tail receives a wider boundary.  A KDE is
    additionally evaluated to record whether the data have a natural bimodal
    split; its valley is retained as diagnostic metadata rather than imposing an
    arbitrary score cutoff.
    """

    raw_scores = np.asarray(scores, dtype=float).reshape(-1)
    finite_scores = raw_scores[np.isfinite(raw_scores)]
    if finite_scores.size == 0:
        raise ValueError("Cannot calibrate confidence thresholds without finite scores.")

    sample_size = int(finite_scores.size)
    if sample_size < 50:
        warnings.warn(
            "Confidence calibration is based on fewer than 50 contigs; "
            "thresholds may be less stable.",
            RuntimeWarning,
            stacklevel=2,
        )

    q1, median, q3 = np.percentile(finite_scores, [25, 50, 75])
    iqr = float(q3 - q1)
    mad = float(np.median(np.abs(finite_scores - median)))
    score_min = float(np.min(finite_scores))
    score_max = float(np.max(finite_scores))
    score_range = score_max - score_min
    tolerance = np.finfo(float).eps * max(1.0, abs(float(median)), score_range)

    is_zero_variance = bool(score_range <= tolerance)
    if is_zero_variance:
        skewness = 0.0
    else:
        skewness = float(stats.skew(finite_scores, bias=False)) if sample_size > 2 else 0.0
        if not np.isfinite(skewness):
            skewness = 0.0

    is_degenerate = bool(iqr <= tolerance)
    kde_is_bimodal, kde_valley = _kde_bimodality(finite_scores)

    if is_degenerate:
        # This mandated MAD fallback remains empirical; clipping preserves score bounds.
        low_threshold = max(score_min, float(median - 1.25 * mad))
        high_threshold = min(score_max, float(median + 1.25 * mad))
        method = "median_mad_fallback"
    else:
        # Positive skew expands the low tail; negative skew expands the high tail.
        alpha = max(skewness, 0.0)
        beta = max(-skewness, 0.0)
        low_threshold = max(score_min, float(q1 - alpha * iqr))
        high_threshold = min(score_max, float(q3 + beta * iqr))
        method = "skew_scaled_iqr"

    return {
        "low_threshold": float(low_threshold),
        "high_threshold": float(high_threshold),
        "q1": float(q1),
        "q3": float(q3),
        "median": float(median),
        "iqr": iqr,
        "mad": mad,
        "skewness": skewness,
        "minimum": score_min,
        "maximum": score_max,
        "sample_size": sample_size,
        "is_bimodal": kde_is_bimodal,
        "kde_valley": kde_valley,
        "is_degenerate": is_degenerate,
        "is_zero_variance": is_zero_variance,
        "calibration_method": method,
    }


def confidence_category(score: float, calibration: Mapping[str, Any]) -> str:
    """Assign a category using the dataset-specific calibrated boundaries."""

    low_threshold = float(calibration["low_threshold"])
    high_threshold = float(calibration["high_threshold"])
    if bool(calibration["is_zero_variance"]):
        return "MEDIUM"
    if score < low_threshold:
        return "LOW"
    if score > high_threshold:
        return "HIGH"
    return "MEDIUM"


def main() -> dict[str, Any]:

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
    if df.empty:
        raise ValueError("Cannot calculate confidence calibration for an empty contig dataset.")

    memberships = build_bin_memberships(df)

    results = []
    raw_confidence_scores: list[float] = []

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
        raw_confidence_scores.append(confidence)

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

            }
        )

    result = pd.DataFrame(results)
    calibration = calculate_adaptive_thresholds(raw_confidence_scores)
    result["confidence_category"] = [
        confidence_category(score, calibration)
        for score in raw_confidence_scores
    ]

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
        OUTPUT,
        index=False
    )

    with CALIBRATION_OUTPUT.open("w", encoding="utf-8") as calibration_file:
        json.dump(calibration, calibration_file, indent=2, sort_keys=True)

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
    print("Calibration metadata:")
    print(CALIBRATION_OUTPUT)

    return calibration


if __name__ == "__main__":
    main()
