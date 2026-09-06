"""
SHAP Analysis for Adaptive Binner Weighting

Explains how length, coverage, and entropy influence the adaptive
weighting logic for MetaBAT2, MaxBin2, and VAMB.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import shap

PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
OUTPUT_DIR = PROJECT_ROOT / "framework/results/explainability"
OUTPUT_CSV = OUTPUT_DIR / "shap_values.csv"

FEATURES = ["length", "coverage", "entropy"]
OUTPUT_WEIGHTS = ["metabat_weight", "maxbin_weight", "vamb_weight"]


def calculate_one_weight(length, coverage, entropy):

    metabat_score = 1.0
    maxbin_score = 1.0
    vamb_score = 1.0

    if length >= 10000:
        metabat_score += 2
    elif length >= 5000:
        metabat_score += 1

    if coverage >= 20:
        maxbin_score += 2
    elif coverage >= 10:
        maxbin_score += 1

    if entropy >= 1.95:
        vamb_score += 2
    elif entropy >= 1.85:
        vamb_score += 1

    total = metabat_score + maxbin_score + vamb_score
    if total == 0:
        return [1 / 3, 1 / 3, 1 / 3]

    return [
        metabat_score / total,
        maxbin_score / total,
        vamb_score / total,
    ]


def predict_weights(X):

    if isinstance(X, pd.DataFrame):
        X = X.to_numpy()

    results = []
    for row in X:
        results.append(calculate_one_weight(row[0], row[1], row[2]))

    return np.array(results)


def run_shap_analysis():

    print("=" * 60)
    print("Running SHAP Analysis for Adaptive Binner Weighting")
    print("=" * 60)

    if not INPUT.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT}")

    df = pd.read_csv(INPUT)
    print(f"Loaded ensemble dataframe: {len(df)} contigs")

    for f in FEATURES:
        if f not in df.columns:
            raise ValueError(f"Missing required feature: {f}")

    X = df[FEATURES].copy()

    # Use a representative background sample for KernelExplainer
    background_size = min(50, len(X))
    bg = shap.sample(X, background_size, random_state=42)

    explainer = shap.KernelExplainer(predict_weights, bg)
    print("Computing SHAP values...")
    shap_vals = explainer.shap_values(X)

    # Convert to standard 3D array (N_samples, N_features, N_outputs) if list returned
    if isinstance(shap_vals, list):
        # List of 3 arrays of shape (N, 3) -> stack along axis 2 -> (N, 3, 3)
        shap_vals = np.stack(shap_vals, axis=2)

    # Construct result DataFrame
    res_df = pd.DataFrame()
    if "contig_id" in df.columns:
        res_df["contig_id"] = df["contig_id"]

    for i, feature in enumerate(FEATURES):
        res_df[feature] = df[feature]

    for o_idx, weight_name in enumerate(OUTPUT_WEIGHTS):
        for f_idx, feature_name in enumerate(FEATURES):
            col_name = f"{weight_name}_shap_{feature_name}"
            res_df[col_name] = shap_vals[:, f_idx, o_idx]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(OUTPUT_CSV, index=False)

    print(f"SHAP values saved to: {OUTPUT_CSV}")
    print(f"Output shape: {res_df.shape}")

    return res_df, shap_vals, X


def main():
    run_shap_analysis()


if __name__ == "__main__":
    main()
