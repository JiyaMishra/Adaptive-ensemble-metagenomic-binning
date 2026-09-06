"""
SHAP Visualizations for Adaptive Binner Weighting

Generates SHAP summary/beeswarm plot and feature importance bar plot
for the adaptive weighting logic.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shap

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXPLAINABILITY_DIR = PROJECT_ROOT / "framework/results/explainability"
SHAP_CSV = EXPLAINABILITY_DIR / "shap_values.csv"
INPUT_ENSEMBLE = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"

FEATURES = ["length", "coverage", "entropy"]
OUTPUT_WEIGHTS = ["metabat_weight", "maxbin_weight", "vamb_weight"]


def visualize_shap():

    print("=" * 60)
    print("Generating SHAP Visualizations")
    print("=" * 60)

    if not SHAP_CSV.exists():
        print("SHAP CSV not found. Running SHAP analysis first...")
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import shap_analysis
        shap_analysis.main()

    df_shap = pd.read_csv(SHAP_CSV)
    df_ensemble = pd.read_csv(INPUT_ENSEMBLE)

    X = df_ensemble[FEATURES].copy()

    # Reconstruct 3D SHAP matrix (N_samples, N_features, N_outputs)
    N = len(df_shap)
    shap_vals = np.zeros((N, len(FEATURES), len(OUTPUT_WEIGHTS)))

    for o_idx, weight_name in enumerate(OUTPUT_WEIGHTS):
        for f_idx, feature_name in enumerate(FEATURES):
            col_name = f"{weight_name}_shap_{feature_name}"
            shap_vals[:, f_idx, o_idx] = df_shap[col_name].values

    EXPLAINABILITY_DIR.mkdir(parents=True, exist_ok=True)

    # 1. SHAP Feature Importance Bar Plot
    # Mean |SHAP| value across all samples and output weights
    mean_abs_shap = np.mean(np.abs(shap_vals), axis=(0, 2))

    plt.figure(figsize=(8, 5))
    bars = plt.barh(FEATURES, mean_abs_shap, color="#2b5c8f", edgecolor="black")
    plt.xlabel("Mean |SHAP Value| (Impact on Adaptive Weights)")
    plt.ylabel("Biological Features")
    plt.title("SHAP Feature Importance for Adaptive Binner Weighting")

    for bar in bars:
        width = bar.get_width()
        plt.text(
            width + 0.0005,
            bar.get_y() + bar.get_height() / 2,
            f"{width:.4f}",
            ha="left",
            va="center",
            fontsize=10,
        )

    plt.tight_layout()
    importance_png = EXPLAINABILITY_DIR / "shap_feature_importance.png"
    plt.savefig(importance_png, dpi=300)
    plt.close()
    print(f"Saved feature importance plot: {importance_png}")

    # 2. SHAP Summary / Beeswarm Plot
    plt.figure(figsize=(10, 6))
    shap_list = [shap_vals[:, :, i] for i in range(len(OUTPUT_WEIGHTS))]
    class_names = ["MetaBAT2 Weight", "MaxBin2 Weight", "VAMB Weight"]

    shap.summary_plot(
        shap_list,
        X,
        feature_names=FEATURES,
        class_names=class_names,
        show=False,
    )
    plt.title("SHAP Summary Plot for Adaptive Weighting Logic", fontsize=12)
    plt.tight_layout()
    summary_png = EXPLAINABILITY_DIR / "shap_summary.png"
    plt.savefig(summary_png, dpi=300)
    plt.close()
    print(f"Saved summary plot: {summary_png}")

    print("\nSHAP Visualizations generated successfully!")


def main():
    visualize_shap()


if __name__ == "__main__":
    main()
