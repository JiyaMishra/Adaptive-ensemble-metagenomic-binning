"""
GMM based metagenomic binning

Input:
    data/processed/normalized_featurematrix.csv

Output:
    results/gmm/gmm_bins.csv
"""


from pathlib import Path
import pandas as pd

from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA


def run_gmm(
    n_bins=10
):

    # Project root
    project_root = Path(__file__).resolve().parents[2]


    input_file = (
        project_root
        / "data"
        / "processed"
        / "normalized_featurematrix.csv"
    )


    output_folder = (
        project_root
        / "results"
        / "gmm"
    )

    output_folder.mkdir(
        parents=True,
        exist_ok=True
    )


    output_file = (
        output_folder
        / "gmm_bins.csv"
    )


    print("="*60)
    print("Loading normalized features")
    print("="*60)


    df = pd.read_csv(input_file)


    # Keep contig names
    contig_ids = df["Contig_ID"]


    # Remove metadata
    X = df.drop(
        columns=[
            "Contig_ID"
        ]
    )


    print("Features used:", X.shape)


    # -------------------------------
    # Optional dimensional reduction
    # -------------------------------

    print("Applying PCA...")


    pca = PCA(
        n_components=10
    )

    X_pca = pca.fit_transform(X)


    print(
        "PCA shape:",
        X_pca.shape
    )


    # -------------------------------
    # GMM clustering
    # -------------------------------

    print("Running Gaussian Mixture Model...")


    gmm = GaussianMixture(
        n_components=n_bins,
        covariance_type="full",
        random_state=42
    )


    clusters = gmm.fit_predict(
        X_pca
    )


    # -------------------------------
    # Save result
    # -------------------------------

    result = pd.DataFrame(
        {
            "Contig_ID": contig_ids,
            "GMM_Bin": clusters
        }
    )


    result.to_csv(
        output_file,
        index=False
    )


    print("\nGMM completed")
    print(
        "Saved:",
        output_file
    )

    print(
        result.head()
    )


if __name__ == "__main__":

    run_gmm(
        n_bins=10
    )