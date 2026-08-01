from pathlib import Path
import pandas as pd


def parse_vamb(vamb_file, output_file):

    df = pd.read_csv(
        vamb_file,
        sep="\t"
    )

    print("Columns detected:")
    print(df.columns.tolist())


    result = pd.DataFrame()

    # VAMB format:
    # clustername   contigname

    result["Contig"] = df["contigname"]

    result["Bin"] = (
        "vamb_bin_" +
        df["clustername"].astype(str)
    )


    result.to_csv(
        output_file,
        index=False
    )

    print(
        f"Saved VAMB assignments: {output_file}"
    )


if __name__ == "__main__":

    PROJECT_ROOT = Path(__file__).resolve().parents[3]


    vamb_file = (
        PROJECT_ROOT /
        "framework/results/vamb/vae_clusters_unsplit.tsv"
    )


    output = (
        PROJECT_ROOT /
        "framework/results/vamb/vamb_assignments.csv"
    )


    parse_vamb(
        vamb_file,
        output
    )