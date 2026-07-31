"""
Convert GMM clustering output
into standard bin assignment format.
"""


from pathlib import Path
import pandas as pd


def parse_gmm():

    project_root = Path(__file__).resolve().parents[2]


    input_file = (
        project_root
        / "results"
        / "gmm"
        / "gmm_bins.csv"
    )


    output_file = (
        project_root
        / "results"
        / "gmm"
        / "gmm_assignments.csv"
    )


    df = pd.read_csv(input_file)


    parsed = pd.DataFrame(
        {
            "Contig_ID": df["Contig_ID"],
            "Bin_ID": "Bin_" + df["GMM_Bin"].astype(str),
            "Source": "GMM"
        }
    )


    parsed.to_csv(
        output_file,
        index=False
    )


    print("GMM parsing completed")
    print("Saved:", output_file)

    print(parsed.head())


if __name__ == "__main__":
    parse_gmm()