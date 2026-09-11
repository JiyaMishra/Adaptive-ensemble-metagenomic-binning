import pandas as pd


def load_depth(depth_file):

    depth = pd.read_csv(
        depth_file,
        sep="\t"
    )

    depth = depth[
        [
            "contigName",
            "totalAvgDepth"
        ]
    ]

    depth.columns = [
        "contig_id",
        "coverage"
    ]

    return depth


if __name__ == "__main__":

    df = load_depth(
        "../data/processed/ERR13958310/depth.txt"
    )

    print(df.head())