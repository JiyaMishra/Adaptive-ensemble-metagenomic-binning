from pathlib import Path


def create_vamb_abundance(depth_file, output_file):

    with open(depth_file, "r") as fin, open(output_file, "w") as fout:

        # skip MetaBAT header
        header = fin.readline()

        # VAMB format
        fout.write("contigname\tsample1\n")

        for line in fin:
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            contig = parts[0]
            coverage = parts[2]

            fout.write(
                f"{contig}\t{coverage}\n"
            )


if __name__ == "__main__":

    PROJECT_ROOT = Path(__file__).resolve().parents[3]

    depth = (
        PROJECT_ROOT /
        "data/processed/metagem_1500/depth.txt"
    )

    output = (
        PROJECT_ROOT /
        "data/processed/metagem_1500/vamb_abundance.tsv"
    )

    create_vamb_abundance(depth, output)

    print("Created:", output)