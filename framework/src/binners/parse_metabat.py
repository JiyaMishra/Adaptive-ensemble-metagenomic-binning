from pathlib import Path
import csv


def parse_metabat(bin_dir, output_file):

    bin_dir = Path(bin_dir)

    assignments = []


    for fasta in sorted(bin_dir.glob("*.fa")):

        bin_name = fasta.stem

        with open(fasta) as f:

            for line in f:

                if line.startswith(">"):

                    contig = line[1:].strip().split()[0]

                    assignments.append(
                        [
                            contig,
                            f"metabat_{bin_name}"
                        ]
                    )


    with open(output_file, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow(
            [
                "Contig",
                "Bin"
            ]
        )

        writer.writerows(assignments)


    print(
        f"Saved MetaBAT assignments: {output_file}"
    )


if __name__ == "__main__":


    PROJECT_ROOT = Path(__file__).resolve().parents[3]


    bin_dir = (
        PROJECT_ROOT /
        "framework/results/metabat2"
    )


    output = (
        PROJECT_ROOT /
        "framework/results/metabat2/metabat_assignments.csv"
    )


    parse_metabat(
        bin_dir,
        output
    )