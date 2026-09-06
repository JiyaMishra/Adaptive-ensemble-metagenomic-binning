from pathlib import Path
import csv


def parse_maxbin(bin_dir, output_file):

    bin_dir = Path(bin_dir)

    assignments = []


    for fasta in sorted(bin_dir.glob("*.fasta")):

        bin_name = fasta.stem


        with open(fasta) as f:

            for line in f:

                if line.startswith(">"):

                    contig = (
                        line[1:]
                        .strip()
                        .split()[0]
                    )


                    assignments.append(
                        [
                            contig,
                            f"maxbin_{bin_name}"
                        ]
                    )


    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
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
        f"Saved MaxBin assignments: {output_file}"
    )



if __name__ == "__main__":


    PROJECT_ROOT = Path(__file__).resolve().parents[3]


    bin_dir = (
        PROJECT_ROOT /
        "framework/results/maxbin2"
    )


    output = (
        PROJECT_ROOT /
        "framework/results/maxbin2/maxbin_assignments.csv"
    )


    parse_maxbin(
        bin_dir,
        output
    )