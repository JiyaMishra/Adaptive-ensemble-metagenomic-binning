from pathlib import Path
import subprocess
import os


def create_vamb_abundance(depth_file, output_file):

    """
    Convert MetaBAT depth file
    into VAMB abundance format
    """

    with open(depth_file) as fin, open(output_file, "w") as fout:

        header = fin.readline().strip().split()

        fout.write("contigname\tsample1\n")

        for line in fin:

            parts = line.strip().split()

            if len(parts) < 3:
                continue

            contig = parts[0]
            depth = parts[2]

            fout.write(
                f"{contig}\t{depth}\n"
            )


def run_vamb():

    PROJECT_ROOT = Path(__file__).resolve().parents[3]


    fasta = (
        PROJECT_ROOT /
        "data/assemblies/ERR13958310_1500/final.contigs.fa"
    )


    depth = (
        PROJECT_ROOT /
        "data/processed/ERR13958310/depth.txt"
    )


    abundance = (
        PROJECT_ROOT /
        "data/processed/ERR13958310_1500/vamb_abundance.tsv"
    )


    output = (
        PROJECT_ROOT /
        "framework/results/vamb"
    )


    if output.exists():
     import shutil
     shutil.rmtree(output)


    create_vamb_abundance(
        depth,
        abundance
    )


    command = [
    "vamb",
    "bin",
    "default",
    "--outdir",
    str(output),
    "--fasta",
    str(fasta),
    "--abundance_tsv",
    str(abundance),
    "-p",
    "1"
    
]


    print("Running VAMB:")
    print(" ".join(command))


    subprocess.run(
        command,
        check=True
    )


if __name__ == "__main__":
    run_vamb()