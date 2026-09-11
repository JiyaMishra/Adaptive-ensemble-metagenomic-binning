"""
Run MaxBin2.
"""

import subprocess
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def create_maxbin_abundance(depth_file, output_file):
    """
    MaxBin2 expects a 2-column abundance file: contig_id<TAB>depth.
    The project's MetaBAT-style depth.txt has 5 columns, so we keep only the
    contig id and the total depth value (column 3) used by the downstream pipeline.
    """
    with open(depth_file) as fin, open(output_file, "w") as fout:
        first = fin.readline()
        if not first:
            return

        for line in fin:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            contig = parts[0]
            depth = parts[2]
            fout.write(f"{contig}\t{depth}\n")


def run_maxbin2():
    output_dir = config.MAXBIN2_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"
    abundance_file = output_dir / "maxbin_abundance.tsv"

    create_maxbin_abundance(config.DEPTH_FILE, abundance_file)

    command = [
        "run_MaxBin.pl",
        "-contig", str(config.FASTA_FILE),
        "-abund", str(abundance_file),
        "-out", str(output_prefix),
        "-thread", str(config.MAXBIN_THREADS),
        "-min_contig_length", str(config.MAXBIN_MIN_CONTIG_LENGTH),
    ]

    print("=" * 60)
    print("Running MaxBin2")
    print("=" * 60)
    print("Using abundance file:", abundance_file)
    print("Command:", " ".join(command))

    subprocess.run(command, check=True)

    print("\nFinished.")
    print(output_dir)


if __name__ == "__main__":
    run_maxbin2()