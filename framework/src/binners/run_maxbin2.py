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
    """Create the two-column abundance format expected by MaxBin2."""
    with open(depth_file) as source, open(output_file, "w") as target:
        source.readline()
        for line in source:
            fields = line.strip().split()
            if len(fields) >= 3:
                target.write(f"{fields[0]}\t{fields[2]}\n")


def run_maxbin2():

    output_dir = config.MAXBIN2_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"
    abundance_file = config.DEPTH_FILE.with_name("maxbin_abundance.txt")
    create_maxbin_abundance(config.DEPTH_FILE, abundance_file)

    command = [
    "run_MaxBin.pl",
    "-contig", str(config.FASTA_FILE),
    "-abund", str(abundance_file),
    "-out", str(output_prefix),
    "-thread", str(config.MAXBIN_THREADS),
]

    print("=" * 60)
    print("Running MaxBin2")
    print("=" * 60)

    subprocess.run(command, check=True)

    print("\nFinished.")
    print(output_dir)


if __name__ == "__main__":
    run_maxbin2()