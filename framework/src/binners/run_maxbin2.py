"""
Run MaxBin2.
"""

import subprocess
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def run_maxbin2():

    output_dir = config.MAXBIN2_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"

    command = [
    "run_MaxBin.pl",
    "-contig", str(config.FASTA_FILE),
    "-abund", str(config.DEPTH_FILE),
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