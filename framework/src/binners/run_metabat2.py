"""
Run MetaBAT2.
"""

import subprocess
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def run_metabat2():

    output_dir = config.METABAT2_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"

    command = [
        "metabat2",
        "-i", str(config.FASTA_FILE),
        "-a", str(config.DEPTH_FILE),
        "-o", str(output_prefix),
        "-m", "1500",
        "-s", "50000",
        "-x", "0",
        "--minCVSum", "0",
        "--maxP", "99",
        "--minS", "2",
        "-t", str(config.THREADS),
        "-v",
    ]

    print("=" * 60)
    print("Running MetaBAT2")
    print("=" * 60)

    subprocess.run(command, check=True)

    print("\nFinished.")
    print(output_dir)


if __name__ == "__main__":
    run_metabat2()