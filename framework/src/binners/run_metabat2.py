"""
Run MetaBAT2 binning.
"""

import subprocess
from pathlib import Path

from pathlib import Path
import sys

# Add framework/src to Python path
SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def run_metabat2():
    """
    Run MetaBAT2 using the assembly FASTA and depth file.
    """

    fasta = config.FASTA_FILE
    depth = config.DEPTH_FILE

    output_dir = config.PROJECT_ROOT / "results" / "metabat2"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"

    command = [
    "metabat2",
    "-i", str(fasta),
    "-o", str(output_prefix),
    "-a", str(depth),
    "-t", "2",
    "-m", "1500",
    "-s", "10000",
    "--minCV",
    "0",
    "--minCVSum",
    "0",
    "--minS",
    "2",
    "--maxP",
    "99",
    "--seed",
    "42",
    "-v"
]
    print("\nRunning MetaBAT2...\n")

    try:
        subprocess.run(command, check=True)

        print("\nMetaBAT2 completed successfully.")
        print(f"Results saved in: {output_dir}")

    except FileNotFoundError:
        print("MetaBAT2 is not installed or not found in PATH.")

    except subprocess.CalledProcessError as e:
        print("MetaBAT2 failed.")
        print(e)


if __name__ == "__main__":
    run_metabat2()