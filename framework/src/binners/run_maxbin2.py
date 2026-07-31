
"""
Run MaxBin2 binning.
"""

import subprocess
from pathlib import Path
import sys

# Add src to Python path
SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def run_maxbin2():
    """
    Run MaxBin2 using the assembly FASTA and depth file.
    """

    fasta = config.FASTA_FILE
    depth = config.DEPTH_FILE

    output_dir = config.PROJECT_ROOT / "results" / "maxbin2"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"

    command = [
        "MaxBin",
        "-fasta", str(fasta),
        "-abund", str(depth),
        "-out", str(output_prefix),
        "-thread", str(config.MAXBIN_THREADS)
    ]

    print("\nRunning MaxBin2...\n")

    try:
        subprocess.run(command, check=True)

        print("\nMaxBin2 completed successfully.")
        print(f"Results saved in: {output_dir}")

    except FileNotFoundError:
        print("MaxBin2 is not installed or not found in PATH.")

    except subprocess.CalledProcessError as e:
        print("MaxBin2 failed.")
        print(e)


if __name__ == "__main__":
    run_maxbin2()