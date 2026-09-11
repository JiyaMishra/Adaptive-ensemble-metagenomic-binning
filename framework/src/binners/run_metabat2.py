"""
Run MetaBAT2.
"""

import subprocess
import sys
import os
import csv
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1]
sys.path.append(str(SRC_DIR))

import config


def validate_inputs():
    """Report contigs that have no matching coverage entry."""
    fasta_headers = {
        line[1:].strip().split()[0]
        for line in config.FASTA_FILE.open()
        if line.startswith(">")
    }
    with config.DEPTH_FILE.open() as depth_file:
        depth_file.readline()
        depth_headers = {
            line.split()[0]
            for line in depth_file
            if line.strip()
        }

    missing = fasta_headers - depth_headers
    if missing:
        print(
            f"WARNING: {len(missing)} FASTA contigs are missing from the depth file."
        )
    else:
        print(f"Validated depth coverage for {len(fasta_headers)} FASTA contigs.")


def write_fallback_assignments(output_dir):
    output_file = output_dir / "metabat_assignments.csv"
    with output_file.open("w", newline="") as assignments_file:
        writer = csv.writer(assignments_file)
        writer.writerow(["Contig", "Bin"])
    print(
        "[WARNING] MetaBAT2 segfaulted on this dataset. "
        f"Generated empty fallback assignments at {output_file}."
    )


def run_metabat2():

    output_dir = config.METABAT2_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / "bin"
    for bin_file in output_dir.glob("bin*.fa"):
        bin_file.unlink()
    validate_inputs()

    command = [
        "metabat2",
        "-i", str(config.FASTA_FILE),
        "-a", str(config.DEPTH_FILE),
        "-o", str(output_prefix),
        "-m", "1500",
        "-s", str(config.METABAT_MIN_CONTIG),
        "-x", "0",
        "--minCVSum", "0",
        "--maxEdges", "500",
        "--minS", "60",
        "--noAdd",
        "-t", str(config.THREADS),
        "-v",
    ]

    print("=" * 60)
    print("Running MetaBAT2")
    print("=" * 60)

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(config.THREADS)
    try:
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
    except subprocess.CalledProcessError as error:
        if error.returncode in (-11, 139):
            write_fallback_assignments(output_dir)
            return
        if error.stdout:
            print(error.stdout)
        if error.stderr:
            print(error.stderr, file=sys.stderr)
        raise

    print("\nFinished.")
    print(output_dir)


if __name__ == "__main__":
    run_metabat2()