from pathlib import Path
import shutil
import subprocess


def create_vamb_abundance(depth_file, output_file):
    """
    Convert the MetaBAT depth file into VAMB's expected 2-column abundance format.
    Output shape: contigname\tdepth
    """
    with open(depth_file) as fin, open(output_file, "w") as fout:
        first = fin.readline()
        if first:
            fout.write("contigname\tsample1\n")

        for line in fin:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            contig = parts[0]
            depth = parts[2]
            fout.write(f"{contig}\t{depth}\n")


def run_vamb():
    PROJECT_ROOT = Path(__file__).resolve().parents[3]

    fasta = PROJECT_ROOT / "data/assemblies/metagem_1500/final.contigs.fa"
    depth = PROJECT_ROOT / "data/processed/metagem_1500/depth.txt"
    abundance = PROJECT_ROOT / "data/processed/metagem_1500/vamb_abundance.tsv"
    output = PROJECT_ROOT / "framework/results/vamb"

    if output.exists():
        shutil.rmtree(output)

    create_vamb_abundance(depth, abundance)

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
    ]

    print("Running VAMB:")
    print(" ".join(command))
    subprocess.run(command, check=True)

    expected = [
        output / "vae_clusters_unsplit.tsv",
        output / "clusters.tsv",
        output / "vae_clusters.tsv",
    ]
    found = [p for p in expected if p.exists()]
    if not found:
        raise FileNotFoundError("VAMB did not generate a cluster TSV file in the output directory.")
    print("VAMB generated:")
    for p in found:
        print("  -", p)


if __name__ == "__main__":
    run_vamb()