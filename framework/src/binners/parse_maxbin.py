
"""
Parse MaxBin2 output bins.

Returns:
    {
        contig_name: bin_name
    }
"""

from pathlib import Path
from Bio import SeqIO


def parse_maxbin(results_dir):
    """
    Parse all MaxBin2 FASTA bin files.
    """

    results_dir = Path(results_dir)

    contig_bins = {}

    fasta_files = sorted(results_dir.glob("*.fasta"))

    if not fasta_files:
        fasta_files = sorted(results_dir.glob("*.fa"))

    if not fasta_files:
        print(f"No MaxBin2 bin files found in {results_dir}")
        return {}

    for fasta in fasta_files:

        bin_name = fasta.stem

        for record in SeqIO.parse(fasta, "fasta"):

            if record.id in contig_bins:
                print(f"Warning: {record.id} appears in multiple MaxBin2 bins.")

            contig_bins[record.id] = bin_name

    return contig_bins


if __name__ == "__main__":

    bins = parse_maxbin("results/maxbin2")

    print(f"\nParsed {len(bins)} contigs.\n")

    for contig, bin_name in list(bins.items())[:10]:
        print(contig, "->", bin_name)