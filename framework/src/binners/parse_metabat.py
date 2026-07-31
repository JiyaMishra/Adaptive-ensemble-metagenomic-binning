from pathlib import Path
from Bio import SeqIO


def parse_metabat(output_dir):
    """
    Parse MetaBAT2 bin FASTA files.

    Returns
    -------
    dict

    {
        "bin.1": ["NODE_1","NODE_2",...],
        "bin.2": [...]
    }

    """

    output_dir = Path(output_dir)

    bins = {}

    fasta_files = sorted(output_dir.glob("*.fa"))

    for fasta in fasta_files:

        bin_name = fasta.stem

        contigs = []

        for record in SeqIO.parse(fasta, "fasta"):
            contigs.append(record.id)

        bins[bin_name] = contigs

    return bins


if __name__ == "__main__":

    import config

    bins = parse_metabat(config.METABAT_OUTPUT_DIR)

    print(f"Total bins: {len(bins)}")

    for name, contigs in bins.items():
        print(name, len(contigs))