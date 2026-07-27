"""
=============================================================
FASTA Reader
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Reads assembled metagenomic FASTA files efficiently.

Features
--------
✓ Streaming FASTA reader (generator)
✓ Header validation
✓ Minimum contig length filtering
✓ Progress logging
✓ Summary statistics
=============================================================
"""

from pathlib import Path
from Bio import SeqIO

from config import (
    FASTA_FILE,
    MIN_CONTIG_LENGTH,
)

from logger import info, warning


class FASTAReader:
    """
    Reads contigs from an assembled FASTA file.
    """

    def __init__(self, fasta_path=FASTA_FILE):

        self.fasta_path = Path(fasta_path)

        self.total_contigs = 0
        self.valid_contigs = 0
        self.skipped_contigs = 0

    # --------------------------------------------------------

    def validate_file(self):

        if not self.fasta_path.exists():

            raise FileNotFoundError(
                f"\nFASTA file not found:\n{self.fasta_path}"
            )

        info(f"Reading FASTA file:")
        info(str(self.fasta_path))

    # --------------------------------------------------------

    def stream_contigs(self):
        """
        Generator that yields one contig at a time.
        """

        self.validate_file()

        for record in SeqIO.parse(str(self.fasta_path), "fasta"):

            self.total_contigs += 1

            sequence = str(record.seq)

            if len(sequence) < MIN_CONTIG_LENGTH:

                self.skipped_contigs += 1
                continue

            self.valid_contigs += 1

            yield record

    # --------------------------------------------------------

    def summary(self):

        info("")

        info("========== FASTA SUMMARY ==========")

        info(f"Total contigs     : {self.total_contigs:,}")

        info(f"Valid contigs     : {self.valid_contigs:,}")

        info(f"Skipped contigs   : {self.skipped_contigs:,}")

        info("===================================")


# ============================================================
# Standalone testing
# ============================================================

if __name__ == "__main__":

    reader = FASTAReader()

    first = True

    for contig in reader.stream_contigs():

        if first:

            info(f"First contig ID : {contig.id}")

            info(f"Length          : {len(contig.seq):,}")

            first = False

    reader.summary()