"""
readfasta.py
-------------
Reads a FASTA assembly file and extracts:

- Contig ID
- Sequence
- Length
- Coverage (from FASTA header)

Example header:
>NODE_4_length_111806_cov_32.509248
"""

import re
from pathlib import Path


class FASTAReader:

    def __init__(self, fasta_path):
        self.fasta_path = Path(fasta_path)

        if not self.fasta_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found:\n{self.fasta_path}"
            )

    def _extract_header_info(self, header):
        """
        Extract contig information from FASTA header.

        Example:
        NODE_4_length_111806_cov_32.509248
        """

        pattern = r"^(.*?)_length_(\d+)_cov_([\d\.]+)"

        match = re.match(pattern, header)

        if match:

            contig_id = match.group(1)

            length = int(match.group(2))

            coverage = float(match.group(3))

        else:
            # fallback for unknown header formats

            contig_id = header.split()[0]

            length = None

            coverage = None

        return contig_id, length, coverage

    def read_contigs(self):
        """
        Reads the FASTA file.

        Returns:
        list of dictionaries.
        """

        contigs = []

        header = None
        sequence_lines = []

        def add_contig(current_header, current_sequence):
            if current_header is None:
                return

            record_id = current_header.split()[0]
            contig_id, length, coverage = self._extract_header_info(record_id)
            sequence = "".join(current_sequence).upper()

            if length is None:
                length = len(sequence)

            contigs.append(
                {
                    "id": contig_id,
                    "length": length,
                    "coverage": coverage,
                    "sequence": sequence,
                }
            )

        with self.fasta_path.open(encoding="utf-8") as fasta_file:
            for line in fasta_file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    add_contig(header, sequence_lines)
                    header = line[1:].strip()
                    sequence_lines = []
                else:
                    sequence_lines.append(line)

        add_contig(header, sequence_lines)

        return contigs


if __name__ == "__main__":

    project_root = Path(__file__).resolve().parents[1]

    fasta = (
        project_root
        / "data"
        / "assemblies"
        / "ERR1018195.fasta"
    )

    reader = FASTAReader(fasta)

    contigs = reader.read_contigs()

    print(f"\nTotal contigs : {len(contigs)}\n")

    if contigs:

        first = contigs[0]

        print("First Contig")
        print("----------------------")
        print("ID        :", first["id"])
        print("Length    :", first["length"])
        print("Coverage  :", first["coverage"])
        print("Seq Start :", first["sequence"][:60], "...")