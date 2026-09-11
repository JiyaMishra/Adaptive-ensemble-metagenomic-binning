"""
Reads a FASTA assembly file and extracts:

- Contig ID
- Sequence
- Length
- Coverage (from FASTA header)
"""

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

        For metagem_1500 headers such as:
            >k119_378

        Coverage is obtained separately from depth.txt,
        so only the contig ID is needed here.
        """

        contig_id = header.split()[0]
        return contig_id, None, None

    def read_contigs(self):
        """
        Reads all contigs from the FASTA file.

        Returns
        -------
        list[dict]
        """

        contigs = []

        header = None
        sequence_lines = []

        def add_contig(current_header, current_sequence):
            if current_header is None:
                return

            contig_id, length, coverage = self._extract_header_info(
                current_header
            )

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

        with self.fasta_path.open("r", encoding="utf-8") as fasta:

            for line in fasta:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    add_contig(header, sequence_lines)

                    header = line[1:]

                    sequence_lines = []

                else:
                    sequence_lines.append(line)

        add_contig(header, sequence_lines)

        return contigs


if __name__ == "__main__":

    project_root = Path(__file__).resolve().parents[2]

    fasta = (
        project_root
        / "data"
        / "assemblies"
        / "metagem_1500"
        / "final.contigs.fa"
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