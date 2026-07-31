import re
from Bio import SeqIO


import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

fasta_file = os.path.join(
    BASE_DIR,
    "data",
    "assemblies",
    "ERR1018195.fasta"
)

output_file = os.path.join(
    BASE_DIR,
    "data",
    "processed",
    "depth.txt"
)

def extract_header_info(header):
    """
    Extract contig name, length and coverage
    Example:
    NODE_1_length_154871_cov_13.594079

    Returns:
    name, length, coverage
    """

    length_match = re.search(r"length_(\d+)", header)
    cov_match = re.search(r"cov_([\d.]+)", header)

    if length_match and cov_match:
        length = int(length_match.group(1))
        coverage = float(cov_match.group(1))

        return header, length, coverage

    else:
        return None, None, None

def create_depth_file():

    count = 0

    with open(output_file, "w") as out:

        # MetaBAT depth header
        out.write(
            "contigName\tcontigLen\ttotalAvgDepth\n"
        )

        for record in SeqIO.parse(fasta_file, "fasta"):

            header = record.id

            contig_name, length, coverage = extract_header_info(header)

            if contig_name:

                # Convert full header:
                # NODE_1_length_154871_cov_13.594079
                # into:
                # NODE_1

                

                out.write(
                    f"{contig_name}\t{length}\t{coverage}\n"
                )

                count += 1

            else:
                print("Skipping:", header)


    print("Depth file created!")
    print("Contigs processed:", count)
    print("Saved as:", output_file)


   


if __name__ == "__main__":
    create_depth_file()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)