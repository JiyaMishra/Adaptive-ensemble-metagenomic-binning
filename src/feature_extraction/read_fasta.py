from pathlib import Path
from Bio import SeqIO

print("Program started")

project_root = Path(__file__).resolve().parents[1]
print("Project root:", project_root)

fasta_path = project_root / "data" / "assemblies" / "ERR1018195.fasta"
print("FASTA path:", fasta_path)

print("File exists:", fasta_path.exists())

if not fasta_path.exists():
    print("❌ FASTA file not found!")
    exit()

contigs = list(SeqIO.parse(str(fasta_path), "fasta"))

print("Finished reading file")

print("Total contigs:", len(contigs))

if len(contigs) > 0:
    print("First contig:", contigs[0].id)