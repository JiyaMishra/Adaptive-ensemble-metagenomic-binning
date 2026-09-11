from pathlib import Path
import gzip
import shutil
import sys


def fastq_to_fasta(fastq_path, fasta_path):
    """Convert paired or single-end FASTQ to FASTA."""
    fastq_path = Path(fastq_path)
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    def iter_fastq_records(path):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                _plus = fh.readline()
                _qual = fh.readline()
                if not seq:
                    continue
                yield header.strip(), seq.strip()

    with open(fasta_path, "w") as out:
        for header, seq in iter_fastq_records(fastq_path):
            header = header[1:] if header.startswith("@") else header
            out.write(f">{header}\n{seq}\n")

    print(f"Converted FASTQ to FASTA: {fastq_path} -> {fasta_path}")


def prepare_dataset(dataset_name, fastq_path, assembly_dir=None, processed_dir=None):
    """Prepare a dataset folder suitable for this pipeline.

    Assumes the user has a FASTQ file and wants to assemble it externally.
    This function only creates the expected project folder structure and can
    optionally convert FASTQ to FASTA if an assembly is not yet present.
    """
    project_root = Path(__file__).resolve().parents[2]
    assembly_dir = Path(assembly_dir) if assembly_dir else project_root / "data" / "assemblies" / dataset_name
    processed_dir = Path(processed_dir) if processed_dir else project_root / "data" / "processed" / dataset_name

    assembly_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = assembly_dir / "final.contigs.fa"
    if fastq_path is not None and not fasta_path.exists():
        fastq_to_fasta(fastq_path, fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"No FASTA assembly found at: {fasta_path}")

    print(f"Dataset ready: {dataset_name}")
    print(f"Assembly: {fasta_path}")
    print(f"Processed dir: {processed_dir}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python prepare_new_dataset.py <dataset_name> <fastq_file>")
        raise SystemExit(1)

    dataset_name = sys.argv[1]
    fastq_path = sys.argv[2]
    prepare_dataset(dataset_name, fastq_path)
