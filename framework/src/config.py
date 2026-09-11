"""
Central configuration for the metagenomic binning framework.
"""

from pathlib import Path


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "framework" / "results"

DATASET_NAME = "metagem_1500"
DATASET_DIR = DATA_DIR / DATASET_NAME

ASSEMBLY_DIR = DATA_DIR / "assemblies"
FEATURE_DIR = DATA_DIR / "featurematrix"
PROCESSED_DIR = DATA_DIR / "processed"

FASTA_FILE = ASSEMBLY_DIR / DATASET_NAME / "final.contigs.fa"
FEATURE_MATRIX = PROCESSED_DIR / DATASET_NAME / "featurematrix.csv"
DEPTH_FILE = PROCESSED_DIR / DATASET_NAME / "depth.txt"
NORMALIZED_MATRIX = PROCESSED_DIR / DATASET_NAME / "normalizedfeatures.csv"


def set_dataset(dataset_name: str):
    """Switch the project to a different dataset folder.

    Example:
        set_dataset("my_big_dataset")
    """
    global DATASET_NAME, DATASET_DIR, FASTA_FILE, FEATURE_MATRIX, DEPTH_FILE, NORMALIZED_MATRIX

    DATASET_NAME = dataset_name
    DATASET_DIR = DATA_DIR / DATASET_NAME

    FASTA_FILE = ASSEMBLY_DIR / DATASET_NAME / "final.contigs.fa"
    FEATURE_MATRIX = PROCESSED_DIR / DATASET_NAME / "featurematrix.csv"
    DEPTH_FILE = PROCESSED_DIR / DATASET_NAME / "depth.txt"
    NORMALIZED_MATRIX = PROCESSED_DIR / DATASET_NAME / "normalizedfeatures.csv"


# ============================================================
# BINNER OUTPUTS
# ============================================================

METABAT2_OUTPUT = RESULTS_DIR / "metabat2"
MAXBIN2_OUTPUT = RESULTS_DIR / "maxbin2"
VAMB_OUTPUT = RESULTS_DIR / "vamb"


# ============================================================
# ENSEMBLE OUTPUT
# ============================================================

ENSEMBLE_OUTPUT = RESULTS_DIR / "ensemble"


# ============================================================
# FEATURE EXTRACTION
# ============================================================

MIN_CONTIG_LENGTH = 1500
KMER_SIZE = 4
EXPORT_FULL_KMER = False


# ============================================================
# PROCESSING
# ============================================================

PROGRESS_INTERVAL = 5000
ROUND_DECIMALS = 4


# ============================================================
# FEATURE FLAGS
# ============================================================

ENABLE_GC_CONTENT = True
ENABLE_GC_SKEW = True
ENABLE_AT_SKEW = True
ENABLE_BASE_COUNTS = True
ENABLE_ENTROPY = True
ENABLE_HOMOPOLYMER = True
ENABLE_N_PERCENTAGE = True
ENABLE_KMER = True


# ============================================================
# BINNING PARAMETERS
# ============================================================

THREADS = 2
METABAT_MIN_CONTIG = 1500
MAXBIN_THREADS = 2
MAXBIN_MIN_CONTIG_LENGTH = 500


# ============================================================
# RANDOM SEED
# ============================================================

RANDOM_SEED = 42
