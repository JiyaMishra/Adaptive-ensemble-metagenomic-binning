"""
config.py
Central configuration for the metagenomic binning framework.
"""

from pathlib import Path

# ============================================================
# PROJECT PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT.parent / "data"

RESULTS_DIR = PROJECT_ROOT / "results"

ASSEMBLY_DIR = DATA_DIR / "assemblies" / "metagem_1500"
PROCESSED_DIR = DATA_DIR / "processed" / "metagem_1500"

FASTA_FILE = ASSEMBLY_DIR / "final.contigs.fa"
DEPTH_FILE = PROCESSED_DIR / "depth.txt"
FEATURE_MATRIX = PROCESSED_DIR / "featurematrix.csv"
NORMALIZED_MATRIX = PROCESSED_DIR / "normalizedfeatures.csv"

METABAT2_OUTPUT = RESULTS_DIR / "metabat2"
MAXBIN2_OUTPUT = RESULTS_DIR / "maxbin2"
SEMIBIN2_OUTPUT = RESULTS_DIR / "semibin2"
CONCOCT_OUTPUT = RESULTS_DIR / "concoct"
VAMB_OUTPUT = RESULTS_DIR / "vamb"
GMM_OUTPUT = RESULTS_DIR / "gmm"
# ============================================================
# OUTPUT ASSEMDIRECTORIES
# ============================================================

METABAT2_OUTPUT = RESULTS_DIR / "metabat2"

MAXBIN2_OUTPUT = RESULTS_DIR / "maxbin2"

SEMIBIN2_OUTPUT = RESULTS_DIR / "semibin2"

CONCOCT_OUTPUT = RESULTS_DIR / "concoct"

VAMB_OUTPUT = RESULTS_DIR / "vamb"

GMM_OUTPUT = RESULTS_DIR / "gmm"


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


# ============================================================
# RANDOM SEED
# ============================================================

RANDOM_SEED = 42