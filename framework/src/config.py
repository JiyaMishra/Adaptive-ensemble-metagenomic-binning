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

ASSEMBLY_DIR = DATA_DIR / "assemblies"
FEATURE_DIR = DATA_DIR / "featurematrix"
PROCESSED_DIR = DATA_DIR / "processed"

FASTA_FILE = ASSEMBLY_DIR / "metagem_1500" / "final.contigs.fa"
FEATURE_MATRIX = PROCESSED_DIR / "metagem_1500" / "featurematrix.csv"
DEPTH_FILE = PROCESSED_DIR / "metagem_1500" / "depth.txt"

NORMALIZED_MATRIX = PROCESSED_DIR / "normalizedfeatures.csv"


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


# ============================================================
# RANDOM SEED
# ============================================================

RANDOM_SEED = 42


# ============================================================
# CLOSED-LOOP XAI
# ============================================================
# These conservative bounds make explanation-derived feedback useful without
# allowing one batch to dominate the following adaptive-ensemble decision.
XAI_BATCH_SIZE = 500
XAI_FEEDBACK_LEARNING_RATE = 0.10
XAI_MINIMUM_FEATURE_WEIGHT = 0.25
XAI_MAXIMUM_FEATURE_WEIGHT = 2.0
