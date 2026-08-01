"""
=============================================================
Configuration File
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Stores all project-wide parameters in one place.
Modify values here instead of changing multiple scripts.
=============================================================
"""

from pathlib import Path

# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"

ASSEMBLY_DIR = DATA_DIR / "assemblies"

FEATURE_DIR = DATA_DIR / "featurematrix"

PROCESSED_DIR = DATA_DIR / "processed"

# ============================================================
# INPUT FILE
# ============================================================

FASTA_FILE = ASSEMBLY_DIR / "ERR1018195.fasta"

# ============================================================
# OUTPUT FILES
# ============================================================

FEATURE_MATRIX = FEATURE_DIR / "featurematrix.csv"

NORMALIZED_MATRIX = PROCESSED_DIR / "normalizedfeatures.csv"

# ============================================================
# FEATURE EXTRACTION PARAMETERS
# ============================================================

# Ignore tiny contigs
MIN_CONTIG_LENGTH = 1000

# k-mer size
KMER_SIZE = 4

# Export every tetranucleotide?
# False = only summary statistics
# True = 256 extra columns
EXPORT_FULL_KMER = False

# ============================================================
# PROCESSING
# ============================================================

# Print progress every N contigs
PROGRESS_INTERVAL = 5000

# Number of decimal places
ROUND_DECIMALS = 4

# ============================================================
# OPTIONAL FEATURES
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
# RANDOM SEED
# ============================================================

RANDOM_SEED = 42