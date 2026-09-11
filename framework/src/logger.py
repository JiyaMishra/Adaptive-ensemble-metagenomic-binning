"""
=============================================================
Logging Utility
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Provides consistent logging across the project.

Features:
✓ Console logging
✓ File logging
✓ Time stamps
✓ Progress updates
=============================================================
"""

import logging
from pathlib import Path
from datetime import datetime

from config import PROJECT_ROOT

# ============================================================
# Create logs folder if it doesn't exist
# ============================================================

LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

# ============================================================
# Log filename
# ============================================================

LOG_FILE = LOG_DIR / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

# ============================================================
# Logger
# ============================================================

logger = logging.getLogger("MetagenomicPipeline")
logger.setLevel(logging.INFO)

# Prevent duplicate handlers
if logger.hasHandlers():
    logger.handlers.clear()

# ============================================================
# Console Handler
# ============================================================

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# ============================================================
# File Handler
# ============================================================

file_handler = logging.FileHandler(LOG_FILE)
file_handler.setLevel(logging.INFO)

# ============================================================
# Format
# ============================================================

formatter = logging.Formatter(
    "[%(asctime)s] %(levelname)s : %(message)s",
    datefmt="%H:%M:%S"
)

console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

logger.addHandler(console_handler)
logger.addHandler(file_handler)

# ============================================================
# Helper Functions
# ============================================================

def info(message):
    logger.info(message)


def warning(message):
    logger.warning(message)


def error(message):
    logger.error(message)


def progress(current, total):

    percentage = (current / total) * 100

    logger.info(
        f"Processed {current:,}/{total:,} contigs "
        f"({percentage:.2f}%)"
    )