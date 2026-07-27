"""
=============================================================
Feature Utility Functions
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Computes biological and statistical features from contigs.
=============================================================
"""

import math
import re
from collections import Counter
from itertools import product

from config import (
    KMER_SIZE,
    ROUND_DECIMALS,
    EXPORT_FULL_KMER
)

# ============================================================
# Header Parser
# ============================================================

HEADER_PATTERN = re.compile(
    r"(NODE_\d+)_length_(\d+)_cov_([\d\.]+)"
)


def parse_header(header):
    """
    Example:
    NODE_1_length_154871_cov_13.594079

    Returns:
        Contig_ID
        Length
        Coverage
    """

    match = HEADER_PATTERN.match(header)

    if match:

        return (
            match.group(1),
            int(match.group(2)),
            float(match.group(3))
        )

    return header, 0, 0.0


# ============================================================
# Base Counts
# ============================================================

def base_counts(sequence):

    sequence = sequence.upper()

    counts = Counter(sequence)

    return {

        "A": counts["A"],

        "T": counts["T"],

        "G": counts["G"],

        "C": counts["C"],

        "N": counts["N"]

    }


# ============================================================
# GC %
# ============================================================

def gc_content(sequence):

    counts = base_counts(sequence)

    total = len(sequence)

    if total == 0:
        return 0

    gc = counts["G"] + counts["C"]

    return round(gc / total * 100, ROUND_DECIMALS)


# ============================================================
# GC Skew
# ============================================================

def gc_skew(sequence):

    counts = base_counts(sequence)

    g = counts["G"]

    c = counts["C"]

    if (g + c) == 0:
        return 0

    return round(
        (g - c) / (g + c),
        ROUND_DECIMALS
    )


# ============================================================
# AT Skew
# ============================================================

def at_skew(sequence):

    counts = base_counts(sequence)

    a = counts["A"]

    t = counts["T"]

    if (a + t) == 0:
        return 0

    return round(
        (a - t) / (a + t),
        ROUND_DECIMALS
    )


# ============================================================
# N Percentage
# ============================================================

def n_percentage(sequence):

    counts = base_counts(sequence)

    if len(sequence) == 0:
        return 0

    return round(
        counts["N"] / len(sequence) * 100,
        ROUND_DECIMALS
    )


# ============================================================
# Shannon Entropy
# ============================================================

def shannon_entropy(sequence):

    sequence = sequence.upper()

    counts = Counter(sequence)

    total = len(sequence)

    entropy = 0

    for nucleotide in ["A", "T", "G", "C"]:

        p = counts[nucleotide] / total

        if p > 0:

            entropy -= p * math.log2(p)

    return round(entropy, ROUND_DECIMALS)


# ============================================================
# Longest Homopolymer
# ============================================================

def longest_homopolymer(sequence):

    longest = 1

    current = 1

    for i in range(1, len(sequence)):

        if sequence[i] == sequence[i-1]:

            current += 1

            longest = max(longest, current)

        else:

            current = 1

    return longest


# ============================================================
# Generate All k-mers
# ============================================================

def generate_kmers(k=KMER_SIZE):

    alphabet = ["A", "T", "G", "C"]

    return [

        "".join(x)

        for x in product(alphabet, repeat=k)

    ]


# ============================================================
# k-mer Frequencies
# ============================================================

def kmer_frequencies(sequence, k=KMER_SIZE):

    sequence = sequence.upper()

    kmers = generate_kmers(k)

    freq = dict.fromkeys(kmers, 0)

    total = len(sequence) - k + 1

    if total <= 0:

        return freq

    for i in range(total):

        kmer = sequence[i:i+k]

        if kmer in freq:

            freq[kmer] += 1

    for key in freq:

        freq[key] /= total

    return freq


# ============================================================
# Summary Statistics of k-mers
# ============================================================

def kmer_statistics(sequence):

    freq = kmer_frequencies(sequence)

    values = list(freq.values())

    mean = sum(values) / len(values)

    variance = sum(
        (x - mean) ** 2 for x in values
    ) / len(values)

    std = math.sqrt(variance)

    summary = {

        "Kmer_Mean":

            round(mean, ROUND_DECIMALS),

        "Kmer_STD":

            round(std, ROUND_DECIMALS)

    }

    if EXPORT_FULL_KMER:

        summary.update(freq)

    return summary