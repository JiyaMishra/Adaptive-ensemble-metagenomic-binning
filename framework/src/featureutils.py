"""
featureutils.py
----------------
Utility functions for extracting biological features
from DNA contigs.
"""

import math
from collections import Counter


# --------------------------------------------------
# Basic Composition
# --------------------------------------------------
def validate_sequence(sequence):
    if len(sequence) == 0:
        return False

    if sequence.count("N") / len(sequence) > 0.20:
        return False

    return True
def base_frequencies(sequence):
    """
    Returns percentage of A, T, G, C and N.
    """

    length = len(sequence)

    counts = Counter(sequence)

    return {
        "A_percent": counts.get("A", 0) / length * 100,
        "T_percent": counts.get("T", 0) / length * 100,
        "G_percent": counts.get("G", 0) / length * 100,
        "C_percent": counts.get("C", 0) / length * 100,
        "N_percent": counts.get("N", 0) / length * 100,
    }


# --------------------------------------------------
# GC Content
# --------------------------------------------------

def gc_content(sequence):

    counts = Counter(sequence)

    gc = counts.get("G", 0) + counts.get("C", 0)

    return (gc / len(sequence)) * 100


# --------------------------------------------------
# GC Skew
# --------------------------------------------------

def gc_skew(sequence):

    counts = Counter(sequence)

    g = counts.get("G", 0)
    c = counts.get("C", 0)

    if g + c == 0:
        return 0

    return (g - c) / (g + c)


# --------------------------------------------------
# AT Skew
# --------------------------------------------------

def at_skew(sequence):

    counts = Counter(sequence)

    a = counts.get("A", 0)
    t = counts.get("T", 0)

    if a + t == 0:
        return 0

    return (a - t) / (a + t)


# --------------------------------------------------
# Shannon Entropy
# --------------------------------------------------

def shannon_entropy(sequence):

    counts = Counter(sequence)

    entropy = 0

    for count in counts.values():

        p = count / len(sequence)

        entropy -= p * math.log2(p)

    return entropy


# --------------------------------------------------
# Longest Homopolymer
# --------------------------------------------------

def longest_homopolymer(sequence):

    longest = 1
    current = 1

    for i in range(1, len(sequence)):

        if sequence[i] == sequence[i - 1]:

            current += 1

            longest = max(longest, current)

        else:

            current = 1

    return longest


# --------------------------------------------------
# Reverse Complement
# --------------------------------------------------

def reverse_complement(sequence):

    table = str.maketrans(
        "ATGC",
        "TACG"
    )

    return sequence.translate(table)[::-1]


# --------------------------------------------------
# Canonical K-mer
# --------------------------------------------------

def canonical_kmer(kmer):

    rev = reverse_complement(kmer)

    return min(kmer, rev)


# --------------------------------------------------
# Canonical Tetranucleotide Frequencies
# --------------------------------------------------

def tetranucleotide_frequencies(sequence, k=4):

    counts = Counter()

    total = 0

    for i in range(len(sequence) - k + 1):

        kmer = sequence[i:i + k]

        if "N" in kmer:
            continue

        canonical = canonical_kmer(kmer)

        counts[canonical] += 1

        total += 1

    frequencies = {}

    if total == 0:
        return frequencies

    for kmer, count in counts.items():

        frequencies[f"TNF_{kmer}"] = count / total

    return frequencies