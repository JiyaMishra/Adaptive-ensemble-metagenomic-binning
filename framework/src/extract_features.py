"""
extract_features.py
-------------------
Extracts biological features from contigs.
"""

from itertools import product

from featureutils import (
    validate_sequence,
    base_frequencies,
    gc_content,
    gc_skew,
    at_skew,
    shannon_entropy,
    longest_homopolymer,
    tetranucleotide_frequencies,
    canonical_kmer
)


# =====================================================
# Generate Canonical TNFs
# =====================================================

def generate_canonical_tnfs():
    """
    Generates the complete list of canonical
    tetranucleotide frequencies.
    """

    bases = ["A", "T", "G", "C"]

    canonical = set()

    for kmer in map("".join, product(bases, repeat=4)):
        canonical.add(canonical_kmer(kmer))

    return sorted(canonical)


CANONICAL_TNFS = generate_canonical_tnfs()


# =====================================================
# Feature Extraction
# =====================================================

def extract_features(contig):
    """
    Extract all biological features from one contig.

    Parameters
    ----------
    contig : dict
        Dictionary returned by FASTAReader.

    Returns
    -------
    dict | None
        Feature dictionary.
        Returns None if contig fails validation.
    """

    sequence = contig["sequence"]

    # ---------------------------------------------
    # Validate sequence
    # ---------------------------------------------

    if not validate_sequence(sequence):
        return None

    features = {}

    # ---------------------------------------------
    # Metadata
    # ---------------------------------------------

    features["contig_id"] = contig["id"]
    features["length"] = contig["length"]
    features["coverage"] = contig.get("coverage", None)

    # ---------------------------------------------
    # Base Composition
    # ---------------------------------------------

    features.update(base_frequencies(sequence))

    # ---------------------------------------------
    # GC Features
    # ---------------------------------------------

    features["gc_content"] = gc_content(sequence)
    features["gc_skew"] = gc_skew(sequence)
    features["at_skew"] = at_skew(sequence)

    # ---------------------------------------------
    # Complexity
    # ---------------------------------------------

    features["entropy"] = shannon_entropy(sequence)
    features["longest_homopolymer"] = longest_homopolymer(sequence)

    # ---------------------------------------------
    # Canonical TNFs
    # ---------------------------------------------

    observed_tnfs = tetranucleotide_frequencies(sequence)

    for tnf in CANONICAL_TNFS:
        features[f"TNF_{tnf}"] = observed_tnfs.get(f"TNF_{tnf}", 0.0)

    return features


# =====================================================
# Batch Extraction
# =====================================================

def extract_all_features(contigs):
    """
    Extract features from every contig.

    Parameters
    ----------
    contigs : list

    Returns
    -------
    list
        List of feature dictionaries.
    """

    feature_matrix = []

    skipped = 0

    for contig in contigs:

        features = extract_features(contig)

        if features is None:
            skipped += 1
            continue

        feature_matrix.append(features)

    print(f"Processed contigs : {len(feature_matrix)}")
    print(f"Skipped contigs   : {skipped}")

    return feature_matrix


# =====================================================
# Standalone Testing
# =====================================================

if __name__ == "__main__":

    from pathlib import Path
    from readfasta import FASTAReader

    project_root = Path(__file__).resolve().parents[1]

    fasta_path = (
        project_root
        / "data"
        / "assemblies"
        / "ERR1018195.fasta"
    )

    reader = FASTAReader(fasta_path)

    contigs = reader.read_contigs()

    features = extract_all_features(contigs)

    print("\nFirst Feature Vector\n")

    for key, value in list(features[0].items())[:15]:
        print(f"{key:25}: {value}")