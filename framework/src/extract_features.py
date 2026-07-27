"""
=============================================================
Feature Extraction Engine
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Reads streamed contigs and extracts biological features.
=============================================================
"""

from logger import info, warning, progress

from config import (
    PROGRESS_INTERVAL,
    MIN_CONTIG_LENGTH
)

from featureutils import (
    parse_header,
    base_counts,
    gc_content,
    gc_skew,
    at_skew,
    n_percentage,
    shannon_entropy,
    longest_homopolymer,
    kmer_statistics
)


class FeatureExtractor:

    def __init__(self):

        self.processed = 0
        self.skipped = 0

    # ---------------------------------------------------------

    def extract(self, contig):

        """
        Extract all features from one contig.
        """

        sequence = str(contig.seq)

        # Skip empty sequences

        if len(sequence) == 0:

            self.skipped += 1
            return None

        # Skip very short contigs

        if len(sequence) < MIN_CONTIG_LENGTH:

            self.skipped += 1
            return None

        # Parse header

        contig_id, header_length, coverage = parse_header(contig.id)

        actual_length = len(sequence)

        # Warn if header and sequence disagree

        if header_length != 0:

            if header_length != actual_length:

                warning(

                    f"{contig_id}: "

                    f"Header length ({header_length}) "

                    f"!= Actual length ({actual_length})"

                )

        # Base counts

        counts = base_counts(sequence)

        # Main feature dictionary

        features = {

            "Contig_ID":

                contig_id,

            "Sequence_Length":

                actual_length,

            "Header_Length":

                header_length,

            "Coverage":

                coverage,

            "GC_Content":

                gc_content(sequence),

            "GC_Skew":

                gc_skew(sequence),

            "AT_Skew":

                at_skew(sequence),

            "Entropy":

                shannon_entropy(sequence),

            "Longest_Homopolymer":

                longest_homopolymer(sequence),

            "N_Percentage":

                n_percentage(sequence),

            "A_Count":

                counts["A"],

            "T_Count":

                counts["T"],

            "G_Count":

                counts["G"],

            "C_Count":

                counts["C"],

            "N_Count":

                counts["N"]

        }

        # Add k-mer statistics

        features.update(

            kmer_statistics(sequence)

        )

        self.processed += 1

        return features

    # ---------------------------------------------------------

    def extract_all(self, reader):

        """
        Extracts features from every contig.
        """

        feature_matrix = []

        for contig in reader.stream_contigs():

            row = self.extract(contig)

            if row is None:

                continue

            feature_matrix.append(row)

            if self.processed % PROGRESS_INTERVAL == 0:

                progress(

                    self.processed,

                    max(reader.valid_contigs, 1)

                )

        info("")

        info("========== FEATURE EXTRACTION ==========")

        info(f"Processed : {self.processed:,}")

        info(f"Skipped   : {self.skipped:,}")

        info("========================================")

        return feature_matrix