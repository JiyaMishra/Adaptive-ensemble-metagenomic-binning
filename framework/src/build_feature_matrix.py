"""
=============================================================
Feature Matrix Builder
Adaptive Explainable Ensemble Framework
-------------------------------------------------------------
Runs the complete feature extraction pipeline and
exports the feature matrix to CSV.
=============================================================
"""

from pathlib import Path
import pandas as pd

from config import (
    FEATURE_DIR,
    FEATURE_MATRIX
)

from logger import info

from readfasta import FASTAReader
from extract_features import FeatureExtractor


class FeatureMatrixBuilder:

    def __init__(self):

        FEATURE_DIR.mkdir(parents=True, exist_ok=True)

        self.reader = FASTAReader()

        self.extractor = FeatureExtractor()

    # ---------------------------------------------------------

    def build(self):

        info("")
        info("========== BUILDING FEATURE MATRIX ==========")

        rows = []

        for contig in self.reader.stream_contigs():

            features = self.extractor.extract(contig)

            if features is not None:

                rows.append(features)

        self.reader.summary()

        df = pd.DataFrame.from_records(rows)

        df.to_csv(FEATURE_MATRIX, index=False)

        info("")
        info("Feature matrix saved successfully!")

        info(f"Location : {FEATURE_MATRIX}")

        info(f"Rows      : {len(df):,}")

        info(f"Columns   : {len(df.columns)}")

        info("============================================")

        return df

    # ---------------------------------------------------------

    def summary(self, dataframe):

        info("")

        info("========== DATASET SUMMARY ==========")

        info(f"Shape : {dataframe.shape}")

        info("")

        info("Missing values")

        missing = dataframe.isnull().sum()

        for column, value in missing.items():

            if value > 0:

                info(f"{column} : {value}")

        info("")

        info("Preview")

        print(dataframe.head())

        info("======================================")


# =============================================================
# MAIN
# =============================================================

if __name__ == "__main__":

    builder = FeatureMatrixBuilder()

    df = builder.build()

    builder.summary(df)