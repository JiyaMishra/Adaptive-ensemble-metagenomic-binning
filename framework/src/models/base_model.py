"""Shared training workflow for Phase 1 classification models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from config import FEATURE_MATRIX, RANDOM_SEED, TRAINED_MODELS_DIR
from logger import info


class BaseClassificationModel(ABC):
    """Base class for loading data, training, evaluating, and saving models."""

    LABEL_CANDIDATES = (
        "label",
        "target",
        "class",
        "bin",
        "bin_id",
        "cluster",
        "genome_id",
        "taxonomy",
    )
    IDENTIFIER_COLUMNS = {"contig_id", "id", "sequence_id", "header"}

    def __init__(
        self,
        label_column: Optional[str] = None,
        test_size: float = 0.2,
        random_state: int = RANDOM_SEED,
    ) -> None:
        """Initialize model settings shared by every classifier."""
        if not 0 < test_size < 1:
            raise ValueError("test_size must be between 0 and 1.")

        self.label_column = label_column
        self.test_size = test_size
        self.random_state = random_state
        self.feature_columns: List[str] = []
        self.model: Any = None
        self.label_encoder: Any = None

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Return the stable artifact name for the concrete model."""

    @abstractmethod
    def _build_model(self) -> Any:
        """Create the concrete classifier instance."""

    def load_data(
        self, file_path: Path = FEATURE_MATRIX
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """Load the matrix and identify features plus its label column."""
        if not file_path.exists():
            raise FileNotFoundError(f"Feature matrix not found: {file_path}")

        dataframe = pd.read_csv(file_path)
        if dataframe.empty:
            raise ValueError(f"Feature matrix is empty: {file_path}")

        label_column = self._resolve_label_column(dataframe)
        self.feature_columns = self._resolve_feature_columns(
            dataframe, label_column
        )
        if not self.feature_columns:
            raise ValueError("No numeric feature columns were detected.")

        features = dataframe[self.feature_columns].copy()
        features = features.fillna(
            features.median(numeric_only=True)
        ).fillna(0)
        labels = dataframe[label_column]

        if labels.isna().any():
            raise ValueError(
                f"Label column '{label_column}' contains missing values."
            )
        if labels.nunique() < 2:
            raise ValueError(
                "At least two unique classes are required for training."
            )

        self.label_column = label_column
        info(
            f"Loaded {len(dataframe):,} rows, "
            f"{len(self.feature_columns)} features, "
            f"and label column '{label_column}'."
        )
        return features, labels

    def train(
        self, features: pd.DataFrame, labels: pd.Series
    ) -> Dict[str, Any]:
        """Split, fit, evaluate, and save the classifier artifact."""
        train_test_split = self._import_train_test_split()
        stratify_labels = labels if labels.value_counts().min() >= 2 else None
        x_train, x_test, y_train, y_test = train_test_split(
            features,
            labels,
            test_size=self.test_size,
            random_state=self.random_state,
            stratify=stratify_labels,
        )

        label_encoder = self._import_label_encoder()
        self.label_encoder = label_encoder()
        encoded_y_train = self.label_encoder.fit_transform(y_train)
        self.model = self._build_model()
        info(f"Training {self.model_name} on {len(x_train):,} samples.")
        self.model.fit(x_train, encoded_y_train)

        encoded_predictions = self.model.predict(x_test)
        predictions = self.label_encoder.inverse_transform(encoded_predictions)
        probabilities = self.model.predict_proba(x_test)
        metrics = self._evaluate(y_test, predictions, probabilities)
        artifact_path = self.save_model()
        metrics["model_path"] = str(artifact_path)
        return metrics

    def run(self, file_path: Path = FEATURE_MATRIX) -> Dict[str, Any]:
        """Execute the complete Phase 1 workflow for one classifier."""
        features, labels = self.load_data(file_path)
        return self.train(features, labels)

    def save_model(self) -> Path:
        """Persist the fitted model and its input schema with joblib."""
        if self.model is None:
            raise RuntimeError("Train the model before saving it.")

        joblib = self._import_joblib()
        TRAINED_MODELS_DIR.mkdir(parents=True, exist_ok=True)
        artifact_path = TRAINED_MODELS_DIR / f"{self.model_name}.joblib"
        joblib.dump(
            {
                "model": self.model,
                "feature_columns": self.feature_columns,
                "label_column": self.label_column,
                "label_encoder": self.label_encoder,
            },
            artifact_path,
        )
        info(f"Saved {self.model_name} model to {artifact_path}.")
        return artifact_path

    def _resolve_label_column(self, dataframe: pd.DataFrame) -> str:
        """Return an explicit or conventionally named label column."""
        if self.label_column is not None:
            if self.label_column not in dataframe.columns:
                raise ValueError(
                    "Requested label column "
                    f"'{self.label_column}' was not found."
                )
            return self.label_column

        normalized_columns = {
            column.lower(): column for column in dataframe.columns
        }
        for candidate in self.LABEL_CANDIDATES:
            if candidate in normalized_columns:
                return normalized_columns[candidate]

        raise ValueError(
            "Could not detect a label column. "
            "Provide label_column explicitly. "
            f"Checked: {', '.join(self.LABEL_CANDIDATES)}."
        )

    def _resolve_feature_columns(
        self, dataframe: pd.DataFrame, label_column: str
    ) -> List[str]:
        """Select numeric columns while excluding labels and identifiers."""
        excluded = self.IDENTIFIER_COLUMNS | {label_column.lower()}
        return [
            column
            for column in dataframe.select_dtypes(include="number").columns
            if column.lower() not in excluded
        ]

    def _evaluate(
        self, labels: pd.Series, predictions: Any, probabilities: Any
    ) -> Dict[str, Any]:
        """Calculate and log standard weighted classification metrics."""
        accuracy_score, precision_recall_fscore_support = (
            self._import_metrics()
        )
        precision, recall, f1_score, _ = precision_recall_fscore_support(
            labels,
            predictions,
            average="weighted",
            zero_division=0,
        )
        metrics = {
            "accuracy": float(accuracy_score(labels, predictions)),
            "precision_weighted": float(precision),
            "recall_weighted": float(recall),
            "f1_weighted": float(f1_score),
            "test_samples": int(len(labels)),
            "probability_shape": tuple(probabilities.shape),
        }
        info(
            f"{self.model_name} evaluation | "
            f"accuracy={metrics['accuracy']:.4f}, "
            f"precision={metrics['precision_weighted']:.4f}, "
            f"recall={metrics['recall_weighted']:.4f}, "
            f"f1={metrics['f1_weighted']:.4f}"
        )
        return metrics

    @staticmethod
    def _import_train_test_split() -> Any:
        """Import scikit-learn splitting functionality on demand."""
        try:
            from sklearn.model_selection import train_test_split
        except ImportError as error:
            raise ImportError(
                "Install scikit-learn to train classification models."
            ) from error
        return train_test_split

    @staticmethod
    def _import_metrics() -> Tuple[Any, Any]:
        """Import scikit-learn metric functions on demand."""
        try:
            from sklearn.metrics import (
                accuracy_score,
                precision_recall_fscore_support,
            )
        except ImportError as error:
            raise ImportError(
                "Install scikit-learn to evaluate classification models."
            ) from error
        return accuracy_score, precision_recall_fscore_support

    @staticmethod
    def _import_label_encoder() -> Any:
        """Import scikit-learn label encoding on demand."""
        try:
            from sklearn.preprocessing import LabelEncoder
        except ImportError as error:
            raise ImportError(
                "Install scikit-learn to encode class labels."
            ) from error
        return LabelEncoder

    @staticmethod
    def _import_joblib() -> Any:
        """Import joblib on demand for model persistence."""
        try:
            import joblib
        except ImportError as error:
            raise ImportError(
                "Install joblib to save trained models."
            ) from error
        return joblib
