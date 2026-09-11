"""LightGBM classifier for the metagenomic feature matrix."""

from __future__ import annotations

from typing import Any

from models.base_model import BaseClassificationModel


class LightGBMModel(BaseClassificationModel):
    """Train and evaluate a LightGBM classifier."""

    @property
    def model_name(self) -> str:
        """Return the filename-safe model name."""
        return "lightgbm"

    def _build_model(self) -> Any:
        """Create the configured LightGBM classifier."""
        try:
            from lightgbm import LGBMClassifier
        except ImportError as error:
            raise ImportError(
                "Install lightgbm to use LightGBMModel."
            ) from error

        return LGBMClassifier(
            n_estimators=300,
            learning_rate=0.05,
            num_leaves=31,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=self.random_state,
            n_jobs=-1,
            verbosity=-1,
        )


if __name__ == "__main__":
    LightGBMModel().run()
