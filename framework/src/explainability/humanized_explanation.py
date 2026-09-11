"""Turn final ensemble and SHAP evidence into reviewer-friendly explanations."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTPUT = PROJECT_ROOT / "framework/results/explanations/humanized_explanations.csv"


def _value(row: pd.Series, name: str, default: str = "not available") -> str:
    value = row.get(name)
    return default if pd.isna(value) or value is None or value == "" else str(value)


def explanation_for_row(row: pd.Series, weight_summary: str = "") -> str:
    """Create a paragraph solely from evidence columns present in ``row``."""

    contig = _value(row, "contig_id", "This contig")
    method = _value(row, "selected_method")
    selected_bin = _value(row, "selected_bin", "no bin")
    category = _value(row, "confidence_category").lower()
    confidence_score = _value(row, "confidence_score")
    route = _value(row, "execution_route").replace("_", " ").lower()
    features = _value(row, "top_shap_features")
    agreement = pd.to_numeric(pd.Series([row.get("tool_agreement_score")]), errors="coerce").iloc[0]
    biology = pd.to_numeric(pd.Series([row.get("biological_consistency_score")]), errors="coerce").iloc[0]

    text = f"Contig {contig} was assigned to {method}, bin {selected_bin}, with {category} confidence (score {confidence_score}) via the {route} route."
    if features != "not available":
        text += f" The strongest SHAP contributors were {features}."
    evidence = []
    if pd.notna(agreement):
        evidence.append("strong" if agreement >= 0.67 else "limited" if agreement < 0.34 else "moderate")
    if pd.notna(biology):
        evidence.append("strong" if biology >= 0.67 else "limited" if biology < 0.34 else "moderate")
    if evidence:
        text += f" Tool agreement and biological consistency were {', '.join(evidence)}, so this assignment is treated as {('well supported' if category == 'high' else 'cautious')}."
    if weight_summary:
        text += f" {weight_summary}"
    return text


def build_humanized_explanations(
    decisions: pd.DataFrame,
    shap_features: pd.DataFrame | None = None,
    feature_feedback: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return the required explanation table without inventing absent evidence."""
    df = decisions.copy()
    if shap_features is not None and not shap_features.empty:
        df = df.merge(shap_features[["contig_id", "top_shap_features"]], on="contig_id", how="left")
    elif "top_shap_features" not in df:
        df["top_shap_features"] = pd.NA

    required = [
        "contig_id", "selected_method", "selected_bin", "confidence_score",
        "confidence_category", "execution_route", "tool_agreement_score",
        "biological_consistency_score", "decision_reason", "top_shap_features",
    ]
    for column in required:
        if column not in df:
            df[column] = pd.NA
    weight_summary = ""
    if feature_feedback is not None and not feature_feedback.empty:
        changed = feature_feedback[feature_feedback["weight_change"].abs() > 1e-12]
        if changed.empty:
            weight_summary = "Feature weights did not change for the next batch."
        else:
            changes = ", ".join(
                f"{row.feature} {row.previous_weight:.3f} to {row.new_weight:.3f}"
                for row in changed.itertuples()
            )
            weight_summary = (
                f"Feature weights changed for the next batch ({changes}); "
                "the adaptive ensemble will use those bounded updates."
            )
    df["humanized_explanation"] = df.apply(
        lambda row: explanation_for_row(row, weight_summary), axis=1
    )
    return df[required + ["humanized_explanation"]]


def save_humanized_explanations(
    decisions: pd.DataFrame,
    shap_features: pd.DataFrame | None = None,
    output: Path = OUTPUT,
    feature_feedback: pd.DataFrame | None = None,
) -> pd.DataFrame:
    result = build_humanized_explanations(decisions, shap_features, feature_feedback)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    return result
