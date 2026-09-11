"""Confidence-aware SHAP feedback for the next adaptive-ensemble batch.

For feature f, a batch score is the SHAP-magnitude-weighted mean of
``confidence_signal + agreement_signal + biological_signal + route_signal``.
HIGH/MEDIUM/LOW confidence map to +1/0/-1; agreement and biological evidence
are centered at 0.5; FAST_PATH is +0.25 and HEAVY_ENSEMBLE is -0.25.  The
bounded update is ``w_new=clip(w_old*(1+learning_rate*score), min, max)``.
Thus a large SHAP magnitude only gives evidence more influence; it never alone
causes a feature penalty.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parents[1]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
import config

MODELS_DIR = SRC_DIR / "models"
if str(MODELS_DIR) not in sys.path:
    sys.path.insert(0, str(MODELS_DIR))
import adaptive_ensemble

from humanized_explanation import save_humanized_explanations


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS = PROJECT_ROOT / "framework/results"
EXPLANATIONS = RESULTS / "explanations"
VISUALIZATIONS = RESULTS / "visualizations"
EVALUATION = RESULTS / "evaluation"


def shap_feature_table(shap_values: pd.DataFrame) -> pd.DataFrame:
    """Collapse multi-output SHAP columns into one magnitude per real feature."""
    result = pd.DataFrame({"contig_id": shap_values["contig_id"]})
    for feature in ("length", "coverage", "entropy"):
        columns = [c for c in shap_values if c.endswith(f"_shap_{feature}")]
        if columns:
            result[f"shap_{feature}"] = shap_values[columns].abs().mean(axis=1)
    return result


def _top_features(feature_table: pd.DataFrame) -> pd.DataFrame:
    feature_columns = [c for c in feature_table if c.startswith("shap_")]
    out = feature_table[["contig_id"]].copy()
    out["top_shap_features"] = feature_table[feature_columns].apply(
        lambda row: ", ".join(
            name.removeprefix("shap_")
            for name in row.sort_values(ascending=False).head(3).index
        ),
        axis=1,
    )
    return out


def _route_signal(value: object) -> float:
    # These are the actual routes emitted by resource_engine.py.
    return {"FAST_PATH": 0.25, "HEAVY_ENSEMBLE": -0.25}.get(str(value), 0.0)


def update_feature_weights(batch: pd.DataFrame, current_weights: dict[str, float]) -> pd.DataFrame:
    """Calculate conservative, evidence-aware feature updates for one batch."""
    feature_columns = [
        feature for feature in ("length", "coverage", "entropy")
        if f"shap_{feature}" in batch
    ]
    if not feature_columns:
        raise ValueError("SHAP feedback requires at least one recognized SHAP feature.")
    category_signal = batch["confidence_category"].map(
        {"HIGH": 1.0, "MEDIUM": 0.0, "LOW": -1.0}
    ).fillna(0.0)
    confidence_score = pd.to_numeric(batch["confidence_score"], errors="coerce").fillna(0.4)
    score_signal = ((confidence_score - 0.4) / 0.2).clip(-1.0, 1.0)
    confidence_signal = 0.5 * category_signal + 0.5 * score_signal
    agreement_signal = 2 * pd.to_numeric(batch["tool_agreement_score"], errors="coerce").fillna(0.5) - 1
    biological_signal = 2 * pd.to_numeric(batch["biological_consistency_score"], errors="coerce").fillna(0.5) - 1
    route_signal = batch["execution_route"].map(_route_signal).fillna(0.0)
    combined = confidence_signal + agreement_signal + biological_signal + route_signal
    rows = []
    for feature in feature_columns:
        magnitude = pd.to_numeric(batch[f"shap_{feature}"], errors="coerce").fillna(0.0)
        total = magnitude.sum()
        weights = magnitude / total if total > 0 else pd.Series(1 / len(batch), index=batch.index)
        score = float((weights * combined).sum())
        previous = float(current_weights.get(feature, 1.0))
        new = float(np.clip(previous * (1 + config.XAI_FEEDBACK_LEARNING_RATE * score), config.XAI_MINIMUM_FEATURE_WEIGHT, config.XAI_MAXIMUM_FEATURE_WEIGHT))
        low_association = float(weights[batch["confidence_category"].eq("LOW")].sum())
        mean_agreement = float((weights * pd.to_numeric(batch["tool_agreement_score"], errors="coerce").fillna(0.0)).sum())
        mean_biology = float((weights * pd.to_numeric(batch["biological_consistency_score"], errors="coerce").fillna(0.0)).sum())
        reason = "reliable" if score > 0.05 else "low-confidence association" if score < -0.05 else "mixed evidence"
        rows.append({"feature": feature, "previous_weight": previous, "mean_abs_shap": float(magnitude.mean()), "low_confidence_association": low_association, "mean_tool_agreement": mean_agreement, "mean_biological_consistency": mean_biology, "feedback_score": score, "new_weight": new, "weight_change": new - previous, "feedback_reason": reason,
                     "confidence_signal": float((weights * confidence_signal).sum()), "agreement_signal": float((weights * agreement_signal).sum()), "biological_signal": float((weights * biological_signal).sum())})
    return pd.DataFrame(rows)


def _plot_weight_evolution(history: pd.DataFrame) -> Path:
    VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    path = VISUALIZATIONS / "feature_weight_evolution.png"
    plt.figure(figsize=(8, 5))
    for feature, values in history.groupby("feature"):
        plt.plot(values["batch_id"], values["new_weight"], marker="o", label=feature)
    plt.axhline(1.0, color="grey", linewidth=0.8, linestyle="--")
    plt.xlabel("Batch")
    plt.ylabel("Feature weight for next batch")
    plt.title("Closed-loop XAI feature-weight evolution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()
    return path


def run_closed_loop(decisions: pd.DataFrame, shap_values: pd.DataFrame, batch_size: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Process batches in order and return feedback, history, explanations."""
    features = shap_feature_table(shap_values)
    data = decisions.merge(features, on="contig_id", how="inner", validate="one_to_one")
    if data.empty:
        raise ValueError("No contigs overlap between final decisions and SHAP values.")
    batch_size = int(batch_size or config.XAI_BATCH_SIZE)
    weights = {feature.removeprefix("shap_"): 1.0 for feature in features.columns if feature.startswith("shap_")}
    all_feedback, history = [], []
    for start in range(0, len(data), batch_size):
        batch_id = start // batch_size + 1
        # Before processing a later batch, run the real adaptive ensemble with
        # the feature weights learned from its predecessor.  The persisted
        # output makes the closed-loop hand-off auditable rather than merely a
        # feedback CSV.
        if batch_id > 1:
            batch_dir = EXPLANATIONS / "batches"
            batch_dir.mkdir(parents=True, exist_ok=True)
            batch_input = batch_dir / f"batch_{batch_id}_ensemble_input.csv"
            batch_output = batch_dir / f"batch_{batch_id}_adaptive_ensemble.csv"
            decision_columns = [
                "contig_id", "metabat_bin", "maxbin_bin", "vamb_bin",
                "length", "coverage", "entropy",
            ]
            available = [column for column in decision_columns if column in data]
            data.iloc[start:start + batch_size][available].to_csv(batch_input, index=False)
            adaptive_ensemble.main(
                feature_weights=weights,
                input_path=batch_input,
                output_path=batch_output,
            )
        feedback = update_feature_weights(data.iloc[start:start + batch_size], weights)
        feedback["batch_id"] = batch_id
        all_feedback.append(feedback)
        history.append(feedback[["batch_id", "feature", "previous_weight", "new_weight", "weight_change", "mean_abs_shap", "confidence_signal", "agreement_signal", "biological_signal", "feedback_score", "feedback_reason"]])
        weights.update(feedback.set_index("feature")["new_weight"].to_dict())
    feedback_all = pd.concat(all_feedback, ignore_index=True)
    history_all = pd.concat(history, ignore_index=True)
    final_feedback = feedback_all.sort_values("batch_id").groupby("feature", as_index=False).tail(1)
    EXPLANATIONS.mkdir(parents=True, exist_ok=True)
    final_feedback.drop(columns="batch_id").to_csv(EXPLANATIONS / "feature_feedback.csv", index=False)
    history_all.to_csv(EXPLANATIONS / "feedback_history.csv", index=False)
    _plot_weight_evolution(history_all)
    explanations = save_humanized_explanations(
        decisions,
        _top_features(features),
        feature_feedback=final_feedback,
    )
    return final_feedback, history_all, explanations


def write_comparison(baseline: pd.DataFrame, closed_loop: pd.DataFrame) -> pd.DataFrame:
    """Report observed internal metrics only; this is not ground-truth accuracy."""
    rows = []
    for name, df in (("baseline", baseline), ("closed_loop", closed_loop)):
        rows.append({"mode": name, "assignment_rate": float(df["selected_bin"].notna().mean()), "mean_confidence": float(df["confidence_score"].mean()), "mean_tool_agreement": float(df["tool_agreement_score"].mean()), "mean_biological_consistency": float(df["biological_consistency_score"].mean()), "ground_truth_accuracy_measured": False})
    result = pd.DataFrame(rows)
    EVALUATION.mkdir(parents=True, exist_ok=True)
    result.to_csv(EVALUATION / "closed_loop_comparison.csv", index=False)
    return result
