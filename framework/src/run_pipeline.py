"""
End-to-end execution runner for the Adaptive Ensemble Metagenomic Binning pipeline.
"""

import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
MODELS_DIR = SRC_DIR / "models"
if str(MODELS_DIR) not in sys.path:
    sys.path.insert(0, str(MODELS_DIR))
BINNERS_DIR = SRC_DIR / "binners"
if str(BINNERS_DIR) not in sys.path:
    sys.path.insert(0, str(BINNERS_DIR))
EVALUATION_DIR = SRC_DIR / "evaluation"
if str(EVALUATION_DIR) not in sys.path:
    sys.path.insert(0, str(EVALUATION_DIR))

EXPLAINABILITY_DIR = SRC_DIR / "explainability"
if str(EXPLAINABILITY_DIR) not in sys.path:
    sys.path.insert(0, str(EXPLAINABILITY_DIR))

import config
from build_feature_matrix import build_feature_matrix
from parse_metabat import parse_metabat
from parse_maxbin import parse_maxbin
from parse_vamb import parse_vamb
from ensemble_utils import build_ensemble_dataframe
import confidence_engine
import resource_engine
import adaptive_decision
import adaptive_ensemble
import adaptive_weight_calculator
import confidence_guided_refinement
import bin_refinement
import evaluate_ensemble
import quality_evaluator
import shap_analysis
import visualize_shap


def run_pipeline():
    print("=" * 60)
    print("STEP 1: Building Feature Matrix")
    print("=" * 60)
    build_feature_matrix()

    print("\n" + "=" * 60)
    print("STEP 2: Parsing Binner Outputs")
    print("=" * 60)
    parse_metabat(config.METABAT2_OUTPUT, config.METABAT2_OUTPUT / "metabat_assignments.csv")
    parse_maxbin(config.MAXBIN2_OUTPUT, config.MAXBIN2_OUTPUT / "maxbin_assignments.csv")
    parse_vamb(config.VAMB_OUTPUT / "vae_clusters_unsplit.tsv", config.VAMB_OUTPUT / "vamb_assignments.csv")

    print("\n" + "=" * 60)
    print("STEP 3: Merging Ensemble Dataframe")
    print("=" * 60)
    df_ensemble = build_ensemble_dataframe()
    print(f"Ensemble Dataframe built with {len(df_ensemble)} contigs.")

    print("\n" + "=" * 60)
    print("STEP 4: Confidence Engine")
    print("=" * 60)
    calibration = confidence_engine.main()
    print(
        "[SELF-CALIBRATION] Calibrated Thresholds: "
        f"Low < {calibration['low_threshold']:.4f}, "
        f"High > {calibration['high_threshold']:.4f}"
    )

    print("\n" + "=" * 60)
    print("STEP 4B: Hardware-Aware Execution Routing")
    print("=" * 60)
    capacity, routed_contigs = resource_engine.prepare_execution_routes()
    has_heavy_workset = (routed_contigs["execution_route"] == "HEAVY_ENSEMBLE").any()
    print(resource_engine.format_startup_telemetry(capacity, routed_contigs))

    print("\n" + "=" * 60)
    print("STEP 5: Adaptive Decisions")
    print("=" * 60)
    adaptive_decision.main()

    print("\n" + "=" * 60)
    print("STEP 6: Adaptive Ensemble Consensus")
    print("=" * 60)
    adaptive_ensemble.main()

    print("\n" + "=" * 60)
    print("STEP 7: Adaptive Weight Calculation")
    print("=" * 60)
    adaptive_weight_calculator.main()

    print("\n" + "=" * 60)
    print("STEP 8: Confidence-Guided Refinement")
    print("=" * 60)
    confidence_guided_refinement.main()

    print("\n" + "=" * 60)
    print("STEP 9: Feature Coherence Bin Refinement")
    print("=" * 60)
    bin_refinement.main()

    print("\n" + "=" * 60)
    print("STEP 10: Evaluating Ensemble")
    print("=" * 60)
    evaluate_ensemble.main()

    print("\n" + "=" * 60)
    print("STEP 11: Quality Evaluation")
    print("=" * 60)
    quality_evaluator.main()

    print("\n" + "=" * 60)
    print("STEP 12: SHAP Explainability Analysis")
    print("=" * 60)
    if has_heavy_workset:
        # KernelSHAP is the pipeline's most memory-intensive matrix operation.
        # Explain only routed heavy-path contigs while retaining all route metadata
        # in the canonical and execution-routes CSVs.
        shap_analysis.main(
            input_path=resource_engine.HEAVY_ENSEMBLE_OUTPUT,
            batch_size=capacity.chunk_size,
        )
    else:
        print("No HEAVY_ENSEMBLE contigs; skipping SHAP computation.")

    print("\n" + "=" * 60)
    print("STEP 13: SHAP Visualizations")
    print("=" * 60)
    if has_heavy_workset:
        visualize_shap.main()
    else:
        print("No HEAVY_ENSEMBLE contigs; skipping SHAP visualizations.")

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 60)


if __name__ == "__main__":
    run_pipeline()
