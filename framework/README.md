# Adaptive Ensemble Metagenomic Binning

## Run the pipeline

From the repository root, create and activate a Python environment, install the
Python dependencies, ensure MetaBAT2, MaxBin2, and VAMB are installed on your
`PATH`, then run:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r framework/requirements.txt
python framework/src/binners/run_metabat2.py
python framework/src/binners/run_maxbin2.py
python framework/src/binners/run_vamb.py
python framework/src/run_pipeline.py
```

The pipeline uses the real assembly at `data/assemblies/metagem_1500/` and its
matching depth file.  It expects binner output directories under
`framework/results/`; run the provided binner wrappers first when those outputs
are absent.  On completion, inspect:

- `framework/results/explanations/humanized_explanations.csv`
- `framework/results/explanations/feature_feedback.csv`
- `framework/results/explanations/feedback_history.csv`
- `framework/results/visualizations/`
- `framework/results/evaluation/closed_loop_comparison.csv`

## XAI layers

SHAP is the technical explanation of the features that influenced adaptive
weighting.  Humanized XAI converts that evidence and the final confidence/route
into a concise explanation for reviewers.  Closed-loop XAI uses SHAP together
with confidence, tool agreement, biological consistency, and the actual
`FAST_PATH`/`HEAVY_ENSEMBLE` execution route to conservatively adjust feature
weights for the next batch.  It does not treat a large SHAP value as bad by
itself, and it does not claim ground-truth accuracy without labels.
