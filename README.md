# Comparing Epigenetic Ageing Clocks (CS-433 ML4Science)

Fellay Lab · Barry Ryan · three-track layout so teammates can work independently.

## Tracks

| Track | Folder | Owner focus |
|-------|--------|-------------|
| A — Data | `01_data_download/`, `03_data_qc/` | GEO cohorts, meta, batch QC, CpG overlap |
| B — Clocks | `02_pyaging/`, `04_clock_predictions/` | pyaging clocks, intermediates, validation |
| C — Analysis | `05_analysis/` | metrics, clinical associations, figures, narrative |

Shared contracts: [`contracts/SCHEMAS.md`](contracts/SCHEMAS.md) · paths: [`shared/paths.py`](shared/paths.py)

```mermaid
flowchart LR
  A[TrackA processed pickles] --> B[TrackB prediction CSVs]
  B --> C[TrackC metrics and figures]
```

## Quick start

```bash
export PY=/home/sommarug/miniconda3/envs/MLproject1/bin/python
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry

# Track B probe lists (feeds Track A overlap)
$PY 04_clock_predictions/export_probe_lists.py

# Track A — synthetic demos when EpiHydra/GEO unavailable
$PY 03_data_qc/make_synthetic_demo.py
$PY 03_data_qc/meta_standardize.py --datasets GSE_DEMO_A GSE_DEMO_B GSE_DEMO_C
$PY 03_data_qc/cpg_overlap.py --datasets GSE_DEMO_A GSE_DEMO_B GSE_DEMO_C
$PY 03_data_qc/batch_correction.py --datasets GSE_DEMO_A GSE_DEMO_B GSE_DEMO_C

# Real GEO (needs network + writable data root)
# EPICLOCK_DATA_ROOT=/path/to/data $PY 01_data_download/geo_data_download.py

# Track B
$PY 04_clock_predictions/ingest_tutorial_gse139307.py
$PY 04_clock_predictions/validate_clocks.py --clocks horvath2013 hannum dnamphenoage grimage2
$PY 04_clock_predictions/run_clocks.py --datasets GSE_DEMO_A GSE_DEMO_B GSE_DEMO_C GSE139307 \
  --clocks horvath2013 hannum dnamphenoage grimage2 --overwrite

# Track C
$PY 05_analysis/compute_metrics.py
$PY 05_analysis/missing_cpg_sensitivity.py --dataset GSE_DEMO_A
$PY 05_analysis/plot_comparisons.py
$PY 05_analysis/strengths_weaknesses.py
```

## Data roots

Default: project-local `data/`. Override with `EPICLOCK_DATA_ROOT` when Bryan's
EpiHydra tree is readable.

## Outcomes mapped

- Re-runnable dataset/clock DB → processed pickles + probe lists + scripts  
- Predictions per clock × dataset → `results/predictions/` (GrimAge/PhenoAge intermediates included)  
- Figures / metric tables / strengths–weaknesses → `results/figures/`, `results/metrics/`
