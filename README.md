# Comparing Epigenetic Ageing Clocks (CS-433 ML4Science)

Fellay Lab · Barry Ryan · three-track layout so teammates can work independently.

## Tracks

| Track | Folder | Owner focus |
|-------|--------|-------------|
| A — Data | `01_data_download/` | GEO cohorts, selection, download, batch QC |
| B — Clocks | `02_pyaging/` | pyaging clocks |
| C — Analysis | (later) | metrics, clinical associations, figures |

## Quick start

```bash
# 1) Overview slide-3 selection (cached 7d)
python 01_data_download/select_datasets.py --force-refresh
python 01_data_download/select_datasets.py --require-age   # MAE subset

# 2) Download
python 01_data_download/geo_data_download.py --selected-only

# 3) Normalise: 01_data_download/normalise_data_candidates.ipynb
```

Config: `shared/paths.py` · Cache: `data/selection/selection_cache.json`  
Overview checks: n>100, clinical outcome, 450K/EPIC (CpG proxy). Beta verified at download.
