# Comparing Epigenetic Ageing Clocks (CS-433 ML4Science)

Fellay Lab · Barry Ryan · three-track layout so teammates can work independently.

## Tracks

| Track | Folder | Owner focus |
|-------|--------|-------------|
| A — Data | `01_data_download/` | GEO cohorts, selection, download, batch QC |
| B — Clocks | `02_pyaging/` | pyaging clocks |
| C — Analysis | (later) | metrics, clinical associations, figures |

## SCITAS (required)

Project path: `/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/`

**Do not run downloads, ComBat, or clocks on the login node.** Login is only for `git`, editing, and submitting jobs.

| Where | Allowed |
|-------|---------|
| Login node | `git`, `sbatch`, `squeue`, `scancel`, light `ls`/`less` |
| Compute (Slurm) | `select_datasets.py`, `geo_data_download.py`, normalise, pyaging |

Submit the data pipeline:

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
mkdir -p logs
sbatch scripts/scitas_select_download.sbatch
squeue -u "$USER"
```

Interactive compute shell (still not the login node):

```bash
Sinteract -c 4 -m 16G -t 02:00:00
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
```

Optional: `export EPICLOCK_DATA_ROOT=/path/to/data` (default: `$PROJECT/data`).

Students on Jed: script uses `#SBATCH --partition=academic --qos=serial`. Change partition/QoS if your account differs (see [SCITAS docs](https://scitas-doc.epfl.ch/user-guide/using-clusters/running-jobs/)).

## Local / after you have a compute shell

```bash
python 01_data_download/select_datasets.py --force-refresh
python 01_data_download/select_datasets.py --require-age
python 01_data_download/geo_data_download.py --selected-only
```

Then open `01_data_download/normalise_data_candidates.ipynb` (on a compute session, not login).

Config: `shared/paths.py` · Cache: `data/selection/selection_cache.json`  
Overview checks: n>100, clinical outcome, 450K/EPIC (CpG proxy). Beta verified at download.
