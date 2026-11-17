# Comparing Epigenetic Ageing Clocks (CS-433 ML4Science)

Fellay Lab · Barry Ryan · three-track layout so teammates can work independently.

## Tracks

| Track | Folder | Owner focus |
|-------|--------|-------------|
| A — Data | `01_data_download/` | GEO cohorts, selection, download, batch QC |
| B — Clocks | `02_pyaging/` | pyaging clocks |
| C — Analysis | (later) | metrics, clinical associations, figures |

## SCITAS

Project: `/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/`

**Never run downloads, ComBat, or clocks on the login node.** Use Slurm.

| Cluster | Who | Guide | Data job script |
|---------|-----|-------|-----------------|
| **Jed** (CPU) | students → `academic`; labs → `standard` | [docs/scitas/jed.md](docs/scitas/jed.md) | `scripts/jed_select_download.sbatch` |
| **Izar** (GPU, teaching) | students / courses | [docs/scitas/izar.md](docs/scitas/izar.md) | `scripts/izar_select_download.sbatch` |
| **Kuma** (GPU, research) | lab accounts | [docs/scitas/kuma.md](docs/scitas/kuma.md) | `scripts/kuma_select_download.sbatch` |

On **Izar** you already hit `invalid partition: academic` — that flag is Jed-only. Use the Izar script instead.

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
mkdir -p logs
sbatch scripts/izar_select_download.sbatch
squeue -u "$USER"
```

## After you are on a compute node / in a job

```bash
python 01_data_download/select_datasets.py --force-refresh
python 01_data_download/geo_data_download.py --selected-only
```

## Jupyter (compute node, not login)

`build` QoS ≈ **1 job/user**. If `sbatch` is already running, `Sinteract -q build` fails with `QOSMaxSubmitJobPerUserLimit` — wait, `scancel <JOBID>`, or use `debug`:

```bash
Sinteract -p gpu -q debug -g gpu:1 -c 4 -m 16G -t 01:00:00
# or, when no other build job: Sinteract -p gpu -q build -c 4 -m 16G -t 02:00:00
```

On the compute host (not `izar1`):

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
conda activate MLproject1
hostname   # e.g. i31
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

Laptop tunnel (replace `USER` / `i31`):

```bash
ssh -N -L 8888:127.0.0.1:8888 USER@i31 -J USER@izar.hpc.epfl.ch
```

Browser: `http://127.0.0.1:8888/lab?token=...` → open `01_data_download/normalise_data_candidates.ipynb`.

Jed/Kuma interactive flags: [jed.md](docs/scitas/jed.md) · [kuma.md](docs/scitas/kuma.md).

Config: `shared/paths.py` · Cache: `data/selection/selection_cache.json`  
Overview: n>100, clinical outcome, 450K/EPIC. Beta checked at download.
