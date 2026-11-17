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
python 01_data_download/geo_data_download.py --selected-only --require-age
```

## Jupyter via VS Code Remote (worked on Izar)

Jupyter runs on the **compute node**. VS Code Remote stays on the **login** node only to forward the port. Do not run kernels on `izar1`.

`build` ≈ **1 job/user**. If a download `sbatch` is running, wait / `scancel` / use `debug` before `Sinteract -q build`.

### Terminal A — compute (prompt like `USER@i31`, not `izar1`)

```bash
Sinteract -p gpu -q build -c 4 -m 16G -t 02:00:00
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
conda activate MLproject1
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

Copy the `token=...` from the log. Leave this running.

### Terminal B — login (`izar1`)

Replace `i31` with your `hostname` from Terminal A. Bind **IPv4** only (plain `-L 1889:...` can fail with `bind [::1]: Cannot assign requested address`):

```bash
ssh -N -L 127.0.0.1:1889:127.0.0.1:8888 i31
```

Leave this running (no prompt = OK).

### VS Code Ports → browser

1. **View → Ports** → **Forward a Port** → `1889`
2. Open in Browser:  
   `http://127.0.0.1:1889/lab?token=<paste-token>`
3. Open `01_data_download/normalise_data_candidates.ipynb`

Optional check on `izar1`: `curl -I --max-time 3 http://127.0.0.1:1889/lab`

Avoid forwarding **8888/8899** on the laptop (Cursor often owns them). `1889` is fine.

Jed/Kuma interactive flags: [jed.md](docs/scitas/jed.md) · [kuma.md](docs/scitas/kuma.md).

Config: `shared/paths.py` · Cache: `data/selection/selection_cache.json`  
Overview: n>100, clinical outcome, 450K/EPIC, prefer age (`--require-age`). Beta checked at download.
