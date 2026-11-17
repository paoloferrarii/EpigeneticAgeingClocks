# SCITAS — Izar (GPU, academic / courses)

Login: `ssh <user>@izar.hpc.epfl.ch`  
Project: `/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/`

**Do not run downloads / ComBat / clocks on the login node.**

Izar is for **students / teaching**. There is **no** `academic` partition here (that is Jed-only). Default partition: `gpu`.

## Data download (no GPU needed)

Use QoS `build` (allows 0 GPU):

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
mkdir -p logs
sbatch scripts/izar_select_download.sbatch
squeue -u "$USER"
```

Script uses: `--partition=gpu --qos=build`

## GPU job (clocks / torch later)

```bash
sbatch --partition=gpu --qos=normal --gres=gpu:1 --time=04:00:00 your_gpu_job.sbatch
```

## Interactive / Jupyter

`build` ≈ 1 job/user. If a download `sbatch` is running, interactive `-q build` hits `QOSMaxSubmitJobPerUserLimit` — wait, `scancel`, or use `debug`:

```bash
Sinteract -p gpu -q debug -g gpu:1 -c 4 -m 16G -t 01:00:00
# or when free: Sinteract -p gpu -q build -c 4 -m 16G -t 02:00:00
```

Then on the compute host: `jupyter lab --no-browser --ip=127.0.0.1 --port=8888`  
Laptop: `ssh -N -L 8888:127.0.0.1:8888 USER@<node> -J USER@izar.hpc.epfl.ch`  
See README for the full notebook steps.

## Docs

https://scitas-doc.epfl.ch/user-guide/using-clusters/slurm-qos-partitions/  
https://scitas-doc.epfl.ch/faq/
