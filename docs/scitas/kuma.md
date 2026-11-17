# SCITAS — Kuma (GPU, research / labs)

Login: `ssh <user>@kuma.hpc.epfl.ch`  
Project: `/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/`

**Do not run downloads / ComBat / clocks on the login node.**

Kuma is for **lab / research** accounts (not the default student teaching cluster).  
**You must set a partition** (`h100`, `l40s`, or `mig`) — there is no default.

## Data download (no GPU needed)

Use QoS `build` (0 GPU, short CPU work):

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
mkdir -p logs
sbatch scripts/kuma_select_download.sbatch
squeue -u "$USER"
```

Script uses: `--partition=l40s --qos=build`

## GPU job (clocks / torch)

H100:

```bash
sbatch --partition=h100 --qos=normal --gres=gpu:1 your_gpu_job.sbatch
```

L40S:

```bash
sbatch --partition=l40s --qos=normal --gres=gpu:1 your_gpu_job.sbatch
```

## Interactive

```bash
Sinteract -p l40s -q build -c 4 -m 16G -t 02:00:00
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
```

With a GPU:

```bash
Sinteract -p l40s -q debug -g gpu:1 -c 8 -m 32G -t 01:00:00
```

## Docs

https://www.epfl.ch/research/facilities/scitas/hardware/kuma/  
https://scitas-doc.epfl.ch/user-guide/using-clusters/slurm-qos-partitions/
