# SCITAS — Jed (CPU)

Login: `ssh <user>@jed.hpc.epfl.ch`  
Project: `/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/`

**Do not run downloads / ComBat / clocks on the login node.**

Best cluster for this project's data download (CPU-only).

## Students (academic partition)

```bash
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
mkdir -p logs
sbatch scripts/jed_select_download.sbatch
squeue -u "$USER"
```

Script uses: `--partition=academic --qos=serial`

## Lab accounts (if you have them)

```bash
sbatch --partition=standard --qos=serial scripts/jed_select_download.sbatch
```

## Interactive compute shell

```bash
Sinteract -p academic -q serial -c 4 -m 16G -t 02:00:00
cd /work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry
```

## Docs

https://scitas-doc.epfl.ch/user-guide/using-clusters/running-jobs/  
https://scitas-doc.epfl.ch/user-guide/using-clusters/slurm-qos-partitions/
