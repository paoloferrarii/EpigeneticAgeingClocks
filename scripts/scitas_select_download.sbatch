#!/bin/bash
#SBATCH --job-name=epiclock-data
#SBATCH --partition=academic
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry}"
cd "$PROJECT_DIR"
mkdir -p logs

# Adjust to your SCITAS conda/env if needed:
# module load gcc python
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate MLproject1

export EPICLOCK_DATA_ROOT="${EPICLOCK_DATA_ROOT:-$PROJECT_DIR/data}"

python 01_data_download/select_datasets.py --force-refresh
python 01_data_download/geo_data_download.py --selected-only
