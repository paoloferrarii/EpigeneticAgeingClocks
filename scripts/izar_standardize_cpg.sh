#!/bin/bash
#SBATCH --job-name=std-cpg
#SBATCH --partition=gpu
#SBATCH --qos=build
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate epiclock
python -u 01_data_download/standardize_cpg.py
