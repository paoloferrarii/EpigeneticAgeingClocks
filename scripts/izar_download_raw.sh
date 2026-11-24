#!/bin/bash
#SBATCH --job-name=download-raw
#SBATCH --partition=gpu
#SBATCH --qos=build
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out

python3 -u 01_data_download/geo_data_download_raw.py
