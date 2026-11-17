"""Paths, catalog, Overview thresholds."""
from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("EPICLOCK_DATA_ROOT", str(PROJECT_ROOT / "data")))

RAW_DIR = DATA_ROOT / "01_raw"
PROCESSED_DIR = DATA_ROOT / "02_processed"
SELECTION_DIR = DATA_ROOT / "selection"
CACHE_FILE = SELECTION_DIR / "selection_cache.json"
TTL_DAYS = float(os.environ.get("EPICLOCK_SELECTION_TTL_DAYS", "7"))

GSE_LIST = [
    "GSE27097", "GSE40279", "GSE41037", "GSE42861", "GSE84727",
    "GSE111629", "GSE157131", "GSE197674", "GSE87571", "GSE55763",
]

# Overview slide 3
MIN_N = 100
PLATFORM_OK = {"GPL13534", "GPL21145", "GPL23976"}  # 450K/EPIC ≈ CpG overlap >50%
OUTCOME_WORDS = (
    "disease", "diagnosis", "status", "mortality", "death", "survival",
    "cancer", "parkinson", "alzheimer", "schizophrenia", "arthritis",
    "ptsd", "trauma", "depression", "hypertension", "blood pressure",
    "systolic", "diastolic", "smoking", "diabetes", "case", "control",
)
