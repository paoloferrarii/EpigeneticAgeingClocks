#!/usr/bin/env python3
"""
standardize_cpg.py

Legge le matrici grezze già scaricate (RAW_DIR/<GSE>/<GSE>_series_matrix.txt.gz),
converte gli M-value in beta, elimina CpG e campioni sotto la soglia di presenza,
imputa i valori mancanti e salva i pickle in OUT_DIR/<GSE>/methylation.pkl
con {"meta": metadati grezzi dall'header, "CpG": matrice}.
"""
import gzip
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# COSTANTI
# --------------------------------------------------------------------------
DATA_DIR = Path("/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/data")
RAW_DIR = DATA_DIR / "01_raw"                    # matrici grezze già scaricate
OUT_DIR = DATA_DIR / "02_standardized"           # cartella affianco per i nuovi pickle
SERIES_MATRIX_NAME = "{gse}_series_matrix.txt.gz"
PICKLE_NAME = "methylation.pkl"

GSE_LIST = [
    "GSE41037", "GSE40279", "GSE42861", "GSE27097",
    "GSE84727", "GSE111629", "GSE157131", "GSE197674",
]

PRESENCE_THRESHOLD = 0.9   # frazione minima di valori presenti per tenere una CpG o un campione


# --------------------------------------------------------------------------
# FUNZIONI
# --------------------------------------------------------------------------
def load_metadata(path: Path) -> pd.DataFrame:
    """Legge i metadati dei campioni dall'header del Series Matrix: una riga per campione.

    Ogni cella "campo: valore" delle righe !Sample_characteristics_ch1 diventa una colonna "campo".
    """
    meta = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            values = [v.strip('"') for v in values]
            if key == "!Sample_geo_accession":
                meta["ID"] = values
            elif key == "!Sample_title":
                meta["title"] = values
            elif key == "!Sample_source_name_ch1":
                meta["source"] = values
            elif key.startswith("!Sample_characteristics"):
                for i, cell in enumerate(values):
                    if ":" in cell:
                        field, value = cell.split(":", 1)
                        meta.setdefault(field.strip(), [None] * len(values))[i] = value.strip()
    return pd.DataFrame(meta)


def load_series_matrix(path: Path) -> pd.DataFrame:
    """Legge la tabella CpG x campioni dal Series Matrix .gz (salta le righe '!' di header)."""
    return pd.read_csv(path, sep="\t", comment="!", index_col=0, compression="gzip")


def mvalue_to_beta(cpg: pd.DataFrame) -> pd.DataFrame:
    """Converte M-value in beta: beta = 2^M / (1 + 2^M)."""
    m = np.exp2(cpg)
    return m / (1 + m)


def remove_below_threshold(cpg: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Elimina le CpG (righe) e poi i campioni (colonne) con frazione di valori presenti < threshold."""
    cpg = cpg.loc[cpg.notna().mean(axis=1) >= threshold]
    cpg = cpg.loc[:, cpg.notna().mean(axis=0) >= threshold]
    return cpg


def impute_median(cpg: pd.DataFrame) -> pd.DataFrame:
    """Sostituisce i valori mancanti con la mediana della stessa CpG su tutti i campioni."""
    values = cpg.to_numpy(dtype=float, copy=True)
    missing = np.isnan(values)
    rows = missing.any(axis=1)                   # median only for CpGs with missing values
    if rows.any():
        sub = values[rows]
        medians = np.nanmedian(sub, axis=1)
        values[rows] = np.where(missing[rows], medians[:, None], sub)
    return pd.DataFrame(values, index=cpg.index, columns=cpg.columns)


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main():
    for gse in GSE_LIST:
        path = RAW_DIR / gse / SERIES_MATRIX_NAME.format(gse=gse)
        cpg = load_series_matrix(path)
        meta = load_metadata(path)

        # M-value se i valori escono da [0, 1] (i beta non sono mai negativi)
        if np.nanpercentile(cpg.values, 1) < -0.5:
            cpg = mvalue_to_beta(cpg)

        cpg = remove_below_threshold(cpg, PRESENCE_THRESHOLD)
        cpg = impute_median(cpg)
        meta = meta[meta["ID"].isin(cpg.columns)]       # solo i campioni rimasti

        out_path = OUT_DIR / gse / PICKLE_NAME
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "wb") as f:
            pickle.dump({"meta": meta, "CpG": cpg}, f)
        print(f"{gse}: {cpg.shape[0]} CpG x {cpg.shape[1]} campioni -> {out_path}")


if __name__ == "__main__":
    main()