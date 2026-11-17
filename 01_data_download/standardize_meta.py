#!/usr/bin/env python3
"""
standardize_meta.py

Carica i pickle prodotti da standardize_cpg.py (OUT_DIR/<GSE>/methylation.pkl).
Per ogni attributo in ATTRIBUTES cerca le colonne dei metadati con i suoi regex,
le stampa con alcuni esempi e chiede da terminale il numero della colonna giusta.
Le colonne scelte sostituiscono "meta" nel pickle e vengono salvate come "std_meta".
"""
import pickle
import re
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# COSTANTI
# --------------------------------------------------------------------------
DATA_DIR = Path("/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/data")
STD_DIR = DATA_DIR / "02_standardized"           # pickle prodotti da standardize_cpg.py
PICKLE_NAME = "methylation.pkl"

GSE_LIST = [
    "GSE41037", "GSE40279", "GSE42861", "GSE27097",
    "GSE84727", "GSE111629", "GSE157131", "GSE197674",
]

N_EXAMPLES = 5   # valori di esempio stampati per ogni colonna

# attributo standard : regex da cercare nei nomi delle colonne
ATTRIBUTES = {
    "age":       [r"^age", r"\bage\b"],
    "sex":       [r"sex", r"gender"],
    "disease":   [r"disease", r"diagnos", r"status", r"condition", r"phenotype", r"subject", r"group"],
    "smoking":   [r"smok", r"tobacco"],
    "ethnicity": [r"race", r"ethnic", r"ancestr"],
    "tissue":    [r"tissue", r"cell.?type", r"source"],
}


# --------------------------------------------------------------------------
# FUNZIONI
# --------------------------------------------------------------------------
def find_columns(meta: pd.DataFrame, patterns: list) -> list:
    """Colonne il cui nome corrisponde ad almeno un regex (maiuscole/minuscole ignorate)."""
    return [c for c in meta.columns
            if any(re.search(p, c, re.IGNORECASE) for p in patterns)]


def print_columns(meta: pd.DataFrame, columns: list) -> None:
    """Stampa le colonne numerate con alcuni valori di esempio."""
    for i, col in enumerate(columns, 1):
        examples = meta[col].dropna().unique()[:N_EXAMPLES].tolist()
        print(f"  {i}) {col}: {examples}")


def ask_column(columns: list):
    """Chiede il numero della colonna. 0 = nessuna."""
    while True:
        answer = input(f"  Numero della colonna [1-{len(columns)}, 0 = nessuna]: ").strip()
        if answer.isdigit() and 0 <= int(answer) <= len(columns):
            return None if answer == "0" else columns[int(answer) - 1]
        print("  Numero non valido.")


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main():
    for gse in GSE_LIST:
        path = STD_DIR / gse / PICKLE_NAME
        with open(path, "rb") as f:
            data = pickle.load(f)
        meta = data["meta"]
        print(f"\n===== {gse} ({len(meta)} campioni) =====")

        std_meta = pd.DataFrame({"ID": meta["ID"]})
        for attr, patterns in ATTRIBUTES.items():
            print(f"\n[{attr}]")
            columns = find_columns(meta, patterns)
            if not columns:
                print("  nessuna colonna trovata")
                std_meta[attr] = None
                continue
            print_columns(meta, columns)
            col = ask_column(columns)
            std_meta[attr] = meta[col].values if col is not None else None

        del data["meta"]
        data["std_meta"] = std_meta
        with open(path, "wb") as f:
            pickle.dump(data, f)
        print(f"\n{gse}: std_meta salvato in {path}")


if __name__ == "__main__":
    main()