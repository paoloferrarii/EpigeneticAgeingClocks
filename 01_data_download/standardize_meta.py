#!/usr/bin/env python3
"""
standardize_meta.py

Carica i pickle prodotti da standardize_cpg.py (OUT_DIR/<GSE>/methylation.pkl).
Per ogni attributo in ATTRIBUTES cerca le colonne dei metadati con i suoi regex,
le stampa con alcuni esempi e chiede da terminale il numero della colonna giusta.
Poi stampa i valori della colonna originale e, sotto, i valori standard numerati
(STANDARD_VALUES): da terminale si scrive un numero per ogni valore originale, nello
stesso ordine. Gli attributi numerici (NUMERIC) vengono solo convertiti in numero.
Le colonne scelte vengono salvate nel pickle come "std_meta", con il nome standard
dell'attributo; "meta" resta invariato. I dataset che hanno già "std_meta" vengono
saltati (vedi OVERWRITE), quindi lo script si può interrompere e rilanciare.
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

OVERWRITE = False   # True = rifà la scelta anche per i dataset che hanno già std_meta
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

NUMERIC = ["age"]   # attributi da convertire in numero, senza valori standard

# attributo standard : valori ammessi (0 da terminale = valore mancante)
STANDARD_VALUES = {
    "sex":       ["male", "female"],
    "disease":   ["control", "schizophrenia", "rheumatoid arthritis",
                  "parkinson", "cancer", "hypertension", "other"],
    "smoking":   ["never", "former", "current", "occasional"],
    "ethnicity": ["european", "african", "asian", "hispanic", "other"],
    "tissue":    ["whole blood", "pbmc", "leukocytes", "other"],
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


def ask_values(values: list, standard: list) -> dict:
    """Stampa valori originali e standard, chiede un numero per valore. 0 = mancante."""
    print("  Valori originali (in ordine):")
    for i, value in enumerate(values, 1):
        print(f"    {i}. {value!r}")
    print("  Valori standard:")
    print("    0) mancante")
    for i, name in enumerate(standard, 1):
        print(f"    {i}) {name}")
    while True:
        answer = input(f"  Scrivi {len(values)} numeri [0-{len(standard)}] separati da spazio, "
                       "uno per valore originale: ").split()
        if len(answer) != len(values):
            print(f"  Servono {len(values)} numeri, ne hai scritti {len(answer)}.")
            continue
        if not all(a.isdecimal() and int(a) <= len(standard) for a in answer):
            print("  Numero non valido.")
            continue
        mapping = {v: standard[int(a) - 1] if int(a) else None for v, a in zip(values, answer)}
        for value, name in mapping.items():
            print(f"    {value!r} -> {name}")
        if input("  Confermi? [s/n]: ").strip().lower() == "s":
            return mapping


def standardize_column(column: pd.Series, attr: str) -> list:
    """Valori standardizzati della colonna scelta, uno per campione."""
    if attr in NUMERIC:
        numbers = pd.to_numeric(column, errors="coerce")
        failed = column[numbers.isna() & column.notna()].unique().tolist()
        if failed:
            print(f"  Attenzione, valori non numerici messi a mancante: {failed}")
        print(f"  Intervallo: {numbers.min()} - {numbers.max()}")
        return numbers.tolist()
    values = sorted(column.dropna().unique().tolist())
    if not values:
        return [None] * len(column)
    mapping = ask_values(values, STANDARD_VALUES[attr])
    return [mapping.get(v) for v in column]


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main():
    for gse in GSE_LIST:
        path = STD_DIR / gse / PICKLE_NAME
        with open(path, "rb") as f:
            data = pickle.load(f)
        if "std_meta" in data and not OVERWRITE:
            print(f"{gse}: std_meta già presente, salto")
            continue
        if "meta" not in data:
            print(f"{gse}: meta grezzo assente, rilancia standardize_cpg.py per questo dataset")
            continue
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
            std_meta[attr] = standardize_column(meta[col], attr) if col is not None else None

        data["std_meta"] = std_meta
        with open(path, "wb") as f:
            pickle.dump(data, f)
        print(f"\n{gse}: std_meta salvato in {path}")


if __name__ == "__main__":
    main()