#!/usr/bin/env python3
"""
harmonize_metadata.py — armonizzazione INTERATTIVA dei metadati GEO.

Per ogni dataset e per ogni attributo definito in ATTRIBUTES:
  1. cerca tra le colonne dei metadati quelle il cui nome corrisponde ai regex
     dell'attributo, e te le mostra con qualche valore di esempio;
  2. aspetta che tu scelga da terminale quale colonna usare (o nessuna);
  3. per gli attributi categorici propone una conversione dei valori
     (es. "F" -> 1, "M" -> 0) e ti chiede conferma o correzione.

Le scelte vengono salvate in CHOICES_FILE: rilanciando lo script non ti
richiede quello che hai già deciso (usa --redo per rifare).

Output:
    PROCESSED_DIR/<GSE>/meta_std.csv    -> ID, Dataset + attributi standard
    PROCESSED_DIR/meta_std_all.csv      -> tutti i dataset in una tabella

Uso:
    python harmonize_metadata.py                      # tutti i dataset
    python harmonize_metadata.py GSE42861 GSE40279    # solo alcuni
    python harmonize_metadata.py --redo GSE42861      # rifà le scelte di un dataset
"""
import json
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

try:
    from shared.paths import GSE_LIST, PROCESSED_DIR
except ImportError:
    GSE_LIST = [
        "GSE41037", "GSE40279", "GSE42861", "GSE27097",
        "GSE84727", "GSE111629", "GSE157131", "GSE197674",
    ]
    PROCESSED_DIR = Path("/work/gr-fe/bryan/EpiHydra/data/02_processed")

CHOICES_FILE = PROCESSED_DIR / "harmonization_choices.json"

# --------------------------------------------------------------------------
# ATTRIBUTI STANDARD
#
#   patterns : regex cercati nel NOME delle colonne (case-insensitive)
#   type     : "numeric"     -> convertito in numero (l'età in mesi/settimane/giorni
#                               diventa anni se l'unità è nel nome o nel valore)
#              "categorical" -> i valori vengono convertiti con "values"
#              "text"        -> copiato così com'è
#   values   : {valore_standard: [regex sui valori grezzi]}  (solo categorical)
#              i valori che non corrispondono a nessun regex te li chiede lui
# --------------------------------------------------------------------------
ATTRIBUTES = {
    "age": {
        "patterns": [r"^age$", r"^age_", r"\bage\b", r"age_at"],
        "type": "numeric",
    },
    "female": {
        "patterns": [r"sex", r"gender", r"female"],
        "type": "categorical",
        "values": {
            1: [r"^f$", r"^female$", r"^wom[ae]n$", r"^girl$"],
            0: [r"^m$", r"^male$", r"^m[ae]n$", r"^boy$"],
        },
    },
    "case": {
        "patterns": [r"disease", r"diagnos", r"status", r"condition", r"group",
                     r"subject", r"phenotype", r"case", r"control", r"patient"],
        "type": "categorical",
        "values": {
            0: [r"control", r"normal", r"healthy", r"^hc$", r"^ctrl", r"unaffected"],
            1: [r"patient", r"case", r"affected$", r"disease"],
        },
    },
    "condition": {
        "patterns": [r"disease", r"diagnos", r"condition", r"phenotype", r"subject"],
        "type": "text",
    },
    "smoking": {
        "patterns": [r"smok", r"tobacco", r"cigar"],
        "type": "categorical",
        "values": {
            "never": [r"never", r"^non", r"^no$", r"^0$"],
            "former": [r"former", r"^ex", r"past", r"quit"],
            "current": [r"current", r"^yes$", r"occasional", r"smoker$"],
        },
    },
    "ethnicity": {
        "patterns": [r"race", r"ethnic", r"ancestr", r"population"],
        "type": "text",
    },
    "tissue": {
        "patterns": [r"tissue", r"cell_type", r"source"],
        "type": "text",
    },
}

N_EXAMPLES = 5          # valori di esempio mostrati per colonna
MAX_CATEGORIES = 25     # oltre questo numero di valori distinti non chiede la mappa a mano


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def load_meta(gse: str) -> pd.DataFrame:
    """Legge meta_raw.csv; se manca, lo estrae dal pickle (lento) e lo crea."""
    csv_path = PROCESSED_DIR / gse / "meta_raw.csv"
    if csv_path.exists():
        return pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    print(f"  {gse}: meta_raw.csv non trovato, lo estraggo dal pickle (può volerci un po')...")
    with open(PROCESSED_DIR / gse / "methylation.pkl", "rb") as f:
        meta = pickle.load(f)["meta"]
    meta.to_csv(csv_path, index=False)
    return pd.read_csv(csv_path, dtype=str, keep_default_na=False)


def load_choices() -> dict:
    if CHOICES_FILE.exists():
        return json.loads(CHOICES_FILE.read_text())
    return {}


def save_choices(choices: dict) -> None:
    CHOICES_FILE.write_text(json.dumps(choices, indent=2, ensure_ascii=False))


# --------------------------------------------------------------------------
# INTERAZIONE
# --------------------------------------------------------------------------
def ask(prompt: str) -> str:
    try:
        return input(prompt).strip()
    except EOFError:
        sys.exit("\nInput terminato: le scelte fatte finora sono salvate.")


def describe_column(meta: pd.DataFrame, col: str) -> str:
    s = meta[col].astype(str).str.strip()
    filled = s[s != ""]
    counts = filled.value_counts()
    examples = ", ".join(f"{v!r}" for v in counts.index[:N_EXAMPLES])
    more = f" … (+{len(counts) - N_EXAMPLES})" if len(counts) > N_EXAMPLES else ""
    return (f"{len(filled)}/{len(s)} compilati, {len(counts)} valori distinti: "
            f"{examples}{more}")


def find_candidates(columns, patterns) -> list:
    regexes = [re.compile(p, re.IGNORECASE) for p in patterns]
    return [c for c in columns if any(r.search(c) for r in regexes)]


SKIP_COLUMNS = {"ID", "Organism", "platform"}


def choose_column(meta: pd.DataFrame, attr: str, patterns: list) -> str | None:
    """Mostra le colonne candidate e chiede quale usare. None = attributo assente."""
    columns = [c for c in meta.columns if c not in SKIP_COLUMNS]
    candidates = find_candidates(columns, patterns)
    show_all = not candidates

    while True:
        options = columns if show_all else candidates
        header = "tutte le colonne" if show_all else "colonne candidate"
        print(f"\n  [{attr}] {header}:")
        for i, col in enumerate(options, 1):
            print(f"    {i:>2}) {col:<30} {describe_column(meta, col)}")
        hint = "a = mostra tutte" if not show_all else "c = solo candidate"
        answer = ask(f"  Quale colonna è '{attr}'? [numero, 0 = nessuna, {hint}]: ").lower()

        if answer == "0":
            return None
        if answer == "a":
            show_all = True
            continue
        if answer == "c" and candidates:
            show_all = False
            continue
        if answer.isdigit() and 1 <= int(answer) <= len(options):
            return options[int(answer) - 1]
        print("  Risposta non valida, riprova.")


def parse_value(raw: str, kind: type):
    raw = raw.strip()
    if raw.lower() in ("", "nan", "na", "none"):
        return None
    if kind is int:
        try:
            return int(raw)
        except ValueError:
            pass
    return raw


def build_value_map(meta: pd.DataFrame, col: str, attr: str, spec: dict) -> dict | None:
    """Propone valore_grezzo -> valore_standard e chiede conferma/correzioni."""
    raw_values = sorted(v for v in meta[col].astype(str).str.strip().unique() if v != "")
    if len(raw_values) > MAX_CATEGORIES:
        print(f"  ATTENZIONE: '{col}' ha {len(raw_values)} valori distinti, "
              f"troppi per un attributo categorico. Lo salto.")
        return None

    allowed = list(spec["values"])
    kind = type(allowed[0])
    value_map = {}
    for v in raw_values:
        value_map[v] = None
        for std, regexes in spec["values"].items():
            if any(re.search(r, v, re.IGNORECASE) for r in regexes):
                value_map[v] = std
                break

    # 1) chiedi una volta i valori che i regex non hanno riconosciuto
    unassigned = [r for r, s in value_map.items() if s is None]
    if unassigned:
        print(f"\n  [{attr}] '{col}': {len(unassigned)} valori non riconosciuti. "
              f"Assegnali tra {allowed} (Invio = lascia mancante):")
        for raw in unassigned:
            value_map[raw] = parse_value(ask(f"    {raw!r} -> "), kind)

    # 2) mostra la tabella finale e chiedi conferma, finché non confermi
    while True:
        print(f"\n  [{attr}] conversione dei valori di '{col}':")
        for i, (raw, std) in enumerate(value_map.items(), 1):
            n = int((meta[col].astype(str).str.strip() == raw).sum())
            label = "(mancante)" if std is None else repr(std)
            print(f"    {i:>2}) {raw!r:<35} ({n:>4} campioni) -> {label}")
        ans = ask("  Confermi? [Invio = sì, numero = correggi quella riga]: ")
        if ans == "":
            return value_map
        if ans.isdigit() and 1 <= int(ans) <= len(value_map):
            raw = list(value_map)[int(ans) - 1]
            value_map[raw] = parse_value(ask(f"    {raw!r} -> "), kind)
        else:
            print("  Risposta non valida.")


# --------------------------------------------------------------------------
# CONVERSIONE
# --------------------------------------------------------------------------
def to_years(col_name: str, raw: str) -> float:
    m = re.search(r"\d+(?:\.\d+)?", str(raw))
    if not m:
        return np.nan
    val = float(m.group())
    unit = f"{col_name} {raw}".lower()
    if re.search(r"month", unit):
        val /= 12
    elif re.search(r"week", unit):
        val /= 52
    elif re.search(r"\bday", unit):
        val /= 365
    return round(val, 3)


def convert(meta: pd.DataFrame, attr: str, choice: dict) -> pd.Series:
    col = choice.get("column")
    if col is None or col not in meta.columns:
        return pd.Series(np.nan, index=meta.index)
    raw = meta[col].astype(str).str.strip()
    kind = ATTRIBUTES[attr]["type"]

    if kind == "numeric":
        return raw.map(lambda v: to_years(col, v) if attr == "age" else
                       pd.to_numeric(v, errors="coerce"))
    if kind == "categorical":
        vmap = choice.get("value_map") or {}
        return raw.map(lambda v: vmap.get(v, np.nan) if vmap.get(v) is not None else np.nan)
    return raw.replace("", np.nan)


# --------------------------------------------------------------------------
# PIPELINE
# --------------------------------------------------------------------------
def harmonize_dataset(gse: str, choices: dict, redo: bool) -> pd.DataFrame:
    meta = load_meta(gse)
    print(f"\n{'=' * 70}\n{gse}: {len(meta)} campioni, {len(meta.columns)} colonne\n{'=' * 70}")

    ds_choices = {} if redo else choices.get(gse, {})
    for attr, spec in ATTRIBUTES.items():
        if attr in ds_choices:
            print(f"  [{attr}] già deciso: colonna = {ds_choices[attr].get('column')}")
            continue
        col = choose_column(meta, attr, spec["patterns"])
        choice = {"column": col}
        if col is not None and spec["type"] == "categorical":
            choice["value_map"] = build_value_map(meta, col, attr, spec)
        ds_choices[attr] = choice
        choices[gse] = ds_choices
        save_choices(choices)          # salva subito: se interrompi non perdi nulla

    std = pd.DataFrame({"ID": meta["ID"], "Dataset": gse})
    for attr in ATTRIBUTES:
        std[attr] = convert(meta, attr, ds_choices[attr]).values

    out = PROCESSED_DIR / gse / "meta_std.csv"
    std.to_csv(out, index=False)
    filled = {a: f"{std[a].notna().mean() * 100:.0f}%" for a in ATTRIBUTES}
    print(f"\n  {gse}: salvato {out}\n  compilazione: {filled}")
    return std


def main(argv: list) -> None:
    redo = "--redo" in argv
    selected = [a for a in argv if a.startswith("GSE")] or GSE_LIST
    choices = load_choices()
    if redo:
        for gse in selected:
            choices.pop(gse, None)

    tables = []
    for gse in selected:
        try:
            tables.append(harmonize_dataset(gse, choices, redo=False))
        except FileNotFoundError as e:
            print(f"\n{gse}: file non trovato ({e}), salto.")

    # Tabella unica con tutti i dataset già armonizzati (anche da run precedenti)
    all_tables = []
    for gse in GSE_LIST:
        path = PROCESSED_DIR / gse / "meta_std.csv"
        if path.exists():
            all_tables.append(pd.read_csv(path))
    if all_tables:
        combined = pd.concat(all_tables, ignore_index=True)
        out = PROCESSED_DIR / "meta_std_all.csv"
        combined.to_csv(out, index=False)
        print(f"\nTabella unica: {combined.shape[0]} campioni da "
              f"{combined['Dataset'].nunique()} dataset -> {out}")
    print(f"Scelte salvate in {CHOICES_FILE}")


if __name__ == "__main__":
    main(sys.argv[1:])