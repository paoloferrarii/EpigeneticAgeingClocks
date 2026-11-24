#!/usr/bin/env python3
"""
geo_data_download_raw.py

Scarica le matrici di metilazione GEO in RAW_DIR/<GSE>/<GSE>_series_matrix.txt.gz,
pronte per standardize_cpg.py. Nessuna imputazione, nessun filtro, nessun pickle.

Per i dataset in SUPPLEMENTARY il Series Matrix di GEO contiene solo i metadati: la
matrice viene scaricata dal file supplementare e accodata all'header, con le colonne
rinominate in GSM. I valori non vengono modificati. L'header originale resta in
<GSE>_header.txt.gz.

I file già pronti vengono saltati e i download interrotti riprendono: se un dataset
fallisce basta rilanciare lo script.
"""
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from operator import itemgetter
from pathlib import Path

# --------------------------------------------------------------------------
# COSTANTI
# --------------------------------------------------------------------------
DATA_DIR = Path("/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/data")
RAW_DIR = DATA_DIR / "01_raw"                    # stessa cartella letta da standardize_cpg.py
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"

GSE_LIST = [
    "GSE41037", "GSE40279", "GSE42861", "GSE27097",
    "GSE84727", "GSE111629", "GSE157131", "GSE197674",
]

# serie con più piattaforme: quale Series Matrix usare
PLATFORM = {"GSE157131": "GPL21145"}

# serie con la matrice solo in un file supplementare
SUPPLEMENTARY = {
    "GSE84727":  "GSE84727_normalisedBetas.csv.gz",
    "GSE111629": "GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz",
    "GSE157131": "GSE157131_Matrix_processed_beta_geo_08252020.txt.gz",
    "GSE197674": "GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_processed.txt.gz",
}

MIN_MATRIX_BYTES = 10 * 1024 * 1024   # sotto: il Series Matrix ha solo metadati
MIN_MATCHED = 0.9                     # frazione minima di campioni da abbinare
ATTEMPTS = 4                          # tentativi per file
WAIT_SECONDS = 60                     # attesa tra i tentativi (60, 120, 180)

CHIP_ID = re.compile(r"\d+_R\d\dC\d\d")                              # es. 9285451020_R01C01
BETA_SUFFIX = re.compile(r"[._ -]?(avg[._ -]?)?beta$", re.IGNORECASE)   # es. S1.AVG_Beta


# --------------------------------------------------------------------------
# DOWNLOAD
# --------------------------------------------------------------------------
def download(url: str, path: Path) -> None:
    """Scarica con curl in un file .download (riprendibile), poi rinomina."""
    if path.exists() and path.stat().st_size > 0:
        print(f"  {path.name}: già presente")
        return
    tmp = path.with_name(path.name + ".download")
    # niente limite di tempo totale: si ferma solo se è bloccato per 5 minuti
    cmd = ["curl", "-fL", "-sS", "-C", "-", "--connect-timeout", "60",
           "--speed-limit", "1000", "--speed-time", "300", "-o", str(tmp), url]
    for attempt in range(1, ATTEMPTS + 1):
        print(f"  {path.name}: download (tentativo {attempt}/{ATTEMPTS})", flush=True)
        code = subprocess.run(cmd).returncode
        if code == 0:
            break
        if attempt == ATTEMPTS:
            raise RuntimeError(f"download fallito (curl {code}): {url}")
        time.sleep(WAIT_SECONDS * attempt)
    try:
        with gzip.open(tmp, "rb") as f:
            f.read(1024)
    except OSError:
        tmp.unlink()
        raise
    tmp.replace(path)


def series_url(gse: str) -> str:
    return f"{GEO_URL}/{gse[:-3]}nnn/{gse}"


def series_matrix_url(gse: str) -> str:
    name = f"{gse}-{PLATFORM[gse]}" if gse in PLATFORM else gse
    return f"{series_url(gse)}/matrix/{name}_series_matrix.txt.gz"


# --------------------------------------------------------------------------
# FILE SUPPLEMENTARE -> SERIES MATRIX
# --------------------------------------------------------------------------
def read_header(path: Path):
    """Righe di header, GSM dei campioni e possibili etichette (una lista per riga)."""
    lines, gsm, labels = [], None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            lines.append(line)
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            values = [v.strip().strip('"').strip() for v in values]
            if key == "!Sample_geo_accession":
                gsm = values
            labels.append(values)
            labels.append([v.split(":", 1)[-1].strip() for v in values])
            # ID del chip, anche dentro i nomi dei file idat
            chips = [CHIP_ID.search(v) for v in values]
            labels.append([m.group() if m else "" for m in chips])
    if not gsm:
        raise ValueError(f"nessun GSM nell'header di {path}")
    return lines, gsm, labels


def map_columns(columns: list, gsm: list, labels: list) -> dict:
    """{posizione della colonna: GSM}, usando la riga di metadati con più corrispondenze."""
    best = {}
    for values in labels:
        lookup, repeated = {}, {""}
        for label, sample in zip(values, gsm):
            if label in lookup:
                repeated.add(label)
            lookup[label] = sample
        for label in repeated:                 # solo etichette univoche
            lookup.pop(label, None)
        found = {}
        for i, col in enumerate(columns):
            base = BETA_SUFFIX.sub("", col)
            for name in (col, base, re.sub(r"^X", "", col), re.sub(r"^X", "", base)):
                if name in lookup:
                    found[i] = lookup[name]
                    break
        if len(set(found.values())) == len(found) > len(best):
            best = found
    return best


def build_series_matrix(header: Path, supplementary: Path, out: Path) -> None:
    """Scrive header + matrice supplementare (colonne GSM), una riga alla volta."""
    lines, gsm, labels = read_header(header)
    tmp = out.with_name(out.name + ".build")
    with gzip.open(supplementary, "rt", encoding="utf-8", errors="replace") as f:
        first, second = f.readline().rstrip("\r\n"), f.readline().rstrip("\r\n")
        sep = max("\t,; ", key=first.count)

        def split(line):
            return [c.strip('"') for c in line.rstrip("\r\n").split(sep)]

        columns, row = split(first), split(second)
        if len(columns) == len(row) - 1:       # header senza nome per la colonna CpG (R)
            columns = [""] + columns
        if len(columns) != len(row):
            raise ValueError(f"header con {len(columns)} colonne, prima riga con {len(row)}")
        mapping = map_columns(columns, gsm, labels)
        print(f"  campioni abbinati: {len(mapping)}/{len(gsm)} (colonne nel file: {len(columns) - 1})")
        if len(mapping) < max(2, MIN_MATCHED * len(gsm)):
            raise ValueError(f"troppi campioni non abbinati. Colonne: {columns[:6]}; "
                             f"prima riga: {row[:6]}; metadati del primo campione: "
                             f"{[v[0] for v in labels if v and v[0]][:60]}")
        positions = list(mapping)
        ids = [i for i, v in enumerate(row) if i not in mapping and re.match(r"(cg|ch\.|rs)", v)]
        if not ids:
            raise ValueError(f"colonna delle CpG non trovata; prima riga: {row[:6]}")
        pick = itemgetter(ids[0], *positions)
        for v in pick(row)[1:]:
            if v not in ("", "NA", "NaN", "null"):
                float(v)                       # errore se la matrice non è numerica

        n = 0
        with gzip.open(tmp, "wt", encoding="utf-8", compresslevel=1, newline="\n") as g:
            g.writelines(lines)
            g.write("!series_matrix_table_begin\n")
            g.write("\t".join(["ID_REF"] + [mapping[i] for i in positions]) + "\n")
            while row:
                if len(row) != len(columns):
                    raise ValueError(f"riga {n + 1}: {len(row)} colonne invece di {len(columns)}")
                g.write("\t".join(pick(row)) + "\n")
                n += 1
                line = f.readline()
                row = split(line) if line.strip() else None
            g.write("!series_matrix_table_end\n")
    tmp.replace(out)
    print(f"  {out.name}: {n} CpG x {len(positions)} campioni")


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def process(gse: str) -> None:
    folder = RAW_DIR / gse
    folder.mkdir(parents=True, exist_ok=True)
    matrix = folder / f"{gse}_series_matrix.txt.gz"
    if matrix.exists() and matrix.stat().st_size >= MIN_MATRIX_BYTES:
        print("  matrice già pronta")
        return
    if not os.access(folder, os.W_OK):
        raise PermissionError(f"non puoi scrivere in {folder}")

    if gse not in SUPPLEMENTARY:
        download(series_matrix_url(gse), matrix)
        if matrix.stat().st_size < MIN_MATRIX_BYTES:
            raise ValueError(f"{matrix.name} contiene solo metadati: aggiungi {gse} a SUPPLEMENTARY")
        return

    header = folder / f"{gse}_header.txt.gz"
    if matrix.exists() and not header.exists():      # Series Matrix di soli metadati già scaricato
        shutil.copyfile(matrix, header)
    download(series_matrix_url(gse), header)
    supplementary = folder / SUPPLEMENTARY[gse]
    download(f"{series_url(gse)}/suppl/{supplementary.name}", supplementary)
    build_series_matrix(header, supplementary, matrix)


def main():
    failed = []
    for gse in GSE_LIST:
        print(f"\n=== {gse} ===", flush=True)
        try:
            process(gse)
        except Exception as e:
            print(f"{gse}: FALLITO - {type(e).__name__}: {e}")
            traceback.print_exc()
            failed.append(gse)
    print(f"\nPronti: {len(GSE_LIST) - len(failed)}/{len(GSE_LIST)}. Falliti: {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
