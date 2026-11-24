#!/usr/bin/env python3
"""
geo_data_download_raw.py

Downloads GEO methylation matrices to RAW_DIR/<GSE>/<GSE>_series_matrix.txt.gz,
ready for standardize_cpg.py. No imputation, no filtering, no pickles.

For the datasets in SUPPLEMENTARY the GEO Series Matrix only holds metadata: the
matrix is downloaded from the supplementary file and appended to the header, with
columns renamed to GSM IDs. Values are not modified. The original header is kept in
<GSE>_header.txt.gz.

Finished files are skipped and interrupted downloads resume: if a dataset fails,
just run the script again.
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
# CONSTANTS
# --------------------------------------------------------------------------
DATA_DIR = Path("/work/scitas-share/FellayMultiOmic/CS433/ML4Science_Barry/data")
RAW_DIR = DATA_DIR / "01_raw"                    # same folder read by standardize_cpg.py
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"

GSE_LIST = [
    "GSE41037", "GSE40279", "GSE42861", "GSE27097",
    "GSE84727", "GSE111629", "GSE157131", "GSE197674",
]

# multi-platform series: which Series Matrix to use
PLATFORM = {"GSE157131": "GPL21145"}

# series whose matrix is only in a supplementary file
SUPPLEMENTARY = {
    "GSE84727":  "GSE84727_normalisedBetas.csv.gz",
    "GSE111629": "GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz",
    "GSE157131": "GSE157131_Matrix_processed_beta_geo_08252020.txt.gz",
    "GSE197674": "GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_processed.txt.gz",
}

MIN_MATRIX_BYTES = 10 * 1024 * 1024   # below: the Series Matrix only holds metadata
MIN_MATCHED = 0.9                     # minimum fraction of samples to match
ATTEMPTS = 4                          # attempts per file
WAIT_SECONDS = 60                     # wait between attempts (60, 120, 180)

CHIP_ID = re.compile(r"\d+_R\d\dC\d\d")                              # e.g. 9285451020_R01C01
BETA_SUFFIX = re.compile(r"[._ -]?(avg[._ -]?)?beta$", re.IGNORECASE)   # e.g. S1.AVG_Beta


# --------------------------------------------------------------------------
# DOWNLOAD
# --------------------------------------------------------------------------
def download(url: str, path: Path) -> None:
    """Download with curl to a resumable .download file, then rename."""
    if path.exists() and path.stat().st_size > 0:
        print(f"  {path.name}: already present")
        return
    tmp = path.with_name(path.name + ".download")
    # no total time limit: stops only if stalled for 5 minutes
    cmd = ["curl", "-fL", "-sS", "-C", "-", "--connect-timeout", "60",
           "--speed-limit", "1000", "--speed-time", "300", "-o", str(tmp), url]
    for attempt in range(1, ATTEMPTS + 1):
        print(f"  {path.name}: download (attempt {attempt}/{ATTEMPTS})", flush=True)
        code = subprocess.run(cmd).returncode
        if code == 0:
            break
        if attempt == ATTEMPTS:
            raise RuntimeError(f"download failed (curl {code}): {url}")
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
# SUPPLEMENTARY FILE -> SERIES MATRIX
# --------------------------------------------------------------------------
def read_header(path: Path):
    """Header lines, sample GSMs and candidate labels (one list per line)."""
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
            # chip ID, also inside idat file names
            chips = [CHIP_ID.search(v) for v in values]
            labels.append([m.group() if m else "" for m in chips])
    if not gsm:
        raise ValueError(f"no GSM in the header of {path}")
    return lines, gsm, labels


def map_columns(columns: list, gsm: list, labels: list) -> dict:
    """{column position: GSM}, using the metadata line with the most matches."""
    best = {}
    for values in labels:
        lookup, repeated = {}, {""}
        for label, sample in zip(values, gsm):
            if label in lookup:
                repeated.add(label)
            lookup[label] = sample
        for label in repeated:                 # unique labels only
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
    """Write header + supplementary matrix (GSM columns), one row at a time."""
    lines, gsm, labels = read_header(header)
    tmp = out.with_name(out.name + ".build")
    with gzip.open(supplementary, "rt", encoding="utf-8", errors="replace") as f:
        first, second = f.readline().rstrip("\r\n"), f.readline().rstrip("\r\n")
        sep = max("\t,; ", key=first.count)

        def split(line):
            return [c.strip('"') for c in line.rstrip("\r\n").split(sep)]

        columns, row = split(first), split(second)
        if len(columns) == len(row) - 1:       # header without a name for the CpG column (R)
            columns = [""] + columns
        if len(columns) != len(row):
            raise ValueError(f"header has {len(columns)} columns, first row has {len(row)}")
        mapping = map_columns(columns, gsm, labels)
        print(f"  matched samples: {len(mapping)}/{len(gsm)} (columns in file: {len(columns) - 1})")
        if len(mapping) < max(2, MIN_MATCHED * len(gsm)):
            raise ValueError(f"too many unmatched samples. Columns: {columns[:6]}; "
                             f"first row: {row[:6]}; metadata of the first sample: "
                             f"{[v[0] for v in labels if v and v[0]][:60]}")
        positions = list(mapping)
        ids = [i for i, v in enumerate(row) if i not in mapping and re.match(r"(cg|ch\.|rs)", v)]
        if not ids:
            raise ValueError(f"CpG column not found; first row: {row[:6]}")
        pick = itemgetter(ids[0], *positions)
        for v in pick(row)[1:]:
            if v not in ("", "NA", "NaN", "null"):
                float(v)                       # error if the matrix is not numeric

        n = 0
        with gzip.open(tmp, "wt", encoding="utf-8", compresslevel=1, newline="\n") as g:
            g.writelines(lines)
            g.write("!series_matrix_table_begin\n")
            g.write("\t".join(["ID_REF"] + [mapping[i] for i in positions]) + "\n")
            while row:
                if len(row) != len(columns):
                    raise ValueError(f"row {n + 1}: {len(row)} columns instead of {len(columns)}")
                g.write("\t".join(pick(row)) + "\n")
                n += 1
                line = f.readline()
                row = split(line) if line.strip() else None
            g.write("!series_matrix_table_end\n")
    tmp.replace(out)
    print(f"  {out.name}: {n} CpGs x {len(positions)} samples")


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def process(gse: str) -> None:
    folder = RAW_DIR / gse
    folder.mkdir(parents=True, exist_ok=True)
    matrix = folder / f"{gse}_series_matrix.txt.gz"
    if matrix.exists() and matrix.stat().st_size >= MIN_MATRIX_BYTES:
        print("  matrix already done")
        return
    if not os.access(folder, os.W_OK):
        raise PermissionError(f"no write permission in {folder}")

    if gse not in SUPPLEMENTARY:
        download(series_matrix_url(gse), matrix)
        if matrix.stat().st_size < MIN_MATRIX_BYTES:
            raise ValueError(f"{matrix.name} only holds metadata: add {gse} to SUPPLEMENTARY")
        return

    header = folder / f"{gse}_header.txt.gz"
    if matrix.exists() and not header.exists():      # metadata-only Series Matrix already downloaded
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
            print(f"{gse}: FAILED - {type(e).__name__}: {e}")
            traceback.print_exc()
            failed.append(gse)
    print(f"\nReady: {len(GSE_LIST) - len(failed)}/{len(GSE_LIST)}. Failed: {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
