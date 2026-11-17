#!/usr/bin/env python3
"""
geo_data_download.py

Scarica i Series Matrix GEO per una lista di GSE, legge la matrice di
metilazione (dal Series Matrix o da un file supplementare) e TUTTI i metadati
dei campioni, e salva un pickle per dataset:

    PROCESSED_DIR/<GSE>/methylation.pkl  ->  {"meta": DataFrame, "CpG": DataFrame}

"meta" contiene:
  - ID, title, source, description, Organism, platform
  - una colonna per ogni caratteristica GEO, con nome completo normalizzato
    (es. "disease state" -> "disease_state", "smoking status" -> "smoking_status")
  - "age": età numerica in anni, ricavata SOLO da campi che si chiamano davvero età

L'armonizzazione (female, case, ...) si fa dopo, con harmonize_metadata.py.
"""
import re
import gzip
import pickle
import subprocess
import sys
import traceback
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin, unquote

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

try:
    from shared.paths import GSE_LIST, PROCESSED_DIR, RAW_DIR
except ImportError:  # fallback standalone
    GSE_LIST = [
        "GSE41037", "GSE40279", "GSE42861", "GSE27097",
        "GSE84727", "GSE111629", "GSE157131", "GSE197674",
    ]
    RAW_DIR = Path("/work/gr-fe/bryan/EpiHydra/data/01_raw")
    PROCESSED_DIR = Path("/work/gr-fe/bryan/EpiHydra/data/02_processed")

# Serie con più piattaforme: quale file del Series Matrix usare
PLATFORM_OVERRIDE = {
    "GSE157131": "GPL21145",
}

# Sotto questa dimensione si assume che il Series Matrix non contenga i beta
SERIES_MATRIX_MIN_BYTES = 10 * 1024 * 1024

# Nome di campo che indica l'età: "age", "age (y)", "age at visit", "age_years", ...
AGE_FIELD_RE = re.compile(r"^age(\b|_|\s|\(|$)", re.IGNORECASE)


# --------------------------------------------------------------------------
# DOWNLOAD
# --------------------------------------------------------------------------
def _gse_prefix(gse: str) -> str:
    """GSE40279 -> GSE40nnn (cartella FTP di GEO)."""
    return gse[:-3] + "nnn"


def series_matrix_url(gse: str) -> str:
    name = f"{gse}_series_matrix.txt.gz"
    if gse in PLATFORM_OVERRIDE:
        name = f"{gse}-{PLATFORM_OVERRIDE[gse]}_series_matrix.txt.gz"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_gse_prefix(gse)}/{gse}/matrix/{name}"


def supplementary_dir_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_gse_prefix(gse)}/{gse}/suppl/"


def curl_download(url: str, filepath: Path) -> None:
    """Scarica con curl in un file .part, valida, poi rinomina."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = filepath.with_name(filepath.name + ".part")
    if tmp_path.exists():
        tmp_path.unlink()

    cmd = [
        "curl", "-fL", "--retry", "5", "--retry-delay", "5", "--retry-connrefused",
        "--connect-timeout", "60", "--max-time", "7200", "--progress-bar",
        "-o", str(tmp_path), url,
    ]
    try:
        subprocess.run(cmd, check=True)
        if not tmp_path.exists() or tmp_path.stat().st_size == 0:
            raise ValueError("Downloaded file is empty.")
        if filepath.name.lower().endswith(".gz"):
            with gzip.open(tmp_path, "rb") as f:
                f.read(1024)
        tmp_path.rename(filepath)
        print(f"Successfully downloaded: {filepath}")
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def download_series_matrix(gse: str, filepath: Path) -> None:
    if filepath.exists() and filepath.stat().st_size > 0:
        print(f"{gse}: series matrix already exists, skipping download.")
        return
    print(f"{gse}: downloading series matrix...")
    curl_download(series_matrix_url(gse), filepath)


def find_supplementary_beta(gse: str) -> tuple[str, str]:
    """Cerca in suppl/ un file con 'beta' o 'normalised/normalized' nel nome.

    Restituisce (filename, url).
    """

    class HrefParser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.hrefs = []

        def handle_starttag(self, tag, attrs):
            if tag.lower() == "a":
                for key, value in attrs:
                    if key.lower() == "href" and value:
                        self.hrefs.append(value)

    supp_url = supplementary_dir_url(gse)
    print(f"{gse}: searching supplementary files at {supp_url}")
    try:
        result = subprocess.run(
            ["curl", "-fLsS", "--connect-timeout", "60", "--max-time", "300", supp_url],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Could not retrieve supplementary directory for {gse}. "
            f"curl exit code: {e.returncode}. Error: {e.stderr}"
        )

    parser = HrefParser()
    parser.feed(result.stdout)

    include_terms = ("beta", "normalised", "normalized")
    exclude_terms = ("detection", "pvalue", "p-value", "intensity", "methylated",
                     "unmethylated", "annotation", "manifest", "control", "idat")
    valid_ext = (".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz")

    candidates = []
    for href in parser.hrefs:
        filename = unquote(href).split("/")[-1]
        lower = filename.lower()
        if not lower.endswith(valid_ext):
            continue
        if any(t in lower for t in exclude_terms) or not any(t in lower for t in include_terms):
            continue
        normalised = "normalised" in lower or "normalized" in lower
        if "beta" in lower and normalised:
            score = 100
        elif normalised:
            score = 80
        else:
            score = 60
        candidates.append({"score": score, "filename": filename, "url": urljoin(supp_url, href)})

    if not candidates:
        raise RuntimeError(f"No beta or normalised/normalized candidate files found for {gse}.")

    candidates.sort(key=lambda x: x["score"], reverse=True)
    print(f"{gse}: supplementary candidates found:")
    for c in candidates:
        print(f"  Score {c['score']}: {c['filename']}")

    best = [c for c in candidates if c["score"] == candidates[0]["score"]]
    if len(best) > 1:
        raise RuntimeError(
            f"Multiple equally ranked supplementary files found for {gse}: "
            f"{[x['filename'] for x in best]}. Please select the correct file manually."
        )
    print(f"{gse}: selected supplementary file: {best[0]['filename']}")
    return best[0]["filename"], best[0]["url"]


def download_supplementary_beta(gse: str, directory: Path) -> Path:
    filename, url = find_supplementary_beta(gse)
    filepath = directory / filename
    if filepath.exists() and filepath.stat().st_size > 0:
        print(f"{gse}: supplementary beta file already exists, skipping.")
        return filepath
    print(f"{gse}: downloading supplementary beta values...")
    curl_download(url, filepath)
    return filepath


# --------------------------------------------------------------------------
# METADATI
# --------------------------------------------------------------------------
def _clean(s: str) -> str:
    return s.strip().strip('"').strip()


def normalise_field_name(name: str) -> str:
    """'Disease State' -> 'disease_state', 'age (y)' -> 'age_y'."""
    name = name.strip().lower()
    name = re.sub(r"[^\w]+", "_", name)
    return name.strip("_") or "unnamed"


def parse_age_years(field: str, value: str) -> float:
    """Estrae l'età in anni da 'age (months)', '45', '45.5 years', ..."""
    match = re.search(r"\d+(?:\.\d+)?", value)
    if not match:
        return np.nan
    years = float(match.group())
    unit_text = f"{field} {value}".lower()
    if re.search(r"\bmonth", unit_text):
        years /= 12
    elif re.search(r"\bweek", unit_text):
        years /= 52
    elif re.search(r"\bday", unit_text):
        years /= 365
    return round(years, 3)


def parse_series_matrix_meta(filepath: Path) -> pd.DataFrame:
    """Legge l'header del Series Matrix e restituisce una riga per campione.

    Le righe Sample_characteristics_ch1 sono lette CAMPIONE PER CAMPIONE:
    ogni cella "campo: valore" va nella colonna del proprio campo, quindi
    righe che mescolano campi diversi tra campioni non creano colonne sbagliate.
    """
    sample_fields = {}           # campo GEO standard -> lista di valori
    characteristics = []         # una lista di celle per ogni riga di caratteristiche

    with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            values = [_clean(v) for v in parts[1:]]
            if key.startswith("Sample_characteristics"):
                characteristics.append(values)
            elif key not in sample_fields:   # tiene la prima occorrenza
                sample_fields[key] = values

    ids = sample_fields["Sample_geo_accession"]
    n = len(ids)

    rows = [{} for _ in range(n)]
    age_fields_seen = set()
    for values in characteristics:
        for i, cell in enumerate(values[:n]):
            if not cell or ":" not in cell:
                continue
            field, value = cell.split(":", 1)
            col = normalise_field_name(field)
            value = value.strip()
            # Se lo stesso campo compare due volte per un campione, non sovrascrivere
            if col in rows[i]:
                col = f"{col}_2"
            rows[i][col] = value
            if AGE_FIELD_RE.match(field.strip()):
                age_fields_seen.add(col)
                if "age" not in rows[i] or pd.isna(rows[i].get("age")):
                    rows[i]["age"] = parse_age_years(field, value)

    chars_df = pd.DataFrame(rows)
    if "age" in chars_df:
        chars_df["age"] = pd.to_numeric(chars_df["age"], errors="coerce").astype(np.float32)
        print(f"  age ricavata dai campi: {sorted(age_fields_seen)}")
    else:
        print("  note: nessun campo età trovato.")

    base = pd.DataFrame({
        "ID": ids,
        "title": sample_fields.get("Sample_title", [""] * n),
        "source": sample_fields.get("Sample_source_name_ch1", [""] * n),
        "description": sample_fields.get("Sample_description", [""] * n),
        "Organism": sample_fields.get("Sample_organism_ch1", [""] * n),
        "platform": sample_fields.get("Sample_platform_id", [""] * n),
    })
    meta = pd.concat([base, chars_df], axis=1)
    print(f"  metadati: {meta.shape[0]} campioni, colonne: {list(meta.columns)}")
    return meta


# --------------------------------------------------------------------------
# MATRICE DI METILAZIONE
# --------------------------------------------------------------------------
def parse_expression_matrix(filepath: Path) -> pd.DataFrame:
    """Matrice CpG x campioni dal Series Matrix."""
    return pd.read_csv(filepath, sep="\t", comment="!", index_col=0, compression="gzip")


def _map_columns_to_gsm(columns, meta: pd.DataFrame) -> dict:
    """Mappa le colonne di un file supplementare agli ID GSM.

    Prova prima la corrispondenza diretta con il GSM, poi con title e description
    (spesso contengono l'ID del chip, es. 9285451020_R01C01, o il nome del campione).
    """
    cols = [str(c) for c in columns]
    mapping = {c: c for c in cols if c in set(meta["ID"])}
    if mapping:
        return mapping

    for field in ("title", "description"):
        lookup = {}
        for gsm, label in zip(meta["ID"], meta[field].astype(str)):
            if label:
                lookup[label] = gsm
        found = {c: lookup[c] for c in cols if c in lookup}
        if len(found) > len(mapping):
            mapping = found

    if not mapping:
        # ultima prova: la colonna contiene il title (es. "Sample_12.AVG_Beta")
        for gsm, label in zip(meta["ID"], meta["title"].astype(str)):
            hits = [c for c in cols if label and label in c and "pval" not in c.lower()]
            if len(hits) == 1:
                mapping[hits[0]] = gsm
    return mapping


def parse_supplementary_beta(filepath: Path, meta: pd.DataFrame) -> pd.DataFrame:
    print(f"Parsing supplementary beta file: {filepath.name}")
    opener = gzip.open if filepath.name.lower().endswith(".gz") else open
    with opener(filepath, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()
    delimiter = "," if header.count(",") > header.count("\t") else "\t"

    df = pd.read_csv(filepath, sep=delimiter, index_col=0, compression="infer", low_memory=False)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    # Campioni sulle colonne (caso tipico) o sulle righe?
    col_map = _map_columns_to_gsm(df.columns, meta)
    row_map = _map_columns_to_gsm(df.index, meta)
    if len(row_map) > len(col_map):
        df, col_map = df.T, row_map

    if not col_map:
        raise ValueError(
            "Could not match supplementary columns to GEO samples (GSM, title, description). "
            f"First columns: {list(df.columns[:5])}; first titles: {list(meta['title'][:5])}"
        )

    expr = df[list(col_map)].rename(columns=col_map)
    missing = [s for s in meta["ID"] if s not in expr.columns]
    if missing:
        print(f"  warning: {len(missing)} campioni GEO non trovati nel file "
              f"supplementare (es. {missing[:5]}); verranno esclusi.")
    expr = expr.loc[:, [s for s in meta["ID"] if s in expr.columns]]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    print(f"Loaded supplementary methylation matrix: {expr.shape[0]} CpGs x {expr.shape[1]} samples")
    return expr


def check_beta_range(expr: pd.DataFrame, gse: str) -> None:
    """Avvisa se i valori non sembrano beta (0-1), es. M-values o percentuali."""
    vmin, vmax = np.nanmin(expr.values), np.nanmax(expr.values)
    print(f"  range valori: {vmin:.3f} - {vmax:.3f}")
    if vmin < -0.01 or vmax > 1.01:
        kind = "percentuali (0-100)?" if vmin >= 0 and vmax <= 100 else "M-values?"
        print(f"  WARNING {gse}: valori fuori da [0, 1] -> probabilmente {kind} "
              "Da convertire prima di usare clock o MethylGPT.")


def impute_missing_median(expr: pd.DataFrame) -> pd.DataFrame:
    """Imputa i NaN con la mediana della stessa CpG su tutti i campioni."""
    n_missing_before = int(expr.isna().sum().sum())
    if n_missing_before == 0:
        print("  no missing values found, skipping imputation.")
        return expr
    row_medians = expr.median(axis=1, skipna=True)
    expr_imputed = expr.T.fillna(row_medians).T
    n_missing_after = int(expr_imputed.isna().sum().sum())
    n_unfillable = int(expr_imputed.isna().any(axis=1).sum())
    print(f"  median imputation: {n_missing_before} missing values found, "
          f"{n_missing_before - n_missing_after} imputed.")
    if n_unfillable:
        print(f"  warning: {n_unfillable} CpG site(s) missing for ALL samples remain NaN.")
    return expr_imputed


# --------------------------------------------------------------------------
# PIPELINE
# --------------------------------------------------------------------------
def process_gse(gse: str) -> None:
    print(f"\n=== Processing {gse} ===")
    raw_dir = RAW_DIR / gse
    raw_filepath = raw_dir / f"{gse}_series_matrix.txt.gz"
    processed_filepath = PROCESSED_DIR / gse / "methylation.pkl"

    if processed_filepath.exists():
        print(f"{gse}: processed pickle already exists at {processed_filepath}, skipping.")
        return

    download_series_matrix(gse, raw_filepath)
    meta = parse_series_matrix_meta(raw_filepath)

    size = raw_filepath.stat().st_size
    print(f"{gse}: series matrix size = {size / 1024**2:.2f} MB")
    if size < SERIES_MATRIX_MIN_BYTES:
        print(f"{gse}: series matrix below threshold, searching supplementary beta values.")
        beta_path = download_supplementary_beta(gse, raw_dir)
        expr = parse_supplementary_beta(beta_path, meta)
    else:
        print(f"{gse}: using series matrix methylation data.")
        expr = parse_expression_matrix(raw_filepath)
        expr = expr.loc[:, [s for s in meta["ID"] if s in expr.columns]]

    # Tieni nei metadati solo i campioni che hanno dati di metilazione
    meta = meta[meta["ID"].isin(expr.columns)].reset_index(drop=True)

    check_beta_range(expr, gse)
    expr = impute_missing_median(expr)

    processed_filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(processed_filepath, "wb") as f:
        pickle.dump({"meta": meta, "CpG": expr}, f)
    # Copia leggera dei metadati: l'armonizzazione la legge senza caricare la matrice CpG
    meta.to_csv(processed_filepath.parent / "meta_raw.csv", index=False)
    print(f"{gse}: saved processed data to {processed_filepath}")


def main():
    succeeded, failed = [], []
    for gse in GSE_LIST:
        try:
            process_gse(gse)
            succeeded.append(gse)
        except Exception as e:
            print(f"\n{gse}: FAILED — {type(e).__name__}: {e}")
            traceback.print_exc()
            failed.append(gse)

    print("\n=== Summary ===")
    print(f"Succeeded ({len(succeeded)}): {succeeded}")
    if failed:
        print(f"Failed ({len(failed)}): {failed}")
        print("Re-run the script to retry only the failed datasets.")


if __name__ == "__main__":
    main()