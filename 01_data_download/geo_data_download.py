#!/usr/bin/env python3
"""
geo_data_download.py

Downloads GEO Series Matrix files for a list of GSE datasets, parses the
expression/methylation matrix plus sample metadata, and pickles the result.

Generalized from a notebook originally written for GSE40279.

Usage:
    python geo_data_download.py

Just edit the GSE_LIST below (and RAW_DIR / PROCESSED_DIR if needed) and run.
"""
import re
import gzip
import pickle
import traceback
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from html.parser import HTMLParser
from urllib.parse import urljoin, unquote
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

# Prefer shared project paths; fall back to the original EpiHydra layout.
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

try:
    from shared.paths import GSE_LIST, PROCESSED_DIR, RAW_DIR
except ImportError:  # pragma: no cover - standalone fallback
    GSE_LIST = [
        "GSE41037",
        "GSE40279",
        "GSE42861",
        "GSE27097",
        "GSE84727",
        "GSE111629",
        "GSE157131",
        "GSE197674",
    ]
    RAW_DIR = Path("/work/gr-fe/bryan/EpiHydra/data/01_raw")
    PROCESSED_DIR = Path("/work/gr-fe/bryan/EpiHydra/data/02_processed")

# Network resilience: NCBI's FTP-over-HTTPS endpoint occasionally drops the
# TLS handshake (SSLEOFError) or resets connections. These settings retry
# transient failures with backoff instead of aborting the whole run.
MAX_RETRIES = 3
RETRY_BACKOFF_SECONDS = 2  # doubles each attempt: 2s, 4s, 8s
REQUEST_TIMEOUT = (30, 300)  # (connect timeout, read timeout) in seconds
REQUEST_HEADERS = {
    "User-Agent": "curl/7.76.1",
    "Accept": "*/*",
}


# def _make_session() -> requests.Session:
#     """Build a requests Session with retry/backoff for transient errors."""
#     session = requests.Session()
#     retry = Retry(
#         total=MAX_RETRIES,
#         backoff_factor=RETRY_BACKOFF_SECONDS,
#         status_forcelist=[429, 500, 502, 503, 504],
#     )
#     adapter = HTTPAdapter(max_retries=retry)
#     session.mount("https://", adapter)
#     return session


# --------------------------------------------------------------------------
# HELPERS
# --------------------------------------------------------------------------
    
def series_matrix_url(gse: str) -> str:
    """Build the NCBI FTP URL for a GSE series matrix file.

    GEO buckets series into folders named after the accession with the
    last 3 digits replaced by 'nnn', e.g. GSE40279 -> GSE40nnn.
    """
    prefix = gse[:-3] + "nnn" 
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"
        f"{gse}_series_matrix.txt.gz"
    )

def curl_download(url: str, filepath: Path) -> None:
    """Download a file using curl, with retries and validation."""

    filepath.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = filepath.with_name(filepath.name + ".part")

    if tmp_path.exists():
        tmp_path.unlink()

    cmd = [
        "curl",
        "-fL",
        "--retry", "5",
        "--retry-delay", "5",
        "--retry-connrefused",
        "--connect-timeout", "60",
        "--max-time", "7200",
        "--progress-bar",
        "-o", str(tmp_path),
        url,
    ]

    try:
        subprocess.run(cmd, check=True)

        if not tmp_path.exists() or tmp_path.stat().st_size == 0:
            raise ValueError("Downloaded file is empty.")

        # Validate compressed files
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
    """Download the GEO series matrix, skipping existing files."""

    if filepath.exists() and filepath.stat().st_size > 0:
        print(f"{gse}: series matrix already exists, skipping download.")
        return

    url = series_matrix_url(gse)

    if gse == "GSE157131":
        url = url.replace(
            "GSE157131_series_matrix.txt.gz",
            "GSE157131-GPL21145_series_matrix.txt.gz"
        )

    print(f"{gse}: downloading series matrix...")
    curl_download(url, filepath)

def supplementary_url(gse: str) -> str:
    """Build the NCBI GEO supplementary directory URL."""

    prefix = gse[:-3] + "nnn"

    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{prefix}/{gse}/suppl/"
    )

def find_supplementary_beta(gse):
    """
    Find a suitable supplementary methylation matrix for a GEO series.

    Searches for filenames containing beta, normalised or normalized.
    Prioritises normalised beta files, then other normalised files,
    followed by beta files.
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

    # Construct GEO supplementary directory URL
    prefix = gse[:-3] + "nnn"

    supp_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{prefix}/{gse}/suppl/"
    )

    print(f"{gse}: searching supplementary files at {supp_url}")

    # Retrieve directory listing
    try:
        result = subprocess.run(
            [
                "curl",
                "-fLsS",
                "--connect-timeout", "60",
                "--max-time", "300",
                supp_url,
            ],
            capture_output=True,
            text=True,
            check=True,
        )

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Could not retrieve supplementary directory for {gse}. "
            f"curl exit code: {e.returncode}. Error: {e.stderr}"
        )

    # Extract filenames from HTML listing
    parser = HrefParser()
    parser.feed(result.stdout)

    include_terms = ("beta", "normalised", "normalized")

    exclude_terms = (
        "detection", "pvalue", "p-value", "intensity",
        "methylated", "unmethylated", "annotation",
        "manifest", "control", "idat",
    )

    valid_extensions = (
        ".txt", ".txt.gz",
        ".tsv", ".tsv.gz",
        ".csv", ".csv.gz",
    )

    candidates = []

    for href in parser.hrefs:
        filename = unquote(href).split("/")[-1]
        lower = filename.lower()

        if not lower.endswith(valid_extensions):
            continue

        if any(term in lower for term in exclude_terms):
            continue

        if not any(term in lower for term in include_terms):
            continue

        # Rank candidate files
        if "beta" in lower and (
            "normalised" in lower or "normalized" in lower
        ):
            score = 100
        elif "normalised" in lower or "normalized" in lower:
            score = 80
        elif "beta" in lower:
            score = 60
        else:
            score = 0

        candidates.append({
            "score": score,
            "filename": filename,
            "url": urljoin(supp_url, href),
        })

    if not candidates:
        raise RuntimeError(
            f"No beta or normalised/normalized candidate files "
            f"found for {gse}."
        )

    candidates.sort(key=lambda x: x["score"], reverse=True)

    print(f"{gse}: supplementary candidates found:")

    for candidate in candidates:
        print(f"  Score {candidate['score']}: {candidate['filename']}")

    # Select highest-ranked candidate
    best_score = candidates[0]["score"]

    best_candidates = [
        candidate for candidate in candidates
        if candidate["score"] == best_score
    ]

    if len(best_candidates) > 1:
        raise RuntimeError(
            f"Multiple equally ranked supplementary files found "
            f"for {gse}: "
            f"{[x['filename'] for x in best_candidates]}. "
            f"Please select the correct file manually."
        )

    selected = best_candidates[0]

    print(
        f"{gse}: selected supplementary file: "
        f"{selected['filename']}"
    )

    return selected["filename"], selected["url"]

def download_supplementary_beta(gse: str, directory: Path) -> Path:
    """Find and download supplementary beta values."""

    filename, url = find_supplementary_beta(gse)

    filepath = directory / filename

    if filepath.exists() and filepath.stat().st_size > 0:
        print(f"{gse}: supplementary beta file already exists, skipping.")
        return filepath

    print(f"{gse}: downloading supplementary beta values...")

    curl_download(url, filepath)

    return filepath

def parse_supplementary_beta(
    filepath: Path,
    sample_ids: list
) -> pd.DataFrame:
    """Parse supplementary beta values and align them with GEO sample IDs."""

    print(f"Parsing supplementary beta file: {filepath.name}")

    # Detect delimiter from the first line
    opener = gzip.open if filepath.name.lower().endswith(".gz") else open

    with opener(filepath, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()

    delimiter = "," if header.count(",") > header.count("\t") else "\t"

    df = pd.read_csv(
        filepath,
        sep=delimiter,
        index_col=0,
        compression="infer",
        low_memory=False,
    )

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    sample_ids = [str(s) for s in sample_ids]

    column_matches = sum(s in df.columns for s in sample_ids)
    index_matches = sum(s in df.index for s in sample_ids)

    if column_matches == 0 and index_matches == 0:
        raise ValueError(
            "Could not identify GEO sample IDs in supplementary beta file. "
            f"First columns: {list(df.columns[:5])}; "
            f"first index values: {list(df.index[:5])}"
        )

    if column_matches >= index_matches:
        # CpGs in rows, samples in columns
        expr = df

    else:
        # Samples in rows, CpGs in columns
        expr = df.T

    missing_samples = [
        s for s in sample_ids if s not in expr.columns
    ]

    if missing_samples:
        raise ValueError(
            f"{len(missing_samples)} GEO samples missing from beta file. "
            f"Examples: {missing_samples[:5]}"
        )

    # Preserve the metadata sample ordering
    expr = expr.loc[:, sample_ids]

    expr = expr.apply(pd.to_numeric, errors="coerce")

    print(
        f"Loaded supplementary methylation matrix: "
        f"{expr.shape[0]} CpGs x {expr.shape[1]} samples"
    )

    return expr

# def download_series_matrix(gse: str, filepath: Path) -> None:
#     """Download the series matrix file, skipping if it already exists."""
#     if filepath.exists() and filepath.stat().st_size > 0:
#         print(f"{gse}: file already exists at {filepath}, skipping download.")
#         return

#     filepath.parent.mkdir(parents=True, exist_ok=True)
#     url = series_matrix_url(gse)
#     if gse == 'GSE157131' : 
#         url = url.replace('GSE157131_series_matrix.txt.gz', 'GSE157131-GPL21145_series_matrix.txt.gzss')
#     tmp_path = filepath.with_suffix(filepath.suffix + ".part")

#     session = _make_session()
#     try:
#         response = session.get(url, headers=REQUEST_HEADERS, stream=True, timeout=REQUEST_TIMEOUT)
#         response.raise_for_status()

#         total_size = int(response.headers.get("content-length", 0))
#         downloaded = 0
#         bar_length = 40

#         with open(tmp_path, "wb") as f:
#             for chunk in response.iter_content(chunk_size=8192):
#                 f.write(chunk)
#                 downloaded += len(chunk)

#                 if total_size > 0:
#                     filled = int(bar_length * downloaded / total_size)
#                     bar = "█" * filled + "░" * (bar_length - filled)
#                     percent = downloaded / total_size * 100
#                     done_mb = downloaded / (1024 ** 2)
#                     total_mb = total_size / (1024 ** 2)
#                     print(
#                         f"\r{gse}: |{bar}| {percent:.1f}% "
#                         f"({done_mb:.1f}/{total_mb:.1f} MB)",
#                         end="",
#                         flush=True,
#                     )
#                 else:
#                     done_mb = downloaded / (1024 ** 2)
#                     print(f"\r{gse}: downloaded {done_mb:.1f} MB", end="", flush=True)

#         tmp_path.rename(filepath)
#         print(f"\nDownloaded to: {filepath}")

#     except Exception:
#         # Don't leave a half-written file behind to be mistaken for
#         # complete on a future run.
#         if tmp_path.exists():
#             tmp_path.unlink()
#         raise
#     finally:
#         session.close()


def parse_expression_matrix(filepath: Path) -> pd.DataFrame:
    """Read the actual data table (probes x samples) from the series matrix."""
    return pd.read_csv(filepath, sep="\t", comment="!", index_col=0, compression="gzip")


def parse_metadata(filepath: Path) -> dict:
    """Parse the '!'-prefixed metadata header lines into a flat dict.

    Sample_characteristics_ch1 lines are special-cased: the key is derived
    from the characteristic name (e.g. 'age', 'gender', 'tissue') rather
    than the literal 'Sample_characteristics_ch1' label, and values for all
    samples are joined with double spaces.
    """
    metadata = {}

    with gzip.open(filepath, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue

            line = line.strip().lstrip("!")
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue

            key, value = parts
            if key == "Sample_characteristics_ch1":
                if 'age' in parts[1] : 
                    key = "age"
                    cleaned_values = []
                    pattern = re.compile(r'(age[^:]*):\s*(\d+(?:\.\d+)?)', re.IGNORECASE)
                    for i in parts[1].split("\t"):
                        match = pattern.search(i)
                        if match:
                            key_part = match.group(1)
                            value = float(match.group(2))
                            if re.search(r'\bmonth', key_part, re.IGNORECASE):
                                value = round(value / 12 , 3)
                            elif re.search(r'\bweek', key_part, re.IGNORECASE):
                                value = round(value / 52 , 3)
                            elif re.search(r'\bday', key_part, re.IGNORECASE):
                                value = round(value / 365 , 3)
                            cleaned_values.append(str(value))
                        else :
                            cleaned_values.append('nan')
                    value = "  ".join(cleaned_values)
    
                else : 
                    key = parts[1].split("\t")[0].split(":")[0].split(" ")[0].strip().replace('"', "")
                    cleaned_values = []
                    for i in parts[1].split("\t"):
                        piece = i.split(":", 1)
                        # Expect "characteristic: value"; if a sample's entry
                        # doesn't contain a colon (blank/malformed field),
                        # fall back to the raw text instead of crashing.
                        v = piece[1] if len(piece) == 2 else piece[0]
                        cleaned_values.append(v.strip().replace('"', ""))
                value = "  ".join(cleaned_values)
            metadata[key] = value.strip('"')

    return metadata


def build_meta_df(metadata: dict) -> pd.DataFrame:
    """Build the sample-level metadata DataFrame, mirroring the notebook logic.

    Only pulls fields that are commonly present. If a dataset lacks one of
    'age'/'gender'/'tissue', that column is skipped with a warning rather
    than raising, so the script degrades gracefully across differing GSEs.
    """
    idx = [i.replace('"', "") for i in metadata["Sample_geo_accession"].split("\t")]
    organism = [i.replace('"', "") for i in metadata["Sample_organism_ch1"].split("\t")]
    platform = [i.replace('"', "") for i in metadata["Sample_platform_id"].split("\t")]

    data = {
        "ID": idx,
        "Organism": organism,
        "platform": platform,
    }

    for field, caster in (("age", np.float32), ("gender", str), ("tissue", str)):
        if field in metadata:
            values = metadata[field].split("  ")
            try:
                data[field] = [caster(v) for v in values]
            except ValueError:
                # Fall back to raw strings if casting (e.g. int(age)) fails
                data[field] = values
        else:
            print(f"  note: metadata field '{field}' not found, skipping column.")

    return pd.DataFrame(data)


def impute_missing_median(expr: pd.DataFrame) -> pd.DataFrame:
    """Median-impute missing values, per CpG (row), across samples.

    Any CpG site with a missing value gets that value replaced by the
    median of the non-missing values for that same CpG across all samples.
    A CpG that is missing for every sample has no median to impute with and
    is left as NaN (reported so it can be reviewed/dropped separately).
    """
    n_missing_before = expr.isna().sum().sum()
    if n_missing_before == 0:
        print("  no missing values found, skipping imputation.")
        return expr

    row_medians = expr.median(axis=1, skipna=True)
    # Note: DataFrame.fillna(Series) aligns the Series on *columns*, not
    # rows, so filling per-CpG (per-row) requires transposing, filling
    # (which then aligns row_medians against the transposed columns, i.e.
    # the original CpG index), and transposing back.
    expr_imputed = expr.T.fillna(row_medians).T

    n_missing_after = expr_imputed.isna().sum().sum()
    n_unfillable = expr_imputed.isna().any(axis=1).sum()

    print(
        f"  median imputation: {n_missing_before} missing values found, "
        f"{n_missing_before - n_missing_after} imputed."
    )
    if n_unfillable:
        print(
            f"  warning: {n_unfillable} CpG site(s) missing for ALL samples "
            f"could not be imputed and remain NaN."
        )

    return expr_imputed


def process_gse(gse: str) -> None:
    print(f"\n=== Processing {gse} ===")

    raw_directory = RAW_DIR / gse

    raw_filepath = raw_directory / f"{gse}_series_matrix.txt.gz"

    processed_filepath = PROCESSED_DIR / gse / "methylation.pkl"

    if processed_filepath.exists():
        print(
            f"{gse}: processed pickle already exists at "
            f"{processed_filepath}, skipping entirely."
        )
        return

    # Download series matrix
    download_series_matrix(gse, raw_filepath)

    # Parse metadata from the series matrix
    metadata = parse_metadata(raw_filepath)
    meta = build_meta_df(metadata)

    # Decide where to obtain methylation measurements
    threshold = 10 * 1024 * 1024

    matrix_size = raw_filepath.stat().st_size

    print(
        f"{gse}: series matrix size = "
        f"{matrix_size / (1024**2):.2f} MB"
    )

    if matrix_size < threshold:

        print(
            f"{gse}: series matrix below 10 MB. "
            "Searching for supplementary beta values."
        )

        beta_filepath = download_supplementary_beta(
            gse,
            raw_directory,
        )

        expr = parse_supplementary_beta(
            beta_filepath,
            meta["ID"].tolist(),
        )

    else:

        print(
            f"{gse}: series matrix above 10 MB. "
            "Using series matrix methylation data."
        )

        df = parse_expression_matrix(raw_filepath)

        expr = df.loc[:, meta["ID"]]

    # Existing processing
    expr = impute_missing_median(expr)

    processed_filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(processed_filepath, "wb") as file:
        pickle.dump(
            {
                "meta": meta,
                "CpG": expr,
            },
            file,
        )

    print(
        f"{gse}: saved processed data to {processed_filepath}"
    )


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main():
    succeeded = []
    failed = []

    for gse in GSE_LIST:
        try:
            process_gse(gse)
            succeeded.append(gse)
        except Exception as e:
            print(f"\n{gse}: FAILED — {type(e).__name__}: {e}")
            traceback.print_exc()
            failed.append(gse)
            continue

    print("\n=== Summary ===")
    print(f"Succeeded ({len(succeeded)}): {succeeded}")
    if failed:
        print(f"Failed ({len(failed)}): {failed}")
        print("Re-run the script to retry only the failed datasets "
              "(completed ones are skipped automatically).")


if __name__ == "__main__":
    main()