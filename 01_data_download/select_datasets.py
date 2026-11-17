#!/usr/bin/env python3
"""Select GEO datasets per Overview slide 3.

Checks: n>100, clinical outcome in metadata, 450K/EPIC platform (CpG proxy).
Beta values are verified later at download — metadata screen stays fast.
Age is reported (needed for MAE) but is not an Overview fail.

Cache: data/selection/selection_cache.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(_ROOT), str(Path(__file__).resolve().parent)]

from shared.paths import (  # noqa: E402
    CACHE_FILE, GSE_LIST, MIN_N, OUTCOME_WORDS, PLATFORM_OK,
    PROCESSED_DIR, RAW_DIR, SELECTION_DIR, TTL_DAYS,
)
import geo_data_download as geo  # noqa: E402

geo.RAW_DIR, geo.PROCESSED_DIR = RAW_DIR, PROCESSED_DIR


def _first(s: str) -> str:
    return s.split("\t", 1)[0].strip().strip('"')


def _header(path: Path) -> tuple[str, str, list[str]]:
    title, platform, keys, seen = "", "", [], set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_title\t"):
                title = _first(line.split("\t", 1)[1])
            elif line.startswith("!Sample_platform_id\t") and not platform:
                platform = _first(line.split("\t", 1)[1])
            elif line.startswith("!Sample_characteristics_ch1\t"):
                k = _first(line.split("\t", 1)[1]).split(":", 1)[0].strip().lower()
                if k and k not in seen:
                    seen.add(k)
                    keys.append(k)
    return title, platform, keys


def _check(gse: str, force_web: bool = False) -> dict:
    d = RAW_DIR / gse
    d.mkdir(parents=True, exist_ok=True)
    matrix = d / f"{gse}_series_matrix.txt.gz"
    huge = matrix.exists() and matrix.stat().st_size >= 50 * 1024 * 1024
    geo.download_series_matrix(gse, matrix, force=force_web and not huge)

    title, platform, keys = _header(matrix)
    meta = geo.build_meta_df(geo.parse_metadata(matrix))
    n = len(meta)
    ages = pd.to_numeric(meta["age"], errors="coerce") if "age" in meta.columns else None
    age_n = int(ages.notna().sum()) if ages is not None else 0

    outcomes = [k for k in keys if any(w in k for w in OUTCOME_WORDS)]
    ok_n, ok_out, ok_plat = n > MIN_N, bool(outcomes), platform in PLATFORM_OK
    # Overview: beta checked at download time (series/suppl), not in this fast pass.
    passed = ok_n and ok_out and ok_plat

    fail = []
    if not ok_n:
        fail.append("n<=100")
    if not ok_out:
        fail.append("no clinical outcome")
    if not ok_plat:
        fail.append(f"platform {platform or '?'} (need 450K/EPIC)")

    return {
        "gse": gse, "n": n, "age_n": age_n, "platform": platform,
        "outcomes": outcomes, "has_age": age_n > 0, "pass": passed,
        "fail": fail, "title": title,
    }


def _fresh(cache: dict | None, need: list[str]) -> bool:
    if not cache or "at" not in cache:
        return False
    try:
        at = datetime.fromisoformat(cache["at"])
    except ValueError:
        return False
    if at.tzinfo is None:
        at = at.replace(tzinfo=timezone.utc)
    if (datetime.now(timezone.utc) - at).total_seconds() > TTL_DAYS * 86400:
        return False
    if set(cache.get("catalog") or []) != set(GSE_LIST):
        return False
    return set(need) <= {r["gse"] for r in cache.get("rows") or []}


def _write(rows: list[dict]) -> dict:
    SELECTION_DIR.mkdir(parents=True, exist_ok=True)
    selected = [r["gse"] for r in rows if r["pass"]]
    with_age = [r["gse"] for r in rows if r["pass"] and r["has_age"]]
    cache = {
        "at": datetime.now(timezone.utc).isoformat(),
        "ttl_days": TTL_DAYS,
        "catalog": list(GSE_LIST),
        "selected": selected,
        "selected_with_age": with_age,  # soft filter for MAE, not Overview
        "rows": rows,
    }
    CACHE_FILE.write_text(json.dumps(cache, indent=2))
    (SELECTION_DIR / "selected_gse.txt").write_text("\n".join(selected) + "\n")
    print(f"cache → {CACHE_FILE}\nPASS ({len(selected)}): {selected}")
    return cache


def get_selected_datasets(
    require_age: bool = False,
    force_refresh: bool = False,
    gse_list: list[str] | None = None,
) -> list[str]:
    need = list(gse_list or GSE_LIST)
    cache = json.loads(CACHE_FILE.read_text()) if CACHE_FILE.exists() else None
    if force_refresh or not _fresh(cache, need):
        print("refreshing from GEO…")
        known = {r["gse"]: r for r in (cache or {}).get("rows") or []}
        for gse in need:
            print(f"  {gse}")
            known[gse] = _check(gse, force_web=True)
        rows = [known[g] for g in GSE_LIST if g in known]
        cache = _write(rows)
    key = "selected_with_age" if require_age else "selected"
    out = list(cache[key])
    return [g for g in out if g in set(need)] if gse_list else out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gse", nargs="*")
    p.add_argument("--force-refresh", action="store_true")
    p.add_argument("--require-age", action="store_true", help="also need age (MAE)")
    args = p.parse_args(argv)
    selected = get_selected_datasets(
        require_age=args.require_age, force_refresh=args.force_refresh, gse_list=args.gse
    )
    for r in json.loads(CACHE_FILE.read_text())["rows"]:
        print(
            f"{r['gse']:12} {'PASS' if r['pass'] else 'FAIL':4} "
            f"n={r['n']:<4} age_n={r['age_n']:<4} {r['platform']:10} {r['fail'] or ''}"
        )
    print("selected:", selected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
