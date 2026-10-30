# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Move Swiss-Ephemeris-derived test data out of the public tree into the private oracle corpus.

Every file that contains outputs of the Swiss Ephemeris (pyswisseph, used as a black-box
oracle) lives in a private corpus that mirrors the repository layout:

    $FE_ORACLE_DIR/<repo-relative path>        (default FE_ORACLE_DIR = <repo>/private/oracle)

`private/` is gitignored and holds every Swiss-derived artefact. Layout (one convention):

    private/oracle/<repo-relative path>   test fixtures and harness reports with Swiss outputs
                                          ($FE_ORACLE_DIR; mirrors the public tree)
    private/stars/                        the optional star-name compatibility dataset (probe,
                                          identifications, star_compat_aliases.tsv; built by
                                          tools/probe_swiss_stars.py, identify_swiss_stars.py,
                                          build_star_compat.py; loaded at runtime through
                                          FRENCHEPHEMERIS_STAR_COMPAT)

Rule: data used only to *test* against the oracle goes to private/oracle/ at its mirrored
path; Swiss-derived data a deployment may *load* lives in a named dataset under private/.

Three kinds of files are handled:

* MOVE  - the whole file is Swiss output: moved to the corpus, removed from the public tree.
* SPLIT - Swiss and non-Swiss data (ERFA, JPL Horizons, NASA, IOTA) share one file: the full
          file goes to the corpus, and a public copy with the Swiss parts stripped stays in place.
* AUTO  - any other JSON under crates/fe-core/tests/fixtures/<area>/ whose area has a gen_*.py
          that imports swisseph is treated as MOVE (with a warning), unless listed in PUBLIC.
          This catches fixtures added by branches merged later.

The script is idempotent: re-run it after merging other branches or after regenerating
fixtures. A newer public file always replaces the corpus copy.

    .venv/bin/python tools/privatize_oracle_fixtures.py            # perform the moves
    .venv/bin/python tools/privatize_oracle_fixtures.py --dry-run  # show what would happen
    .venv/bin/python tools/privatize_oracle_fixtures.py --check    # verify the public tree only

Generators import `oracle_file` / `write_split` from here so they write straight to the corpus.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIX = "crates/fe-core/tests/fixtures/"


def private_dir() -> Path:
    return REPO / "private"


def oracle_dir() -> Path:
    d = os.environ.get("FE_ORACLE_DIR")
    return Path(d) if d else REPO / "private" / "oracle"


def rel(path) -> str:
    return str(Path(path).resolve().relative_to(REPO))


def oracle_file(public_path) -> str:
    """Corpus path mirroring `public_path` (a path inside the repo); parent dirs are created."""
    p = oracle_dir() / rel(public_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p)


# ------------------------------------------------------------------------------- the inventory
# Whole files that are Swiss output.
MOVE = [
    FIX + "ayanamsa/ayanamsa.json",
    FIX + "ayanamsa_all/ayanamsa_all.json",
    FIX + "eclipses/eclipses.json",          # NASA canon columns and Swiss columns share each row
    FIX + "eclipses/visibility.json",        # Swiss how / when_loc flags (gen_eclipse_visibility.py)
    FIX + "events/events.json",
    FIX + "fictitious/fictitious.json",
    FIX + "heliacal/heliacal.json",
    FIX + "houses/gauquelin.json",
    FIX + "houses/house_pos.json",
    FIX + "houses/house_pos_ext.json",
    FIX + "houses/house_speeds.json",
    FIX + "houses/houses.json",
    FIX + "houses/houses_ext.json",
    FIX + "houses/houses_sidereal.json",
    FIX + "occultations/occultations.json",  # Swiss columns + IOTA predictions (redistribution unconfirmed)
    FIX + "points/points.json",
    FIX + "quickwins/quickwins.json",
    FIX + "sidtime/sidtime.json",
    FIX + "stars/swiss.json",
    FIX + "stars/swiss_names.json",
    "harness/swiss_report.json",
    "harness/swiss_report_compat.json",
    "tools/star_names_swiss.tsv",            # names/strings returned by swe.fixstar2_ut
    "tools/star_identifications.tsv",        # separations from Swiss's star places
]

# MOVE entries that belong to a private dataset rather than the mirrored oracle corpus:
# public path -> destination relative to private/.
DATASETS = {
    "tools/star_names_swiss.tsv": "stars/star_names_swiss.tsv",
    "tools/star_identifications.tsv": "stars/star_identifications.tsv",
}


def private_dest(r: str) -> Path:
    return private_dir() / DATASETS[r] if r in DATASETS else oracle_dir() / r

# Files mixing Swiss and non-Swiss data: what to strip from the public copy.
#   drop: top-level keys; meta: keys removed from the top-level "meta" dict;
#   row_fields: {list key: [fields removed from every row]}; set: top-level keys overwritten;
#   deep: remove, at any depth, dict keys matching this regex.
SPLIT = {
    FIX + "frames/frames.json": dict(drop=["swiss"], meta=["swisseph"]),
    FIX + "horizons/horizons.json": dict(drop=["swiss_jplhora", "swiss_jplhor"]),
    FIX + "obliquity/obliquity.json": dict(drop=["swiss"], meta=["swisseph"]),
    FIX + "time/time.json": dict(drop=["julday", "swe_deltat"], meta=["swisseph", "tid_acc"]),
    FIX + "eclipses/extras.json": dict(drop=["swiss_solar", "swiss_lunar"]),
    FIX + "events/horizons.json": dict(row_fields={"cases": ["t_swiss"]}),
    "harness/horizons_report.json": dict(deep=r"(?i)^swe|swiss"),
    "harness/minor_report.json": dict(deep=r"(?i)^swe|swiss"),
    "docs/CHIRON_STUDY.json": dict(deep=r"(?i)^swe|swiss"),
}

# Fixtures in an area with a Swiss generator that are nevertheless independent of Swiss.
PUBLIC = {
    FIX + "stars/erfa.json",                 # ERFA (pyerfa) only, gen_erfa.py
}

# JSON keys naming Swiss values, or Swiss named as the source of the data.
SWISS_MARKERS = re.compile(r'"[A-Za-z0-9_]*(?:swiss|swisseph|swe_)[A-Za-z0-9_]*"\s*:|pyswisseph|Swiss Ephemeris \d', re.I)


def strip(data, spec):
    data = json.loads(json.dumps(data))
    for k in spec.get("drop", []):
        data.pop(k, None)
    if "meta" in spec and isinstance(data.get("meta"), dict):
        for k in spec["meta"]:
            data["meta"].pop(k, None)
    for key, fields in spec.get("row_fields", {}).items():
        for row in data.get(key, []):
            for fld in fields:
                row.pop(fld, None)
    data.update(spec.get("set", {}))
    if "deep" in spec:
        rx = re.compile(spec["deep"])

        def walk(x):
            if isinstance(x, dict):
                return {k: walk(v) for k, v in x.items() if not rx.search(k)}
            if isinstance(x, list):
                return [walk(v) for v in x]
            return x
        data = walk(data)
    return data


def dump_like(data, original_text: str) -> str:
    if "\n " in original_text[:200]:
        return json.dumps(data, indent=1) + "\n"
    if '": ' in original_text[:2000] or '", "' in original_text[:2000]:
        return json.dumps(data)
    return json.dumps(data, separators=(",", ":"))


def write_split(public_path, full_data, indent=None) -> None:
    """For a SPLIT file: write `full_data` to the corpus and the stripped copy to `public_path`."""
    spec = SPLIT[rel(public_path)]
    text = json.dumps(full_data, indent=indent) + "\n" if indent else json.dumps(full_data, separators=(",", ":"))
    Path(oracle_file(public_path)).write_text(text)
    Path(public_path).write_text(dump_like(strip(full_data, spec), text))
    print(f"wrote {oracle_file(public_path)} (full) and {public_path} (Swiss parts stripped)")


def save_report(path, data, indent=1, private=False) -> None:
    """Write a harness report: SPLIT paths are split, MOVE paths (or `private=True`, or any
    repo path whose content names Swiss values) go to the oracle corpus only."""
    path = Path(path).resolve()
    try:
        r = rel(path)
    except ValueError:  # outside the repo: the caller's choice
        path.write_text(json.dumps(data, indent=indent, default=float))
        return
    if r in SPLIT:
        write_split(path, json.loads(json.dumps(data, default=float)), indent=indent)
        return
    text = json.dumps(data, indent=indent, default=float)
    if private or r in MOVE or SWISS_MARKERS.search(text):
        o = oracle_file(path)
        Path(o).write_text(text)
        print(f"wrote {o} (Swiss-derived: private oracle corpus)")
    else:
        path.write_text(text)


# ------------------------------------------------------------------------------- discovery
def swiss_areas() -> set[str]:
    out = set()
    for gen in (REPO / FIX).glob("*/*.py"):
        if re.search(r"^\s*import swisseph", gen.read_text(), re.M):
            out.add(gen.parent.name)
    return out


def auto_moves() -> list[str]:
    known = set(MOVE) | set(SPLIT) | PUBLIC
    out = []
    for area in sorted(swiss_areas()):
        for p in sorted((REPO / FIX / area).glob("*.json")):
            r = rel(p)
            if r not in known:
                out.append(r)
    return out


def public_files():
    skip = ("private/", "target/", ".venv/", "data/", ".git/")
    for p in REPO.rglob("*"):
        r = str(p.relative_to(REPO))
        if p.is_file() and not r.startswith(skip) and "__pycache__" not in r:
            yield r, p


def check() -> int:
    """Swiss-output markers left in public data files (JSON/TSV). Code and prose are not data."""
    bad = []
    for r, p in public_files():
        if p.suffix not in (".json", ".tsv", ".txt", ".csv"):
            continue
        try:
            txt = p.read_text(errors="replace")
        except OSError:
            continue
        stripped_differs = False
        if r in SPLIT:
            try:
                stripped_differs = strip(json.loads(txt), SPLIT[r]) != json.loads(txt)
            except ValueError:
                stripped_differs = True
        if r in MOVE or stripped_differs or SWISS_MARKERS.search(txt) or (
                p.suffix == ".tsv" and re.search(r"(?i)swiss", txt[:2000])):
            bad.append(r)
    for r in auto_moves():
        if (REPO / r).exists() and r not in bad:
            bad.append(r)
    for r in bad:
        print(f"SWISS DATA IN PUBLIC TREE: {r}")
    print("check:", "FAIL" if bad else "ok (no Swiss-output data files in the public tree)")
    return 1 if bad else 0


def unguarded_tests(moved: list[str]) -> None:
    """Warn about Rust tests that read moved fixtures by their public path."""
    tests = REPO / "crates/fe-core/tests"
    for t in tests.glob("*.rs"):
        src = t.read_text()
        for r in moved:
            name = r[len(FIX):] if r.startswith(FIX) else None
            if name and f"tests/fixtures/{name}" in src.replace('"/tests', 'tests'):
                print(f"WARNING: {t.relative_to(REPO)} reads {name} by its public path; use "
                      f"common::oracle_path(\"{name}\") + require_oracle!")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true", help="only verify the public tree")
    a = ap.parse_args()
    if a.check:
        return check()
    od = oracle_dir()
    print(f"oracle corpus: {od}")
    moved, nbytes = [], 0
    auto = auto_moves()
    for r in auto:
        print(f"WARNING: unclassified fixture in a Swiss area, treated as MOVE: {r} "
              f"(add it to MOVE, SPLIT or PUBLIC in {rel(__file__)})")
    for r in MOVE + auto:
        src, dst = REPO / r, private_dest(r)
        if not src.exists():
            print(f"  ok      {r}" + ("" if dst.exists() else "  (absent in both trees)"))
            continue
        nbytes += src.stat().st_size
        print(f"  move    {r} ({src.stat().st_size:,} B)")
        if not a.dry_run:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
        moved.append(r)
    for r, spec in SPLIT.items():
        src, dst = REPO / r, od / r
        if not src.exists():
            print(f"  missing {r}")
            continue
        text = src.read_text()
        data = json.loads(text)
        public = strip(data, spec)
        if public == data:
            print(f"  ok      {r} (public copy already stripped" + ("" if dst.exists() else "; NO oracle copy") + ")")
            continue
        nbytes += len(text.encode()) - len(dump_like(public, text).encode())
        print(f"  split   {r}")
        if not a.dry_run:
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text(text)
            src.write_text(dump_like(public, text))
        moved.append(r)
    print(f"{len(moved)} file(s) {'would be ' if a.dry_run else ''}privatized, "
          f"{nbytes / 1e6:.1f} MB of Swiss-derived data removed from the public tree")
    unguarded_tests([r for r in MOVE + auto])
    gi = REPO / ".gitignore"
    if "/private/" not in gi.read_text().split():
        print("WARNING: /private/ is not in .gitignore")
    if not a.dry_run:
        return check()
    return 0


if __name__ == "__main__":
    sys.exit(main())
