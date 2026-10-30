#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Fetch the JPL source kernels for the minor-body pack and write the provenance ledger.

Sources (all JPL, no Swiss data):

* Ceres, Pallas, Juno, Vesta: `sb441-n16.bsp`, the 16 asteroid perturbers
  integrated together with DE441 (JPL IOM 392R-21-005), SPK type 2, centre
  Sun, JD -1200525.5 .. 5008242.5.  Read remotely with HTTP Range requests
  (tools/spk.py `HttpRangeFile`); only the needed records are cached.
* Chiron, Pholus, Eris, Eros: JPL Horizons small-body SPK files (type 21,
  centre Sun), generated through the Horizons API with EPHEM_TYPE='SPK'.
  Horizons refuses SPK and vector requests for these bodies before
  1599-12-11 TDB (its small-body integration limit), so that is where their
  coverage starts.

    python tools/fetch_minor.py            # download missing files, write ledger
    python tools/fetch_minor.py --ledger-only

The ledger (`docs/MINOR_PROVENANCE.json`, committed) records request
parameters, retrieval time, file sha256, Horizons solution ID (JPL#...),
solution date, integration source tags and licence/attribution notes.
`tools/fit_minor.py` copies it into the pack manifest.
"""
from __future__ import annotations

import argparse
import base64
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spk import SPK  # noqa: E402

DEFAULT_DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(REPO, "docs", "MINOR_PROVENANCE.json")

HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"
SB441_URL = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp"
SB441_DOC = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf"

# Horizons small-body SPK request span. Start = the earliest epoch Horizons
# accepts for these bodies ("No ephemeris ... prior to A.D. 1599-DEC-11").
HZ_START = "1599-12-11"
HZ_STOP = "2401-04-01"

HORIZONS_BODIES = [(2060, "chiron"), (5145, "pholus"), (136199, "eris"), (433, "eros")]
SB441_BODIES = [(1, "ceres"), (2, "pallas"), (3, "juno"), (4, "vesta")]

LICENSE_NOTE = (
    "JPL Solar System Dynamics data (NASA/JPL-Caltech). No licence terms are attached to the "
    "files; JPL asks for credit to the Horizons On-Line Ephemeris System / SSD group. The "
    "Horizons SPK comments state the files are informal releases, not for spacecraft-safety use. "
    "Confirm the attribution wording with JPL SSD before public redistribution of derived packs."
)


def hz_spk_path(minor_dir: str, number: int, name: str) -> str:
    return os.path.join(minor_dir, f"hz_{number}_{name}.bsp")


def horizons_spk(number: int, out: str) -> dict:
    params = {"format": "json", "COMMAND": f"'{number};'", "EPHEM_TYPE": "'SPK'", "OBJ_DATA": "'YES'",
              "START_TIME": f"'{HZ_START}'", "STOP_TIME": f"'{HZ_STOP}'"}
    url = HORIZONS_API + "?" + urllib.parse.urlencode(params)
    r = subprocess.run(["curl", "-sSf", "--max-time", "600", url], check=True, capture_output=True)
    d = json.loads(r.stdout)
    if "spk" not in d:
        raise RuntimeError(f"Horizons refused SPK for {number}: {d.get('error')}\n{d.get('result', '')[-1500:]}")
    tmp = out + ".tmp"
    with open(tmp, "wb") as f:
        f.write(base64.b64decode(d["spk"]))
    os.replace(tmp, out)
    with open(out + ".result.txt", "w") as f:
        f.write(d["result"])
    with open(out + ".request.json", "w") as f:
        json.dump({"api": HORIZONS_API, "params": params,
                   "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                   "spk_file_id": d.get("spk_file_id")}, f, indent=1)
    return d


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _grab(pattern: str, text: str) -> str | None:
    m = re.search(pattern, text, re.M)
    return m.group(1).strip() if m else None


def horizons_entry(path: str, number: int, name: str) -> dict:
    spk = SPK.open(path)
    com = spk.comments()
    res = open(path + ".result.txt").read() if os.path.exists(path + ".result.txt") else ""
    req = {}
    if os.path.exists(path + ".request.json"):
        req = json.load(open(path + ".request.json"))
    else:
        # files fetched before this script existed: same parameters, time from the comments
        req = {"api": HORIZONS_API,
               "params": {"COMMAND": f"'{number};'", "EPHEM_TYPE": "'SPK'", "OBJ_DATA": "'YES'",
                          "START_TIME": f"'{HZ_START}'", "STOP_TIME": f"'{HZ_STOP}'"},
               "retrieved_utc": None}
    segs = spk.segments
    return {
        "body": name.capitalize(), "number": number,
        "file": os.path.basename(path), "size_bytes": os.path.getsize(path), "sha256": sha256(path),
        "kind": "JPL Horizons small-body SPK (type 21)",
        "request": req,
        "execution_date_pasadena": _grab(r"Execution date\s*:\s*(.+)$", com),
        "target_line": _grab(r"Target body\s*:\s*(.+)$", com),
        "solution_id": _grab(r"source:\s*Horizons_SPK:(JPL#\d+)", com),
        "solution_date": _grab(r"Soln\.date:\s*(\S+)", res),
        "observations": _grab(r"# obs:\s*(.+)$", res),
        "spk_target_id": segs[0].target, "spk_center_id": segs[0].center, "spk_frame": segs[0].frame,
        "coord_center": _grab(r"Coord center\s*:\s*(.+)$", com),
        "small_body_perturbers": _grab(r"Small-body perts:\s*(.+)$", com),
        "start_time": _grab(r"Start time\s*:\s*(.+)$", com),
        "stop_time": _grab(r"Stop  time\s*:\s*(.+)$", com),
        "coverage_jd_tdb": [min(s.start_jd for s in segs), max(s.end_jd for s in segs)],
        "n_segments": len(segs),
        "force_model": ("planets, Moon, 16 largest asteroids (SB441-N16), PPN solar relativity, Earth "
                        "oblateness, non-gravitational terms of the solution (from the SPK comments)"),
        "horizons_span_limit": ("Horizons refuses ephemerides for this body before A.D. 1599-DEC-11 TDB "
                                "('No ephemeris for target ... prior to'), so no JPL data exists before"),
        "license": LICENSE_NOTE,
    }


def sb441_entry(minor_dir: str) -> dict:
    cache = os.path.join(minor_dir, "sb441-n16_cache")
    spk = SPK.open(SB441_URL, cache)
    f = spk.daf.f
    bodies = {}
    for num, name in SB441_BODIES:
        segs = spk.segments_for(2000000 + num, 10)
        bodies[name.capitalize()] = {
            "spk_target_id": 2000000 + num, "spk_center_id": 10,
            "segments_jd_tdb": [[s.start_jd, s.end_jd] for s in segs],
            "record_days": segs[0].interval_days, "n_coeff": segs[0].ncoef, "type": segs[0].data_type}
    return {
        "file": "sb441-n16.bsp", "url": SB441_URL, "size_bytes": f.size,
        "last_modified": getattr(f, "last_modified", ""),
        "how_read": "HTTP Range requests, 1 MiB chunk cache (data/sources/minor/sb441-n16_cache)",
        "sha256": None,
        "integrity": ("not hashed (range-read); DAF/type-2 directory checks in tools/spk.py, record grid "
                      "checked, record-boundary continuity <= 3e-7 km"),
        "kind": "JPL DE441 small-body perturber SPK (type 2, Chebyshev, 32-day records)",
        "documentation": SB441_DOC,
        "reference": "Park, Folkner, Williams, Boggs (2021), JPL IOM 392R-21-005; DE440/DE441 paper AJ 161:105",
        "ifname": spk.daf.ifname,
        "bodies": bodies,
        "note": ("Asteroid trajectories integrated simultaneously with the DE441 planets (the same ones "
                 "used as perturbers in DE441 and in every Horizons small-body integration)."),
        "license": LICENSE_NOTE,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--ledger-only", action="store_true")
    ap.add_argument("--ledger", default=LEDGER)
    args = ap.parse_args()
    minor_dir = os.path.join(args.data, "sources", "minor")
    os.makedirs(minor_dir, exist_ok=True)
    entries = []
    for num, name in HORIZONS_BODIES:
        p = hz_spk_path(minor_dir, num, name)
        if not os.path.exists(p):
            if args.ledger_only:
                raise SystemExit(f"missing {p}")
            print(f"fetching Horizons SPK for {num} {name}", flush=True)
            horizons_spk(num, p)
        entries.append(horizons_entry(p, num, name))
    ledger = {
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tool": "tools/fetch_minor.py",
        "clean_room": "All minor-body data come from JPL. Swiss Ephemeris is used only as a black-box oracle "
                      "in harness/check_minor.py, never as a data source.",
        "sb441": sb441_entry(minor_dir),
        "horizons": entries,
    }
    with open(args.ledger, "w") as f:
        json.dump(ledger, f, indent=1)
    print(f"wrote {args.ledger}")
    for e in entries:
        print(f"  {e['body']:<7} {e['solution_id']} soln {e['solution_date']} "
              f"JD {e['coverage_jd_tdb'][0]:.1f}..{e['coverage_jd_tdb'][1]:.1f} {e['file']}")


if __name__ == "__main__":
    main()
