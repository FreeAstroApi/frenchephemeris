#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Extended packs vs JPL Horizons (DE441) geometric vectors.

Horizons serves DE441, but only for 9999 BC Mar 21 … AD 9999 Dec 31 (both
observer and vector tables refuse earlier/later epochs), so this covers
−9998 … −1 and 2401 … 9999; the rest of the extended range (−13200 … −9999 and
10000 … 17191) is verified against DE441 itself inside the fitter (every record,
held-out points; see the pack manifest `fit.bodies`).

For random TDB epochs, fetches the geometric geocentric ICRF position of each
body (VEC_CORR NONE, CENTER 500@399, km) and compares with the pack evaluated
through the Rust reader (`pack_eval` example, main + extended packs chained):
3-D difference (km) and the angle seen from the geocentre (arcsec).

    .venv/bin/python harness/horizons_ext.py [--per-era 30] [--json out.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
import tempfile
import urllib.parse
from pathlib import Path

import numpy as np

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
REPO = Path(__file__).resolve().parents[1]
CACHE = Path(__file__).parent / "cache" / "horizons_ext"
DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
PACKS = ",".join(p for p in (os.path.join(DATA, "packs", n) for n in
                 ("fepk-de441-de442-v1.bin", "fepk-de441-ext-past-v1.bin", "fepk-de441-ext-future-v1.bin"))
                 if os.path.exists(p))
# (name, our NAIF target for pack_eval, Horizons command)
TARGETS = [("Sun", 10, "10"), ("Moon", 301, "301"), ("Mercury", 1, "1"), ("Venus", 2, "2"),
           ("Mars", 4, "4"), ("Jupiter", 5, "5"), ("Saturn", 6, "6"), ("Uranus", 7, "7"),
           ("Neptune", 8, "8"), ("Pluto", 9, "9"), ("EMB-Earth", 3, "3")]
# TDB JD ranges inside Horizons' limits
ERAS = {"-9998..-1": (-1930100.0, 1720600.0), "2401..9999": (2598200.0, 5373100.0)}


def fetch(cmd: str, jds: list[float]) -> np.ndarray:
    params = {"format": "json", "COMMAND": f"'{cmd}'", "EPHEM_TYPE": "'VECTORS'", "CENTER": "'500@399'",
              "TLIST_TYPE": "'JD'", "TIME_TYPE": "'TDB'", "TLIST": "'" + " ".join(f"{j:.6f}" for j in jds) + "'",
              "VEC_TABLE": "'1'", "OUT_UNITS": "'KM-S'", "CSV_FORMAT": "'YES'", "REF_PLANE": "'FRAME'",
              "REF_SYSTEM": "'ICRF'", "VEC_CORR": "'NONE'", "OBJ_DATA": "'NO'"}
    key = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:24]
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{cmd}-{key}.json"
    if path.exists():
        data = json.loads(path.read_text())
    else:
        url = API + "?" + urllib.parse.urlencode(params)
        out = subprocess.run(["curl", "-sSf", "--max-time", "300", url], check=True, capture_output=True)
        data = json.loads(out.stdout)
        path.write_text(json.dumps(data))
    body = data["result"].split("$$SOE")[1].split("$$EOE")[0]
    rows = []
    for line in body.strip().splitlines():
        c = [x.strip() for x in line.split(",")]
        rows.append([float(c[0]), float(c[2]), float(c[3]), float(c[4])])
    a = np.array(rows)
    assert len(a) == len(jds) and np.allclose(a[:, 0], jds, atol=1e-5)
    return a[:, 1:]


def ours(target: int, jds: list[float]) -> np.ndarray:
    with tempfile.TemporaryDirectory() as d:
        fi, fo = os.path.join(d, "t.bin"), os.path.join(d, "o.bin")
        np.asarray(jds, "<f8").tofile(fi)
        subprocess.run(["cargo", "run", "-q", "-p", "fe-core", "--example", "pack_eval", "--", PACKS,
                        str(target), "399", fi, fo], check=True, cwd=REPO)
        return np.fromfile(fo, "<f8").reshape(len(jds), 6)[:, :3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-era", type=int, default=30)
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--json")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    rep = {}
    for name, tgt, cmd in TARGETS:
        rep[name] = {}
        for era, (lo, hi) in ERAS.items():
            jds = sorted(round(rng.uniform(lo, hi), 6) for _ in range(a.per_era))
            h = fetch(cmd, jds)
            if name == "EMB-Earth":
                tgt_eval = 3  # EMB wrt Earth = −(Earth wrt EMB): tests the Earth alias
            else:
                tgt_eval = tgt
            o = ours(tgt_eval, jds)
            if np.isnan(o).all():
                print(f"{name:<10} {era:<11} (pack not present)")
                continue
            d = np.linalg.norm(o - h, axis=1)
            ang = d / np.linalg.norm(h, axis=1) * 206264.806
            rep[name][era] = {"n": len(jds), "max_km": float(d.max()), "p50_km": float(np.median(d)),
                              "max_arcsec": float(ang.max()), "min_dist_km": float(np.linalg.norm(h, axis=1).min())}
            print(f"{name:<10} {era:<11} max {d.max():9.4f} km  p50 {np.median(d):9.4f} km  "
                  f"max {ang.max():.6f}\"", flush=True)
    if a.json:
        Path(a.json).write_text(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
