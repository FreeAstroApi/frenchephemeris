#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Generate pack reader fixtures with jplephem (an independent SPK reader).

Writes `pack_states.json` beside this script: for many random TDB epochs and
body pairs, position (km) and velocity (km/day) evaluated directly from the
JPL kernels:

* DE442 (`data/sources/de442.bsp`) for JD >= 2287184.5 (DE442 start = the seam),
* DE441 part-1 before that, read through a *sparse* local copy assembled from
  the HTTP-range cache in `data/sources/de441_cache/` (run tools/build_pack.py
  for the full pack first so the needed ranges are cached).

Chained pairs (e.g. Moon wrt Earth, Sun wrt Earth) are formed here by adding
the stored JPL segments with numpy, independently of the Rust chaining code.

Usage:  python gen_pack_fixtures.py [--data <repo>/data]
"""
import argparse
import json
import os
import sys

import numpy as np
from jplephem.spk import SPK as JplSPK

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tools"))
from spk import HttpRangeFile  # noqa: E402  (only to materialise the sparse DE441 copy)

SEAM = 2287184.5
DE441_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp"

# path to SSB for each body, as (target, center) segments
PATH = {0: [], 3: [(3, 0)], 301: [(301, 3), (3, 0)], 399: [(399, 3), (3, 0)]}
for b in (1, 2, 4, 5, 6, 7, 8, 9, 10):
    PATH[b] = [(b, 0)]

PAIRS = ([(t, 0) for t in range(1, 11)] + [(301, 3), (399, 3)]
         + [(399, 0), (301, 0), (301, 399), (10, 399), (4, 399), (5, 301), (399, 10)])

# (label, jd0, jd1, n per pair)
REGIONS = [
    ("de441", 1720700.0, SEAM, 60),
    ("de442_early", SEAM, 2414956.5, 30),
    ("dev_1900_2100", 2415100.0, 2488400.0, 80),
    ("de442_late", 2488500.0, 2598000.0, 30),
]


def sparse_de441(data_dir):
    cache = os.path.join(data_dir, "sources", "de441_cache")
    out = os.path.join(data_dir, "sources", "de441_part-1.sparse.bsp")
    if not os.path.isdir(cache):
        return None
    hf = HttpRangeFile(DE441_URL, cache)
    if not os.path.exists(out) or os.path.getmtime(out) < max(
            os.path.getmtime(os.path.join(cache, n)) for n in os.listdir(cache)):
        hf.write_sparse_copy(out)
    return out


def segs(kernel):
    return {(s.target, s.center): s for s in kernel.segments}


def state(seg_map, target, center, jd):
    p = np.zeros(3)
    v = np.zeros(3)
    pt, pc = PATH[target], PATH[center]
    # drop common tail like the reader does (keeps the same arithmetic grouping)
    while pt and pc and pt[-1] == pc[-1]:
        pt, pc = pt[:-1], pc[:-1]
    for sgn, path in ((1.0, pt), (-1.0, pc)):
        for key in path:
            # whole/fraction split keeps jplephem's time arithmetic exact
            w = float(np.floor(jd))
            pos, vel = seg_map[key].compute_and_differentiate(w, jd - w)
            p += sgn * pos
            v += sgn * vel
    return p, v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../..", "data"))
    ap.add_argument("--seed", type=int, default=20261028)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    k442 = JplSPK.open(os.path.join(args.data, "sources", "de442.bsp"))
    s442 = segs(k442)
    sp = sparse_de441(args.data)
    s441 = segs(JplSPK.open(sp)) if sp else None

    rows = []
    for t, c in PAIRS:
        times = []
        for label, a, b, n in REGIONS:
            if label == "de441" and s441 is None:
                continue
            times += [(label, float(x)) for x in rng.uniform(a, b, n)]
        # record boundaries, the seam itself and its neighbourhood
        times += [("boundary", 2451545.0 + 0.5 + 4 * k) for k in range(3)]
        times += [("seam", SEAM), ("seam", SEAM + 1e-6), ("seam", SEAM + 3.7)]
        if s441 is not None:
            times += [("seam", SEAM - 1e-6), ("seam", SEAM - 2.3)]
        for label, jd in times:
            src = s441 if jd < SEAM else s442
            p, v = state(src, t, c, jd)
            if not np.all(np.isfinite(p)) or np.all(p == 0):
                raise SystemExit(f"bad JPL state {t}/{c} at {jd} (sparse DE441 range missing?)")
            # jd also as its IEEE-754 bit pattern: serde_json (without the
            # float_roundtrip feature) may parse decimal floats 1 ulp off, and one
            # JD ulp moves Mercury by ~1e-3 km.
            bits = int(np.float64(jd).view(np.uint64))
            rows.append([t, c, jd, *map(float, p), *map(float, v), label, bits])

    out = os.path.join(HERE, "pack_states.json")
    with open(out, "w") as f:
        json.dump({"generator": "gen_pack_fixtures.py (jplephem " + __import__("jplephem").__version__ + ")",
                   "units": "km, km/day, JD TDB",
                   "columns": ["target", "center", "jd_tdb", "px", "py", "pz", "vx", "vy", "vz", "region",
                               "jd_tdb_f64_bits"],
                   "rows": rows}, f, separators=(",", ":"))
    print(f"wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
