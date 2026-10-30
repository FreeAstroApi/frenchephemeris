#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""DE441 -> DE442 seam report for an FEPK pack.

At the seam epoch (DE442 start, a shared record boundary) compare:
* geocentric geometric directions (body wrt SSB minus Earth wrt SSB, no light
  time) from DE441 vs DE442, in arcseconds, plus the geocentric distance change;
* the jump of each stored series (position km, velocity km/day) between the
  last DE441 record evaluated at its end and the first DE442 record at its start,
  read from the pack itself.

Prints a Markdown section (used in docs/DATA_REPORT.md).
"""
import argparse
import os
import sys

import numpy as np
from jplephem.spk import SPK as JplSPK

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_pack_vs_jpl import DATA, FepkPy, jpl_state, sparse_de441  # noqa: E402

NAMES = {1: "Mercury (bary)", 2: "Venus (bary)", 3: "EMB", 4: "Mars (bary)", 5: "Jupiter (bary)",
         6: "Saturn (bary)", 7: "Uranus (bary)", 8: "Neptune (bary)", 9: "Pluto (bary)", 10: "Sun",
         301: "Moon", 399: "Earth"}
ARCSEC = np.degrees(1.0) * 3600.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default=os.path.join(DATA, "packs", "fepk-de441-de442-v1.bin"))
    args = ap.parse_args()
    pk = FepkPy(args.pack)
    seam = min(s["start"] for s in pk.series if pk.manifest["sources"][s["source_index"]]["name"] == "de442.bsp")
    k = {"de441": {(s.target, s.center): s for s in JplSPK.open(sparse_de441()).segments},
         "de442": {(s.target, s.center): s for s in JplSPK.open(os.path.join(DATA, "sources", "de442.bsp")).segments}}
    jd = np.array([seam])
    print(f"### Geocentric directions at the seam, JD {seam} TDB (DE441 vs DE442)\n")
    print("| body | angle DE441 vs DE442 (arcsec) | geocentric distance change (km) |")
    print("|---|---:|---:|")
    for b in (10, 301, 1, 2, 4, 5, 6, 7, 8, 9):
        r1, _ = jpl_state(k, np.inf, b, 399, jd)      # seam=inf -> DE441 everywhere
        r2, _ = jpl_state(k, -np.inf, b, 399, jd)     # seam=-inf -> DE442 everywhere
        r1, r2 = r1[0], r2[0]
        ang = np.arctan2(np.linalg.norm(np.cross(r1, r2)), np.dot(r1, r2)) * ARCSEC
        print(f"| {NAMES[b]} | {ang:.6f} | {np.linalg.norm(r2) - np.linalg.norm(r1):+.3f} |")

    print(f"\n### Jump of each stored series across the seam (pack: DE441 at end of its last record "
          f"vs DE442 at start of its first)\n")
    print("| series | position jump (km) | velocity jump (km/day) |")
    print("|---|---:|---:|")
    before = np.array([np.nextafter(seam, -np.inf)])
    for t, c in sorted({(s["target"], s["center"]) for s in pk.series}):
        # evaluate DE441's last record exactly at its end: FepkPy picks DE442 at t == seam,
        # so evaluate the DE441 series directly at u = +1
        s441 = [s for s in pk.pair_series(t, c) if s["end"] == seam][0]
        C = s441["coeffs"][-1]
        p1 = C.sum(axis=1)  # T_k(1) = 1
        n = s441["n_coeff"]
        v1 = (C * (np.arange(n) ** 2)).sum(axis=1) * 2.0 / s441["interval"]  # T_k'(1) = k^2
        p2, v2 = pk.state(t, c, jd)
        pb, _ = pk.state(t, c, before)
        dp = np.linalg.norm(p2[0] - p1)
        dv = np.linalg.norm(v2[0] - v1)
        assert np.linalg.norm(pb[0] - p1) < 1e-3 + dp, "pack before-seam value inconsistent"
        print(f"| {t} wrt {c} ({NAMES[t]}) | {dp:.3f} | {dv:.3e} |")


if __name__ == "__main__":
    main()
