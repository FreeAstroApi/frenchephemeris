# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Fixtures for fe-core/tests/stars.rs: space motion and astrometric->apparent
chain computed with ERFA (pyerfa), independently of the Rust code.

(a) erfa.pmsafe from the catalogue epoch to dates spanning 0..2400 CE;
(b) rigorous GCRS apparent direction: erfa.starpv position at the date,
    minus the observer (erfa.epv00 Earth), erfa.ldsun, erfa.ab.

Run: .venv/bin/python -I crates/fe-core/tests/fixtures/stars/gen_erfa.py
"""
import json
import math
import sys
from pathlib import Path

import erfa
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from _catalog import load  # noqa: E402

MAS = math.pi / 180 / 3600 / 1000
C_AU_DAY = 299792.458 * 86400 / 149597870.7

epoch, stars, _ = load()
# Every 8th star plus the high proper-motion / nearby / ayanamsa stars.
keep = {"Spica", "Asellus Australis", "Rigil Kentaurus", "Toliman", "Sirius", "Procyon", "Arcturus",
        "Aldebaran", "Regulus", "Antares", "Polaris", "Achird", "Barnard's Star", "Deneb", "Rigel"}
sel = [dict(s, idx=i) for i, s in enumerate(stars) if i % 8 == 0 or s["name"] in keep]
tdbs = [1721424.0, 1830000.5, 2000000.5, 2299160.5, 2415020.0, 2451545.0, 2460000.5, 2524594.0, 2597641.5]

rows = []
for s in sel:
    ra, dec = s["ra"], s["dec"]
    pmr = s["pmra"] * MAS / math.cos(dec)
    pmd = s["pmde"] * MAS
    px = s["plx"] / 1000.0
    for t in tdbs:
        ra2, dec2, pmr2, pmd2, px2, rv2 = erfa.pmsafe(ra, dec, pmr, pmd, px, s["rv"], epoch, 0.0, t, 0.0)
        pv = erfa.starpv(ra2, dec2, pmr2, pmd2, px2, rv2)
        p = np.array(pv["p"], dtype=float)
        pvh, pvb = erfa.epv00(t, 0.0)
        eh = np.array(pvh["p"], dtype=float)
        pob, vob = np.array(pvb["p"], dtype=float), np.array(pvb["v"], dtype=float)
        sun = pob - eh
        d = p - pob
        dist = float(np.linalg.norm(d))
        pn = d / dist
        em = float(np.linalg.norm(eh))
        p1 = erfa.ldsun(pn, eh / em, em)
        v = vob / C_AU_DAY
        bm1 = math.sqrt(1 - float(v @ v))
        p2 = erfa.ab(p1, v, em, bm1)
        rows.append({
            "idx": s["idx"], "tdb": t,
            "pm": [float(ra2), float(dec2), float(pmr2), float(pmd2), float(px2), float(rv2)],
            "obs": [*map(float, pob), *map(float, vob), *map(float, sun)],
            "astrometric": [*map(float, pn)], "apparent": [*map(float, p2)], "dist_au": dist,
        })

out = Path(__file__).with_name("erfa.json")
out.write_text(json.dumps({"catalog_epoch": epoch, "rows": rows}, separators=(",", ":")))
print(f"wrote {out.name}: {len(rows)} rows ({len(sel)} stars x {len(tdbs)} dates)")
