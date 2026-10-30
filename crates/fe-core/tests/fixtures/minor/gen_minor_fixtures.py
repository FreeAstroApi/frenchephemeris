# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Fixtures for crates/fe-core/tests/minor.rs: heliocentric ICRF positions (km)
of the eight minor bodies straight from the JPL source SPKs (tools/spk.py:
sb441-n16 type 2, Horizons type 21), at random epochs inside each body's
supported range plus both range ends.

    .venv/bin/python crates/fe-core/tests/fixtures/minor/gen_minor_fixtures.py

Rows: [naif_id, name, jd_tdb, x, y, z, jd_tdb_bits, source kind]. Also records each body's
coverage (from the pack builder's constants) for the range tests.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../../../tools"))
import fit_minor as F  # noqa: E402

from fepk import Fepk  # noqa: E402

rng = np.random.default_rng(2060)
rows, coverage = [], {}
F.Source.main = F.IntegratedSource.main = Fepk(os.path.join(F.DEFAULT_DATA, "packs", "fepk-de441-de442-v1.bin"))
for name, number, api_id, kind in F.BODIES:
    srcs = [F.Source(name, number, kind, F.DEFAULT_DATA)]
    if kind == "horizons":
        # our own integration before JPL's first epoch (tools/integrate_minor.py)
        srcs.insert(0, F.IntegratedSource(name, number, F.DEFAULT_DATA))
    for src in srcs:
        jds = np.concatenate([[src.start, src.end], np.round(rng.uniform(src.start, src.end, 30), 6)])
        p = src.helio(jds)
        for jd, v in zip(jds, p):
            rows.append([F.naif_id(number), name, float(jd), *map(float, v), int(np.float64(jd).view(np.uint64)),
                         src.kind])
    coverage[name] = {"naif": F.naif_id(number), "api_id": api_id, "jd_tdb": [srcs[0].start, srcs[-1].end],
                      "jpl_jd_tdb": [srcs[-1].start, srcs[-1].end]}
out = {"generator": "gen_minor_fixtures.py", "units": "km, ICRF, wrt Sun (10), JD TDB",
       "source": ("JPL sb441-n16.bsp / Horizons SPK, and before 1599-12-11 our integration from the JPL "
                  "state (row[7] = source kind; see docs/MINOR_PROVENANCE.json)"),
       "coverage": coverage, "rows": rows}
with open(os.path.join(HERE, "minor_states.json"), "w") as f:
    json.dump(out, f, indent=0)
print(len(rows), "rows")
