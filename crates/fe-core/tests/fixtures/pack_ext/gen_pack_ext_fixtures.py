#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Fixtures for tests/pack_ext.rs: JPL Horizons (DE441) geometric geocentric ICRF
vectors (km, TDB) at random epochs in −9998 … −1 and 2401 … 9999 (Horizons serves
DE441 only inside 9999 BC … AD 9999), for the Sun, Moon, planet barycentres and
the EMB. Uses the cached fetcher of harness/horizons_ext.py.

    .venv/bin/python crates/fe-core/tests/fixtures/pack_ext/gen_pack_ext_fixtures.py
"""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "harness"))
import numpy as np  # noqa: E402
from horizons_ext import ERAS, TARGETS, fetch  # noqa: E402

rng = random.Random(20261029)
rows = []
for name, naif, cmd in TARGETS:
    for era, (lo, hi) in ERAS.items():
        jds = sorted(round(rng.uniform(lo, hi), 6) for _ in range(12))
        h = fetch(cmd, jds)
        for jd, v in zip(jds, h):
            rows.append([naif, 399, jd, *map(float, v), int(np.float64(jd).view(np.uint64)), era])
json.dump({"generator": "gen_pack_ext_fixtures.py (JPL Horizons API, DE441, VEC_CORR NONE)",
           "columns": ["target", "center", "jd_tdb", "x_km", "y_km", "z_km", "jd_bits", "era"],
           "rows": rows}, open(os.path.join(HERE, "horizons_vectors.json"), "w"), separators=(",", ":"))
print(len(rows), "rows")
