#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""ERFA (pyerfa) frames at extreme epochs, −13200 … +17191, for tests/frames_ext.rs:
Vondrák 2011 long-term precession (ltpb, valid ±200 kyr), IAU 2000A nutation,
Earth rotation angle, IAU 2006 obliquity polynomial (for reference only)."""
import json, os
import numpy as np
import erfa

J2000 = 2451545.0
rng = np.random.default_rng(441)
years = [-13200, -13000, -10000, -7000, -4000, -1000, 3000, 6000, 9000, 12000, 15000, 17191]
jds = [J2000 + (y - 2000) * 365.25 + float(rng.uniform(0, 1)) for y in years]
jds += [J2000 + float(u) * 365.25 for u in rng.uniform(-15200, 15190, 20)]
rows = []
for jd in jds:
    rows.append(dict(jd_bits=int(np.float64(jd).view(np.uint64)), ltpb=erfa.ltpb(2000.0 + (jd - J2000) / 365.25).tolist(),
                     ltp=erfa.ltp(2000.0 + (jd - J2000) / 365.25).tolist(), nut00a=list(erfa.nut00a(jd, 0.0)), era00=erfa.era00(jd, 0.0)))
json.dump({"generator": "gen_frames_ext.py (pyerfa " + erfa.__version__ + ")", "rows": rows},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "frames_ext.json"), "w"))
print(len(rows))
