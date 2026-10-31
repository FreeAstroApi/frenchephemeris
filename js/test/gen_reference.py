# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Reference values for js/test/engine.test.js, from the Python package with the same mobile
packs: the npm build must give the same numbers. Run: .venv/bin/python -I js/test/gen_reference.py"""
import json, sys
import frenchephemeris as fe
from frenchephemeris import _native
from pathlib import Path
P = str(Path(__file__).resolve().parents[2] / "data" / "packs") + "/"
eng = fe.Engine(P + "fepk-de442-mobile-1800-2400-v1.bin", True, [P + "fepk-minor-mobile-1800-2400-v1.bin"])
jd = _native.julday(1990, 5, 17, 14.5, True)
ids = {"sun": 0, "moon": 1, "mercury": 2, "mars": 4, "pluto": 9, "trueNode": 11, "chiron": 15, "ceres": 17}
out = {"jd": jd,
       "bodies": {n: list(eng.calc_ut(jd, i, 256)) for n, i in ids.items()},
       "sidereal": {n: list(eng.calc_ut(jd, i, 256 | 65536, "lahiri")) for n, i in ids.items()},
       "equatorial": list(eng.calc_ut(jd, 4, 256 | 2048)),
       "topo_moon": list(eng.calc_ut(jd, 1, 256 | 32768, None, (2.3522, 48.8566, 35.0))),
       "ayanamsa": eng.ayanamsa_ut(jd, "lahiri", False),
       "deltat_s": eng.deltat(jd) * 86400}
c, a = _native.houses_ex2(eng, jd, 48.8566, 2.3522, "P")[:2]
out["houses"] = {"cusps": list(c)[:12], "ascmc": list(a)[:8]}
out["spica"] = list(_native.fixstar_ut(eng, "Spica", jd, 256)[0])
from frenchephemeris.ext import space
s = space.Site("mars", 77.45, 18.44)
out["mars_earth"] = space.position(eng, jd, "earth", s)
out["mars_houses"] = space.houses(eng, jd, s, "P")
json.dump(out, open(Path(__file__).with_name("reference.json"), "w"))
print("ok")
