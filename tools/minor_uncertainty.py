#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Position uncertainty tables for the integrated minor bodies (Chiron, Pholus, Eris, Eros).

    python tools/minor_uncertainty.py   -> python/frenchephemeris/minor_uncertainty.json

The backward half (covariance epoch -> start of each body's range) is the clone propagation
already recorded by `tools/integrate_minor.py build` in docs/MINOR_INTEGRATION.json. This tool
adds the forward half: the same six covariance clones of JPL's orbit solution (SBDB API),
integrated with the same model from the covariance epoch to 2401. Each row is the largest
1-sigma axis of the position covariance, in km; `frenchephemeris.minor_uncertainty` turns it
into arcseconds at the body's geocentric distance.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import integrate_minor as im  # noqa: E402
from fepk import Fepk  # noqa: E402

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(REPO, "python", "frenchephemeris", "minor_uncertainty.json")
END_JD = 2598000.5  # 2401
# Scale applied so the table is never below JPL Horizons' own RA/Dec 3-sigma (QUANTITIES 36,
# queried 2026-10-30 at 1700, 1900, 2000, 2100, 2300 and 2400): the clone propagation gave
# 0.27-0.88 x Horizons for Chiron and 0.88-2.6 x for Pholus; Eris and Eros were at or above it.
SCALE = {"Chiron": 3.7, "Pholus": 1.14}
STEP_DAYS = 4.0


def forward(name: str, main_path: str) -> list[tuple[float, float]]:
    main = Fepk(main_path)
    src = im.source(name)
    _t_start, _h, em = im.TARGETS[name]
    pert = im.perturbers(main, em)
    ep, _nom, states, _info = im.clone_states(name, main, src, pert[0][1])
    n = int(round((END_JD - ep) / STEP_DAYS))
    h = (END_JD - ep) / n
    every = max(1, int(round(3650.0 / h)))
    cn, cl, _ = im.run(main_path, ep, END_JD, h, states, pert, every, im.scratch(name))
    tc = cl[0][:, 0]
    idx = np.round((tc - ep) / ((END_JD - ep) / (len(cn) - 1))).astype(int)
    nomc = cn[idx, 1:4]
    d = np.stack([c[:, 1:4] - nomc for c in cl[1:]])
    return [(float(t), float(np.sqrt(np.linalg.eigvalsh(d[:, j, :].T @ d[:, j, :]).max())))
            for j, t in enumerate(tc)]


def main() -> int:
    main_path = os.path.join(im.F.DEFAULT_DATA, "packs", "fepk-de441-de442-v1.bin")
    integ = json.load(open(os.path.join(REPO, "docs", "MINOR_INTEGRATION.json")))["bodies"]
    out = {"description": "1-sigma position uncertainty (largest axis, km) from JPL's orbit covariance",
           "tool": "tools/minor_uncertainty.py", "bodies": {}}
    for name, number in im.NUMBER.items():
        back = [(r["jd_tdb"], r["km"] / 3.0) for r in integ[name]["integration"]["uncertainty_3sigma"]]
        fwd = forward(name, main_path)
        k = SCALE.get(name, 1.0)
        rows = sorted({round(t, 3): s * k for t, s in back + fwd}.items())
        cov = integ[name]["integration"]["covariance"]
        out["bodies"][str(number)] = {"name": name, "orbit_id": cov["orbit_id"], "soln_date": cov["soln_date"],
                                      "scale_to_horizons": k,
                                      "jd_tdb": [t for t, _ in rows], "sigma_km": [s for _, s in rows]}
        print(f"{name}: {len(rows)} rows, 1-sigma {min(s for _, s in rows):.3g} .. {max(s for _, s in rows):.3g} km")
    with open(OUT, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print(f"wrote {OUT} ({os.path.getsize(OUT)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
