#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Why does our Chiron differ from the Swiss oracle before ~1000 AD?

    CARGO_TARGET_DIR=... python tools/chiron_study.py   -> docs/CHIRON_STUDY.json

Controlled experiments with our integrator (tools/integrate_minor.py); every
variant is compared with the production trajectory (JPL#171 state at
1599-12-11, full model) as a geocentric angle (observer at the EMB; geometric
directions; arcsec):

(b) dynamical model, same start:
    no_asteroids  - Ceres/Pallas/Juno/Vesta removed (Swiss's docs: it includes
                    Ceres, Pallas, Vesta);
    no_gr         - Sun's 1PN term removed;
    de200_masses  - planetary masses of DE200 (Standish 1990; the docs say the
                    Swiss asteroid integrator was built to reproduce DE200);
    saturn_1as    - Saturn barycentre displaced by 6,500 km (~1" at 9 AU), a
                    proxy for older planetary ephemerides near the 720 AD
                    encounter; uranus_2200km likewise (the DE441/DE442 size).
(a) initial orbit, same model, integrated from the solution epoch 2026 (JD
    2461200.5) instead of 1599:
    jpl_2026      - JPL#171 state from the Horizons SPK at that epoch;
    mpc_2026      - the MPC's independent orbit (orbfit, DE431, 17 asteroids;
                    data.minorplanetcenter.net get-orb, Cartesian ecliptic J2000).
(c) chaos: the JPL covariance clones of docs/MINOR_INTEGRATION.json (3-sigma,
    scaled x4 to Horizons' Chiron 3-sigma).

The FE - Swiss column is the black-box comparison made by
harness/check_minor.py (section e) at the same 10-year samples. Nothing here
is fitted to Swiss output; Swiss only appears as the column being explained.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import integrate_minor as I  # noqa: E402
from privatize_oracle_fixtures import oracle_file, write_split  # noqa: E402,F401  (private Swiss-oracle corpus)
from fepk import Fepk  # noqa: E402

MAIN = os.path.join(I.F.DEFAULT_DATA, "packs", "fepk-de441-de442-v1.bin")
OUT = os.path.join(REPO, "docs", "CHIRON_STUDY.json")
EPOCH_2026 = 2461200.5
T_START = 1967601.5  # 675-01-01
H = 4.0
YEARS = [(2305426.5 - 2451545.0) / 365.25 + 2000, 1500, 1400, 1300, 1200, 1100, 1000, 950, 900, 850, 800, 760, 740, 730, 725, 721, 715, 710,
         700, 690, 680, 676]

# DE200 planetary mass ratios Sun/planet (Standish 1990, A&A 233, 252; also the IAU 1976 values
# for the outer planets used by DE200). Pluto as a system.
DE200_RATIO = {1: 6023600.0, 2: 408523.5, 3: 328900.55, 4: 3098710.0, 5: 1047.350, 6: 3498.0,
               7: 22960.0, 8: 19314.0, 9: 130000000.0}


def mpc_state(main: Fepk) -> np.ndarray:
    p = os.path.join(I.F.DEFAULT_DATA, "sources", "minor", "mpc_2060.json")
    if not os.path.exists(p):
        r = subprocess.run(["curl", "-sSf", "--max-time", "120", "-X", "GET",
                            "https://data.minorplanetcenter.net/api/get-orb",
                            "-H", "Content-Type: application/json", "-d", '{"desig":"2060"}'],
                           check=True, capture_output=True)
        open(p, "wb").write(r.stdout)
    o = json.load(open(p))[0]["mpc_orb"][0]
    assert o["epoch_data"]["epoch"] == 61200.0 and o["system_data"]["refsys"] == "Ecliptic"
    x = np.array(o["CAR"]["coefficient_values"])  # au, au/day, heliocentric ecliptic J2000
    eps = np.radians(float(o["system_data"]["EclipticObliquityArcseconds"]) / 3600.0)
    R = np.array([[1, 0, 0], [0, np.cos(eps), -np.sin(eps)], [0, np.sin(eps), np.cos(eps)]])
    au = I.AU_KM
    sp, sv = main.state(10, 0, [EPOCH_2026])
    meta = {k: o[k] for k in ("epoch_data", "system_data", "software_data")}
    meta["orbit_fit_statistics"] = {k: o["orbit_fit_statistics"][k] for k in ("arc_length_total", "nobs_total_sel",
                                                                             "normalized_RMS")}
    return np.r_[R @ x[:3] * au + sp[0], R @ x[3:] * au + sv[0]], meta


def angles(main: Fepk, a: np.ndarray, b: np.ndarray, t: np.ndarray) -> np.ndarray:
    e = main.pos(3, 0, t)
    u, v = a - e, b - e
    c = np.einsum("ij,ij->i", u, v) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600.0


def at_years(traj: np.ndarray, years) -> np.ndarray:
    """Rows (t, x, y, z) at the exact epochs of `years`, by quintic Hermite
    interpolation of the integration grid (t, x, v, a)."""
    g = traj[np.argsort(traj[:, 0])]
    jd = np.array([2451545.0 + (y - 2000) * 365.25 for y in years])
    jd = np.clip(jd, g[0, 0], g[-1, 0])
    obj = I.F.IntegratedSource.__new__(I.F.IntegratedSource)
    obj.g = g
    return np.c_[jd, obj.bary(jd)]


def main():
    main_pk = Fepk(MAIN)
    src = I.source("Chiron")
    pert = I.perturbers(main_pk, False)
    x0, v0 = I.jpl_bary_state(src, main_pk, I.T0)
    s1599 = np.r_[x0[0], v0[0]]
    tmp = I.scratch("chiron_study")

    def back(state, t0, p=pert, gr=True, offsets=()):
        nom, _, _ = I.run(MAIN, t0, T_START, H, [state], p, 1, tmp, gr, offsets=offsets)
        return nom

    runs = {}
    base = back(s1599, I.T0)
    runs["no_asteroids"] = back(s1599, I.T0, p=[q for q in pert if q[0] < 1000000])
    runs["no_gr"] = back(s1599, I.T0, gr=False)
    gms = pert[0][1]
    runs["de200_masses"] = back(s1599, I.T0, p=[(t, gms / DE200_RATIO[t]) if t in DE200_RATIO else (t, g)
                                                for t, g in pert])
    # displace along Saturn's orbital velocity at 720 AD (an along-track ephemeris error)
    jd720 = 2451545.0 + (720.7 - 2000) * 365.25
    _, sv = main_pk.state(6, 0, [jd720])
    u = sv[0] / np.linalg.norm(sv[0])
    runs["saturn_1as"] = back(s1599, I.T0, offsets=[(6, *(6500.0 * u))])
    _, uv = main_pk.state(7, 0, [jd720])
    u7 = uv[0] / np.linalg.norm(uv[0])
    runs["uranus_2200km"] = back(s1599, I.T0, offsets=[(7, *(2200.0 * u7))])
    runs["h_2d"] = I.run(MAIN, I.T0, T_START, H / 2, [s1599], pert, 1, tmp)[0]
    # initial-orbit variants from the 2026 solution epoch
    xj, vj = I.jpl_bary_state(src, main_pk, EPOCH_2026)
    jpl26 = back(np.r_[xj[0], vj[0]], EPOCH_2026)
    mpc, mpc_meta = mpc_state(main_pk)
    mpc26 = back(mpc, EPOCH_2026)
    start_diff_km = float(np.linalg.norm(mpc[:3] - np.r_[xj[0]]))

    b = at_years(base, YEARS)
    t = b[:, 0]
    out = {"years": YEARS, "jd_tdb": t.tolist(), "variants_vs_production_arcsec": {}}
    for k, r in runs.items():
        out["variants_vs_production_arcsec"][k] = angles(main_pk, at_years(r, YEARS)[:, 1:4], b[:, 1:4], t).tolist()
    j = at_years(jpl26, YEARS)
    m = at_years(mpc26, YEARS)
    out["variants_vs_production_arcsec"]["jpl_2026_our_model"] = angles(main_pk, j[:, 1:4], b[:, 1:4], t).tolist()
    out["mpc_minus_jpl_both_from_2026_arcsec"] = angles(main_pk, m[:, 1:4], j[:, 1:4], t).tolist()
    out["mpc_vs_jpl_start_km"] = start_diff_km
    out["mpc_orbit"] = mpc_meta
    # chaos: clone 3-sigma (x4) from the production integration summary
    integ = json.load(open(os.path.join(REPO, "docs", "MINOR_INTEGRATION.json")))["bodies"]["Chiron"]["integration"]
    u3 = integ["uncertainty_3sigma"]
    out["jpl_3sigma_x4_arcsec"] = [4 * min(u3, key=lambda r: abs(r["year"] - y))["arcsec"] for y in YEARS]
    # FE - Swiss from the harness (black box), nearest 10-year sample
    rep = json.load(open(oracle_file(os.path.join(REPO, "harness", "minor_report.json"))))  # full report (private)
    rows = rep["e_integrated_span_vs_swiss"]["Chiron"]
    out["fe_minus_swiss_arcsec"] = [min(rows, key=lambda r: abs(r["year"] - y))["fe_minus_swiss_arcsec"]
                                    if min(abs(r["year"] - y) for r in rows) < 6 else None for y in YEARS]
    # close approaches to Saturn in the production run
    out["saturn_encounters_lt_1au"] = [c for c in I.close_approaches(main_pk, base[::2, 0], base[::2, 1:4], 1.0)]
    write_split(OUT, out, indent=1)  # public copy without the Swiss column
    cols = ["fe_minus_swiss", "3sig_x4", "mpc-jpl"] + list(out["variants_vs_production_arcsec"])
    print("year " + " ".join(f"{c[:12]:>12}" for c in cols))
    for i, y in enumerate(YEARS):
        vals = [out["fe_minus_swiss_arcsec"][i], out["jpl_3sigma_x4_arcsec"][i], out["mpc_minus_jpl_both_from_2026_arcsec"][i]]
        vals += [out["variants_vs_production_arcsec"][k][i] for k in out["variants_vs_production_arcsec"]]
        print(f"{y:>7} " + " ".join(f"{v:12.3f}" if v is not None else f"{'-':>12}" for v in vals))
    print("saturn encounters:", out["saturn_encounters_lt_1au"])


if __name__ == "__main__":
    main()
