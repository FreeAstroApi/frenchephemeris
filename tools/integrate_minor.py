#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Extend Chiron, Pholus, Eris and Eros before JPL's 1599-12-11 limit by our own
numerical integration ("documented orbit" path; no Swiss data).

    python tools/integrate_minor.py validate Chiron [--h 2]   # forward 1600->2401 vs JPL
    python tools/integrate_minor.py build [Chiron,...]         # backward + clones -> products

Model (crates/fe-core/examples/integrate_minor.rs): massless particle,
barycentric ICRF; point masses Sun, Mercury and Venus barycentres, Earth and
Moon (Eros) or the EMB (others), Mars..Pluto system barycentres from the main
pack (DE441 before 1550, DE442 after); 1PN Sun; GMs from the main pack
manifest (DE442 au^3/d^2 table x AU^3). Gauss-Legendre 6-stage (order 12)
fixed step.

Initial state: the JPL Horizons SPK state (heliocentric, + the main pack's
Sun) at the first JPL epoch, JD 2305426.5 (1599-12-11 TDB), integrated
backwards.  The joint is therefore continuous in position and velocity.

Uncertainty: JPL's orbit covariance (SBDB API, elements e q tp node peri i at
the solution epoch) is mapped to Cartesian (two-body Jacobian at that epoch),
its Cholesky columns are added to the JPL state at the solution epoch, and the
6 clones plus the nominal are integrated back to the target start with the
same model.  Position covariance P(t) = sum_k d_k d_k^T (d_k = clone_k - nominal)
gives 1-sigma axes; reported as 3-sigma km and arcsec (at the geocentric
distance) per epoch. Check: the mapped 3-sigma at 1600-1800 vs Horizons'
own RA/Dec 3-sigma (QUANTITIES 36) at the same epochs.

Products (gitignored): data/sources/minor/int_<number>_<name>.npy  (nominal grid:
t x y z vx vy vz ax ay az, barycentric, km, km/day, km/day^2) and a JSON
summary in docs/MINOR_INTEGRATION.json (committed).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import fit_minor as F  # noqa: E402
from fepk import Fepk  # noqa: E402

AU_KM = 149597870.7
C_KM_DAY = 299792.458 * 86400.0
T0 = F.HZ_START  # 2305426.5, first JPL epoch
T_END = F.MAIN_END
SUMMARY = os.path.join(REPO, "docs", "MINOR_INTEGRATION.json")
BIN = os.path.join(os.environ.get("CARGO_TARGET_DIR", os.path.join(REPO, "target")), "debug", "examples",
                   "integrate_minor")

# body -> (target start JD TDB, step days, use Earth+Moon separately)
# Starts match the Swiss oracle's first dates (black-box probing): Chiron
# 675-01-01, Pholus before year 0 (-> main pack start), Eris/Eros 1500-01-01 (Julian).
TARGETS = {
    "Chiron": (1967601.5, 4.0, False),
    "Pholus": (F.MAIN_START, 4.0, False),
    "Eris": (2268932.5, 8.0, False),
    "Eros": (2268932.5, 0.5, True),
}
NUMBER = {"Chiron": 2060, "Pholus": 5145, "Eris": 136199, "Eros": 433}


def gm_table(main: Fepk) -> dict:
    c = main.manifest["constants"]
    au = c.get("AU_km", AU_KM)
    g = {k: v * au ** 3 for k, v in c["GM_au3_d2"].items()}  # km^3/day^2
    return g


# Big asteroids available in the minor pack (sb441-n16, 0..2400) and their GMs,
# km^3/s^2, DE440 values (Park et al. 2021, AJ 161:105) -- rounded; any residual
# central-mass error is absorbed by the fitted solar GM scale (see fit_state).
ASTEROIDS = [(2000001, 62.6284), (2000002, 13.6659), (2000004, 17.2882), (2000003, 1.82)]
# Only Ceres-Vesta (identical, JPL sb441-n16 fits) are read from it.
MINOR_PACK = next((p for p in (os.path.join(F.DEFAULT_DATA, "packs", f"fepk-minor-v{v}.bin") for v in (2, 1))
                   if os.path.exists(p)), os.path.join(F.DEFAULT_DATA, "packs", "fepk-minor-v2.bin"))


def perturbers(main: Fepk, earth_moon: bool, asteroids: bool = True) -> list:
    g = gm_table(main)
    p = [(10, g["GMS"]), (1, g["GM1"]), (2, g["GM2"])]
    if earth_moon:
        p += [(399, g["GMB"] - g["GMM"]), (301, g["GMM"])]
    else:
        p += [(3, g["GMB"])]
    p += [(k, g[f"GM{k}"]) for k in range(4, 10)]
    if asteroids:
        p += [(t, gm * 86400.0 ** 2) for t, gm in ASTEROIDS]
    return p


def scaled(pert: list, dmu: float) -> list:
    """Solar GM scaled by (1 + dmu)."""
    return [(pert[0][0], pert[0][1] * (1.0 + dmu))] + pert[1:]


def source(name: str) -> F.Source:
    return F.Source(name, NUMBER[name], "horizons", F.DEFAULT_DATA)


def jpl_bary_state(src: F.Source, main: Fepk, jd) -> tuple[np.ndarray, np.ndarray]:
    jd = np.atleast_1d(np.asarray(jd, float))
    p, v = src.spk.state_et(src.target, 10, (jd - F.J2000_JD) * 86400.0)
    sp, sv = main.state(10, 0, jd)
    return p + sp, v * 86400.0 + sv


def run(main_path: str, t0: float, t1: float, h: float, states, pert, every: int, out: str, gr=True,
        extra_packs=(), offsets=()):
    n = max(1, int(round(abs(t1 - t0) / h)))
    h = abs(t1 - t0) / n  # exact: Python and Rust agree on the step count
    cfg = {"packs": [main_path, MINOR_PACK, *extra_packs], "t0": t0, "t1": t1, "h": h, "states": [list(map(float, s)) for s in states],
           "perturbers": [[int(a), float(b)] for a, b in pert], "gm_sun": float(pert[0][1]),
           "c_km_day": C_KM_DAY, "gr": gr, "out": out, "every": every,
           "offsets": [[int(o[0]), float(o[1]), float(o[2]), float(o[3])] for o in offsets]}
    cpath = out + ".json"
    with open(cpath, "w") as f:
        json.dump(cfg, f)
    if not os.path.exists(BIN):
        subprocess.run(["cargo", "build", "-p", "fe-core", "--example", "integrate_minor"], cwd=REPO, check=True)
    # (CARGO_TARGET_DIR is honoured by both cargo and BIN)
    r = subprocess.run([BIN, cpath], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr)
    os.remove(cpath)
    raw = np.fromfile(out, dtype="<f8")
    os.remove(out)
    nom = raw[: (n + 1) * 10].reshape(n + 1, 10)
    rest = raw[(n + 1) * 10:]
    nc = len(states) - 1
    clones = rest.reshape(nc, -1, 7) if nc else None
    return nom, clones, r.stderr.strip()


def scratch(name: str) -> str:
    d = os.path.join(F.DEFAULT_DATA, "sources", "minor")
    return os.path.join(d, f"_run_{name}.bin")


def geo_dist(main: Fepk, jd, pos):
    e = main.pos(3, 0, jd)  # EMB is close enough for the angular scale
    return np.linalg.norm(pos - e, axis=1)


YEAR = 365.25


def jpl_errors(src, main, t, pos):
    jp, _ = jpl_bary_state(src, main, t)
    err = np.linalg.norm(pos - jp, axis=1)
    return err, err / geo_dist(main, t, jp) * 206264.806


def fit_state(name: str, main_path: str, wa: float, wb: float, iters: int = 4):
    """Least-squares fit of the barycentric state at `wa` and a solar GM scale
    (1 + dmu) to the JPL positions over [wa, wb] (samples every ~10 d).
    dmu stands in for the central mass our model lacks (the other 12 SB441
    asteroids, rounded asteroid GMs) and for small model differences.
    Returns (state, dmu, rms_km, max_km)."""
    main = Fepk(main_path)
    src = source(name)
    _, h, em = TARGETS[name]
    pert = perturbers(main, em)
    x, v = jpl_bary_state(src, main, wa)
    p = np.r_[x[0], v[0], 0.0]
    every = max(1, int(round(10.0 / h)))
    eps = np.r_[np.full(3, 10.0), np.full(3, 1e-3)]  # km, km/day
    for it in range(iters):
        st = [p[:6]] + [p[:6] + eps[k] * np.eye(6)[k] for k in range(6)]
        nom, cl, _ = run(main_path, wa, wb, h, st, scaled(pert, p[6]), every, scratch(name))
        t = cl[0][:, 0]
        idx = np.round((t - wa) / ((wb - wa) / (len(nom) - 1))).astype(int)
        base = nom[idx, 1:4]
        dmu_eps = 1e-9
        nd, _, _ = run(main_path, wa, wb, h, [p[:6]], scaled(pert, p[6] + dmu_eps), 1, scratch(name))
        cols = [(c[:, 1:4] - base) / eps[k] for k, c in enumerate(cl)] + [(nd[idx, 1:4] - base) / dmu_eps]
        J = np.stack([c.ravel() for c in cols], axis=1)
        jp, _ = jpl_bary_state(src, main, t)
        r = (jp - base).ravel()
        dp, *_ = np.linalg.lstsq(J, r, rcond=None)
        p = p + dp
        res = np.linalg.norm((jp - base), axis=1)
    nom, _, _ = run(main_path, wa, wb, h, [p[:6]], scaled(pert, p[6]), 1, scratch(name))
    err, ang = jpl_errors(src, main, nom[:, 0][::every], nom[::every, 1:4])
    return p[:6], float(p[6]), float(np.sqrt(np.mean(err ** 2))), float(err.max()), float(ang.max())


def close_approaches(main: Fepk, t, pos, limit_au: float = 1.0) -> list:
    """Minima of the distance to the planet-system barycentres 4..9 below limit_au."""
    out = []
    for k, nm in [(4, "Mars"), (5, "Jupiter"), (6, "Saturn"), (7, "Uranus"), (8, "Neptune"), (9, "Pluto"), (3, "EMB")]:
        d = np.linalg.norm(pos - main.pos(k, 0, t), axis=1) / AU_KM
        i = np.flatnonzero((d[1:-1] < d[:-2]) & (d[1:-1] <= d[2:]) & (d[1:-1] < limit_au)) + 1
        out += [{"planet": nm, "jd_tdb": float(t[j]), "year": float((t[j] - 2451545.0) / YEAR + 2000),
                 "dist_au": float(d[j])} for j in i]
    return sorted(out, key=lambda x: x["jd_tdb"])


def validate(name: str, main_path: str) -> dict:
    """Measure our model against JPL over the JPL span:
    (1) plain: JPL state at 1600, dmu = 0, forward to 2401;
    (2) fit state + solar GM scale dmu over 2000-2401;
    (3) production analog: JPL state at 2401 (exact, as at 1600 in production),
        dmu from (2), integrated backwards 800 years to 1600 (and with dmu = 0)."""
    main = Fepk(main_path)
    src = source(name)
    _, h, em = TARGETS[name]
    pert = perturbers(main, em)
    out = {"step_days": h}
    k = max(1, int(round(10.0 / h)))
    x0, v0 = jpl_bary_state(src, main, T0)
    nom, _, log = run(main_path, T0, T_END, h, [np.r_[x0[0], v0[0]]], pert, 1, scratch(name))
    e, a = jpl_errors(src, main, nom[::k, 0], nom[::k, 1:4])
    out["plain_forward_1600_2401"] = {"max_km": float(e.max()), "max_arcsec": float(a.max()), "log": log}
    jt = np.arange(T0, T_END, 2.0 if name != "Eros" else 0.25)
    jp, _ = jpl_bary_state(src, main, jt)
    out["jpl_close_approaches_lt_1au"] = close_approaches(main, jt, jp)
    mid = 2451545.0
    st, dmu, rms, mx, _ = fit_state(name, main_path, mid, T_END)
    out["fit_2000_2401"] = {"dmu": dmu, "fit_rms_km": rms, "fit_max_km": mx}
    xe, ve = jpl_bary_state(src, main, T_END)
    for label, d in (("dmu_fit", dmu), ("dmu_0", 0.0)):
        back, _, _ = run(main_path, T_END, T0, h, [np.r_[xe[0], ve[0]]], scaled(pert, d), 1, scratch(name))
        e, a = jpl_errors(src, main, back[::k, 0], back[::k, 1:4])
        yrs = (back[::k, 0] - 2451545.0) / YEAR + 2000
        out[f"from_jpl_2401_backwards_{label}"] = {
            str(n): {"max_km": float(e[yrs >= 2401 - n].max()), "max_arcsec": float(a[yrs >= 2401 - n].max())}
            for n in (100, 200, 400, 600, 800)}
    out["dmu"] = dmu
    return out


# ----------------------------------------------------------- covariance
def elements_to_state(el: np.ndarray, gm: float) -> np.ndarray:
    """Heliocentric ecliptic J2000 (IAU76 obliquity) elements e,q(au),tp(JD),node,peri,i (deg)
    at epoch el[6] -> ICRF state km, km/day (two-body)."""
    e, q, tp, node, peri, inc, epoch = el
    q *= AU_KM
    a = q / (1 - e)
    if e < 1:
        n = np.sqrt(gm / a ** 3)
        M = n * (epoch - tp)
        E = M
        for _ in range(50):
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        x = a * (np.cos(E) - e)
        y = a * np.sqrt(1 - e * e) * np.sin(E)
        r = a * (1 - e * np.cos(E))
        vx = -np.sqrt(gm * a) / r * np.sin(E)
        vy = np.sqrt(gm * a * (1 - e * e)) / r * np.cos(E)
    else:
        raise ValueError("hyperbolic")
    O, w, i = np.radians([node, peri, inc])

    def rot(px, py):
        cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
        X = (cO * cw - sO * sw * ci) * px + (-cO * sw - sO * cw * ci) * py
        Y = (sO * cw + cO * sw * ci) * px + (-sO * sw + cO * cw * ci) * py
        Z = (sw * si) * px + (cw * si) * py
        return np.array([X, Y, Z])
    eps = np.radians(84381.448 / 3600.0)
    R = np.array([[1, 0, 0], [0, np.cos(eps), -np.sin(eps)], [0, np.sin(eps), np.cos(eps)]])
    return np.r_[R @ rot(x, y), R @ rot(vx, vy)]


def sbdb(number: int) -> dict:
    d = os.path.join(F.DEFAULT_DATA, "sources", "minor")
    p = os.path.join(d, f"sbdb_{number}.json")
    if not os.path.exists(p):
        url = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={number}&cov=mat&full-prec=1"
        r = subprocess.run(["curl", "-sSf", "--max-time", "120", url], check=True, capture_output=True)
        with open(p, "wb") as f:
            f.write(r.stdout)
    return json.load(open(p))


def clone_states(name: str, main: Fepk, src: F.Source, gm_sun: float):
    """Nominal JPL state + 6 covariance clones at the covariance epoch."""
    d = sbdb(NUMBER[name])
    cov = d["orbit"]["covariance"]
    labels = cov["labels"]
    assert labels == ["e", "q", "tp", "node", "peri", "i"], labels
    C = np.array(cov["data"], dtype=float)
    el = {e["name"]: float(e["value"]) for e in cov["elements"]}
    ep = float(cov["epoch"])
    E0 = np.array([el["e"], el["q"], el["tp"], el["om"], el["w"], el["i"], ep])
    # Jacobian d state / d elements by central differences of the two-body map
    J = np.zeros((6, 6))
    for k in range(6):
        s = np.sqrt(C[k, k]) * 1e-3 or 1e-9
        a, b = E0.copy(), E0.copy()
        a[k] += s
        b[k] -= s
        J[:, k] = (elements_to_state(a, gm_sun) - elements_to_state(b, gm_sun)) / (2 * s)
    Cx = J @ C @ J.T
    L = np.linalg.cholesky(Cx + 1e-30 * np.eye(6))
    x, v = jpl_bary_state(src, main, ep)
    nom = np.r_[x[0], v[0]]
    return ep, nom, [nom] + [nom + L[:, k] for k in range(6)], {
        "orbit_id": d["orbit"]["orbit_id"], "soln_date": d["orbit"]["soln_date"], "epoch": ep,
        "labels": labels, "sigma": [float(np.sqrt(C[k, k])) for k in range(6)],
        "position_sigma_km_at_epoch": float(np.sqrt(np.trace(Cx[:3, :3])))}


def build(name: str, main_path: str, dmu: float) -> dict:
    main = Fepk(main_path)
    src = source(name)
    t_start, h, em = TARGETS[name]
    pert = scaled(perturbers(main, em), dmu)
    gm_sun = pert[0][1]
    # 1) nominal: JPL state at T0, backwards to the target start. Use an exact
    # step h = span / n so that the h/2 check grid nests (2n steps).
    n = int(round(abs(T0 - t_start) / h))
    h = abs(T0 - t_start) / n
    x0, v0 = jpl_bary_state(src, main, T0)
    nom, _, log = run(main_path, T0, t_start, h, [np.r_[x0[0], v0[0]]], pert, 1, scratch(name))
    # 2) step-size check: same with h/2
    nom2, _, _ = run(main_path, T0, t_start, h / 2, [np.r_[x0[0], v0[0]]], pert, 1, scratch(name))
    step_err = np.linalg.norm(nom[:, 1:4] - nom2[::2, 1:4], axis=1)
    # 3) round trip: from the end state back to T0
    back, _, _ = run(main_path, t_start, T0, h, [nom[-1, 1:7]], pert, 1, scratch(name))
    rt = float(np.linalg.norm(back[-1, 1:4] - x0[0]))
    # 4) clones from the covariance epoch to the target start
    ep, cnom, states, covinfo = clone_states(name, main, src, gm_sun)
    every = max(1, int(round(3650.0 / h)))  # ~10-year sampling
    cn, cl, _ = run(main_path, ep, t_start, h, states, pert, every, scratch(name))
    tc = cl[0][:, 0]
    nomc = cn[::every][: len(tc)] if len(cn[::every]) >= len(tc) else None
    # nominal at clone sample times (the last sample is the final step)
    idx = np.round((tc - ep) / ((t_start - ep) / (len(cn) - 1))).astype(int)
    nomc = cn[idx, 1:4]
    d = np.stack([c[:, 1:4] - nomc for c in cl[1:]])  # (6, nt, 3)
    sig3_km, sig3_as = [], []
    egeo = main.pos(3, 0, tc)
    for j in range(len(tc)):
        P = d[:, j, :].T @ d[:, j, :]
        los = nomc[j] - egeo[j]
        dist = np.linalg.norm(los)
        u = los / dist
        Pp = (np.eye(3) - np.outer(u, u)) @ P @ (np.eye(3) - np.outer(u, u))  # transverse part
        sig3_km.append(3 * float(np.sqrt(np.linalg.eigvalsh(P).max())))
        sig3_as.append(3 * float(np.sqrt(max(np.linalg.eigvalsh(Pp).max(), 0.0))) / dist * 206264.806)
    # our nominal from the covariance epoch vs the JPL-joined nominal (model error over the JPL span
    # plus everything before T0)
    jd_years = (tc - 2451545.0) / 365.25 + 2000
    # product
    prod = os.path.join(F.DEFAULT_DATA, "sources", "minor", f"int_{NUMBER[name]}_{name.lower()}.npy")
    np.save(prod, nom[::-1].copy())  # increasing time
    info = {
        "body": name, "number": NUMBER[name], "product": os.path.basename(prod),
        "span_jd_tdb": [t_start, T0], "step_days": h, "steps": len(nom) - 1, "log": log,
        "initial_state": {"jd_tdb": T0, "source": f"Horizons SPK {src.file} ({F.naif_id(NUMBER[name])})",
                          "bary_km_kmday": list(map(float, np.r_[x0[0], v0[0]]))},
        "perturbers": [[a, b] for a, b in pert], "gr": True, "dmu": dmu,
        "close_approaches_lt_3au": close_approaches(main, nom[::max(1, int(4 / h)), 0], nom[::max(1, int(4 / h)), 1:4], 3.0),
        "step_check_max_km": float(step_err.max()),
        "round_trip_km": rt,
        "covariance": covinfo,
        "uncertainty_3sigma": [{"year": float(y), "jd_tdb": float(t), "km": k, "arcsec": a}
                               for y, t, k, a in zip(jd_years, tc, sig3_km, sig3_as)],
    }
    return info


def write_ledger(summary: dict):
    """Record the integrations as sources in docs/MINOR_PROVENANCE.json."""
    path = os.path.join(REPO, "docs", "MINOR_PROVENANCE.json")
    led = json.load(open(path))
    hz = {e["number"]: e for e in led["horizons"]}
    ent = {}
    for name, d in summary["bodies"].items():
        i = d.get("integration")
        if not i:
            continue
        sol = hz[i["number"]]["solution_id"]
        ent[name] = {
            "source": f"FrenchEphemeris integration from {sol}",
            "kind": "own numerical integration (model extrapolation before JPL's first epoch)",
            "product": i["product"], "span_jd_tdb": i["span_jd_tdb"],
            "initial_state": i["initial_state"], "jpl_solution": sol,
            "covariance_from": f"SBDB API orbit {i['covariance']['orbit_id']} ({i['covariance']['soln_date']})",
            "tool": "tools/integrate_minor.py + crates/fe-core/examples/integrate_minor.rs",
            "model": "point masses Sun, planet-system barycentres (Earth+Moon for Eros), Ceres/Pallas/Juno/Vesta; "
                     "1PN Sun; main pack DE441 (<1550) / DE442; Gauss-Legendre 6-stage",
            "step_days": i["step_days"],
            "checks": {"step_halving_max_km": i["step_check_max_km"], "round_trip_km": i["round_trip_km"]},
            "details": "docs/MINOR_INTEGRATION.json",
        }
    led["integrations"] = ent
    with open(path, "w") as f:
        json.dump(led, f, indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["validate", "build", "ledger"])
    ap.add_argument("bodies", nargs="?", default="Chiron,Pholus,Eris,Eros")
    ap.add_argument("--main-pack", default=os.path.join(F.DEFAULT_DATA, "packs", "fepk-de441-de442-v1.bin"))
    args = ap.parse_args()
    names = args.bodies.split(",")
    summary = json.load(open(SUMMARY)) if os.path.exists(SUMMARY) else {"bodies": {}}
    for name in names if args.cmd != "ledger" else []:
        if args.cmd == "validate":
            r = validate(name, args.main_pack)
            print(name, json.dumps(r, indent=1), flush=True)
            summary["bodies"].setdefault(name, {})["validation"] = r
        else:
            # The fitted solar-GM scale (validate) differs in sign between bodies and makes
            # Eros' backward extrapolation worse (0.1" -> 8.8" after 100 years), so it is
            # not physical; production uses the plain model (dmu = 0).
            info = build(name, args.main_pack, 0.0)
            u = info["uncertainty_3sigma"]
            print(f"{name}: steps {info['steps']} step-check {info['step_check_max_km']:.3g} km round-trip "
                  f"{info['round_trip_km']:.3g} km; 3sigma at start {u[-1]['km']:.3g} km {u[-1]['arcsec']:.3g}\"",
                  flush=True)
            summary["bodies"].setdefault(name, {})["integration"] = info
    if args.cmd in ("build", "ledger"):
        write_ledger(summary)
    summary["generated_utc"] = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    summary["tool"] = "tools/integrate_minor.py + crates/fe-core/examples/integrate_minor.rs"
    with open(SUMMARY, "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
