# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Emit the truncated ELP 2000-82B series used by the natural (interpolated)
lunar apogee and perigee (crates/fe-core/src/points/elp_main.rs).

Source (no Swiss Ephemeris numbers are used):
  * Theory: ELP 2000-82B, M. Chapront-Touzé & J. Chapront, "The lunar
    ephemeris ELP 2000", A&A 124, 50 (1983) and "ELP 2000-85: a
    semi-analytical lunar ephemeris adequate for historical times", A&A 190,
    342 (1988); constants fitted to DE200/LE200 (notice "Lunar solution
    ELP 2000-82B", G. Francou, Bureau des Longitudes, Feb 1996).
  * Files: the authors' 36 series files ELP1..ELP36 as distributed by IMCCE -
    Observatoire de Paris, ftp://ftp.imcce.fr/pub/ephem/moon/elp82b/ (also
    CDS catalogue VI/79). They are downloaded into a cache and NOT shipped;
    the SHA-256 of each file is printed and must match the table in
    docs/FICTITIOUS_BODIES.md.
  * Evaluation conventions (arguments, constant corrections, units, linear
    arguments of the perturbation series, the a0/ath distance factor) follow
    the authors' Fortran routine ELP82B_2 (file elp82b_2, same directory).
  * Truncation: |A|·20^k > CUT rad (distance: |A|/384747.98 km), k = power of T.

    python tools/gen_intp_series.py [--cache data/sources/elp82b] [--cut 1e-7] \
        --out crates/fe-core/src/points/elp_main.rs
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess

import numpy as np

URL = "ftp://ftp.imcce.fr/pub/ephem/moon/elp82b/"
FILES = [f"ELP{i}" for i in range(1, 37)]

RAD = 648000.0 / np.pi  # arcsec per radian
DEG = np.pi / 180.0
ATH = 384747.9806743165
A0 = 384747.9806448954


def dms(d, m, s):
    return (d + m / 60.0 + s / 3600.0) * DEG


# ---- ELP 2000-82B constants (Fortran ELP82B_2) -------------------------------
AM = 0.074801329518
ALPHA = 0.002571881335
DTASM = 2.0 * ALPHA / (3.0 * AM)
W = [  # W1, W2, W3: constant (rad), T1..T4 (rad)
    (dms(218, 18, 59.95571), [1732559343.73604 / RAD, -5.8883 / RAD, 0.6604e-2 / RAD, -0.3169e-4 / RAD]),
    (dms(83, 21, 11.67475), [14643420.2632 / RAD, -38.2776 / RAD, -0.45047e-1 / RAD, 0.21301e-3 / RAD]),
    (dms(125, 2, 40.39816), [-6967919.3622 / RAD, 6.3622 / RAD, 0.7625e-2 / RAD, -0.3586e-4 / RAD]),
]
EART = (dms(100, 27, 59.22059), [129597742.2758 / RAD, -0.0202 / RAD, 0.9e-5 / RAD, 0.15e-6 / RAD])
PERI = (dms(102, 56, 14.42753), [1161.2283 / RAD, 0.5327 / RAD, -0.138e-3 / RAD, 0.0])
PRECESS = 5029.0966 / RAD
PLANETS = [  # Me, V, EMB (= eart), Ma, J, S, U, N: constant, T1
    (dms(252, 15, 3.25986), 538101628.68898 / RAD),
    (dms(181, 58, 47.28305), 210664136.43355 / RAD),
    (EART[0], EART[1][0]),
    (dms(355, 25, 59.78866), 68905077.59284 / RAD),
    (dms(34, 21, 5.34212), 10925660.42861 / RAD),
    (dms(50, 4, 38.89694), 4399609.65932 / RAD),
    (dms(314, 3, 18.01841), 1542481.19393 / RAD),
    (dms(304, 20, 55.19575), 786550.32074 / RAD),
]
W11 = W[0][1][0]
DELNU = 0.55604 / RAD / W11
DELE = 0.01789 / RAD
DELG = -0.08066 / RAD
DELNP = -0.06424 / RAD / W11
DELEP = -0.12879 / RAD


def comb(*terms):
    return (sum(s * t[0] for s, t in terms), [sum(s * t[1][k] for s, t in terms) for k in range(4)])


# Delaunay arguments; the authors' order del(1..4) is D, l', l, F
D = comb((1, W[0]), (-1, EART))
D = (D[0] + np.pi, D[1])
LP = comb((1, EART), (-1, PERI))
L = comb((1, W[0]), (-1, W[1]))
F = comb((1, W[0]), (-1, W[2]))
# Output argument order: D, F, l, l', Me, V, EMB, Ma, J, S, U, N, zeta
ARGS = [D, F, L, LP] + [(c0, [c1, 0.0, 0.0, 0.0]) for c0, c1 in PLANETS] + [(W[0][0], [W11 + PRECESS, 0.0, 0.0, 0.0])]
DEL_IDX = [0, 3, 2, 1]  # del(1..4) = D, l', l, F -> output index


def fetch(cache):
    os.makedirs(cache, exist_ok=True)
    sums = {}
    for f in FILES:
        p = os.path.join(cache, f)
        if not os.path.exists(p):
            subprocess.run(["curl", "-sf", "-o", p, URL + f], check=True)
        sums[f] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return sums


def read(cache, f):
    with open(os.path.join(cache, f)) as fh:
        return fh.read().splitlines()[1:]


def parse(cache):
    """Terms per component: (mult[13], amp (rad or km), phase rad, T power, main?)."""
    out = {1: [], 2: [], 3: []}
    for ific in range(1, 37):
        iv = (ific - 1) % 3 + 1
        itab = (ific + 2) // 3
        k = 1 if itab in (3, 5, 7, 9) else 2 if itab == 12 else 0
        for line in read(cache, f"ELP{ific}"):
            if not line.strip():
                continue
            m = [0] * 13
            if ific <= 3:  # main problem: (4i3,2x,f13.5,6(2x,f10.2))
                ilu = [int(line[3 * i:3 * i + 3]) for i in range(4)]
                coef = [float(line[14:27])] + [float(line[29 + 12 * j:39 + 12 * j]) for j in range(6)]
                if ific == 3:
                    coef[0] = coef[0] - 2.0 * coef[0] * DELNU / 3.0
                tgv = coef[1] + DTASM * coef[5]
                x = coef[0] + tgv * (DELNP - AM * DELNU) + coef[2] * DELG + coef[3] * DELE + coef[4] * DELEP
                for i in range(4):
                    m[DEL_IDX[i]] += ilu[i]
                ph = np.pi / 2 if iv == 3 else 0.0
                main = True
            elif ific <= 9 or ific >= 22:  # (5i3,1x,f9.5,1x,f9.5): iz, ilu(4), phase deg, amp
                iz = int(line[0:3])
                ilu = [int(line[3 + 3 * i:6 + 3 * i]) for i in range(4)]
                ph = float(line[16:25]) * DEG
                x = float(line[26:35])
                m[12] += iz
                for i in range(4):
                    m[DEL_IDX[i]] += ilu[i]
                main = False
            else:  # planetary: (11i3,1x,f9.5,1x,f9.5)
                ipla = [int(line[3 * i:3 * i + 3]) for i in range(11)]
                ph = float(line[34:43]) * DEG
                x = float(line[44:53])
                if ific < 16:  # table 1: 8 planets, then D, l, F
                    for i in range(8):
                        m[4 + i] += ipla[i]
                    m[0] += ipla[8]
                    m[2] += ipla[9]
                    m[1] += ipla[10]
                else:  # table 2: 7 planets, then D, l', l, F
                    for i in range(7):
                        m[4 + i] += ipla[i]
                    for i in range(4):
                        m[DEL_IDX[i]] += ipla[7 + i]
                main = False
            amp = x * A0 / ATH if iv == 3 else x / RAD
            out[iv].append((m, amp, ph, k, main))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="data/sources/elp82b")
    ap.add_argument("--cut", type=float, default=1e-7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sums = fetch(a.cache)
    for f in FILES:
        print(f"{sums[f]}  {f}")
    terms = parse(a.cache)
    keep = {}
    for iv, rows in terms.items():
        sc = 1.0 / ATH if iv == 3 else 1.0
        keep[iv] = [r for r in rows if abs(r[1]) * sc * 20.0 ** r[3] > a.cut]
    n = sum(len(v) for v in keep.values())

    def poly(p):
        return f"({float(p[0])!r}, [{', '.join(repr(float(v)) for v in p[1])}])"

    with open(a.out, "w") as f:
        f.write(f"""// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
// Data source: ELP 2000-82B (Chapront-Touzé & Chapront 1983, 1988), series files
// ELP1..ELP36, IMCCE - Observatoire de Paris, ftp://ftp.imcce.fr/pub/ephem/moon/elp82b/
// (checksums in docs/FICTITIOUS_BODIES.md).
// @generated by tools/gen_intp_series.py — do not edit.
//! ELP 2000-82B series (main problem; Earth and Moon figure, planetary, tidal,
//! relativistic and solar-eccentricity perturbations), truncated at
//! |A|·20^k > {a.cut:g} rad (distance /384747.98 km): {n} terms. Each term is
//! A·T^k·sin(Σ mᵢ·argᵢ + φ) over the 13 arguments (D, F, l, l', Me, V, EMB, Ma,
//! J, S, U, N, ζ); longitude and latitude in radians, distance in km (the
//! authors' a0/ath factor applied). Main-problem terms (`true`) take the full
//! polynomial arguments [`ARGS`], perturbation terms their linear part, as in
//! the authors' routine ELP82B_2. Polynomials: constant (rad) and T¹..T⁴ (rad),
//! T in Julian centuries TDB from J2000.

// Measured coefficients, not mathematical constants.
#![allow(clippy::approx_constant)]

pub const W1: (f64, [f64; 4]) = {poly(W[0])};
pub const W2: (f64, [f64; 4]) = {poly(W[1])};
pub const ARGS: [(f64, [f64; 4]); 13] = [
""")
        for p in ARGS:
            f.write(f"    {poly(p)},\n")
        f.write("];\n")
        for iv, name in ((1, "LON"), (2, "LAT"), (3, "DIST")):
            rows = keep[iv]
            f.write(f"\npub const {name}: [([i8; 13], f64, f64, u8, bool); {len(rows)}] = [\n")
            for m, A, p, k, mn in rows:
                f.write(f"    ([{', '.join(map(str, m))}], {float(A)!r}, {float(p)!r}, {k}, {str(mn).lower()}),\n")
            f.write("];\n")
    print("wrote", a.out, {k: len(v) for k, v in keep.items()})


if __name__ == "__main__":
    main()
