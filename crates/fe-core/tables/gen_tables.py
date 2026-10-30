#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Generate the Rust constant tables used by `time` and `frames`.

Outputs (overwritten):
    crates/fe-core/src/time/tables.rs
    crates/fe-core/src/frames/tables.rs

Inputs
------
* Delta T (shipped in `tables/sources/`, small text files, fetched once):
  - SMH2016-Table-S15.txt   Stephenson, Morrison & Hohenkerk (2016), Proc. R. Soc. A 472:20160404,
                            supplementary "ASCII Data Text Files", figshare doi:10.6084/m9.figshare.4290866.v1
                            (https://ndownloader.figshare.com/files/6998699, file Table-S15.txt)
  - SMH2020-Table-S15.2020.txt  Morrison, Stephenson, Hohenkerk & Zawilski (2021) "Addendum 2020 to
                            'Measurement of the Earth's rotation: 720 BC to AD 2015'", Proc. R. Soc. A 477:20200776,
                            figshare doi:10.6084/m9.figshare.13885863.v2 (https://ndownloader.figshare.com/files/26513569)
  - USNO-historic_deltat.data  https://maia.usno.navy.mil/ser7/historic_deltat.data (semi-annual, 1657-1984.5)
  - USNO-deltat.data        https://maia.usno.navy.mil/ser7/deltat.data (monthly, observed, 1973-02 .. latest)
* IAU series (nutation 2000A/2000B, s06, EE complementary terms, Fairhead-Bretagnon TDB-TT):
  numeric coefficient tables from ERFA (BSD-3, derived from IAU SOFA), pinned commit
  ERFA_SHA below, downloaded to a cache directory on first run. The coefficients are the
  published IERS Conventions / MHB2000 / Fairhead & Bretagnon (1990) values.

Usage:  python gen_tables.py [--erfa-cache DIR]
"""
import argparse
import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "sources")
CRATE = os.path.dirname(HERE)
ERFA_SHA = "1a8044cde5b7763295d472a6443387239127c6c8"
ERFA_URL = "https://raw.githubusercontent.com/liberfa/erfa/{sha}/src/{name}.c"


# ----------------------------------------------------------------------------- ERFA parsing
def erfa_text(cache, name):
    os.makedirs(cache, exist_ok=True)
    p = os.path.join(cache, f"{ERFA_SHA[:12]}_{name}.c")
    if not os.path.exists(p):
        with urllib.request.urlopen(ERFA_URL.format(sha=ERFA_SHA, name=name), timeout=60) as r:
            data = r.read()
        with open(p, "wb") as f:
            f.write(data)
    return open(p, encoding="utf-8", errors="replace").read()


def c_block(txt, start):
    i = txt.index(start)
    j = txt.index("};", i)
    body = txt[i + len(start):j]
    return re.sub(r"/\*.*?\*/", "", body, flags=re.S)


NUM = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")


def rows(txt, start, n):
    v = [float(x) for x in NUM.findall(c_block(txt, start))]
    assert len(v) % n == 0, (start, len(v))
    return [v[i:i + n] for i in range(0, len(v), n)]


def fnum(x):
    """Shortest round-trip float literal."""
    s = repr(float(x))
    if "e" in s or "E" in s:
        m, e = s.split("e")
        if "." not in m:
            m += ".0"
        return f"{m}e{e}"
    if "." not in s:
        s += ".0"
    return s


def ints(v):
    return "[" + ", ".join(str(int(round(x))) for x in v) + "]"


# ----------------------------------------------------------------------------- Delta T inputs
def smh_table(path):
    out = []
    for line in open(path):
        p = line.split()
        if len(p) == 7 and p[0].isdigit():
            out.append([float(x) for x in p[1:]])
    return out


def jd_greg(y, m, d):
    # Gregorian calendar date at 0h UT -> JD (Fliegel & Van Flandern style, integer arithmetic).
    a = (14 - m) // 12
    yy = y + 4800 - a
    mm = m + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    return jdn - 0.5


def decimal_year_to_jd(yf):
    y = int(yf // 1)
    frac = yf - y
    j0 = jd_greg(y, 1, 1)
    j1 = jd_greg(y + 1, 1, 1)
    return j0 + frac * (j1 - j0)


def observed_table():
    """(jd_ut, delta_t) pairs: USNO historic semi-annual 1955.0..1973.0 then monthly 1973-02-01.."""
    out = []
    for line in open(os.path.join(SRC, "USNO-historic_deltat.data")).readlines()[2:]:
        p = line.split()
        if len(p) < 2:
            continue
        y = float(p[0])
        if 1955.0 <= y <= 1973.0:
            out.append((decimal_year_to_jd(y), float(p[1])))
    for line in open(os.path.join(SRC, "USNO-deltat.data")):
        p = line.split()
        if len(p) == 4:
            out.append((jd_greg(int(p[0]), int(p[1]), int(p[2])), float(p[3])))
    jds = [j for j, _ in out]
    assert jds == sorted(jds) and len(set(jds)) == len(jds)
    return out


# ----------------------------------------------------------------------------- writers
SPDX = ("// SPDX-License-Identifier: (AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial) AND BSD-3-Clause\n"
        "// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI\n")
DATA_TIME = ("// Data: Delta T from Stephenson, Morrison & Hohenkerk 2016/2020 (Table S15) and USNO observed values; "
             "TDB-TT series from ERFA dtdb.c (BSD-3-Clause, Copyright (C) 2013-2021 NumFOCUS Foundation; derived from IAU SOFA). See NOTICE.\n")
DATA_FRAMES = "// Data: coefficients transcribed from ERFA (BSD-3-Clause, Copyright (C) 2013-2021 NumFOCUS Foundation; derived from IAU SOFA). See NOTICE.\n"
HEADER = "// @generated by crates/fe-core/tables/gen_tables.py -- do not edit by hand.\n"


def write_time(cache):
    s16 = smh_table(os.path.join(SRC, "SMH2016-Table-S15.txt"))
    s20 = smh_table(os.path.join(SRC, "SMH2020-Table-S15.2020.txt"))
    obs = observed_table()
    dtdb = erfa_text(cache, "dtdb")
    fair = rows(dtdb, "fairhd[787][3] = {", 3)
    assert len(fair) == 787
    power = [0] * 474 + [1] * 205 + [2] * 85 + [3] * 20 + [4] * 3
    # Truncation: keep terms whose amplitude * 2^power >= 3e-8 s (|t| <= 2 millennia over
    # years 0..2400). Measured error vs the full series: max 0.42 us, rms 0.09 us.
    trunc = [(a, f, p, n) for (a, f, p), n in zip(fair, power) if a * 2.0 ** n >= 3e-8]

    o = [SPDX, DATA_TIME, HEADER, "//! Numeric tables for `time` (Delta T splines, observed Delta T, TDB-TT series).\n\n"]
    o.append("/// Spline row: (K_i, K_{i+1}, a0, a1, a2, a3); DT = a0 + a1 t + a2 t^2 + a3 t^3, t = (Y-K_i)/(K_{i+1}-K_i).\n")
    o.append("pub type SplineRow = (f64, f64, f64, f64, f64, f64);\n\n")
    for name, tab, cite in [
        ("SMH2016", s16, "Stephenson, Morrison & Hohenkerk 2016, Table S15 (doi:10.1098/rspa.2016.0404)"),
        ("SMH2020", s20, "Morrison, Stephenson, Hohenkerk & Zawilski 2021, Addendum 2020, Table S15.2020 (doi:10.1098/rspa.2020.0776)"),
    ]:
        o.append(f"/// {cite}.\npub const {name}: [SplineRow; {len(tab)}] = [\n")
        for r in tab:
            o.append("    (" + ", ".join(fnum(x) for x in r) + "),\n")
        o.append("];\n\n")
    o.append(
        "/// Observed Delta T (TT-UT1, s) at JD(UT1): USNO historic_deltat.data semi-annual values 1955.0..1973.0,\n"
        "/// then USNO deltat.data monthly values (https://maia.usno.navy.mil/ser7/).\n"
        f"pub const OBSERVED: [(f64, f64); {len(obs)}] = [\n")
    for j, v in obs:
        o.append(f"    ({fnum(j)}, {fnum(v)}),\n")
    o.append("];\n\n")
    o.append("/// Fairhead & Bretagnon (1990) TDB-TT series, as in IERS Conventions / SOFA dtdb:\n"
             "/// (amplitude s, frequency rad per Julian millennium, phase rad, power of T).\n")
    o.append(f"pub static FB_FULL: [(f64, f64, f64, u8); {len(fair)}] = [\n")
    for (a, f, p), n in zip(fair, power):
        o.append(f"    ({fnum(a)}, {fnum(f)}, {fnum(p)}, {n}),\n")
    o.append("];\n\n")
    o.append("/// Truncated subset of FB_FULL (amplitude * 2^power >= 3e-8 s): max error 0.42 us over years 0..2400.\n")
    o.append(f"pub const FB_TRUNC: [(f64, f64, f64, u8); {len(trunc)}] = [\n")
    for a, f, p, n in trunc:
        o.append(f"    ({fnum(a)}, {fnum(f)}, {fnum(p)}, {n}),\n")
    o.append("];\n")
    path = os.path.join(CRATE, "src", "time", "tables.rs")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("".join(o))
    print("wrote", path, len(obs), "observed rows;", len(trunc), "TDB terms kept")


def write_frames(cache):
    n00a = erfa_text(cache, "nut00a")
    ls = rows(n00a, "} xls[] = {", 11)
    pl = rows(n00a, "} xpl[] = {", 17)
    n00b = rows(erfa_text(cache, "nut00b"), "} x[] = {", 11)
    s06t = erfa_text(cache, "s06")
    sp = [float(x) for x in NUM.findall(c_block(s06t, "sp[] = {"))]
    s = [rows(s06t, f"TERM s{k}[] = {{", 10) for k in range(5)]
    eect = erfa_text(cache, "eect00")
    e0 = rows(eect, "TERM e0[] = {", 10)
    e1 = rows(eect, "TERM e1[] = {", 10)
    assert (len(ls), len(pl), len(n00b)) == (678, 687, 77)
    assert [len(x) for x in s] == [33, 3, 25, 4, 1] and (len(e0), len(e1)) == (33, 1)

    o = [SPDX, DATA_FRAMES, HEADER, "//! IAU nutation / CIO-locator / equation-of-equinoxes series coefficients\n"
         "//! (IERS Conventions 2003/2010, MHB2000; transcribed from ERFA, BSD-3, NumFOCUS, derived from IAU SOFA).\n\n"]
    o.append("/// Luni-solar term: multipliers of (l, l', F, D, Om); dpsi = (sp + spt t) sin + cp cos;\n"
             "/// deps = (ce + cet t) cos + se sin. Units 0.1 uas.\n"
             "pub type LsTerm = ([i8; 5], f64, f64, f64, f64, f64, f64);\n")
    o.append("/// Planetary term: multipliers of (l, F, D, Om, Me, Ve, E, Ma, Ju, Sa, Ur, Ne, pA);\n"
             "/// dpsi = sp sin + cp cos; deps = se sin + ce cos. Units 0.1 uas.\n"
             "pub type PlTerm = ([i8; 13], f64, f64, f64, f64);\n")
    o.append("/// Series term with multipliers of (l, l', F, D, Om, LVe, LE, pA) and (sin, cos) coefficients (arcsec).\n"
             "pub type FaTerm = ([i8; 8], f64, f64);\n\n")

    def ls_rows(tab, name, doc):
        o.append(f"/// {doc}\npub {'static' if name == 'NUT00A_LS' else 'const'} {name}: [LsTerm; {len(tab)}] = [\n")
        for r in tab:
            o.append(f"    ({ints(r[:5])}, " + ", ".join(fnum(x) for x in r[5:]) + "),\n")
        o.append("];\n\n")

    ls_rows(ls, "NUT00A_LS", "IAU 2000A (MHB2000) luni-solar series, 678 terms.")
    o.append(f"/// IAU 2000A planetary series, 687 terms.\npub static NUT00A_PL: [PlTerm; {len(pl)}] = [\n")
    for r in pl:
        o.append(f"    ({ints(r[:13])}, " + ", ".join(fnum(x) for x in r[13:]) + "),\n")
    o.append("];\n\n")
    ls_rows(n00b, "NUT00B_LS", "IAU 2000B luni-solar series, 77 terms (McCarthy & Luzum 2003).")
    o.append(f"/// s + XY/2 polynomial part (arcsec), IAU 2006/2000A.\npub const S06_POLY: [f64; 6] = [{', '.join(fnum(x) for x in sp)}];\n\n")
    for k in range(5):
        o.append(f"/// s + XY/2 series, terms of order t^{k}.\npub const S06_T{k}: [FaTerm; {len(s[k])}] = [\n")
        for r in s[k]:
            o.append(f"    ({ints(r[:8])}, {fnum(r[8])}, {fnum(r[9])}),\n")
        o.append("];\n\n")
    for name, tab, k in [("EECT_T0", e0, 0), ("EECT_T1", e1, 1)]:
        o.append(f"/// Equation-of-the-equinoxes complementary terms (IERS 2003), order t^{k}.\npub const {name}: [FaTerm; {len(tab)}] = [\n")
        for r in tab:
            o.append(f"    ({ints(r[:8])}, {fnum(r[8])}, {fnum(r[9])}),\n")
        o.append("];\n\n")
    path = os.path.join(CRATE, "src", "frames", "tables.rs")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("".join(o))
    print("wrote", path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--erfa-cache", default=os.path.join(HERE, ".erfa_cache"))
    a = ap.parse_args()
    write_time(a.erfa_cache)
    write_frames(a.erfa_cache)
