#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Fit our own Chebyshev series for minor bodies and write `fepk-minor-v2.bin`.

    python tools/fit_minor.py [--tol-km 0.1] [--bodies Chiron,Eros] [--out PATH]

For each body the source SPK (tools/spk.py: type 2 from sb441-n16, type 21
from Horizons) is sampled and re-fitted as heliocentric (target wrt the Sun,
NAIF 10) ICRF positions in FEPK v1 series (docs/FORMAT.md).  The engine
chains them to the SSB through the main pack's Sun (10 wrt 0).

Targets use NAIF small-body ids 2000000 + number (Ceres 2000001, Chiron
2002060, Eris 2136199, ...).

Fit
---
* Each record is the Chebyshev interpolant of degree N-1 through the source
  at the N Chebyshev-Gauss-Lobatto points s_k = cos(pi k/(N-1)).  The end
  points s = +-1 are nodes, so adjacent records agree with the source - and
  with each other - exactly at their common boundary (position continuity up
  to rounding; velocity jumps are of the order of the fit error / record).
* Centre: Sun (10) or SSB (0), whichever needs fewer coefficients.  Heliocentric
  coordinates carry the Sun's reflex motion (Mercury ~10 km at 88 d, Venus
  ~270 km, EMB ~450 km), so the barycentric form usually wins.  Barycentric =
  JPL heliocentric + the main pack's Sun wrt SSB, which the engine subtracts
  again exactly; barycentric series restart at the main pack's DE441 -> DE442
  seam (the Sun jumps by ~1 km there).
* Record length (power of two days, 1..1024) and N (4..MAX_N) are chosen per
  body to minimise size subject to the tolerance on a sample of records,
  then the whole span is fitted and checked; N is raised until every record
  passes.
* Check ("out of sample"): in every record, the midpoints between adjacent
  nodes (where interpolation error peaks) plus 8 uniformly random points.
  None of these are fit nodes.

Tolerance: max 3-D error <= tol_km everywhere.  The spec allows the looser of
1 km and 1e-5 arcsec as seen from Earth; 1e-5" equals 1 km only at 138 AU,
so for all eight bodies (geocentric distance < 100 AU) 1 km is the looser
bound.  We fit to --tol-km (default 0.1 km, 10x margin) and record the
achieved max/p99/rms and an upper bound of the angular error seen from Earth
(error / (r_helio_min - 1.0167 AU)) in the manifest.

The span ends at JD 2598072.5 (common end of the main pack) for every body and
starts at JD 1720624.5 (main pack start) for the sb441 bodies and at the
Horizons limit 2305426.5 (1599-12-11) for the others.  When a span (or a
piece between seams) is not a whole number of records, a one-record series
of the same pair covers the remainder (FEPK allows contiguous series per
pair).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spk import SPK, J2000_JD, S_PER_DAY  # noqa: E402
from fepk import Fepk  # noqa: E402
import fetch_minor  # noqa: E402

FIT_VERSION = "1.0.0"
MAGIC = b"FEPK\0\0\0\x01"
HEADER_LEN = 64
SERIES_LEN = 64
AU_KM = 149597870.7
MAX_N = 32
DEFAULT_DATA = fetch_minor.DEFAULT_DATA

MAIN_START = 1720624.5   # main pack (fepk-de441-de442-v1) earliest series start
MAIN_END = 2598072.5     # main pack common coverage end (Moon/Earth; Sun ends 2598080.5)
HZ_START = 2305426.5     # 1599-12-11 0h TDB, Horizons small-body limit

# name, MPC number, Swiss/API id, source
BODIES = [
    ("Ceres", 1, 17, "sb441"),
    ("Pallas", 2, 18, "sb441"),
    ("Juno", 3, 19, "sb441"),
    ("Vesta", 4, 20, "sb441"),
    ("Chiron", 2060, 15, "horizons"),
    ("Pholus", 5145, 16, "horizons"),
    ("Eris", 136199, 10000 + 136199, "horizons"),
    ("Eros", 433, 10000 + 433, "horizons"),
]


def naif_id(number: int) -> int:
    return 2000000 + number


class Source:
    """ICRF position (km) of one body from its JPL SPK, wrt the Sun (center 10)
    or, with `center = 0`, wrt the SSB (source heliocentric + the main pack's
    Sun wrt SSB, i.e. exactly what the engine's chain subtracts again)."""

    center = 10
    main: Fepk | None = None
    # [a, b] of the piece being fitted: the main-pack Sun is taken from the
    # series covering it, so a piece ending on the DE441/DE442 seam uses the
    # DE441 side at its last node (the engine's own rule picks DE442 there,
    # but also the next piece of this pack, which was fitted with DE442).
    piece: tuple[float, float] | None = None

    def __init__(self, name: str, number: int, kind: str, data: str):
        minor = os.path.join(data, "sources", "minor")
        self.name, self.number, self.kind = name, number, kind
        if kind == "sb441":
            self.spk = SPK.open(fetch_minor.SB441_URL, os.path.join(minor, "sb441-n16_cache"))
            self.target = 2000000 + number
            self.file = "sb441-n16.bsp"
            self.start, self.end = MAIN_START, MAIN_END
            self._preload()
        else:
            nm = dict(fetch_minor.HORIZONS_BODIES)[number]
            path = fetch_minor.hz_spk_path(minor, number, nm)
            self.spk = SPK.open(path)
            self.target = 20000000 + number  # Horizons SPK id convention
            self.file = os.path.basename(path)
            self.start, self.end = HZ_START, MAIN_END
        segs = self.spk.segments_for(self.target, 10)
        lo, hi = min(s.start_jd for s in segs), max(s.end_jd for s in segs)
        if not (lo <= self.start and self.end <= hi):
            raise ValueError(f"{name}: source covers JD {lo}..{hi}, need {self.start}..{self.end}")
        self.n_evals = 0

    def _preload(self):
        """Read every type-2 record of the span once (bit-exact) into memory."""
        self._t2 = []
        for s in self.spk.segments_for(self.target, 10):
            a, b = max(self.start, s.start_jd), min(self.end, s.end_jd)
            if a >= b:
                continue
            i0, recs = s.records_for_jd(a, b)
            self._t2.append((s, i0, recs))

    def pos(self, jd: np.ndarray) -> np.ndarray:
        jd = np.asarray(jd, dtype=np.float64).ravel()
        h = self.helio(jd)
        if self.center == 0:
            return h + self.main.pos(10, 0, jd, within=self.piece)
        return h

    def helio(self, jd: np.ndarray) -> np.ndarray:
        self.n_evals += len(jd)
        et = (jd - J2000_JD) * S_PER_DAY
        if self.kind != "sb441":
            return self.spk.state_et(self.target, 10, et)[0]
        out = np.full((len(jd), 3), np.nan)
        done = np.zeros(len(jd), bool)
        for s, i0, recs in reversed(self._t2):  # later segment wins on a shared boundary
            m = (~done) & (et >= s.start_et) & (et <= s.end_et)
            if not m.any():
                continue
            idx = np.clip(((et[m] - s.init) // s.intlen).astype(np.int64), 0, s.n - 1) - i0
            idx = np.clip(idx, 0, len(recs) - 1)
            r = recs[idx]
            x = (et[m] - r[:, 0]) / r[:, 1]
            c = r[:, 2:].reshape(-1, 3, s.ncoef)
            # Clenshaw-free direct sum (ncoef = 8): T_k(x) table
            t = np.empty((len(x), s.ncoef))
            t[:, 0] = 1.0
            t[:, 1] = x
            for k in range(2, s.ncoef):
                t[:, k] = 2 * x * t[:, k - 1] - t[:, k - 2]
            out[m] = np.einsum("nck,nk->nc", c, t)
            done |= m
        if not done.all():
            raise ValueError("epoch outside source")
        return out


class IntegratedSource:
    """Our own backward integration (tools/integrate_minor.py) of a body before
    the JPL span: a barycentric grid t, x, v, a (km, km/day, km/day^2) read from
    data/sources/minor/int_<number>_<name>.npy, evaluated by quintic Hermite
    interpolation between grid points (error ~ h^6 x^(6) / 46080, far below the
    integration's own error)."""

    kind = "integration"
    main: Fepk | None = None

    def __init__(self, name: str, number: int, data: str):
        self.name, self.number = name, number
        nm = dict(fetch_minor.HORIZONS_BODIES)[number]
        self.file = f"int_{number}_{nm}.npy"
        self.g = np.load(os.path.join(data, "sources", "minor", self.file))
        self.start, self.end = float(self.g[0, 0]), float(self.g[-1, 0])
        if self.end != HZ_START:
            raise ValueError(f"{name}: integration ends at {self.end}, expected the JPL start {HZ_START}")
        self.center = 0
        self.piece = None
        self.n_evals = 0

    def bary(self, jd: np.ndarray) -> np.ndarray:
        jd = np.asarray(jd, dtype=np.float64).ravel()
        g = self.g
        t = g[:, 0]
        if jd.min() < t[0] or jd.max() > t[-1]:
            raise ValueError("epoch outside integration")
        i = np.clip(np.searchsorted(t, jd, side="right") - 1, 0, len(t) - 2)
        h = (t[i + 1] - t[i])[:, None]
        s = ((jd - t[i]) / h[:, 0])[:, None]
        p0, v0, a0 = g[i, 1:4], g[i, 4:7] * h, g[i, 7:10] * h * h
        p1, v1, a1 = g[i + 1, 1:4], g[i + 1, 4:7] * h, g[i + 1, 7:10] * h * h
        s2, s3 = s * s, s * s * s
        s4, s5 = s3 * s, s3 * s2
        h0 = 1 - 10 * s3 + 15 * s4 - 6 * s5
        h1 = s - 6 * s3 + 8 * s4 - 3 * s5
        h2 = 0.5 * (s2 - 3 * s3 + 3 * s4 - s5)
        h3 = 10 * s3 - 15 * s4 + 6 * s5
        h4 = -4 * s3 + 7 * s4 - 3 * s5
        h5 = 0.5 * (s3 - 2 * s4 + s5)
        return h0 * p0 + h1 * v0 + h2 * a0 + h3 * p1 + h4 * v1 + h5 * a1

    def helio(self, jd: np.ndarray) -> np.ndarray:
        jd = np.asarray(jd, dtype=np.float64).ravel()
        self.n_evals += len(jd)
        return self.bary(jd) - self.main.pos(10, 0, jd)

    def pos(self, jd: np.ndarray) -> np.ndarray:
        jd = np.asarray(jd, dtype=np.float64).ravel()
        self.n_evals += len(jd)
        b = self.bary(jd)
        if self.center == 0:
            return b
        return b - self.main.pos(10, 0, jd, within=self.piece)


# --------------------------------------------------------------------------
# Chebyshev-Lobatto interpolation
# --------------------------------------------------------------------------
def lobatto_nodes(n: int) -> np.ndarray:
    return np.cos(np.pi * np.arange(n) / (n - 1))  # +1 .. -1


def lobatto_matrix(n: int) -> np.ndarray:
    """M (n x n) with coeffs = M @ values-at-nodes (DCT-I)."""
    k = np.arange(n)
    c = np.cos(np.pi * np.outer(k, k) / (n - 1))  # [j, k]
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    m = (2.0 / (n - 1)) * c * w[None, :]
    m[0] *= 0.5
    m[-1] *= 0.5
    return m


def cheb_eval(coef: np.ndarray, s: np.ndarray) -> np.ndarray:
    """coef (nrec, 3, n), s (nrec, m) -> (nrec, m, 3)."""
    n = coef.shape[2]
    t = np.empty(s.shape + (n,))
    t[..., 0] = 1.0
    t[..., 1] = s
    for k in range(2, n):
        t[..., k] = 2 * s * t[..., k - 1] - t[..., k - 2]
    return np.einsum("rck,rmk->rmc", coef, t)


def fit_records(src: Source, start: float, interval: float, idx: np.ndarray, n: int) -> np.ndarray:
    src.piece = (start + idx.min() * interval, start + (idx.max() + 1) * interval)
    nodes = lobatto_nodes(n)
    t0 = start + idx * interval
    jd = t0[:, None] + (nodes[None, :] + 1.0) * 0.5 * interval
    vals = src.pos(jd).reshape(len(idx), n, 3)
    m = lobatto_matrix(n)
    return np.einsum("jk,rkc->rcj", m, vals)  # (nrec, 3, n)


def check_points(n: int, nrec: int, rng: np.random.Generator) -> np.ndarray:
    nodes = lobatto_nodes(n)
    mids = 0.5 * (nodes[:-1] + nodes[1:])
    rnd = rng.uniform(-1, 1, (nrec, 8))
    return np.concatenate([np.broadcast_to(mids, (nrec, len(mids))), rnd], axis=1)


def record_errors(src, start, interval, idx, coef, rng) -> np.ndarray:
    """Max 3-D error (km) per record on the out-of-sample check points; also returns all errors."""
    src.piece = (start + idx.min() * interval, start + (idx.max() + 1) * interval)
    s = check_points(coef.shape[2], len(idx), rng)
    jd = (start + idx * interval)[:, None] + (s + 1.0) * 0.5 * interval
    truth = src.pos(jd).reshape(len(idx), s.shape[1], 3)
    err = np.linalg.norm(cheb_eval(coef, s) - truth, axis=2)
    return err


def layout(src: Source, interval: float) -> list[dict]:
    """Series covering [src.start, src.end]: whole records of `interval` days plus a
    one-record tail, restarted at every seam of the main pack's Sun (10 wrt 0)
    when fitting barycentrically -- the DE441 -> DE442 Sun has a ~1 km jump there,
    which no smooth record may straddle."""
    cuts = [src.start, src.end]
    # Barycentric JPL fits add the main-pack Sun; integrated spans subtract it for a
    # heliocentric centre: either way the Sun's seam jump enters the fitted data.
    if src.center == 0 or getattr(src, "kind", "") == "integration":
        ss = [x for x in Source.main.series if x["target"] == 10 and x["center"] == 0]
        cuts += [x["start"] for x in ss if src.start < x["start"] < src.end]
    cuts = sorted(set(cuts))
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        nrec = int((b - a) // interval)
        if nrec:
            out.append(dict(start=a, interval=interval, n_records=nrec))
        t = a + nrec * interval
        if b - t > 1e-9:
            out.append(dict(start=t, interval=b - t, n_records=1))
    return out


def choose(src: Source, tol: float, rng, sample: int = 160):
    """(interval, n) minimising coefficients per day with sampled max error <= tol/2."""
    best = None
    for p in range(0, 11):
        interval = float(2 ** p)
        main = [x for x in layout(src, interval) if x["interval"] == interval]
        if not main:
            break
        for n in range(4, MAX_N + 1):
            e = 0.0
            for x in main:
                k = max(2, int(sample * x["n_records"] * interval / (src.end - src.start)))
                idx = np.unique(np.linspace(0, x["n_records"] - 1, min(k, x["n_records"])).round().astype(np.int64))
                coef = fit_records(src, x["start"], interval, idx, n)
                e = max(e, record_errors(src, x["start"], interval, idx, coef, rng).max())
            if e <= tol * 0.5:  # margin for records the sample misses
                cost = n / interval
                if best is None or cost < best[0]:
                    best = (cost, interval, n, e)
                break
    if best is None:
        raise RuntimeError(f"{src.name}: no (interval, n) reaches {tol} km")
    return best[1], best[2]


def refine_series(src, x: dict, n: int, tol: float, rng, depth: int = 0):
    """Fit grid series `x`; records failing `tol` are replaced, run by run, by
    series of half-length records (recursively). Returns (series list, errors)."""
    idx = np.arange(x["n_records"])
    coef = fit_records(src, x["start"], x["interval"], idx, n)
    err = record_errors(src, x["start"], x["interval"], idx, coef, rng)
    bad = err.max(axis=1) > tol
    if not bad.any():
        return [dict(x, coef=coef)], [err.ravel()]
    if depth > 12:
        raise RuntimeError(f"{src.name}: cannot reach {tol} km near JD {x['start']}")
    out, errs = [], []
    k = 0
    while k < len(idx):
        j = k
        while j < len(idx) and bad[j] == bad[k]:
            j += 1
        sub = dict(start=x["start"] + k * x["interval"], interval=x["interval"], n_records=j - k)
        if bad[k]:
            half = dict(start=sub["start"], interval=x["interval"] / 2, n_records=2 * (j - k))
            s2, e2 = refine_series(src, half, n, tol, rng, depth + 1)
            out += s2
            errs += e2
        else:
            out.append(dict(sub, coef=coef[k:j]))
            errs.append(err[k:j].ravel())
        k = j
    return out, errs


def fit_body(src: Source, tol: float, seed: int = 7, center: int | None = None, refine: bool = False):
    """`refine`: keep the chosen degree and halve the record length locally where
    records fail (close planetary encounters in integrated spans)."""
    rng = np.random.default_rng(seed)
    # Heliocentric coordinates of a distant body carry the Sun's reflex motion
    # (Mercury ~10 km at 88 d, Venus ~270 km, EMB ~450 km): barycentric ones are
    # far smoother there. Pick whichever centre needs fewer coefficients
    # (or use `center`, e.g. for a body's second span: one centre per pair).
    picks = {}
    for center in ((10, 0) if center is None else (center,)):
        src.center = center
        interval, n = choose(src, tol, rng)
        picks[center] = (n / interval, interval, n)
    src.center = min(picks, key=lambda c: picks[c][0])
    _, interval, n = picks[src.center]
    lay = layout(src, interval)
    while refine:
        series, errs, vj = [], [], 0.0
        for x in lay:
            s2, e2 = refine_series(src, x, n, tol, rng)
            series += s2
            errs += e2
        for x in series:
            vj = max(vj, velocity_jumps(x["coef"], x["interval"]))
        allerr = np.concatenate(errs)
        break
    while not refine:
        series, errs, vj = [], [], 0.0
        for x in lay:
            idx = np.arange(x["n_records"])
            coef = fit_records(src, x["start"], x["interval"], idx, n)
            err = record_errors(src, x["start"], x["interval"], idx, coef, rng)
            series.append(dict(x, coef=coef))
            errs.append(err.ravel())
            vj = max(vj, velocity_jumps(coef, x["interval"]))
        allerr = np.concatenate(errs)
        if allerr.max() <= tol or n == MAX_N:
            break
        n += 1
    if allerr.max() > tol:
        raise RuntimeError(f"{src.name}: max error {allerr.max():.3f} km > {tol} km at N={MAX_N}")
    # heliocentric distance range on a coarse grid (for the angular bound)
    src.piece = None
    grid = np.linspace(src.start, src.end, 4001)
    r = np.linalg.norm(src.helio(grid), axis=1)
    rmin = r.min()
    geo_min_km = max(rmin - 1.0167 * AU_KM, 0.01 * AU_KM)
    stats = {
        "center": src.center,
        "center_choice_coeffs_per_day": {str(c): v[0] for c, v in picks.items()},
        "interval_days": interval, "n_coeff": n,
        "min_interval_days": float(min(x["interval"] for x in series)),
        "n_records": int(sum(x["n_records"] for x in series)), "n_series": len(series),
        "bytes": int(sum(x["n_records"] for x in series) * 3 * n * 8),
        "check_points": int(allerr.size),
        "max_km": float(allerr.max()), "p99_km": float(np.percentile(allerr, 99)),
        "rms_km": float(np.sqrt(np.mean(allerr ** 2))),
        "helio_dist_au": [float(rmin / AU_KM), float(r.max() / AU_KM)],
        "max_arcsec_from_earth_bound": float(allerr.max() / geo_min_km * 206264.806),
        "boundary_velocity_jump_km_day_max": float(vj),
        "source_evaluations": src.n_evals,
    }
    return series, stats


def velocity_jumps(coef: np.ndarray, interval: float) -> float:
    n = coef.shape[2]
    k = np.arange(n)
    d_end = k ** 2              # T_k'(1)
    d_start = (-1.0) ** (k + 1) * k ** 2  # T_k'(-1)
    v_end = coef @ d_end * 2 / interval
    v_start = coef @ d_start * 2 / interval
    if len(coef) < 2:
        return 0.0
    return float(np.linalg.norm(v_end[:-1] - v_start[1:], axis=1).max())


# --------------------------------------------------------------------------
# FEPK writer (docs/FORMAT.md v1)
# --------------------------------------------------------------------------
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def write_pack(path: str, rows: list[dict], manifest: dict):
    rows = sorted(rows, key=lambda r: (r["target"], r["center"], r["start"]))
    for r in rows:
        r["nbytes"] = r["n_records"] * 3 * r["n_coeff"] * 8
    table_meta = [{k: r[k] for k in ("target", "center", "n_coeff", "n_records", "start", "interval",
                                     "source_index")} | {"end": r["start"] + r["n_records"] * r["interval"]}
                  for r in rows]
    manifest["series"] = table_meta
    man = json.dumps(manifest, indent=1).encode()
    table_off = HEADER_LEN
    man_off = table_off + len(rows) * SERIES_LEN
    data_off = (man_off + len(man) + 7) // 8 * 8
    total = data_off + sum(r["nbytes"] for r in rows)
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(MAGIC + struct.pack("<II", 1, len(rows)))
        f.write(struct.pack("<QQQQ", table_off, man_off, len(man), data_off) + bytes(16))
        off = data_off
        for r in rows:
            end = r["start"] + r["n_records"] * r["interval"]
            f.write(struct.pack("<iiII", r["target"], r["center"], r["n_coeff"], r["n_records"]))
            f.write(struct.pack("<dddQI", r["start"], r["interval"], end, off, r["source_index"]))
            f.write(bytes(12))
            off += r["nbytes"]
        f.write(man)
        f.write(bytes(data_off - man_off - len(man)))
        for r in rows:
            # [record][component][k], little-endian f64
            f.write(np.ascontiguousarray(r["coef"], dtype="<f8").tobytes())
        assert f.tell() == total
    os.replace(tmp, path)
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--tol-km", type=float, default=0.1)
    ap.add_argument("--bodies", help="comma-separated subset (default: all)")
    ap.add_argument("--out")
    ap.add_argument("--main-pack", help="main FEPK pack (Sun wrt SSB for barycentric fits)")
    ap.add_argument("--no-integration", action="store_true",
                    help="JPL spans only (no integrated extension before 1599-12-11)")
    args = ap.parse_args()
    main_pack = args.main_pack or os.path.join(args.data, "packs", "fepk-de441-de442-v1.bin")
    Source.main = Fepk(main_pack)
    main_meta = {"name": Source.main.manifest.get("name"), "file": os.path.basename(main_pack),
                 "created_utc": Source.main.manifest.get("created_utc")}
    out = args.out or os.path.join(args.data, "packs", "fepk-minor-v2.bin")
    want = set(args.bodies.split(",")) if args.bodies else None
    ledger = json.load(open(fetch_minor.LEDGER))
    ipath = os.path.join(fetch_minor.REPO, "docs", "MINOR_INTEGRATION.json")
    integrations = json.load(open(ipath))["bodies"] if os.path.exists(ipath) else {}
    IntegratedSource.main = Source.main
    hz_by_num = {e["number"]: e for e in ledger["horizons"]}

    sources = []   # manifest.sources
    rows = []
    bodies = {}
    for name, number, api_id, kind in BODIES:
        if want and name not in want:
            continue
        src = Source(name, number, kind, args.data)
        series, stats = fit_body(src, args.tol_km)
        if kind == "sb441":
            if not any(s["name"] == "sb441-n16.bsp" for s in sources):
                sb = ledger["sb441"]
                sources.append({"name": "sb441-n16.bsp", "url": sb["url"], "size_bytes": sb["size_bytes"],
                                "last_modified": sb["last_modified"], "kind": sb["kind"],
                                "reference": sb["reference"], "license": sb["license"]})
            si = next(i for i, s in enumerate(sources) if s["name"] == "sb441-n16.bsp")
        else:
            e = hz_by_num[number]
            sources.append({"name": e["file"], "sha256": e["sha256"], "size_bytes": e["size_bytes"],
                            "kind": e["kind"], "solution_id": e["solution_id"],
                            "solution_date": e["solution_date"], "observations": e["observations"],
                            "execution_date_pasadena": e["execution_date_pasadena"],
                            "request": e["request"], "license": e["license"]})
            si = len(sources) - 1
        tgt = naif_id(number)
        for s in series:
            rows.append(dict(target=tgt, center=stats["center"], n_coeff=s["coef"].shape[2], n_records=s["n_records"],
                             start=s["start"], interval=s["interval"], source_index=si, coef=s["coef"]))
        bodies[name] = {"naif_id": tgt, "center": stats["center"], "mpc_number": number, "api_id": api_id,
                        "source_index": si, "source": sources[si]["name"],
                        "solution_id": sources[si].get("solution_id", "SB441-N16 (DE441)"),
                        "coverage_jd_tdb": [src.start, src.end], "fit": stats}
        integ = integrations.get(name, {}).get("integration")
        if integ and not args.no_integration:
            isrc = IntegratedSource(name, number, args.data)
            iseries, istats = fit_body(isrc, args.tol_km, center=stats["center"], refine=True)
            sources.append({"name": isrc.file,
                            "kind": f"FrenchEphemeris integration from {sources[si]['solution_id']}",
                            "solution_id": sources[si]["solution_id"],
                            "initial_state": integ["initial_state"], "dmu": integ["dmu"],
                            "step_days": integ["step_days"], "span_jd_tdb": integ["span_jd_tdb"],
                            "model": "tools/integrate_minor.py (see docs/MINOR_INTEGRATION.json)",
                            "note": ("model extrapolation before JPL's first epoch; physical uncertainty grows "
                                     "backwards (docs/MINOR_BODIES_REPORT.md)")})
            ii = len(sources) - 1
            for s in iseries:
                rows.append(dict(target=tgt, center=istats["center"], n_coeff=s["coef"].shape[2],
                                 n_records=s["n_records"], start=s["start"], interval=s["interval"],
                                 source_index=ii, coef=s["coef"]))
            bodies[name]["coverage_jd_tdb"] = [isrc.start, src.end]
            bodies[name]["integrated_span"] = {"jd_tdb": [isrc.start, isrc.end], "source_index": ii,
                                               "source": isrc.file, "fit": istats}
            print(f"{name:<7} integrated span: {istats['n_records']} rec x {istats['interval_days']:g} d x "
                  f"N={istats['n_coeff']}  max {istats['max_km']:.4f} km", flush=True)
        print(f"{name:<7} {tgt}: {stats['n_records']} rec x {stats['interval_days']:g} d x N={stats['n_coeff']}"
              f"  max {stats['max_km']:.4f} km  p99 {stats['p99_km']:.4f}  rms {stats['rms_km']:.4f}"
              f"  <= {stats['max_arcsec_from_earth_bound']:.2e}\"  ({stats['check_points']} checks)", flush=True)

    manifest = {
        "name": "fepk-minor",
        "format_version": 1,
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builder": {"tool": "tools/fit_minor.py", "version": FIT_VERSION, "argv": sys.argv[1:]},
        "kind": ("own Chebyshev fit (Chebyshev-Gauss-Lobatto interpolation) of JPL small-body SPKs and, before "
                 "JPL's 1599-12-11 limit, of our own integration from the JPL state (sources[].kind)"),
        "tolerance_km": args.tol_km,
        "tolerance_note": ("spec: max error <= max(1 km, 1e-5 arcsec x geocentric distance); 1 km is the looser "
                           "bound for every body here; fitted to tolerance_km"),
        "frame": "ICRF (SPK frame 1), positions km, time JD TDB; center per body: Sun (10, chained to the "
                 "SSB through the main pack's 10 wrt 0) or SSB (0, = JPL heliocentric + main-pack Sun)",
        "requires": "fepk-de441-de442-v1 (Sun wrt SSB) for chaining",
        "barycentric_sun_from": main_meta,
        "bodies": bodies,
        "sources": sources,
        "provenance_ledger": "docs/MINOR_PROVENANCE.json",
        "license": fetch_minor.LICENSE_NOTE,
    }
    total = write_pack(out, rows, manifest)
    side = dict(manifest)
    side.pop("series", None)
    side["file"] = {"path": os.path.basename(out), "size_bytes": total, "sha256": sha256_file(out)}
    with open(out + ".manifest.json", "w") as f:
        json.dump(side, f, indent=1)
    print(f"wrote {out}: {total} bytes ({total / 1e6:.2f} MB), {len(rows)} series")


if __name__ == "__main__":
    main()
