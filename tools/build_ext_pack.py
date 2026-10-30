#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Build the extended-range FEPK v2 packs by refitting DE441 (docs/FORMAT.md, v2).

    python tools/build_ext_pack.py --era past     # DE441 start (-13200) .. main pack start (-1)
    python tools/build_ext_pack.py --era future   # main pack end (2401) .. DE441 end (17191)
    python tools/build_ext_pack.py --era future --max-chunks 1 --assemble-partial \
        --out data/packs/fepk-de441-ext-near-future-v1.bin   # 2401 .. ~2580 only (~ few MB)
    python tools/build_ext_pack.py --era compact --start-year 1800 --end-year 2400
        # compact refit of the local DE442 file for a year range (small download)
    python tools/build_ext_pack.py --era compact --quant --max-n 64 --tol-scale 10 \
        --record-days 1:128,301:64,2:256,3:256,4:512,10:512 \
        --out data/packs/fepk-de442-mobile-1800-2400-v1.bin          # mobile planets (1.74 MB)
    python tools/build_ext_pack.py --era compact --quant --max-n 64 \
        --source-pack data/packs/fepk-minor-v2.bin --record-days \
        2000001:1024,2000002:1024,2000003:1024,2000004:1024,2002060:2048,2005145:2048,2000433:512 \
        --out data/packs/fepk-minor-mobile-1800-2400-v1.bin         # mobile minor bodies (0.29 MB)

Outputs `fepk-de441-ext-past-v1.bin` / `fepk-de441-ext-future-v1.bin` in
data/packs/, separate from the main pack; `Pack::open_many` chains their series
with the main pack's (contiguous at the seams, which are the main pack's own
per-pair coverage limits).

Why refit: a bit-identical copy of DE441 over 30,000 years is ~3 GB. These
packs hold our own Chebyshev fits to DE441, sized by an explicit error budget:

* every record is the Chebyshev interpolant at N Chebyshev-Gauss-Lobatto nodes
  (end points included, so a record matches DE441 - and its neighbours - at
  its boundaries);
* per body, the record length L is fixed (chosen from a size/accuracy study,
  see docs/EXTENDED_RANGE_REPORT.md) and N is chosen *per chunk* of 65,536 days
  (~179 years, one FEPK series): the smallest N whose fit passes the check;
* coefficients of degree >= n_f64 are stored as float32 (FEPK v2 mixed
  precision). n_f64 is the smallest value with the rounding bound
  sum_{k>=n_f64} max|c_k| * 2^-24 <= tol/4 for the chunk;
* check (held out, every record, after float32 rounding): the N-1 midpoints
  between adjacent nodes (where interpolation error peaks) plus 3N uniform
  random epochs; accept if max |r_fit - r_DE441| <= 0.8 tol;
* Earth wrt EMB is not stored: a v2 "alias" series defines it as
  -1/EMRAT x (Moon wrt EMB) (DE441: identical to 5e-13 km in its coefficients).

DE441 is streamed by HTTP Range requests (tools/spk.py StreamRangeFile): one
body x one chunk at a time, held in memory only, prefetched one chunk ahead.
Nothing of the source is written to disk. The output is written data-first
(header patched at the end), so no second copy is made.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import shutil
import struct
import sys
import threading
import time

import numpy as np
from numpy.polynomial import chebyshev as C

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spk import SPK, LocalFile, StreamRangeFile  # noqa: E402

BUILDER_VERSION = "ext-1.0.0"
MAGIC = b"FEPK\0\0\0\x01"
FORMAT_VERSION = 2
NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
SSD_DE441 = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de441.bsp"
DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
MAIN_PACK = os.path.join(DATA, "packs", "fepk-de441-de442-v1.bin")
CHUNK_DAYS = 65536
EMRAT_DE441 = 81.300568221497215  # printed in the DE441 kernel comments
# DE442: GMB/GMM - 1 from its GM table (EMRAT is not printed in its comments), as in the main pack
EMRAT_DE442 = 81.30056846634166
DE442_LOCAL = os.path.join(DATA, "sources", "de442.bsp")
MIN_FREE = 2 << 30

# body: (center, record days L, starting N, tolerance km)
# Tolerances: 0.001" geocentric at the body's closest approach to Earth, split
# between the body and the EMB (0.04 km); planets capped at 1 km (0.9 used);
# Moon 0.1 km (= 0.058" at perigee; 0.001" would need 1.7 m, see report).
BODIES = {
    301: (3, 32, 27, 0.1),
    3: (0, 128, 27, 0.04),
    1: (0, 32, 17, 0.3),
    2: (0, 128, 15, 0.15),
    4: (0, 256, 17, 0.2),
    10: (0, 256, 24, 0.6),
    5: (0, 1024, 15, 0.9),
    6: (0, 1024, 9, 0.9),
    7: (0, 1024, 6, 0.9),
    8: (0, 1024, 5, 0.9),
    9: (0, 1024, 5, 0.9),
}
MAX_N = 32  # --max-n

# Minor bodies for --era compact --source-pack (refit of fepk-minor-v2.bin): (center, record days,
# starting N, tolerance km). Budgets ~0.01" at each body's closest approach to Earth.
MINOR_BODIES = {
    2000001: (0, 256, 12, 8.0),     # Ceres
    2000002: (0, 256, 12, 8.0),     # Pallas
    2000003: (0, 256, 12, 8.0),     # Juno
    2000004: (0, 256, 12, 8.0),     # Vesta
    2002060: (0, 512, 10, 40.0),    # Chiron
    2005145: (0, 512, 10, 40.0),    # Pholus
    2136199: (0, 1024, 6, 200.0),   # Eris
    2000433: (10, 64, 12, 1.0),     # Eros (wrt Sun)
}


class PackTruth:
    """Truth callable over one FEPK pair (for refits of our own minor-body pack)."""

    def __init__(self, pack, target, center):
        self.pack, self.target, self.center = pack, target, center

    def __call__(self, jd, vel=False):
        p, v = self.pack.state(self.target, self.center, np.asarray(jd, dtype=np.float64).ravel())
        return (p, v) if vel else p


def free_bytes(path):
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def main_pack_coverage():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from fepk import Fepk
    pk = Fepk(MAIN_PACK)
    cov = {}
    for s in pk.series:
        k = (s["target"], s["center"])
        lo, hi = cov.get(k, (np.inf, -np.inf))
        cov[k] = (min(lo, s["start"]), max(hi, s["end"]))
    return cov


class Truth:
    """DE441 position/velocity from a block of type-2 records (vectorised)."""

    def __init__(self, seg, i0, recs):
        self.t0 = seg.init_jd + i0 * seg.interval_days
        self.iv = seg.interval_days
        self.n = seg.ncoef
        self.co = recs[:, 2:].reshape(-1, 3, self.n)

    def __call__(self, jd, vel=False):
        x = (jd - self.t0) / self.iv
        i = np.clip(np.floor(x).astype(np.int64), 0, len(self.co) - 1)
        u = 2.0 * (x - i) - 1.0
        T = C.chebvander(u, self.n - 1)
        p = np.einsum("mck,mk->mc", self.co[i], T)
        if not vel:
            return p
        D = np.zeros_like(T)
        D[:, 1] = 1.0
        for k in range(2, self.n):
            D[:, k] = 2 * u * D[:, k - 1] + 2 * T[:, k - 1] - D[:, k - 2]
        return p, np.einsum("mck,mk->mc", self.co[i], D) * (2.0 / self.iv)


def lobatto(N):
    return np.cos(np.pi * np.arange(N) / (N - 1))


QSTATS = []
QUANT = False  # --quant: KIND_QUANT records (integers of step q), set in main()


QUANT_C = 0.6  # --quant-c: q = C*tol/sqrt(N) (0 = worst-case step 0.4*tol/N)


def quant_step(tol, N):
    """Integer step (km): rounding adds at most N*q/2 = 0.2 tol per component."""
    if QUANT_C:
        return QUANT_C * tol / np.sqrt(N)
    return 0.4 * tol / N


def fit_chunk(truth, start, L, nrec, N, tol, rng):
    """Fit nrec records of L days from `start` with N coefficients. Returns
    (ok, coef (nrec,3,N) f64 after rounding, n64, stats)."""
    s = lobatto(N)
    starts = start + L * np.arange(nrec)
    vals = truth((starts[:, None] + (s[None, :] + 1.0) * (0.5 * L)).ravel()).reshape(nrec, N, 3)
    V = C.chebvander(s, N - 1)
    coef = np.linalg.solve(V, vals.transpose(1, 0, 2).reshape(N, -1)).reshape(N, nrec, 3)
    mag = np.abs(coef).max(axis=(1, 2))
    n64 = N
    for k in range(1, N + 1):
        if mag[k:].sum() * 2.0**-24 <= 0.25 * tol:
            n64 = k
            break
    if QUANT:
        q = quant_step(tol, N)
        coef = np.round(coef / q) * q
        n64 = N
    else:
        coef[n64:] = coef[n64:].astype(np.float32).astype(np.float64)
    mids = np.cos(np.pi * (np.arange(N - 1) + 0.5) / (N - 1))
    u = np.concatenate([mids[None, :].repeat(nrec, 0), rng.uniform(-1, 1, (nrec, 3 * N))], axis=1)
    jd = starts[:, None] + (u + 1.0) * (0.5 * L)
    tp, tv = truth(jd.ravel(), vel=True)
    T = C.chebvander(u.ravel(), N - 1).reshape(nrec, u.shape[1], N)
    D = np.zeros_like(T)
    D[..., 1] = 1.0
    for k in range(2, N):
        D[..., k] = 2 * u * D[..., k - 1] + 2 * T[..., k - 1] - D[..., k - 2]
    cr = coef.transpose(1, 0, 2)  # (nrec, N, 3)
    fp = np.einsum("ruk,rkc->ruc", T, cr).reshape(-1, 3)
    fv = np.einsum("ruk,rkc->ruc", D, cr).reshape(-1, 3) * (2.0 / L)
    dp = np.linalg.norm(fp - tp, axis=1)
    dv = np.linalg.norm(fv - tv, axis=1)
    stats = dict(max_km=float(dp.max()), rms_km=float(np.sqrt((dp**2).mean())),
                 max_vel_km_d=float(dv.max()), n_check=int(dp.size), r_min_km=float(np.linalg.norm(tp, axis=1).min()))
    return dp.max() <= 0.8 * tol, cr, n64, stats


def encode(cr, n64):
    """(nrec, N, 3) -> bytes in FEPK v2 record layout."""
    nrec, N, _ = cr.shape
    c = cr.transpose(0, 2, 1)  # (nrec, 3, N)
    hi = np.ascontiguousarray(c[:, :, :n64]).astype("<f8").reshape(nrec, -1)
    lo = c[:, :, n64:].astype("<f4").reshape(nrec, -1)
    pad = (-lo.shape[1]) % 2
    if pad:
        lo = np.concatenate([lo, np.zeros((nrec, pad), "<f4")], axis=1)
    return np.concatenate([hi.view(np.uint8), lo.view(np.uint8)], axis=1).tobytes()


def encode_quant(cr, tol):
    """(nrec, N, 3) multiples of q -> KIND_QUANT block (header, bit-packed records, padding)."""
    nrec, N, _ = cr.shape
    q = quant_step(tol, N)
    ints = np.rint(cr.transpose(0, 2, 1) / q).astype(np.int64)  # (nrec, 3, N)
    hi, lo = ints.max(axis=0), ints.min(axis=0)
    widths = np.zeros((3, N), dtype=np.int64)
    for a in range(3):
        for k in range(N):
            m = max(int(hi[a, k]), -int(lo[a, k]) - 1, 0)
            widths[a, k] = 0 if hi[a, k] == 0 and lo[a, k] == 0 else m.bit_length() + 1
    w = widths.ravel()
    rec_bits = int(w.sum())
    if os.environ.get("FE_QUANT_STATS"):
        # Rice-code estimate: per (component, degree) parameter from the mean |v|
        z = np.abs(ints).reshape(nrec, -1).astype(float)
        mean = z.mean(axis=0)
        kpar = np.maximum(np.floor(np.log2(np.maximum(mean, 1.0) * 0.69)), 0)
        rice = (np.floor(z / 2.0 ** kpar) + 1 + kpar + 1).sum()  # unary + stop + k bits + sign
        QSTATS.append((nrec * rec_bits, float(rice) + 32 * nrec))
    out = bytearray(struct.pack("<dII", q, rec_bits, N))
    out += bytes(w.astype(np.uint8)) + bytes((-3 * N) % 8)
    acc, accbits = 0, 0
    flat = ints.reshape(nrec, -1)
    for r in range(nrec):
        rec, off = 0, 0
        for v, b in zip(flat[r].tolist(), w.tolist()):
            if b:
                rec |= (v & ((1 << b) - 1)) << off
                off += b
        acc |= rec << accbits
        accbits += rec_bits
        if accbits >= 1 << 15:
            nb = accbits // 8
            out += (acc & ((1 << (8 * nb)) - 1)).to_bytes(nb, "little")
            acc >>= 8 * nb
            accbits -= 8 * nb
    nb = (accbits + 7) // 8
    out += acc.to_bytes(nb, "little") + bytes(16)
    out += bytes((-len(out)) % 8)
    return bytes(out)


def record_bytes(N, n64):
    return 24 * n64 + 8 * (-(-3 * (N - n64) // 2))


class Prefetcher:
    """Fetch the next chunk's DE441 records in a background thread."""

    def __init__(self, seg, jobs):
        self.seg, self.jobs, self.k = seg, jobs, 0
        self.next = None
        self._start()

    def _start(self):
        if self.k >= len(self.jobs):
            self.thread = None
            return
        c0, c1 = self.jobs[self.k][:2]
        seg = self.seg
        i0 = max(0, int(np.floor((c0 - seg.init_jd) / seg.interval_days)))
        i1 = min(seg.n, int(np.ceil((c1 - seg.init_jd) / seg.interval_days)))
        box = {}

        def run():
            box["v"] = (i0, seg.records(i0, i1))
        self.box = box
        self.thread = threading.Thread(target=run, daemon=True)
        self.thread.start()

    def get(self):
        self.thread.join()
        if "v" not in self.box:
            raise IOError("prefetch failed")
        v = self.box["v"]
        self.k += 1
        self._start()
        return v


def year_jd(year):
    """JD (TDB) of 1 January 00:00 of a Gregorian year."""
    a = (14 - 1) // 12
    y, m = year + 4800 - a, 1 + 12 * a - 3
    return (1 + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045) - 0.5


def build(args):
    era = args.era
    source_pack = None
    if era == "compact":
        url = args.source_pack or args.de442
        src = LocalFile(url)
        src.bytes_fetched, src.last_modified = 0, None
        compact = (year_jd(args.start_year), year_jd(args.end_year))
        if args.source_pack:
            from fepk import Fepk
            source_pack = Fepk(args.source_pack)
            BODIES.clear()
            BODIES.update(MINOR_BODIES)
    else:
        # ssd.jpl.nasa.gov serves the same DE441 as one file (both parts' segments) and was
        # several times faster than naif.jpl.nasa.gov here; --source naif uses the parts.
        if args.source == "ssd":
            url = SSD_DE441
        else:
            url = NAIF + ("de441_part-1.bsp" if era == "past" else "de441_part-2.bsp")
        src = StreamRangeFile(url, parallel=args.parallel)
    spk = None if source_pack else SPK(src)
    cov = None if era == "compact" else main_pack_coverage()
    default_name = (f"fepk-minor-compact-{args.start_year}-{args.end_year}-v1.bin" if args.source_pack
                    else f"fepk-de442-compact-{args.start_year}-{args.end_year}-v1.bin" if era == "compact"
                    else f"fepk-de441-ext-{era}-v1.bin")
    out = args.out or os.path.join(DATA, "packs", default_name)
    rng = np.random.default_rng(args.seed)
    bodies = [int(b) for b in args.bodies.split(",")] if args.bodies else list(BODIES)

    parts = out + ".parts"
    os.makedirs(parts, exist_ok=True)
    t_start = time.time()
    for b in bodies:
        if free_bytes(os.path.dirname(out)) < MIN_FREE:
            raise SystemExit(f"free disk below 2 GB ({free_bytes(os.path.dirname(out))/1e9:.2f} GB): stopping")
        ctr, L, N0, tol = BODIES[b]
        if source_pack:
            lo_main, hi_main = compact
            nrec_total = int(np.floor((hi_main - lo_main) / L))
            e = lo_main + nrec_total * L
            edges = [lo_main + k * CHUNK_DAYS for k in range(0, int((e - lo_main) // CHUNK_DAYS) + 1)]
            if edges[-1] < e:
                edges.append(e)
            jobs = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
            seg = None
        elif era == "compact":
            lo_main, hi_main = compact
            probe = lo_main + 1.0
        else:
            lo_main, hi_main = cov[(b, ctr)]
            probe = lo_main - 1.0 if era == "past" else hi_main + 1.0
        segs = [None] if source_pack else [
            g for g in spk.segments if g.target == b and g.center == ctr and g.start_jd <= probe <= g.end_jd]
        if len(segs) != 1:
            raise SystemExit(f"expected one source segment for {b}/{ctr} covering JD {probe}")
        seg = segs[0]
        if source_pack:
            pass  # jobs computed above
        elif era == "compact":
            if not (seg.start_jd <= lo_main and hi_main <= seg.end_jd):
                raise SystemExit(f"{b}/{ctr}: DE442 covers JD {seg.start_jd}..{seg.end_jd} only")
            nrec_total = int(np.floor((hi_main - lo_main) / L))
            e = lo_main + nrec_total * L
            edges = [lo_main + k * CHUNK_DAYS for k in range(0, int((e - lo_main) // CHUNK_DAYS) + 1)]
            if edges[-1] < e:
                edges.append(e)
            jobs = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
        elif era == "past":
            seam = lo_main
            nrec_total = int(np.floor((seam - seg.start_jd) / L))
            a = seam - nrec_total * L
            edges = [a] + [seam - k * CHUNK_DAYS for k in range(int((seam - a) // CHUNK_DAYS), 0, -1)
                           if seam - k * CHUNK_DAYS > a] + [seam]
            # process from the seam outwards (most used first)
            jobs = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)][::-1]
        else:
            seam = hi_main
            nrec_total = int(np.floor((seg.end_jd - seam) / L))
            e = seam + nrec_total * L
            edges = [seam + k * CHUNK_DAYS for k in range(0, int((e - seam) // CHUNK_DAYS) + 1)]
            if edges[-1] < e:
                edges.append(e)
            jobs = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
        if args.max_chunks:
            jobs = jobs[:args.max_chunks]
        # resume: completed chunks of this body
        jpath, bpath = os.path.join(parts, f"{b}.json"), os.path.join(parts, f"{b}.bin")
        state = json.load(open(jpath)) if os.path.exists(jpath) else {"rows": [], "stats": []}
        done = {r["start_jd_tdb"] for r in state["rows"]}
        size_ok = sum(r.get("nbytes", r["n_records"] * record_bytes(r["n_coeff"], r["n_coeff_f64"]))
                      for r in state["rows"])
        with open(bpath, "ab") as g:
            g.truncate(size_ok)
        todo = [j for j in jobs if j[0] not in done]
        if not todo:
            print(f"body {b}: complete ({len(state['rows'])} chunks)", flush=True)
            continue
        pf = None if source_pack else Prefetcher(seg, todo)
        N = max([N0] + [r["n_coeff"] for r in state["rows"][-1:]])
        for k, (c0, c1) in enumerate(todo):
            if free_bytes(os.path.dirname(out)) < MIN_FREE:
                raise SystemExit("free disk below 2 GB: stopping (resumable)")
            if source_pack:
                recs, truth = None, PackTruth(source_pack, b, ctr)
            else:
                i0, recs = pf.get()
                truth = Truth(seg, i0, recs)
            nrec = int(round((c1 - c0) / L))
            assert abs(c0 + nrec * L - c1) < 1e-9
            # smallest passing N, starting two below the previous chunk's
            N = max(5, N - 2)
            while True:
                ok, cr, n64, st = fit_chunk(truth, c0, L, nrec, N, tol, rng)
                if ok:
                    break
                N += 1
                if N > MAX_N:
                    raise SystemExit(f"body {b}: no N <= {MAX_N} reaches {tol} km in chunk {c0}..{c1}")
            blob = encode_quant(cr, tol) if QUANT else encode(cr, n64)
            assert QUANT or len(blob) == nrec * record_bytes(N, n64)
            with open(bpath, "ab") as g:
                off = g.tell()
                g.write(blob)
            state["rows"].append(dict(target=b, center=ctr, n_coeff=N, n_records=nrec, start_jd_tdb=c0,
                                      interval_days=float(L), end_jd_tdb=c0 + nrec * L, rel_offset=off,
                                      n_coeff_f64=0 if QUANT else n64, kind=2 if QUANT else 0,
                                      alias_target=0, nbytes=len(blob)))
            st.update(jd=[c0, c1], N=N, n_f64=n64)
            state["stats"].append(st)
            with open(jpath + ".tmp", "w") as g:
                json.dump(state, g)
            os.replace(jpath + ".tmp", jpath)
            del recs, truth
            el = time.time() - t_start
            print(f"  {b:>3}/{ctr} chunk {k+1}/{len(todo)} JD {c0:.1f}..{c1:.1f} N={N} n64={n64} "
                  f"max {st['max_km']:.3g} km (tol {tol}) streamed {src.bytes_fetched/1e6:.0f} MB "
                  f"{el:.0f} s", flush=True)
    if (args.max_chunks and not args.assemble_partial) or args.no_assemble:
        print("parts written; not assembling (--max-chunks/--no-assemble)")
        return
    assemble(args, era, out, parts, bodies, url, src, t_start)


def assemble(args, era, out, parts, bodies, url, src, t_start):
    tmp = out + ".tmp"
    table, stats_all = [], {}
    f = open(tmp, "wb")
    f.write(bytes(64))
    sizes = {}
    for b in bodies:
        st = json.load(open(os.path.join(parts, f"{b}.json")))
        sizes[b] = sum(r.get("nbytes", r["n_records"] * record_bytes(r["n_coeff"], r["n_coeff_f64"]))
                       for r in st["rows"])
        base = f.tell()
        with open(os.path.join(parts, f"{b}.bin"), "rb") as g:
            shutil.copyfileobj(g, f, 1 << 20)
        for r in st["rows"]:
            r = dict(r)
            r.pop("nbytes", None)
            r["coeff_offset"] = base + r.pop("rel_offset")
            table.append(r)
        bs = st["stats"]
        ctr, L, _n0, tol = BODIES[b]
        rows = [r for r in table if r["target"] == b]
        stats_all[f"{b}/{ctr}"] = dict(
            tol_km=tol, record_days=L, accept_km=0.8 * tol,
            max_km=max(s["max_km"] for s in bs), max_vel_km_d=max(s["max_vel_km_d"] for s in bs),
            rms_km=float(np.sqrt(np.mean([s["rms_km"] ** 2 for s in bs]))),
            n_checked=sum(s["n_check"] for s in bs),
            N_range=[min(s["N"] for s in bs), max(s["N"] for s in bs)],
            n_f64_range=[min(s["n_f64"] for s in bs), max(s["n_f64"] for s in bs)],
            r_min_km=min(s["r_min_km"] for s in bs),
            coverage_jd_tdb=[min(r["start_jd_tdb"] for r in rows), max(r["end_jd_tdb"] for r in rows)],
            n_series=len(rows),
            bytes=sizes[b],
            worst_chunks=sorted(bs, key=lambda s: -s["max_km"])[:3])
        while f.tell() % 8:
            f.write(b"\0")
    emrat = EMRAT_DE442 if era == "compact" else EMRAT_DE441
    # Earth wrt EMB: alias of the Moon, scaled by -1/EMRAT, over the Moon's span
    if 301 in bodies:
        moon = sorted([r for r in table if r["target"] == 301], key=lambda r: r["start_jd_tdb"])
        lo, hi = moon[0]["start_jd_tdb"], moon[-1]["end_jd_tdb"]
        off = f.tell()
        f.write(struct.pack("<d", -1.0 / emrat))
        table.append(dict(target=399, center=3, n_coeff=0, n_records=1, start_jd_tdb=lo,
                          interval_days=hi - lo, end_jd_tdb=hi, coeff_offset=off, n_coeff_f64=0,
                          kind=1, alias_target=301))
    table.sort(key=lambda r: (r["target"], r["center"], r["start_jd_tdb"]))
    data_end = f.tell()

    pair_cov = {}
    for r in table:
        k = f"{r['target']}/{r['center']}"
        a, b2 = pair_cov.get(k, (np.inf, -np.inf))
        pair_cov[k] = (min(a, r["start_jd_tdb"]), max(b2, r["end_jd_tdb"]))
    manifest = {
        "name": (f"fepk-de442-compact-{args.start_year}-{args.end_year}" if era == "compact"
                 else f"fepk-de441-ext-{era}"),
        "format_version": FORMAT_VERSION,
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builder": {"tool": "tools/build_ext_pack.py", "version": BUILDER_VERSION, "argv": sys.argv[1:]},
        "kind": "refit",
        "coverage_jd_tdb": [max(v[0] for v in pair_cov.values()), min(v[1] for v in pair_cov.values())],
        "pair_coverage_jd_tdb": {k: list(v) for k, v in pair_cov.items()},
        "sources": [{"name": os.path.basename(url), "url": url, "size_bytes": src.size,
                     "last_modified": src.last_modified, "bytes_streamed": src.bytes_fetched,
                     "md5": None, "md5_note": "streamed by HTTP Range; not stored, md5 not computed"}],
        "fit": {"method": "Chebyshev interpolation at Gauss-Lobatto nodes per record; float32 for degree >= n_coeff_f64",
                "chunk_days": CHUNK_DAYS, "accept": "max held-out 3-D error <= 0.8 tol after float32 rounding",
                "bodies": stats_all},
        "earth": {"alias": "399 wrt 3 = -1/EMRAT x (301 wrt 3)", "EMRAT": emrat},
        "constants": {"EMRAT": emrat, "AU_km": 149597870.7, "CLIGHT_km_s": 299792.458,
                      "source": "DE442 GM table (GMB/GMM - 1)" if era == "compact" else "DE441 kernel comments"},
        "series": [{k: r[k] for k in ("target", "center", "n_coeff", "n_coeff_f64", "n_records",
                                      "start_jd_tdb", "interval_days", "kind")} for r in table],
        "license": f"Derived from JPL {'DE442' if era == 'compact' else 'DE441'} (public): https://naif.jpl.nasa.gov/naif/rules.html",
    }
    man = json.dumps(manifest, indent=1).encode()
    pad = (-data_end) % 8
    f.write(bytes(pad))
    table_off = f.tell()
    for r in table:
        f.write(struct.pack("<iiII", r["target"], r["center"], r["n_coeff"], r["n_records"]))
        f.write(struct.pack("<dddQ", r["start_jd_tdb"], r["interval_days"], r["end_jd_tdb"], r["coeff_offset"]))
        f.write(struct.pack("<IIIi", 0, r["n_coeff_f64"], r["kind"], r["alias_target"]))
    man_off = f.tell()
    f.write(man)
    total = f.tell()
    f.seek(0)
    f.write(MAGIC + struct.pack("<II", FORMAT_VERSION, len(table)))
    f.write(struct.pack("<QQQQ", table_off, man_off, len(man), 64) + bytes(16))
    f.close()
    os.replace(tmp, out)
    if not args.keep_parts:
        shutil.rmtree(parts)
    h = hashlib.sha256()
    with open(out, "rb") as g:
        for blk in iter(lambda: g.read(1 << 20), b""):
            h.update(blk)
    side = dict(manifest, file={"path": os.path.basename(out), "size_bytes": total, "sha256": h.hexdigest()})
    with open(out + ".manifest.json", "w") as g:
        json.dump(side, g, indent=1)
    if QSTATS:
        fixed, rice = (sum(v) / 8e6 for v in zip(*QSTATS))
        print(f"coefficient bits: fixed-width {fixed:.2f} MB, Rice-coded estimate {rice:.2f} MB")
    print(f"wrote {out}: {total} bytes ({total/1e6:.1f} MB), {len(table)} series, streamed "
          f"{src.bytes_fetched/1e6:.0f} MB in {time.time()-t_start:.0f} s")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--era", choices=["past", "future", "compact"], required=True)
    ap.add_argument("--start-year", type=int, default=1800, help="compact: first year")
    ap.add_argument("--end-year", type=int, default=2400, help="compact: end year (exclusive)")
    ap.add_argument("--de442", default=DE442_LOCAL, help="compact: local de442.bsp")
    ap.add_argument("--source-pack", help="compact: refit the minor bodies of this FEPK pack instead")
    ap.add_argument("--record-days", help="override record lengths, e.g. 1:64,301:16 (powers of 2)")
    ap.add_argument("--quant", action="store_true", help="bit-packed integer coefficients (KIND_QUANT)")
    ap.add_argument("--quant-c", type=float, default=0.6, help="--quant step q = C*tol/sqrt(N)")
    ap.add_argument("--max-n", type=int, default=32, help="largest coefficient count per record (<= 64)")
    ap.add_argument("--tol-scale", type=float, default=1.0, help="multiply every error budget")
    ap.add_argument("--bodies", help="comma list (default all)")
    ap.add_argument("--max-chunks", type=int, help="test: only the chunks nearest the seam")
    ap.add_argument("--parallel", type=int, default=12)
    ap.add_argument("--source", choices=["ssd", "naif"], default="ssd")
    ap.add_argument("--seed", type=int, default=441)
    ap.add_argument("--out")
    ap.add_argument("--no-assemble", action="store_true", help="only fit (resumable parts)")
    ap.add_argument("--keep-parts", action="store_true")
    ap.add_argument("--assemble-partial", action="store_true",
                    help="with --max-chunks: assemble a pack of just those chunks (e.g. 2401-2580)")
    args = ap.parse_args()
    global QUANT, QUANT_C, MAX_N
    QUANT, QUANT_C, MAX_N = args.quant, args.quant_c, args.max_n
    for table in (BODIES, MINOR_BODIES):
        for body, (ctr, L, n0, tol) in list(table.items()):
            table[body] = (ctr, L, n0, tol * args.tol_scale)
    for item in (args.record_days or "").split(","):
        if item:
            body, days = (int(v) for v in item.split(":"))
            table = MINOR_BODIES if body in MINOR_BODIES else BODIES
            ctr, _L, n0, tol = table[body]
            table[body] = (ctr, days, n0, tol)
    build(args)


if __name__ == "__main__":
    main()
