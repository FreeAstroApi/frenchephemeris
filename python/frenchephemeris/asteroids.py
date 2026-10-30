# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Numbered asteroids on demand: JPL Horizons small-body SPK -> our Chebyshev fit -> FEPK pack.

    import frenchephemeris as fe
    from frenchephemeris import asteroids
    path = asteroids.fetch(1181)                       # Lilith, 1599-12-11 .. 2401
    eng = fe.open_default()                            # picks up the cache directory
    eng.calc_ut(2451545.0, fe.AST_OFFSET + 1181)

`fetch` downloads the body's SPK (type 21) from the Horizons API, fits it
(barycentric ICRF, Chebyshev-Gauss-Lobatto records, max error <= `tol_km` on
out-of-sample points, record length halved locally where needed), writes
`ast_<number>.bin` (FEPK v1, docs/FORMAT.md) plus a JSON sidecar with the
JPL solution id into the cache directory, and deletes the SPK. The cache
directory is `cache_dir`, else `$FRENCHEPHEMERIS_ASTEROID_CACHE`, else
`~/.cache/frenchephemeris/asteroids`. Network access happens only here, in
Python, never in the Rust core.

This module needs only numpy (it does not import the native extension), so
the offline pack builder (tools/build_asteroids.py) uses the same code.
"""
from __future__ import annotations

import base64
import datetime as _dt
import hashlib
import json
import os
import re
import struct
import subprocess
import tempfile
import urllib.parse
import urllib.request

import numpy as np

HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"
# Horizons small-body ephemerides start at A.D. 1599-12-11 TDB.
DEFAULT_START = "1599-12-11"
DEFAULT_STOP = "2401-04-01"
FALLBACK_SPANS = [("1700-01-01", "2300-01-01"), ("1800-01-01", "2200-01-01"), ("1900-01-01", "2100-01-01")]
J2000 = 2451545.0
SPD = 86400.0
GM_SUN_AU3_D2 = 2.9591220828411951e-04
AU_KM = 149597870.7
MAX_N = 32
AST_OFFSET = 10000
BUILDER = "frenchephemeris.asteroids 1.0.0"


def naif_id(number: int) -> int:
    """NAIF small-body id: 2000000 + n (n < 1e6), else 20000000 + n."""
    return 2_000_000 + number if number < 1_000_000 else 20_000_000 + number


def default_cache_dir() -> str:
    d = os.environ.get("FRENCHEPHEMERIS_ASTEROID_CACHE")
    return d if d else os.path.join(os.path.expanduser("~"), ".cache", "frenchephemeris", "asteroids")


# --------------------------------------------------------------------------- network
def _get(url: str, timeout: float = 600.0) -> bytes:
    """HTTP GET: curl when available (system trust store), else urllib."""
    import shutil

    if shutil.which("curl"):
        r = subprocess.run(["curl", "-sSf", "--connect-timeout", "30", "--max-time", str(int(timeout)), url],
                           capture_output=True)
        if r.returncode == 0:
            return r.stdout
        raise RuntimeError(f"curl failed ({r.returncode}): {r.stderr.decode(errors='replace')[:300]}")
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def horizons_spk(number: int, out_path: str, start: str = DEFAULT_START, stop: str = DEFAULT_STOP) -> dict:
    """Download the Horizons SPK of a numbered asteroid; returns request metadata."""
    params = {"format": "json", "COMMAND": f"'{number};'", "EPHEM_TYPE": "'SPK'", "OBJ_DATA": "'YES'",
              "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'"}
    d = json.loads(_get(HORIZONS_API + "?" + urllib.parse.urlencode(params)))
    if "spk" not in d:
        raise RuntimeError(f"Horizons refused SPK for {number} ({start}..{stop}): {d.get('error')}")
    with open(out_path, "wb") as f:
        f.write(base64.b64decode(d["spk"]))
    res = d.get("result", "")

    def grab(p):
        m = re.search(p, res, re.M)
        return m.group(1).strip() if m else None
    return {"api": HORIZONS_API, "params": params, "spk_file_id": d.get("spk_file_id"),
            "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "target": grab(r"JPL/HORIZONS\s+(.+?)\s{2,}"), "solution_date": grab(r"Soln\.date:\s*(\S+)"),
            "observations": grab(r"# obs:\s*(.+)$"), "soln_ref": grab(r"soln ref\.=\s*([^,\n]+)")}


# --------------------------------------------------------------------------- SPK type 21
class SPK21:
    """Minimal DAF/SPK reader for Horizons small-body files (type 21 segments).
    Same algorithm as tools/spk.py (validated there against Horizons vectors),
    vectorised over epochs."""

    def __init__(self, path: str):
        raw = open(path, "rb").read()
        self.raw = raw
        if not raw[:7] == b"DAF/SPK":
            raise ValueError(f"{path}: not an SPK")
        e = {"LTL-IEEE": "<", "BIG-IEEE": ">"}[raw[88:96].decode()]
        nd, ni = struct.unpack(e + "ii", raw[8:16])
        fward = struct.unpack(e + "i", raw[76:80])[0]
        ss = nd + (ni + 1) // 2
        self.comments = raw[1024:(fward - 1) * 1024].replace(b"\0", b"\n").split(b"\x04")[0].decode("ascii", "replace")
        segs = []
        rec = fward
        while rec:
            o = (rec - 1) * 1024
            nxt, _, nsum = struct.unpack(e + "ddd", raw[o:o + 24])
            for k in range(int(nsum)):
                b = o + 24 + k * ss * 8
                t0, t1 = struct.unpack(e + "2d", raw[b:b + 16])
                tgt, ctr, frm, typ, ia, ib = struct.unpack(e + "6i", raw[b + 16:b + 40])
                segs.append((t0, t1, tgt, ctr, frm, typ, ia, ib))
            rec = int(nxt)
        lines, epochs, t0s = [], [], []
        md_all = set()
        for t0, t1, tgt, ctr, frm, typ, ia, ib in sorted(segs):
            if typ != 21 or ctr != 10 or frm != 1:
                raise ValueError(f"unsupported segment type {typ} centre {ctr} frame {frm}")
            md, n = struct.unpack(e + "2d", raw[(ib - 2) * 8:ib * 8])
            md, n = int(md), int(n)
            dl = 4 * md + 11
            if ia + n * dl + n + n // 100 + 2 - 1 != ib:
                raise ValueError("inconsistent type-21 segment")
            a = np.frombuffer(raw[(ia - 1) * 8:(ia - 1 + n * dl + n) * 8], dtype=e + "f8").astype("<f8")
            lines.append(a[:n * dl].reshape(n, dl))
            epochs.append(a[n * dl:])
            md_all.add(md)
            t0s.append(t0)
            self.target = tgt
        if len(md_all) != 1:
            raise ValueError("mixed MAXDIM")
        self.md = md_all.pop()
        # segments are contiguous in time; concatenate lines in time order (later wins on ties)
        self.lines = np.concatenate(lines)
        self.epochs = np.concatenate(epochs)
        order = np.argsort(self.epochs, kind="stable")
        self.lines, self.epochs = self.lines[order], self.epochs[order]
        self.start_et = min(s[0] for s in segs)
        self.end_et = max(s[1] for s in segs)

    @property
    def coverage_jd(self):
        return J2000 + self.start_et / SPD, J2000 + self.end_et / SPD

    def pos(self, jd) -> np.ndarray:
        """Heliocentric ICRF position (km) at JD TDB (array)."""
        et = (np.asarray(jd, dtype=np.float64).ravel() - J2000) * SPD
        if et.min() < self.start_et - 1e-3 or et.max() > self.end_et + 1e-3:
            raise ValueError("epoch outside the SPK coverage")
        li = np.minimum(np.searchsorted(self.epochs, et, side="left"), len(self.epochs) - 1)
        md = self.md
        L = self.lines[li]
        K = L[:, 7 + 4 * md].astype(int)
        out = np.empty((len(et), 3))
        for k in np.unique(K):
            m = K == k
            out[m] = self._eval(L[m], et[m], int(k))
        return out

    def _eval(self, L, et, kqmax1):
        md = self.md
        tl = L[:, 0]
        g = L[:, 1:1 + md]                      # G[1..md] -> g[:, j-1]
        rv = L[:, 1 + md:7 + md]
        refpos, refvel = rv[:, 0::2], rv[:, 1::2]
        dt = L[:, 7 + md:7 + 4 * md].reshape(-1, 3, md)
        kq = L[:, 8 + 4 * md:11 + 4 * md].astype(int)
        npt = len(et)
        delta = et - tl
        mq2 = kqmax1 - 2
        fc = np.zeros((kqmax1 + 2, npt))
        wc = np.zeros((kqmax1 + 2, npt))
        tp = delta.copy()
        for j in range(1, mq2 + 1):
            fc[j + 1] = tp / g[:, j - 1]
            wc[j] = delta / g[:, j - 1]
            tp = delta + g[:, j - 1]
        w = np.zeros((2 * kqmax1 + 2, npt))
        for j in range(1, kqmax1 + 1):
            w[j] = 1.0 / j
        ks, ks1, jx = kqmax1 - 1, kqmax1 - 2, 0
        while ks >= 2:
            jx += 1
            for j in range(1, jx + 1):
                w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks]
            ks, ks1 = ks1, ks1 - 1
        jj = np.arange(1, md + 1)
        mask = jj[None, None, :] <= kq[:, :, None]   # (npt, 3, md)
        wk = w[1 + ks:1 + ks + md].T if 1 + ks + md <= w.shape[0] else np.pad(w[1 + ks:], ((0, md), (0, 0)))[:md].T
        acc = np.einsum("pcj,pj->pc", np.where(mask, dt, 0.0), wk)
        return refpos + delta[:, None] * (refvel + delta[:, None] * acc)


# --------------------------------------------------------------------------- FEPK
class Fepk:
    """Read-only FEPK v1 (only what the fitter needs: one stored pair)."""

    def __init__(self, path: str):
        self.mm = np.memmap(path, dtype="<f8", mode="r")
        with open(path, "rb") as f:
            raw = f.read(64)
            if raw[:8] != b"FEPK\0\0\0\x01":
                raise ValueError(f"{path}: not an FEPK v1 file")
            n = struct.unpack("<I", raw[12:16])[0]
            toff, moff, mlen, _ = struct.unpack("<QQQQ", raw[16:48])
            f.seek(toff)
            tab = f.read(n * 64)
            f.seek(moff)
            self.manifest = json.loads(f.read(mlen))
        self.series = []
        for k in range(n):
            o = k * 64
            tgt, ctr, nc, nr = struct.unpack("<iiII", tab[o:o + 16])
            st, iv, en, off = struct.unpack("<dddQ", tab[o + 16:o + 48])
            self.series.append(dict(target=tgt, center=ctr, n_coeff=nc, n_records=nr, start=st, interval=iv,
                                    end=en, offset=off))

    def targets(self):
        return sorted({s["target"] for s in self.series})

    def pos(self, target, center, jd):
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        ss = sorted((s for s in self.series if s["target"] == target and s["center"] == center),
                    key=lambda s: s["start"])
        out = np.full((len(jd), 3), np.nan)
        done = np.zeros(len(jd), bool)
        for s in reversed(ss):
            m = (~done) & (jd >= s["start"]) & (jd <= s["end"])
            if not m.any():
                continue
            n, nr = s["n_coeff"], s["n_records"]
            i = np.minimum(((jd[m] - s["start"]) / s["interval"]).astype(np.int64), nr - 1)
            u = 2.0 * (jd[m] - (s["start"] + i * s["interval"])) / s["interval"] - 1.0
            c = np.asarray(self.mm[s["offset"] // 8:s["offset"] // 8 + nr * 3 * n]).reshape(nr, 3, n)[i]
            out[m] = np.einsum("rck,rk->rc", c, _cheb(u, n))
            done |= m
        if not done.all():
            raise ValueError(f"epoch outside {target}/{center}")
        return out


def _cheb(u, n):
    t = np.empty(np.shape(u) + (n,))
    t[..., 0] = 1.0
    t[..., 1] = u
    for k in range(2, n):
        t[..., k] = 2 * u * t[..., k - 1] - t[..., k - 2]
    return t


def write_fepk(path: str, rows: list, manifest: dict) -> int:
    """rows: dicts target, center, start, interval, coef (n_rec, 3, n), source_index. Sorted here."""
    rows = sorted(rows, key=lambda r: (r["target"], r["center"], r["start"]))
    man = json.dumps(manifest, indent=1).encode()
    table_off, nrow = 64, len(rows)
    man_off = table_off + nrow * 64
    data_off = (man_off + len(man) + 7) // 8 * 8
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(b"FEPK\0\0\0\x01" + struct.pack("<II", 1, nrow))
        f.write(struct.pack("<QQQQ", table_off, man_off, len(man), data_off) + bytes(16))
        off = data_off
        for r in rows:
            nr, _, n = r["coef"].shape
            end = r["start"] + nr * r["interval"]
            f.write(struct.pack("<iiII", r["target"], r["center"], n, nr))
            f.write(struct.pack("<dddQI", r["start"], r["interval"], end, off, r.get("source_index", 0)))
            f.write(bytes(12))
            off += nr * 3 * n * 8
        f.write(man)
        f.write(bytes(data_off - man_off - len(man)))
        for r in rows:
            f.write(np.ascontiguousarray(r["coef"], dtype="<f8").tobytes())
        total = f.tell()
    os.replace(tmp, path)
    return total


# --------------------------------------------------------------------------- fit
def _lobatto(n):
    x = np.cos(np.pi * np.arange(n) / (n - 1))
    k = np.arange(n)
    c = np.cos(np.pi * np.outer(k, k) / (n - 1))
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    m = (2.0 / (n - 1)) * c * w[None, :]
    m[0] *= 0.5
    m[-1] *= 0.5
    return x, m


def fit_body(pos_fn, start: float, end: float, interval: float, n: int = 20, tol_km: float = 0.5,
             seed: int = 0, min_interval: float = 0.25):
    """Fit pos_fn (JD -> (m,3) km) on [start, end]. Records of `interval` days
    (last one shorter), halved recursively where the out-of-sample error
    exceeds tol_km. Returns (series list, error stats)."""
    rng = np.random.default_rng(seed)
    nodes, M = _lobatto(n)
    mids = 0.5 * (nodes[:-1] + nodes[1:])

    def fit_grid(t0s, h):
        jd = t0s[:, None] + (nodes[None, :] + 1) * 0.5 * h
        vals = pos_fn(jd.ravel()).reshape(len(t0s), n, 3)
        coef = np.einsum("jk,rkc->rcj", M, vals)
        s = np.concatenate([np.broadcast_to(mids, (len(t0s), n - 1)), rng.uniform(-1, 1, (len(t0s), 4))], axis=1)
        jc = t0s[:, None] + (s + 1) * 0.5 * h
        truth = pos_fn(jc.ravel()).reshape(len(t0s), s.shape[1], 3)
        err = np.linalg.norm(np.einsum("rck,rmk->rmc", coef, _cheb(s, n)) - truth, axis=2)
        return coef, err

    series, errs = [], []
    # pieces: (start, interval, n_records)
    nrec = int((end - start) // interval)
    pieces = []
    if nrec:
        pieces.append((start, interval, nrec))
    t = start + nrec * interval
    if end - t > 1e-9:
        pieces.append((t, end - t, 1))
    while pieces:
        a, h, nr = pieces.pop(0)
        t0s = a + h * np.arange(nr)
        coef, err = fit_grid(t0s, h)
        bad = err.max(axis=1) > tol_km
        if h / 2 < min_interval:
            bad[:] = False
        k = 0
        while k < nr:
            j = k
            while j < nr and bad[j] == bad[k]:
                j += 1
            if bad[k]:
                pieces.append((a + k * h, h / 2, 2 * (j - k)))
            else:
                series.append(dict(start=a + k * h, interval=h, coef=coef[k:j]))
                errs.append(err[k:j].ravel())
            k = j
    series.sort(key=lambda s: s["start"])
    e = np.concatenate(errs)
    stats = {"n_coeff": n, "base_interval_days": interval,
             "min_interval_days": float(min(s["interval"] for s in series)),
             "n_series": len(series), "n_records": int(sum(len(s["coef"]) for s in series)),
             "bytes": int(sum(s["coef"].size * 8 for s in series)),
             "check_points": int(e.size), "max_km": float(e.max()),
             "p99_km": float(np.percentile(e, 99)), "rms_km": float(np.sqrt(np.mean(e ** 2)))}
    return series, stats


def base_interval(a_au: float) -> float:
    """Initial record length from the orbital period scale (refined locally by fit_body)."""
    tau = 58.13 * max(a_au, 0.3) ** 1.5  # days, 1/n for a circular orbit
    p = int(np.floor(np.log2(max(tau, 16.0))))
    return float(2 ** min(max(p, 4), 10))


def fit_spk(spk: SPK21, main: Fepk, tol_km: float = 0.5, start=None, end=None, n: int | None = None):
    """Barycentric fit (Sun from the main pack) of an SPK over [start, end] (JD TDB)."""
    lo, hi = spk.coverage_jd
    sun = [s for s in main.series if s["target"] == 10 and s["center"] == 0]
    lo = max(lo, min(s["start"] for s in sun)) if start is None else max(start, lo)
    hi = min(hi, max(s["end"] for s in sun), 2598072.5) if end is None else min(end, hi)
    # snap the span start to a whole day boundary (.5) inside coverage
    lo = np.ceil(lo - 0.5) + 0.5
    hi = np.floor(hi - 0.5) + 0.5

    def bary(jd):
        return spk.pos(jd) + main.pos(10, 0, jd)
    grid = np.linspace(lo, hi, 2001)
    r = np.linalg.norm(spk.pos(grid), axis=1)
    a_au = 0.5 * (r.min() + r.max()) / AU_KM
    iv0 = base_interval(a_au)
    # Horizons files can be discontinuous at a few line boundaries (e.g. 2101
    # Adonis: 232 km in 2036, at a close planetary approach). Records must not
    # straddle such a jump, so the span is cut there and the pack reproduces it.
    et = spk.epochs[:-1]
    eb = 1e-3  # s
    jump = np.linalg.norm(spk.pos((et - eb) / SPD + J2000) - spk.pos((et + eb) / SPD + J2000), axis=1)
    cut_jd = [float(t / SPD + J2000) for t, j in zip(et, jump) if j > 0.2 * tol_km]
    cut_jd = [t for t in cut_jd if lo < t < hi]
    edges = [lo, *cut_jd, hi]
    best, closest = None, None
    # smallest pack over a small (record length, degree) grid; each trial is a
    # complete fit with local refinement, so every candidate meets tol_km
    for iv in sorted({min(iv0 * f, 1024.0) for f in (0.5, 1, 2, 4, 8, 16)}):
        for nn in ((n,) if n else (8, 12, 16, 20, 24, 28)):
            parts = [fit_body(bary, a, b, iv, n=nn, tol_km=tol_km) for a, b in zip(edges[:-1], edges[1:])]
            ser = [x for p in parts for x in p[0]]
            st = dict(parts[0][1])
            for k in ("n_series", "n_records", "bytes", "check_points"):
                st[k] = int(sum(p[1][k] for p in parts))
            st["max_km"] = float(max(p[1]["max_km"] for p in parts))
            st["p99_km"] = float(max(p[1]["p99_km"] for p in parts))
            st["rms_km"] = float(np.sqrt(sum(p[1]["rms_km"] ** 2 * p[1]["check_points"] for p in parts)
                                         / st["check_points"]))
            st["min_interval_days"] = float(min(p[1]["min_interval_days"] for p in parts))
            if st["max_km"] <= tol_km and (best is None or st["bytes"] < best[1]["bytes"]):
                best = (ser, st)
            if closest is None or st["max_km"] < closest[1]["max_km"]:
                closest = (ser, st)
    if best is None:
        raise RuntimeError(f"fit cannot reach {tol_km} km (best max error {closest[1]['max_km']:.3f} km, "
                           f"records down to {closest[1]['min_interval_days']} d)")
    series, stats = best
    stats.update({"source_discontinuities": [{"jd_tdb": t, "jump_km": float(j)}
                                             for t, j in zip((et / SPD + J2000).tolist(), jump) if j > 0.2 * tol_km
                                             and lo < t / 1 < hi],
                  "helio_dist_au": [float(r.min() / AU_KM), float(r.max() / AU_KM)],
                  "coverage_jd_tdb": [float(lo), float(hi)], "center": 0})
    return series, stats


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def spk_metadata(spk: SPK21) -> dict:
    c = spk.comments

    def grab(p):
        m = re.search(p, c, re.M)
        return m.group(1).strip() if m else None
    return {"target_line": grab(r"Target body\s*:\s*(.+)$"), "solution_id": grab(r"Horizons_SPK:(JPL#\d+)"),
            "start_time": grab(r"Start time\s*:\s*(.+)$"), "stop_time": grab(r"Stop  time\s*:\s*(.+)$"),
            "small_body_perturbers": grab(r"Small-body perts:\s*(.+)$"),
            "execution_date_pasadena": grab(r"Execution date\s*:\s*(.+)$"), "spk_target_id": spk.target}


def default_main_pack() -> str:
    env = os.environ.get("FRENCHEPHEMERIS_PACK")
    if env:
        return env
    from pathlib import Path
    p = Path(__file__).resolve().parents[2] / "data" / "packs" / "fepk-de441-de442-v1.bin"
    return str(p)


def fetch(mpc_number: int, start: str | None = None, end: str | None = None, cache_dir: str | None = None,
          main_pack: str | None = None, tol_km: float = 0.5, force: bool = False) -> str:
    """Fetch, fit and cache a numbered asteroid; returns the pack path.

    `start`/`end`: Horizons date strings (default 1599-12-11 .. 2401-04-01, the
    longest span Horizons serves; narrower fallbacks are tried if refused).
    The pack holds NAIF id 2000000 + n (wrt the SSB); the engine id is 10000 + n."""
    n = int(mpc_number)
    if n < 1:
        raise ValueError("MPC number must be >= 1")
    cache = cache_dir or default_cache_dir()
    os.makedirs(cache, exist_ok=True)
    out = os.path.join(cache, f"ast_{n}.bin")
    if os.path.exists(out) and not force:
        return out
    main = Fepk(main_pack or default_main_pack())
    spans = [(start or DEFAULT_START, end or DEFAULT_STOP)] + ([] if start or end else FALLBACK_SPANS)
    with tempfile.TemporaryDirectory(dir=cache) as td:
        sp = os.path.join(td, f"{n}.bsp")
        req, err = None, None
        for a, b in spans:
            try:
                req = horizons_spk(n, sp, a, b)
                break
            except RuntimeError as e:
                err = e
        if req is None:
            raise err
        spk = SPK21(sp)
        meta = spk_metadata(spk)
        meta.update({"sha256": _sha256(sp), "size_bytes": os.path.getsize(sp)})
        series, stats = fit_spk(spk, main, tol_km=tol_km)
        # reference states straight from the JPL SPK (heliocentric ICRF km) for tests
        lo, hi = stats["coverage_jd_tdb"]
        jd = np.r_[lo, hi, np.round(np.random.default_rng(n).uniform(lo, hi, 6), 6)]
        samples = [[float(t), *map(float, p)] for t, p in zip(jd, spk.pos(jd))]
    tgt = naif_id(n)
    manifest = pack_manifest(f"ast_{n}", {n: {"request": req, "spk": meta, "fit": stats,
                                              "jpl_helio_samples_km": samples}}, tol_km, main)
    write_fepk(out, [dict(target=tgt, center=0, start=s["start"], interval=s["interval"], coef=s["coef"])
                     for s in series], manifest)
    with open(out + ".json", "w") as f:
        json.dump(manifest, f, indent=1)
    return out


def pack_manifest(name: str, bodies: dict, tol_km: float, main: Fepk) -> dict:
    return {
        "name": name, "format_version": 1, "builder": BUILDER,
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "kind": "own Chebyshev-Lobatto fit of JPL Horizons small-body SPKs (barycentric, centre 0)",
        "tolerance_km": tol_km,
        "barycentric_sun_from": main.manifest.get("name"),
        "ids": "target = NAIF 2000000 + MPC number; engine/Swiss id = 10000 + MPC number",
        "license": ("JPL Solar System Dynamics data (NASA/JPL-Caltech), credit JPL Horizons; "
                    "files are informal releases (see the Horizons SPK comments)"),
        "bodies": {str(k): v for k, v in bodies.items()},
    }


def pack_numbers(path: str) -> set:
    """MPC numbers of the small bodies stored in an FEPK file."""
    f = Fepk(path)
    out = set()
    for t in f.targets():
        if 2_000_000 <= t < 3_000_000:
            out.add(t - 2_000_000)
        elif t >= 20_000_000:
            out.add(t - 20_000_000)
    return out


def cached_packs(cache_dir: str | None = None, exclude: set | None = None) -> list:
    """Cached on-demand packs, skipping bodies already provided by other packs."""
    d = cache_dir or default_cache_dir()
    if not os.path.isdir(d):
        return []
    out = []
    for name in sorted(os.listdir(d)):
        m = re.fullmatch(r"ast_(\d+)\.bin", name)
        if m and not (exclude and int(m.group(1)) in exclude):
            out.append(os.path.join(d, name))
    return out
