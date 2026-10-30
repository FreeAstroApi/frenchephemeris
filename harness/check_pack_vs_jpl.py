#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Check an FEPK pack against direct JPL evaluation (jplephem) over its full range.

1. Bit-identity: every coefficient of every series is compared with the source
   SPK record words read by jplephem's DAF reader (DE442 local file; DE441 via a
   sparse local copy of the HTTP-range cache, see tools/spk.py).
2. Dense random evaluation: for each stored pair, N random epochs over the pair
   coverage, evaluated (a) with a small numpy FEPK evaluator in this script and
   (b) optionally with the Rust reader (`--rust`, runs the `pack_eval` example),
   compared with jplephem position/velocity of the source used at that epoch.
   With `--rust`, chained pairs (geocentric Sun/Moon/planets...) are checked too.

Usage:
  python harness/check_pack_vs_jpl.py [--pack PATH] [--n 20000] [--rust]
"""
import argparse
import json
import os
import struct
import subprocess
import sys
import tempfile

import numpy as np
from jplephem.spk import SPK as JplSPK

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "tools"))
from spk import HttpRangeFile  # noqa: E402

DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
DE441_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp"
PATH = {0: [], 3: [(3, 0)], 301: [(301, 3), (3, 0)], 399: [(399, 3), (3, 0)]}
for b in (1, 2, 4, 5, 6, 7, 8, 9, 10):
    PATH[b] = [(b, 0)]
CHAINED = [(10, 399), (301, 399), (4, 399), (5, 399), (1, 399), (399, 0), (301, 0), (9, 301)]


class FepkPy:
    """Minimal numpy reader of FEPK v1 (independent of the Rust code)."""

    def __init__(self, path):
        self.buf = np.memmap(path, dtype=np.uint8, mode="r")
        b = self.buf
        assert bytes(b[:8]) == b"FEPK\0\0\0\x01"
        ver, n = struct.unpack("<II", bytes(b[8:16]))
        toff, moff, mlen, doff = struct.unpack("<QQQQ", bytes(b[16:48]))
        self.manifest = json.loads(bytes(b[moff:moff + mlen]))
        self.series = []
        for k in range(n):
            o = toff + 64 * k
            t, c, nc, nr = struct.unpack("<iiII", bytes(b[o:o + 16]))
            st, iv, en, co, si = struct.unpack("<dddQI", bytes(b[o + 16:o + 52]))
            coeffs = np.frombuffer(b, dtype="<f8", count=nr * 3 * nc, offset=co).reshape(nr, 3, nc)
            self.series.append(dict(target=t, center=c, n_coeff=nc, n_records=nr, start=st,
                                    interval=iv, end=en, coeffs=coeffs, source_index=si))

    def pair_series(self, t, c):
        return [s for s in self.series if s["target"] == t and s["center"] == c]

    def state(self, t, c, jd):
        jd = np.asarray(jd, dtype=float)
        pos = np.full((len(jd), 3), np.nan)
        vel = np.full((len(jd), 3), np.nan)
        for s in self.pair_series(t, c):
            m = (jd >= s["start"]) & (jd <= s["end"])
            if not m.any():
                continue
            x = jd[m]
            i = np.minimum(((x - s["start"]) / s["interval"]).astype(np.int64), s["n_records"] - 1)
            u = 2.0 * (x - (s["start"] + i * s["interval"])) / s["interval"] - 1.0
            n = s["n_coeff"]
            T = np.zeros((n, len(x)))
            D = np.zeros((n, len(x)))
            T[0] = 1.0
            T[1] = u
            D[1] = 1.0
            for k in range(2, n):
                T[k] = 2 * u * T[k - 1] - T[k - 2]
                D[k] = 2 * u * D[k - 1] + 2 * T[k - 1] - D[k - 2]
            C = s["coeffs"][i]  # (m, 3, n)
            pos[m] = np.einsum("mck,km->mc", C, T)
            vel[m] = np.einsum("mck,km->mc", C, D) * (2.0 / s["interval"])
        return pos, vel


def sparse_de441():
    cache = os.path.join(DATA, "sources", "de441_cache")
    out = os.path.join(DATA, "sources", "de441_part-1.sparse.bsp")
    hf = HttpRangeFile(DE441_URL, cache)
    newest = max(os.path.getmtime(os.path.join(cache, n)) for n in os.listdir(cache))
    if not os.path.exists(out) or os.path.getmtime(out) < newest:
        hf.write_sparse_copy(out)
    return out


def jpl_state(kernels, seam, t, c, jd):
    """Direct jplephem evaluation, choosing DE441 before the seam, DE442 after."""
    p = np.zeros((len(jd), 3))
    v = np.zeros((len(jd), 3))
    pt, pc = PATH[t], PATH[c]
    while pt and pc and pt[-1] == pc[-1]:
        pt, pc = pt[:-1], pc[:-1]
    for name, mask in (("de441", jd < seam), ("de442", jd >= seam)):
        if not mask.any():
            continue
        segs = kernels[name]
        for sgn, path in ((1, pt), (-1, pc)):
            for key in path:
                # split whole/fraction days so jplephem keeps full time precision
                # (a single JD float costs ~2e-6 s in its seconds-offset arithmetic)
                w = np.floor(jd[mask])
                pos, vel = segs[key].compute_and_differentiate(w, jd[mask] - w)
                p[mask] += sgn * pos.T
                v[mask] += sgn * vel.T
    return p, v


def run_rust(pack_path, t, c, jd):
    with tempfile.TemporaryDirectory() as d:
        fin, fout = os.path.join(d, "in.bin"), os.path.join(d, "out.bin")
        np.asarray(jd, dtype="<f8").tofile(fin)
        subprocess.run(["cargo", "run", "-q", "--release", "-p", "fe-core", "--example", "pack_eval", "--",
                        pack_path, str(t), str(c), fin, fout], check=True, cwd=REPO)
        out = np.fromfile(fout, dtype="<f8").reshape(len(jd), 6)
    return out[:, :3], out[:, 3:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default=os.path.join(DATA, "packs", "fepk-de441-de442-v1.bin"))
    ap.add_argument("--n", type=int, default=20000, help="random epochs per pair")
    ap.add_argument("--rust", action="store_true", help="also evaluate through the Rust reader")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    pk = FepkPy(args.pack)
    names = [s["name"] for s in pk.manifest["sources"]]
    kernels, jpl_files = {}, {}
    jpl_files["de442"] = JplSPK.open(os.path.join(DATA, "sources", "de442.bsp"))
    if "de441_part-1.bsp" in names:
        jpl_files["de441"] = JplSPK.open(sparse_de441())
    for k, f in jpl_files.items():
        kernels[k] = {(s.target, s.center): s for s in f.segments}
    seam = min(s["start"] for s in pk.series if names[s["source_index"]] == "de442.bsp")
    print(f"pack {os.path.basename(args.pack)}  sources {names}  DE442 from JD {seam}")

    # ---- 1. bit identity of coefficients ----
    total = 0
    for s, row in zip(pk.series, pk.manifest["series"]):
        src = "de441" if names[s["source_index"]] == "de441_part-1.bsp" else "de442"
        seg = kernels[src][(s["target"], s["center"])]
        init, intlen, rsize, n = seg.daf.read_array(seg.end_i - 3, seg.end_i)
        raw = seg.daf.map_array(seg.start_i, seg.end_i - 4).reshape(int(n), int(rsize))
        i0 = row["first_source_record"]
        ref = raw[i0:i0 + s["n_records"], 2:].reshape(s["n_records"], 3, s["n_coeff"])
        start = 2451545.0 + (init + i0 * intlen) / 86400.0
        assert start == s["start"] and intlen / 86400.0 == s["interval"], "grid mismatch"
        same = np.array_equal(ref.view(np.uint64), np.ascontiguousarray(s["coeffs"]).view(np.uint64))
        if not same:
            raise SystemExit(f"coefficients differ for {s['target']}/{s['center']} from {src}")
        total += s["n_records"]
    print(f"bit-identity: OK, {total} records in {len(pk.series)} series identical to source SPK words")

    # ---- 2. dense random evaluation ----
    pairs = sorted({(s["target"], s["center"]) for s in pk.series})
    summary = {}
    print(f"{'pair':>8} {'src':>6} {'N':>7} {'max|dp| km':>11} {'max rel':>9} {'max|dv| km/d':>12}"
          + ("  rust max|dp|  rust max|dv|" if args.rust else ""))
    for t, c in pairs + (CHAINED if args.rust else []):
        if (t, c) in pairs:
            lo = min(s["start"] for s in pk.pair_series(t, c))
            hi = max(s["end"] for s in pk.pair_series(t, c))
        else:
            ch = [k for k in PATH[t] + PATH[c]]
            lo = max(min(s["start"] for s in pk.pair_series(*k)) for k in ch)
            hi = min(max(s["end"] for s in pk.pair_series(*k)) for k in ch)
        jd = np.sort(np.concatenate([rng.uniform(lo, hi, args.n), [lo, hi, seam, np.nextafter(seam, 0)]]))
        jd = jd[(jd >= lo) & (jd <= hi)]
        rp, rv = jpl_state(kernels, seam, t, c, jd)
        if args.rust:
            up, uv = run_rust(args.pack, t, c, jd)
        for src, m in (("de441", jd < seam), ("de442", jd >= seam)):
            if not m.any():
                continue
            line = f"{t:>4}/{c:<3} {src:>6} {m.sum():>7}"
            rec = {}
            if (t, c) in pairs:
                pp, pv = pk.state(t, c, jd[m])
                dp = np.linalg.norm(pp - rp[m], axis=1)
                dv = np.linalg.norm(pv - rv[m], axis=1)
                rel = dp / np.linalg.norm(rp[m], axis=1)
                line += f" {dp.max():11.3e} {rel.max():9.2e} {dv.max():12.3e}"
                rec.update(py_dp=float(dp.max()), py_rel=float(rel.max()), py_dv=float(dv.max()))
            else:
                line += f" {'(chained)':>11} {'':>9} {'':>12}"
            if args.rust:
                dp = np.linalg.norm(up[m] - rp[m], axis=1)
                dv = np.linalg.norm(uv[m] - rv[m], axis=1)
                line += f"  {dp.max():12.3e}  {dv.max():12.3e}"
                rec.update(rust_dp=float(dp.max()), rust_dv=float(dv.max()))
            summary[f"{t}/{c}/{src}"] = rec
            print(line, flush=True)
    worst = {k: max(v.get(k, 0) for v in summary.values()) for k in ("py_dp", "py_dv", "rust_dp", "rust_dv")}
    print("overall max:", {k: f"{v:.3e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
