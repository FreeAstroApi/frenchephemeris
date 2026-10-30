#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Build an FEPK v1 coefficient pack (docs/FORMAT.md) from JPL DE441 + DE442.

Full pack (default):
    python tools/build_pack.py
        -> data/packs/fepk-de441-de442-v1.bin covering astronomical year -1
           (JD 1720692.5) .. end of 2400 (JD 2598007.5) plus a margin.
        DE441 (de441_part-1.bsp, read remotely through HTTP Range requests with
        an on-disk chunk cache) supplies records before DE442's start, DE442
        (local de442.bsp) the rest.

Dev pack (DE442 only):
    python tools/build_pack.py --range 1900 2100
        -> data/packs/fepk-de442-1900-2100-v1.bin

Coefficients are copied bit-for-bit from the source type-2 records (the MID and
RADIUS words are dropped; they are implied by the uniform record grid and are
checked against it record by record).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import struct
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spk import SPK, S_PER_DAY, J2000_JD  # noqa: E402

BUILDER_VERSION = "1.0.0"
MAGIC = b"FEPK\0\0\0\x01"
FORMAT_VERSION = 1
HEADER_LEN = 64
SERIES_LEN = 64

NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
DEFAULT_DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))

PAIRS = [(t, 0) for t in range(1, 11)] + [(301, 3), (399, 3)]

FULL_START_JD = 1720692.5  # 1 Jan -1 (Julian calendar), astronomical year numbering
FULL_END_JD = 2598007.5    # 1 Jan 2401 (Gregorian) = end of 2400


def calendar_to_jd(year: int, month: int = 1, day: float = 1.0) -> float:
    """JD at 0h of a civil date: Gregorian from 1582-10-15, Julian before (Meeus 7.1)."""
    y, m = (year - 1, month + 12) if month <= 2 else (year, month)
    gregorian = (year, month, day) >= (1582, 10, 15)
    b = 0
    if gregorian:
        a = y // 100
        b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + day + b - 1524.5


def md5_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def published_md5(sources_dir: str, name: str) -> str | None:
    p = os.path.join(sources_dir, "aa_checksums.txt")
    if not os.path.exists(p):
        return None
    for line in open(p):
        parts = line.split()
        if len(parts) == 2 and parts[1] == name:
            return parts[0]
    return None


def parse_constants(text: str) -> dict:
    """Pull AU, EMRAT, CLIGHT and planetary GMs from SPK comments.

    The 'Mass parameter (GM)' table gives each GM in au^3/day^2 (16 digits) and
    km^3/s^2 (rounded). The 'constants used for integration' block (present in
    DE441's comments, absent from DE442's) prints AU, EMRAT and CLIGHT; when it
    is missing, EMRAT is derived from the au^3/day^2 GMB and GMM columns.
    """
    out: dict = {}
    gm_km, gm_au = {}, {}
    for m in re.finditer(r"^[ \t]*(GM[0-9SMB])[ \t]+(\S+)[ \t]+(\S+)[ \t]+(\S+)[ \t]*$", text, re.M):
        gm_au[m.group(1)] = float(m.group(2))
        gm_km[m.group(1)] = float(m.group(4))
    if gm_km:
        out["GM_km3_s2"] = gm_km
        out["GM_au3_d2"] = gm_au
        out["GM_source"] = "kernel comments, 'Mass parameter (GM)' table"
    for key in ("AU", "EMRAT", "CLIGHT"):
        m = re.search(rf"^{key}[ \t]+(\S+)[ \t]*$", text, re.M)
        if m:
            out[{"AU": "AU_km", "CLIGHT": "CLIGHT_km_s"}.get(key, key)] = float(m.group(1).replace("D", "E"))
    if "EMRAT" not in out and {"GMB", "GMM"} <= gm_au.keys():
        # GMB = GM(Earth) + GM(Moon)  ->  EMRAT = GM(Earth)/GM(Moon)
        out["EMRAT"] = gm_au["GMB"] / gm_au["GMM"] - 1.0
        out["EMRAT_note"] = "derived as GMB/GMM - 1 from the au^3/day^2 GM table (not printed in comments)"
    if "AU_km" not in out and "GMS" in gm_au:
        # GMS[au^3/d^2] * AU^3 / 86400^2 = GMS[km^3/s^2]
        au = (gm_km["GMS"] * S_PER_DAY**2 / gm_au["GMS"]) ** (1.0 / 3.0)
        out["AU_km"] = 149597870.7
        out["AU_note"] = ("not printed in these comments; DE44x use the IAU 2012 value; the two GMS columns "
                          f"imply {au:.4f} km")
    return out


def exact(x: float) -> Fraction:
    return Fraction(x)


def check_grid(seg, i0: int, recs: np.ndarray):
    """Assert MID/RADIUS of records i0.. match the uniform grid (exact in seconds)."""
    idx = np.arange(i0, i0 + len(recs), dtype=np.float64)
    mid = seg.init + (idx + 0.5) * seg.intlen
    if not np.array_equal(recs[:, 1], np.full(len(recs), seg.intlen / 2)):
        raise ValueError(f"{seg}: RADIUS != INTLEN/2 in some record")
    err = np.max(np.abs(recs[:, 0] - mid))
    if err > 1e-6:  # seconds; grid epochs are integers so this should be exactly 0
        raise ValueError(f"{seg}: MID off-grid by {err} s")
    return float(err)


def seam_epoch(s441, s442) -> Fraction:
    """Smallest epoch >= DE442 start lying on both record grids (seconds past J2000)."""
    a0, da = exact(s441.init), exact(s441.intlen)
    b0, db = exact(s442.init), exact(s442.intlen)
    # search a bounded number of DE442 records for one that is also a DE441 boundary
    for j in range(0, 4096):
        e = b0 + j * db
        if ((e - a0) / da).denominator == 1:
            if e > exact(s441.init) + s441.n * da:
                break
            return e
    raise ValueError(f"record grids of {s441} and {s442} never meet near the DE442 start")


def build(args):
    sources_dir = os.path.join(args.data, "sources")
    packs_dir = os.path.join(args.data, "packs")
    os.makedirs(packs_dir, exist_ok=True)

    de442_path = os.path.join(sources_dir, "de442.bsp")
    de442 = SPK.open(de442_path)
    if args.range:
        y0, y1 = args.range
        start_jd = calendar_to_jd(y0) - args.margin_days
        end_jd = calendar_to_jd(y1 + 1) - 0.0 + args.margin_days
        name = f"fepk-de442-{y0}-{y1}"
        use441 = False
    else:
        start_jd = FULL_START_JD - args.margin_days
        end_jd = FULL_END_JD + args.margin_days
        name = "fepk-de441-de442"
        use441 = True
    out_path = args.out or os.path.join(packs_dir, f"{name}-v1.bin")

    de441 = None
    if use441:
        de441 = SPK.open(NAIF + "de441_part-1.bsp", os.path.join(sources_dir, "de441_cache"))

    print(f"building {out_path}: JD {start_jd} .. {end_jd}", flush=True)

    # ---------------- plan series ----------------
    plan = []  # (target, center, source_index, seg, i0, i1)
    seams = []
    t_start, t_end = exact(jd_to_et_exact(start_jd)), exact(jd_to_et_exact(end_jd))
    for tgt, ctr in PAIRS:
        s442 = de442.segment(tgt, ctr)
        b0, db = exact(s442.init), exact(s442.intlen)
        j1 = int(-((b0 - t_end) // db))  # ceil
        if j1 > s442.n:
            raise ValueError(f"requested end beyond DE442 for {tgt}")
        if t_start >= b0:
            j0 = int((t_start - b0) // db)
            plan.append((tgt, ctr, 0, s442, j0, j1))
            continue
        if de441 is None:
            raise ValueError("range starts before DE442: DE441 required (do not use --range so early)")
        s441 = de441.segment(tgt, ctr)
        a0, da = exact(s441.init), exact(s441.intlen)
        e = seam_epoch(s441, s442)
        i0 = int((t_start - a0) // da)
        i1 = int((e - a0) / da)
        j0 = int((e - b0) / db)
        if i0 < 0:
            raise ValueError("requested start before DE441 part-1")
        plan.append((tgt, ctr, 1, s441, i0, i1))
        plan.append((tgt, ctr, 0, s442, j0, j1))
        seams.append({"target": tgt, "center": ctr, "jd_tdb": float(J2000_JD + e / 86400),
                      "de441_records_end": i1, "de442_record_start": j0,
                      "grids_aligned": e == b0})

    # sources[] ordering: 0 = de442 always present; 1 = de441 when used.
    # Manifest lists them in time order, so remap: time-ordered index.
    src_order = [1, 0] if use441 else [0]
    src_index = {s: k for k, s in enumerate(src_order)}
    plan.sort(key=lambda p: (p[0], p[1], p[3].init + p[4] * p[3].intlen))

    # ---------------- layout ----------------
    table = []
    coeff_bytes_total = 0
    for tgt, ctr, src, seg, i0, i1 in plan:
        n = i1 - i0
        start = J2000_JD + (seg.init + i0 * seg.intlen) / S_PER_DAY
        interval = seg.intlen / S_PER_DAY
        end = start + n * interval
        table.append(dict(target=tgt, center=ctr, n_coeff=seg.ncoef, n_records=n, start_jd_tdb=start,
                          interval_days=interval, end_jd_tdb=end, source_index=src_index[src],
                          source=["de442.bsp", "de441_part-1.bsp"][src], first_source_record=i0))
        coeff_bytes_total += n * 3 * seg.ncoef * 8

    # pair coverage
    cov = {}
    for row in table:
        k = (row["target"], row["center"])
        lo, hi = cov.get(k, (row["start_jd_tdb"], row["end_jd_tdb"]))
        cov[k] = (min(lo, row["start_jd_tdb"]), max(hi, row["end_jd_tdb"]))
    common = (max(v[0] for v in cov.values()), min(v[1] for v in cov.values()))

    # sources
    sources = []
    for s in src_order:
        rows = [r for r in table if r["source_index"] == src_index[s]]
        used = [min(r["start_jd_tdb"] for r in rows), max(r["end_jd_tdb"] for r in rows)]
        if s == 0:
            md5 = md5_file(de442_path)
            pub = published_md5(sources_dir, "de442.bsp")
            if pub and pub != md5:
                raise ValueError(f"de442.bsp md5 {md5} != published {pub}")
            sources.append({"name": "de442.bsp", "url": NAIF + "de442.bsp", "md5": md5,
                            "md5_verified": pub == md5, "size_bytes": os.path.getsize(de442_path),
                            "coverage_jd_tdb": [de442.segments[0].start_jd, de442.segments[0].end_jd],
                            "used_jd_tdb": used, "constants": parse_constants(de442.comments())})
        else:
            f = de441.daf.f
            sources.append({"name": "de441_part-1.bsp", "url": NAIF + "de441_part-1.bsp",
                            "md5": published_md5(sources_dir, "de441_part-1.bsp"),
                            "md5_verified": False,
                            "md5_note": ("published value from aa_checksums.txt; not recomputed because only "
                                         "the needed byte ranges were fetched (HTTP Range)"),
                            "size_bytes": f.size, "last_modified": getattr(f, "last_modified", ""),
                            "coverage_jd_tdb": [de441.segments[0].start_jd, de441.segments[0].end_jd],
                            "used_jd_tdb": used, "constants": parse_constants(de441.comments())})
    de442_consts = sources[src_index[0]]["constants"]
    constants = dict(de442_consts)
    constants["source"] = "de442.bsp kernel comments (see sources[].constants for each kernel)"
    if use441:
        # AU and CLIGHT are printed only in the DE441 comments; DE44x share them.
        c441 = sources[src_index[1]]["constants"]
        for k in ("AU_km", "CLIGHT_km_s"):
            if k in c441:
                constants[k] = c441[k]
                constants[k + "_source"] = "de441_part-1.bsp kernel comments (same value in all DE44x)"

    manifest = {
        "name": name,
        "format_version": FORMAT_VERSION,
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builder": {"tool": "tools/build_pack.py", "version": BUILDER_VERSION,
                    "argv": sys.argv[1:]},
        "requested_jd_tdb": [start_jd, end_jd],
        "margin_days": args.margin_days,
        "coverage_jd_tdb": list(common),
        "pair_coverage_jd_tdb": {f"{t}/{c}": list(v) for (t, c), v in sorted(cov.items())},
        "sources": sources,
        "seams": seams,
        "series": table,
        "constants": constants,
        "units": {"position": "km", "time": "JD TDB", "frame": "ICRF (J2000 frame id 1)"},
        "license": "JPL ephemeris data: public, see https://naif.jpl.nasa.gov/naif/rules.html",
    }
    man_bytes = json.dumps(manifest, indent=1).encode()

    table_off = HEADER_LEN
    man_off = table_off + len(table) * SERIES_LEN
    data_off = (man_off + len(man_bytes) + 7) // 8 * 8
    total = data_off + coeff_bytes_total

    free = _free_bytes(packs_dir)
    if free < total + (512 << 20):
        raise SystemExit(f"not enough disk: need {total/1e6:.0f} MB + 512 MB slack, have {free/1e6:.0f} MB")

    # ---------------- write ----------------
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as f:
        hdr = MAGIC + struct.pack("<II", FORMAT_VERSION, len(table))
        hdr += struct.pack("<QQQQ", table_off, man_off, len(man_bytes), data_off) + bytes(16)
        assert len(hdr) == HEADER_LEN
        f.write(hdr)
        off = data_off
        for row in table:
            row["coeff_offset"] = off
            f.write(struct.pack("<iiII", row["target"], row["center"], row["n_coeff"], row["n_records"]))
            f.write(struct.pack("<dddQI", row["start_jd_tdb"], row["interval_days"], row["end_jd_tdb"],
                                off, row["source_index"]))
            f.write(bytes(12))
            off += row["n_records"] * 3 * row["n_coeff"] * 8
        f.write(man_bytes)
        f.write(bytes(data_off - man_off - len(man_bytes)))
        assert f.tell() == data_off
        for (tgt, ctr, src, seg, i0, i1), row in zip(plan, table):
            assert f.tell() == row["coeff_offset"]
            step = max(1, (8 << 20) // (seg.rsize * 8))  # ~8 MB per read
            for a in range(i0, i1, step):
                b = min(a + step, i1)
                recs = seg.records(a, b)
                check_grid(seg, a, recs)
                # bit-identical copy of the X/Y/Z coefficient words
                f.write(np.ascontiguousarray(recs[:, 2:]).astype("<f8", copy=False).tobytes())
            print(f"  {tgt:>3} wrt {ctr}: {row['source']:<16} records {i0}..{i1} "
                  f"({row['n_records']} x {row['n_coeff']} coeffs, {row['interval_days']:g} d) "
                  f"JD {row['start_jd_tdb']}..{row['end_jd_tdb']}", flush=True)
        assert f.tell() == total
    os.replace(tmp, out_path)

    side = dict(manifest)
    side["file"] = {"path": os.path.basename(out_path), "size_bytes": total, "sha256": sha256_file(out_path)}
    if use441:
        side["de441_bytes_fetched_this_run"] = de441.daf.f.bytes_fetched
    with open(out_path + ".manifest.json", "w") as f:
        json.dump(side, f, indent=1)
    print(f"wrote {out_path}: {total} bytes ({total/1e6:.1f} MB), {len(table)} series, "
          f"coverage JD {common[0]} .. {common[1]}")
    for s in seams:
        print(f"  seam {s['target']}/{s['center']} at JD {s['jd_tdb']} aligned={s['grids_aligned']}")
    return out_path


def jd_to_et_exact(jd: float) -> Fraction:
    return (Fraction(jd) - Fraction(J2000_JD)) * 86400


def _free_bytes(path: str) -> int:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--range", nargs=2, type=int, metavar=("Y0", "Y1"),
                    help="dev pack from DE442 only covering years Y0..Y1 inclusive")
    ap.add_argument("--margin-days", type=float, default=64.0)
    ap.add_argument("--data", default=DEFAULT_DATA, help="data dir with sources/ and packs/")
    ap.add_argument("--out", help="output path")
    build(ap.parse_args())


if __name__ == "__main__":
    main()
