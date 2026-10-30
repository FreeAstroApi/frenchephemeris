# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Minimal vectorised FEPK v1 reader (docs/FORMAT.md) for the offline tools.

    pk = Fepk(path); pk.pos(10, 0, jd_array) -> (n, 3) km

Only stored pairs are evaluated (no chaining); the seam rule of the format
applies (an epoch on a series boundary uses the later series).
"""
from __future__ import annotations

import json
import struct

import numpy as np


class Fepk:
    def __init__(self, path: str):
        self.path = path
        self.mm = np.memmap(path, dtype="<f8", mode="r")
        raw = open(path, "rb").read(64)
        if raw[:8] != b"FEPK\0\0\0\x01":
            raise ValueError(f"{path}: not an FEPK v1 file")
        (_ver, n) = struct.unpack("<II", raw[8:16])
        table_off, man_off, man_len, _data_off = struct.unpack("<QQQQ", raw[16:48])
        with open(path, "rb") as f:
            f.seek(table_off)
            table = f.read(n * 64)
            f.seek(man_off)
            self.manifest = json.loads(f.read(man_len))
        self.series = []
        for k in range(n):
            o = k * 64
            tgt, ctr, nc, nr = struct.unpack("<iiII", table[o:o + 16])
            start, interval, end, off, src = struct.unpack("<dddQI", table[o + 16:o + 52])
            self.series.append(dict(target=tgt, center=ctr, n_coeff=nc, n_records=nr, start=start,
                                    interval=interval, end=end, offset=off, source_index=src))

    def pairs(self):
        return sorted({(s["target"], s["center"]) for s in self.series})

    def coverage(self, target: int, center: int):
        ss = [s for s in self.series if s["target"] == target and s["center"] == center]
        return min(s["start"] for s in ss), max(s["end"] for s in ss)

    def state(self, target: int, center: int, jd) -> tuple[np.ndarray, np.ndarray]:
        """(pos km, vel km/day), shape (n, 3), from one stored pair (later series on a seam)."""
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        ss = sorted((s for s in self.series if s["target"] == target and s["center"] == center),
                    key=lambda s: s["start"])
        pos = np.full((len(jd), 3), np.nan)
        vel = np.full((len(jd), 3), np.nan)
        done = np.zeros(len(jd), bool)
        for s in reversed(ss):
            m = (~done) & (jd >= s["start"]) & (jd <= s["end"])
            if not m.any():
                continue
            n, nr = s["n_coeff"], s["n_records"]
            i = np.minimum(((jd[m] - s["start"]) / s["interval"]).astype(np.int64), nr - 1)
            u = 2.0 * (jd[m] - (s["start"] + i * s["interval"])) / s["interval"] - 1.0
            base = s["offset"] // 8
            c = np.asarray(self.mm[base:base + nr * 3 * n]).reshape(nr, 3, n)[i]
            t = np.zeros((len(u), n))
            d = np.zeros((len(u), n))
            t[:, 0] = 1.0
            t[:, 1] = u
            d[:, 1] = 1.0
            for k in range(2, n):
                t[:, k] = 2 * u * t[:, k - 1] - t[:, k - 2]
                d[:, k] = 2 * u * d[:, k - 1] + 2 * t[:, k - 1] - d[:, k - 2]
            pos[m] = np.einsum("rck,rk->rc", c, t)
            vel[m] = np.einsum("rck,rk->rc", c, d) * 2.0 / s["interval"]
            done |= m
        if not done.all():
            raise ValueError(f"epoch outside coverage of {target}/{center}")
        return pos, vel

    def pos(self, target: int, center: int, jd, within: tuple[float, float] | None = None) -> np.ndarray:
        """Positions (n, 3) km. `within=(a, b)` restricts evaluation to the one series
        covering [a, b] (used to take a seam epoch from the earlier side)."""
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        ss = sorted((s for s in self.series if s["target"] == target and s["center"] == center),
                    key=lambda s: s["start"])
        if within is not None:
            ss = [s for s in ss if s["start"] <= within[0] and within[1] <= s["end"]]
        if not ss:
            raise KeyError(f"pair {target}/{center} not stored")
        out = np.full((len(jd), 3), np.nan)
        done = np.zeros(len(jd), bool)
        for s in reversed(ss):
            m = (~done) & (jd >= s["start"]) & (jd <= s["end"])
            if not m.any():
                continue
            n, nr = s["n_coeff"], s["n_records"]
            i = np.minimum(((jd[m] - s["start"]) / s["interval"]).astype(np.int64), nr - 1)
            t0 = s["start"] + i * s["interval"]
            u = 2.0 * (jd[m] - t0) / s["interval"] - 1.0
            base = s["offset"] // 8
            c = np.asarray(self.mm[base:base + nr * 3 * n]).reshape(nr, 3, n)[i]
            t = np.empty((len(u), n))
            t[:, 0] = 1.0
            t[:, 1] = u
            for k in range(2, n):
                t[:, k] = 2 * u * t[:, k - 1] - t[:, k - 2]
            out[m] = np.einsum("rck,rk->rc", c, t)
            done |= m
        if not done.all():
            raise ValueError(f"epoch outside coverage of {target}/{center}")
        return out
