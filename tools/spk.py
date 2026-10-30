# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Minimal DAF / SPK reader: type 2 (Chebyshev) and type 21 (extended
modified difference arrays), independent of jplephem.

Works on a local file or on an HTTP-range-backed file object (`HttpRangeFile`),
so that a handful of segments can be extracted from a multi-GB kernel without
downloading it.

References: NAIF "DAF Required Reading" and "SPK Required Reading" (type 2:
Chebyshev, position only).

Layout reminders
----------------
* File record (record 1, 1024 bytes): LOCIDW[8] ND NI LOCIFN[60] FWARD BWARD
  FREE LOCFMT[8] ...  Records are 1024 bytes, addresses are 1-based 8-byte
  words.
* Summary record: NEXT PREV NSUM (as doubles) then NSUM summaries of
  ND doubles + NI ints (packed, padded to whole doubles). The following record
  holds the segment names (NSUM x 8*SS chars).
* SPK: ND=2 (start_et, end_et, seconds past J2000 TDB), NI=6
  (target, center, frame, type, begin_addr, end_addr).
* Type 2 segment: N records of RSIZE doubles [MID, RADIUS, X[n], Y[n], Z[n]]
  followed by the directory INIT INTLEN RSIZE N.
* Type 21 segment (JPL small-body integrator output, as written by Horizons):
  N difference lines of DLSIZE = 4*MAXDIM + 11 doubles, then the N final
  epochs (seconds TDB, increasing), then floor(N/100) directory epochs
  (every 100th final epoch), then MAXDIM and N.  A difference line is
      TL, G[MAXDIM], (REFPOS_k, REFVEL_k) k=x,y,z interleaved (6),
      DT[MAXDIM][3] (column-major: all MAXDIM of x, then y, then z),
      KQMAX1, KQ[3].
  TL is the line's reference epoch; the line is used for epochs in
  (previous final epoch, final epoch] (first line whose final epoch >= t).
  In Horizons files TL usually equals the previous line's final epoch.  Evaluation is Krogh's
  modified divided-difference interpolation; see `Type21` for the algorithm.
"""
from __future__ import annotations

import os
import struct
import subprocess
from dataclasses import dataclass

import numpy as np

S_PER_DAY = 86400.0
J2000_JD = 2451545.0
RECLEN = 1024


def _groups(idx: np.ndarray):
    """Yield (value, index array) for each distinct value of an int array (O(n log n))."""
    order = np.argsort(idx, kind="stable")
    sidx = idx[order]
    cuts = np.flatnonzero(np.diff(sidx)) + 1
    for part in np.split(order, cuts):
        if len(part):
            yield int(idx[part[0]]), part


def et_to_jd(et: float) -> float:
    return J2000_JD + et / S_PER_DAY


def jd_to_et(jd: float) -> float:
    return (jd - J2000_JD) * S_PER_DAY


# --------------------------------------------------------------------------
# File-like backends exposing read_at(offset, n) -> bytes
# --------------------------------------------------------------------------
class LocalFile:
    def __init__(self, path: str):
        self.path = path
        self._f = open(path, "rb")
        self.size = os.fstat(self._f.fileno()).st_size

    def read_at(self, offset: int, n: int) -> bytes:
        self._f.seek(offset)
        b = self._f.read(n)
        if len(b) != n:
            raise EOFError(f"short read at {offset}+{n} in {self.path}")
        return b

    def close(self):
        self._f.close()


def _curl_range(url: str, lo: int, hi: int) -> bytes:
    """Bytes lo..hi inclusive of `url` (curl: system TLS trust store), with retries."""
    import time as _time
    for attempt in range(12):
        # abort stalled transfers (< 1 kB/s for 90 s) and retry
        r = subprocess.run(["curl", "-sSf", "--max-time", "1800", "--speed-limit", "1000",
                            "--speed-time", "90", "-r", f"{lo}-{hi}", url], capture_output=True)
        if r.returncode == 0 and len(r.stdout) == hi - lo + 1:
            return r.stdout
        _time.sleep(min(60, 5 * (attempt + 1)))
        if attempt == 11:
            raise IOError(f"range fetch {lo}-{hi} failed: rc={r.returncode} len={len(r.stdout)} "
                          f"{r.stderr.decode(errors='replace')}")
    raise AssertionError


class StreamRangeFile:
    """Remote file read through HTTP Range requests, nothing written to disk.

    Small reads (file record, summaries, segment directories) are served from an
    in-memory cache of 64 KiB blocks; large reads are fetched directly, split
    into `parallel` concurrent requests. Used to stream multi-GB kernels.
    """

    BLOCK = 1 << 16

    def __init__(self, url: str, parallel: int = 6):
        self.url = url
        self.parallel = parallel
        hdr = subprocess.run(["curl", "-sSIL", "--max-time", "60", url], check=True,
                             capture_output=True, text=True).stdout
        heads = {}
        for line in hdr.splitlines():
            if ":" in line:
                k, v = line.split(":", 1)
                heads[k.strip().lower()] = v.strip()
        self.size = int(heads["content-length"])
        self.last_modified = heads.get("last-modified", "")
        self.bytes_fetched = 0
        self._blocks: dict[int, bytes] = {}

    def read_at(self, offset: int, n: int) -> bytes:
        if offset < 0 or offset + n > self.size:
            raise EOFError("read beyond end of remote file")
        if n == 0:
            return b""
        if n <= 4 * self.BLOCK:
            out = bytearray()
            for i in range(offset // self.BLOCK, (offset + n - 1) // self.BLOCK + 1):
                if i not in self._blocks:
                    lo = i * self.BLOCK
                    hi = min(lo + self.BLOCK, self.size) - 1
                    self._blocks[i] = _curl_range(self.url, lo, hi)
                    self.bytes_fetched += hi - lo + 1
                blk = self._blocks[i]
                a = max(offset - i * self.BLOCK, 0)
                b = min(offset + n - i * self.BLOCK, len(blk))
                out += blk[a:b]
            return bytes(out)
        from concurrent.futures import ThreadPoolExecutor
        k = max(1, min(self.parallel, n // (1 << 17)))
        cuts = [offset + (n * j) // k for j in range(k + 1)]
        with ThreadPoolExecutor(k) as ex:
            parts = list(ex.map(lambda j: _curl_range(self.url, cuts[j], cuts[j + 1] - 1), range(k)))
        self.bytes_fetched += n
        return b"".join(parts)

    def close(self):
        self._blocks.clear()


class HttpRangeFile:
    """Read-only remote file using HTTP Range requests with an on-disk chunk cache.

    The file is split into fixed `chunk` byte blocks; each block fetched is
    stored as `<cache_dir>/<chunk_index>.bin`.  Runs of missing blocks are
    fetched with a single request.  `bytes_fetched` counts network bytes.
    """

    def __init__(self, url: str, cache_dir: str, chunk: int = 1 << 20, size: int | None = None):
        self.url = url
        self.cache_dir = cache_dir
        self.chunk = chunk
        os.makedirs(cache_dir, exist_ok=True)
        meta = os.path.join(cache_dir, "meta.txt")
        if size is None and os.path.exists(meta):
            with open(meta) as f:
                kv = dict(line.strip().split("=", 1) for line in f if "=" in line)
            if kv.get("url") == url and int(kv.get("chunk", 0)) == chunk:
                size = int(kv["size"])
                self.last_modified = kv.get("last_modified", "")
        if size is None:
            hdr = subprocess.run(["curl", "-sSIL", "--max-time", "60", url], check=True,
                                 capture_output=True, text=True).stdout
            heads = {}
            for line in hdr.splitlines():
                if ":" in line:
                    k, v = line.split(":", 1)
                    heads[k.strip().lower()] = v.strip()
            size = int(heads["content-length"])
            self.last_modified = heads.get("last-modified", "")
            with open(meta, "w") as f:
                f.write(f"url={url}\nsize={size}\nchunk={chunk}\nlast_modified={self.last_modified}\n")
        self.size = size
        self.last_modified = getattr(self, "last_modified", "")
        self.bytes_fetched = 0

    def _path(self, i: int) -> str:
        return os.path.join(self.cache_dir, f"{i:06d}.bin")

    def _chunk_len(self, i: int) -> int:
        return min(self.chunk, self.size - i * self.chunk)

    def _fetch(self, i0: int, i1: int):
        """Fetch chunks i0..i1 inclusive in one request and cache them."""
        lo = i0 * self.chunk
        hi = min((i1 + 1) * self.chunk, self.size) - 1
        # curl (system TLS trust store) rather than urllib: python builds without
        # a CA bundle fail certificate verification against naif.jpl.nasa.gov.
        for attempt in range(5):
            r = subprocess.run(["curl", "-sSf", "--max-time", "600", "-r", f"{lo}-{hi}", self.url],
                               capture_output=True)
            data = r.stdout
            if r.returncode == 0 and len(data) == hi - lo + 1:
                break
            if attempt == 4:
                raise IOError(f"range fetch {lo}-{hi} failed: rc={r.returncode} len={len(data)} "
                              f"{r.stderr.decode(errors='replace')}")
        if len(data) != hi - lo + 1:
            raise IOError(f"short range response {len(data)} != {hi - lo + 1}")
        self.bytes_fetched += len(data)
        for i in range(i0, i1 + 1):
            a = (i - i0) * self.chunk
            blob = data[a:a + self._chunk_len(i)]
            tmp = self._path(i) + ".tmp"
            with open(tmp, "wb") as f:
                f.write(blob)
            os.replace(tmp, self._path(i))

    def prefetch(self, offset: int, n: int, max_run: int = 16):
        i0, i1 = offset // self.chunk, (offset + n - 1) // self.chunk
        run_start = None
        for i in range(i0, i1 + 2):
            missing = i <= i1 and not os.path.exists(self._path(i))
            if missing and run_start is None:
                run_start = i
            if run_start is not None and (not missing or i - run_start >= max_run):
                self._fetch(run_start, i - 1)
                run_start = i if missing else None

    def read_at(self, offset: int, n: int) -> bytes:
        if offset < 0 or offset + n > self.size:
            raise EOFError("read beyond end of remote file")
        if n == 0:
            return b""
        self.prefetch(offset, n)
        out = bytearray()
        i0, i1 = offset // self.chunk, (offset + n - 1) // self.chunk
        for i in range(i0, i1 + 1):
            with open(self._path(i), "rb") as f:
                blob = f.read()
            a = max(offset - i * self.chunk, 0)
            b = min(offset + n - i * self.chunk, len(blob))
            out += blob[a:b]
        return bytes(out)

    def cached_chunks(self) -> list[int]:
        return sorted(int(n[:-4]) for n in os.listdir(self.cache_dir) if n.endswith(".bin"))

    def write_sparse_copy(self, path: str) -> int:
        """Write a sparse local file of the full remote size holding every cached chunk
        at its true offset (holes read as zeros). Lets an independent reader such as
        jplephem open the kernel normally while only the fetched ranges use disk.
        Returns the number of bytes written."""
        written = 0
        with open(path, "wb") as f:
            f.truncate(self.size)
            for i in self.cached_chunks():
                with open(self._path(i), "rb") as c:
                    blob = c.read()
                f.seek(i * self.chunk)
                f.write(blob)
                written += len(blob)
        return written

    def close(self):
        pass


# --------------------------------------------------------------------------
# DAF / SPK
# --------------------------------------------------------------------------
@dataclass
class Segment:
    daf: "DAF"
    name: str
    start_et: float
    end_et: float
    target: int
    center: int
    frame: int
    data_type: int
    start_i: int  # 1-based word address of first double
    end_i: int    # 1-based word address of last double
    # type 2 directory
    init: float = 0.0
    intlen: float = 0.0
    rsize: int = 0
    n: int = 0

    @property
    def ncoef(self) -> int:
        return (self.rsize - 2) // 3

    @property
    def start_jd(self) -> float:
        return et_to_jd(self.start_et)

    @property
    def end_jd(self) -> float:
        return et_to_jd(self.end_et)

    @property
    def init_jd(self) -> float:
        return et_to_jd(self.init)

    @property
    def interval_days(self) -> float:
        return self.intlen / S_PER_DAY

    def record_epoch_et(self, i: int) -> float:
        """Start epoch (seconds) of record i."""
        return self.init + i * self.intlen

    def record_byte_range(self, i0: int, i1: int) -> tuple[int, int]:
        """(offset, nbytes) of records [i0, i1)."""
        off = (self.start_i - 1) * 8 + i0 * self.rsize * 8
        return off, (i1 - i0) * self.rsize * 8

    def records(self, i0: int, i1: int) -> np.ndarray:
        """Raw records [i0, i1) as an (i1-i0, rsize) float64 array (bit-exact)."""
        if not (0 <= i0 <= i1 <= self.n):
            raise IndexError(f"record slice {i0}:{i1} outside 0:{self.n}")
        off, nb = self.record_byte_range(i0, i1)
        raw = self.daf.f.read_at(off, nb)
        return np.frombuffer(raw, dtype=self.daf.dtype).astype("<f8", copy=False).reshape(i1 - i0, self.rsize)

    def record_index_range(self, jd0: float, jd1: float) -> tuple[int, int]:
        """Smallest [i0, i1) of records whose union covers [jd0, jd1]."""
        t0, t1 = jd_to_et(jd0), jd_to_et(jd1)
        if t0 < self.init or t1 > self.init + self.n * self.intlen:
            raise ValueError(f"JD range {jd0}..{jd1} outside segment {self.start_jd}..{self.end_jd}")
        i0 = int((t0 - self.init) // self.intlen)
        i1 = int(-((self.init - t1) // self.intlen))  # ceil
        i1 = max(i1, i0 + 1)
        return i0, min(i1, self.n)

    def records_for_jd(self, jd0: float, jd1: float) -> tuple[int, np.ndarray]:
        i0, i1 = self.record_index_range(jd0, jd1)
        return i0, self.records(i0, i1)

    # Reference evaluation (used for self-checks; vectorised over one record)
    def position(self, jd: float) -> np.ndarray:
        t = jd_to_et(jd)
        i = min(int((t - self.init) // self.intlen), self.n - 1)
        rec = self.records(i, i + 1)[0]
        mid, rad = rec[0], rec[1]
        s = (t - mid) / rad
        c = rec[2:].reshape(3, self.ncoef)
        return np.polynomial.chebyshev.chebval(s, c.T)

    def cheb_state(self, et: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised type-2 (pos km, vel km/s) at epochs `et` (seconds)."""
        idx = np.clip(((et - self.init) // self.intlen).astype(np.int64), 0, self.n - 1)
        pos = np.empty((len(et), 3))
        vel = np.empty((len(et), 3))
        for i, m in _groups(idx):
            rec = self.records(int(i), int(i) + 1)[0]
            mid, rad = rec[0], rec[1]
            x = (et[m] - mid) / rad
            c = rec[2:].reshape(3, self.ncoef)
            pos[m] = np.polynomial.chebyshev.chebval(x, c.T).T
            dc = np.polynomial.chebyshev.chebder(c.T)
            vel[m] = np.polynomial.chebyshev.chebval(x, dc).T / rad
        return pos, vel

    def __repr__(self):
        return (f"Segment({self.target} wrt {self.center} frame={self.frame} type={self.data_type} "
                f"jd {self.start_jd:.1f}..{self.end_jd:.1f} n={self.n} ncoef={self.ncoef} "
                f"interval={self.interval_days:g}d)")


class Type21:
    """SPK type 21 (extended modified difference arrays) evaluator.

    Implemented from the algorithm description of the NAIF SPK Required
    Reading (types 1/21) and Krogh's variable-order Adams integrator output
    form: with  delta = t - TL  and the step-size vector G,

        FC[j+1] = tp/G[j],  WC[j] = delta/G[j],  tp <- delta + G[j]   (j < KQMAX1-2)
        W[j]    = 1/j                                          (j = 1..KQMAX1)
        repeated integration of the W vector (KQMAX1-2 passes), then
        pos_i = REFPOS_i + delta*(REFVEL_i + delta * sum_j DT[j,i] W[j+1])
        one more pass, then
        vel_i = REFVEL_i + delta * sum_j DT[j,i] W[j].

    Arrays are 1-based in the formulas; the code below keeps that indexing
    (index 0 unused) to stay checkable against them.
    """

    def __init__(self, seg: "Segment"):
        f, e = seg.daf.f, seg.daf.endian
        self.seg = seg
        maxdim, n = struct.unpack(e + "2d", f.read_at((seg.end_i - 2) * 8, 16))
        self.maxdim, self.n = int(maxdim), int(n)
        self.dlsize = 4 * self.maxdim + 11
        ndir = self.n // 100
        expect = seg.start_i + self.n * self.dlsize + self.n + ndir + 2 - 1
        if expect != seg.end_i:
            raise ValueError(f"inconsistent type-21 segment: maxdim={self.maxdim} n={self.n} "
                             f"words {seg.start_i}..{seg.end_i} (expected end {expect})")
        off = (seg.start_i - 1) * 8
        raw = f.read_at(off, (self.n * self.dlsize + self.n) * 8)
        arr = np.frombuffer(raw, dtype=seg.daf.dtype).astype("<f8")
        self.lines = arr[: self.n * self.dlsize].reshape(self.n, self.dlsize)
        self.epochs = arr[self.n * self.dlsize:].copy()
        if not np.all(np.diff(self.epochs) > 0):
            raise ValueError("type-21 epochs not increasing")
        if self.epochs[-1] < seg.end_et:
            raise ValueError("type-21 lines end before the segment end")
        # Informational: in Horizons files most lines start where the previous
        # one ended (TL[i+1] == final epoch[i]).
        self.tiled_fraction = float(np.mean(self.lines[1:, 0] == self.epochs[:-1])) if self.n > 1 else 1.0

    def record_index(self, et: np.ndarray) -> np.ndarray:
        # First record whose final epoch TL >= et.
        return np.minimum(np.searchsorted(self.epochs, et, side="left"), self.n - 1)

    def state(self, et) -> tuple[np.ndarray, np.ndarray]:
        et = np.atleast_1d(np.asarray(et, dtype=np.float64))
        idx = self.record_index(et)
        pos = np.empty((len(et), 3))
        vel = np.empty((len(et), 3))
        for i, m in _groups(idx):
            pos[m], vel[m] = self._eval_line(self.lines[i], et[m])
        return pos, vel

    def _eval_line(self, rec: np.ndarray, et: np.ndarray):
        md = self.maxdim
        tl = rec[0]
        g = np.concatenate(([0.0], rec[1:1 + md]))             # G[1..md]
        rv = rec[1 + md:7 + md]
        refpos = rv[0::2]
        refvel = rv[1::2]
        dt = rec[7 + md:7 + 4 * md].reshape(3, md)              # DT[i][j-1]
        kqmax1 = int(rec[7 + 4 * md])
        kq = rec[8 + 4 * md:11 + 4 * md].astype(int)
        m = len(et)
        delta = et - tl
        mq2 = kqmax1 - 2
        fc = np.zeros((kqmax1 + 2, m))
        wc = np.zeros((kqmax1 + 2, m))
        tp = delta.copy()
        for j in range(1, mq2 + 1):
            fc[j + 1] = tp / g[j]
            wc[j] = delta / g[j]
            tp = delta + g[j]
        w = np.zeros((2 * kqmax1 + 2, m))
        for j in range(1, kqmax1 + 1):
            w[j] = 1.0 / j
        ks = kqmax1 - 1
        ks1 = ks - 1
        jx = 0
        while ks >= 2:
            jx += 1
            for j in range(1, jx + 1):
                w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks]
            ks = ks1
            ks1 -= 1
        pos = np.empty((m, 3))
        for i in range(3):
            acc = np.zeros(m)
            for j in range(kq[i], 0, -1):
                acc += dt[i, j - 1] * w[j + ks]
            pos[:, i] = refpos[i] + delta * (refvel[i] + delta * acc)
        for j in range(1, jx + 1):
            w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks]
        ks -= 1
        vel = np.empty((m, 3))
        for i in range(3):
            acc = np.zeros(m)
            for j in range(kq[i], 0, -1):
                acc += dt[i, j - 1] * w[j + ks]
            vel[:, i] = refvel[i] + delta * acc
        return pos, vel


class DAF:
    def __init__(self, f):
        self.f = f
        rec = f.read_at(0, RECLEN)
        self.idword = rec[:8].decode("ascii", "replace")
        if not self.idword.startswith("DAF/"):
            raise ValueError(f"not a DAF file: {self.idword!r}")
        locfmt = rec[88:96].decode("ascii", "replace")
        if locfmt == "LTL-IEEE":
            e = "<"
        elif locfmt == "BIG-IEEE":
            e = ">"
        else:
            raise ValueError(f"unsupported binary format {locfmt!r}")
        self.endian = e
        self.dtype = np.dtype(e + "f8")
        self.nd, self.ni = struct.unpack(e + "ii", rec[8:16])
        self.ifname = rec[16:76].decode("ascii", "replace").strip()
        self.fward, self.bward, self.free = struct.unpack(e + "iii", rec[76:88])
        self.ss = self.nd + (self.ni + 1) // 2  # summary size in doubles

    def comments(self) -> str:
        """Comment area (records 2..FWARD-1): 1000 chars/record, NUL = newline, EOT ends."""
        nrec = self.fward - 2
        if nrec <= 0:
            return ""
        raw = self.f.read_at(RECLEN, nrec * RECLEN)
        parts = [raw[i * RECLEN:i * RECLEN + 1000] for i in range(nrec)]
        text = b"".join(parts)
        end = text.find(b"\x04")
        if end >= 0:
            text = text[:end]
        return text.replace(b"\x00", b"\n").decode("ascii", "replace")

    def summaries(self):
        e = self.endian
        recno = self.fward
        nd, ni, ss = self.nd, self.ni, self.ss
        while recno:
            srec = self.f.read_at((recno - 1) * RECLEN, 2 * RECLEN)
            nrec = srec[RECLEN:]
            nxt, _prev, nsum = struct.unpack(e + "ddd", srec[:24])
            for k in range(int(nsum)):
                base = 24 + k * ss * 8
                dvals = struct.unpack(e + f"{nd}d", srec[base:base + nd * 8])
                ivals = struct.unpack(e + f"{ni}i", srec[base + nd * 8:base + nd * 8 + ni * 4])
                name = nrec[k * ss * 8:(k + 1) * ss * 8].decode("ascii", "replace").strip()
                yield name, dvals, ivals
            recno = int(nxt)


class SPK:
    def __init__(self, f):
        self.daf = DAF(f)
        if not self.daf.idword.startswith("DAF/SPK"):
            raise ValueError(f"not an SPK: {self.daf.idword!r}")
        self.segments: list[Segment] = []
        for name, (t0, t1), (tgt, ctr, frm, typ, a, b) in self.daf.summaries():
            seg = Segment(self.daf, name, t0, t1, tgt, ctr, frm, typ, a, b)
            if typ == 2:
                e = self.daf.endian
                init, intlen, rsize, n = struct.unpack(e + "4d", f.read_at((b - 4) * 8, 32))
                seg.init, seg.intlen, seg.rsize, seg.n = init, intlen, int(rsize), int(n)
                if seg.start_i + seg.rsize * seg.n + 4 - 1 != seg.end_i:
                    raise ValueError(f"inconsistent type-2 segment {seg}")
            elif typ == 21:
                seg.t21 = Type21(seg)
            self.segments.append(seg)

    @classmethod
    def open(cls, path_or_url: str, cache_dir: str | None = None) -> "SPK":
        if path_or_url.startswith(("http://", "https://")):
            if cache_dir is None:
                raise ValueError("cache_dir required for remote kernels")
            return cls(HttpRangeFile(path_or_url, cache_dir))
        return cls(LocalFile(path_or_url))

    def segments_for(self, target: int, center: int) -> list[Segment]:
        """All segments of a pair, sorted by start epoch."""
        return sorted((s for s in self.segments if s.target == target and s.center == center),
                      key=lambda s: s.start_et)

    def state_et(self, target: int, center: int, et) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised (pos km, vel km/s), shape (n, 3), from type-2 or type-21 segments.

        When epochs fall on a boundary shared by two segments the later
        segment (greatest start <= et) is used, matching CSPICE's preference
        for the last-loaded applicable segment only up to boundary choice --
        both segments are valid there.
        """
        et = np.atleast_1d(np.asarray(et, dtype=np.float64))
        segs = self.segments_for(target, center)
        if not segs:
            raise KeyError(f"no segment {target} wrt {center}")
        pos = np.full((len(et), 3), np.nan)
        vel = np.full((len(et), 3), np.nan)
        done = np.zeros(len(et), bool)
        for s in reversed(segs):
            m = (~done) & (et >= s.start_et) & (et <= s.end_et)
            if not m.any():
                continue
            if s.data_type == 21:
                p, v = s.t21.state(et[m])
            elif s.data_type == 2:
                p, v = s.cheb_state(et[m])
            else:
                raise NotImplementedError(f"SPK type {s.data_type}")
            pos[m], vel[m] = p, v
            done |= m
        if not done.all():
            bad = et[~done][0]
            raise ValueError(f"epoch {bad} s ({et_to_jd(bad)} JD) not covered for {target} wrt {center}")
        return pos, vel

    def segment(self, target: int, center: int) -> Segment:
        found = [s for s in self.segments if s.target == target and s.center == center]
        if len(found) != 1:
            raise KeyError(f"expected one segment {target} wrt {center}, found {len(found)}")
        return found[0]

    def comments(self) -> str:
        return self.daf.comments()


if __name__ == "__main__":
    import sys

    cache = sys.argv[2] if len(sys.argv) > 2 else None
    spk = SPK.open(sys.argv[1], cache)
    print(spk.daf.idword, spk.daf.ifname, "nd/ni", spk.daf.nd, spk.daf.ni, "fward", spk.daf.fward)
    for s in spk.segments:
        extra = f"init_jd={s.init_jd}" if s.data_type == 2 else ""
        if s.data_type == 21:
            extra = f"maxdim={s.t21.maxdim} n={s.t21.n}"
        print(s, extra)
