# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Positions per second through `swe_compat` (one body per call, with speeds; tropical and
Lahiri sidereal) on 1 and N threads, and through the batch API.

    python tools/bench_throughput.py [--threads 8] [--swiss]

`--swiss` measures pyswisseph the same way (needs `SE_EPHE_PATH`); its settings are
process-wide, so its sidereal runs set the mode once before timing.
"""

from __future__ import annotations

import argparse
import os
import random
import threading
import time

import numpy as np

BODIES = list(range(10))
N_CALLS = 200_000


def _jds(n, seed=1):
    rnd = random.Random(seed)
    return [2415020.5 + rnd.random() * 73049.0 for _ in range(n)]


def per_call_rate(mod, flags, threads, sidereal):
    jds = _jds(N_CALLS // len(BODIES))
    per = len(jds) // threads
    barrier = threading.Barrier(threads + 1)

    def work(chunk):
        if sidereal:
            mod.set_sid_mode(mod.SIDM_LAHIRI)
        barrier.wait()
        for jd in chunk:
            for b in BODIES:
                mod.calc_ut(jd, b, flags)

    ts = [threading.Thread(target=work, args=(jds[i * per:(i + 1) * per],)) for i in range(threads)]
    for t in ts:
        t.start()
    barrier.wait()
    t0 = time.perf_counter()
    for t in ts:
        t.join()
    return per * threads * len(BODIES) / (time.perf_counter() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--swiss", action="store_true")
    a = ap.parse_args()

    import frenchephemeris as fe
    import frenchephemeris.swe_compat as fswe

    engines = [("FrenchEphemeris", fswe)]
    if a.swiss:
        import swisseph as sw

        sw.set_ephe_path(os.environ["SE_EPHE_PATH"])
        engines.append(("Swiss Ephemeris", sw))
    for name, mod in engines:
        for sidereal in (False, True):
            flags = mod.FLG_SWIEPH * (name != "FrenchEphemeris") | mod.FLG_SPEED
            flags |= mod.FLG_SIDEREAL if sidereal else 0
            label = "sidereal" if sidereal else "tropical"
            for th in (1, a.threads):
                mod.calc_ut(2451545.0, 0, flags)  # warm-up
                print(f"{name:16} {label:8} {th} thread(s): {per_call_rate(mod, flags, th, sidereal):>10,.0f} positions/s")

    eng = fe.open_default()
    times = np.array(_jds(10_000, seed=2))
    eng.calc_batch_ut(times[:10], BODIES)
    t0 = time.perf_counter()
    eng.calc_batch_ut(times, BODIES)
    print(f"{'FrenchEphemeris':16} batch    {len(times) * len(BODIES) / (time.perf_counter() - t0):>10,.0f} positions/s")


if __name__ == "__main__":
    main()
