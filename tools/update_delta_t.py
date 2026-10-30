#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Refresh the Delta T sources from USNO / IERS before a release.

    python tools/update_delta_t.py            # download, write the sources
    python crates/fe-core/tables/gen_tables.py   # then regenerate time/tables.rs

Writes, in crates/fe-core/tables/sources/:
  - USNO-deltat.data            observed monthly Delta T (https://maia.usno.navy.mil/ser7/deltat.data)
  - IERS-predicted-deltat.data  monthly Delta T predicted by IERS Bulletin A for the months after the
                                last observed value: Delta T = 32.184 s + (TAI - UTC) - (UT1 - UTC),
                                UT1 - UTC from the "P" (predicted) rows of finals2000A.data
                                (https://maia.usno.navy.mil/ser7/finals2000A.data)

After the predicted months the engine extrapolates (time.rs). Re-run this every release: the
predictions become observations, and the forecast of record is kept in docs/DELTA_T_REPORT.md.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import urllib.request

SER7 = "https://maia.usno.navy.mil/ser7/"
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "crates", "fe-core", "tables", "sources")
TT_TAI = 32.184
TAI_UTC = 37.0  # since 2017-01-01; utc.rs carries the leap-second table


def fetch(name: str) -> str:
    req = urllib.request.Request(SER7 + name, headers={"User-Agent": "frenchephemeris"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def predicted_rows(finals: str, after: tuple[int, int, int]) -> list[tuple[int, int, int, float]]:
    """First-of-month Delta T from the predicted UT1 - UTC rows, strictly after `after`."""
    out = []
    for line in finals.splitlines():
        if len(line) < 68 or line[57] != "P" or line[4:6].strip() != "1":
            continue
        yy, mm = int(line[0:2]), int(line[2:4])
        year = 2000 + yy if yy < 70 else 1900 + yy
        if (year, mm, 1) > after:
            out.append((year, mm, 1, TT_TAI + TAI_UTC - float(line[58:68])))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", default=SRC)
    args = ap.parse_args()
    observed = fetch("deltat.data")
    last = [line.split() for line in observed.splitlines() if len(line.split()) == 4][-1]
    after = (int(last[0]), int(last[1]), int(last[2]))
    rows = predicted_rows(fetch("finals2000A.data"), after)
    if not rows:
        raise SystemExit("no predicted months after the last observed value")
    with open(os.path.join(args.src, "USNO-deltat.data"), "w") as f:
        f.write(observed)
    today = _dt.date.today().isoformat()
    with open(os.path.join(args.src, "IERS-predicted-deltat.data"), "w") as f:
        f.write(f"# IERS Bulletin A predictions (finals2000A.data, retrieved {today}):\n")
        f.write(f"# Delta T = {TT_TAI} + {TAI_UTC} - (UT1-UTC), first of each month after the last USNO value\n")
        for y, m, d, v in rows:
            f.write(f" {y} {m:2d} {d:2d}  {v:.4f}\n")
    print(f"observed to {after[0]}-{after[1]:02d}-{after[2]:02d} ({last[3]} s); "
          f"{len(rows)} predicted months to {rows[-1][0]}-{rows[-1][1]:02d} ({rows[-1][3]:.4f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
