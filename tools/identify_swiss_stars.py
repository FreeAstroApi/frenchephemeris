#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Identify every object behind the names Swiss Ephemeris accepts, independently.

Input: private/stars/swiss_probe.json (from tools/probe_swiss_stars.py).
For each distinct returned "Name,designation" string we translate the
designation into a SIMBAD identifier (",al-1Cap" -> "* alf01 Cap",
",80UMa" -> "* 80 UMa", ",HR6630" -> "HR 6630", ",M31" -> "M 31", ...) or, for
named reference points, a SIMBAD "NAME ..." identifier, and resolve it with
the SIMBAD TAP service (CfA mirror). Swiss's J2000 place is used only to
*check* the identification (< 10″) and, when the designation does not
resolve or disagrees, to pick the nearest non-planet SIMBAD object within 3′
(method "cone"). Nothing from Swiss enters the catalogue except names.

Output (private, gitignored): private/stars/star_identifications.tsv
  returned string, SIMBAD main_id, HIP (or 0), method, separation from Swiss at J2000 (″)
Objects that cannot be identified are listed with main_id "-".

Run: .venv/bin/python -I tools/identify_swiss_stars.py
"""
from __future__ import annotations

import csv
import io
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "sources" / "stars"
PRIVATE = ROOT / "private" / "stars"  # gitignored private corpus (oracle-derived data)
OUT = PRIVATE / "star_identifications.tsv"
SIMBAD_TAP = "https://simbad.harvard.edu/simbad/sim-tap/sync"

GREEK3 = {"al": "alf", "be": "bet", "ga": "gam", "de": "del", "ep": "eps", "ze": "zet", "et": "eta", "th": "tet",
          "io": "iot", "ka": "kap", "la": "lam", "mu": "mu.", "nu": "nu.", "xi": "ksi", "omi": "omi", "pi": "pi.",
          "rh": "rho", "si": "sig", "ta": "tau", "up": "ups", "ph": "phi", "ch": "chi", "ps": "psi", "ome": "ome"}
# Named reference points: SIMBAD identifiers tried for a Swiss display name.
NAMED = {
    "Galactic Center": ["NAME Sgr A*"], "Gal. Center": ["NAME Sgr A*"],
    "Great Attractor": ["NAME Great Attractor"], "Virgo Cluster": ["NAME Virgo Cluster", "NAME Vir Cluster"],
    "Andromeda Galaxy": ["M 31"], "Praesepe Cluster": ["M 44"], "Orion Nebula": ["M 42"],
}


# Swiss reference points that are not catalogue objects but conventional,
# published round-number definitions (reproduced exactly; see STARS_REPORT.md).
DEFINED = {
    "Great Attractor": "DEF Great Attractor (l=325.3, b=-7.2)",
    "Virgo Cluster": "DEF Virgo Cluster centre (B1950 12h24m +13d00m)",
}
# Swiss reference points for which no published definition could be identified.
NO_SOURCE = {"Apex", "Gal.Equ"}


def tap(query: str) -> list[dict]:
    cmd = ["curl", "-sS", "-f", "-m", "180", SIMBAD_TAP, "--data-urlencode", "REQUEST=doQuery",
           "--data-urlencode", "LANG=ADQL", "--data-urlencode", "FORMAT=csv", "--data-urlencode", f"QUERY={query}"]
    for attempt in range(4):  # the TAP mirror occasionally returns transient HTTP errors
        try:
            out = subprocess.run(cmd, check=True, capture_output=True).stdout.decode("utf-8")
            break
        except subprocess.CalledProcessError:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    return list(csv.DictReader(io.StringIO(out)))


def simbad_ids(desig: str) -> list[str]:
    m = re.fullmatch(r"(al|be|ga|de|ep|ze|et|th|io|ka|la|mu|nu|xi|omi|pi|rh|si|ta|up|ph|ch|ps|ome)(?:-(\d))?([A-Z]\w+)", desig)
    if m:
        g, k, con = m.groups()
        base = GREEK3[g]
        return [f"* {base}0{k} {con}"] if k else [f"* {base} {con}"]
    m = re.fullmatch(r"(\d+)([A-Z]\w+)", desig)
    if m:
        n = m.group(1)
        return [f"* {n} {m.group(2)}", f"* {n:>3} {m.group(2)}", f"*  {n} {m.group(2)}"]
    m = re.fullmatch(r"([a-zA-Z])([A-Z][A-Za-z]{1,2})", desig)
    if m:
        return [f"* {m.group(1)} {m.group(2)}", f"V* {m.group(1).upper()} {m.group(2)}"]
    m = re.fullmatch(r"(V\d+|[R-Z]{1,2}|[A-Q][A-Z]) (\w+)", desig)
    if m:
        return [f"V* {m.group(1)} {m.group(2)}"]
    for pre, width in (("HR", 5), ("HD", 6), ("NGC", 5), ("IC", 5), ("M", 3)):
        m = re.fullmatch(pre + r"(\d+)", desig)
        if m:  # SIMBAD pads catalogue numbers ("M  31", "NGC   869", "HR  6630")
            n = m.group(1)
            return sorted({f"{pre} {n}", f"{pre} {n:>{width}}"})
    m = re.fullmatch(r"PSR([BJ][\d+\-]+)", desig)
    if m:
        return [f"PSR {m.group(1)}"]
    return []


COLS = "b.oid, b.main_id, b.ra, b.dec, b.otype"


def sep_arcsec(ra1, de1, ra2, de2) -> float:
    r = math.radians
    c = (math.sin(r(de1)) * math.sin(r(de2)) + math.cos(r(de1)) * math.cos(r(de2)) * math.cos(r(ra1 - ra2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600


def main() -> None:
    probe = json.loads((PRIVATE / "swiss_probe.json").read_text())
    entries = probe["entries"]  # returned string -> [ra, dec, dist] (J2000 barycentric ICRS)
    want: dict[str, list[str]] = {}
    for ret in entries:
        name, desig = ret.split(",", 1)
        want[ret] = NAMED.get(name, []) + simbad_ids(desig)
    ids = sorted({i for v in want.values() for i in v})
    found: dict[str, dict] = {}
    for k in range(0, len(ids), 150):
        q = ",".join("'" + i.replace("'", "''") + "'" for i in ids[k:k + 150])
        for r in tap(f"SELECT i.id, {COLS} FROM ident AS i JOIN basic AS b ON i.oidref = b.oid WHERE i.id IN ({q})"):
            found[r["id"]] = r

    result: dict[str, dict] = {}
    for ret, v in entries.items():
        ra, de = v[0], v[1]
        name = ret.split(",", 1)[0]
        if name in DEFINED:
            result[ret] = ({"main_id": DEFINED[name], "oid": "-1"}, 0.0, "defined")
            continue
        if name in NO_SOURCE:
            result[ret] = None
            continue
        best = None
        for i in want[ret]:
            if i in found and found[i]["ra"]:
                r = found[i]
                d = sep_arcsec(ra, de, float(r["ra"]), float(r["dec"]))
                if best is None or d < best[1]:
                    best = (r, d, "designation")
        if best is None or best[1] > 10.0:
            # Nearest non-planet SIMBAD object within 3' of Swiss's place; used only if it
            # agrees within 10" (Swiss's label points elsewhere) or the designation failed.
            rows = tap(f"SELECT TOP 20 {COLS} FROM basic AS b "
                       f"WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {ra}, {de}, 0.05)) = 1")
            rows = [r for r in rows if r["ra"] and r["otype"] not in ("Pl", "Pl?")]
            if rows:
                r = min(rows, key=lambda r: sep_arcsec(ra, de, float(r["ra"]), float(r["dec"])))
                d = sep_arcsec(ra, de, float(r["ra"]), float(r["dec"]))
                if best is None or d < 10.0:
                    best = (r, d, "cone")
        result[ret] = best

    oids = sorted({b[0]["oid"] for b in result.values() if b and b[0]["oid"] != "-1"})
    hip: dict[str, int] = {}
    for k in range(0, len(oids), 300):
        q = ",".join(oids[k:k + 300])
        for r in tap(f"SELECT oidref, id FROM ident WHERE oidref IN ({q}) AND id LIKE 'HIP %'"):
            if r["id"].split()[1].isdigit():  # skip component ids such as "HIP 84345A"
                hip[r["oidref"]] = int(r["id"].split()[1])

    # Multiple-star objects without a HIP number of their own: use component A's HIP
    # (e.g. "* tet01 Eri" -> "* tet01 Eri A" = HIP 13847) so the name maps to the
    # Hipparcos star when Swiss's place is not the system place (the test checks it).
    sys_ids = sorted({b[0]["main_id"] for b in result.values()
                      if b and b[0]["oid"] != "-1" and b[0]["oid"] not in hip and b[0]["main_id"].startswith("* ")})
    comp_hip: dict[str, int] = {}
    for k in range(0, len(sys_ids), 150):
        q = ",".join("'" + (i + " A").replace("'", "''") + "'" for i in sys_ids[k:k + 150])
        rows = tap(f"SELECT i.id, h.id AS hip FROM ident AS i JOIN ident AS h ON h.oidref = i.oidref "
                   f"WHERE i.id IN ({q}) AND h.id LIKE 'HIP %'")
        for r in rows:
            if r["hip"].split()[1].isdigit():
                comp_hip[r["id"][:-2]] = int(r["hip"].split()[1])
    for ret, b in result.items():
        # Only where Swiss's place is not the system's own catalogue place (> 1").
        if b and b[0]["oid"] not in hip and b[0]["main_id"] in comp_hip and b[1] > 1.0:
            hip[b[0]["oid"]] = comp_hip[b[0]["main_id"]]

    lines = ["# Independent identification of every object behind a Swiss-accepted fixed-star name",
             "# (tools/identify_swiss_stars.py). Columns: returned string, SIMBAD main_id, HIP (0 = none),",
             "# method (designation | cone | defined | no-source), separation of SIMBAD's J2000 place from Swiss's (arcsec)."]
    n_bad = 0
    for ret in sorted(result, key=str.lower):
        b = result[ret]
        if b is None:
            n_bad += 1
            lines.append(f"{ret}\t-\t0\tno-source\t-")
        else:
            r, d, how = b
            lines.append(f"{ret}\t{r['main_id']}\t{hip.get(r['oid'], 0)}\t{how}\t{d:.3f}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    seps = sorted(b[1] for b in result.values() if b)
    print(f"{len(result)} entries; unidentified {n_bad}; cone {sum(1 for b in result.values() if b and b[2] == 'cone')}; "
          f"sep median {seps[len(seps) // 2]:.3f}\" p95 {seps[int(len(seps) * .95)]:.3f}\" max {seps[-1]:.1f}\"",
          file=sys.stderr)


if __name__ == "__main__":
    main()
