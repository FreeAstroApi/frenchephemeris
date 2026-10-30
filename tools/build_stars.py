#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Build the embedded fixed-star table `crates/fe-core/src/stars/catalog.rs`
(the public default catalogue) and `docs/STAR_NAMES.md`.

Every name and designation in the default catalogue comes from a public source:

* IAU Catalog of Star Names (IAU-CSN, WGSN), machine-readable edition by
  E. Mamajek (CC BY 4.0, IAU): proper names, Bayer/Flamsteed IDs, HIP numbers.
* SIMBAD (CDS) identifiers of each object: "NAME ..." proper names, Bayer and
  Flamsteed designations ("* alf Vir"), variable-star names, HR, HD, HIP,
  Messier/NGC/IC numbers.
* EXTRA_ALIASES: traditional names as given in R. H. Allen, *Star Names: Their
  Lore and Meaning* (1899; public domain) and common astrological usage.
* DEFINED: names of two conventional reference points.

Objects: the IAU-named stars, the alias targets, and `tools/star_extra_objects.tsv`
(HIP numbers / SIMBAD identifiers only, no names), plus SYSTEMS and DEFINED.

Astrometry: Hipparcos, the New Reduction (van Leeuwen 2007, A&A 474, 653; CDS
catalogue I/311, read from a VizieR mirror): ICRS, epoch J1991.25 (TT). Objects
without a Hipparcos 2007 solution use SIMBAD's ICRS J2000 place, proper motion,
parallax and RV (bibcodes in `rv_ref`), propagated to J1991.25 with ERFA
`pmsafe`; open clusters with Gaia-derived SIMBAD values use Dias et al. 2002
(CDS B/ocl). RV and Johnson V: SIMBAD (`basic`, `flux`, `mesVelocities`),
Hipparcos 1997 V (CDS I/239, field H5) as V fallback.

No Gaia-derived value (CC BY-SA 3.0 IGO) may enter the catalogue: see
GAIA_BIBCODES / is_gaia(); replaced values are logged in
tools/star_gaia_replacements.tsv.

The optional compatibility alias table (classic-software spellings) is built
separately by tools/build_star_compat.py into the private corpus; nothing here
reads it.

Usage:
    python tools/build_stars.py fetch   # download into data/sources/stars/ (small)
    python tools/build_stars.py build   # write catalog.rs and docs/STAR_NAMES.md
    python tools/build_stars.py         # both

Downloads are cached (incrementally), so `build` is deterministic for a given
cache. Python's urllib fails TLS on the dev machine: downloads use `curl`.
"""
from __future__ import annotations

import csv
import io
import json
import math
import re
import subprocess
import urllib.parse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "sources" / "stars"
OUT = ROOT / "crates" / "fe-core" / "src" / "stars" / "catalog.rs"
DOC = ROOT / "docs" / "STAR_NAMES.md"
EXTRA_OBJECTS_TSV = ROOT / "tools" / "star_extra_objects.tsv"

CSN_URL = "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt"
# ---------------------------------------------------------------------------
# No Gaia-derived values. Gaia data (DR1/DR2/EDR3/DR3 and their derived products)
# are released under CC BY-SA 3.0 IGO; none may enter the catalogue. A value is
# Gaia-derived when its SIMBAD/VizieR bibcode is a Gaia release (VizieR I/337,
# I/345-I/358), is listed below, or the publication title names Gaia or TGAS.
GAIA_VIZIER = re.compile(r"^\d{4}yCat\.13(37|4[5-9]|5[0-8])\.")
GAIA_BIBCODES = {
    "2016A&A...595A...1G", "2016A&A...595A...2G", "2018A&A...616A...1G", "2021A&A...649A...1G",
    "2023A&A...674A...1G", "2023A&A...674A...5K", "2019A&A...622A.205K", "2018A&A...616A...5S",
    "2018A&A...616A...6K", "2018A&A...616A...7S", "2018A&A...616A..10G", "2018A&A...618A..93C",
    "2020A&A...640A...1C", "2021A&A...647A..19T", "2019MNRAS.482.5138B", "2019MNRAS.483.5026L",
    "2019ApJ...870...32K", "2021MNRAS.504..356D", "2021MNRAS.505.5978V", "2023A&A...674A..37G",
    "2023A&A...673A.114H", "2021A&A...649A...6G", "2021A&A...649A...7G",
    # Excluded out of caution: catalogues that merge Gaia values with their own.
    "2021MNRAS.506..150B",  # GALAH DR3
    "2018MNRAS.473..849D",  # cluster distances partly based on Gaia
}
_GAIA_TITLE = re.compile(r"\bGaia\b|\bTGAS\b", re.I)


def is_gaia(bibcode: str, titles: dict) -> bool:
    """True if the value with this bibcode is (or is derived from) Gaia data."""
    if not bibcode:
        return False
    b = bibcode.strip()
    return b in GAIA_BIBCODES or bool(GAIA_VIZIER.match(b)) or bool(_GAIA_TITLE.search(titles.get(b) or ""))


# Preferred radial-velocity compilations (non-Gaia) when SIMBAD's adopted value is Gaia.
RV_PREFERRED = ["2006AstL...32..759G", "2007AN....328..889K", "2005A&A...430..165F", "2004A&A...418..989N",
                "2004A&A...424..727P", "1979IAUS...30...57E", "1953GCRV..C......0W"]
DIAS_OCL = "2002A&A...389..871D"  # Dias et al. 2002 (CDS B/ocl, version 3.5, 2015)

# VizieR (CDS) mirror serving the original CDS catalogues (I/311 Hipparcos 2007,
# I/239 Hipparcos 1997, B/ocl open clusters) in FITS binary form at full precision.
VIZIER = "http://vizier.idia.ac.za/viz-bin"
SIMBAD_TAP = "https://simbad.harvard.edu/simbad/sim-tap/sync"

# Hipparcos catalogue epoch J1991.25 (TT); SIMBAD places are epoch J2000.
HIP_EPOCH_JD = 2448349.0625
J2000 = 2451545.0
# Magnitude returned for objects without a V magnitude (the classic API's sentinel).
NO_MAG = 999.99


# Record key of the alpha Cen AB barycentre (no HIP number; hip = 0 in the table).
ALF_CEN_AB = "SYS alf Cen AB"

# Traditional names (R. H. Allen, Star Names: Their Lore and Meaning, 1899, public
# domain; astrological usage), lower precedence than IAU WGSN names.
# (alias, record key: HIP int, SIMBAD main_id or ALF_CEN_AB).
EXTRA_ALIASES = [
    ("Agena", 68702),            # beta Cen (WGSN: Hadar)
    ("Rigil Kentaurus", ALF_CEN_AB),
    ("Toliman", ALF_CEN_AB),     # astrological usage: alpha Cen; WGSN gives Toliman to B
    ("Bungula", ALF_CEN_AB),
    ("Pushya", 42911),           # delta Cnc, yogatara of the Pushya nakshatra
    ("Benetnasch", 67301),       # eta UMa = Alkaid
    ("Alnath", 25428),           # beta Tau = Elnath
    ("Wega", 91262),             # Vega
    ("Gemma", 76267),            # alpha CrB = Alphecca
    ("Sirrah", 677),             # alpha And = Alpheratz
    ("Graffias", 78820),         # beta Sco (astrological usage)
    ("Zuben Elschemali", 74785),
    ("Kitalpha", 104987),        # alpha Equ
    ("Cor Scorpii", 80763),      # Antares
    ("Cor Leonis", 49669),       # Regulus
    ("Cor Hydrae", 46390),       # Alphard
    ("Algenib Persei", 15863),   # alpha Per = Mirfak
    ("Mirphak", 15863),
    ("Deneb Adige", 102098),     # alpha Cyg
    ("Arided", 102098),
    ("Kaus Meridionalis", 89931),
    ("Galactic Center", "NAME Sgr A*"),
    ("Sgr A*", "NAME Sgr A*"),
]

# alpha Cen A and B: Hipparcos proper motions are instantaneous and include
# ~0.2"/yr of orbital motion (P = 79.9 yr), a ~450" error after 2000 years.
# Both components get the system motion of the barycentre record, so each
# keeps its J1991.25 offset from the barycentre (orbital amplitude <= ~10").
_SYS = {"pm_ra": -3608.0, "pm_de": 686.0, "plx": 750.81, "rv": -22.3,
        "rv_ref": "alpha Cen AB system motion: SIMBAD pm/RV (1979IAUS...30...57E), 2021AJ....162...14A plx"}
OVERRIDES = {71683: _SYS, 71681: _SYS}

# Star systems represented by their barycentre: place at J1991.25 = mass-weighted
# mean of the Hipparcos 2007 component places, motion = the system's.
# alpha Cen AB: masses 1.1055 / 0.9373 Msun (Kervella et al. 2016, A&A 594, A107);
# system proper motion and RV from SIMBAD "alf Cen" (RV 1979IAUS...30...57E);
# parallax Akeson et al. 2021 (AJ 162, 14); V of the pair combined from A and B.
SYSTEMS = {
    ALF_CEN_AB: {"components": (71683, 71681), "masses": (1.1055, 0.9373), "nom": "alCen",
                 "desig": "alf Cen AB", "vmag": -0.27, "vmag_ref": "combined A (V=-0.01) + B (V=1.33)",
                 "pm_ra": -3608.0, "pm_de": 686.0, "plx": 750.81, "rv": -22.3,
                 "rv_ref": "alpha Cen AB barycentre: Kervella+2016 masses, SIMBAD pm/RV "
                 "(1979IAUS...30...57E), 2021AJ....162...14A plx"},
}
# SIMBAD objects that are represented by a SYSTEMS record.
SYSTEM_OF = {"* alf Cen": ALF_CEN_AB, "* alf Cen A": ALF_CEN_AB}

# Conventional reference points (no catalogue object): defining coordinates.
#  Great Attractor: centre l = 325.3 deg, b = -7.2 deg (the value conventionally
#    quoted for the GA overdensity, cf. Kraan-Korteweg et al. 1996, Nature 379, 519);
#    distance 4350 km/s (Lynden-Bell et al. 1988, ApJ 326, 19) / H0 = 70.
#  Virgo Cluster centre: 12h24m, +13 deg 00' (FK4, B1950), the round-number centre used
#    in Local Supercluster work; distance 16.5 Mpc (Mei et al. 2007, ApJ 655, 144).
DEFINED = {
    "DEF Great Attractor (l=325.3, b=-7.2)": {"name": "Great Attractor", "gal": (325.3, -7.2), "dist_mpc": 4350 / 70,
        "ref": "defined point: GA centre l=325.3 b=-7.2 (cf. 1996Natur.379..519K); dist 1988ApJ...326...19L, H0=70"},
    "DEF Virgo Cluster centre (B1950 12h24m +13d00m)": {"name": "Virgo Cluster", "fk4_b1950": (186.0, 13.0),
        "dist_mpc": 16.5,
        "ref": "defined point: Virgo centre 12h24m +13d00m (FK4 B1950); dist 2007ApJ...655..144M"},
}

# Sgr A* (Galactic Center): SIMBAD place; apparent proper motion of Reid &
# Brunthaler 2020 (ApJ 892, 39): mu_l = -6.411, mu_b = -0.219 mas/yr (converted to
# ICRS); distance R0 = 8.178 kpc (GRAVITY Collaboration 2019, A&A 625, L10).
SGR_A = "NAME Sgr A*"
SGR_A_PM_GAL = (-6.411, -0.219)
SGR_A_PLX = 1.0 / 8.178  # mas


def curl(url: str, data: dict | None = None, timeout: int = 180) -> bytes:
    cmd = ["curl", "-sS", "-f", "-m", str(timeout), "-L", url]
    if data is not None:
        for k, v in data.items():
            cmd += ["--data-urlencode", f"{k}={v}"]
    for attempt in range(4):  # TAP mirrors occasionally return transient HTTP errors
        try:
            return subprocess.run(cmd, check=True, capture_output=True).stdout
        except subprocess.CalledProcessError:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    raise AssertionError


def tap(url: str, query: str) -> list[dict]:
    raw = curl(url, {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query})
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def parse_csn(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "$")):
            continue
        rest = line[81:].split()
        if rest[2] == "_":
            continue
        ident, con, comp = line[49:55].strip(), line[61:65].strip(), line[65:70].strip()
        rows.append({
            "name": line[0:18].strip(),
            "id": "" if ident == "_" else ident,
            "con": "" if con == "_" else con,
            "comp": "" if comp == "_" else comp,
            "hip": int(rest[2]),
        })
    return rows


def parse_tsv(path: Path) -> list[list[str]]:
    return [l.rstrip("\n").split("\t") for l in path.read_text(encoding="utf-8").splitlines()
            if l.strip() and not l.startswith("#")]


def extra_objects() -> list:
    """tools/star_extra_objects.tsv: HIP numbers (int) and SIMBAD main_ids (str)."""
    out = []
    for r in parse_tsv(EXTRA_OBJECTS_TSV):
        out.append(int(r[1]) if r[0] == "HIP" else r[1])
    return out


def wanted() -> tuple[list[int], list[str]]:
    hips = {r["hip"] for r in parse_csn(CACHE / "IAU-CSN.txt")}
    hips |= {h for _, h in EXTRA_ALIASES if isinstance(h, int)}
    for sy in SYSTEMS.values():
        hips |= set(sy["components"])
    objs = {h for _, h in EXTRA_ALIASES if isinstance(h, str) and h not in SYSTEMS}
    for k in extra_objects():
        (hips if isinstance(k, int) else objs).add(k)
    return sorted(hips), sorted(objs)


SIMBAD_COLS = ("b.main_id, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value, b.rvz_radvel, b.coo_bibcode, "
               "b.pm_bibcode, b.plx_bibcode, b.rvz_bibcode, b.otype")


def fetch_incremental(path: Path, keys: list, key_of, query_chunk, chunk: int) -> None:
    """Fetch rows for `keys` not yet cached. Keys with no row get a {"_missing": key} marker."""
    rows = json.loads(path.read_text()) if path.exists() else []
    have = {r["_missing"] if "_missing" in r else key_of(r) for r in rows}
    todo = [k for k in keys if k not in have]
    for i in range(0, len(todo), chunk):
        got = query_chunk(todo[i:i + chunk])
        found = {key_of(r) for r in got}
        rows += got + [{"_missing": k} for k in todo[i:i + chunk] if k not in found]
    if todo or not path.exists():
        path.write_text(json.dumps(rows, indent=0))


def load_rows(name: str) -> list[dict]:
    return [r for r in json.loads((CACHE / name).read_text()) if "_missing" not in r]


def sq(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def vizier_fits(source: str, cols: list[str], key: str, values: list) -> list[dict]:
    """Rows of a CDS catalogue from the VizieR mirror (FITS binary table: full precision)."""
    import struct
    url = (f"{VIZIER}/asu-binfits?-source={source}&-out={','.join(cols)}&-out.max=unlimited"
           f"&{key}={','.join(map(str, values))}")
    d = curl(url)
    start = d.find(b"XTENSION")
    if start < 0:
        return []
    h = d[start:]
    kv, n = {}, 0
    while True:
        c = h[n * 80:(n + 1) * 80].decode("ascii")
        if c.startswith("END "):
            break
        if c[8:10] == "= ":
            kv[c[:8].strip()] = c[10:].split("/")[0].strip().strip("'").strip()
        n += 1
    hdr = ((n + 1) * 80 + 2879) // 2880 * 2880
    nrow, width, nf = int(kv["NAXIS2"]), int(kv["NAXIS1"]), int(kv["TFIELDS"])
    names = [kv[f"TTYPE{i}"] for i in range(1, nf + 1)]
    codes = {"D": "d", "E": "f", "J": "i", "I": "h", "K": "q"}
    fmt = ">" + "".join(codes[kv[f"TFORM{i}"][-1]] for i in range(1, nf + 1))
    body = h[hdr:hdr + nrow * width]
    return [dict(zip(names, struct.unpack(fmt, body[r * width:(r + 1) * width]))) for r in range(nrow)]


def fetch() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    csn = CACHE / "IAU-CSN.txt"
    if not csn.exists():
        csn.write_bytes(curl(CSN_URL))
    hips, objs = wanted()
    # Hipparcos 2007 (CDS I/311/hip2). VizieR serves RArad/DErad in degrees as doubles;
    # the catalogue values are radians to 1e-10, restored by rounding.
    fetch_incremental(
        CACHE / "hip2.json", hips, lambda r: int(r["hip"]),
        lambda c: [{"hip": str(r["HIP"]), "ra_rad": f"{math.radians(r['RArad']):.10f}",
                    "de_rad": f"{math.radians(r['DErad']):.10f}", "plx": f"{r['Plx']:.2f}",
                    "pm_ra": f"{r['pmRA']:.2f}", "pm_de": f"{r['pmDE']:.2f}"}
                   for r in vizier_fits("I/311/hip2", ["HIP", "RArad", "DErad", "Plx", "pmRA", "pmDE"], "HIP", c)],
        200)
    # SIMBAD basic data for every HIP star (RV; fallback astrometry) and every non-HIP object.
    fetch_incremental(
        CACHE / "simbad_hip.json", hips, lambda r: int(r["id"].split()[1]),
        lambda c: tap(SIMBAD_TAP, f"SELECT i.id, b.oid, {SIMBAD_COLS} FROM ident AS i JOIN basic AS b "
                                  f"ON i.oidref = b.oid WHERE i.id IN ({','.join(sq(f'HIP {h}') for h in c)})"), 200)
    fetch_incremental(
        CACHE / "simbad_obj.json", objs, lambda r: r["main_id"],
        lambda c: tap(SIMBAD_TAP, f"SELECT b.oid, {SIMBAD_COLS} FROM basic AS b "
                                  f"WHERE b.main_id IN ({','.join(sq(o) for o in c)})"), 200)
    # Hipparcos Catalogue (ESA 1997) V magnitude (field H5): fallback where SIMBAD has no V.
    fetch_incremental(
        CACHE / "hip1_v.json", hips, lambda r: int(r["hip"]),
        lambda c: [{"hip": str(r["HIP"]), "vmag": f"{r['Vmag']:.2f}"}
                   for r in vizier_fits("I/239/hip_main", ["HIP", "Vmag"], "HIP", c)], 200)
    # Multiple-star objects SIMBAD lists without motion (e.g. "* 61 Cyg"): motion of component A.
    need = sorted(r["main_id"] for r in load_rows("simbad_obj.json") if not r["pmra"] and r["otype"].startswith("**"))
    fetch_incremental(
        CACHE / "simbad_comp.json", need, lambda r: r["system"],
        lambda c: [dict(r, system=r["id"][:-2]) for r in
                   tap(SIMBAD_TAP, f"SELECT i.id, b.oid, {SIMBAD_COLS} FROM ident AS i JOIN basic AS b "
                                   f"ON i.oidref = b.oid WHERE i.id IN ({','.join(sq(o + ' A') for o in c)})")], 200)
    # SIMBAD identifiers (names, Bayer/Flamsteed, variable, HR, HD, HIP, catalogue numbers), by oid.
    oids_all = sorted({int(r["oid"]) for f in ("simbad_hip.json", "simbad_obj.json") for r in load_rows(f)})
    fetch_incremental(
        CACHE / "simbad_ids.json", oids_all, lambda r: int(r["oidref"]),
        lambda c: tap(SIMBAD_TAP, "SELECT oidref, id FROM ident WHERE oidref IN "
                                  f"({','.join(map(str, c))}) AND (id LIKE 'NAME %' OR id LIKE '* %' OR id LIKE 'V* %' "
                                  "OR id LIKE 'HR %' OR id LIKE 'HD %' OR id LIKE 'HIP %' OR id LIKE 'M %' "
                                  "OR id LIKE 'NGC %' OR id LIKE 'IC %')"), 150)
    # Johnson V magnitudes (SIMBAD flux table), by SIMBAD oid.
    oids = sorted({int(r["oid"]) for f in ("simbad_hip.json", "simbad_obj.json", "simbad_comp.json")
                   for r in load_rows(f)})
    fetch_incremental(
        CACHE / "simbad_v.json", oids, lambda r: int(r["oidref"]),
        lambda c: tap(SIMBAD_TAP, "SELECT oidref, flux, bibcode FROM flux WHERE filter = 'V' AND oidref IN "
                                  f"({','.join(map(str, c))})"), 300)
    fetch_provenance()


def all_bibcodes() -> set[str]:
    out = set()
    for f, cols in (("simbad_hip.json", ("coo_bibcode", "pm_bibcode", "plx_bibcode", "rvz_bibcode")),
                    ("simbad_obj.json", ("coo_bibcode", "pm_bibcode", "plx_bibcode", "rvz_bibcode")),
                    ("simbad_comp.json", ("coo_bibcode", "pm_bibcode", "plx_bibcode", "rvz_bibcode")),
                    ("simbad_v.json", ("bibcode",)), ("mes_vel.json", ("bibcode",)), ("mes_pm.json", ("bibcode",)),
                    ("mes_plx.json", ("bibcode",)), ("mes_dist.json", ("bibcode",))):
        if (CACHE / f).exists():
            out |= {r[c].strip() for r in load_rows(f) for c in cols if r.get(c)}
    return out


def fetch_titles() -> dict:
    fetch_incremental(
        CACHE / "ref_titles.json", sorted(all_bibcodes()), lambda r: r["bibcode"],
        lambda c: tap(SIMBAD_TAP, f"SELECT bibcode, title FROM ref WHERE bibcode IN ({','.join(sq(b) for b in c)})"),
        150)
    return {r["bibcode"]: r["title"] for r in load_rows("ref_titles.json")}


def fetch_provenance() -> None:
    """Publication titles of every bibcode (for the Gaia rule) and the non-Gaia
    alternatives (SIMBAD measurement tables, Dias 2002 clusters) where needed."""
    titles = fetch_titles()
    g = lambda b: is_gaia(b, titles)  # noqa: E731
    rows = load_rows("simbad_hip.json") + load_rows("simbad_obj.json") + load_rows("simbad_comp.json")
    need_rv = sorted({int(r["oid"]) for r in rows if g(r["rvz_bibcode"])})
    need_ast = sorted({int(r["oid"]) for r in rows if g(r["pm_bibcode"]) or g(r["plx_bibcode"])})
    fetch_incremental(
        CACHE / "mes_vel.json", need_rv, lambda r: int(r["oidref"]),
        lambda c: tap(SIMBAD_TAP, "SELECT oidref, veltype, velvalue, quality, nbmes, bibcode FROM mesVelocities "
                                  f"WHERE oidref IN ({','.join(map(str, c))})"), 100)
    for name, cols, table in (("mes_pm.json", "pmra, pmde, pmra_err, bibcode", "mesPM"),
                              ("mes_plx.json", "plx, plx_err, bibcode", "mesPLX"),
                              ("mes_dist.json", "dist, unit, qual, bibcode", "mesDistance")):
        fetch_incremental(
            CACHE / name, need_ast, lambda r: int(r["oidref"]),
            lambda c, cols=cols, table=table: tap(SIMBAD_TAP, f"SELECT oidref, {cols} FROM {table} "
                                                              f"WHERE oidref IN ({','.join(map(str, c))})"), 100)
    # Open clusters whose SIMBAD place/motion is Gaia-derived: Dias et al. 2002 (B/ocl).
    clusters = {r["main_id"]: r for r in rows if r["otype"].strip() in ("OpC", "Cl*") and g(r["coo_bibcode"])}

    def dias_rows(m: str) -> list[dict]:
        r = clusters[m]
        dec = float(r["dec"])
        near = vizier_tsv("B/ocl/clusters", "Cluster,_RAJ2000,_DEJ2000,Dist,pmRA,pmDE,r_pmRA,RV,r_RV",
                          {"-c": f"{float(r['ra'])}{'+' if dec >= 0 else ''}{dec}", "-c.rm": "30"})
        return [dict(x, main_id=m) for x in near if " ".join(x["Cluster"].split()) == " ".join(m.split())]

    fetch_incremental(CACHE / "dias_ocl.json", sorted(clusters), lambda r: r["main_id"],
                      lambda c: [x for m in c for x in dias_rows(m)], 1)
    if not (CACHE / "dias_refs.json").exists():
        (CACHE / "dias_refs.json").write_text(json.dumps(vizier_tsv("B/ocl/refs", "Ref,BibCode", {}), indent=0))
    fetch_titles()


def vizier_tsv(source: str, cols: str, constraints: dict) -> list[dict]:
    q = "&".join(f"{k}={urllib.parse.quote(v)}" for k, v in constraints.items())
    raw = curl(f"{VIZIER}/asu-tsv?-source={source}&-oc.form=d&-out={cols}&-out.max=unlimited&{q}").decode("utf-8")
    lines = [l for l in raw.splitlines() if l and not l.startswith("#")]
    if len(lines) < 3:
        return []
    head = [h.strip() for h in lines[0].split("\t")]
    return [dict(zip(head, (x.strip() for x in l.split("\t")))) for l in lines[3:]]


def nomenclature(ident: str, con: str) -> tuple[str, str]:
    """("", IAU 'alf Vir' form) for a WGSN ID + constellation (Flamsteed numbers without zeros)."""
    if not ident or not con:
        return "", ""
    if ident.isdigit():  # Flamsteed number
        return "", f"{int(ident)} {con}"
    return "", f"{ident} {con}"


def fold(s: str) -> str:
    """Lookup key: lower case, whitespace removed (as the classic API matches names)."""
    return "".join(s.lower().split())


def rs_str(s: str) -> str:
    return json.dumps(s, ensure_ascii=False)


def rs_f64(x: float) -> str:
    return repr(float(x))


def unit(ra: float, dec: float) -> list[float]:
    return [math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)]


def num(x) -> float | None:
    return float(x) if x not in ("", None) else None


MAS_RAD = math.pi / 180 / 3.6e6


def gal_pm_to_icrs(ra: float, dec: float, mu_l: float, mu_b: float) -> tuple[float, float]:
    """Galactic proper motion (mas/yr) -> (mu_alpha*, mu_delta) mas/yr at an ICRS place (rad)."""
    import erfa
    l, b = erfa.icrs2g(ra, dec)
    h = 1000.0  # years: a 1000-yr step keeps the difference well above rounding
    r2, d2 = erfa.g2icrs(l + mu_l * MAS_RAD * h / math.cos(b), b + mu_b * MAS_RAD * h)
    dra = (r2 - ra + math.pi) % (2 * math.pi) - math.pi
    return (dra * math.cos(dec) / MAS_RAD / h, (d2 - dec) / MAS_RAD / h)


def to_hip_epoch(ra, dec, pmra, pmde, plx, rv) -> dict:
    """ICRS J2000 place + motion -> the same star at J1991.25 (ERFA pmsafe)."""
    import erfa
    pmr = pmra * MAS_RAD / math.cos(dec)
    pmd = pmde * MAS_RAD
    if pmra or pmde:
        ra, dec, pmr, pmd, _, rv = erfa.pmsafe(ra, dec, pmr, pmd, plx / 1000, rv, J2000, 0.0, HIP_EPOCH_JD, 0.0)
    return {"ra_rad": float(ra) % (2 * math.pi), "de_rad": float(dec),
            "pm_ra": float(pmr * math.cos(dec) / MAS_RAD), "pm_de": float(pmd / MAS_RAD), "plx": plx, "rv": float(rv)}


def simbad_astrometry(r: dict, comp: dict | None = None) -> tuple[dict, str]:
    ra, dec = math.radians(float(r["ra"])), math.radians(float(r["dec"]))
    if comp is not None and not r["pmra"]:
        # The system's place with the motion of its component A (common proper motion pair).
        r = dict(r, **{k: comp[k] for k in ("pmra", "pmdec", "plx_value", "rvz_radvel", "pm_bibcode",
                                             "plx_bibcode", "rvz_bibcode")})
        r["coo_bibcode"] = f"{r['coo_bibcode']} (motion of {comp['main_id'].strip()})"
    pmra, pmde = num(r["pmra"]) or 0.0, num(r["pmdec"]) or 0.0
    plx, rv = num(r["plx_value"]), num(r["rvz_radvel"]) or 0.0
    note = f"SIMBAD {r['main_id'].strip()}: place {r['coo_bibcode'] or '?'}"
    if r["main_id"] == SGR_A:
        pmra, pmde = gal_pm_to_icrs(ra, dec, *SGR_A_PM_GAL)
        plx = SGR_A_PLX
        note += "; pm 2020ApJ...892...39R; dist 2019A&A...625L..10G"
    else:
        if r["pm_bibcode"] and (pmra or pmde):
            note += f"; pm {r['pm_bibcode']}"
        if plx is not None and plx > 0:
            note += f"; plx {r['plx_bibcode']}"
    if plx is None or plx <= 0:
        plx = 0.0
        note += "; no parallax"
    if r["rvz_bibcode"]:
        note += f"; RV {r['rvz_bibcode']}"
    return to_hip_epoch(ra, dec, pmra, pmde, plx, rv), note


def defined_astrometry(d: dict) -> dict:
    import erfa
    if "gal" in d:
        ra, dec = erfa.g2icrs(math.radians(d["gal"][0]), math.radians(d["gal"][1]))
    else:
        ra, dec = erfa.fk45z(math.radians(d["fk4_b1950"][0]), math.radians(d["fk4_b1950"][1]), 1950.0)
    return {"ra_rad": float(ra) % (2 * math.pi), "de_rad": float(dec), "pm_ra": 0.0, "pm_de": 0.0,
            "plx": 1e-3 / d["dist_mpc"], "rv": 0.0}


def record_line(s: dict, h: int, a: dict, rv_kms: float, ref: str) -> str:
    return (
        "    StarRecord {{ id: {}, name: {}, nomenclature: {}, designation: {}, hip: {}, "
        "ra_rad: {}, dec_rad: {}, pmra_mas_yr: {}, pmdec_mas_yr: {}, plx_mas: {}, "
        "rv_km_s: {}, vmag: {}, rv_ref: {} }},".format(
            rs_str(s["id"]), rs_str(s["name"]), rs_str(s["nom"]), rs_str(s["desig"]), h,
            rs_f64(a["ra_rad"]), rs_f64(a["de_rad"]), rs_f64(a["pm_ra"]), rs_f64(a["pm_de"]),
            rs_f64(a["plx"]), rs_f64(rv_kms), rs_f64(s["vmag"]), rs_str(ref),
        )
    )


GREEK_UNICODE = {
    "alf": "α", "bet": "β", "gam": "γ", "del": "δ", "eps": "ε", "zet": "ζ", "eta": "η", "tet": "θ",
    "iot": "ι", "kap": "κ", "lam": "λ", "mu": "μ", "nu": "ν", "ksi": "ξ", "omi": "ο", "pi": "π",
    "rho": "ρ", "sig": "σ", "tau": "τ", "ups": "υ", "phi": "φ", "chi": "χ", "psi": "ψ", "ome": "ω",
}
_GREEK_RE = re.compile(r"^(alf|bet|gam|del|eps|zet|eta|tet|iot|kap|lam|mu|nu|ksi|omi|pi|rho|sig|tau|ups|phi|chi|psi|ome)"
                       r"(\d*) (.+)$")


def norm_id(i: str) -> str:
    """SIMBAD identifier -> standard form: "* alf Vir" -> "alf Vir", "*  61 Cyg" -> "61 Cyg",
    "* mu. Cas" -> "mu Cas", "V* V645 Cen" -> "V645 Cen", "M  31" -> "M 31", "NAME Sgr A*" -> "Sgr A*"."""
    i = " ".join(i.split())
    for pre in ("* ", "V* ", "NAME "):
        if i.startswith(pre):
            i = i[len(pre):]
            break
    return re.sub(r"\b(mu|nu|pi)\.", r"\1", i)


def desig_variants(d: str) -> list[str]:
    """Accepted spellings of a designation: as is, Greek letter in Unicode, superscript
    index without leading zeros ("alf01 Cap" -> "alf1 Cap", "α1 Cap")."""
    out = [d]
    m = _GREEK_RE.match(d)
    if m:
        g, n, rest = m.groups()
        n2 = n.lstrip("0")
        for nn in {n, n2}:
            out += [f"{g}{nn} {rest}", f"{GREEK_UNICODE[g]}{nn} {rest}"]
    return list(dict.fromkeys(out))


def id_rank(i: str) -> int:
    """Preference of a SIMBAD identifier as the displayed designation."""
    n = norm_id(i)
    if i.startswith("* ") and _GREEK_RE.match(n):
        return 0
    if i.startswith("* ") and n.split()[0].isdigit():
        return 1
    if i.startswith("* "):
        return 2
    if i.startswith("V* "):
        return 3
    if i.startswith(("M ", "NGC ", "IC ")):
        return 4
    return {"HR": 5, "HD": 6, "HIP": 7}.get(i.split()[0], 9)


def record_id(k) -> str:
    """Stable identifier of a record (used by the optional compatibility table)."""
    if isinstance(k, int):
        return f"HIP {k}"
    if k in SYSTEMS or k in DEFINED:
        return k
    return "SIMBAD " + " ".join(k.split())


def build() -> None:
    csn = parse_csn(CACHE / "IAU-CSN.txt")
    hip2 = {int(r["hip"]): r for r in load_rows("hip2.json")}
    sim_hip = {int(r["id"].split()[1]): r for r in load_rows("simbad_hip.json")}
    sim_obj = {r["main_id"]: r for r in load_rows("simbad_obj.json")}
    vflux = {int(r["oidref"]): (float(r["flux"]), r["bibcode"]) for r in load_rows("simbad_v.json") if r["flux"]}
    comps = {r["system"]: r for r in load_rows("simbad_comp.json")}
    hip1_v = {int(r["hip"]): float(r["vmag"]) for r in load_rows("hip1_v.json") if r["vmag"]}
    titles = {r["bibcode"]: r["title"] for r in load_rows("ref_titles.json")}
    g = lambda b: is_gaia(b, titles)  # noqa: E731
    mes_vel: dict[int, list] = {}
    for r in load_rows("mes_vel.json"):
        mes_vel.setdefault(int(r["oidref"]), []).append(r)
    mes_pm: dict[int, list] = {}
    for r in load_rows("mes_pm.json"):
        mes_pm.setdefault(int(r["oidref"]), []).append(r)
    mes_plx: dict[int, list] = {}
    for r in load_rows("mes_plx.json"):
        mes_plx.setdefault(int(r["oidref"]), []).append(r)
    mes_dist: dict[int, list] = {}
    for r in load_rows("mes_dist.json"):
        mes_dist.setdefault(int(r["oidref"]), []).append(r)
    dias = {r["main_id"]: r for r in load_rows("dias_ocl.json")}
    dias_refs = {r["Ref"].strip(): r["BibCode"].strip() for r in json.loads((CACHE / "dias_refs.json").read_text())}
    changes: list[tuple[str, str, str, str]] = []  # (object, field, old bibcode, new source or "dropped")

    def year(b: str) -> int:
        return int(b[:4]) if b[:4].isdigit() else 0

    def replace_rv(r: dict) -> dict:
        if not g(r.get("rvz_bibcode")):
            return r
        qrank = {q: i for i, q in enumerate("ABCDE")}
        cands = [m for m in mes_vel.get(int(r["oid"]), []) if m["veltype"].strip().lower() == "v"
                 and m["velvalue"] and not g(m["bibcode"]) and len(m["bibcode"].strip()) == 19]
        cands.sort(key=lambda m: (qrank.get(m["quality"].strip(), 9),
                                  RV_PREFERRED.index(m["bibcode"]) if m["bibcode"] in RV_PREFERRED else 99,
                                  -year(m["bibcode"])))
        old = r["rvz_bibcode"]
        if cands:
            m = cands[0]
            changes.append((r["main_id"].strip(), "rv", old, m["bibcode"]))
            return dict(r, rvz_radvel=f"{float(m['velvalue']):.3f}".rstrip("0").rstrip("."), rvz_bibcode=m["bibcode"])
        changes.append((r["main_id"].strip(), "rv", old, "dropped"))
        return dict(r, rvz_radvel="", rvz_bibcode="")

    def replace_astrometry(r: dict, motion_only: bool = False) -> dict:
        """Non-Gaia place/pm/plx for a SIMBAD row (clusters: Dias 2002). `motion_only`:
        the row's place is not used (component rows lend only their motion)."""
        r = dict(r)
        oid = int(r["oid"])
        if r["main_id"] in dias and (g(r["coo_bibcode"]) or g(r["pm_bibcode"]) or g(r["plx_bibcode"])):
            d = dias[r["main_id"]]
            for f, old in (("place", r["coo_bibcode"]), ("pm", r["pm_bibcode"]), ("plx", r["plx_bibcode"])):
                if old:
                    changes.append((r["main_id"].strip(), f, old, DIAS_OCL))
            r.update(ra=d["_RAJ2000"], dec=d["_DEJ2000"], coo_bibcode=DIAS_OCL,
                     pmra=d["pmRA"], pmdec=d["pmDE"], pm_bibcode=dias_refs.get(d["r_pmRA"], DIAS_OCL),
                     plx_value=f"{1000.0 / float(d['Dist']):.4f}" if d["Dist"] else "", plx_bibcode=DIAS_OCL)
            if g(r["rvz_bibcode"]) and d["RV"]:
                changes.append((r["main_id"].strip(), "rv", r["rvz_bibcode"], dias_refs.get(d["r_RV"], DIAS_OCL)))
                r.update(rvz_radvel=d["RV"], rvz_bibcode=dias_refs.get(d["r_RV"], DIAS_OCL))
            return r
        if motion_only:
            r.update(ra="", dec="", coo_bibcode="")
        assert not g(r["coo_bibcode"]), f"{r['main_id']}: Gaia-derived place and no non-Gaia source"
        if g(r["pm_bibcode"]):
            c = sorted((m for m in mes_pm.get(oid, []) if not g(m["bibcode"]) and m["pmra"] and m["pmde"]),
                       key=lambda m: (float(m["pmra_err"]) if m["pmra_err"] else 99.0, -year(m["bibcode"])))
            new = c[0]["bibcode"] if c else "dropped"
            changes.append((r["main_id"].strip(), "pm", r["pm_bibcode"], new))
            r.update(pmra=c[0]["pmra"] if c else "", pmdec=c[0]["pmde"] if c else "",
                     pm_bibcode=c[0]["bibcode"] if c else "")
        if g(r["plx_bibcode"]):
            c = sorted((m for m in mes_plx.get(oid, []) if not g(m["bibcode"]) and m["plx"] and float(m["plx"]) > 0),
                       key=lambda m: (float(m["plx_err"]) if m["plx_err"] else 99.0, -year(m["bibcode"])))
            if c:
                changes.append((r["main_id"].strip(), "plx", r["plx_bibcode"], c[0]["bibcode"]))
                r.update(plx_value=c[0]["plx"], plx_bibcode=c[0]["bibcode"])
            else:
                scale = {"pc": 1.0, "kpc": 1e3, "Mpc": 1e6}
                d = sorted((m for m in mes_dist.get(oid, []) if not g(m["bibcode"]) and m["dist"]
                            and m["unit"].strip() in scale), key=lambda m: -year(m["bibcode"]))
                if d:
                    pc = float(d[0]["dist"]) * scale[d[0]["unit"].strip()]
                    changes.append((r["main_id"].strip(), "plx", r["plx_bibcode"], d[0]["bibcode"] + " (distance)"))
                    r.update(plx_value=f"{1000.0 / pc:.6g}", plx_bibcode=d[0]["bibcode"])
                else:
                    changes.append((r["main_id"].strip(), "plx", r["plx_bibcode"], "dropped"))
                    r.update(plx_value="", plx_bibcode="")
        return replace_rv(r)

    sim_hip = {k: replace_rv(v) for k, v in sim_hip.items()}
    wanted_objs = set(wanted()[1])
    sim_obj = {k: replace_astrometry(v) for k, v in sim_obj.items() if k in wanted_objs}
    # Component-A motion for systems (e.g. 61 Cyg A): the most precise non-Gaia measurement.
    comps = {k: replace_astrometry(c, motion_only=True) for k, c in comps.items()}
    ids_by_oid: dict[int, list[str]] = {}
    for r in load_rows("simbad_ids.json"):
        ids_by_oid.setdefault(int(r["oidref"]), []).append(r["id"])

    # Records keyed by HIP (int), SIMBAD main_id or DEFINED key (str), or a SYSTEMS key.
    recs: dict = {}

    def ensure(key) -> bool:
        if key in recs:
            return True
        ok = key in SYSTEMS or key in DEFINED or (key in hip2 or key in sim_hip if isinstance(key, int)
                                                  else key in sim_obj)
        if not ok:
            print(f"warning: {key}: no Hipparcos 2007 / SIMBAD data, skipped", file=sys.stderr)
            return False
        recs[key] = {"names": [], "desigs": []}
        return True

    wgsn = []
    for r in csn:
        if not ensure(r["hip"]):
            continue
        rec = recs[r["hip"]]
        if "wgsn_desig" not in rec:
            rec["wgsn_desig"] = nomenclature(r["id"], r["con"])[1]
            rec["comp"] = r["comp"]
        rec["names"].append((r["name"], "IAU-CSN"))
        wgsn.append((r["name"], r["hip"]))
    for a, k in EXTRA_ALIASES:
        if ensure(k):
            recs[k]["names"].append((a, "Allen 1899 / traditional"))
    for k in extra_objects():
        ensure(k)
    for k in SYSTEMS:
        ensure(k)
        recs[k]["desigs"].append(SYSTEMS[k]["desig"])
    for k, d in DEFINED.items():
        ensure(k)
        recs[k]["names"].append((d["name"], "conventional name"))

    # SIMBAD identifiers: names (lowest precedence) and designations.
    for k, rec in recs.items():
        row = sim_hip.get(k) if isinstance(k, int) else sim_obj.get(k)
        sids = ids_by_oid.get(int(row["oid"]), []) if row else []
        if row and row["main_id"] not in sids:
            sids = [row["main_id"]] + sids
        rec["simbad_ids"] = sids
        for i in sids:
            if i.startswith("NAME "):
                rec["names"].append((norm_id(i), "SIMBAD"))

    # V magnitudes (also used to rank shared keys: the brightest object wins).
    for k, rec in recs.items():
        if k in SYSTEMS:
            rec["vmag"], rec["vref"] = SYSTEMS[k]["vmag"], SYSTEMS[k]["vmag_ref"]
            continue
        row = sim_hip.get(k) if isinstance(k, int) else sim_obj.get(k)
        v = vflux.get(int(row["oid"])) if row else None
        if v is None and k in comps:
            v = vflux.get(int(comps[k]["oid"]))
            v = (v[0], f"{v[1]} (component A)") if v else None
        if v is None and isinstance(k, int) and k in hip1_v:
            v = (hip1_v[k], "Hipparcos Catalogue 1997 field H5")
        rec["vmag"], rec["vref"] = v if v else (NO_MAG, "")
    bright = sorted(recs, key=lambda k: (recs[k]["vmag"], record_id(k)))

    # Names: IAU-CSN, then traditional aliases and conventional names, then SIMBAD
    # names (brightest object first). The first registration of a key wins.
    name_keys: dict[str, object] = {}
    spell: dict[str, str] = {}
    for prio in ("IAU-CSN", "Allen 1899 / traditional", "conventional name", "SIMBAD"):
        for k in bright:
            for n, src in recs[k]["names"]:
                if src == prio and fold(n) not in name_keys:
                    name_keys[fold(n)] = k
                    spell[fold(n)] = n

    # Designation in standard form, and every accepted designation spelling.
    desig_count: dict[str, int] = {}
    for rec in recs.values():
        if rec.get("wgsn_desig"):
            desig_count[rec["wgsn_desig"]] = desig_count.get(rec["wgsn_desig"], 0) + 1
    for k, rec in recs.items():
        sids = sorted((i for i in rec["simbad_ids"] if not i.startswith("NAME ")), key=lambda i: (id_rank(i), i))
        wd = rec.get("wgsn_desig", "")
        if wd and desig_count[wd] > 1 and rec.get("comp") in ("A", "B", "C"):
            wd = f"{wd} {rec['comp']}"
        if rec["desigs"]:
            disp = rec["desigs"][0]
        elif wd:
            disp = wd
        elif sids:
            disp = norm_id(sids[0])
        elif isinstance(k, str) and k not in DEFINED:
            disp = norm_id(k)
        else:
            disp = ""
        rec["desig"] = disp
        keys = [disp, rec.get("wgsn_desig", "")] + [norm_id(i) for i in sids]
        if isinstance(k, str) and k not in SYSTEMS and k not in DEFINED:
            keys.append(norm_id(k))
        if isinstance(k, int):
            keys.append(f"HIP {k}")
        rec["desig_keys"] = [v for d in keys if d for v in desig_variants(d)]
    desig_keys: dict[str, object] = {}
    for k in bright:
        for d in recs[k]["desig_keys"]:
            if fold(d) not in desig_keys:
                desig_keys[fold(d)] = k
                spell.setdefault("," + fold(d), d)

    for k, rec in recs.items():
        assert not rec["desig"] or desig_keys.get(fold(rec["desig"])) == k, \
            f"{record_id(k)}: designation {rec['desig']!r} belongs to another record"

    # Display name: the record's first name that resolves to it; returned string
    # "Name,designation" (just the name for designation-less reference points).
    for k, rec in recs.items():
        rec["name"] = next((n for n, _ in rec["names"] if name_keys.get(fold(n)) == k), "")
        rec["full"] = f"{rec['name']},{rec['desig']}" if rec["desig"] else rec["name"]

    # Astrometry.
    for k, rec in recs.items():
        if k in SYSTEMS:
            sy = SYSTEMS[k]
            (c1, c2), (m1, m2) = sy["components"], sy["masses"]
            u1 = unit(float(hip2[c1]["ra_rad"]), float(hip2[c1]["de_rad"]))
            u2 = unit(float(hip2[c2]["ra_rad"]), float(hip2[c2]["de_rad"]))
            w = m2 / (m1 + m2)
            v = [(1 - w) * x1 + w * x2 for x1, x2 in zip(u1, u2)]
            rec["a"] = {"ra_rad": math.atan2(v[1], v[0]) % (2 * math.pi),
                        "de_rad": math.atan2(v[2], math.hypot(v[0], v[1])),
                        "pm_ra": sy["pm_ra"], "pm_de": sy["pm_de"], "plx": sy["plx"]}
            rec["rv"], rec["ref"], rec["hip"] = sy["rv"], sy["rv_ref"], 0
        elif k in DEFINED:
            rec["a"] = defined_astrometry(DEFINED[k])
            rec["rv"], rec["ref"], rec["hip"] = 0.0, DEFINED[k]["ref"], 0
        elif isinstance(k, int) and k in hip2:
            a = dict(hip2[k])
            r = sim_hip.get(k, {})
            rec["rv"], rec["ref"] = (num(r.get("rvz_radvel")) or 0.0), (r.get("rvz_bibcode") or "")
            if k in OVERRIDES:
                o = OVERRIDES[k]
                a.update({x: o[x] for x in ("pm_ra", "pm_de", "plx")})
                rec["rv"], rec["ref"] = o["rv"], o["rv_ref"]
            rec["a"], rec["hip"] = a, k
        else:
            r = sim_hip[k] if isinstance(k, int) else sim_obj[k]
            rec["a"], rec["ref"] = simbad_astrometry(r, comps.get(k))
            rec["rv"], rec["hip"] = rec["a"]["rv"], (k if isinstance(k, int) else 0)

    order = sorted(recs, key=lambda k: (fold(recs[k]["full"]).lstrip(","), record_id(k)))
    index = {k: i for i, k in enumerate(order)}
    assert len(order) < 65536

    lines = [
        "// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial",
        "// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI",
        "// Data: Hipparcos New Reduction (ESA; van Leeuwen 2007), SIMBAD (CDS), IAU WGSN star names (CC BY 4.0), traditional names (R. H. Allen 1899). See NOTICE.",
        "// @generated by tools/build_stars.py — do not edit by hand.",
        "//",
        "// Names: IAU Catalog of Star Names (WGSN, CC BY 4.0), SIMBAD (CDS) identifiers, traditional names",
        "//   (R. H. Allen 1899). Designations in standard form; see docs/STAR_NAMES.md.",
        "// Astrometry: Hipparcos New Reduction (van Leeuwen 2007, CDS I/311),",
        "//   ICRS, epoch J1991.25 TT; non-Hipparcos objects from SIMBAD (CDS), propagated to J1991.25.",
        "// RV and V magnitudes: SIMBAD (CDS). Per-record references in `rv_ref`. See docs/STARS_REPORT.md.",
        "",
        "// Catalogue values, not mathematical constants.",
        "#![allow(clippy::approx_constant)]",
        "",
        "use super::StarRecord;",
        "",
        "/// Catalogue epoch (Hipparcos J1991.25), JD TT.",
        f"pub const CATALOG_EPOCH_JD: f64 = {HIP_EPOCH_JD:.4f};",
        "",
        f"pub static STARS: [StarRecord; {len(order)}] = [",
    ]
    for k in order:
        rec = recs[k]
        s = {"id": record_id(k), "name": rec["name"], "nom": rec["desig"], "desig": rec["desig"], "vmag": rec["vmag"]}
        ref = rec["ref"] + (f"; V {rec['vref']}" if rec["vref"] else "; no V")
        lines.append(record_line(s, rec["hip"], rec["a"], rec["rv"], ref))
    lines += ["];", ""]
    lines.append("/// Default returned string of each record (same order as `STARS`).")
    lines.append(f"pub static FULL: [&str; {len(order)}] = [")
    lines += [f"    {rs_str(recs[k]['full'])}," for k in order]
    lines += ["];", ""]
    for title, table in (("NAMES", name_keys), ("DESIGNATIONS", desig_keys)):
        rows = sorted((key, recs[k]["full"], index[k]) for key, k in table.items())
        lines.append("/// (folded key, returned \"Name,designation\", star index), sorted by key.")
        lines.append(f"pub static {title}: [(&str, &str, u16); {len(rows)}] = [")
        lines += [f"    ({rs_str(key)}, {rs_str(full)}, {i})," for key, full, i in rows]
        lines += ["];", ""]
    # No Gaia-derived value may remain: every bibcode cited by a record is checked.
    bib_re = re.compile(r"(?<![\w.])(\d{4}[^\s;,()]{15})(?![\w.])")
    for k in order:
        for b in bib_re.findall(recs[k]["ref"] + " " + recs[k]["vref"]):
            assert not g(b), f"{record_id(k)}: Gaia-derived value ({b}) in the catalogue"
    used = {record_id(k) for k in order}
    log = ["# Values replaced because their source was Gaia-derived (tools/build_stars.py).",
           "# object\tfield\tGaia-derived source\treplacement (bibcode, or 'dropped')"]
    log += [f"{o}\t{f}\t{old}\t{new}" for o, f, old, new in sorted(set(changes))]
    (ROOT / "tools" / "star_gaia_replacements.tsv").write_text("\n".join(log) + "\n", encoding="utf-8")
    summary = {}
    for _, f, _, new in set(changes):
        summary.setdefault(f, [0, 0])[new == "dropped"] += 1
    print("Gaia-derived values replaced/dropped by field:", summary, file=sys.stderr)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    write_names_doc(recs, order, name_keys, desig_keys, spell)
    n_v = sum(1 for k in recs if recs[k]["vmag"] != NO_MAG)
    n_named = sum(1 for k in recs if recs[k]["name"])
    print(f"wrote {OUT.relative_to(ROOT)}: {len(order)} records ({n_named} named, {n_v} with V), "
          f"{len(name_keys)} names, {len(desig_keys)} designations")


def write_names_doc(recs, order, name_keys, desig_keys, spell) -> None:
    """docs/STAR_NAMES.md: every object with its returned string, names (with source) and designations."""
    by_rec_names: dict = {}
    for key, k in name_keys.items():
        by_rec_names.setdefault(k, []).append(spell[key])
    by_rec_desigs: dict = {}
    for key, k in desig_keys.items():
        by_rec_desigs.setdefault(k, []).append(spell.get("," + key, key))
    src_of = {}
    for k, rec in recs.items():
        for n, src in rec["names"]:
            src_of.setdefault((k, n), src)
    tag = {"IAU-CSN": "IAU", "Allen 1899 / traditional": "trad.", "conventional name": "conv.", "SIMBAD": "SIMBAD"}
    md = ["# Fixed-star names",
          "",
          "Generated by `tools/build_stars.py`. Name sources: **IAU** = IAU Catalog of Star Names (WGSN,",
          "CC BY 4.0); **SIMBAD** = SIMBAD (CDS) `NAME` identifiers; **trad.** = traditional names as given",
          "in R. H. Allen, *Star Names: Their Lore and Meaning* (1899, public domain) and common astrological",
          "usage; **conv.** = conventional names of defined reference points. Designations are SIMBAD/IAU",
          "standard forms (Bayer `alf Vir` / `α Vir`, Flamsteed `61 Cyg`, variable stars, HR, HD, HIP,",
          "Messier/NGC/IC). See `docs/STARS_REPORT.md` for the astrometric sources.",
          "",
          "**Matching.** Case-insensitive, spaces ignored (`\"asellus australis\"` = `\"Asellus Australis\"`).",
          "`\"Spica,anything\"` uses the name part; `\",alf Vir\"` looks up a designation. When a key is",
          "shared, IAU names win over traditional names, which win over SIMBAD names; shared",
          "designations go to the brightest object. Queries starting with a digit are not star numbers.",
          "",
          "**Returned string.** `\"Name,designation\"` as in the first column (just the name for the two",
          "defined reference points); objects without a public name return `\",designation\"`.",
          "",
          "**Optional compatibility table.** Deployments that must reproduce the names and",
          "`\"Name,nomenclature\"` strings of classic astrology software can load a private alias",
          "table (`Engine::with_star_compat`, Python `Engine(..., star_compat=path)`, env",
          "`FRENCHEPHEMERIS_STAR_COMPAT`). It is not part of this repository; with it loaded, lookups",
          "and returned strings follow that table, and every name below remains accepted.",
          "",
          f"## Objects ({len(recs)}; {sum(1 for k in recs if recs[k]['name'])} with a public name)",
          "",
          "| Returned | Names (source) | Designations | HIP | V |",
          "|---|---|---|---:|---:|"]
    for k in order:
        r = recs[k]
        names = ", ".join(f"{n} ({tag[src_of.get((k, n), 'SIMBAD')]})"
                          for n in sorted(by_rec_names.get(k, []), key=str.lower))
        des = ", ".join(sorted(set(by_rec_desigs.get(k, [])), key=str.lower))
        v = "" if r["vmag"] == NO_MAG else f"{r['vmag']:.2f}"
        full = r["full"].replace("|", "\\|")
        md.append(f"| `{full}` | {names} | {des} | {r['hip'] or ''} | {v} |")
    DOC.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("fetch", "all"):
        fetch()
    if cmd in ("build", "all"):
        build()
