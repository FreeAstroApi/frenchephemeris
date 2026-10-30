# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Read the generated Rust catalogue so fixtures use exactly the shipped numbers."""
import re
from pathlib import Path

CATALOG = Path(__file__).resolve().parents[3] / "src" / "stars" / "catalog.rs"
_REC = re.compile(r'StarRecord \{ id: "[^"]*", name: "(?P<name>[^"]*)", nomenclature: "(?P<nom>[^"]*)", designation: "[^"]*", '
                  r'hip: (?P<hip>\d+), ra_rad: (?P<ra>[-\d.e]+), dec_rad: (?P<dec>[-\d.e]+), '
                  r'pmra_mas_yr: (?P<pmra>[-\d.e]+), pmdec_mas_yr: (?P<pmde>[-\d.e]+), plx_mas: (?P<plx>[-\d.e]+), '
                  r'rv_km_s: (?P<rv>[-\d.e]+),')
_NAME = re.compile(r'^    \("([^"]*)", "([^"]*)", (\d+)\),$', re.M)
_EPOCH = re.compile(r"CATALOG_EPOCH_JD: f64 = ([\d.]+);")


def load():
    text = CATALOG.read_text(encoding="utf-8")
    stars = []
    for m in _REC.finditer(text):
        d = m.groupdict()
        stars.append({"name": d["name"], "nom": d["nom"], "hip": int(d["hip"]),
                      **{k: float(d[k]) for k in ("ra", "dec", "pmra", "pmde", "plx", "rv")}})
    sec = text[text.index("pub static NAMES"):text.index("pub static DESIGNATIONS")]
    names = [(key, full, int(i)) for key, full, i in _NAME.findall(sec)]
    return float(_EPOCH.search(text).group(1)), stars, names
