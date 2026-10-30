#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Write docs/ASTEROID_LIST.json: the numbered asteroids bundled in
fepk-asteroids-v1.bin, each with the reasons it is included.

Sources of the selection (names are filled in from the JPL SPK at build time):
  swiss_docs   - listed in the Swiss Ephemeris documentation's asteroid
                 introduction (swisseph.htm section on asteroids: families,
                 "transplutonians", love/wisdom/etc. asteroids), read as
                 documentation; the section heading is kept as the reason;
  requested    - named in the FreeAstroAPI asteroid request (astrological use);
  classical    - MPC numbers 1-200: the classical, mostly mythologically named
                 asteroids that asteroid-astrology programs ship by default;
  tno_centaur  - dwarf planets, large trans-Neptunian objects and centaurs;
  neo_myth     - further mythologically named near-Earth/Mars-crossing asteroids
                 often used in astrology.
Excluded: the 8 bodies of the minor pack (Ceres, Pallas, Juno, Vesta, Chiron,
Pholus, Eris, Eros), which keep their own ids and longer coverage.

    python tools/asteroid_list.py [swisseph_htm_path]
"""
from __future__ import annotations

import html
import json
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "ASTEROID_LIST.json")
MINOR_PACK = {1, 2, 3, 4, 2060, 5145, 136199, 433}

REQUESTED = {10: "Hygiea", 5: "Astraea", 16: "Psyche", 1181: "Lilith", 80: "Sappho", 1221: "Amor",
             944: "Hidalgo", 7066: "Nessus", 10199: "Chariklo", 8405: "Asbolus", 90377: "Sedna",
             50000: "Quaoar", 90482: "Orcus", 136108: "Haumea", 136472: "Makemake", 20000: "Varuna",
             28978: "Ixion", 225088: "Gonggong", 120347: "Salacia", 174567: "Varda", 6: "Hebe", 7: "Iris",
             8: "Flora", 9: "Metis", 114: "Kassandra", 30: "Urania", 1930: "Lucifer", 2063: "Bacchus",
             1862: "Apollo", 1566: "Icarus", 1685: "Toro", 55: "Pandora", 43: "Ariadne", 75: "Eurydike",
             1170: "Siva", 149: "Medusa", 399: "Persephone", 100: "Hekate", 128: "Nemesis", 349: "Dembowska"}

TNO_CENTAUR = {38628: "Huya", 65489: "Ceto", 42355: "Typhon", 19521: "Chaos", 58534: "Logos",
               55565: "2002 AW197", 84522: "2002 TC302", 52872: "Okyrhoe", 54598: "Bienor",
               32532: "Thereus", 31824: "Elatus", 10370: "Hylonome", 49036: "Pelion", 52975: "Cyllarus",
               60558: "Echeclus", 55576: "Amycus", 83982: "Crantor", 471143: "Dziewanna", 79360: "Sila-Nunam",
               38083: "Rhadamanthus", 15760: "Albion", 47171: "Lempo", 26375: "1999 DE9", 148780: "Altjira",
               20461: "Dioretsa", 5335: "Damocles", 944: "Hidalgo", 3552: "Don Quixote"}

NEO_MYTH = {2101: "Adonis", 2212: "Hephaistos", 3200: "Phaethon", 1864: "Daedalus", 1915: "Quetzalcoatl",
            3671: "Dionysus", 2340: "Hathor", 5011: "Ptah", 1943: "Anteros", 1036: "Ganymed", 4179: "Toutatis",
            4660: "Nereus", 99942: "Apophis", 1620: "Geographos", 1866: "Sisyphus", 1580: "Betulia",
            5786: "Talos", 3362: "Khufu", 2062: "Aten", 2100: "Ra-Shalom", 3753: "Cruithne", 1627: "Ivar",
            1980: "Tezcatlipoca", 3122: "Florence", 1685: "Toro", 4183: "Cuno", 2329: "Orthos"}


def swiss_docs(path: str) -> dict:
    t = open(path, encoding="latin-1").read()
    t = re.sub(r"<br[^>]*>|</p>|</tr>|</div>", "\n", t)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t)).replace("\xa0", " ").replace("Â", " ")
    i = t.find("Centaurs:")
    block = t[i - 15000:i + 25000]
    out, sec = {}, ""
    for line in (l.strip() for l in block.splitlines() if l.strip()):
        m = re.match(r"^(\d{1,6})\s+([A-Z][^#]*?)\s*(#.*)?$", line)
        if m:
            out.setdefault(int(m.group(1)), (m.group(2).strip(), sec))
        elif line.startswith("#") and ":" in line and len(line) < 80:
            sec = line.strip("# ").strip()
    return out


def main():
    htm = sys.argv[1] if len(sys.argv) > 1 else None
    entries: dict[int, dict] = {}

    def add(n, name, reason):
        if n in MINOR_PACK:
            return
        e = entries.setdefault(n, {"number": n, "name": name, "reasons": []})
        if name and not e["name"]:
            e["name"] = name
        if reason not in e["reasons"]:
            e["reasons"].append(reason)
    if htm:
        for n, (name, sec) in swiss_docs(htm).items():
            add(n, name, f"swiss_docs: {sec.rstrip(':')[:60] or 'main list'}")
    for n, name in REQUESTED.items():
        add(n, name, "requested")
    for n in range(1, 201):
        add(n, "", "classical")
    for n, name in TNO_CENTAUR.items():
        add(n, name, "tno_centaur")
    for n, name in NEO_MYTH.items():
        add(n, name, "neo_myth")
    rows = sorted(entries.values(), key=lambda e: e["number"])
    json.dump({"generator": "tools/asteroid_list.py", "excluded_minor_pack": sorted(MINOR_PACK),
               "count": len(rows), "bodies": rows}, open(OUT, "w"), indent=0)
    print(len(rows), "bodies ->", OUT)


if __name__ == "__main__":
    main()
