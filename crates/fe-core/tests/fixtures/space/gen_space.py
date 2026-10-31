# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Generate space.json: JPL Horizons places for observers on other worlds.

For each site (planetodetic latitude, east longitude on the IAU 2015 ellipsoid) and TT date,
Horizons returns the astrometric ICRF RA/Dec (quantity 1) and the apparent airless azimuth and
elevation (quantity 4) of a few targets. Horizons takes west longitudes on IAU west-positive
bodies, so east longitudes are sent negated there; the header line it echoes is stored to show
how the site was read.

Usage: python gen_space.py  (needs `curl`; writes space.json next to this file)
"""
import json
import subprocess
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent

# world, NAIF centre, east lon, lat, height km
SITES = [
    ("mars", 499, 77.45, 18.44, 0.0),     # Jezero crater
    ("mars", 499, 200.0, -60.0, 1.0),
    ("moon", 301, 23.47, 0.67, 0.0),      # Tranquility Base
    ("mercury", 199, 10.0, -30.0, 0.0),
    ("venus", 299, 30.0, 40.0, 0.0),
    ("jupiter", 599, 100.0, 20.0, 0.0),
    ("saturn", 699, 250.0, -10.0, 0.0),
    ("uranus", 799, 40.0, -30.0, 0.0),
    ("neptune", 899, 300.0, 50.0, 0.0),
    ("pluto", 999, 180.0, 25.0, 0.0),
    ("earth", 399, 2.35, 48.85, 0.0),     # Paris, control
]
JDS_TT = [2433282.5, 2448029.104167, 2460000.5, 2473459.75]
# Pack ids: outer planets as system barycentres, like the engine.
TARGETS = {"sun": 10, "earth": 399, "moon": 301, "venus": 299, "mars": 4, "jupiter": 5, "saturn": 6}


def horizons(target, centre, lon_field, lat, h):
    params = {
        "format": "text", "COMMAND": f"'{target}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'", "CENTER": f"'coord@{centre}'", "COORD_TYPE": "'GEODETIC'",
        "SITE_COORD": f"'{lon_field},{lat},{h}'", "TLIST": "'" + " ".join(str(j) for j in JDS_TT) + "'",
        "TLIST_TYPE": "'JD'", "TIME_TYPE": "'TT'", "QUANTITIES": "'1,4'", "ANG_FORMAT": "'DEG'",
        "EXTRA_PREC": "'YES'", "CSV_FORMAT": "'YES'", "CAL_FORMAT": "'JD'", "APPARENT": "'AIRLESS'",
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)
    return subprocess.run(["curl", "-s", "--max-time", "120", url], capture_output=True, text=True, check=True).stdout


def fetch(target, centre, lon, lat, h):
    text = horizons(target, centre, lon, lat, h)
    if "Input east-longitude as negative" in text:
        text = horizons(target, centre, -lon, lat, h)
    if "$$SOE" not in text:
        raise RuntimeError(text[-800:])
    site_line = next(l for l in text.splitlines() if l.startswith("Center geodetic"))
    rows = []
    for row in text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        c = [x.strip() for x in row.split(",")]
        nums = [float(x) for x in c[1:] if x and x.lstrip("-").replace(".", "", 1).isdigit()]
        rows.append({"jd_tt": float(c[0]), "ra": nums[0], "dec": nums[1], "az": nums[2], "el": nums[3]})
    return site_line.split(":", 1)[1].split("{")[0].strip(), rows


def main():
    out = []
    for world, centre, lon, lat, h in SITES:
        for name, tid in TARGETS.items():
            if name == world:
                continue
            site_line, rows = fetch(tid, centre, lon, lat, h)
            out.append({"world": world, "lon": lon, "lat": lat, "height_km": h, "target": name,
                        "horizons_site": site_line, "rows": rows})
            print(world, name, site_line)
    (HERE / "space.json").write_text(json.dumps({
        "source": "JPL Horizons API, astrometric ICRF RA/Dec (q1) and apparent airless Az/El (q4), TT",
        "cases": out}, indent=1) + "\n")


if __name__ == "__main__":
    main()
