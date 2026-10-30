# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""House variants built on a point other than the Ascendant, as offered by chart software
beyond the classic house systems:

    from frenchephemeris import house_variants as hv
    cusps, ascmc = hv.houses(jd_ut, lat, lon, "sunrise")       # Sun on the 1st cusp, equal houses
    cusps, ascmc = hv.houses(jd_ut, lat, lon, "whole_sign_leo")  # 1st house = 0° Leo

Variants (cusps are ecliptic longitudes in degrees, cusps[0] = 1st house):

- `sunrise`, `moonrise`, `north_node`, `south_node`: equal houses, the body on the 1st cusp
  (the "solar chart" and its lunar and nodal counterparts);
- `sunset`, `moonset`: equal houses, the body on the 7th cusp;
- `equal_vehlow_mc`: equal houses with the MC in the middle of the 10th house (the MC-based
  counterpart of Vehlow's equal houses, which put the Ascendant in the middle of the 1st);
- `whole_sign_<sign>` (`whole_sign_aries` … `whole_sign_pisces`): whole-sign houses whose 1st
  house is that sign, whatever the Ascendant;
- `asc_mc`: no intermediate cusps; only the angles (cusps 1, 4, 7, 10 = ASC, IC, DSC, MC; the
  other cusps are None).

`ascmc` is always the chart's real angles (houses_ex with Placidus is used only for them), so
the Ascendant and MC can still be drawn. `flags` takes swe_compat.FLG_SIDEREAL for sidereal
cusps (with the current thread's sidereal mode); `node` picks TRUE_NODE or MEAN_NODE.
"""

from __future__ import annotations

from . import swe_compat as swe

SIGNS = ("aries", "taurus", "gemini", "cancer", "leo", "virgo",
         "libra", "scorpio", "sagittarius", "capricorn", "aquarius", "pisces")
_BODY_ON_FIRST = {"sunrise": swe.SUN, "moonrise": swe.MOON, "north_node": None, "south_node": None}
_BODY_ON_SEVENTH = {"sunset": swe.SUN, "moonset": swe.MOON}
VARIANTS = (tuple(_BODY_ON_FIRST) + tuple(_BODY_ON_SEVENTH) + ("equal_vehlow_mc", "asc_mc")
            + tuple(f"whole_sign_{s}" for s in SIGNS))


def _equal(first: float) -> tuple[float, ...]:
    return tuple((first + 30.0 * k) % 360.0 for k in range(12))


def houses(jd_ut: float, lat: float, lon: float, variant: str, flags: int = 0,
           node: int = swe.TRUE_NODE):
    """(cusps[12], ascmc[8]) for one of VARIANTS."""
    if variant not in VARIANTS:
        raise ValueError(f"unknown house variant {variant!r}; expected one of {', '.join(VARIANTS)}")
    side = flags & swe.FLG_SIDEREAL
    _cusps, ascmc = swe.houses_ex(jd_ut, lat, lon, b"P", side)
    asc, mc = ascmc[0], ascmc[1]

    def body_lon(body: int) -> float:
        return swe.calc_ut(jd_ut, body, side)[0][0]

    if variant in _BODY_ON_FIRST:
        if variant in ("north_node", "south_node"):
            first = body_lon(node) + (180.0 if variant == "south_node" else 0.0)
        else:
            first = body_lon(_BODY_ON_FIRST[variant])
        return _equal(first), ascmc
    if variant in _BODY_ON_SEVENTH:
        return _equal(body_lon(_BODY_ON_SEVENTH[variant]) + 180.0), ascmc
    if variant == "equal_vehlow_mc":
        # MC at 15° of the 10th house: 10th cusp = MC − 15°, 1st cusp = 10th cusp + 90°
        return _equal(mc - 15.0 + 90.0), ascmc
    if variant == "asc_mc":
        cusps = [None] * 12
        cusps[0], cusps[3], cusps[6], cusps[9] = asc, (mc + 180.0) % 360.0, (asc + 180.0) % 360.0, mc
        return tuple(cusps), ascmc
    sign = SIGNS.index(variant.removeprefix("whole_sign_"))
    return _equal(30.0 * sign), ascmc
