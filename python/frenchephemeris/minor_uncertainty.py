# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""How well a minor body's position is known, from JPL's orbit covariance.

    import frenchephemeris.swe_compat as swe
    from frenchephemeris import minor_uncertainty as mu
    mu.sigma_arcsec(swe.CHIRON, jd_ut)      # 1-sigma error of the geocentric direction, arcsec
    mu.is_uncertain(swe.CHIRON, jd_ut)      # True when 3 sigma exceeds the 0.001 deg display step
    xx, sigma, uncertain = mu.calc_ut(jd_ut, swe.CHIRON)

Covers Chiron, Pholus, Eris and Eros, whose positions before JPL's tabulation come from our own
integration: the six covariance clones of JPL's orbit solution are propagated with the same
model (tools/minor_uncertainty.py), and the largest 1-sigma axis of the position covariance is
turned into an angle at the body's distance from Earth. Chiron is known to a fraction of an
arcsecond after 1000 but only to about a degree before its close Saturn encounter of 720;
other bodies return None (their errors are far below the display step).
"""

from __future__ import annotations

import functools
import json
import math
import os

from . import swe_compat as swe

AU_KM = 149597870.7
ARCSEC_PER_RAD = 206264.80624709636
DISPLAY_STEP_ARCSEC = 3.6  # 0.001 deg
_IDS = {swe.CHIRON: 2060, swe.PHOLUS: 5145, swe.AST_OFFSET + 136199: 136199, swe.AST_OFFSET + 433: 433}


@functools.lru_cache(maxsize=1)
def _tables():
    with open(os.path.join(os.path.dirname(__file__), "minor_uncertainty.json")) as f:
        data = json.load(f)["bodies"]
    return {int(k): (v["jd_tdb"], [math.log(s) for s in v["sigma_km"]]) for k, v in data.items()}


def sigma_km(body: int, jd: float) -> float | None:
    """Largest 1-sigma axis of the position error, km (log-interpolated; clamped at the ends)."""
    number = _IDS.get(body)
    if number is None:
        return None
    jds, logs = _tables()[number]
    if jd <= jds[0]:
        return math.exp(logs[0])
    if jd >= jds[-1]:
        return math.exp(logs[-1])
    lo, hi = 0, len(jds) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if jds[mid] <= jd else (lo, mid)
    w = (jd - jds[lo]) / (jds[hi] - jds[lo])
    return math.exp(logs[lo] + w * (logs[hi] - logs[lo]))


def sigma_arcsec(body: int, jd_ut: float, dist_au: float | None = None) -> float | None:
    """1-sigma error of the geocentric direction, arcsec; None for bodies without a table."""
    s = sigma_km(body, jd_ut)
    if s is None:
        return None
    if dist_au is None:
        dist_au = swe.calc_ut(jd_ut, body, 0)[0][2]
    return min(s / (dist_au * AU_KM), math.pi) * ARCSEC_PER_RAD


def is_uncertain(body: int, jd_ut: float, n_sigma: float = 3.0,
                 limit_arcsec: float = DISPLAY_STEP_ARCSEC) -> bool:
    """True when n_sigma x the 1-sigma error exceeds limit_arcsec (default: the 0.001 deg step)."""
    s = sigma_arcsec(body, jd_ut)
    return s is not None and n_sigma * s > limit_arcsec


def calc_ut(jd_ut: float, body: int, flags: int = swe.FLG_SWIEPH | swe.FLG_SPEED):
    """(xx, sigma_arcsec, uncertain): swe_compat.calc_ut's position plus its uncertainty."""
    xx, _ = swe.calc_ut(jd_ut, body, flags)
    dist = xx[2] if not flags & (swe.FLG_HELCTR | swe.FLG_BARYCTR) else None
    s = sigma_arcsec(body, jd_ut, dist)
    return xx, s, s is not None and 3.0 * s > DISPLAY_STEP_ARCSEC
