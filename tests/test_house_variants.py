# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
import pytest

import frenchephemeris.swe_compat as swe
from frenchephemeris import house_variants as hv

JD = swe.julday(1990, 5, 17, 14.5)
PARIS = (48.8566, 2.3522)


def close(a, b, tol=1e-9):
    return abs(((a - b) + 180.0) % 360.0 - 180.0) < tol


def test_body_on_first_and_seventh_cusp():
    sun = swe.calc_ut(JD, swe.SUN, 0)[0][0]
    moon = swe.calc_ut(JD, swe.MOON, 0)[0][0]
    node = swe.calc_ut(JD, swe.TRUE_NODE, 0)[0][0]
    assert close(hv.houses(JD, *PARIS, "sunrise")[0][0], sun)
    assert close(hv.houses(JD, *PARIS, "moonrise")[0][0], moon)
    assert close(hv.houses(JD, *PARIS, "north_node")[0][0], node)
    assert close(hv.houses(JD, *PARIS, "south_node")[0][0], node + 180.0)
    assert close(hv.houses(JD, *PARIS, "sunset")[0][6], sun)
    assert close(hv.houses(JD, *PARIS, "moonset")[0][6], moon)
    cusps = hv.houses(JD, *PARIS, "sunrise")[0]
    assert all(close(cusps[k], sun + 30.0 * k) for k in range(12))


def test_vehlow_mc_whole_sign_and_angles():
    cusps, ascmc = hv.houses(JD, *PARIS, "equal_vehlow_mc")
    assert close(cusps[9] + 15.0, ascmc[1])
    assert hv.houses(JD, *PARIS, "whole_sign_leo")[0][0] == 120.0
    assert hv.houses(JD, *PARIS, "whole_sign_aries")[0][11] == 330.0
    cusps, ascmc = hv.houses(JD, *PARIS, "asc_mc")
    assert cusps[0] == ascmc[0] and cusps[9] == ascmc[1] and cusps[1] is None
    _, real = swe.houses_ex(JD, *PARIS, b"P")
    assert ascmc == real


def test_sidereal_and_unknown():
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    sun = swe.calc_ut(JD, swe.SUN, swe.FLG_SIDEREAL)[0][0]
    assert close(hv.houses(JD, *PARIS, "sunrise", swe.FLG_SIDEREAL)[0][0], sun)
    with pytest.raises(ValueError):
        hv.houses(JD, *PARIS, "nonsense")
