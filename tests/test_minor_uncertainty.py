# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
import frenchephemeris.swe_compat as swe
from frenchephemeris import minor_uncertainty as mu


def test_chiron_before_saturn_encounter_is_flagged():
    # The table itself needs no pack; give the distance so any installed pack works.
    early, modern = swe.julday(700, 1, 1, 0), swe.julday(2000, 1, 1, 0)
    assert mu.sigma_arcsec(swe.CHIRON, early, dist_au=12.0) > 100
    assert mu.sigma_km(swe.CHIRON, early) > 1e6
    assert mu.sigma_arcsec(swe.CHIRON, modern) < 0.5
    assert not mu.is_uncertain(swe.CHIRON, modern)


def test_calc_ut_and_other_bodies():
    xx, sigma, uncertain = mu.calc_ut(swe.julday(1990, 5, 17, 14.5), swe.CHIRON)
    assert len(xx) == 6 and 0 < sigma < 1 and uncertain is False
    assert mu.sigma_arcsec(swe.MARS, 2451545.0) is None
    assert mu.is_uncertain(swe.PHOLUS, swe.julday(1901, 6, 1, 0))
