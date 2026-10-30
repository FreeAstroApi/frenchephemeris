# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""swe_compat settings must be scoped per thread/task, never process-global."""

import threading

import frenchephemeris.swe_compat as swe


def test_sid_mode_is_thread_scoped():
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    seen = {}

    def worker():
        seen["before"] = swe._sid_mode.get()[0]
        swe.set_sid_mode(swe.SIDM_RAMAN)
        seen["after"] = swe._sid_mode.get()[0]

    t = threading.Thread(target=worker)
    t.start()
    t.join()
    assert seen == {"before": swe.SIDM_FAGAN_BRADLEY, "after": swe.SIDM_RAMAN}
    assert swe._sid_mode.get()[0] == swe.SIDM_LAHIRI


def test_topo_is_thread_scoped():
    swe.set_topo(2.35, 48.85, 35)
    out = []
    t = threading.Thread(target=lambda: out.append(swe._topo.get()))
    t.start()
    t.join()
    assert out == [None]


def test_lunar_points_honour_equatorial_flag():
    """Nodes and apogees are ecliptic points: their RA/Dec must be the ecliptic
    position rotated by the true obliquity, like any other body."""
    jd = 2448028.98
    eps = swe.calc(jd, swe.ECL_NUT, 0)[0][0]
    for body in (swe.MEAN_NODE, swe.TRUE_NODE, swe.MEAN_APOG, swe.OSCU_APOG, swe.INTP_APOG, swe.INTP_PERG):
        lon, lat = swe.calc(jd, body, swe.FLG_SWIEPH)[0][:2]
        ra, dec = swe.calc(jd, body, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL)[0][:2]
        exp_ra, exp_dec, _ = swe.cotrans((lon, lat, 1.0), -eps)
        assert abs((ra - exp_ra + 180) % 360 - 180) < 1e-6, body
        assert abs(dec - exp_dec) < 1e-6, body
        assert abs(dec) > 1.0, body  # not left on the ecliptic


def test_true_star_ayanamsa_follows_the_body_corrections():
    """True Citra is Spica at 180 deg. A body computed without aberration and light
    deflection (TRUEPOS) is measured against Spica computed the same way, so its
    sidereal longitude stays star-relative."""
    jd = 2440667.2
    swe.set_sid_mode(swe.SIDM_TRUE_CITRA)
    for flags in (swe.FLG_TRUEPOS, swe.FLG_NOABERR | swe.FLG_NOGDEFL):
        tropical = swe.calc_ut(jd, swe.MOON, flags)[0][0]
        sidereal = swe.calc_ut(jd, swe.MOON, flags | swe.FLG_SIDEREAL)[0][0]
        spica = swe.fixstar2_ut("Spica", jd, swe.FLG_NOABERR | swe.FLG_NOGDEFL)[0][0]
        expected = (tropical - (spica - 180.0)) % 360.0
        assert abs((sidereal - expected + 180.0) % 360.0 - 180.0) * 3600 < 0.01, flags
    swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)


def test_sidereal_and_tropical_paths_through_every_wrapper():
    """calc_ut, calc, houses_ex, get_ayanamsa*, fixstar2_ut with and without FLG_SIDEREAL,
    including option bits and a user-defined mode."""
    jd = 2451545.0
    for mode in (swe.SIDM_LAHIRI, swe.SIDM_TRUE_CITRA, swe.SIDM_LAHIRI | swe.SIDBIT_ECL_DATE):
        swe.set_sid_mode(mode)
        trop = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)[0][0]
        sid = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SIDEREAL)[0][0]
        assert abs(sid - swe.calc(jd + swe.deltat(jd), swe.SUN, swe.FLG_SWIEPH | swe.FLG_SIDEREAL)[0][0]) < 1e-9
        assert 20 < (trop - sid) % 360 < 26, mode
        cusps = swe.houses_ex(jd, 48.85, 2.35, b"P", swe.FLG_SIDEREAL)[0]
        assert abs(((swe.houses_ex(jd, 48.85, 2.35, b"P")[0][0] - cusps[0]) % 360) - swe.get_ayanamsa_ut(jd)) < 0.01
        assert swe.get_ayanamsa_ex_ut(jd)[1] > 0 and swe.get_ayanamsa(jd) > 0
        assert swe.fixstar2_ut("Spica", jd, swe.FLG_SWIEPH | swe.FLG_SIDEREAL)[0][0] < 360
    swe.set_sid_mode(swe.SIDM_USER, 2451545.0, 23.5)
    assert abs(swe.get_ayanamsa(2451545.0) - 23.5) < 1e-9
    assert swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SIDEREAL)[0][0] < 360
    swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)
