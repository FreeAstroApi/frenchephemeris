# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""swe_compat quick wins (gap analysis step 1) against the pyswisseph oracle, when installed.

Delta T settings must be scoped per thread/task like set_sid_mode."""

import os
import threading

import pytest

import frenchephemeris.swe_compat as swe

try:
    import swisseph as sw

    sw.set_ephe_path(os.environ.get("SE_EPHE_PATH", ""))
except ImportError:  # pragma: no cover
    sw = None

needs_oracle = pytest.mark.skipif(sw is None or not os.environ.get("SE_EPHE_PATH"),
                                  reason="pyswisseph and its data files (SE_EPHE_PATH) needed for the oracle")
JD = 2460000.3


def test_delta_t_settings_are_context_scoped():
    base = swe.deltat(JD)
    swe.set_delta_t_userdef(100 / 86400)
    seen = {}

    def worker():
        seen["dt"] = swe.deltat(JD)
        seen["tid"] = swe.get_tid_acc()
        swe.set_tid_acc(swe.TIDAL_DE200)
        seen["tid_after"] = swe.get_tid_acc()

    t = threading.Thread(target=worker)
    t.start()
    t.join()
    assert seen["dt"] == base
    assert seen["tid"] == -25.936 and seen["tid_after"] == swe.TIDAL_DE200
    assert swe.deltat(JD) == pytest.approx(100 / 86400, abs=1e-15)
    assert swe.get_tid_acc() == -25.936
    moon_ut = swe.calc_ut(JD, swe.MOON)[0][0]
    moon_tt = swe.calc(JD + 100 / 86400, swe.MOON)[0][0]
    assert moon_ut == pytest.approx(moon_tt, abs=1e-9)
    swe.set_delta_t_userdef(swe.DELTAT_AUTOMATIC)
    assert swe.deltat(JD) == base


@needs_oracle
def test_ecl_nut_and_output_flags():
    ours = swe.calc(JD, swe.ECL_NUT, 0)[0]
    ref = sw.calc(JD, sw.ECL_NUT, 0)[0]
    assert all(abs(a - b) < 1e-5 for a, b in zip(ours, ref))
    for fl in (swe.FLG_XYZ, swe.FLG_RADIANS, swe.FLG_J2000 | swe.FLG_EQUATORIAL | swe.FLG_ICRS, swe.FLG_ASTROMETRIC):
        a = swe.calc(JD, swe.MARS, fl | swe.FLG_SPEED)[0]
        b = sw.calc(JD, sw.MARS, fl | sw.FLG_SPEED)[0]
        assert all(abs(x - y) < 1e-6 for x, y in zip(a[:3], b[:3])), fl


@needs_oracle
def test_calc_pctr_crossings_and_local_lunar_eclipse():
    a = swe.calc_pctr(JD, swe.MARS, swe.JUPITER)[0]
    b = sw.calc_pctr(JD, sw.MARS, sw.JUPITER, sw.FLG_SPEED)[0]
    assert abs(a[0] - b[0]) < 1e-4
    t = swe.mooncross_node(JD)
    assert abs(t[0] - sw.mooncross_node(JD)[0]) * 86400 < 1
    h = swe.helio_cross(swe.MARS, 100.0, JD)
    assert abs(h - sw.helio_cross(sw.MARS, 100.0, JD)) * 86400 < 60
    geo = (2.35, 48.85, 35.0)
    f, tret, attr = swe.lun_eclipse_when_loc(JD, geo)
    g, sret, sattr = sw.lun_eclipse_when_loc(JD, geo)
    assert f == g
    assert all((x == 0) == (y == 0) and abs(x - y) * 86400 < 60 for x, y in zip(tret, sret))


@needs_oracle
def test_utc_and_calendar():
    assert swe.utc_to_jd(2016, 12, 31, 23, 59, 60.5)[0] == pytest.approx(sw.utc_to_jd(2016, 12, 31, 23, 59, 60.5)[0], abs=1e-9)
    with pytest.raises(swe.Error):
        swe.utc_to_jd(2015, 6, 29, 23, 59, 60.0)
    y, m, d, h, mi, s = swe.jdet_to_utc(swe.utc_to_jd(2016, 12, 31, 23, 59, 60.5)[0])
    assert (y, m, d, h, mi) == (2016, 12, 31, 23, 59) and abs(s - 60.5) < 1e-4
    assert swe.utc_time_zone(2023, 1, 1, 0, 30, 15.5, 5.5)[:5] == sw.utc_time_zone(2023, 1, 1, 0, 30, 15.5, 5.5)[:5]
    assert swe.date_conversion(2023, 2, 30, 12.0, b"g") == sw.date_conversion(2023, 2, 30, 12.0, b"g")
    assert swe.day_of_week(2451545.0) == sw.day_of_week(2451545.0)
    assert abs(swe.sidtime0(JD, 23.44, 0.003) - sw.sidtime0(JD, 23.44, 0.003)) * 3600 < 0.01


@needs_oracle
def test_helpers_and_names():
    for v in (12.3456789, 29.9999998, 359.9999999, -12.3456789):
        for fl in (0, 1, 2, 4, 8, 9, 25, 1024, 1025, 1041):
            a, b = swe.split_deg(v, fl), sw.split_deg(v, fl)
            assert a[:3] == b[:3] and a[4] == b[4] and abs(a[3] - b[3]) < 1e-6, (v, fl)
    assert swe.cs2lonlatstr(123456789, b"E", b"W") == sw.cs2lonlatstr(123456789, b"E", b"W")
    assert swe.cs2timestr(123456789, b":", False) == sw.cs2timestr(123456789, b":", False)
    assert swe.cs2degstr(123456789) == "12°56'07"
    for x in (-1, 129599999, 10799950):
        assert swe.csnorm(x) == sw.csnorm(x) and swe.csroundsec(x) == sw.csroundsec(x)
    assert swe.difdeg2n(10, 350) == sw.difdeg2n(10, 350) and swe.degnorm(-30) == sw.degnorm(-30)
    assert [swe.get_planet_name(i) for i in range(23)] == [sw.get_planet_name(i) for i in range(23)]
    assert swe.get_ayanamsa_name(swe.SIDM_LAHIRI) == sw.get_ayanamsa_name(sw.SIDM_LAHIRI)
    path, start, end, denum = swe.get_current_file_data(0)
    assert path.endswith(".bin") and start < 2451545.0 < end and denum == 441
    assert swe.get_library_path().endswith((".so", ".dylib", ".pyd"))
    assert swe.version.startswith("FrenchEphemeris")
