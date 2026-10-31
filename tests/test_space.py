# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""frenchephemeris.ext.space: charts for observers away from Earth.

The Rust tests (crates/fe-core/tests/space.rs) compare positions and horizons against JPL
Horizons; these check the Python API and published physical values."""

import pytest

import frenchephemeris
from frenchephemeris.ext import space

JD = 2460000.5


@pytest.fixture(scope="module")
def eng():
    return frenchephemeris.open_default()


# Axial tilt against the world's own orbit (published values, degrees).
@pytest.mark.parametrize("world, tilt", [
    ("mercury", 0.03), ("venus", 177.36), ("earth", 23.44), ("moon", 1.54), ("mars", 25.19),
    ("jupiter", 3.13), ("saturn", 26.73), ("uranus", 97.77), ("neptune", 28.32), ("pluto", 119.61),
])
def test_native_obliquity_is_the_axial_tilt(eng, world, tilt):
    h = space.houses(eng, JD, space.Site(world, 10.0, 0.0), "R", "native")
    assert h["obliquity"] == pytest.approx(tilt, abs=0.1)


def test_chart_on_mars_has_earth_and_no_mars(eng):
    c = space.chart(eng, JD, space.Site("mars", 77.45, 18.44))
    assert "mars" not in c["bodies"] and "earth" in c["bodies"]
    assert len(c["houses"]["cusps"]) == 12
    # Earth and Moon are 0.1° apart or less seen from Mars.
    d = abs((c["bodies"]["earth"]["lon"] - c["bodies"]["moon"]["lon"] + 180) % 360 - 180)
    assert d < 0.5


def test_earth_site_matches_ordinary_houses(eng):
    import frenchephemeris.swe_compat as swe

    cusps, ascmc = swe.houses(JD, 48.85, 2.35, b"P")
    h = space.houses(eng, JD, space.Site("earth", 2.35, 48.85), "P")
    assert h["asc"] == pytest.approx(ascmc[0], abs=0.02)
    assert h["mc"] == pytest.approx(ascmc[1], abs=0.02)


def test_sun_is_up_at_mars_local_noon(eng):
    # Local mean noon at Jezero on some sol: the Sun is near the meridian, above the horizon.
    lon = 77.45
    jd_tt = space.jd_from_mars_time(53021, 12.0, lon)
    jd_ut = jd_tt - frenchephemeris.deltat(jd_tt)
    az, el = space.horizon(eng, jd_ut, "sun", space.Site("mars", lon, 18.44))
    assert el > 30.0
    assert space.mars_time(jd_tt, lon)["local_mean_solar_time"] == pytest.approx(12.0, abs=1e-6)


def test_point_in_space_and_native_needs_a_world(eng):
    far = space.Point("sun", z=1e9)  # 1e9 km from the Sun towards the ICRF north pole
    p = space.position(eng, JD, "earth", far)
    au = 1.495978707e8
    # Earth stays within ~1.02 AU of the Sun, so it is 1e9 km ± 1.02 AU away.
    assert (1e9 / au - 1.02) < p["distance_au"] < (1e9 / au + 1.02)
    with pytest.raises(ValueError):
        space.position(eng, JD, "earth", far, zodiac="native")
    with pytest.raises(ValueError):
        space.position(eng, JD, "mars", space.Center("mars"))


def test_placidus_is_polar_on_uranus(eng):
    with pytest.raises(frenchephemeris.PolarHousesError):
        space.houses(eng, JD, space.Site("uranus", 0.0, 30.0), "P")
    assert len(space.houses(eng, JD, space.Site("uranus", 0.0, 30.0), "W")["cusps"]) == 12


@pytest.mark.no_packs
def test_mars_time_round_trip():
    jd = space.jd_from_mars_time(53000, 6.5, 137.4)
    t = space.mars_time(jd, 137.4)
    assert t["local_mean_solar_time"] == pytest.approx(6.5, abs=1e-6)
    assert space.WORLDS[0] == "mercury" and "pluto" in space.WORLDS
