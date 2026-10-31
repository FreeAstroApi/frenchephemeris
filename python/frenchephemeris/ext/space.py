# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Charts for an observer away from Earth.

Cast a chart as if someone were born on Mars, on the Moon, on any planet out to Pluto, at the
centre of a world, or at any point in space::

    import frenchephemeris
    from frenchephemeris.ext import space

    eng = frenchephemeris.open_default()
    jezero = space.Site("mars", lon=77.45, lat=18.44)           # east longitude, latitude
    chart = space.chart(eng, 2460000.5, jezero, house_system="P")
    chart["bodies"]["earth"]["lon"], chart["houses"]["asc"]

Times are Julian days in UT (Earth clock time), like the core ``calc_ut``. On Mars,
:func:`mars_time` and :func:`jd_from_mars_time` convert to and from Mars Sol Date and local
mean solar time.

Zodiac:

- ``"earth"`` (default): longitudes in Earth's tropical zodiac (true ecliptic and equinox of
  date), so signs mean what they mean in an ordinary chart.
- ``"native"``: the world's own zodiac: its orbital plane, with 0° at its own vernal
  equinox (the ascending node of its equator on that plane).

Houses use the world's equator and rotation (IAU 2015 rotation models; Earth uses the
engine's own Earth rotation), so the Ascendant is the rising point of the zodiac on that world,
also on Venus, Uranus and Pluto, which spin backwards. Sites use east longitudes and
planetodetic latitudes on the IAU reference ellipsoid. Inside a world's polar circle, Placidus
and Koch are undefined, as on Earth (``PolarHousesError``); Uranus's polar circle reaches 8°
from its equator. An observer at a world's centre or at a point in space has no horizon, so
its chart has positions only.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Union

from .. import _native

__all__ = [
    "WORLDS", "BODIES", "Site", "Center", "Point", "position", "houses", "horizon", "chart",
    "orientation", "radii", "mars_sol_date", "mars_time", "jd_from_mars_time",
]

#: Worlds an observer can stand on.
WORLDS: tuple = tuple(_native.space_worlds())

#: Bodies placed in a chart, with their classic ids.
BODIES: dict = {
    "sun": 0, "moon": 1, "mercury": 2, "venus": 3, "earth": 14, "mars": 4,
    "jupiter": 5, "saturn": 6, "uranus": 7, "neptune": 8, "pluto": 9,
}

# NAIF ids for the centre of a Point observer (outer planets: system barycentres).
_CENTERS = {
    "ssb": 0, "barycenter": 0, "sun": 10, "mercury": 1, "venus": 2, "earth": 399, "moon": 301,
    "mars": 4, "jupiter": 5, "saturn": 6, "uranus": 7, "neptune": 8, "pluto": 9,
}


@dataclass(frozen=True)
class Site:
    """A place on a world's surface: east longitude and planetodetic latitude in degrees,
    height in km above the reference ellipsoid."""

    world: str
    lon: float
    lat: float
    height_km: float = 0.0


@dataclass(frozen=True)
class Center:
    """The centre of a world (no horizon, so no houses)."""

    world: str


@dataclass(frozen=True)
class Point:
    """Any point in space: an offset in km (ICRF axes) from ``center`` (a body name, ``"ssb"``
    for the solar system barycentre, or a NAIF id), moving with that body plus ``velocity``
    (km/day). No horizon, so no houses."""

    center: Union[str, int] = "sun"
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    velocity: Sequence[float] = (0.0, 0.0, 0.0)


Observer = Union[Site, Center, Point]


def _body_id(body: Union[str, int]) -> int:
    if isinstance(body, int):
        return body
    try:
        return BODIES[body.lower()]
    except KeyError:
        raise ValueError(f"unknown body {body!r}; use one of {sorted(BODIES)} or a body id") from None


def _observer_args(obs: Observer) -> dict:
    if isinstance(obs, Site):
        return {"kind": "surface", "world": obs.world, "lon": obs.lon, "lat": obs.lat, "height_km": obs.height_km}
    if isinstance(obs, Center):
        return {"kind": "center", "world": obs.world}
    if isinstance(obs, Point):
        c = obs.center
        if isinstance(c, str):
            try:
                c = _CENTERS[c.lower()]
            except KeyError:
                raise ValueError(f"unknown centre {obs.center!r}; use one of {sorted(_CENTERS)} or a NAIF id") from None
        return {"kind": "point", "center": int(c), "offset": [obs.x, obs.y, obs.z], "velocity": list(obs.velocity)}
    raise TypeError("observer must be a Site, Center or Point")


def _world_of(obs: Observer) -> Optional[str]:
    return obs.world.lower() if isinstance(obs, (Site, Center)) else None


def position(engine, jd_ut: float, body: Union[str, int], observer: Observer,
             zodiac: str = "earth", apparent: bool = True) -> dict:
    """Place of ``body`` seen by ``observer``: longitude, latitude (degrees), distance (AU)
    and their daily rates, plus ``retrograde``. ``apparent=False`` gives the astrometric
    place (no aberration, no light deflection)."""
    lon, lat, dist, dlon, dlat, ddist = _native.space_calc_ut(
        engine, jd_ut, _body_id(body), zodiac=zodiac, apparent=apparent, **_observer_args(observer))
    return {"lon": lon, "lat": lat, "distance_au": dist, "lon_speed": dlon, "lat_speed": dlat,
            "distance_speed": ddist, "retrograde": dlon < 0}


def houses(engine, jd_ut: float, site: Site, house_system: str = "P", zodiac: str = "earth") -> dict:
    """House cusps and angles for a site on a world. ``sidereal_angle`` is the local
    meridian's angle along the world's equator (the ARMC analogue), ``obliquity`` the angle
    between the spin axis and the zodiac pole (above 90° when the world spins against it)."""
    cusps, ascmc, sidereal_angle, spin_lat, obliquity = _native.space_houses_ut(
        engine, jd_ut, site.world, site.lon, site.lat, site.height_km, house_system, zodiac)
    return {
        "system": house_system, "cusps": list(cusps), "asc": ascmc[0], "mc": ascmc[1],
        "vertex": ascmc[3], "equatorial_asc": ascmc[4], "sidereal_angle": sidereal_angle,
        "spin_latitude": spin_lat, "obliquity": obliquity,
    }


def horizon(engine, jd_ut: float, body: Union[str, int], site: Site) -> tuple:
    """Azimuth (from the world's north towards east) and elevation of ``body`` at ``site``,
    degrees, apparent and without refraction."""
    return _native.space_azalt_ut(engine, jd_ut, _body_id(body), site.world, site.lon, site.lat, site.height_km)


def chart(engine, jd_ut: float, observer: Observer, house_system: Optional[str] = "P",
          zodiac: str = "earth", bodies: Optional[Sequence[str]] = None) -> dict:
    """A whole chart: every body (except the observer's own world) and, for a surface site
    with ``house_system`` set, houses and angles."""
    home = _world_of(observer)
    names = [b for b in (bodies or BODIES) if b.lower() != home]
    out = {
        "observer": observer, "zodiac": zodiac, "jd_ut": jd_ut,
        "bodies": {b: position(engine, jd_ut, b, observer, zodiac) for b in names},
        "houses": None,
    }
    if isinstance(observer, Site) and house_system:
        out["houses"] = houses(engine, jd_ut, observer, house_system, zodiac)
    return out


def orientation(world: str, jd_tdb: float) -> dict:
    """IAU 2015 orientation: pole right ascension and declination, prime meridian angle W
    (degrees) and rotation rate (degrees per day, negative for retrograde rotators)."""
    ra, dec, w, rate = _native.space_orientation(world, jd_tdb)
    return {"pole_ra": ra, "pole_dec": dec, "prime_meridian": w, "rotation_rate": rate}


def radii(world: str) -> tuple:
    """(equatorial, polar) radius of the world's reference ellipsoid, km."""
    return _native.space_radii(world)


def mars_sol_date(jd_tt: float) -> float:
    """Mars Sol Date (Mars24 convention) for a JD in TT."""
    return _native.mars_sol_date(jd_tt)


def mars_time(jd_tt: float, lon: float = 0.0) -> dict:
    """Mars Sol Date, Coordinated Mars Time and local mean solar time at east longitude
    ``lon`` (hours) for a JD in TT."""
    return {
        "mars_sol_date": _native.mars_sol_date(jd_tt),
        "coordinated_mars_time": _native.mars_coordinated_time(jd_tt),
        "local_mean_solar_time": _native.mars_local_mean_time(jd_tt, lon),
    }


def jd_from_mars_time(sol: int, hours: float, lon: float = 0.0) -> float:
    """JD TT for local mean solar time ``hours`` on Mars Sol Date number ``sol`` at east
    longitude ``lon``. Convert to UT with the engine before calling the chart functions."""
    msd = sol + ((hours - lon / 15.0) / 24.0)
    return _native.jd_tt_from_mars_sol_date(msd)
