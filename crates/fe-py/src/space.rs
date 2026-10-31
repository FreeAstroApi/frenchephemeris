// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Bindings for the `space` extension: observers on other worlds or anywhere in space.
//! Low-level functions taking the `Engine`; `frenchephemeris.ext.space` builds the public
//! Python API on top of them.

use super::{to_py, Engine};
use fe_core::context::Body;
use fe_core::ephem::Corrections;
use fe_core::space::{self, Site, SpaceObserver, World, Zodiac};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn world(name: &str) -> PyResult<World> {
    World::from_name(name).map_err(to_py)
}

fn zodiac(name: &str) -> PyResult<Zodiac> {
    match name.trim().to_ascii_lowercase().as_str() {
        "earth" | "tropical" => Ok(Zodiac::Earth),
        "native" => Ok(Zodiac::Native),
        _ => Err(PyValueError::new_err(format!("unknown zodiac {name:?}: use 'earth' or 'native'"))),
    }
}

fn site(world_name: &str, lon: f64, lat: f64, height_km: f64) -> PyResult<Site> {
    Ok(Site { world: world(world_name)?, lon_deg: lon, lat_deg: lat, height_km })
}

/// Observer from flat arguments: kind "surface" (world, lon, lat, height_km), "center"
/// (world) or "point" (center NAIF id, offset km, velocity km/day, ICRF axes).
#[allow(clippy::too_many_arguments)]
fn observer(
    kind: &str,
    world_name: Option<&str>,
    lon: f64,
    lat: f64,
    height_km: f64,
    center: i32,
    offset: [f64; 3],
    velocity: [f64; 3],
) -> PyResult<SpaceObserver> {
    let need_world =
        || world_name.ok_or_else(|| PyValueError::new_err(format!("observer kind {kind:?} needs a world")));
    match kind {
        "surface" => Ok(SpaceObserver::Surface(site(need_world()?, lon, lat, height_km)?)),
        "center" => Ok(SpaceObserver::Center(world(need_world()?)?)),
        "point" => Ok(SpaceObserver::Point { center, offset_km: offset, velocity_km_day: velocity }),
        _ => Err(PyValueError::new_err(format!("unknown observer kind {kind:?}"))),
    }
}

/// Position of `body_id` seen by the observer, at JD UT: (lon, lat, dist AU, lon speed,
/// lat speed, dist speed) in the chosen zodiac. `apparent=False` gives the astrometric place.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body_id, kind, world=None, lon=0.0, lat=0.0, height_km=0.0,
    center=0, offset=[0.0; 3], velocity=[0.0; 3], zodiac="earth", apparent=true))]
#[allow(clippy::too_many_arguments)]
fn space_calc_ut(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body_id: i32,
    kind: &str,
    world: Option<&str>,
    lon: f64,
    lat: f64,
    height_km: f64,
    center: i32,
    offset: [f64; 3],
    velocity: [f64; 3],
    zodiac: &str,
    apparent: bool,
) -> PyResult<[f64; 6]> {
    let obs = observer(kind, world, lon, lat, height_km, center, offset, velocity)?;
    let (body, z) = (Body::from_id(body_id).map_err(to_py)?, self::zodiac(zodiac)?);
    let corr = if apparent { Corrections::APPARENT } else { Corrections::ASTROMETRIC };
    let p = py.detach(|| engine.inner.space_position(body, jd_ut, &obs, z, corr)).map_err(to_py)?;
    Ok([p.lon, p.lat, p.dist_au, p.lon_speed, p.lat_speed, p.dist_speed])
}

/// (cusps[12], ascmc[8], sidereal angle, spin latitude, obliquity).
type SpaceHousesTuple = ([f64; 12], [f64; 8], f64, f64, f64);

/// Houses on another world: (cusps[12], ascmc[8], sidereal angle, spin latitude, obliquity).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, world, lon, lat, height_km=0.0, hsys="P", zodiac="earth"))]
#[allow(clippy::too_many_arguments)]
fn space_houses_ut(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    world: &str,
    lon: f64,
    lat: f64,
    height_km: f64,
    hsys: &str,
    zodiac: &str,
) -> PyResult<SpaceHousesTuple> {
    let s = site(world, lon, lat, height_km)?;
    let mut it = hsys.chars();
    let system = match (it.next(), it.next()) {
        (Some(c), None) => c,
        _ => return Err(PyValueError::new_err("hsys must be exactly one character")),
    };
    let z = self::zodiac(zodiac)?;
    let h = py.detach(|| engine.inner.space_houses(jd_ut, &s, system, z)).map_err(to_py)?;
    Ok((h.houses.cusps, h.houses.ascmc(), h.sidereal_angle, h.spin_latitude, h.obliquity))
}

/// Azimuth (from the world's north towards east) and elevation of `body_id`, degrees.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body_id, world, lon, lat, height_km=0.0))]
#[allow(clippy::too_many_arguments)]
fn space_azalt_ut(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body_id: i32,
    world: &str,
    lon: f64,
    lat: f64,
    height_km: f64,
) -> PyResult<(f64, f64)> {
    let s = site(world, lon, lat, height_km)?;
    let body = Body::from_id(body_id).map_err(to_py)?;
    let h = py.detach(|| engine.inner.space_horizontal(body, jd_ut, &s)).map_err(to_py)?;
    Ok((h.azimuth, h.elevation))
}

/// IAU 2015 orientation of a world at JD TDB: (pole RA, pole Dec, prime meridian W,
/// rotation rate deg/day).
#[pyfunction]
fn space_orientation(world_name: &str, jd_tdb: f64) -> PyResult<(f64, f64, f64, f64)> {
    let o = world(world_name)?.orientation(jd_tdb);
    Ok((o.pole_ra, o.pole_dec, o.prime_meridian, o.rotation_rate))
}

/// (equatorial, polar) radius of a world's reference ellipsoid, km.
#[pyfunction]
fn space_radii(world_name: &str) -> PyResult<(f64, f64)> {
    Ok(world(world_name)?.radii_km())
}

/// Names of the worlds an observer can stand on.
#[pyfunction]
fn space_worlds() -> Vec<&'static str> {
    World::ALL.iter().map(|w| w.name()).collect()
}

#[pyfunction]
fn mars_sol_date(jd_tt: f64) -> f64 {
    space::mars_sol_date(jd_tt)
}

#[pyfunction]
fn jd_tt_from_mars_sol_date(msd: f64) -> f64 {
    space::jd_tt_from_mars_sol_date(msd)
}

#[pyfunction]
fn mars_coordinated_time(jd_tt: f64) -> f64 {
    space::mars_coordinated_time(jd_tt)
}

#[pyfunction]
fn mars_local_mean_time(jd_tt: f64, lon: f64) -> f64 {
    space::mars_local_mean_time(jd_tt, lon)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(space_calc_ut, m)?)?;
    m.add_function(wrap_pyfunction!(space_houses_ut, m)?)?;
    m.add_function(wrap_pyfunction!(space_azalt_ut, m)?)?;
    m.add_function(wrap_pyfunction!(space_orientation, m)?)?;
    m.add_function(wrap_pyfunction!(space_radii, m)?)?;
    m.add_function(wrap_pyfunction!(space_worlds, m)?)?;
    m.add_function(wrap_pyfunction!(mars_sol_date, m)?)?;
    m.add_function(wrap_pyfunction!(jd_tt_from_mars_sol_date, m)?)?;
    m.add_function(wrap_pyfunction!(mars_coordinated_time, m)?)?;
    m.add_function(wrap_pyfunction!(mars_local_mean_time, m)?)?;
    Ok(())
}
