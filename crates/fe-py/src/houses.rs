// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Python bindings for the house systems (every Swiss code incl. 'J'), speeds, house positions, house names and Gauquelin sectors.
//! Module-level functions (the dated ones take the `Engine`), shapes as in
//! the classic API so `frenchephemeris.swe_compat` forwards them unchanged.

use crate::{options, parse_sidereal, to_py, Engine, FLG_SIDEREAL};
use fe_core::context::Body;
use fe_core::ephem::Observer;
use fe_core::events::Atmosphere;
use fe_core::houses::{self, GauquelinMethod, HouseSpeeds};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn system(hsys: &str) -> PyResult<char> {
    hsys.chars().next().ok_or_else(|| PyValueError::new_err("empty hsys"))
}

type Ex2 = (Vec<f64>, [f64; 8], Vec<f64>, [f64; 8]);

fn ex2(h: HouseSpeeds) -> Ex2 {
    (h.houses.cusps, h.houses.ascmc, h.cusp_speeds, h.ascmc_speeds)
}

/// Houses from ARMC (deg), geographic latitude and obliquity: (cusps, ascmc[8]);
/// 12 cusps, 36 for 'G'. `ascmc9` is the Sun's declination for Sunshine ('I', 'i').
/// Placidus/Koch/Gauquelin inside the polar circle (and 'i' with a circumpolar Sun)
/// raise PolarHousesError unless `swiss_fallback=True` (Porphyry, as the C library).
#[pyfunction]
#[pyo3(signature = (armc, lat, eps, hsys="P", ascmc9=0.0, swiss_fallback=false))]
fn houses_armc(
    armc: f64,
    lat: f64,
    eps: f64,
    hsys: &str,
    ascmc9: f64,
    swiss_fallback: bool,
) -> PyResult<(Vec<f64>, [f64; 8])> {
    let s = system(hsys)?;
    let h = if swiss_fallback {
        houses::house_cusps_swiss_compat(s, armc, lat, eps, Some(ascmc9)).map(|x| x.0)
    } else {
        houses::house_cusps(s, armc, lat, eps, Some(ascmc9))
    }
    .map_err(to_py)?;
    Ok((h.cusps, h.ascmc))
}

/// `houses_armc` plus speeds (deg/day): (cusps, ascmc, cusp_speeds, ascmc_speeds).
#[pyfunction]
#[pyo3(signature = (armc, lat, eps, hsys="P", ascmc9=0.0, swiss_fallback=false))]
fn houses_armc_ex2(
    armc: f64,
    lat: f64,
    eps: f64,
    hsys: &str,
    ascmc9: f64,
    swiss_fallback: bool,
) -> PyResult<Ex2> {
    let s = system(hsys)?;
    let h = if swiss_fallback {
        houses::house_cusps_speeds_swiss_compat(s, armc, lat, eps, Some(ascmc9)).map(|x| x.0)
    } else {
        houses::house_cusps_speeds(s, armc, lat, eps, Some(ascmc9))
    }
    .map_err(to_py)?;
    Ok(ex2(h))
}

/// Houses with speeds at JD UT and place: (cusps, ascmc, cusp_speeds, ascmc_speeds).
/// Sunshine uses the Sun's apparent declination of the date. `sidereal`: a
/// sidereal spec as for `Engine.calc_ut` (name, Swiss number, `"lahiri|256"`);
/// projection modes intersect the houses with the plane of t0 (oracle semantics).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, lat, lon, hsys="P", sidereal=None, swiss_fallback=false))]
#[allow(clippy::too_many_arguments)]
fn houses_ex2(
    engine: &Engine,
    jd_ut: f64,
    lat: f64,
    lon: f64,
    hsys: &str,
    sidereal: Option<&str>,
    swiss_fallback: bool,
) -> PyResult<Ex2> {
    let s = system(hsys)?;
    let sid = sidereal.map(parse_sidereal).transpose()?;
    // A few microseconds: keep the GIL (releasing it costs more and makes threads contend).
    let (h, _) =
        engine.inner.house_cusps_ut_bits(jd_ut, lat, lon, s, sid, swiss_fallback, true).map_err(to_py)?;
    Ok(ex2(h))
}

/// Continuous house position in [1, 13) (Gauquelin: [1, 37), clockwise) of an
/// ecliptic point. `sun_dec` is required for Sunshine ('I', 'i').
#[pyfunction]
#[pyo3(signature = (armc, lat, eps, lon, lat_ecl=0.0, hsys="P", sun_dec=None))]
#[allow(clippy::too_many_arguments)]
fn house_pos(
    armc: f64,
    lat: f64,
    eps: f64,
    lon: f64,
    lat_ecl: f64,
    hsys: &str,
    sun_dec: Option<f64>,
) -> PyResult<f64> {
    houses::house_position_ex(system(hsys)?, armc, lat, eps, lon, lat_ecl, sun_dec).map_err(to_py)
}

/// Name of a house system ("" if unknown).
#[pyfunction]
fn house_name(hsys: &str) -> PyResult<&'static str> {
    Ok(houses::house_name(system(hsys)?))
}

/// Gauquelin sector (1 ≤ s < 37) of a body at JD UT; `method` 0–5 as the classic
/// `gauquelin_sector`; `flags` select the place for methods 0/1 (FLG_TOPOCTR uses
/// `geopos` as the observer). Raises for bodies that do not rise and set (2–5).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body, method, geopos, atpress=0.0, attemp=0.0, flags=crate::FLG_TOPOCTR))]
#[allow(clippy::too_many_arguments)]
fn gauquelin_sector(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: i32,
    method: i32,
    geopos: (f64, f64, f64),
    atpress: f64,
    attemp: f64,
    flags: i32,
) -> PyResult<f64> {
    let body = Body::from_id(body).map_err(to_py)?;
    let m = GauquelinMethod::from_code(method).map_err(to_py)?;
    let opt = options(flags & !FLG_SIDEREAL, None, Some(geopos))?;
    let obs = Observer { lon_deg: geopos.0, lat_deg: geopos.1, height_m: geopos.2 };
    let atm = Atmosphere { pressure_hpa: atpress, temperature_c: attemp };
    py.detach(|| engine.inner.gauquelin_sector(jd_ut, body, m, &obs, atm, &opt)).map_err(to_py)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(houses_armc, m)?)?;
    m.add_function(wrap_pyfunction!(houses_armc_ex2, m)?)?;
    m.add_function(wrap_pyfunction!(houses_ex2, m)?)?;
    m.add_function(wrap_pyfunction!(house_pos, m)?)?;
    m.add_function(wrap_pyfunction!(house_name, m)?)?;
    m.add_function(wrap_pyfunction!(gauquelin_sector, m)?)?;
    Ok(())
}
