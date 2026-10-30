// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Python bindings for fe_core::events: module-level functions taking the
//! `Engine` (no second `#[pymethods]` block, so no PyO3 feature change).
//! Argument and return shapes follow the classic API so
//! `frenchephemeris.swe_compat` can forward them unchanged.

use crate::{options, to_py, Engine};
use fe_core::context::Body;
use fe_core::ephem::Observer;
use fe_core::events::{
    self, Atmosphere, AzaltInput, DiscPoint, RefracDir, RefractionModel, RiseSet, RiseSetKind, RiseSetOptions,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub const CALC_RISE: i32 = 1;
pub const CALC_SET: i32 = 2;
pub const CALC_MTRANSIT: i32 = 4;
pub const CALC_ITRANSIT: i32 = 8;
pub const BIT_GEOCTR_NO_ECL_LAT: i32 = 128;
pub const BIT_DISC_CENTER: i32 = 256;
pub const BIT_NO_REFRACTION: i32 = 512;
pub const BIT_CIVIL_TWILIGHT: i32 = 1024;
pub const BIT_NAUTIC_TWILIGHT: i32 = 2048;
pub const BIT_ASTRO_TWILIGHT: i32 = 4096;
pub const BIT_DISC_BOTTOM: i32 = 8192;
pub const BIT_FIXED_DISC_SIZE: i32 = 16384;
pub const BIT_FORCE_SLOW_METHOD: i32 = 32768;
pub const BIT_HINDU_RISING: i32 = BIT_DISC_CENTER | BIT_NO_REFRACTION | BIT_GEOCTR_NO_ECL_LAT;
pub const ECL2HOR: i32 = 0;
pub const EQU2HOR: i32 = 1;
pub const TRUE_TO_APP: i32 = 0;
pub const APP_TO_TRUE: i32 = 1;

fn observer(geopos: (f64, f64, f64)) -> Observer {
    Observer { lon_deg: geopos.0, lat_deg: geopos.1, height_m: geopos.2 }
}

fn rise_options(
    rsmi: i32,
    atpress: f64,
    attemp: f64,
    horhgt: f64,
) -> PyResult<(RiseSetKind, RiseSetOptions)> {
    let kind = match rsmi & 15 {
        CALC_RISE => RiseSetKind::Rise,
        CALC_SET => RiseSetKind::Set,
        CALC_MTRANSIT => RiseSetKind::UpperTransit,
        CALC_ITRANSIT => RiseSetKind::LowerTransit,
        k => {
            return Err(PyValueError::new_err(format!(
                "rsmi must select exactly one of rise/set/transit, got {k}"
            )))
        }
    };
    let mut o = RiseSetOptions {
        refraction: Some(Atmosphere { pressure_hpa: atpress, temperature_c: attemp }),
        horizon_deg: horhgt,
        ..Default::default()
    };
    if rsmi & BIT_DISC_CENTER != 0 {
        o.disc = DiscPoint::Center;
    }
    if rsmi & BIT_DISC_BOTTOM != 0 {
        o.disc = DiscPoint::LowerLimb;
    }
    if rsmi & BIT_NO_REFRACTION != 0 {
        o.refraction = None;
    }
    o.geocentric_no_lat = rsmi & BIT_GEOCTR_NO_ECL_LAT != 0;
    o.fixed_disc_size = rsmi & BIT_FIXED_DISC_SIZE != 0;
    // Twilight: disc centre at a true altitude of −6/−12/−18°, no refraction (measured).
    for (bit, h) in [(BIT_CIVIL_TWILIGHT, -6.0), (BIT_NAUTIC_TWILIGHT, -12.0), (BIT_ASTRO_TWILIGHT, -18.0)] {
        if rsmi & bit != 0 {
            o.disc = DiscPoint::Center;
            o.refraction = None;
            o.horizon_deg = h;
        }
    }
    Ok((kind, o))
}

/// Next rise / set / transit after `jd_ut`.
/// Returns (res, jd_ut, category): res 0 with the event JD, or −2 with jd 0.0
/// and category "always_above" / "always_below" (circumpolar, polar day/night).
/// Ephemeris flags are not needed: the place is always the apparent
/// topocentric one (or geocentric for BIT_GEOCTR_NO_ECL_LAT).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body, rsmi, geopos, atpress=0.0, attemp=0.0, horhgt=0.0))]
#[allow(clippy::too_many_arguments)]
fn rise_trans(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: i32,
    rsmi: i32,
    geopos: (f64, f64, f64),
    atpress: f64,
    attemp: f64,
    horhgt: f64,
) -> PyResult<(i32, f64, &'static str)> {
    let body = Body::from_id(body).map_err(to_py)?;
    let (kind, o) = rise_options(rsmi, atpress, attemp, horhgt)?;
    let obs = observer(geopos);
    let r = py.detach(|| engine.inner.rise_set(body, jd_ut, kind, &obs, &o)).map_err(to_py)?;
    Ok(match r {
        RiseSet::Event(t) => (0, t, "event"),
        RiseSet::AlwaysAbove => (-2, 0.0, "always_above"),
        RiseSet::AlwaysBelow => (-2, 0.0, "always_below"),
    })
}

/// Planetary phenomena at JD UT: 20 floats, [phase angle, phase, elongation,
/// apparent diameter, magnitude, horizontal parallax (Moon), 0...].
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body, flags=0, topo=None))]
fn pheno_ut(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: i32,
    flags: i32,
    topo: Option<(f64, f64, f64)>,
) -> PyResult<[f64; 20]> {
    let opt = options(flags & !crate::FLG_SIDEREAL, None, topo)?;
    let body = Body::from_id(body).map_err(to_py)?;
    py.detach(|| engine.inner.pheno_ut(body, jd_ut, &opt)).map(|p| p.to_array()).map_err(to_py)
}

/// `pheno_ut` with JD in TT.
#[pyfunction]
#[pyo3(signature = (engine, jd_tt, body, flags=0, topo=None))]
fn pheno(
    py: Python<'_>,
    engine: &Engine,
    jd_tt: f64,
    body: i32,
    flags: i32,
    topo: Option<(f64, f64, f64)>,
) -> PyResult<[f64; 20]> {
    let opt = options(flags & !crate::FLG_SIDEREAL, None, topo)?;
    let body = Body::from_id(body).map_err(to_py)?;
    py.detach(|| engine.inner.pheno_tt(body, jd_tt, &opt)).map(|p| p.to_array()).map_err(to_py)
}

fn cross(
    py: Python<'_>,
    engine: &Engine,
    body: Body,
    x: f64,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<f64> {
    let opt = options(flags, sidereal, None)?;
    py.detach(|| engine.inner.longitude_cross_ut(body, x, jd_ut, &opt)).map_err(to_py)
}

/// Next JD (UT) after `jd_ut` when the Sun's longitude is `x2cross`.
#[pyfunction]
#[pyo3(signature = (engine, x2cross, jd_ut, flags=0, sidereal=None))]
fn solcross_ut(
    py: Python<'_>,
    engine: &Engine,
    x2cross: f64,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<f64> {
    cross(py, engine, Body::Sun, x2cross, jd_ut, flags, sidereal)
}

/// TT variant of `solcross_ut` (input and result in TT).
#[pyfunction]
#[pyo3(signature = (engine, x2cross, jd_tt, flags=0, sidereal=None))]
fn solcross(
    py: Python<'_>,
    engine: &Engine,
    x2cross: f64,
    jd_tt: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<f64> {
    cross(py, engine, Body::Sun, x2cross, engine.inner.tt_to_ut(jd_tt), flags, sidereal)
        .map(|t| engine.inner.ut_to_tt(t))
}

/// Next JD (UT) after `jd_ut` when the Moon's longitude is `x2cross`.
#[pyfunction]
#[pyo3(signature = (engine, x2cross, jd_ut, flags=0, sidereal=None))]
fn mooncross_ut(
    py: Python<'_>,
    engine: &Engine,
    x2cross: f64,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<f64> {
    cross(py, engine, Body::Moon, x2cross, jd_ut, flags, sidereal)
}

/// TT variant of `mooncross_ut`.
#[pyfunction]
#[pyo3(signature = (engine, x2cross, jd_tt, flags=0, sidereal=None))]
fn mooncross(
    py: Python<'_>,
    engine: &Engine,
    x2cross: f64,
    jd_tt: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<f64> {
    cross(py, engine, Body::Moon, x2cross, engine.inner.tt_to_ut(jd_tt), flags, sidereal)
        .map(|t| engine.inner.ut_to_tt(t))
}

/// Equation of time (days): local apparent − local mean time.
#[pyfunction]
fn time_equ(engine: &Engine, jd_ut: f64) -> PyResult<f64> {
    engine.inner.time_equ(jd_ut).map_err(to_py)
}

#[pyfunction]
fn lmt_to_lat(engine: &Engine, jd_lmt: f64, geolon: f64) -> PyResult<f64> {
    engine.inner.lmt_to_lat(jd_lmt, geolon).map_err(to_py)
}

#[pyfunction]
fn lat_to_lmt(engine: &Engine, jd_lat: f64, geolon: f64) -> PyResult<f64> {
    engine.inner.lat_to_lmt(jd_lat, geolon).map_err(to_py)
}

/// True solar (sundial) time as a JD-like value: UT + lon/360 + E(UT).
#[pyfunction]
fn true_solar_time(engine: &Engine, jd_ut: f64, geolon: f64) -> PyResult<f64> {
    engine.inner.true_solar_time(jd_ut, geolon).map_err(to_py)
}

/// Ecliptic ↔ equatorial rotation: eps < 0 ecliptic→equatorial, eps > 0 back.
#[pyfunction]
fn cotrans(coord: (f64, f64, f64), eps: f64) -> (f64, f64, f64) {
    let r = events::cotrans([coord.0, coord.1, coord.2], eps);
    (r[0], r[1], r[2])
}

/// `cotrans` for position and speed (6 values).
#[pyfunction]
fn cotrans_sp(coord: [f64; 6], eps: f64) -> [f64; 6] {
    events::cotrans_sp(coord, eps)
}

/// (azimuth from south westwards, true altitude, apparent altitude), degrees.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, flag, geopos, atpress, attemp, xin))]
fn azalt(
    engine: &Engine,
    jd_ut: f64,
    flag: i32,
    geopos: (f64, f64, f64),
    atpress: f64,
    attemp: f64,
    xin: Vec<f64>,
) -> PyResult<(f64, f64, f64)> {
    let input = match flag {
        ECL2HOR => AzaltInput::Ecliptic,
        EQU2HOR => AzaltInput::Equatorial,
        f => return Err(PyValueError::new_err(format!("azalt flag must be ECL2HOR or EQU2HOR, got {f}"))),
    };
    if xin.len() < 2 {
        return Err(PyValueError::new_err("xin needs (lon, lat[, dist])"));
    }
    let atm = Atmosphere { pressure_hpa: atpress, temperature_c: attemp };
    let r = engine.inner.azalt(
        jd_ut,
        input,
        &observer(geopos),
        &atm,
        [xin[0], xin[1], *xin.get(2).unwrap_or(&1.0)],
    );
    Ok((r[0], r[1], r[2]))
}

fn refraction_model(name: &str) -> PyResult<RefractionModel> {
    Ok(match name.to_ascii_lowercase().replace(['-', ' '], "_").as_str() {
        "sinclair" | "sinclair1982" => RefractionModel::Sinclair1982,
        "saemundsson" => RefractionModel::Saemundsson,
        "bennett" => RefractionModel::Bennett,
        "sinclair_simple" => RefractionModel::SinclairSimple,
        "almanac34" | "aa" => RefractionModel::Almanac34,
        other => return Err(PyValueError::new_err(format!("unknown refraction model {other:?}"))),
    })
}

/// Refraction conversion (TRUE_TO_APP / APP_TO_TRUE), degrees. `model`:
/// "sinclair1982" (default; Explanatory Supplement 1992 §3.28), "saemundsson",
/// "bennett", "sinclair_simple", "almanac34".
#[pyfunction]
#[pyo3(signature = (alt, atpress, attemp, flag, model="sinclair1982"))]
fn refrac(alt: f64, atpress: f64, attemp: f64, flag: i32, model: &str) -> PyResult<f64> {
    let dir = match flag {
        TRUE_TO_APP => RefracDir::TrueToApparent,
        APP_TO_TRUE => RefracDir::ApparentToTrue,
        f => {
            return Err(PyValueError::new_err(format!(
                "refrac flag must be TRUE_TO_APP or APP_TO_TRUE, got {f}"
            )))
        }
    };
    Ok(events::refrac_with(refraction_model(model)?, alt, atpress, attemp, dir))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rise_trans, m)?)?;
    m.add_function(wrap_pyfunction!(pheno_ut, m)?)?;
    m.add_function(wrap_pyfunction!(pheno, m)?)?;
    m.add_function(wrap_pyfunction!(solcross_ut, m)?)?;
    m.add_function(wrap_pyfunction!(solcross, m)?)?;
    m.add_function(wrap_pyfunction!(mooncross_ut, m)?)?;
    m.add_function(wrap_pyfunction!(mooncross, m)?)?;
    m.add_function(wrap_pyfunction!(time_equ, m)?)?;
    m.add_function(wrap_pyfunction!(lmt_to_lat, m)?)?;
    m.add_function(wrap_pyfunction!(lat_to_lmt, m)?)?;
    m.add_function(wrap_pyfunction!(true_solar_time, m)?)?;
    m.add_function(wrap_pyfunction!(cotrans, m)?)?;
    m.add_function(wrap_pyfunction!(cotrans_sp, m)?)?;
    m.add_function(wrap_pyfunction!(azalt, m)?)?;
    m.add_function(wrap_pyfunction!(refrac, m)?)?;
    for (k, v) in [
        ("CALC_RISE", CALC_RISE),
        ("CALC_SET", CALC_SET),
        ("CALC_MTRANSIT", CALC_MTRANSIT),
        ("CALC_ITRANSIT", CALC_ITRANSIT),
        ("BIT_GEOCTR_NO_ECL_LAT", BIT_GEOCTR_NO_ECL_LAT),
        ("BIT_DISC_CENTER", BIT_DISC_CENTER),
        ("BIT_NO_REFRACTION", BIT_NO_REFRACTION),
        ("BIT_CIVIL_TWILIGHT", BIT_CIVIL_TWILIGHT),
        ("BIT_NAUTIC_TWILIGHT", BIT_NAUTIC_TWILIGHT),
        ("BIT_ASTRO_TWILIGHT", BIT_ASTRO_TWILIGHT),
        ("BIT_DISC_BOTTOM", BIT_DISC_BOTTOM),
        ("BIT_FIXED_DISC_SIZE", BIT_FIXED_DISC_SIZE),
        ("BIT_FORCE_SLOW_METHOD", BIT_FORCE_SLOW_METHOD),
        ("BIT_HINDU_RISING", BIT_HINDU_RISING),
        ("ECL2HOR", ECL2HOR),
        ("EQU2HOR", EQU2HOR),
        ("TRUE_TO_APP", TRUE_TO_APP),
        ("APP_TO_TRUE", APP_TO_TRUE),
    ] {
        m.add(k, v)?;
    }
    m.add("RISE_SET_WINDOW_DAYS", events::RISE_SET_WINDOW_DAYS)?;
    Ok(())
}
