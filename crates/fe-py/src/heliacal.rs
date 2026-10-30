// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Heliacal-event bindings: module-level functions taking the `Engine`, with the
//! classic argument shapes (`geopos`, `datm[4]`, `dobs[6]`, object name, event
//! type, `HELFLAG_*` bits) so `frenchephemeris.swe_compat` forwards them unchanged.

use super::{to_py, Engine};
use fe_core::ephem::Observer;
use fe_core::heliacal::{HelAtmosphere, HelEvent, HelFlags, HelObserver, HelTarget};
use pyo3::prelude::*;

pub const HELFLAG_LONG_SEARCH: i32 = 128;
pub const HELFLAG_HIGH_PRECISION: i32 = 256;
pub const HELFLAG_OPTICAL_PARAMS: i32 = 512;
pub const HELFLAG_NO_DETAILS: i32 = 1024;
pub const HELFLAG_SEARCH_1_PERIOD: i32 = 2048;
pub const HELFLAG_VISLIM_DARK: i32 = 4096;
pub const HELFLAG_VISLIM_NOMOON: i32 = 8192;
pub const HELFLAG_VISLIM_PHOTOPIC: i32 = 16384;
pub const HELFLAG_AV: i32 = 65536;

fn observer(g: (f64, f64, f64)) -> Observer {
    Observer { lon_deg: g.0, lat_deg: g.1, height_m: g.2 }
}

fn atmosphere(a: &[f64]) -> HelAtmosphere {
    let g = |i: usize| a.get(i).copied().unwrap_or(0.0);
    HelAtmosphere { pressure_hpa: g(0), temperature_c: g(1), humidity_pct: g(2), vr_or_k: g(3) }
}

fn hel_observer(o: &[f64]) -> HelObserver {
    let g = |i: usize| o.get(i).copied().unwrap_or(0.0);
    HelObserver {
        age: g(0),
        snellen: g(1),
        binocular: g(2) != 0.0,
        magnification: g(3),
        aperture_mm: g(4),
        transmission: g(5),
    }
    .resolved()
}

fn hel_flags(f: i32) -> HelFlags {
    HelFlags {
        optical: f & HELFLAG_OPTICAL_PARAMS != 0,
        dark: f & HELFLAG_VISLIM_DARK != 0,
        no_moon: f & HELFLAG_VISLIM_NOMOON != 0,
        photopic: f & HELFLAG_VISLIM_PHOTOPIC != 0,
    }
}

/// Next heliacal event: (first visible, best, last visible) JD UT of the event
/// day (`HELFLAG_NO_DETAILS`: (first visible, 0, 0)). Event types: 1 morning
/// first (heliacal rising), 2 evening last (heliacal setting), 3 evening first,
/// 4 morning last (3/4: Mercury, Venus, Moon; the Moon has only 3/4).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, geopos, atmo, observer_, objname, eventtype, flags=0))]
#[allow(clippy::too_many_arguments)]
fn heliacal_ut(
    py: Python<'_>,
    engine: &Bound<'_, Engine>,
    jd_ut: f64,
    geopos: (f64, f64, f64),
    atmo: Vec<f64>,
    observer_: Vec<f64>,
    objname: &str,
    eventtype: i32,
    flags: i32,
) -> PyResult<(f64, f64, f64)> {
    let target = HelTarget::parse(objname).map_err(to_py)?;
    let ev = HelEvent::from_id(eventtype).map_err(to_py)?;
    let eng = &engine.get().inner;
    let (obs, atm, ho, fl) =
        (observer(geopos), atmosphere(&atmo), hel_observer(&observer_), hel_flags(flags));
    let r = py.detach(|| eng.heliacal_ut(jd_ut, &obs, &atm, &ho, target, ev, fl)).map_err(to_py)?;
    if flags & HELFLAG_NO_DETAILS != 0 {
        return Ok((r[0], 0.0, 0.0));
    }
    Ok((r[0], r[1], r[2]))
}

/// Heliacal phenomena at an instant: 50 floats (see fe_core::heliacal).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, geopos, atmo, observer_, objname, eventtype, flags=0))]
#[allow(clippy::too_many_arguments)]
fn heliacal_pheno_ut(
    py: Python<'_>,
    engine: &Bound<'_, Engine>,
    jd_ut: f64,
    geopos: (f64, f64, f64),
    atmo: Vec<f64>,
    observer_: Vec<f64>,
    objname: &str,
    eventtype: i32,
    flags: i32,
) -> PyResult<Vec<f64>> {
    let target = HelTarget::parse(objname).map_err(to_py)?;
    let ev = HelEvent::from_id(eventtype).map_err(to_py)?;
    let eng = &engine.get().inner;
    let (obs, atm, ho, fl) =
        (observer(geopos), atmosphere(&atmo), hel_observer(&observer_), hel_flags(flags));
    let r = py.detach(|| eng.heliacal_pheno_ut(jd_ut, &obs, &atm, &ho, target, ev, fl)).map_err(to_py)?;
    Ok(r.to_vec())
}

/// Limiting visual magnitude at the object's position: (res, dret[10]);
/// res −2 below the horizon, 0 photopic, 1 scotopic.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, geopos, atmo, observer_, objname, flags=0))]
#[allow(clippy::too_many_arguments)]
fn vis_limit_mag(
    py: Python<'_>,
    engine: &Bound<'_, Engine>,
    jd_ut: f64,
    geopos: (f64, f64, f64),
    atmo: Vec<f64>,
    observer_: Vec<f64>,
    objname: &str,
    flags: i32,
) -> PyResult<(i32, [f64; 10])> {
    let target = HelTarget::parse(objname).map_err(to_py)?;
    let eng = &engine.get().inner;
    let (obs, atm, ho, fl) =
        (observer(geopos), atmosphere(&atmo), hel_observer(&observer_), hel_flags(flags));
    let r = py.detach(|| eng.vis_limit_mag(jd_ut, &obs, &atm, &ho, target, fl)).map_err(to_py)?;
    Ok((r.res, r.dret))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(heliacal_ut, m)?)?;
    m.add_function(wrap_pyfunction!(heliacal_pheno_ut, m)?)?;
    m.add_function(wrap_pyfunction!(vis_limit_mag, m)?)?;
    for (k, v) in [
        ("HELIACAL_RISING", 1),
        ("MORNING_FIRST", 1),
        ("HELIACAL_SETTING", 2),
        ("EVENING_LAST", 2),
        ("EVENING_FIRST", 3),
        ("MORNING_LAST", 4),
        ("ACRONYCHAL_RISING", 5),
        ("ACRONYCHAL_SETTING", 6),
        ("COSMICAL_SETTING", 6),
        ("HELFLAG_LONG_SEARCH", HELFLAG_LONG_SEARCH),
        ("HELFLAG_HIGH_PRECISION", HELFLAG_HIGH_PRECISION),
        ("HELFLAG_OPTICAL_PARAMS", HELFLAG_OPTICAL_PARAMS),
        ("HELFLAG_NO_DETAILS", HELFLAG_NO_DETAILS),
        ("HELFLAG_SEARCH_1_PERIOD", HELFLAG_SEARCH_1_PERIOD),
        ("HELFLAG_VISLIM_DARK", HELFLAG_VISLIM_DARK),
        ("HELFLAG_VISLIM_NOMOON", HELFLAG_VISLIM_NOMOON),
        ("HELFLAG_VISLIM_PHOTOPIC", HELFLAG_VISLIM_PHOTOPIC),
        ("HELFLAG_AV", HELFLAG_AV),
        ("PHOTOPIC_FLAG", 0),
        ("SCOTOPIC_FLAG", 1),
        ("MIXEDOPIC_FLAG", 2),
    ] {
        m.add(k, v)?;
    }
    Ok(())
}
