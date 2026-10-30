// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Fixed-star bindings: module-level functions taking the `Engine` handle
//! (a second `#[pymethods]` block would need PyO3's `multiple-pymethods`).
//!
//! Return shapes follow the classic `fixstar2_ut` convention used by the
//! FreeAstroAPI adapter: `(xx[6], "Name,designation")`. Names and returned strings
//! follow the engine's naming: the public default catalogue, or the optional
//! compatibility alias table given as `Engine(..., star_compat=path)`.

use super::{options, to_py, Engine, FLG_SPEED};
use fe_core::context::Position;
use fe_core::stars;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// (xx, name) for one star at JD UT. `xx` = lon/RA, lat/Dec (deg), distance
/// (AU), and their daily rates; `name` = "Name,designation", e.g. "Spica,alf Vir".
/// Names are matched case-insensitively, ignoring spaces; ",alf Vir", "α Vir",
/// "HR 5056", "HD 116658" and "HIP 65474" are accepted too. ValueError if unknown.
#[pyfunction]
#[pyo3(signature = (engine, name, jd_ut, flags=FLG_SPEED, sidereal=None, topo=None))]
fn fixstar_ut(
    py: Python<'_>,
    engine: &Bound<'_, Engine>,
    name: &str,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
    topo: Option<(f64, f64, f64)>,
) -> PyResult<([f64; 6], String)> {
    let opt = options(flags, sidereal, topo)?;
    let eng = &engine.get().inner;
    let (n, p) = py
        .detach(|| match super::horizons_mode(flags) {
            Some(m) => {
                let (star, full) = eng.fixstar_lookup(name)?;
                let prep = eng.prepare_horizons(jd_ut, true, opt.speed, m);
                Ok((full, eng.fixstar_prepared(&star, &prep, &opt)?))
            }
            None => eng.fixstar_ut(name, jd_ut, &opt),
        })
        .map_err(to_py)?;
    Ok((p.to_array(), n))
}

/// Many stars at one JD UT (instant prepared once, runs without the GIL).
/// Returns a list aligned with `names`: `(xx, name)` or `(None, error message)`.
#[pyfunction]
#[pyo3(signature = (engine, names, jd_ut, flags=FLG_SPEED, sidereal=None, topo=None))]
fn fixstars_ut(
    py: Python<'_>,
    engine: &Bound<'_, Engine>,
    names: Vec<String>,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
    topo: Option<(f64, f64, f64)>,
) -> PyResult<Vec<(Option<[f64; 6]>, String)>> {
    let opt = options(flags, sidereal, topo)?;
    let eng = &engine.get().inner;
    let refs: Vec<&str> = names.iter().map(String::as_str).collect();
    let out = py.detach(|| eng.fixstars_ut(&refs, jd_ut, &opt));
    Ok(out
        .into_iter()
        .map(|r| match r {
            Ok((n, p)) => (Some(Position::to_array(p)), n),
            Err(e) => (None, e.to_string()),
        })
        .collect())
}

/// (magnitude, "Name,designation") for a star, like the classic `fixstar2_mag`:
/// Johnson V (SIMBAD), 999.99 where no V magnitude is published. ValueError if unknown.
#[pyfunction]
fn fixstar_mag(engine: &Bound<'_, Engine>, name: &str) -> PyResult<(f64, String)> {
    let (n, m) = engine.get().inner.fixstar_mag(name).map_err(to_py)?;
    Ok((m, n))
}

/// Many names: list of `(mag, name)` or `(None, error message)`.
#[pyfunction]
fn fixstars_mag(engine: &Bound<'_, Engine>, names: Vec<String>) -> Vec<(Option<f64>, String)> {
    let refs: Vec<&str> = names.iter().map(String::as_str).collect();
    engine
        .get()
        .inner
        .fixstars_mag(&refs)
        .into_iter()
        .map(|r| match r {
            Ok((n, m)) => (Some(m), n),
            Err(e) => (None, e.to_string()),
        })
        .collect()
}

/// Every accepted star name (traditional / IAU names and aliases).
#[pyfunction]
fn fixstar_names() -> Vec<&'static str> {
    let mut v: Vec<&'static str> = stars::names().map(|(n, _)| n).filter(|n| !n.is_empty()).collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Catalogue data for a star: name, nomenclature, designation, hip,
/// ra_deg/dec_deg (ICRS, epoch J1991.25), pm (mas/yr), parallax (mas), rv (km/s), vmag.
#[pyfunction]
fn fixstar_info<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let s = stars::lookup(name).map_err(to_py)?;
    let r = s.record;
    let d = PyDict::new(py);
    d.set_item("name", s.full_name())?;
    d.set_item("primary_name", r.name)?;
    d.set_item("nomenclature", r.nomenclature)?;
    d.set_item("designation", r.designation)?;
    d.set_item("hip", r.hip)?;
    d.set_item("ra_deg", r.ra_rad.to_degrees())?;
    d.set_item("dec_deg", r.dec_rad.to_degrees())?;
    d.set_item("pmra_mas_yr", r.pmra_mas_yr)?;
    d.set_item("pmdec_mas_yr", r.pmdec_mas_yr)?;
    d.set_item("plx_mas", r.plx_mas)?;
    d.set_item("rv_km_s", r.rv_km_s)?;
    d.set_item("vmag", r.vmag)?;
    d.set_item("epoch_jd_tt", stars::CATALOG_EPOCH_JD)?;
    Ok(d)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fixstar_ut, m)?)?;
    m.add_function(wrap_pyfunction!(fixstars_ut, m)?)?;
    m.add_function(wrap_pyfunction!(fixstar_names, m)?)?;
    m.add_function(wrap_pyfunction!(fixstar_info, m)?)?;
    m.add_function(wrap_pyfunction!(fixstar_mag, m)?)?;
    m.add_function(wrap_pyfunction!(fixstars_mag, m)?)?;
    Ok(())
}
