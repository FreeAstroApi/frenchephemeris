// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Bindings for the gap-analysis "quick wins": planet-centric positions, heliocentric and
//! node crossings, local lunar eclipses, UTC / leap seconds, Delta T controls (as derived
//! engines, never global state), calendar and formatting helpers, names. Module-level
//! functions taking the `Engine` where they need one; return shapes follow the classic API
//! so `swe_compat` only re-packs arguments.

use super::{options, to_py, Engine, FLG_BARYCTR, FLG_HELCTR, FLG_NOABERR, FLG_NOGDEFL, FLG_TOPOCTR};
use fe_core::context::Body;
use fe_core::ephem::Observer;
use fe_core::time::Calendar;
use fe_core::{lunar_local, utc, util};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub const ECL_PARTBEG_VISIBLE: i32 = 512;
pub const ECL_TOTBEG_VISIBLE: i32 = 1024;
pub const ECL_TOTEND_VISIBLE: i32 = 2048;
pub const ECL_PARTEND_VISIBLE: i32 = 4096;
pub const ECL_PENUMBBEG_VISIBLE: i32 = 8192;
pub const ECL_PENUMBEND_VISIBLE: i32 = 16384;
const NO_SAROS: f64 = -99_999_999.0;

fn cal(gregorian: bool) -> Calendar {
    if gregorian {
        Calendar::Gregorian
    } else {
        Calendar::Julian
    }
}

fn body(id: i32) -> PyResult<Body> {
    Body::from_id(id).map_err(to_py)
}

fn one_char(s: &str, what: &str) -> PyResult<char> {
    let mut it = s.chars();
    match (it.next(), it.next()) {
        (Some(c), None) => Ok(c),
        _ => Err(PyValueError::new_err(format!("{what} must be exactly one character"))),
    }
}

/// Planet-centric position (JD TT) of `target` seen from the centre of `center`:
/// (x[6]). Centre flags (HELCTR/BARYCTR/TOPOCTR) are ignored as by the original; HELCTR or
/// BARYCTR switch to an astrometric place (no aberration, no deflection), as it does.
#[pyfunction]
#[pyo3(signature = (engine, jd_tt, target, center, flags=256, sidereal=None))]
fn calc_pctr(
    py: Python<'_>,
    engine: &Engine,
    jd_tt: f64,
    target: i32,
    center: i32,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<[f64; 6]> {
    let mut f = flags & !(FLG_TOPOCTR | FLG_HELCTR | FLG_BARYCTR);
    if flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        f |= FLG_NOABERR | FLG_NOGDEFL;
    }
    let opt = options(f, sidereal, None)?;
    let (t, c) = (body(target)?, body(center)?);
    py.detach(|| engine.inner.calc_pctr(t, c, jd_tt, &opt)).map(|p| p.to_array()).map_err(to_py)
}

/// Next (previous) heliocentric longitude crossing, JD TT in and out.
#[pyfunction]
#[pyo3(signature = (engine, body_id, x2cross, jd_tt, flags=0, sidereal=None, backward=false))]
#[allow(clippy::too_many_arguments)]
fn helio_cross(
    py: Python<'_>,
    engine: &Engine,
    body_id: i32,
    x2cross: f64,
    jd_tt: f64,
    flags: i32,
    sidereal: Option<&str>,
    backward: bool,
) -> PyResult<f64> {
    let opt = options(flags & !(FLG_TOPOCTR | FLG_HELCTR | FLG_BARYCTR), sidereal, None)?;
    let b = body(body_id)?;
    py.detach(|| engine.inner.helio_cross(b, x2cross, jd_tt, &opt, backward)).map_err(to_py)
}

/// `helio_cross` with JD UT in and out.
#[pyfunction]
#[pyo3(signature = (engine, body_id, x2cross, jd_ut, flags=0, sidereal=None, backward=false))]
#[allow(clippy::too_many_arguments)]
fn helio_cross_ut(
    py: Python<'_>,
    engine: &Engine,
    body_id: i32,
    x2cross: f64,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
    backward: bool,
) -> PyResult<f64> {
    let opt = options(flags & !(FLG_TOPOCTR | FLG_HELCTR | FLG_BARYCTR), sidereal, None)?;
    let b = body(body_id)?;
    py.detach(|| engine.inner.helio_cross_ut(b, x2cross, jd_ut, &opt, backward)).map_err(to_py)
}

/// Next Moon crossing of the ecliptic: (JD TT, lon, lat).
#[pyfunction]
#[pyo3(signature = (engine, jd_tt, flags=0, sidereal=None))]
fn mooncross_node(
    py: Python<'_>,
    engine: &Engine,
    jd_tt: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<(f64, f64, f64)> {
    let opt = options(flags & !(FLG_TOPOCTR | FLG_HELCTR | FLG_BARYCTR), sidereal, None)?;
    py.detach(|| engine.inner.mooncross_node(jd_tt, &opt)).map_err(to_py)
}

/// `mooncross_node` with JD UT in and out: (JD UT, lon, lat).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, flags=0, sidereal=None))]
fn mooncross_node_ut(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    flags: i32,
    sidereal: Option<&str>,
) -> PyResult<(f64, f64, f64)> {
    let opt = options(flags & !(FLG_TOPOCTR | FLG_HELCTR | FLG_BARYCTR), sidereal, None)?;
    py.detach(|| engine.inner.mooncross_node_ut(jd_ut, &opt)).map_err(to_py)
}

/// Next/previous lunar eclipse visible at a place: (flags, tret[10], attr[20]) with
/// tret = max (or the visible instant nearest to it), 0, U1, U4, U2, U3, P1, P4, moonrise,
/// moonset (JD UT; 0 when absent or not visible) and attr as `lun_eclipse_how` at tret[0].
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, geopos, backward=false))]
fn lun_eclipse_when_loc(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    geopos: (f64, f64, f64),
    backward: bool,
) -> PyResult<(i32, [f64; 10], [f64; 20])> {
    use crate::eclipses::{ECL_MAX_VISIBLE, ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL, ECL_VISIBLE};
    use fe_core::eclipses::LunarKind;
    let obs = Observer { lon_deg: geopos.0, lat_deg: geopos.1, height_m: geopos.2 };
    let e = py
        .detach(|| lunar_local::lunar_eclipse_local_search(&engine.inner, jd_ut, &obs, backward))
        .map_err(to_py)?;
    let mut flags = match e.kind {
        LunarKind::Total => ECL_TOTAL,
        LunarKind::Partial => ECL_PARTIAL,
        LunarKind::Penumbral => ECL_PENUMBRAL,
    } | ECL_VISIBLE;
    for (t, bit) in [
        (e.max_visible.then_some(e.max), ECL_MAX_VISIBLE),
        (e.u1, ECL_PARTBEG_VISIBLE),
        (e.u2, ECL_TOTBEG_VISIBLE),
        (e.u3, ECL_TOTEND_VISIBLE),
        (e.u4, ECL_PARTEND_VISIBLE),
        (e.p1, ECL_PENUMBBEG_VISIBLE),
        (e.p4, ECL_PENUMBEND_VISIBLE),
    ] {
        if t.is_some() {
            flags |= bit;
        }
    }
    let z = |t: Option<f64>| t.unwrap_or(0.0);
    let tret =
        [e.max, 0.0, z(e.u1), z(e.u4), z(e.u2), z(e.u3), z(e.p1), z(e.p4), z(e.moonrise), z(e.moonset)];
    let a = &e.attrs;
    let (um, pm) = (a.shadow.umbral_mag.max(0.0), a.shadow.penumbral_mag.max(0.0));
    let mut attr = [0.0; 20];
    attr[..11].copy_from_slice(&[
        um,
        pm,
        0.0,
        0.0,
        (a.azimuth_deg + 180.0).rem_euclid(360.0),
        a.true_alt_deg,
        a.apparent_alt_deg,
        a.opposition_dist_deg,
        um,
        NO_SAROS,
        NO_SAROS,
    ]);
    Ok((flags, tret, attr))
}

/// UTC → (JD TT, JD UT1) with this engine's Delta T; ValueError for invalid date/time.
#[pyfunction]
#[pyo3(signature = (engine, year, month, day, hour, minute, second, gregorian=true))]
#[allow(clippy::too_many_arguments)]
fn utc_to_jd(
    engine: &Engine,
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    gregorian: bool,
) -> PyResult<(f64, f64)> {
    engine.inner.utc_to_jd(year, month, day, hour, minute, second, cal(gregorian)).map_err(to_py)
}

/// JD TT → UTC (y, m, d, h, min, sec); sec ≥ 60 inside a leap second.
#[pyfunction]
#[pyo3(signature = (engine, jd_tt, gregorian=true))]
fn jdet_to_utc(engine: &Engine, jd_tt: f64, gregorian: bool) -> (i32, i32, i32, i32, i32, f64) {
    engine.inner.jdet_to_utc(jd_tt, cal(gregorian))
}

/// JD UT1 → UTC.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, gregorian=true))]
fn jdut1_to_utc(engine: &Engine, jd_ut: f64, gregorian: bool) -> (i32, i32, i32, i32, i32, f64) {
    engine.inner.jdut1_to_utc(jd_ut, cal(gregorian))
}

/// Local time ↔ UTC: subtracts `offset` hours (east positive).
#[pyfunction]
fn utc_time_zone(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    offset: f64,
) -> (i32, i32, i32, i32, i32, f64) {
    utc::utc_time_zone(year, month, day, hour, minute, second, offset)
}

/// (valid, jd, (y, m, d, hour)).
#[pyfunction]
#[pyo3(signature = (year, month, day, hour=12.0, gregorian=true))]
fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    gregorian: bool,
) -> (bool, f64, (i32, i32, i32, f64)) {
    utc::date_conversion(year, month, day, hour, cal(gregorian))
}

/// Monday = 0 … Sunday = 6.
#[pyfunction]
fn day_of_week(jd: f64) -> i32 {
    utc::day_of_week(jd)
}

/// A new engine sharing this one's pack and frames, with a fixed Delta T (`user_days`,
/// days) and/or another lunar tidal acceleration (″/cy²). `None` keeps the current value.
#[pyfunction]
#[pyo3(signature = (engine, user_days=None, tid_acc=None))]
fn engine_with_delta_t(engine: &Engine, user_days: Option<f64>, tid_acc: Option<f64>) -> Engine {
    let mut m = engine.inner.delta_t_model();
    if let Some(d) = user_days {
        m.user_seconds = Some(d * 86400.0);
    }
    if let Some(a) = tid_acc {
        m.lunar_ndot = a;
    }
    Engine { inner: engine.inner.with_delta_t(m) }
}

/// Lunar tidal acceleration (″/cy²) of the engine's Delta T model.
#[pyfunction]
fn tid_acc(engine: &Engine) -> f64 {
    engine.inner.tidal_acceleration()
}

/// Sidereal time (hours) from a given obliquity and nutation (degrees).
#[pyfunction]
fn sidtime0(engine: &Engine, jd_ut: f64, eps: f64, nut: f64) -> f64 {
    engine.inner.sidtime0(jd_ut, eps, nut)
}

/// TDB coverage (start, end) of a pack body (NAIF id) relative to the barycentre.
#[pyfunction]
fn coverage(engine: &Engine, naif: i32) -> Option<(f64, f64)> {
    engine.inner.pack().chain_coverage(naif, 0)
}

/// (deg, min, sec, secfr, sign).
#[pyfunction]
fn split_deg(ddeg: f64, roundflag: i32) -> (i32, i32, i32, f64, i32) {
    let s = util::split_deg(ddeg, roundflag);
    (s.deg, s.min, s.sec, s.secfr, s.sign)
}

#[pyfunction]
fn cs2degstr(cs: i32) -> String {
    util::cs2degstr(cs)
}

#[pyfunction]
fn cs2lonlatstr(cs: i32, plus: &str, minus: &str) -> PyResult<String> {
    Ok(util::cs2lonlatstr(cs, one_char(plus, "plus")?, one_char(minus, "minus")?))
}

#[pyfunction]
#[pyo3(signature = (cs, sep, suppresszero=false))]
fn cs2timestr(cs: i32, sep: &str, suppresszero: bool) -> PyResult<String> {
    Ok(util::cs2timestr(cs, one_char(sep, "sep")?, suppresszero))
}

#[pyfunction]
fn csnorm(cs: i32) -> i32 {
    util::csnorm(cs)
}

#[pyfunction]
fn csroundsec(cs: i32) -> i32 {
    util::csroundsec(cs)
}

#[pyfunction]
fn d2l(x: f64) -> i32 {
    util::d2l(x)
}

#[pyfunction]
fn degnorm(x: f64) -> f64 {
    util::degnorm(x)
}

#[pyfunction]
fn radnorm(x: f64) -> f64 {
    util::radnorm(x)
}

#[pyfunction]
fn difdegn(p1: f64, p2: f64) -> f64 {
    util::difdegn(p1, p2)
}

#[pyfunction]
fn difdeg2n(p1: f64, p2: f64) -> f64 {
    util::difdeg2n(p1, p2)
}

#[pyfunction]
fn difrad2n(p1: f64, p2: f64) -> f64 {
    util::difrad2n(p1, p2)
}

#[pyfunction]
fn difcsn(p1: i32, p2: i32) -> i32 {
    util::difcsn(p1, p2)
}

#[pyfunction]
fn difcs2n(p1: i32, p2: i32) -> i32 {
    util::difcs2n(p1, p2)
}

#[pyfunction]
fn planet_name(id: i32) -> String {
    util::planet_name(id)
}

#[pyfunction]
fn ayanamsa_name(sid_mode: i32) -> &'static str {
    util::ayanamsa_name(sid_mode)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_pctr, m)?)?;
    m.add_function(wrap_pyfunction!(helio_cross, m)?)?;
    m.add_function(wrap_pyfunction!(helio_cross_ut, m)?)?;
    m.add_function(wrap_pyfunction!(mooncross_node, m)?)?;
    m.add_function(wrap_pyfunction!(mooncross_node_ut, m)?)?;
    m.add_function(wrap_pyfunction!(lun_eclipse_when_loc, m)?)?;
    m.add_function(wrap_pyfunction!(utc_to_jd, m)?)?;
    m.add_function(wrap_pyfunction!(jdet_to_utc, m)?)?;
    m.add_function(wrap_pyfunction!(jdut1_to_utc, m)?)?;
    m.add_function(wrap_pyfunction!(utc_time_zone, m)?)?;
    m.add_function(wrap_pyfunction!(date_conversion, m)?)?;
    m.add_function(wrap_pyfunction!(day_of_week, m)?)?;
    m.add_function(wrap_pyfunction!(engine_with_delta_t, m)?)?;
    m.add_function(wrap_pyfunction!(tid_acc, m)?)?;
    m.add_function(wrap_pyfunction!(sidtime0, m)?)?;
    m.add_function(wrap_pyfunction!(coverage, m)?)?;
    m.add_function(wrap_pyfunction!(split_deg, m)?)?;
    m.add_function(wrap_pyfunction!(cs2degstr, m)?)?;
    m.add_function(wrap_pyfunction!(cs2lonlatstr, m)?)?;
    m.add_function(wrap_pyfunction!(cs2timestr, m)?)?;
    m.add_function(wrap_pyfunction!(csnorm, m)?)?;
    m.add_function(wrap_pyfunction!(csroundsec, m)?)?;
    m.add_function(wrap_pyfunction!(d2l, m)?)?;
    m.add_function(wrap_pyfunction!(degnorm, m)?)?;
    m.add_function(wrap_pyfunction!(radnorm, m)?)?;
    m.add_function(wrap_pyfunction!(difdegn, m)?)?;
    m.add_function(wrap_pyfunction!(difdeg2n, m)?)?;
    m.add_function(wrap_pyfunction!(difrad2n, m)?)?;
    m.add_function(wrap_pyfunction!(difcsn, m)?)?;
    m.add_function(wrap_pyfunction!(difcs2n, m)?)?;
    m.add_function(wrap_pyfunction!(planet_name, m)?)?;
    m.add_function(wrap_pyfunction!(ayanamsa_name, m)?)?;
    m.add("LEAP_SECONDS_SOURCE", utc::LEAP_SECONDS_SOURCE)?;
    m.add("VERSION", env!("CARGO_PKG_VERSION"))?;
    for (k, v) in [
        ("ECL_PARTBEG_VISIBLE", ECL_PARTBEG_VISIBLE),
        ("ECL_TOTBEG_VISIBLE", ECL_TOTBEG_VISIBLE),
        ("ECL_TOTEND_VISIBLE", ECL_TOTEND_VISIBLE),
        ("ECL_PARTEND_VISIBLE", ECL_PARTEND_VISIBLE),
        ("ECL_PENUMBBEG_VISIBLE", ECL_PENUMBBEG_VISIBLE),
        ("ECL_PENUMBEND_VISIBLE", ECL_PENUMBEND_VISIBLE),
        ("SPLIT_DEG_ROUND_SEC", util::SPLIT_DEG_ROUND_SEC),
        ("SPLIT_DEG_ROUND_MIN", util::SPLIT_DEG_ROUND_MIN),
        ("SPLIT_DEG_ROUND_DEG", util::SPLIT_DEG_ROUND_DEG),
        ("SPLIT_DEG_ZODIACAL", util::SPLIT_DEG_ZODIACAL),
        ("SPLIT_DEG_KEEP_SIGN", util::SPLIT_DEG_KEEP_SIGN),
        ("SPLIT_DEG_KEEP_DEG", util::SPLIT_DEG_KEEP_DEG),
        ("SPLIT_DEG_NAKSHATRA", util::SPLIT_DEG_NAKSHATRA),
        ("ECL_NUT", fe_core::extras::ECL_NUT),
    ] {
        m.add(k, v)?;
    }
    Ok(())
}
