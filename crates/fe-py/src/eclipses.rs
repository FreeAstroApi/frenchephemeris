// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Eclipse bindings: module-level functions taking the `Engine` (no
//! `multiple-pymethods` feature needed). Return tuples follow the classic
//! `swe_*eclipse*` layouts (times JD UT, flags with the `ECL_*` bit values),
//! so `swe_compat` only re-packs arguments. Everything runs without the GIL.

use super::{to_py, Engine};
use fe_core::context::Body;
use fe_core::eclipses::{self as ecl, LunarKind, SolarKind};
use fe_core::ephem::Observer;
use fe_core::occultations::{self as occ, Occulted, Search};
use fe_core::stars;
use pyo3::prelude::*;

pub const ECL_CENTRAL: i32 = 1;
pub const ECL_NONCENTRAL: i32 = 2;
pub const ECL_TOTAL: i32 = 4;
pub const ECL_ANNULAR: i32 = 8;
pub const ECL_PARTIAL: i32 = 16;
pub const ECL_ANNULAR_TOTAL: i32 = 32;
pub const ECL_PENUMBRAL: i32 = 64;
pub const ECL_VISIBLE: i32 = 128;
pub const ECL_MAX_VISIBLE: i32 = 256;
pub const ECL_1ST_VISIBLE: i32 = 512;
pub const ECL_2ND_VISIBLE: i32 = 1024;
pub const ECL_3RD_VISIBLE: i32 = 2048;
pub const ECL_4TH_VISIBLE: i32 = 4096;
pub const ECL_OCC_BEG_DAYLIGHT: i32 = 8192;
pub const ECL_OCC_END_DAYLIGHT: i32 = 16384;
pub const ECL_ONE_TRY: i32 = 32768;
/// Saros fields are not computed (the classic API's "not available" value).
const NO_SAROS: f64 = -99_999_999.0;

fn saros_pair(jd_tt: f64, solar: bool) -> (f64, f64) {
    ecl::saros(jd_tt, solar).map_or((NO_SAROS, NO_SAROS), |(s, m)| (s as f64, m as f64))
}

/// TT → UT with the engine's own Delta T model (not the library default).
fn ut(engine: &Engine, t: f64) -> f64 {
    engine.inner.tt_to_ut(t)
}

fn opt_ut(engine: &Engine, t: Option<f64>) -> f64 {
    t.map_or(0.0, |t| ut(engine, t))
}

fn observer(geopos: (f64, f64, f64)) -> Observer {
    Observer { lon_deg: geopos.0, lat_deg: geopos.1, height_m: geopos.2 }
}

/// Azimuth from north (east positive) → classic convention (from south, westward).
fn az_south(az_north: f64) -> f64 {
    (az_north + 180.0).rem_euclid(360.0)
}

fn lunar_flags(k: LunarKind) -> i32 {
    match k {
        LunarKind::Total => ECL_TOTAL,
        LunarKind::Partial => ECL_PARTIAL,
        LunarKind::Penumbral => ECL_PENUMBRAL,
    }
}

fn solar_type_flags(k: SolarKind) -> i32 {
    match k {
        SolarKind::Total => ECL_TOTAL,
        SolarKind::Annular => ECL_ANNULAR,
        SolarKind::Hybrid => ECL_ANNULAR_TOTAL,
        SolarKind::Partial => ECL_PARTIAL,
    }
}

/// Next/previous global lunar eclipse: (flags, tret[10]) with tret =
/// max, 0, U1, U4, U2, U3, P1, P4, 0, 0 (JD UT; 0 when a phase does not occur).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, ecltype=0, backward=false))]
pub fn lun_eclipse_when(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    ecltype: i32,
    backward: bool,
) -> PyResult<(i32, [f64; 10])> {
    let e = py
        .detach(|| {
            let mut t = engine.inner.ut_to_tt(jd_ut);
            loop {
                let e = ecl::lunar_eclipse_search(&engine.inner, t, backward)?;
                if ecltype == 0 || lunar_flags(e.kind) & ecltype != 0 {
                    return Ok(e);
                }
                t = e.max;
            }
        })
        .map_err(to_py)?;
    let tret = [
        ut(engine, e.max),
        0.0,
        opt_ut(engine, e.u1),
        opt_ut(engine, e.u4),
        opt_ut(engine, e.u2),
        opt_ut(engine, e.u3),
        ut(engine, e.p1),
        ut(engine, e.p4),
        0.0,
        0.0,
    ];
    Ok((lunar_flags(e.kind), tret))
}

/// Lunar-eclipse attributes at an instant and place: (flags, attr[20]) with
/// attr = umbral mag, penumbral mag, 0, 0, azimuth (from south), true alt,
/// apparent alt, distance from opposition, umbral mag, saros, saros member.
#[pyfunction]
pub fn lun_eclipse_how(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    geopos: (f64, f64, f64),
) -> PyResult<(i32, [f64; 20])> {
    let l = py
        .detach(|| ecl::lunar_local(&engine.inner, engine.inner.ut_to_tt(jd_ut), &observer(geopos)))
        .map_err(to_py)?;
    let (sar, mem) = if l.shadow.penumbral_mag > 0.0 {
        saros_pair(engine.inner.ut_to_tt(jd_ut), false)
    } else {
        (NO_SAROS, NO_SAROS)
    };
    let (um, pm) = (l.shadow.umbral_mag.max(0.0), l.shadow.penumbral_mag.max(0.0));
    let eclipse = pm > 0.0;
    // As the classic function: no flags while the Moon's centre is below the
    // apparent horizon (the magnitudes are still given); no eclipse at all gives
    // zero magnitudes and elongation.
    let flags = ecl::lunar_seen_kind(&l).map_or(0, |k| match k {
        LunarKind::Total => ECL_TOTAL,
        LunarKind::Partial => ECL_PARTIAL,
        LunarKind::Penumbral => ECL_PENUMBRAL,
    });
    let mut attr = [0.0; 20];
    attr[..11].copy_from_slice(&[
        um,
        pm,
        0.0,
        0.0,
        az_south(l.azimuth_deg),
        l.true_alt_deg,
        l.apparent_alt_deg,
        if eclipse { l.opposition_dist_deg } else { 0.0 },
        um,
        sar,
        mem,
    ]);
    Ok((flags, attr))
}

/// Next/previous global solar eclipse: (flags, tret[10]) with tret = max,
/// eclipse at local apparent noon (0 if never), P1, P4, U1, U4, C1, C2,
/// hybrid→total, hybrid→annular (JD UT).
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, ecltype=0, backward=false))]
pub fn sol_eclipse_when_glob(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    ecltype: i32,
    backward: bool,
) -> PyResult<(i32, [f64; 10])> {
    let flags_of = |e: &ecl::SolarEclipse| {
        solar_type_flags(e.kind) | if e.central { ECL_CENTRAL } else { ECL_NONCENTRAL }
    };
    let e = py
        .detach(|| {
            let mut t = engine.inner.ut_to_tt(jd_ut);
            loop {
                let e = ecl::solar_eclipse_search(&engine.inner, t, backward)?;
                let f = flags_of(&e);
                let types = ecltype & !(ECL_CENTRAL | ECL_NONCENTRAL);
                let geom = ecltype & (ECL_CENTRAL | ECL_NONCENTRAL);
                if (types == 0 || f & types != 0) && (geom == 0 || f & geom != 0) {
                    return Ok::<_, fe_core::Error>((e, ecl::solar_local_noon(&engine.inner, &e)?));
                }
                t = e.max;
            }
        })
        .map_err(to_py)?;
    let (e, noon) = e;
    let tret = [
        ut(engine, e.max),
        opt_ut(engine, noon),
        ut(engine, e.p1),
        ut(engine, e.p4),
        opt_ut(engine, e.u1),
        opt_ut(engine, e.u4),
        opt_ut(engine, e.c1),
        opt_ut(engine, e.c2),
        opt_ut(engine, e.h1),
        opt_ut(engine, e.h2),
    ];
    Ok((flags_of(&e), tret))
}

/// Classic attr layout. Fully covered discs report obscuration = (diameter
/// ratio)², which exceeds 1 (the classic convention; callers clamp).
fn solar_attr_array(a: &ecl::SolarAttrs, saros: (f64, f64)) -> [f64; 20] {
    let mut attr = [0.0; 20];
    let obsc = if a.in_umbra && a.ratio > 1.0 { a.ratio * a.ratio } else { a.obscuration };
    attr[..11].copy_from_slice(&[
        a.fraction_diameter,
        a.ratio,
        obsc,
        a.core_shadow_km,
        az_south(a.azimuth_deg),
        a.true_alt_deg,
        a.apparent_alt_deg,
        a.separation_deg,
        a.nasa_magnitude,
        saros.0,
        saros.1,
    ]);
    attr
}

fn geopos_of(w: &ecl::SolarWhere) -> [f64; 10] {
    let mut g = [0.0; 10];
    g[0] = w.lon_deg;
    g[1] = w.lat_deg;
    for (i, l) in w.limits.iter().enumerate() {
        if let Some((lon, lat)) = l {
            g[2 + 2 * i] = *lon;
            g[3 + 2 * i] = *lat;
        }
    }
    g
}

fn where_flags(w: &ecl::SolarWhere) -> i32 {
    let mut flags = local_flags(&w.attrs);
    if flags & (ECL_TOTAL | ECL_ANNULAR) != 0 {
        flags |= if w.central { ECL_CENTRAL } else { ECL_NONCENTRAL };
    }
    flags
}

/// Local type flags from the attributes at one instant.
fn local_flags(a: &ecl::SolarAttrs) -> i32 {
    if a.fraction_diameter <= 0.0 {
        0
    } else if a.in_umbra && (a.ratio >= 1.0 || a.ratio == 0.0) {
        ECL_TOTAL
    } else if a.in_umbra {
        ECL_ANNULAR
    } else {
        ECL_PARTIAL
    }
}

/// Where the eclipse is central (or greatest) at an instant: (flags,
/// geopos[10], attr[20]); geopos = lon, lat of the central line, then the
/// northern/southern limits of the umbra and of the penumbra (0 where a limit
/// does not touch the Earth).
#[pyfunction]
pub fn sol_eclipse_where(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
) -> PyResult<(i32, [f64; 10], [f64; 20])> {
    let w = py
        .detach(|| ecl::solar_where_ut(&engine.inner, engine.inner.ut_to_tt(jd_ut), jd_ut))
        .map_err(to_py)?;
    let flags = where_flags(&w);
    let sar = if flags != 0 { saros_pair(engine.inner.ut_to_tt(jd_ut), true) } else { (NO_SAROS, NO_SAROS) };
    Ok((flags, geopos_of(&w), solar_attr_array(&w.attrs, sar)))
}

/// Solar-eclipse attributes at an instant and place: (flags, attr[20]).
#[pyfunction]
pub fn sol_eclipse_how(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    geopos: (f64, f64, f64),
) -> PyResult<(i32, [f64; 20])> {
    let a = py
        .detach(|| ecl::solar_how(&engine.inner, engine.inner.ut_to_tt(jd_ut), &observer(geopos)))
        .map_err(to_py)?;
    // As the classic function: an eclipse is reported only while the Sun's centre
    // is above the apparent horizon; then VISIBLE plus CENTRAL / NONCENTRAL
    // (whether the shadow axis meets the Earth at that instant). Otherwise the
    // flags are 0 and the eclipse attributes are zero (direction fields kept).
    let mut flags = ecl::solar_seen_kind(&a).map_or(0, solar_type_flags);
    if flags == 0 {
        let mut attr = solar_attr_array(&a, (0.0, 0.0));
        for i in [0, 1, 2, 3, 8] {
            attr[i] = 0.0;
        }
        return Ok((0, attr));
    }
    flags |= ECL_VISIBLE | if a.axis_on_earth { ECL_CENTRAL } else { ECL_NONCENTRAL };
    Ok((flags, solar_attr_array(&a, saros_pair(engine.inner.ut_to_tt(jd_ut), true))))
}

/// Next/previous solar eclipse visible at a place: (flags, tret[10],
/// attr[20]) with tret = max, C1, C2, C3, C4, sunrise, sunset (JD UT; 0 if
/// absent). `attr` refers to the (visible) maximum.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, geopos, backward=false))]
pub fn sol_eclipse_when_loc(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    geopos: (f64, f64, f64),
    backward: bool,
) -> PyResult<(i32, [f64; 10], [f64; 20])> {
    let (l, g) = py
        .detach(|| {
            ecl::solar_local_search(&engine.inner, engine.inner.ut_to_tt(jd_ut), &observer(geopos), backward)
        })
        .map_err(to_py)?;
    let mut flags = solar_type_flags(l.kind);
    // The classic local search sets NONCENTRAL (never CENTRAL) when the shadow
    // axis misses the Earth at the reported maximum.
    if !l.attrs.axis_on_earth {
        flags |= ECL_NONCENTRAL;
    }
    for (on, bit) in [
        (l.visible, ECL_VISIBLE),
        (l.max_visible, ECL_MAX_VISIBLE),
        (l.c1_visible, ECL_1ST_VISIBLE),
        (l.c2_visible, ECL_2ND_VISIBLE),
        (l.c3_visible, ECL_3RD_VISIBLE),
        (l.c4_visible, ECL_4TH_VISIBLE),
    ] {
        if on {
            flags |= bit;
        }
    }
    let tret = [
        ut(engine, l.max),
        ut(engine, l.c1),
        opt_ut(engine, l.c2),
        opt_ut(engine, l.c3),
        ut(engine, l.c4),
        opt_ut(engine, l.sunrise),
        opt_ut(engine, l.sunset),
        0.0,
        0.0,
        0.0,
    ];
    Ok((flags, tret, solar_attr_array(&l.attrs, saros_pair(g.max, true))))
}

// --- occultations --------------------------------------------------------------------------

/// Planet id (0 = Sun, 2–9) or a star name.
fn occulted(body: &Bound<'_, PyAny>) -> PyResult<Occulted> {
    if let Ok(name) = body.extract::<String>() {
        return Ok(Occulted::Star(stars::lookup(&name).map_err(to_py)?));
    }
    let id: i32 = body.extract()?;
    Ok(Occulted::Body(Body::from_id(id).map_err(to_py)?))
}

/// Next/previous occultation of a planet (id) or star (name) globally:
/// (flags, tret[10]) like `sol_eclipse_when_glob`. With `ECL_ONE_TRY` in
/// `flags` only the next conjunction is tested; if it gives no occultation the
/// result is (0, tret) with tret[0] = time of that conjunction.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body, flags=0, ecltype=0, backward=false))]
pub fn lun_occult_when_glob(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: &Bound<'_, PyAny>,
    flags: i32,
    ecltype: i32,
    backward: bool,
) -> PyResult<(i32, [f64; 10])> {
    let obj = occulted(body)?;
    let one_try = flags & ECL_ONE_TRY != 0;
    let flags_of = |e: &ecl::SolarEclipse| {
        solar_type_flags(e.kind) | if e.central { ECL_CENTRAL } else { ECL_NONCENTRAL }
    };
    let r = py
        .detach(|| {
            let mut t = engine.inner.ut_to_tt(jd_ut);
            loop {
                match occ::occultation_search(&engine.inner, &obj, t, backward, one_try)? {
                    Search::Found(e) => {
                        let f = flags_of(&e);
                        let types = ecltype & !(ECL_CENTRAL | ECL_NONCENTRAL);
                        let geom = ecltype & (ECL_CENTRAL | ECL_NONCENTRAL);
                        if one_try || ((types == 0 || f & types != 0) && (geom == 0 || f & geom != 0)) {
                            return Ok::<_, fe_core::Error>(Search::Found((
                                e,
                                occ::local_noon(&engine.inner, &obj, &e)?,
                            )));
                        }
                        t = e.max;
                    }
                    Search::NoneAt(c) => return Ok(Search::NoneAt(c)),
                }
            }
        })
        .map_err(to_py)?;
    Ok(match r {
        Search::Found((e, noon)) => (
            flags_of(&e),
            [
                ut(engine, e.max),
                opt_ut(engine, noon),
                ut(engine, e.p1),
                ut(engine, e.p4),
                opt_ut(engine, e.u1),
                opt_ut(engine, e.u4),
                opt_ut(engine, e.c1),
                opt_ut(engine, e.c2),
                0.0,
                0.0,
            ],
        ),
        Search::NoneAt(c) => (0, [ut(engine, c), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    })
}

/// Where a body is centred on (or closest to the centre of) the lunar disc at
/// an instant: (flags, geopos[10], attr[20]) like `sol_eclipse_where`.
#[pyfunction]
pub fn lun_occult_where(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: &Bound<'_, PyAny>,
) -> PyResult<(i32, [f64; 10], [f64; 20])> {
    let obj = occulted(body)?;
    let w = py
        .detach(|| occ::occultation_where_ut(&engine.inner, &obj, engine.inner.ut_to_tt(jd_ut), jd_ut))
        .map_err(to_py)?;
    Ok((where_flags(&w), geopos_of(&w), solar_attr_array(&w.attrs, (NO_SAROS, NO_SAROS))))
}

/// Next/previous occultation visible at a place: (flags, tret[10],
/// attr[20]); tret = max, C1, C2, C3, C4, rise and set of the body between C1
/// and C4 (0 if none). Flags add `ECL_OCC_BEG_DAYLIGHT` / `ECL_OCC_END_DAYLIGHT`
/// when the Sun is up at C1 / C4. `ECL_ONE_TRY` as in `lun_occult_when_glob`.
#[pyfunction]
#[pyo3(signature = (engine, jd_ut, body, geopos, flags=0, backward=false))]
pub fn lun_occult_when_loc(
    py: Python<'_>,
    engine: &Engine,
    jd_ut: f64,
    body: &Bound<'_, PyAny>,
    geopos: (f64, f64, f64),
    flags: i32,
    backward: bool,
) -> PyResult<(i32, [f64; 10], [f64; 20])> {
    let obj = occulted(body)?;
    let one_try = flags & ECL_ONE_TRY != 0;
    let r = py
        .detach(|| {
            occ::occultation_local_search(
                &engine.inner,
                &obj,
                engine.inner.ut_to_tt(jd_ut),
                &observer(geopos),
                backward,
                one_try,
            )
        })
        .map_err(to_py)?;
    let (l, _) = match r {
        Search::Found(x) => x,
        Search::NoneAt(c) => {
            return Ok((0, [ut(engine, c), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0; 20]))
        }
    };
    let mut f = solar_type_flags(l.kind);
    for (on, bit) in [
        (l.visible, ECL_VISIBLE),
        (l.max_visible, ECL_MAX_VISIBLE),
        (l.c1_visible, ECL_1ST_VISIBLE),
        (l.c2_visible, ECL_2ND_VISIBLE),
        (l.c3_visible, ECL_3RD_VISIBLE),
        (l.c4_visible, ECL_4TH_VISIBLE),
        (l.c1_daylight, ECL_OCC_BEG_DAYLIGHT),
        (l.c4_daylight, ECL_OCC_END_DAYLIGHT),
    ] {
        if on {
            f |= bit;
        }
    }
    let tret = [
        ut(engine, l.max),
        ut(engine, l.c1),
        opt_ut(engine, l.c2),
        opt_ut(engine, l.c3),
        ut(engine, l.c4),
        opt_ut(engine, l.rise),
        opt_ut(engine, l.set),
        0.0,
        0.0,
        0.0,
    ];
    let mut attr = solar_attr_array(&l.attrs, (NO_SAROS, NO_SAROS));
    // Local occultation attributes are covered fractions, at most 1.
    attr[0] = attr[0].min(1.0);
    attr[2] = l.attrs.obscuration.min(1.0);
    Ok((f, tret, attr))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(lun_eclipse_when, m)?)?;
    m.add_function(wrap_pyfunction!(lun_eclipse_how, m)?)?;
    m.add_function(wrap_pyfunction!(sol_eclipse_when_glob, m)?)?;
    m.add_function(wrap_pyfunction!(sol_eclipse_where, m)?)?;
    m.add_function(wrap_pyfunction!(sol_eclipse_how, m)?)?;
    m.add_function(wrap_pyfunction!(sol_eclipse_when_loc, m)?)?;
    m.add_function(wrap_pyfunction!(lun_occult_when_glob, m)?)?;
    m.add_function(wrap_pyfunction!(lun_occult_when_loc, m)?)?;
    m.add_function(wrap_pyfunction!(lun_occult_where, m)?)?;
    for (k, v) in [
        ("ECL_CENTRAL", ECL_CENTRAL),
        ("ECL_NONCENTRAL", ECL_NONCENTRAL),
        ("ECL_TOTAL", ECL_TOTAL),
        ("ECL_ANNULAR", ECL_ANNULAR),
        ("ECL_PARTIAL", ECL_PARTIAL),
        ("ECL_ANNULAR_TOTAL", ECL_ANNULAR_TOTAL),
        ("ECL_HYBRID", ECL_ANNULAR_TOTAL),
        ("ECL_PENUMBRAL", ECL_PENUMBRAL),
        (
            "ECL_ALLTYPES_SOLAR",
            ECL_CENTRAL | ECL_NONCENTRAL | ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL,
        ),
        ("ECL_ALLTYPES_LUNAR", ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL),
        ("ECL_VISIBLE", ECL_VISIBLE),
        ("ECL_MAX_VISIBLE", ECL_MAX_VISIBLE),
        ("ECL_1ST_VISIBLE", ECL_1ST_VISIBLE),
        ("ECL_2ND_VISIBLE", ECL_2ND_VISIBLE),
        ("ECL_3RD_VISIBLE", ECL_3RD_VISIBLE),
        ("ECL_4TH_VISIBLE", ECL_4TH_VISIBLE),
        ("ECL_OCC_BEG_DAYLIGHT", ECL_OCC_BEG_DAYLIGHT),
        ("ECL_OCC_END_DAYLIGHT", ECL_OCC_END_DAYLIGHT),
        ("ECL_ONE_TRY", ECL_ONE_TRY),
    ] {
        m.add(k, v)?;
    }
    Ok(())
}
