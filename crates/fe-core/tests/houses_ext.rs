// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Additional house systems, speeds, house positions and Gauquelin sectors vs
//! the Swiss oracle (fixtures from tests/fixtures/houses/gen_houses_ext.py).
//! Run with `-- --nocapture` for the per-system max-delta tables quoted in
//! docs/GEOMETRY_REPORT.md.

#[macro_use]
mod common;

use fe_core::context::{Body, CalcOptions, Engine, EngineConfig};
use fe_core::ephem::{Corrections, Observer, Origin};
use fe_core::events::Atmosphere;
use fe_core::houses::*;
use fe_core::math::diff_deg;
use fe_core::Error;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn load(name: &str) -> Value {
    let p = common::oracle_path(&format!("houses/{name}"));
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// null (the oracle's NaN) -> NaN
fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

fn sys_of(v: &Value) -> char {
    v.as_str().unwrap().chars().next().unwrap()
}

/// Ecliptic (nearly) on the horizon: see tests/houses.rs.
fn near_degenerate(armc: f64, lat: f64, eps: f64) -> bool {
    let c = (armc.to_radians()).cos().abs();
    ((lat.abs() - (90.0 - eps)).abs() < 1e-6) && c < 1e-3
}

/// tan φ tan δ of the Sun (Sunshine), with the oracle's latitude clamp.
fn sun_t(lat: f64, dec: f64) -> f64 {
    lat.clamp(-90.0 + 1e-10, 90.0 - 1e-10).to_radians().tan() * dec.to_radians().tan()
}

fn class(lat: f64, eps: f64, err: bool) -> &'static str {
    if lat.abs() >= 90.0 {
        // exact poles: clamped to 90 − 1e-10°, tan φ ≈ 6e11 amplifies
        // 1e-12 differences of the Asc declination (B, U, Y, I, Q); reported only
        "pole"
    } else if err {
        "fallback"
    } else if lat.abs() >= 90.0 - eps {
        "polar"
    } else if lat.abs() > 60.0 {
        "60..polar"
    } else {
        "<=60"
    }
}

/// Cusps the oracle's 'i' (Makransky) computes differently from the
/// geometric Sunshine construction: when 2/3 of a semi-arc of the Sun
/// exceeds 90° (Makransky's published formulas are not available to us).
fn makransky_unreproduced(lat: f64, dec: f64) -> Vec<usize> {
    let t = sun_t(lat, dec);
    let dsa = 90.0 + t.clamp(-1.0, 1.0).asin().to_degrees();
    let mut v = vec![];
    if 2.0 * dsa / 3.0 > 90.0 {
        v.extend([7, 11]);
    }
    if 2.0 * (180.0 - dsa) / 3.0 > 90.0 {
        v.extend([1, 5]);
    }
    v
}

#[test]
fn extra_systems_match_swiss() {
    require_oracle!("houses/houses_ext.json");
    let rows = load("houses_ext.json")["rows"].as_array().unwrap().clone();
    let mut worst: BTreeMap<(char, &str), (f64, String)> = BTreeMap::new();
    let (mut n, mut skipped_i, mut nan_rows) = (0, 0, 0);
    for r in &rows {
        let r = r.as_array().unwrap();
        let sys = sys_of(&r[0]);
        let (armc, lat, eps, dec) = (f(&r[1]), f(&r[2]), f(&r[3]), f(&r[4]));
        let err = r[5].as_bool().unwrap();
        let exp: Vec<f64> = r[6..].iter().map(f).collect();
        if exp.iter().any(|x| x.is_nan()) {
            nan_rows += 1; // oracle NaN cusps (Sunshine, circumpolar Sun)
        }
        if sys == 'I' || sys == 'i' {
            // the oracle's Sunshine fails numerically in a ~1e-8 band around a
            // just-circumpolar Sun (ARMC dependent); not reproduced
            if (sun_t(lat, dec).abs() - 1.0).abs() < 1e-6 {
                continue;
            }
        }
        let res = house_cusps(sys, armc, lat, eps, Some(dec));
        if err {
            assert!(
                matches!(res, Err(Error::PolarHouses { .. })),
                "{sys} {armc} {lat} {eps} {dec}: expected PolarHouses, got {res:?}"
            );
        } else {
            assert!(res.is_ok(), "{sys} {armc} {lat} {eps} {dec}: unexpected {res:?}");
        }
        let (h, fb) = house_cusps_swiss_compat(sys, armc, lat, eps, Some(dec)).unwrap();
        assert_eq!(fb, err);
        if near_degenerate(armc, lat, eps) {
            continue;
        }
        // the oracle's 'i' has a singular branch for a Sun declination of
        // exactly 0 (cusps 2/8 and 5/11 change half-plane); any |δ| ≥ 1e-9 agrees
        // Horizon on the equator with the equinox at the zenith: every vertical
        // meets the ecliptic at the equinox (degenerate, measure zero)
        if sys == 'H' && lat == 0.0 && (armc == 0.0 || armc == 180.0) {
            continue;
        }
        if sys == 'i' && (dec == 0.0 || lat == 0.0) {
            continue; // (also exactly on the equator)
        }
        let skip = if sys == 'i' && !err { makransky_unreproduced(lat, dec) } else { vec![] };
        if !skip.is_empty() {
            skipped_i += 1;
        }
        let got: Vec<f64> = h.cusps.iter().chain(h.ascmc.iter()).copied().collect();
        assert_eq!(got.len(), exp.len(), "{sys}");
        let d = got
            .iter()
            .zip(&exp)
            .enumerate()
            // co-Asc (Munkasey) = asc(ARMC, 90° − φ) is degenerate for φ = 0
            // exactly (pole height 90°); the oracle returns 0 or 180 depending
            // on the house system
            .filter(|(i, _)| !(*i == got.len() - 2 && lat == 0.0))
            .filter(|(i, (_, b))| !skip.contains(i) && !b.is_nan())
            .map(|(_, (a, b))| diff_deg(*a, *b).abs())
            .fold(0.0, f64::max);
        let c = class(lat, eps, err);
        let e = worst.entry((sys, c)).or_insert((0.0, String::new()));
        if d > e.0 {
            *e = (d, format!("armc={armc} lat={lat} eps={eps} dec={dec}"));
        }
        n += 1;
    }
    println!("compared {n} rows ('i' rows with Makransky-specific cusps excluded from 2 cusps: {skipped_i}; rows with oracle NaN cusps: {nan_rows})");
    for ((s, c), (d, w)) in &worst {
        println!("{s} {c:10} max |Δ| = {d:.2e}°  ({w})");
    }
    for ((s, c), (d, w)) in &worst {
        let tol = match (*s, *c) {
            // Gauquelin = Placidus construction: the oracle's Placidus
            // iteration error (tests/houses.rs)
            ('G', "<=60") => 2e-7,
            ('G', "60..polar") => 5e-6,
            // Sunshine next to a circumpolar Sun: steep in the declination
            (_, "pole") => f64::INFINITY,
            ('I', _) | ('i', _) => 1e-6,
            _ => 1e-7,
        };
        assert!(*d <= tol, "{s} {c}: {d:e} > {tol:e} at {w}");
    }
}

#[test]
fn speeds_match_swiss() {
    require_oracle!("houses/house_speeds.json");
    let rows = load("house_speeds.json")["rows"].as_array().unwrap().clone();
    let mut vs_fd: BTreeMap<char, (f64, String)> = BTreeMap::new();
    let mut vs_oracle: BTreeMap<char, (f64, String)> = BTreeMap::new();
    let mut defects: BTreeMap<(char, usize), usize> = BTreeMap::new();
    let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1.0);
    for r in &rows {
        let r = r.as_array().unwrap();
        let sys = sys_of(&r[0]);
        let (armc, lat, eps, dec) = (f(&r[1]), f(&r[2]), f(&r[3]), f(&r[4]));
        let err = r[5].as_bool().unwrap();
        if (sys == 'I' || sys == 'i') && (sun_t(lat, dec).abs() - 1.0).abs() < 1e-3 {
            continue;
        }
        let spd: Vec<f64> = r[6].as_array().unwrap().iter().map(f).collect();
        let fd: Vec<f64> = r[7].as_array().unwrap().iter().map(f).collect();
        let (h, fb) = house_cusps_speeds_swiss_compat(sys, armc, lat, eps, Some(dec)).unwrap();
        assert_eq!(fb, err, "{sys} {armc} {lat}");
        let got: Vec<f64> = h.cusp_speeds.iter().chain(h.ascmc_speeds.iter()).copied().collect();
        let skip = if sys == 'i' && !err { makransky_unreproduced(lat, dec) } else { vec![] };
        for (i, ((g, e), d)) in got.iter().zip(&spd).zip(&fd).enumerate() {
            if skip.contains(&i) || e.is_nan() || d.is_nan() {
                continue;
            }
            let tag =
                format!("[{i}] armc={armc} lat={lat} eps={eps} dec={dec} got={g} oracle={e} oracle_fd={d}");
            // ours vs the oracle's own cusps differentiated (definition check);
            // a cusp jumping inside ±1e-3° (W/N boundaries) is not comparable
            if rel(*d, *g) > 1e-4 && (d.abs() > 1e4 || g.abs() < 1e-12) {
                continue;
            }
            let w = vs_fd.entry(sys).or_insert((0.0, String::new()));
            if rel(*g, *d) > w.0 {
                *w = (rel(*g, *d), tag.clone());
            }
            // oracle speed inconsistent with its own cusps: an oracle defect
            if rel(*e, *d) > 1e-5 {
                *defects.entry((sys, i)).or_default() += 1;
                continue;
            }
            let w = vs_oracle.entry(sys).or_insert((0.0, String::new()));
            if rel(*g, *e) > w.0 {
                *w = (rel(*g, *e), tag);
            }
        }
    }
    for (s, (d, w)) in &vs_fd {
        println!("{s} vs oracle finite diff: max rel {d:.2e} ({w})");
    }
    for (s, (d, w)) in &vs_oracle {
        println!("{s} vs oracle speed     : max rel {d:.2e} ({w})");
    }
    let mut by_sys: BTreeMap<char, Vec<String>> = BTreeMap::new();
    for ((s, i), n) in &defects {
        let name = if *i < 36 { format!("c{}:{n}", i + 1) } else { format!("a{}:{n}", i - 36) };
        by_sys.entry(*s).or_default().push(name);
    }
    for (s, v) in &by_sys {
        println!("oracle speed != its own finite difference: {s} {}", v.join(" "));
    }
    for (s, (d, w)) in vs_fd.iter().chain(vs_oracle.iter()) {
        assert!(*d <= 2e-5, "{s}: {d:e} at {w}");
    }
}

#[test]
fn house_pos_ext_matches_swiss() {
    require_oracle!("houses/house_pos_ext.json");
    let rows = load("house_pos_ext.json")["rows"].as_array().unwrap().clone();
    let mut worst: BTreeMap<char, (f64, String)> = BTreeMap::new();
    let mut quirks: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &rows {
        let r = r.as_array().unwrap();
        let sys = sys_of(&r[0]);
        let v: Vec<f64> = r[1..8].iter().map(f).collect();
        let (armc, lat, eps, lon, blat, dec, exp) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
        let got = house_position_ex(sys, armc, lat, eps, lon, blat, Some(dec)).unwrap();
        if exp.is_nan() {
            *quirks.entry("oracle NaN").or_default() += 1;
            continue;
        }
        let n = if sys == 'G' { 36.0 } else { 12.0 };
        assert!((1.0..n + 1.0).contains(&got), "{sys} {got}");
        // Sripati: the oracle returns 1.0 for every point whose Porphyry
        // position + 0.5 is ≥ 12 (Sripati house 12 and the first half of
        // house 1 collapse to 1.0); we return the geometric value.
        if sys == 'S' && exp == 1.0 && (got >= 12.0 || got < 1.5) {
            *quirks.entry("S: oracle 1.0 for Sripati house 12").or_default() += 1;
            continue;
        }
        // Horizon, southern hemisphere: the oracle's cusps are the mirror of
        // its own house_pos azimuth numbering, and it snaps points on a cusp
        // to that cusp's number; a point on a cusp is therefore inconsistent.
        if sys == 'H' && exp == exp.round() && (got - got.round()).abs() < 1e-6 {
            *quirks.entry("H: point on a cusp numbered inconsistently by the oracle").or_default() += 1;
            continue;
        }
        // the oracle's Krusinski position is ill-conditioned on the equator
        // (it moves 3e-5 house between φ = 0 and 1e-9; ours is continuous)
        if sys == 'U' && lat.abs() < 1e-6 {
            *quirks.entry("U: equator").or_default() += 1;
            continue;
        }
        // Sunshine with a circumpolar Sun: the oracle's positions degenerate
        if (sys == 'I' || sys == 'i') && sun_t(lat, dec).abs() >= 1.0 {
            *quirks.entry("I/i: circumpolar Sun").or_default() += 1;
            continue;
        }
        let d = diff_deg(got * 360.0 / n, exp * 360.0 / n).abs() * n / 12.0; // 30 per house or sector
        let w = worst.entry(sys).or_insert((0.0, String::new()));
        if d > w.0 {
            *w = (
                d,
                format!(
                    "armc={armc} lat={lat} eps={eps} lon={lon} blat={blat} dec={dec} got={got} exp={exp}"
                ),
            );
        }
    }
    for (s, (d, w)) in &worst {
        println!("{s} max |Δ| = {d:.2e} house-deg ({w})");
    }
    for (k, n) in &quirks {
        println!("oracle quirk skipped: {k}: {n}");
    }
    for (s, (d, w)) in &worst {
        // 1 mas offset/snap of the oracle (see tests/houses.rs), in longitude
        // for the linear systems (A V D: narrow houses near the poles)
        let tol = match s {
            'A' | 'V' | 'D' | 'L' | 'Q' | 'N' | 'S' => 3e-6,
            'Y' => 2e-6, // the oracle's APC position on the equator (φ = 0): 1.1e-6
            'G' => 1e-6, // the oracle's 1 mas Placidus house-degree offset, ×3 in sector units
            _ => 2.8e-7,
        };
        assert!(*d <= tol, "{s}: {d:e} > {tol:e} at {w}");
    }
}

#[test]
fn names() {
    assert_eq!(house_name('P'), "Placidus");
    assert_eq!(house_name('i'), "Sunshine/alt.");
    assert_eq!(house_name('I'), "Sunshine");
    assert_eq!(house_name('g'), "Gauquelin sectors");
    assert_eq!(house_name('Z'), "");
    for s in SYSTEMS {
        if s != 'J' {
            assert!(!house_name(s).is_empty(), "{s}");
        }
    }
}

fn pack_path() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FE_PACK") {
        return Some(PathBuf::from(p)).filter(|p| p.exists());
    }
    let name = "fepk-de441-de442-v1.bin";
    [
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs").join(name),
        PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join(name),
    ]
    .into_iter()
    .find(|p| p.exists())
}

#[test]
fn gauquelin_sector_and_dated_houses() {
    require_oracle!("houses/gauquelin.json");
    let Some(p) = pack_path() else {
        eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
        return;
    };
    let eng = Engine::with_config(fe_core::pack::Pack::open(&p).unwrap(), EngineConfig::SWISS_TIME);
    let fx = load("gauquelin.json");
    let mut worst: BTreeMap<i64, (f64, String)> = BTreeMap::new();
    for r in fx["sectors"].as_array().unwrap() {
        let r = r.as_array().unwrap();
        let t = f(&r[0]);
        let body = Body::from_id(r[1].as_i64().unwrap() as i32).unwrap();
        let meth = r[2].as_i64().unwrap();
        let g: Vec<f64> = r[3].as_array().unwrap().iter().map(f).collect();
        let obs = Observer { lon_deg: g[0], lat_deg: g[1], height_m: g[2] };
        let atm = Atmosphere { pressure_hpa: f(&r[4]), temperature_c: f(&r[5]) };
        let topo = r[6].as_bool().unwrap();
        let opt = CalcOptions {
            origin: if topo { Origin::Topocentric(obs) } else { Origin::Geocentric },
            corrections: Corrections::APPARENT,
            ..Default::default()
        };
        let res =
            eng.gauquelin_sector(t, body, GauquelinMethod::from_code(meth as i32).unwrap(), &obs, atm, &opt);
        match (&r[7], res) {
            (Value::Null, Err(_)) => {}
            (Value::Null, Ok(s)) => panic!("expected error, got {s} for {r:?}"),
            (e, Err(err)) => panic!("expected {e}, got {err} for {r:?}"),
            (e, Ok(s)) => {
                let d = diff_deg(s * 10.0, f(e) * 10.0).abs() / 10.0;
                let w = worst.entry(meth).or_insert((0.0, String::new()));
                if d > w.0 {
                    *w = (d, format!("{r:?} got {s}"));
                }
            }
        }
    }
    for (m, (d, w)) in &worst {
        println!("gauquelin method {m}: max |Δ| = {d:.2e} sector ({w})");
    }
    for (m, (d, w)) in &worst {
        // methods 0/1: planet places (JPL vs the oracle's files, ~1e-6 deg);
        // 2–5: rise/set instants (seconds) over a ~0.5-day arc
        let tol = if *m <= 1 { 1e-5 } else { 2e-4 };
        assert!(*d <= tol, "method {m}: {d:e} at {w}");
    }
    let mut wd: f64 = 0.0;
    for r in fx["dated_houses"].as_array().unwrap() {
        let r = r.as_array().unwrap();
        let (t, lat, lon, s) = (f(&r[0]), f(&r[1]), f(&r[2]), sys_of(&r[3]));
        let res = eng.house_cusps_ut(t, lat, lon, s, None, false);
        if r[4].is_null() {
            assert!(res.is_err(), "{s} {t} {lat}");
            continue;
        }
        let (h, _) = res.unwrap();
        let exp: Vec<f64> = r[4].as_array().unwrap().iter().chain(r[5].as_array().unwrap()).map(f).collect();
        let got: Vec<f64> = h.houses.cusps.iter().chain(h.houses.ascmc.iter()).copied().collect();
        let skip =
            if s == 'i' { makransky_unreproduced(lat, eng.sun_declination_ut(t).unwrap()) } else { vec![] };
        for (i, (a, b)) in got.iter().zip(&exp).enumerate() {
            if skip.contains(&i) {
                continue;
            }
            let d = diff_deg(*a, *b).abs();
            if d > 1e-5 {
                println!("dated {s} t={t} lat={lat} lon={lon}: got {a} exp {b}");
            }
            wd = wd.max(d);
        }
        if s == 'I' {
            let hh = eng.houses_ut(t, lat, lon, 'I', None).unwrap();
            assert_eq!(hh.cusps.to_vec(), h.houses.cusps);
        }
    }
    println!("dated houses (I i G P Y H): max |Δ| = {wd:.2e}°");
    assert!(wd < 1e-4, "{wd}");
}

/// Sidereal houses, traditional and projection modes (`SIDBIT_ECL_T0`,
/// `SIDBIT_SSY_PLANE`, J2000/J1900/B1950), vs `houses_ex2(FLG_SIDEREAL)`.
#[test]
fn sidereal_projection_houses() {
    require_oracle!("houses/houses_sidereal.json");
    use fe_core::ayanamsa::{Ayanamsa, SidBits};
    let Some(p) = pack_path() else {
        eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
        return;
    };
    let eng = Engine::with_config(fe_core::pack::Pack::open(&p).unwrap(), EngineConfig::SWISS_TIME);
    let fx = load("houses_sidereal.json");
    let mut worst: BTreeMap<i64, (f64, String)> = BTreeMap::new();
    for r in fx["rows"].as_array().unwrap() {
        let r = r.as_array().unwrap();
        let mode = r[0].as_i64().unwrap();
        let (t, lat, lon, s) = (f(&r[1]), f(&r[2]), f(&r[3]), sys_of(&r[4]));
        let m = Ayanamsa::from_swiss((mode & 0xff) as i32).unwrap();
        let bits = SidBits((mode as u32) & !0xff);
        let res = eng.house_cusps_ut_bits(t, lat, lon, s, Some((m, bits)), false, false);
        if r[5].is_null() {
            assert!(res.is_err(), "{mode} {s} {t} {lat}");
            continue;
        }
        let (h, _) = res.unwrap();
        let exp: Vec<f64> = r[5].as_array().unwrap().iter().chain(r[6].as_array().unwrap()).map(f).collect();
        let got: Vec<f64> = h.houses.cusps.iter().chain(h.houses.ascmc.iter()).copied().collect();
        let d = got.iter().zip(&exp).map(|(a, b)| diff_deg(*a, *b).abs()).fold(0.0, f64::max);
        if std::env::var("FE_DEBUG_SID").is_ok() && mode > 255 {
            let dv: Vec<String> =
                got.iter().zip(&exp).map(|(a, b)| format!("{:.1e}", diff_deg(*a, *b))).collect();
            println!("{mode} {s} t={t:.1} lat={lat:.1}: {}", dv.join(" "));
        }
        let w = worst.entry(mode).or_insert((0.0, String::new()));
        if d > w.0 {
            *w = (d, format!("{s} t={t} lat={lat} lon={lon}"));
        }
    }
    for (m, (d, w)) in &worst {
        println!("sidereal mode {m}: max |Δ| = {d:.2e}° ({w})");
    }
    for (m, (d, w)) in &worst {
        // SSY_PLANE: the invariable plane's orientation differs from the
        // oracle's by ~1e-5° (0.042″ in body latitudes, ayanamsa report), which
        // moves Ω and is amplified in high-latitude cusps
        let tol = if m & 512 != 0 { 2e-4 } else { 1e-5 };
        assert!(*d < tol, "mode {m}: {d:e} at {w}");
    }
}
