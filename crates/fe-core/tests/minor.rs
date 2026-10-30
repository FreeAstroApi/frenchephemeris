// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Minor bodies: ids, the fitted minor-body pack vs the JPL source SPKs,
//! supported ranges, and apparent places vs JPL Horizons.
//!
//! Fixtures: tests/fixtures/minor/gen_minor_fixtures.py (JPL source states)
//! and minor_horizons.json (written by harness/check_minor.py from cached
//! Horizons OBSERVER responses). Packs are gitignored build products, looked
//! up in `$FE_DATA_DIR/packs`, `<repo>/data/packs` and the main checkout's
//! data dir; tests needing them skip with a message when absent.

use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::ephem::Corrections;
use fe_core::pack::{Pack, MINOR_BODIES};
use fe_core::Error;
use std::path::PathBuf;

const MAIN: &str = "fepk-de441-de442-v1.bin";
const MINOR: &str = "fepk-minor-v2.bin";

fn find(name: &str) -> Option<PathBuf> {
    let mut dirs = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        dirs.push(PathBuf::from(d).join("packs"));
    }
    dirs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    dirs.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    dirs.into_iter().map(|d| d.join(name)).find(|p| p.exists())
}

fn packs() -> Option<(PathBuf, PathBuf)> {
    match (find(MAIN), find(MINOR)) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => {
            eprintln!(
                "SKIP: {MAIN} and/or {MINOR} not found (build with tools/build_pack.py, tools/fit_minor.py)"
            );
            None
        }
    }
}

fn fixture(name: &str) -> serde_json::Value {
    let path = format!("{}/tests/fixtures/minor/{name}", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(path).expect("fixture")).unwrap()
}

const ALL: [(Body, i32, i32); 8] = [
    (Body::Chiron, 15, 2002060),
    (Body::Pholus, 16, 2005145),
    (Body::Ceres, 17, 2000001),
    (Body::Pallas, 18, 2000002),
    (Body::Juno, 19, 2000003),
    (Body::Vesta, 20, 2000004),
    (Body::Eris, 146199, 2136199),
    (Body::Eros, 10433, 2000433),
];

#[test]
fn ids_follow_the_api_convention() {
    for (b, id, naif) in ALL {
        assert_eq!(Body::from_id(id).unwrap(), b);
        assert_eq!(b.id(), id);
        assert_eq!(b.naif(), Some(naif));
        assert!(MINOR_BODIES.contains(&naif));
    }
    // AST_OFFSET + MPC number aliases
    assert_eq!(Body::from_id(10001).unwrap(), Body::Ceres);
    assert_eq!(Body::from_id(12060).unwrap(), Body::Chiron);
    assert_eq!(Body::from_id(15145).unwrap(), Body::Pholus);
    // 5145 (an id the API's chart package passes for Pholus) is not a body
    assert!(matches!(Body::from_id(5145), Err(Error::UnknownBody(5145))));
    // 21/22 are the natural apogee/perigee (Swiss numbering); 23 is no body
    assert!(matches!(Body::from_id(23), Err(Error::UnknownBody(23))));
}

#[test]
fn fitted_pack_matches_sources() {
    let Some((main, minor)) = packs() else { return };
    let pack = Pack::open_many(&[main, minor]).unwrap();
    let fx = fixture("minor_states.json");
    let mut worst = 0.0f64;
    for r in fx["rows"].as_array().unwrap() {
        let naif = r[0].as_i64().unwrap() as i32;
        let jd = f64::from_bits(r[6].as_u64().unwrap());
        let want = [r[3].as_f64().unwrap(), r[4].as_f64().unwrap(), r[5].as_f64().unwrap()];
        let got = pack.position(naif, 10, jd).unwrap();
        let d = ((got[0] - want[0]).powi(2) + (got[1] - want[1]).powi(2) + (got[2] - want[2]).powi(2)).sqrt();
        assert!(d < 0.1, "{} at JD {jd}: {d} km from its source ({})", r[1], r[7]);
        worst = worst.max(d);
    }
    eprintln!("max |pack - JPL source| = {worst:.4} km over {} rows", fx["rows"].as_array().unwrap().len());
}

#[test]
fn supported_ranges_are_explicit() {
    let Some((main, minor)) = packs() else { return };
    let pack = Pack::open_many(&[main.clone(), minor.clone()]).unwrap();
    let fx = fixture("minor_states.json");
    for (name, c) in fx["coverage"].as_object().unwrap() {
        let naif = c["naif"].as_i64().unwrap() as i32;
        let (a, b) = (c["jd_tdb"][0].as_f64().unwrap(), c["jd_tdb"][1].as_f64().unwrap());
        let (lo, hi) = pack.chain_coverage(naif, 0).unwrap();
        assert_eq!((lo, hi), (a, b), "{name} coverage");
        assert!(pack.position(naif, 0, a).is_ok() && pack.position(naif, 0, b).is_ok());
        // geocentric: the intersection with the main pack's Earth (starts 1720628.5)
        let (glo, ghi) = pack.chain_coverage(naif, 399).unwrap();
        assert_eq!((glo, ghi), (a.max(1720628.5), b));
        for t in [a - 1e-6, b + 1e-6, a - 365.25 * 100.0] {
            assert!(
                matches!(pack.position(naif, 0, t), Err(Error::OutOfRange { body, .. }) if body == naif),
                "{name} at {t} must be OutOfRange"
            );
        }
    }
    // Engine level, UT. JPL data for Chiron/Pholus/Eris/Eros start 1599-12-11; before
    // that our integration extends them to the Swiss oracle's first dates: Chiron
    // 675-01-01, Eris/Eros 1500-01-01 (Julian), Pholus the main pack start.
    let eng = Engine::new(Pack::open_many(&[main.clone(), minor]).unwrap());
    let opt = CalcOptions::default();
    let jd_0001 = 1721423.5;
    let jd_1550 = 2287184.5;
    for (b, first_ok, first_bad) in [
        (Body::Chiron, 1967602.5, 1967600.5),
        (Body::Eris, 2268933.5, 2268931.5),
        (Body::Eros, 2268933.5, 2268931.5),
        (Body::Pholus, jd_0001, 1720620.0),
    ] {
        assert!(eng.calc_ut(b, first_ok, &opt).is_ok(), "{b:?} at {first_ok}");
        assert!(
            matches!(eng.calc_ut(b, first_bad, &opt), Err(Error::OutOfRange { .. })),
            "{b:?} at {first_bad}"
        );
        assert!(eng.calc_ut(b, jd_1550, &opt).is_ok());
        assert!(eng.calc_ut(b, 2305426.5, &opt).is_ok()); // the JPL joint
        assert!(eng.calc_ut(b, 2451545.0, &opt).is_ok());
        assert!(eng.calc_ut(b, 2597640.5, &opt).is_ok()); // 2400-01-01
    }
    for b in [Body::Ceres, Body::Pallas, Body::Juno, Body::Vesta] {
        assert!(eng.calc_ut(b, jd_0001, &opt).is_ok(), "{b:?} in year 1");
        assert!(eng.calc_ut(b, 2597640.5, &opt).is_ok());
    }
    let jd_600 = 1940299.5;
    // Same call twice: same error (no state).
    assert_eq!(eng.calc_ut(Body::Chiron, jd_600, &opt), eng.calc_ut(Body::Chiron, jd_600, &opt));
    assert!(matches!(eng.calc_ut(Body::Chiron, jd_600, &opt), Err(Error::OutOfRange { .. })));
    // Without the minor pack the bodies are unknown, not silently wrong.
    let bare = Engine::new(Pack::open(&main).unwrap());
    assert!(matches!(bare.calc_ut(Body::Ceres, 2451545.0, &opt), Err(Error::UnknownBody(2000001))));
}

#[test]
fn integrated_extension_joins_jpl_continuously() {
    let Some((main, minor)) = packs() else { return };
    let pack = Pack::open_many(&[main, minor]).unwrap();
    let t: f64 = 2305426.5; // first JPL epoch = end of our integration
    for naif in [2002060, 2005145, 2136199, 2000433] {
        let before = f64::from_bits(t.to_bits() - 1); // previous double: integrated side
        let (p0, v0) = pack.state(naif, 0, before).unwrap();
        let (p1, v1) = pack.state(naif, 0, t).unwrap(); // JPL side (later series wins)
        let dp = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2) + (p1[2] - p0[2]).powi(2)).sqrt();
        let dv = ((v1[0] - v0[0]).powi(2) + (v1[1] - v0[1]).powi(2) + (v1[2] - v0[2]).powi(2)).sqrt();
        assert!(dp < 0.05 && dv < 0.5, "{naif}: jump {dp} km, {dv} km/day at the JPL joint");
    }
}

#[test]
fn astrometric_places_match_horizons() {
    let Some((main, minor)) = packs() else { return };
    let fx = fixture("minor_horizons.json");
    let eng = Engine::new(Pack::open_many(&[main, minor]).unwrap());
    // Horizons QUANTITIES 1: light-time corrected, no aberration/deflection, ICRF.
    let opt = CalcOptions {
        corrections: Corrections { light_time: true, deflection: false, aberration: false },
        equatorial: true,
        j2000: true,
        icrs: true,
        speed: false,
        ..Default::default()
    };
    let lim = fx["limit_arcsec"].as_object().unwrap();
    for r in fx["rows"].as_array().unwrap() {
        let name = r["body"].as_str().unwrap();
        let b = Body::from_id(r["id"].as_i64().unwrap() as i32).unwrap();
        let jd_tt = r["jd_tt"].as_f64().unwrap();
        let p = eng.calc_tt(b, jd_tt, &opt).unwrap();
        let (ra, de) = (r["ra"].as_f64().unwrap().to_radians(), r["dec"].as_f64().unwrap().to_radians());
        let (ra1, de1) = (p.lon.to_radians(), p.lat.to_radians());
        let h = ((de1 - de) / 2.0).sin().powi(2) + de.cos() * de1.cos() * ((ra1 - ra) / 2.0).sin().powi(2);
        let sep = 2.0 * h.sqrt().asin().to_degrees() * 3600.0;
        let tol = lim[name].as_f64().unwrap();
        assert!(sep < tol, "{name} JD_TT {jd_tt}: {sep:.4}\" from Horizons (tol {tol}\")");
    }
}
