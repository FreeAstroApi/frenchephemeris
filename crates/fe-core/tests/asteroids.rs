// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Numbered asteroids beyond the minor pack: ids, dynamic chaining of any
//! number of small bodies, and the bundled pack `fepk-asteroids-v1.bin` vs
//! JPL SPK states recorded at fit time (tests/fixtures/asteroids/, written by
//! tools/build_asteroids.py). Pack tests skip when the packs are absent.

use fe_core::context::{Body, CalcOptions, Engine, AST_OFFSET};
use fe_core::pack::Pack;
use fe_core::Error;
use std::path::PathBuf;

fn find(name: &str) -> Option<PathBuf> {
    let mut dirs = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        dirs.push(PathBuf::from(d).join("packs"));
    }
    dirs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    dirs.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    dirs.into_iter().map(|d| d.join(name)).find(|p| p.exists())
}

#[test]
fn asteroid_ids() {
    let lilith = Body::from_id(AST_OFFSET + 1181).unwrap();
    assert_eq!(lilith, Body::Asteroid(1181));
    assert_eq!(lilith.id(), 11181);
    assert_eq!(lilith.naif(), Some(2001181));
    let sedna = Body::from_id(AST_OFFSET + 90377).unwrap();
    assert_eq!(sedna.naif(), Some(2090377));
    // the minor-pack bodies keep their own variants and ids
    assert_eq!(Body::from_id(10433).unwrap(), Body::Eros);
    assert_eq!(Body::from_id(10001).unwrap(), Body::Ceres);
    assert_eq!(Body::Eros.id(), 10433);
    assert_eq!(Body::Chiron.id(), 15);
    assert_eq!(Body::Cupido.id(), 40);
    assert!(Body::from_id(AST_OFFSET).is_err());
    // never-loaded asteroid: UnknownBody (NAIF id), not a wrong value
    if let Some(main) = find("fepk-de441-de442-v1.bin") {
        let eng = Engine::new(Pack::open(&main).unwrap());
        let r = eng.calc_ut(Body::Asteroid(123456), 2451545.0, &CalcOptions::default());
        assert!(matches!(r, Err(Error::UnknownBody(2123456))), "{r:?}");
    }
}

#[test]
fn bundled_asteroids_match_jpl() {
    let (Some(main), Some(ast)) = (find("fepk-de441-de442-v1.bin"), find("fepk-asteroids-v1.bin")) else {
        eprintln!("SKIP: main or asteroid pack not found");
        return;
    };
    let mut paths = vec![main.clone()];
    if let Some(minor) = find("fepk-minor-v2.bin") {
        paths.push(minor);
    }
    paths.push(ast);
    let pack = Pack::open_many(&paths).unwrap();
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/asteroids/asteroid_states.json");
    let fx: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let tol = fx["tolerance_km"].as_f64().unwrap() * 1.0001;
    let rows = fx["rows"].as_array().unwrap();
    let mut worst = 0.0f64;
    for r in rows {
        let naif = r[0].as_i64().unwrap() as i32;
        let jd = r[2].as_f64().unwrap();
        let want = [r[3].as_f64().unwrap(), r[4].as_f64().unwrap(), r[5].as_f64().unwrap()];
        let got = pack.position(naif, 10, jd).unwrap();
        let d = ((got[0] - want[0]).powi(2) + (got[1] - want[1]).powi(2) + (got[2] - want[2]).powi(2)).sqrt();
        assert!(d <= tol, "{naif} at {jd}: {d} km");
        worst = worst.max(d);
    }
    eprintln!("{} JPL states, max {worst:.4} km", rows.len());
    // supported range is explicit per body
    let eng = Engine::new(Pack::open_many(&paths).unwrap());
    let opt = CalcOptions::default();
    for (id, cov) in fx["coverage"].as_object().unwrap() {
        let b = Body::from_id(id.parse().unwrap()).unwrap();
        let (lo, hi) = (cov[0].as_f64().unwrap(), cov[1].as_f64().unwrap());
        assert!(eng.calc_tt(b, 0.5 * (lo + hi), &opt).is_ok(), "{id}");
        assert!(matches!(eng.calc_tt(b, lo - 1.0, &opt), Err(Error::OutOfRange { .. })), "{id} before {lo}");
    }
}
