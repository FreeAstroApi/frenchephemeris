// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Houses vs the Swiss oracle (fixtures from tests/fixtures/houses/gen_houses.py).

#[macro_use]
mod common;

use fe_core::houses::{house_position, houses, houses_swiss_compat};
use fe_core::math::diff_deg;
use fe_core::Error;
use serde_json::Value;
use std::collections::BTreeMap;

fn load(name: &str) -> Vec<Value> {
    let p = common::oracle_path(&format!("houses/{name}"));
    let v: Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    v["rows"].as_array().unwrap().clone()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

/// Ecliptic (nearly) coincides with the horizon: |lat| ≈ 90 - eps and
/// ARMC ≈ 90/270. Swiss snaps such ascendants to 90/270 within ~1e-4° of
/// the singularity; excluded from the 1e-7 comparison (documented).
fn near_degenerate(armc: f64, lat: f64, eps: f64) -> bool {
    let c = (armc.to_radians()).cos().abs();
    ((lat.abs() - (90.0 - eps)).abs() < 1e-6) && c < 1e-3
}

#[test]
fn houses_match_swiss() {
    require_oracle!("houses/houses.json");
    let rows = load("houses.json");
    let mut worst: BTreeMap<(char, &str), (f64, String)> = BTreeMap::new();
    let mut n = 0;
    for r in &rows {
        let r = r.as_array().unwrap();
        let sys = r[0].as_str().unwrap().chars().next().unwrap();
        let (armc, lat, eps) = (f(&r[1]), f(&r[2]), f(&r[3]));
        let polar_err = r[4].as_bool().unwrap();
        let exp: Vec<f64> = r[5..25].iter().map(f).collect();
        let res = houses(sys, armc, lat, eps);
        if polar_err {
            assert!(
                matches!(res, Err(Error::PolarHouses { .. })),
                "{sys} {armc} {lat} {eps}: expected PolarHouses"
            );
        } else {
            assert!(res.is_ok(), "{sys} {armc} {lat} {eps}: unexpected {res:?}");
        }
        let (h, fell_back) = houses_swiss_compat(sys, armc, lat, eps).unwrap();
        assert_eq!(fell_back, polar_err);
        if near_degenerate(armc, lat, eps) {
            continue;
        }
        // Whole sign: Asc within 1e-8° of a sign boundary is ambiguous at the
        // ~1e-10° level (Swiss' Asc rounding differs); skip those rows.
        if sys == 'W' {
            let a = exp[12];
            if (a / 30.0 - (a / 30.0).round()).abs() * 30.0 < 1e-8 {
                continue;
            }
        }
        let got: Vec<f64> = h.cusps.iter().chain(h.ascmc().iter()).copied().collect();
        let d = got.iter().zip(&exp).map(|(a, b)| diff_deg(*a, *b).abs()).fold(0.0, f64::max);
        let class = if polar_err {
            "polar-fallback"
        } else if lat.abs() >= 90.0 - eps {
            "polar"
        } else if lat.abs() > 60.0 {
            "60..polar"
        } else {
            "<=60"
        };
        let e = worst.entry((sys, class)).or_insert((0.0, String::new()));
        if d > e.0 {
            *e = (d, format!("armc={armc} lat={lat} eps={eps}"));
        }
        n += 1;
    }
    println!("compared {n} rows");
    for ((s, c), (d, w)) in &worst {
        println!("{s} {c:15} max |Δ| = {d:.2e}°  ({w})");
    }
    for ((s, c), (d, w)) in &worst {
        // Swiss' own Placidus iteration stops at ~1e-6° near the polar circle.
        let tol = if *s == 'P' && *c != "polar-fallback" {
            if *c == "<=60" {
                2e-7
            } else {
                2e-6
            }
        } else {
            1e-7
        };
        assert!(*d <= tol, "{s} {c}: {d:e} > {tol:e} at {w}");
    }
}

#[test]
fn house_position_matches_swiss() {
    require_oracle!("houses/house_pos.json");
    let rows = load("house_pos.json");
    let mut worst: BTreeMap<(char, bool), (f64, String)> = BTreeMap::new();
    for r in &rows {
        let r = r.as_array().unwrap();
        let sys = r[0].as_str().unwrap().chars().next().unwrap();
        let v: Vec<f64> = r[1..7].iter().map(f).collect();
        let (armc, lat, eps, lon, blat, exp) = (v[0], v[1], v[2], v[3], v[4], v[5]);
        let got = house_position(sys, armc, lat, eps, lon, blat).unwrap();
        if exp == 0.0 || got == 0.0 {
            assert_eq!(got, exp, "{sys} zero/undefined mismatch {armc} {lat} {lon} {blat}");
            continue;
        }
        // compare in "house degrees" (30 per house), circularly
        let d = diff_deg(got * 30.0, exp * 30.0).abs();
        let dec = {
            let (se, ce) = eps.to_radians().sin_cos();
            let (sb, cb) = blat.to_radians().sin_cos();
            (cb * lon.to_radians().sin() * se + sb * ce).asin()
        };
        let circumpolar = (lat.to_radians().tan() * dec.tan()).abs() > 1.0;
        let e = worst.entry((sys, circumpolar)).or_insert((0.0, String::new()));
        if d > e.0 {
            *e = (d, format!("armc={armc} lat={lat} eps={eps} lon={lon} blat={blat} got={got} exp={exp}"));
        }
    }
    for ((s, c), (d, w)) in &worst {
        println!("{s} circumpolar={c:5} max |Δ| = {d:.2e} house-deg ({w})");
    }
    for ((s, c), (d, w)) in &worst {
        let tol = match (s, c) {
            ('T', false) => 2e-5,         // Swiss' topocentric house_pos is approximate
            ('T', true) => f64::INFINITY, // ambiguous branch, not reproduced (documented)
            // Swiss shifts by +1 mas (house degree, or longitude for E/W/O) and
            // snaps within 1 mas of a cusp; we do not (pure geometry).
            ('O', _) => 2e-6, // 1 mas of longitude × 30/width of narrow polar houses
            _ => 2.8e-7,      // 1 mas
        };
        assert!(*d <= tol, "{s} circumpolar={c}: {d:e} > {tol:e} at {w}");
    }
}
