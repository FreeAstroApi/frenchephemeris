// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Natural apsides (21/22) and fictitious bodies (40+) vs the Swiss 2.10
//! black-box oracle (fixtures: tests/fixtures/fictitious/gen_fictitious.py).
//! Skips when no pack is present. Run with `--nocapture` for the tables.

#[macro_use]
mod common;

use fe_core::context::{Body, CalcOptions, Engine, EngineConfig};
use fe_core::ephem::{Corrections, Origin};
use fe_core::math::diff_deg;
use fe_core::pack::Pack;
use serde_json::Value;
use std::path::PathBuf;

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

fn engine() -> Option<Engine> {
    match pack_path() {
        Some(p) => Some(Engine::with_config(Pack::open(&p).expect("pack"), EngineConfig::SWISS_DELTA_T)),
        None => {
            eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
            None
        }
    }
}

fn rows() -> Vec<Value> {
    let p = common::oracle_path("fictitious/fictitious.json");
    let v: Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    v["rows"].as_array().unwrap().clone()
}

fn era(jd: f64) -> usize {
    let y = 2000.0 + (jd - 2451545.0) / 365.25;
    if y < 1000.0 {
        0
    } else if y < 1800.0 {
        1
    } else if y <= 2100.0 {
        2
    } else {
        3
    }
}
const ERAS: [&str; 4] = ["0–1000", "1000–1800", "1800–2100", "2100–2400"];

/// max |Δ| per era of [lon °, lat °, dist rel, lon speed °/d, lat speed °/d, dist speed rel/d]
type Stats = [[f64; 6]; 4];

fn accumulate(st: &mut Stats, k: usize, ours: &[f64; 6], sw: &[f64]) {
    let d = [
        diff_deg(ours[0], sw[0]).abs(),
        (ours[1] - sw[1]).abs(),
        (ours[2] / sw[2] - 1.0).abs(),
        (ours[3] - sw[3]).abs(),
        (ours[4] - sw[4]).abs(),
        ((ours[5] - sw[5]) / sw[2]).abs(),
    ];
    for j in 0..6 {
        st[k][j] = st[k][j].max(d[j]);
    }
}

fn print(name: &str, st: &Stats) {
    for k in 0..4 {
        println!(
            "{name:12} {:10} lon {:9.5}° lat {:9.5}° dist {:8.2e}  vlon {:8.2e} vlat {:8.2e} vdist {:8.2e}",
            ERAS[k], st[k][0], st[k][1], st[k][2], st[k][3], st[k][4], st[k][5]
        );
    }
}

fn arr(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

#[test]
fn natural_apsides_vs_swiss() {
    require_oracle!("fictitious/fictitious.json");
    let Some(e) = engine() else { return };
    let opt = CalcOptions::default();
    // max |Δlon| (deg) per era, measured and rounded up: [apogee, perigee]
    let tol_lon: [[f64; 2]; 4] = [[0.016, 0.024], [0.020, 0.018], [0.016, 0.020], [0.013, 0.018]];
    for (id, body) in [(21, Body::IntpApogee), (22, Body::IntpPerigee)] {
        let mut st: Stats = Default::default();
        let (mut ss, mut n) = ([0.0f64; 4], [0usize; 4]);
        for r in rows() {
            let jd = r["jd_ut"].as_f64().unwrap();
            let sw = arr(&r[id.to_string()]["default"]);
            let p = e.calc_ut(body, jd, &opt).unwrap().to_array();
            accumulate(&mut st, era(jd), &p, &sw);
            ss[era(jd)] += diff_deg(p[0], sw[0]).powi(2);
            n[era(jd)] += 1;
        }
        print(&format!("intp {id}"), &st);
        for k in 0..4 {
            println!("intp {id} {:10} lon rms {:.4}°", ERAS[k], (ss[k] / n[k] as f64).sqrt());
        }
        for k in 0..4 {
            assert!(st[k][0] <= tol_lon[k][(id - 21) as usize], "{id} {}: {}", ERAS[k], st[k][0]);
            assert!(st[k][1] <= 0.002, "{id} lat {}", st[k][1]);
            assert!(st[k][2] <= 1.5e-5, "{id} dist {}", st[k][2]);
        }
    }
}

#[allow(dead_code)]
fn opts(name: &str) -> CalcOptions {
    let mut o = CalcOptions::default();
    match name {
        "helio" => o.origin = Origin::Heliocentric,
        "nonut_truepos" => {
            o.mean_equinox = true;
            o.corrections = Corrections::GEOMETRIC;
        }
        "j2000" => {
            o.j2000 = true;
            o.mean_equinox = true;
        }
        _ => {}
    }
    o
}

/// Fictitious bodies with Swiss built-in elements (40–54) vs the oracle, all
/// flag variants of the fixture. Tolerances: measured max |Δlon| (deg) per
/// body over 0–2400, rounded up (docs/FICTITIOUS_BODIES.md explains them).
#[test]
fn fictitious_vs_swiss() {
    require_oracle!("fictitious/fictitious.json");
    let Some(e) = engine() else { return };
    let cases: [(i32, f64, f64); 13] = [
        (40, 0.002, 1e-6),
        (41, 0.002, 1e-6),
        (42, 0.002, 1e-6),
        (43, 0.09, 7e-5), // Kronos: Neely's a = 64.81690; Swiss's orbit is ~0.003 AU larger
        (44, 0.002, 1e-6),
        (45, 0.002, 1e-6),
        (46, 0.002, 1e-6),
        (47, 0.002, 1e-6),
        (50, 0.002, 1e-6),
        (51, 0.6, 1.5e-3), // Le Verrier: CR 1846 elements; Swiss uses another variant
        (52, 0.4, 1.2e-3), // Adams hypothesis II
        (53, 0.5, 2.5e-3), // Lowell X1 (via Bower 1931)
        (0, 0.0, 0.0),
    ];
    for (id, tol_lon, tol_dist) in cases {
        if id == 0 {
            continue;
        }
        let body = Body::from_id(id).unwrap();
        for flag in ["default", "helio", "nonut_truepos", "j2000"] {
            let opt = opts(flag);
            let mut st: Stats = Default::default();
            for r in rows() {
                let jd = r["jd_ut"].as_f64().unwrap();
                let sw = arr(&r[id.to_string()][flag]);
                let p = e.calc_ut(body, jd, &opt).unwrap().to_array();
                accumulate(&mut st, era(jd), &p, &sw);
            }
            print(&format!("{id} {flag}"), &st);
            for k in 0..4 {
                assert!(st[k][0] <= tol_lon, "{id} {flag} {}: lon {}", ERAS[k], st[k][0]);
                assert!(st[k][2] <= tol_dist, "{id} {flag} {}: dist {}", ERAS[k], st[k][2]);
            }
        }
    }
}

#[test]
fn unprovided_bodies_are_unknown() {
    for id in [48, 49, 54, 55, 57] {
        assert!(Body::from_id(id).is_err(), "{id}");
    }
}
