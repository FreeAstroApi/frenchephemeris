// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! IAU 1976/1980 frames vs ERFA, and the JPL Horizons modes (`SEFLG_JPLHOR`,
//! `SEFLG_JPLHOR_APPROX`) vs JPL Horizons itself and vs the Swiss oracle.
//! Fixtures: tests/fixtures/horizons/gen_horizons.py (the Swiss part is in the private oracle
//! corpus; without it the Swiss comparison is skipped and the Horizons/ERFA checks still run).

#[macro_use]
mod common;

use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::frames::iau1980::*;
use fe_core::math::*;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::OnceLock;

fn fixture() -> &'static Value {
    static F: OnceLock<Value> = OnceLock::new();
    F.get_or_init(|| {
        let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/horizons/horizons.json");
        serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
    })
}

fn engine() -> Option<&'static Engine> {
    static E: OnceLock<Option<Engine>> = OnceLock::new();
    E.get_or_init(|| {
        let mut c = vec![];
        if let Ok(p) = std::env::var("FE_PACK") {
            c.push(PathBuf::from(p));
        }
        c.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs/fepk-de441-de442-v1.bin"));
        c.push(
            PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join("fepk-de441-de442-v1.bin"),
        );
        c.into_iter().find(|p| p.exists()).map(|p| Engine::open(&p).expect("open pack"))
    })
    .as_ref()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn maxabs(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
}

#[test]
fn iau1980_matches_erfa() {
    let (mut dn, mut dobl, mut dp, mut dg, mut de) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for r in fixture()["erfa"].as_array().unwrap() {
        let jd = f(&r["jd"]);
        let (a, b) = nut80(jd);
        dn = dn.max((a - f(&r["nut80"][0])).abs()).max((b - f(&r["nut80"][1])).abs());
        dobl = dobl.max((obl80(jd) - f(&r["obl80"])).abs());
        let m = pmat76(jd);
        for i in 0..3 {
            for j in 0..3 {
                dp = dp.max((m[i][j] - f(&r["pmat76"][i][j])).abs());
            }
        }
        dg = dg.max(diff_deg(gmst82(jd) / DEG, f(&r["gmst82"]) / DEG).abs() * DEG);
        let (dpsi, _) = nut80(jd);
        de = de.max((eqeq94(jd, dpsi, obl80(jd)) - f(&r["eqeq94"])).abs());
    }
    println!(
        "vs ERFA (rad): nut80 {dn:.1e}, obl80 {dobl:.1e}, pmat76 {dp:.1e}, gmst82 {dg:.1e}, eqeq94 {de:.1e}"
    );
    assert!(dn < 1e-14 && dobl < 1e-14 && dp < 1e-14 && dg < 1e-11 && de < 1e-14);
}

/// (ecliptic lon, lat, RA, Dec) of date in degrees.
fn place(e: &Engine, body: Body, jd_tt: f64, mode: Option<PoleOffsets>) -> [f64; 4] {
    let mut out = [0.0; 4];
    for (k, equatorial) in [(0, false), (2, true)] {
        let opt = CalcOptions { equatorial, speed: false, ..Default::default() };
        let p = match mode {
            Some(m) => e.calc_prepared(body, &e.prepare_horizons(jd_tt, false, false, m), &opt).unwrap(),
            None => e.calc_tt(body, jd_tt, &opt).unwrap(),
        };
        out[k] = p.lon;
        out[k + 1] = p.lat;
    }
    out
}

#[test]
fn horizons_modes_vs_horizons_and_swiss() {
    let Some(e) = engine() else {
        eprintln!("SKIP: no pack");
        return;
    };
    let fx = fixture();
    // Swiss-oracle values (private corpus); `None` skips the Swiss comparison only.
    let swiss = common::oracle_json("horizons/horizons.json");
    if swiss.is_none() {
        eprintln!(
            "SKIP (Swiss part only): oracle fixture horizons/horizons.json not found in {}",
            common::oracle_dir().display()
        );
    }
    let bodies = [
        ("Sun", Body::Sun),
        ("Moon", Body::Moon),
        ("Mars", Body::Mars),
        ("Jupiter", Body::Jupiter),
        ("Saturn", Body::Saturn),
    ];
    // [mode][era] → (ecliptic, equatorial) maximum errors in mas; eras: <1700, ≥1700.
    let names = ["default (Vondrák+2000A)", "JPLHOR (EOP C04)", "JPLHOR_APPROX (model)"];
    let modes = [None, Some(PoleOffsets::EopC04), Some(PoleOffsets::Model)];
    let mut err = [
        [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
        [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
        [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
    ];
    let mut astro = [Vec::new(), Vec::new()];
    let (mut sw_ecl, mut sw_eq) = (vec![], vec![]);
    let mut inner_modern = vec![];
    for (name, body) in bodies {
        let hz = fx["horizons"][name].as_array().unwrap();
        let sw = swiss.as_ref().map(|s| s["swiss_jplhora"][name].as_array().unwrap().clone());
        for (i, row) in hz.iter().enumerate() {
            let jd = f(&row[0]);
            let era = (jd >= 2_341_972.5) as usize; // 1700-01-01
            let (ra0, de0, ra, de, lon, lat) =
                (f(&row[1]), f(&row[2]), f(&row[3]), f(&row[4]), f(&row[5]), f(&row[6]));
            // Astrometric ICRF: pure ephemeris difference (DE441 vs our pack).
            let opt = CalcOptions {
                equatorial: true,
                j2000: true,
                icrs: true,
                speed: false,
                corrections: fe_core::ephem::Corrections {
                    light_time: true,
                    deflection: false,
                    aberration: false,
                },
                ..Default::default()
            };
            let a = e.calc_tt(body, jd, &opt).unwrap();
            astro[era].push((diff_deg(a.lon, ra0) * (de0 * DEG).cos()).hypot(a.lat - de0) * 3.6e6);
            let ast_err = (diff_deg(a.lon, ra0) * (de0 * DEG).cos(), a.lat - de0);
            for (m, mode) in modes.iter().enumerate() {
                let p = place(e, body, jd, *mode);
                let d_ecl = (diff_deg(p[0], lon) * (lat * DEG).cos()).hypot(p[1] - lat) * 3.6e6;
                let d_eq = (diff_deg(p[2], ra) * (de * DEG).cos()).hypot(p[3] - de) * 3.6e6;
                err[m][era][0].push(d_ecl);
                err[m][era][1].push(d_eq);
                if m == 1 && era == 1 {
                    // Remove the ephemeris (DE442 vs Horizons' DE441) difference: the
                    // apparent-minus-astrometric error, RA·cosδ and Dec components (the
                    // basis turns by ≤ 6° in 1700–2400).
                    let app_err = (diff_deg(p[2], ra) * (de * DEG).cos(), p[3] - de);
                    inner_modern.push((app_err.0 - ast_err.0).hypot(app_err.1 - ast_err.1) * 3.6e6);
                }
                if let (2, Some(sw)) = (m, &sw) {
                    let s: Vec<f64> = sw[i].as_array().unwrap().iter().map(f).collect();
                    sw_ecl.push((diff_deg(p[0], s[0]) * (s[1] * DEG).cos()).hypot(p[1] - s[1]) * 3.6e6);
                    sw_eq.push((diff_deg(p[2], s[2]) * (s[3] * DEG).cos()).hypot(p[3] - s[3]) * 3.6e6);
                }
            }
        }
    }
    println!("vs JPL Horizons, apparent geocentric of date, 5 bodies × 40 epochs (max, mas):");
    println!(
        "  astrometric ICRF (ephemeris only)       <1700 {:>8.2}   ≥1700 {:>8.2}",
        maxabs(&astro[0]),
        maxabs(&astro[1])
    );
    for (m, n) in names.iter().enumerate() {
        println!(
            "  {n:<38}  <1700 ecl {:>8.2} eq {:>8.2}   ≥1700 ecl {:>8.2} eq {:>8.2}",
            maxabs(&err[m][0][0]),
            maxabs(&err[m][0][1]),
            maxabs(&err[m][1][0]),
            maxabs(&err[m][1][1])
        );
    }
    println!(
        "  JPLHOR (EOP) RA/Dec, ephemeris difference removed, 1700–2400: max {:.2} mas",
        maxabs(&inner_modern)
    );
    // pyswisseph with its own ephemeris files ignores the Horizons flags (they need an
    // original JPL file), so these values are its default frames.
    println!(
        "  JPLHOR_APPROX vs Swiss (flag ignored by the oracle): ecl max {:.2} mas, eq max {:.2} mas",
        maxabs(&sw_ecl),
        maxabs(&sw_eq)
    );
    assert!(maxabs(&inner_modern) < 15.0);
}
