// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Ayanamsas (published definitions) vs the Swiss oracle
//! (fixtures from tests/fixtures/ayanamsa/gen_ayanamsa.py).
//!
//! Tolerances are the measured "published model vs Swiss" deltas; they are
//! not accuracy claims about either side.

#[macro_use]
mod common;

use fe_core::ayanamsa::{ayanamsa_mean_deg, true_pushya_deg, Ayanamsa};
use fe_core::math::diff_deg;
use serde_json::Value;

fn rows() -> Vec<Value> {
    let p = common::oracle_path("ayanamsa/ayanamsa.json");
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

fn print(name: &str, e: &[f64; 4], signed: &[(f64, f64); 4]) {
    let cells: Vec<String> = (0..4)
        .map(|i| {
            format!(
                "{}: {:.2e}° [{:+.3}\"..{:+.3}\"]",
                ERAS[i],
                e[i],
                signed[i].0 * 3600.0,
                signed[i].1 * 3600.0
            )
        })
        .collect();
    println!("{name:13} {}", cells.join("  "));
}

#[test]
fn fixed_epoch_modes_vs_swiss() {
    require_oracle!("ayanamsa/ayanamsa.json");
    // (mode, max |Δ| over 0–2400 measured, rounded up)
    let modes = [
        ("Lahiri", Ayanamsa::Lahiri, 2.5e-6),
        ("Raman", Ayanamsa::Raman, 2e-8),
        ("Krishnamurti", Ayanamsa::Krishnamurti, 2e-8),
        ("FaganBradley", Ayanamsa::FaganBradley, 2e-8),
        ("Yukteshwar", Ayanamsa::Yukteshwar, 2e-8),
    ];
    let rows = rows();
    for (name, mode, tol) in modes {
        let mut e = [0.0f64; 4];
        let mut sg = [(f64::MAX, f64::MIN); 4];
        for r in &rows {
            let jd = r["jd_tt"].as_f64().unwrap();
            let d = diff_deg(ayanamsa_mean_deg(mode, jd), r[name][0].as_f64().unwrap());
            let k = era(jd);
            e[k] = e[k].max(d.abs());
            sg[k] = (sg[k].0.min(d), sg[k].1.max(d));
            // true = mean + Δψ holds identically on both sides
            let dt = diff_deg(
                ayanamsa_mean_deg(mode, jd) + r["dpsi"].as_f64().unwrap(),
                r[name][1].as_f64().unwrap(),
            );
            assert!((dt - d).abs() < 1e-9);
        }
        print(name, &e, &sg);
        for v in e {
            assert!(v <= tol, "{name}: {v:e} > {tol:e}");
        }
    }
}

fn v3(v: &Value, off: usize) -> [f64; 3] {
    [v[off].as_f64().unwrap(), v[off + 1].as_f64().unwrap(), v[off + 2].as_f64().unwrap()]
}

#[test]
fn true_pushya_vs_swiss() {
    require_oracle!("ayanamsa/ayanamsa.json");
    use fe_core::ayanamsa::star::delta_cnc_apparent_mean_from_states;
    let rows = rows();
    let mut defn = 0.0f64;
    let (mut eb, mut es) = ([0.0f64; 4], [0.0f64; 4]);
    let (mut sb, mut ss) = ([(f64::MAX, f64::MIN); 4], [(f64::MAX, f64::MIN); 4]);
    for r in &rows {
        let jd = r["jd_tt"].as_f64().unwrap();
        let mean = r["TruePushya"][0].as_f64().unwrap();
        let tru = r["TruePushya"][1].as_f64().unwrap();
        // the definition itself (star lon − 106°) is exact
        defn = defn.max(diff_deg(true_pushya_deg(r["delta_cnc"][0].as_f64().unwrap()), mean).abs());
        defn = defn.max(diff_deg(true_pushya_deg(r["delta_cnc"][2].as_f64().unwrap()), tru).abs());
        // skip the hours when δ Cnc is behind the solar disc (Swiss attenuates deflection)
        let sep = diff_deg(r["sun_lon"].as_f64().unwrap(), r["delta_cnc"][0].as_f64().unwrap()).abs();
        if sep < 0.3 {
            continue;
        }
        let k = era(jd);
        let d = diff_deg(ayanamsa_mean_deg(Ayanamsa::TruePushya, jd), mean);
        assert_eq!(ayanamsa_mean_deg(Ayanamsa::Pushya, jd), ayanamsa_mean_deg(Ayanamsa::TruePushya, jd));
        eb[k] = eb[k].max(d.abs());
        sb[k] = (sb[k].0.min(d), sb[k].1.max(d));
        let e = &r["earth_bary"];
        let (lon, _) = delta_cnc_apparent_mean_from_states(jd, v3(e, 0), v3(e, 3), v3(&r["sun_bary"], 0));
        let d = diff_deg(true_pushya_deg(lon), mean);
        es[k] = es[k].max(d.abs());
        ss[k] = (ss[k].0.min(d), ss[k].1.max(d));
    }
    println!("TruePushya definition (star lon − 106°): {defn:.2e}°");
    print("TP built-in", &eb, &sb);
    print("TP states", &es, &ss);
    assert!(defn < 1e-11);
    let tol = [6.0e-6, 6.0e-6, 6.0e-6, 6.0e-6];
    for k in 0..4 {
        // built-in low-precision Earth vs exact (oracle) Earth/Sun states
        assert!(eb[k] <= tol[k] && es[k] <= 1.5e-6, "{}: {:e} {:e}", ERAS[k], eb[k], es[k]);
    }
}
