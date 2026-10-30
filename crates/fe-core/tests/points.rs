// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Lunar nodes / Lilith vs the Swiss oracle (fixtures from tests/fixtures/points/gen_points.py).

#[macro_use]
mod common;

use fe_core::math::diff_deg;
use fe_core::points::*;
use serde_json::Value;

fn rows() -> Vec<Value> {
    let p = common::oracle_path("points/points.json");
    let v: Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    v["rows"].as_array().unwrap().clone()
}

fn g(v: &Value, i: usize) -> f64 {
    v[i].as_f64().unwrap()
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

/// Published ELP/MPP02 + IAU 2006 p_A mean elements vs the oracle. The
/// tolerances are the measured deltas (arcsec), per era.
#[test]
fn mean_points_vs_swiss() {
    require_oracle!("points/points.json");
    // [node lon, apogee lon, apogee lat] max |Δ| in arcsec per era (measured, rounded up)
    let tol: [[f64; 3]; 4] =
        [[1.21, 2.91, 0.26], [0.15, 2.69, 0.25], [0.035, 0.55, 0.05], [0.15, 0.36, 0.03]];
    let mut e = [[0.0f64; 3]; 4];
    let mut node_v = 0.0f64;
    for r in rows() {
        let jd = r["jd_tt"].as_f64().unwrap();
        let k = era(jd);
        e[k][0] = e[k][0].max(diff_deg(mean_node_deg(jd), g(&r["mean_node"], 0)).abs() * 3600.0);
        let (l, b) = mean_apogee_deg(jd);
        e[k][1] = e[k][1].max(diff_deg(l, g(&r["mean_apogee"], 0)).abs() * 3600.0);
        e[k][2] = e[k][2].max((b - g(&r["mean_apogee"], 1)).abs() * 3600.0);
        node_v = node_v.max((mean_node_speed_deg_per_day(jd) - g(&r["mean_node"], 1)).abs());
        // conventional distances: 384 400 km and 384 400 km × (1 + 0.054900489)
        assert!((g(&r["mean_apogee"], 2) - MEAN_APOGEE_DIST_AU).abs() < 1e-15);
    }
    for k in 0..4 {
        println!(
            "{:10} mean node {:.3}\"  mean apogee lon {:.3}\"  lat {:.3}\"",
            ERAS[k], e[k][0], e[k][1], e[k][2]
        );
        for j in 0..3 {
            assert!(e[k][j] <= tol[k][j], "{} {j}: {} > {}", ERAS[k], e[k][j], tol[k][j]);
        }
    }
    println!("mean node speed max |Δ| {node_v:.2e}°/day");
    assert!(node_v < 1e-6, "{node_v:e}");
}

#[test]
fn osculating_points_match_swiss() {
    require_oracle!("points/points.json");
    // GM(Earth+Moon) the oracle uses for its osculating elements (measured,
    // test-only — the engine uses its ephemeris' GM).
    let gm = 3.98600448e5 * (1.0 + 1.0 / 81.30056) * 86400.0 * 86400.0;
    let (mut node, mut node_r, mut apo_l, mut apo_b, mut apo_r) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for r in rows() {
        let p = [g(&r["moon_pos_km"], 0), g(&r["moon_pos_km"], 1), g(&r["moon_pos_km"], 2)];
        let v = [g(&r["moon_vel_km_day"], 0), g(&r["moon_vel_km_day"], 1), g(&r["moon_vel_km_day"], 2)];
        let n = osculating_node_deg(p, v);
        node = node.max(diff_deg(n, g(&r["true_node"], 0)).abs());
        let (n2, dist) = osculating_node_full(p, v, gm);
        assert_eq!(n, n2);
        node_r = node_r.max((dist / g(&r["true_node"], 1) - 1.0).abs());
        let (l, b, d) = osculating_apogee_deg(p, v, gm);
        apo_l = apo_l.max(diff_deg(l, g(&r["oscu_apogee"], 0)).abs());
        apo_b = apo_b.max((b - g(&r["oscu_apogee"], 1)).abs());
        apo_r = apo_r.max((d / g(&r["oscu_apogee"], 2) - 1.0).abs());
    }
    println!("true node     max |Δlon| {node:.2e}°  rel |Δdist| {node_r:.2e}");
    println!("oscu apogee   max |Δlon| {apo_l:.2e}°  |Δlat| {apo_b:.2e}°  rel |Δdist| {apo_r:.2e}");
    assert!(node < 1e-9, "{node:e}");
    assert!(node_r < 1e-9, "{node_r:e}");
    assert!(apo_l < 1e-6 && apo_b < 1e-7, "{apo_l:e} {apo_b:e}");
    assert!(apo_r < 1e-9, "{apo_r:e}");
}
