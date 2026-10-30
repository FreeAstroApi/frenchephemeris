// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Mean obliquity: Vondrák, Capitaine & Wallace (2011) published ε_A series vs the Swiss
//! Ephemeris 2.10 oracle (SE_ECL_NUT, sidtime, Placidus houses) over −13000 … +17000, and vs
//! ERFA `obl06` near J2000. Fixture: tests/fixtures/obliquity/gen_obliquity.py.
//!
//! `cargo test -p fe-core --test obliquity -- --nocapture`

#[macro_use]
mod common;

use fe_core::context::API_FRAME_MODEL;
use fe_core::frames::{self, FrameModel, FrameSet, Obliquity, SiderealTime};
use fe_core::houses;
use fe_core::math::{diff_deg, ARCSEC, DEG};
use serde_json::Value;

fn fixture() -> Value {
    let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/obliquity/obliquity.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// Full fixture including the Swiss-oracle part (private corpus; see tests/common).
fn oracle() -> Value {
    serde_json::from_str(&common::read_oracle("obliquity/obliquity.json")).unwrap()
}
fn bits(v: &Value) -> f64 {
    f64::from_bits(v.as_u64().unwrap())
}
fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

#[test]
fn vondrak_obliquity_matches_published_value_at_j2000() {
    // Eq. (10) + Table 3 at T = 0 must give the IAU 2006 J2000 obliquity 84381.406″ (to the
    // 1e-6″ rounding of the published coefficients).
    let e = frames::vondrak_obliquity(2451545.0) / ARCSEC;
    assert!((e - 84381.406).abs() < 2e-6, "{e}");
}

#[test]
fn vondrak_obliquity_vs_iau2006_near_j2000() {
    let fx = fixture();
    let mut mx: f64 = 0.0;
    for c in fx["obl06"].as_array().unwrap() {
        let jd = bits(&c["jd_bits"]);
        let d = (frames::vondrak_obliquity(jd) - f(&c["obl06"])) / ARCSEC;
        assert!((frames::obl06(jd) - f(&c["obl06"])).abs() < 1e-15);
        mx = mx.max(d.abs());
    }
    println!("max |ε_A(Vondrák) − obl06| over 0 … 3000: {mx:.4}″");
    assert!(mx < 0.025, "{mx}");
}

#[test]
fn obliquity_vs_swiss() {
    require_oracle!("obliquity/obliquity.json");
    let fx = oracle();
    let swiss_model = API_FRAME_MODEL;
    assert_eq!(swiss_model.obliquity, Obliquity::Vondrak2011);
    assert_eq!(FrameModel::SWISS_COMPAT.obliquity, Obliquity::Vondrak2011);
    let poles = FrameModel { obliquity: Obliquity::Consistent, ..API_FRAME_MODEL };
    let (mut m_mean, mut m_true, mut m_poles, mut m_cusp, mut m_cusp_poles): (f64, f64, f64, f64, f64) =
        (0.0, 0.0, 0.0, 0.0, 0.0);
    println!(
        "{:>8} {:>14} {:>14} {:>14} {:>12}",
        "year", "ε_V − swe ″", "ε_poles − swe ″", "true − swe ″", "Plac. cusp ″"
    );
    for c in fx["swiss"].as_array().unwrap() {
        let jd_tt = bits(&c["jd_bits"]);
        let y = f(&c["year"]);
        let fs = FrameSet::with_model(jd_tt, swiss_model);
        let fp = FrameSet::with_model(jd_tt, poles);
        let d_mean = (fs.eps_mean / DEG - f(&c["eps_mean"])) * 3600.0;
        let d_true = (fs.eps_true / DEG - f(&c["eps_true"])) * 3600.0;
        let d_poles = (fp.eps_mean / DEG - f(&c["eps_mean"])) * 3600.0;
        m_mean = m_mean.max(d_mean.abs());
        m_true = m_true.max(d_true.abs());
        m_poles = m_poles.max(d_poles.abs());

        // Houses: Swiss's own ARMC, our true obliquity at Swiss's TT (= jd_ut + Swiss ΔT), so
        // only the obliquity (and nutation) enter.
        let jd_ut = jd_tt;
        let tt = jd_ut + f(&c["dt_s"]) / 86400.0;
        let cusp_err = |m: FrameModel| -> f64 {
            let e = FrameSet::with_model(tt, m).eps_true / DEG;
            let (h, _) = houses::houses_swiss_compat('P', f(&c["armc"]), 47.0, e).unwrap();
            let mut w: f64 = diff_deg(h.asc, f(&c["asc"])).abs();
            for i in 0..12 {
                w = w.max(diff_deg(h.cusps[i], f(&c["cusps"][i])).abs());
            }
            w * 3600.0
        };
        let ce = cusp_err(swiss_model);
        let cp = cusp_err(poles);
        m_cusp = m_cusp.max(ce);
        m_cusp_poles = m_cusp_poles.max(cp);
        if (y % 2000.0).abs() < 1e-9 || [-1000.0, 0.0, 1000.0, 3000.0, 14000.0].contains(&y) {
            println!("{y:8.0} {d_mean:14.7} {d_poles:14.4} {d_true:14.7} {ce:8.4} (poles {cp:.2})");
        }
    }
    println!("max |ε_mean − swe|: Vondrák series {m_mean:.2e}″, Vondrák poles {m_poles:.3}″");
    println!("max |ε_true − swe|: {m_true:.2e}″;  max Placidus cusp/Asc error: {m_cusp:.2e}″ (poles: {m_cusp_poles:.2}″)");
    assert!(m_mean < 1e-6, "{m_mean}");
    assert!(m_true < 1e-6, "{m_true}");
    assert!(m_cusp < 1e-4, "{m_cusp}");
}

#[test]
fn long_term_sidereal_time_extended_range() {
    require_oracle!("obliquity/obliquity.json");
    let fx = oracle();
    let model = FrameModel { sidereal: SiderealTime::SwissLongTerm, ..API_FRAME_MODEL };
    let own = frames::LONG_TERM_ST_DELTA_T;
    let mut rows = vec![];
    for c in fx["swiss"].as_array().unwrap() {
        let jd = bits(&c["jd_bits"]);
        let dt = f(&c["dt_s"]);
        let fs = FrameSet::with_model(jd + dt / 86400.0, model);
        let sw = f(&c["st_h"]) * 15.0 * DEG;
        let t = std::f64::consts::TAU;
        let mut r = ((fs.gast(jd) - sw + t / 2.0).rem_euclid(t) - t / 2.0) / ARCSEC;
        if jd < frames::LONG_TERM_ST_LIMITS.0 || jd >= frames::LONG_TERM_ST_LIMITS.1 {
            r -= 0.041_067 * (fe_core::time::delta_t_seconds_with(jd, &own) - dt);
        }
        rows.push((f(&c["year"]), r));
    }
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut mx: f64 = 0.0;
    for (y, r) in &rows {
        mx = mx.max(r.abs());
        if (y % 2000.0).abs() < 1e-9 {
            println!("sidtime {y:8.0}: GAST − swe = {r:+.4}″");
        }
    }
    println!("max |GAST(long-term) − swe.sidtime| over −13000 … 17000: {mx:.4}″");
    // With the Vondrák ε_A in the mean-Sun projection the residual is a constant on each side
    // (Swiss's unpublished continuity constants); with obl06 it drifted to 23″.
    let side = |lo: f64, hi: f64| {
        let v: Vec<f64> = rows.iter().filter(|r| r.0 >= lo && r.0 < hi).map(|r| r.1).collect();
        let mn = v.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    };
    let (b0, b1) = side(-14000.0, 1849.0);
    let (a0, a1) = side(2051.0, 18000.0);
    println!("before 1850: {b0:+.4} … {b1:+.4}″; after 2050: {a0:+.4} … {a1:+.4}″");
    // After 2050 the residual also carries the first-order correction's error for the
    // difference between our Delta T forecast (FE-DT/2, IERS + 10-year rate) and Swiss's.
    assert!(b1 - b0 < 0.03 && a1 - a0 < 0.04, "residual not constant");
    assert!(mx < 2.05, "{mx}");
}
