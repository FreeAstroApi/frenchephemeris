// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Long-term sidereal time (`SiderealTime::SwissLongTerm`) vs the Swiss Ephemeris 2.10
//! black-box oracle (fixture: tests/fixtures/sidtime/gen_sidtime.py). The oracle is used to
//! report residuals only; no constant is derived from it.
//!
//! `cargo test -p fe-core --test sidtime -- --nocapture`

#[macro_use]
mod common;

use fe_core::context::{Engine, EngineConfig};
use fe_core::frames::{self, FrameModel, FrameSet, SiderealTime};
use fe_core::math::ARCSEC;
use serde_json::Value;
use std::path::PathBuf;

fn fixture() -> Value {
    let p = common::oracle_path("sidtime/sidtime.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

fn bits(v: &Value) -> f64 {
    f64::from_bits(v.as_u64().unwrap())
}

fn wrap_as(x_rad: f64) -> f64 {
    let t = std::f64::consts::TAU;
    ((x_rad + t / 2.0).rem_euclid(t) - t / 2.0) / ARCSEC
}

fn wrap_deg(x: f64) -> f64 {
    (x + 180.0).rem_euclid(360.0) - 180.0
}

const LT: FrameModel =
    FrameModel { sidereal: SiderealTime::SwissLongTerm, ..fe_core::context::API_FRAME_MODEL };
const IAU: FrameModel = fe_core::context::API_FRAME_MODEL;

/// Residual (″) of a model's GAST against Swiss `sidtime`, both at Swiss's own ΔT so that only
/// the sidereal-time definition is compared.
fn residual(c: &Value, model: FrameModel) -> (f64, f64) {
    let jd = bits(&c["jd_ut_bits"]);
    let jd_tt = jd + c["dt_s"].as_f64().unwrap() / 86400.0;
    let fs = FrameSet::with_model(jd_tt, model);
    let sw = c["st_h"].as_f64().unwrap() * 15.0 * fe_core::math::DEG;
    let mut r = wrap_as(fs.gast(jd) - sw);
    if model.sidereal == SiderealTime::SwissLongTerm
        && (jd < frames::LONG_TERM_ST_LIMITS.0 || jd >= frames::LONG_TERM_ST_LIMITS.1)
    {
        // The long-term model takes its own Delta T (`LONG_TERM_ST_DELTA_T`); remove the
        // known ΔT-extrapolation difference (0.0411″ per second, mean solar motion) so that
        // only the sidereal-time definition is compared.
        let own = fe_core::time::delta_t_seconds_with(jd, &frames::LONG_TERM_ST_DELTA_T);
        r -= 0.041_067 * (own - c["dt_s"].as_f64().unwrap());
    }
    (jd, r)
}

fn year(jd: f64) -> f64 {
    2000.0 + (jd - 2451545.0) / 365.2425
}

#[test]
fn long_term_sidereal_time_vs_swiss() {
    require_oracle!("sidtime/sidtime.json");
    let fx = fixture();
    let cases = fx["cases"].as_array().unwrap();
    let (lo, hi) = frames::LONG_TERM_ST_LIMITS;
    let off = frames::long_term_offsets();
    println!(
        "continuity constants: before 1850 {:+.4}″, from 2050 {:+.4}″",
        off[0] / ARCSEC,
        off[1] / ARCSEC
    );
    // bands: (label, from year, to year, max |LT − swe|, max |IAU − swe|, n)
    let mut bands: Vec<(&str, f64, f64, f64, f64, usize)> = vec![
        ("-600..0", -601.0, 0.0, 0.0, 0.0, 0),
        ("0..1000", 0.0, 1000.0, 0.0, 0.0, 0),
        ("1000..1850", 1000.0, 1850.0, 0.0, 0.0, 0),
        ("1850..2050", 1850.0, 2050.0, 0.0, 0.0, 0),
        ("2050..2400", 2050.0, 2401.0, 0.0, 0.0, 0),
    ];
    let mut before = vec![];
    let mut after = vec![];
    for c in cases {
        let (jd, r) = residual(c, LT);
        let (_, ri) = residual(c, IAU);
        let y = year(jd);
        for b in bands.iter_mut() {
            if y >= b.1 && y < b.2 {
                b.3 = b.3.max(r.abs());
                b.4 = b.4.max(ri.abs());
                b.5 += 1;
            }
        }
        if jd < lo {
            before.push(r);
        } else if jd >= hi {
            after.push(r);
        } else {
            // Inside 1850–2050 both models are the IAU one; Swiss is gst06a to ~1 mas and the
            // API frames (Vondrák + 77-term nutation) differ from gst06a by a few mas.
            assert!(r.abs() < 0.01, "inside {y}: {r}″");
        }
    }
    println!("max |GAST − swe.sidtime| (″), both at Swiss's ΔT:");
    println!("  {:>12} {:>5} {:>12} {:>12}", "years", "n", "LongTerm", "IAU 2006");
    for b in &bands {
        println!("  {:>12} {:>5} {:>12.4} {:>12.3}", b.0, b.5, b.3, b.4);
    }
    let spread = |v: &[f64]| {
        let mn = v.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    };
    let (b0, b1) = spread(&before);
    let (a0, a1) = spread(&after);
    println!(
        "  before 1850: residual {b0:+.4} … {b1:+.4}″  (constant part = continuity constant difference)"
    );
    println!("  from 2050:   residual {a0:+.4} … {a1:+.4}″");
    // The published method reproduces the shape of Swiss's function to ≪ 1″: the residual on
    // each side is a constant (the two continuity constants are Swiss's own, unpublished).
    assert!(b1 - b0 < 0.01, "before-1850 residual not constant: {b0} … {b1}");
    assert!(a1 - a0 < 0.01, "after-2050 residual not constant: {a0} … {a1}");
    assert!(b0.abs().max(b1.abs()) < 0.15, "before 1850: {b0} … {b1}");
    assert!(a0.abs().max(a1.abs()) < 2.5, "after 2050: {a0} … {a1}");
}

#[test]
fn long_term_is_continuous_at_limits() {
    let (lo, hi) = frames::LONG_TERM_ST_LIMITS;
    let m = fe_core::time::DeltaTModel::SWISS_COMPAT;
    for b in [lo, hi] {
        let e = 1e-7;
        let g = |jd: f64| {
            let fs = FrameSet::with_model(fe_core::time::ut_to_tt_with(jd, &m), LT);
            fs.gast(jd)
        };
        let jump = wrap_as(g(b + e) - g(b - e)) - 2.0 * e * 1.002_737_909 * 1_296_000.0;
        assert!(jump.abs() < 0.005, "jump at {b}: {jump}″");
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

/// End to end through `Engine::houses_ut` with each engine's own ΔT.
#[test]
fn houses_long_term_vs_swiss() {
    require_oracle!("sidtime/sidtime.json");
    let Some(path) = pack_path() else {
        eprintln!("SKIP: no pack");
        return;
    };
    let fx = fixture();
    let cfgs = [
        ("SWISS_DELTA_T (IAU ST)", EngineConfig::SWISS_DELTA_T),
        ("SWISS_TIME (long-term ST)", EngineConfig::SWISS_TIME),
    ];
    let mut worst = vec![];
    for (name, cfg) in cfgs {
        let eng = Engine::with_config(fe_core::pack::Pack::open(&path).unwrap(), cfg);
        let (mut armc_max, mut cusp_max, mut armc_post, mut cusp_post) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for h in fx["houses"].as_array().unwrap() {
            let jd = bits(&h["jd_ut_bits"]);
            let (lat, lon) = (h["lat"].as_f64().unwrap(), h["lon"].as_f64().unwrap());
            let hs = eng.houses_ut(jd, lat, lon, 'P', None).unwrap();
            let da = wrap_deg(hs.armc - h["armc"].as_f64().unwrap()).abs() * 3600.0;
            let mut dc = wrap_deg(hs.asc - h["asc"].as_f64().unwrap())
                .abs()
                .max(wrap_deg(hs.mc - h["mc"].as_f64().unwrap()).abs());
            for (i, c) in h["cusps"].as_array().unwrap().iter().enumerate().take(12) {
                dc = dc.max(wrap_deg(hs.cusps[i] - c.as_f64().unwrap()).abs());
            }
            if jd >= frames::LONG_TERM_ST_LIMITS.1 {
                armc_post = armc_post.max(da);
                cusp_post = cusp_post.max(dc);
            } else {
                armc_max = armc_max.max(da);
                cusp_max = cusp_max.max(dc);
            }
        }
        println!(
            "{name}: before 2050 max ΔARMC {armc_max:.3}″, max Δcusp/ASC/MC {cusp_max:.5}°; after 2050 {armc_post:.3}″, {cusp_post:.5}°"
        );
        worst.push((armc_max, cusp_max));
    }
    // Before 2050 the long-term model removes the sidereal-time difference; what remains is
    // ΔT (SWISS_COMPAT = Swiss to ~0.006 s before 1950) and the continuity constant (0.1″).
    assert!(worst[1].0 < 0.5, "ARMC before 2050 {}″", worst[1].0);
    assert!(worst[1].1 < 0.001, "cusps before 2050 {}°", worst[1].1);
}
