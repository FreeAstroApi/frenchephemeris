// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Events vs the Swiss Ephemeris 2.10 black-box oracle (fixtures:
//! tests/fixtures/events/gen_events.py) and vs JPL Horizons airless altitudes
//! (tests/fixtures/events/gen_horizons.py). Skips when no pack is present.

#[macro_use]
mod common;

use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::ephem::Observer;
use fe_core::events::*;
use fe_core::math::diff_deg;
use serde_json::Value;
use std::collections::BTreeMap;
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
        Some(p) => Some(Engine::open(&p).expect("open pack")),
        None => {
            eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
            None
        }
    }
}

fn fixture(name: &str) -> Value {
    let p = common::oracle_path(&format!("events/{name}"));
    serde_json::from_str(&std::fs::read_to_string(p).expect("fixture")).unwrap()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

/// Per-group statistics of absolute errors.
#[derive(Default)]
struct Stats(BTreeMap<String, Vec<f64>>);
impl Stats {
    fn add(&mut self, k: &str, v: f64) {
        self.0.entry(k.to_string()).or_default().push(v.abs());
    }
    fn report(&self, unit: &str) -> BTreeMap<String, (usize, f64, f64, f64)> {
        let mut out = BTreeMap::new();
        for (k, v) in &self.0 {
            let mut s = v.clone();
            s.sort_by(f64::total_cmp);
            let n = s.len();
            let p = |q: f64| s[((n as f64 - 1.0) * q).round() as usize];
            eprintln!(
                "{k:<36} n={n:<5} p50 {:>9.4} p95 {:>9.4} p99 {:>9.4} max {:>9.4} {unit}",
                p(0.5),
                p(0.95),
                p(0.99),
                s[n - 1]
            );
            out.insert(k.clone(), (n, p(0.5), p(0.99), s[n - 1]));
        }
        out
    }
}

fn rise_opts(bits: &str, c: &Value) -> RiseSetOptions {
    let atm = Atmosphere { pressure_hpa: f(&c["P"]), temperature_c: f(&c["T"]) };
    let mut o = RiseSetOptions { refraction: Some(atm), horizon_deg: f(&c["hor"]), ..Default::default() };
    match bits {
        "default" => {}
        "center" => o.disc = DiscPoint::Center,
        "bottom" => o.disc = DiscPoint::LowerLimb,
        "norefr" => o.refraction = None,
        "center_norefr" => {
            o.disc = DiscPoint::Center;
            o.refraction = None
        }
        "hindu" => o = RiseSetOptions { horizon_deg: o.horizon_deg, ..RiseSetOptions::hindu() },
        "fixed" => o.fixed_disc_size = true,
        b => panic!("bits {b}"),
    }
    o
}

fn kind(k: &str) -> RiseSetKind {
    match k {
        "rise" => RiseSetKind::Rise,
        "set" => RiseSetKind::Set,
        "mtransit" => RiseSetKind::UpperTransit,
        "itransit" => RiseSetKind::LowerTransit,
        _ => panic!(),
    }
}

#[test]
fn rise_set_vs_swiss() {
    require_oracle!("events/events.json");
    let Some(e) = engine() else { return };
    let fx = fixture("events.json");
    let mut st = Stats::default();
    let (mut mismatch, mut n) = (Vec::new(), 0usize);
    let t0 = std::time::Instant::now();
    for c in fx["rise"].as_array().unwrap() {
        let body = Body::from_id(c["body"].as_i64().unwrap() as i32).unwrap();
        let k = c["kind"].as_str().unwrap();
        let bits = c["bits"].as_str().unwrap();
        let obs = Observer { lon_deg: f(&c["lon"]), lat_deg: f(&c["lat"]), height_m: f(&c["h"]) };
        let r = e.rise_set(body, f(&c["jd"]), kind(k), &obs, &rise_opts(bits, c)).unwrap();
        n += 1;
        let hor = f(&c["hor"]);
        let group = match c["group"].as_str().unwrap() {
            "atmosphere" if hor < 0.0 && f(&c["h"]) > 0.0 => "atmosphere_neg_horizon_elevated",
            "atmosphere" if hor != 0.0 => "atmosphere_nonzero_horizon",
            g => g,
        };
        let swiss_found = c["ret"].as_i64().unwrap() == 0;
        match (r.jd(), swiss_found) {
            (Some(t), true) => {
                let d = (t - f(&c["t"])) * 86400.0;
                let lbl = match body {
                    Body::Sun => "sun",
                    Body::Moon => "moon",
                    _ => "planet",
                };
                st.add(&format!("{group}/{lbl}"), d);
                // Refracted rise/set by latitude band (the published-refraction residual).
                let refracted =
                    !k.contains("transit") && !matches!(bits, "norefr" | "center_norefr" | "hindu");
                let modern = group != "long_range" && !group.starts_with("atmosphere");
                if modern {
                    let a = obs.lat_deg.abs();
                    let band = if a < 30.0 {
                        "00-30"
                    } else if a < 50.0 {
                        "30-50"
                    } else if a < 60.0 {
                        "50-60"
                    } else if a < 70.0 {
                        "60-70"
                    } else {
                        "70-89"
                    };
                    let kind = if refracted { "refracted" } else { "unrefracted+transits" };
                    st.add(&format!("band {kind} |lat| {band}"), d);
                }
                if group == "long_range" {
                    // Remove the sidereal-time model difference (Swiss sidtime vs IAU ERA-based GAST).
                    let ts = f(&c["t"]);
                    let gast = e.prepare_ut(ts, false).epoch.gast().to_degrees();
                    let off_s = diff_deg(f(&c["st"]) * 15.0, gast) / 15.0 * 3600.0;
                    st.add(&format!("long_range_st_corrected/{lbl}"), d - off_s * 0.99727);
                    // ... and the Delta T model difference (Moon: ~13.2°/d RA drift vs 347.8°/d hour angle).
                    let ddt = fe_core::time::delta_t_seconds(f(&c["jd"])) - f(&c["dt"]) * 86400.0;
                    let k = if body == Body::Moon { 13.18 / 347.8 } else { 0.0 };
                    st.add(&format!("long_range_st_dt_corrected/{lbl}"), d - off_s * 0.99727 - k * ddt);
                }
                if d.abs() > 1.0 && std::env::var("FE_DEBUG").is_ok() {
                    eprintln!("BIG {d:.2}s {group} {body:?} {k} {bits} jd {} lon {:.2} lat {:.2} h {} P {} T {} hor {} dt {}", f(&c["jd"]), obs.lon_deg, obs.lat_deg, obs.height_m, c["P"], c["T"], c["hor"], c["dt"]);
                }
                if !k.contains("transit") && !group.starts_with("atmosphere") && group != "long_range" {
                    st.add(&format!("bits/{bits}"), d);
                }
            }
            (None, false) => {}
            _ => mismatch.push(format!(
                "{group} {body:?} {k} {bits} jd {} lat {:.2}: fe {r:?} swiss ret {} t {}",
                f(&c["jd"]),
                obs.lat_deg,
                c["ret"],
                c["t"]
            )),
        }
    }
    let el = t0.elapsed();
    eprintln!(
        "rise/set: {n} cases in {:.1} ms ({:.1} us/case, debug build)",
        el.as_secs_f64() * 1e3,
        el.as_secs_f64() * 1e6 / n as f64
    );
    let rep = st.report("s");
    for m in &mismatch {
        eprintln!("MISMATCH {m}");
    }
    eprintln!("classification mismatches: {} / {n}", mismatch.len());
    // Tolerances (s): ≤ 1 s everywhere in 1900–2100, refracted or not (the
    // refraction model is Sinclair 1982 / Explanatory Supplement 1992 §3.28).
    for (k, (_, _, _, max)) in &rep {
        let tol = if k.starts_with("long_range") {
            60.0
        } else if k.contains("neg_horizon_elevated") {
            f64::INFINITY // dip-aware refraction for elevated observers: not modelled
        } else {
            1.0
        };
        assert!(*max <= tol, "{k}: max {max} s > {tol}");
    }
    assert!(mismatch.is_empty(), "classification mismatches");
}

#[test]
fn pheno_vs_swiss() {
    require_oracle!("events/events.json");
    let Some(e) = engine() else { return };
    let fx = fixture("events.json");
    let names = ["phase_angle", "phase", "elongation", "diameter", "magnitude", "hor_parallax"];
    let mut st = Stats::default();
    let opt = CalcOptions::default();
    let t0 = std::time::Instant::now();
    for c in fx["pheno"].as_array().unwrap() {
        let id = c["body"].as_i64().unwrap() as i32;
        let body = Body::from_id(id).unwrap();
        let jd = f(&c["jd"]);
        // Outside 1900-2030 the Delta T model difference moves the Moon by up to ~100″.
        let era = if jd < 2415020.0 {
            "0-1900"
        } else if jd < 2462502.5 {
            "1900-2030"
        } else {
            "2030-2100"
        };
        let p = e.pheno_ut(body, jd, &opt).unwrap().to_array();
        let a = c["attr"].as_array().unwrap();
        for i in 0..6 {
            if i == 5 && id != 1 {
                continue;
            }
            let d = p[i] - f(&a[i]);
            // angles in arcsec, phase (illuminated fraction) in ppm, magnitude in mag
            let (d, unit) = match i {
                0 | 2 | 3 | 5 => (d * 3600.0, "\""),
                1 => (d * 1e6, "ppm"),
                _ => (d, "mag"),
            };
            let tag = if i == 4 && id == 1 && f(&a[0]) > 150.0 { " (alpha>150)" } else { "" };
            st.add(&format!("{}/{:<7}{tag} {era} {unit}", names[i], format!("{body:?}")), d);
        }
    }
    eprintln!(
        "pheno: {:.1} us/case",
        t0.elapsed().as_secs_f64() * 1e6 / fx["pheno"].as_array().unwrap().len() as f64
    );
    let rep = st.report("");
    for (k, (_, _, _, max)) in rep {
        if !k.contains("1900-2030") || k.contains("alpha>150") {
            continue; // reported only
        }
        let moon = k.contains("Moon");
        let tol = match k.split('/').next().unwrap() {
            "phase_angle" => {
                if moon {
                    1.5
                } else {
                    120.0
                }
            }
            "elongation" => {
                if moon {
                    1.5
                } else {
                    25.0
                }
            }
            "phase" => {
                if moon {
                    10.0
                } else {
                    150.0
                }
            }
            "diameter" => {
                if moon {
                    0.2
                } else {
                    0.05
                }
            } // Moon radius 1737.4 (WGCCRE) vs 1737.5
            "hor_parallax" => 1.0,
            _ => {
                if k.contains("Sun") {
                    0.13
                } else {
                    0.06
                }
            } // magnitude: model differences
        };
        assert!(max <= tol, "{k}: {max} > {tol}");
    }
}

#[test]
fn crossings_vs_swiss() {
    require_oracle!("events/events.json");
    let Some(e) = engine() else { return };
    let fx = fixture("events.json");
    let mut st = Stats::default();
    let opt = CalcOptions::default();
    let t0 = std::time::Instant::now();
    for c in fx["cross"].as_array().unwrap() {
        let body = if c["body"].as_i64() == Some(0) { Body::Sun } else { Body::Moon };
        let jd = f(&c["jd"]);
        let t = e.longitude_cross_ut(body, f(&c["x"]), jd, &opt).unwrap();
        assert!(t > jd);
        let era = if jd > 2415020.0 { "1900-2100" } else { "0-1900" };
        let d = (t - f(&c["t"])) * 86400.0;
        st.add(&format!("{body:?} UT {era}"), d);
        // Same instant in TT: removes the Delta T model difference.
        let ddt = fe_core::time::delta_t_seconds(t) - f(&c["dt"]) * 86400.0;
        st.add(&format!("{body:?} TT {era}"), d + ddt);
    }
    eprintln!(
        "crossings: {:.1} us/case",
        t0.elapsed().as_secs_f64() * 1e6 / fx["cross"].as_array().unwrap().len() as f64
    );
    for (k, (_, _, _, max)) in st.report("s") {
        let tol = match (k.contains(" TT "), k.contains("0-1900")) {
            (true, false) => 0.5,
            (true, true) => 5.0,
            (false, false) => 120.0, // UT: Delta T extrapolation after 2026 differs
            (false, true) => 600.0,
        };
        assert!(max <= tol, "{k}: {max}");
    }
}

#[test]
fn equation_of_time_vs_swiss() {
    require_oracle!("events/events.json");
    let Some(e) = engine() else { return };
    let fx = fixture("events.json");
    let mut st = Stats::default();
    for c in fx["eot"].as_array().unwrap() {
        let jd = f(&c["jd"]);
        let lon = f(&c["lon"]);
        let era = if jd > 2415020.0 { "1900-2100" } else { "0-1900" };
        let d = (e.time_equ(jd).unwrap() - f(&c["e"])) * 86400.0;
        st.add(&format!("time_equ {era}"), d);
        // Remove the sidereal-time model difference (E = GAST − α☉ − UT + 12h).
        let gast = e.prepare_ut(jd, false).epoch.gast().to_degrees();
        let off_s = diff_deg(f(&c["st"]) * 15.0, gast) / 15.0 * 3600.0;
        // ... and the Delta T model difference (α☉ drifts ~1/365.25 of a day per day).
        let ddt = fe_core::time::delta_t_seconds(jd) - f(&c["dt"]) * 86400.0;
        st.add(&format!("time_equ ST+DT-corrected {era}"), d + off_s + ddt / 365.25);
        st.add(
            &format!("lmt_to_lat {era}"),
            (e.lmt_to_lat(jd, lon).unwrap() - f(&c["lmt_to_lat"])) * 86400.0,
        );
        st.add(
            &format!("lat_to_lmt {era}"),
            (e.lat_to_lmt(jd, lon).unwrap() - f(&c["lat_to_lmt"])) * 86400.0,
        );
        let back = e.lmt_to_lat(e.lat_to_lmt(jd, lon).unwrap(), lon).unwrap();
        assert!((back - jd).abs() * 86400.0 < 1e-4);
    }
    for (k, (_, _, _, max)) in st.report("s") {
        let tol = if k.contains("corrected") {
            if k.contains("0-1900") {
                0.2
            } else {
                0.05
            }
        } else if k.contains("0-1900") {
            30.0
        } else {
            1.0
        };
        assert!(max <= tol, "{k}: {max}");
    }
}

#[test]
fn cotrans_and_azalt_vs_swiss() {
    require_oracle!("events/events.json");
    let fx = fixture("events.json");
    let mut st = Stats::default();
    for c in fx["cotrans"].as_array().unwrap() {
        let i: Vec<f64> = c["in"].as_array().unwrap().iter().map(f).collect();
        let o = cotrans([i[0], i[1], i[2]], f(&c["eps"]));
        let ref_o: Vec<f64> = c["out"].as_array().unwrap().iter().map(f).collect();
        st.add("cotrans lon/lat deg", diff_deg(o[0], ref_o[0]).max((o[1] - ref_o[1]).abs()));
        let s = cotrans_sp([i[0], i[1], i[2], i[3], i[4], i[5]], f(&c["eps"]));
        let rs: Vec<f64> = c["out_sp"].as_array().unwrap().iter().map(f).collect();
        st.add("cotrans_sp speeds deg/d", (s[3] - rs[3]).abs().max((s[4] - rs[4]).abs()));
    }
    let rep = st.report("");
    assert!(rep["cotrans lon/lat deg"].3 < 1e-9);
    assert!(rep["cotrans_sp speeds deg/d"].3 < 1e-9);

    let Some(e) = engine() else { return };
    let mut st = Stats::default();
    for c in fx["azalt"].as_array().unwrap() {
        let obs = Observer { lon_deg: f(&c["lon"]), lat_deg: f(&c["lat"]), height_m: f(&c["h"]) };
        let x: Vec<f64> = c["x"].as_array().unwrap().iter().map(f).collect();
        let input = if c["equ"].as_bool().unwrap() { AzaltInput::Equatorial } else { AzaltInput::Ecliptic };
        let atm = Atmosphere { pressure_hpa: f(&c["P"]), temperature_c: f(&c["T"]) };
        let r = e.azalt(f(&c["jd"]), input, &obs, &atm, [x[0], x[1], x[2]]);
        let o: Vec<f64> = c["out"].as_array().unwrap().iter().map(f).collect();
        st.add("azalt az \"", diff_deg(r[0], o[0]) * 3600.0);
        st.add("azalt true alt \"", (r[1] - o[1]) * 3600.0);
        let site = if obs.height_m > 0.0 { "elevated" } else { "sea level" };
        let band = if o[1] < 0.0 {
            "<0"
        } else if o[1] < 15.0 {
            "0-15"
        } else {
            ">=15"
        };
        st.add(&format!("azalt app alt {site} true alt {band} \""), (r[2] - o[2]) * 3600.0);
    }
    let rep = st.report("");
    // Residual = the Swiss sidereal-time offset after ~2050 (0.12 s ≈ 1.8″).
    assert!(rep["azalt true alt \""].3 < 2.0);
    // Apparent altitude above the horizon: Sinclair 1982, same as the oracle (residual = ST).
    for (k, (_, _, _, max)) in &rep {
        if k.starts_with("azalt app alt") && !k.contains("<0") {
            assert!(*max < 2.0, "{k}: {max}");
        }
    }
}

#[test]
fn rise_set_vs_horizons() {
    let Some(e) = engine() else { return };
    let p = format!("{}/tests/fixtures/events/horizons.json", env!("CARGO_MANIFEST_DIR"));
    let Ok(text) = std::fs::read_to_string(p) else {
        eprintln!("SKIP: no horizons.json");
        return;
    };
    let fx: Value = serde_json::from_str(&text).unwrap();
    let mut st = Stats::default();
    for c in fx["cases"].as_array().unwrap() {
        let body = Body::from_id(c["body"].as_i64().unwrap() as i32).unwrap();
        let obs = Observer { lon_deg: f(&c["lon"]), lat_deg: f(&c["lat"]), height_m: f(&c["h"]) };
        // Same fixed horizon on both sides: upper limb at true altitude −34′ (AA convention).
        let o = RiseSetOptions { refraction_model: HorizonRefraction::Almanac34, ..Default::default() };
        let k = if c["kind"].as_str() == Some("rise") { RiseSetKind::Rise } else { RiseSetKind::Set };
        let t = e.rise_set(body, f(&c["jd"]), k, &obs, &o).unwrap().jd().unwrap();
        let d = (t - f(&c["t_horizons"])) * 86400.0;
        st.add(&format!("horizons {body:?}"), d);
    }
    for (k, (_, _, _, max)) in st.report("s") {
        assert!(max <= 1.0, "{k}: {max}");
    }
}

/// Long-range rise/set/transit (years 0–1900) with the Swiss-documented long-term sidereal time
/// and the Swiss-compatible ΔT (`EngineConfig::SWISS_TIME`): no ST or ΔT correction applied.
#[test]
fn long_range_rise_set_swiss_time() {
    require_oracle!("events/events.json");
    let Some(path) = pack_path() else { return };
    let e = Engine::with_config(
        fe_core::pack::Pack::open(&path).unwrap(),
        fe_core::context::EngineConfig::SWISS_TIME,
    );
    let fx = fixture("events.json");
    let mut st = Stats::default();
    for c in fx["rise"].as_array().unwrap() {
        if c["group"].as_str().unwrap() != "long_range" || c["ret"].as_i64().unwrap() != 0 {
            continue;
        }
        let body = Body::from_id(c["body"].as_i64().unwrap() as i32).unwrap();
        let k = c["kind"].as_str().unwrap();
        let bits = c["bits"].as_str().unwrap();
        let obs = Observer { lon_deg: f(&c["lon"]), lat_deg: f(&c["lat"]), height_m: f(&c["h"]) };
        let Some(t) = e.rise_set(body, f(&c["jd"]), kind(k), &obs, &rise_opts(bits, c)).unwrap().jd() else {
            continue;
        };
        let d = (t - f(&c["t"])) * 86400.0;
        let geom = k.contains("transit") || matches!(bits, "norefr" | "center_norefr" | "hindu");
        let lbl = match body {
            Body::Sun => "sun",
            Body::Moon => "moon",
            _ => "planet",
        };
        st.add(&format!("swiss_time/{}/{lbl}", if geom { "geometric" } else { "refracted" }), d);
    }
    let rep = st.report("s");
    for (k, (_, _, _, max)) in &rep {
        // Sun/planets: geometry only. Moon: also the lunar ephemeris (DE441 vs Swiss's
        // DE431-based files) at ancient dates, ~1″ in RA.
        if k.contains("geometric") {
            let tol = if k.ends_with("moon") { 3.0 } else { 1.0 };
            assert!(*max < tol, "{k}: {max} s");
        }
    }
}

/// Refraction formulas vs the oracle's extended (`refrac_extended`, geoalt 0) and
/// simple (`refrac`) conversions, apparent/true altitudes −20…90°.
#[test]
fn refraction_vs_swiss() {
    require_oracle!("events/events.json");
    let fx = fixture("events.json");
    let mut st = Stats::default();
    for c in fx["refrac"].as_array().unwrap() {
        let (h, p, t) = (f(&c["h"]), f(&c["P"]), f(&c["T"]));
        let band = if h < 0.0 {
            "h<0"
        } else if h < 17.9 {
            "0-17.9"
        } else {
            "17.9-90"
        };
        let a2t = refrac(h, p, t, RefracDir::ApparentToTrue);
        let t2a = refrac(h, p, t, RefracDir::TrueToApparent);
        st.add(&format!("Sinclair1982 vs extended A2T {band}"), (a2t - f(&c["ext_a2t"])) * 3600.0);
        st.add(&format!("Sinclair1982 vs extended T2A {band}"), (t2a - f(&c["ext_t2a"])) * 3600.0);
        st.add(&format!("Sinclair1982 vs simple  T2A {band}"), (t2a - f(&c["simple_t2a"])) * 3600.0);
        let s2a = refrac_with(RefractionModel::Saemundsson, h, p, t, RefracDir::TrueToApparent);
        st.add(&format!("Saemundsson  vs simple  T2A {band}"), (s2a - f(&c["simple_t2a"])) * 3600.0);
    }
    let rep = st.report("\"");
    for (k, (_, _, _, max)) in &rep {
        if k.contains("extended") {
            // A2T exact; T2A ≤ 0.1″ (the oracle's own inversion tolerance near 0°).
            assert!(*max < 0.1, "{k}: {max}\"");
        }
    }
}
