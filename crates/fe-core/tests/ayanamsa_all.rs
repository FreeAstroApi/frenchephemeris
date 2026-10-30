// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! All Swiss sidereal modes (0–46), SE_SIDM_USER and the SIDBIT options:
//! published definitions vs the Swiss oracle
//! (fixtures: tests/fixtures/ayanamsa_all/gen_ayanamsa_all.py).
//!
//! `cargo test -p fe-core --test ayanamsa_all -- --nocapture` prints the
//! per-mode delta tables quoted in docs/GEOMETRY_REPORT.md. Tolerances are
//! the measured "published definition vs Swiss" deltas (rounded up), not
//! accuracy claims.

#[macro_use]
mod common;

use fe_core::ayanamsa::{self, Ayanamsa, Kind, SidBits};
use fe_core::context::{Body, CalcOptions, Engine, EngineConfig};
use fe_core::math::diff_deg;
use fe_core::pack::Pack;
use serde_json::Value;
use std::path::PathBuf;

fn fx() -> Value {
    let p = common::oracle_path("ayanamsa_all/ayanamsa_all.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

fn engine() -> Option<Engine> {
    let name = "fepk-de441-de442-v1.bin";
    let dirs = [
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"),
        PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()),
    ];
    let p = std::env::var("FE_PACK").map(PathBuf::from).ok().filter(|p| p.exists());
    let p = p.or_else(|| dirs.iter().map(|d| d.join(name)).find(|p| p.exists()))?;
    Some(Engine::with_config(Pack::open(&p).ok()?, EngineConfig::SWISS_TIME))
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

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

/// Max |Δ| per era in arcsec, and signed range over all eras.
#[derive(Default, Clone, Copy)]
struct Acc {
    e: [f64; 4],
    lo: f64,
    hi: f64,
}
impl Acc {
    fn add(&mut self, jd: f64, d_deg: f64) {
        let d = d_deg * 3600.0;
        let k = era(jd);
        self.e[k] = self.e[k].max(d.abs());
        self.lo = self.lo.min(d);
        self.hi = self.hi.max(d);
    }
    fn row(&self) -> String {
        format!(
            "{:>9.4} {:>9.4} {:>9.4} {:>9.4}   [{:+.4}..{:+.4}]",
            self.e[0], self.e[1], self.e[2], self.e[3], self.lo, self.hi
        )
    }
}

/// Measured max |Δ| (arcsec, any era) of the published definition vs Swiss,
/// built-in (no pack) path. NaN = not implemented (no published definition).
fn tolerance(m: Ayanamsa) -> f64 {
    use Ayanamsa::*;
    match m {
        Lahiri => 0.01,
        LahiriIcrc => 0.07,
        BabylEtaPiscium => 0.94,
        Hipparchos | BabylKugler1 | BabylKugler2 | BabylKugler3 | BabylHuber => 0.015,
        UshaShashi => 0.005,
        DjwhalKhul => 5.0,
        TrueCitra | TrueRevati | TrueMula | TruePushya => 0.04,
        GalCent0Sag | GalCentCochrane | GalCentRgBrand | GalCentMulaWilhelm => 1.0,
        GalEquIau1958 => 0.32,
        GalEquTrue | GalEquMula => 0.12,
        BabylAldebaran15Tau => 0.025,
        JnBhasin | TrueSheoran | ValensMoon => f64::NAN,
        _ => 0.002,
    }
}

#[test]
fn all_modes_builtin_vs_swiss() {
    require_oracle!("ayanamsa_all/ayanamsa_all.json");
    let fx = fx();
    let rows = fx["rows"].as_array().unwrap();
    println!(
        "mode                               0–1000   1000–1800 1800–2100 2100–2400 (max |Δ| arcsec)  signed"
    );
    for (i, &m) in Ayanamsa::all().iter().enumerate() {
        let mut acc = Acc { lo: f64::MAX, hi: f64::MIN, ..Default::default() };
        for r in rows {
            let jd = f(&r["jd_tt"]);
            let ours = ayanamsa::ayanamsa_mean_deg(m, jd);
            if ours.is_nan() {
                break;
            }
            acc.add(jd, diff_deg(ours, f(&r["mean"][i])));
        }
        let tol = tolerance(m);
        if tol.is_nan() {
            println!("{i:>2} {:32} not implemented (no published definition available)", m.name());
            assert!(ayanamsa::ayanamsa_mean_deg(m, 2451545.0).is_nan());
            continue;
        }
        println!("{i:>2} {:32} {}", m.name(), acc.row());
        for e in acc.e {
            assert!(e <= tol, "{}: {e}\" > {tol}\"", m.name());
        }
    }
}

#[test]
fn star_and_galactic_modes_engine_vs_swiss() {
    require_oracle!("ayanamsa_all/ayanamsa_all.json");
    let Some(eng) = engine() else { return };
    let fx = fx();
    let rows = fx["rows"].as_array().unwrap();
    for (i, &m) in Ayanamsa::all().iter().enumerate() {
        if !matches!(
            ayanamsa::kind(m, SidBits(0)),
            Kind::Star(..) | Kind::DeltaCnc(..) | Kind::GalCenter(..) | Kind::GalCenterPolar(..)
        ) {
            continue;
        }
        let mut acc = Acc { lo: f64::MAX, hi: f64::MIN, ..Default::default() };
        let mut acc_true = 0.0f64;
        for r in rows {
            let jd = f(&r["jd_tt"]);
            let ep = eng.prepare_tt(jd, false).epoch;
            let Ok(a) = eng.ayanamsa_deg(m, &ep, false) else { continue };
            acc.add(jd, diff_deg(a, f(&r["mean"][i])));
            let t = eng.ayanamsa_deg(m, &ep, true).unwrap();
            acc_true =
                acc_true.max((diff_deg(t, f(&r["true"][i])) - diff_deg(a, f(&r["mean"][i]))).abs() * 3600.0);
        }
        println!("engine {i:>2} {:32} {}  (true−mean consistency {acc_true:.1e}\")", m.name(), acc.row());
        assert!(acc_true < 0.002);
        for e in acc.e {
            assert!(e <= tolerance(m), "{}: {e}", m.name());
        }
    }
}

#[test]
fn option_bits_vs_swiss() {
    require_oracle!("ayanamsa_all/ayanamsa_all.json");
    let fx = fx();
    let combos: Vec<(i32, u32)> = fx["bits_combos"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| (c[0].as_i64().unwrap() as i32, c[1].as_u64().unwrap() as u32))
        .collect();
    for (k, &(n, b)) in combos.iter().enumerate() {
        let m = Ayanamsa::from_swiss(n).unwrap();
        let mut acc = Acc { lo: f64::MAX, hi: f64::MIN, ..Default::default() };
        for r in fx["bits_rows"].as_array().unwrap() {
            let jd = f(&r["jd_tt"]);
            acc.add(jd, diff_deg(ayanamsa::ayanamsa_mean_deg_bits(m, jd, SidBits(b)), f(&r["v"][k])));
        }
        println!("bits {n:>2}|{b:<5} {:22} {}", m.name(), acc.row());
        // ECL_DATE and NO_PREC_OFFSET reproduce Swiss to the base-mode residual
        // (Lahiri −0.008", ICRC see report); PREC_ORIG (historical precession
        // models) is a documented test option, measured only.
        let tol = match (n, b) {
            (46, 2048) => 0.07,
            (46, 4096) => 0.07,
            (_, 8192) => 0.8,
            _ => 0.01,
        };
        assert!(acc.e.iter().all(|&e| e <= tol), "{n}|{b}");
    }
    // USER mode (t0 TT and t0 UT)
    let mut e = [0.0f64; 2];
    for u in fx["user"].as_array().unwrap() {
        let jd = f(&u[0]);
        let a =
            ayanamsa::ayanamsa_mean_deg(Ayanamsa::User { t0: 2440000.5, ayan_t0: 22.5, t0_is_ut: false }, jd);
        let b = ayanamsa::ayanamsa_mean_deg_bits(
            Ayanamsa::User { t0: 2440000.5, ayan_t0: 22.5, t0_is_ut: false },
            jd,
            SidBits(SidBits::USER_UT),
        );
        e[0] = e[0].max(diff_deg(a, f(&u[1])).abs() * 3600.0);
        e[1] = e[1].max(diff_deg(b, f(&u[2])).abs() * 3600.0);
    }
    println!("user TT t0 {:.5}\"  user UT t0 {:.5}\"", e[0], e[1]);
    assert!(e[0] < 1e-4 && e[1] < 1e-3);
}

#[test]
fn projected_positions_vs_swiss() {
    require_oracle!("ayanamsa_all/ayanamsa_all.json");
    let Some(eng) = engine() else { return };
    let fx = fx();
    let mut worst: std::collections::BTreeMap<(i64, i64), (f64, f64, f64, f64)> = Default::default();
    for p in fx["proj"].as_array().unwrap() {
        let (n, b) = (p["mode"].as_i64().unwrap(), p["bits"].as_i64().unwrap());
        let flags = p["flags"].as_i64().unwrap();
        let opt = CalcOptions {
            sidereal: Some(Ayanamsa::from_swiss(n as i32).unwrap()),
            sid_bits: SidBits(b as u32),
            mean_equinox: flags & 64 != 0,
            ..Default::default()
        };
        let body = Body::from_id(p["body"].as_i64().unwrap() as i32).unwrap();
        let jd = f(&p["jd_ut"]);
        // Compare the sidereal *transformation* (sidereal − tropical, same flags),
        // so that planetary-ephemeris and Delta T differences cancel.
        let r = eng.calc_ut(body, jd, &opt).unwrap();
        let t = eng.calc_ut(body, jd, &CalcOptions { sidereal: None, ..opt }).unwrap();
        let tr = &p["trop"];
        let w = worst.entry((n, b)).or_insert((0.0, 0.0, 0.0, 0.0));
        let dl = diff_deg(diff_deg(r.lon, t.lon), diff_deg(f(&p["lon"]), f(&tr[0])));
        w.0 = w.0.max(dl.abs() * 3600.0);
        w.1 = w.1.max(((r.lat - t.lat) - (f(&p["lat"]) - f(&tr[1]))).abs() * 3600.0);
        w.2 = w.2.max(((r.lon_speed - t.lon_speed) - (f(&p["lon_speed"]) - f(&tr[2]))).abs() * 3600.0);
        let eopt = CalcOptions { equatorial: true, ..opt };
        let re = eng.calc_ut(body, jd, &eopt).unwrap();
        let te = eng.calc_ut(body, jd, &CalcOptions { sidereal: None, ..eopt }).unwrap();
        let tq = &p["trop_eq"];
        let da = diff_deg(diff_deg(re.lon, te.lon), diff_deg(f(&p["ra"]), f(&tq[0])));
        let dd = (re.lat - te.lat) - (f(&p["dec"]) - f(&tq[1]));
        // (lunar points: the engine's tropical EQUATORIAL output for them is
        // ecliptic — a pre-existing gap outside this module — so skip RA/Dec)
        if !matches!(body, Body::MeanNode) {
            w.3 = w.3.max(da.abs().max(dd.abs()) * 3600.0);
        }
        if std::env::var("FE_DEBUG").is_ok() && da.abs().max(dd.abs()) * 3600.0 > 0.01 {
            println!(
                "    debug mode {n} body {body:?} flags {flags} jd {jd} dRA {:.4} dDec {:.4}",
                da * 3600.0,
                dd * 3600.0
            );
        }
    }
    println!("sidereal − tropical transformation, ours vs Swiss (max |Δ|, 1800–2100):");
    for ((n, b), (l, bt, s, e)) in &worst {
        println!(
            "  mode {n:>2} bits {b:>4}: lon {l:.4}\" lat {bt:.4}\" lon speed {s:.4}\"/d  ra/dec {e:.4}\""
        );
        let tol_lon = match n {
            1 => 0.01,
            27 | 31 => 0.2,
            _ => 0.002,
        };
        let tol_lat = if *b == 512 { 0.06 } else { 0.003 };
        assert!(*l < tol_lon && *bt < tol_lat && *e < 0.003 && *s < 0.01, "mode {n} bits {b}");
    }
}

#[test]
fn sidereal_fixed_stars_vs_swiss() {
    require_oracle!("ayanamsa_all/ayanamsa_all.json");
    let Some(eng) = engine() else { return };
    let fx = fx();
    let mut worst: std::collections::BTreeMap<(i64, i64), f64> = Default::default();
    for r in fx["stars"].as_array().unwrap() {
        let (n, b) = (r["mode"].as_i64().unwrap(), r["bits"].as_i64().unwrap());
        let jd = f(&r["jd_ut"]);
        let name = r["star"].as_str().unwrap();
        // Ecliptic only: for fixed stars the oracle's SIDEREAL|EQUATORIAL output
        // is not the mean-equator RA/Dec it returns for planets (RA shifted by
        // ≈ the ayanamsa); we keep the planet rule for stars.
        {
            let (key, eq) = ("ecl", false);
            let opt = CalcOptions {
                sidereal: Some(Ayanamsa::from_swiss(n as i32).unwrap()),
                sid_bits: SidBits(b as u32),
                equatorial: eq,
                speed: false,
                ..Default::default()
            };
            let s = eng.fixstar_ut(name, jd, &opt).unwrap().1;
            let t = eng.fixstar_ut(name, jd, &CalcOptions { sidereal: None, ..opt }).unwrap().1;
            let x = &r[key];
            // compare the sidereal transformation (star catalogues differ)
            let dl = diff_deg(diff_deg(s.lon, t.lon), diff_deg(f(&x[0]), f(&x[2])));
            let db = (s.lat - t.lat) - (f(&x[1]) - f(&x[3]));
            let w = worst.entry((n, b)).or_insert(0.0);
            *w = w.max(dl.abs().max(db.abs()) * 3600.0);
        }
    }
    for ((n, b), w) in &worst {
        println!("fixstar sidereal mode {n:>2} bits {b:>3}: {w:.4}\" (ecliptic lon/lat)");
        let tol = if *b == 512 {
            0.06
        } else if *n == 1 {
            0.01
        } else {
            0.003
        };
        assert!(*w < tol, "mode {n} bits {b}: {w}");
    }
}

#[test]
fn mardyks_epoch_is_the_1998_equinox() {
    let Some(eng) = engine() else { return };
    let opt = CalcOptions { speed: false, ..Default::default() };
    let sun = eng.calc_tt(Body::Sun, ayanamsa::MARDYKS_T0_TT, &opt).unwrap();
    // 1e-6 day ≈ 0.1 s ≈ 4e-8° of solar motion
    assert!(diff_deg(sun.lon, 180.0).abs() < 1e-6, "{}", diff_deg(sun.lon, 180.0));
}
