// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Time scales vs Swiss (black-box julday/revjul/deltat) and ERFA (dtdb).
//! Fixtures: tests/fixtures/time/gen_time.py.

#[macro_use]
mod common;

use fe_core::time::*;
use serde_json::Value;

fn fixture() -> Value {
    let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/time/time.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// Full fixture including the Swiss-oracle part (private corpus; see tests/common).
fn oracle() -> Value {
    serde_json::from_str(&common::read_oracle("time/time.json")).unwrap()
}

fn bits(v: &Value) -> f64 {
    f64::from_bits(v.as_u64().unwrap())
}

fn cal(v: &Value) -> Calendar {
    if v["cal"] == "G" {
        Calendar::Gregorian
    } else {
        Calendar::Julian
    }
}

#[test]
fn julday_matches_swiss_exactly() {
    require_oracle!("time/time.json");
    let f = oracle();
    let mut n = 0;
    for c in f["julday"].as_array().unwrap() {
        let jd = julday(
            c["y"].as_i64().unwrap() as i32,
            c["m"].as_u64().unwrap() as u32,
            c["d"].as_u64().unwrap() as u32,
            bits(&c["h_bits"]),
            cal(c),
        );
        assert_eq!(jd.to_bits(), bits(&c["jd_bits"]).to_bits(), "{c}");
        n += 1;
    }
    for c in f["known"].as_array().unwrap() {
        let jd = julday(
            c["y"].as_i64().unwrap() as i32,
            c["m"].as_u64().unwrap() as u32,
            c["d"].as_u64().unwrap() as u32,
            c["h"].as_f64().unwrap(),
            cal(c),
        );
        assert!((jd - c["jd"].as_f64().unwrap()).abs() < 1e-9, "{c} -> {jd}");
    }
    assert!(n > 3000);
}

#[test]
fn revjul_matches_swiss() {
    require_oracle!("time/time.json");
    let f = oracle();
    let mut max_dh: f64 = 0.0;
    for c in f["julday"].as_array().unwrap() {
        let (y, m, d, h) = revjul(bits(&c["jd_bits"]), cal(c));
        assert_eq!(
            (y as i64, m as u64, d as u64),
            (c["ry"].as_i64().unwrap(), c["rm"].as_u64().unwrap(), c["rd"].as_u64().unwrap()),
            "{c}"
        );
        max_dh = max_dh.max((h - bits(&c["rh_bits"])).abs());
    }
    eprintln!("revjul max |hour - swe| = {max_dh:e} h");
    assert_eq!(max_dh, 0.0, "revjul hour must be bit-identical to swe");
}

#[test]
fn tdb_vs_erfa_dtdb() {
    let f = fixture();
    let (mut full, mut trunc): (f64, f64) = (0.0, 0.0);
    for c in f["dtdb"].as_array().unwrap() {
        let jd = c["jd"].as_f64().unwrap();
        let v = c["v"].as_f64().unwrap();
        full = full.max((tdb_minus_tt_seconds_full(jd) - v).abs());
        trunc = trunc.max((tdb_minus_tt_seconds(jd) - v).abs());
    }
    eprintln!("TDB-TT max error: full {full:e} s, truncated {trunc:e} s");
    assert!(full < 1e-12);
    assert!(trunc < 0.5e-6);
}

#[test]
fn delta_t_vs_swiss() {
    require_oracle!("time/time.json");
    let f = oracle();
    let mut worst_compat: f64 = 0.0;
    eprintln!("year  swe       default   d(default)  d(compat)");
    for c in f["swe_deltat"].as_array().unwrap() {
        let (y, jd, s) = (c["year"].as_f64().unwrap(), c["jd"].as_f64().unwrap(), c["swe"].as_f64().unwrap());
        let d = delta_t_seconds(jd);
        let k = delta_t_seconds_with(jd, &DeltaTModel::SWISS_COMPAT);
        if (y as i64) % 100 == 99 || (y as i64) % 100 == 9 {
            eprintln!("{y:6.0} {s:9.3} {d:9.3} {:+9.3} {:+9.3}", d - s, k - s);
        }
        // The fixture holds Swiss's "cold" swe_deltat (ṅ −25.80). Inside a chart before
        // 1800 Swiss has its DE441-based file loaded and corrects to ṅ −25.936; expect that.
        let swiss_in_chart =
            if y < 1800.0 { s + fe_core::time::tidal_correction_seconds(jd, -25.80, -25.936) } else { s };
        if y < 1800.0 {
            worst_compat = worst_compat.max((k - swiss_in_chart).abs());
        } else if y < 1950.0 {
            assert!((k - s).abs() < 0.35, "1800-1950 compat within 0.35 s, year {y}");
        }
    }
    eprintln!("SWISS_COMPAT max |dT - swe(in chart)| before 1800: {worst_compat:.4} s");
    assert!(worst_compat < 0.02);
}

#[test]
fn delta_t_sanity() {
    // Observed era: USNO values (J2000: 63.83 s; 2020-01-01: 69.3612 s).
    let j = julday(2020, 1, 1, 0.0, Calendar::Gregorian);
    assert!((delta_t_seconds(j) - 69.3612).abs() < 1e-9);
    assert!((delta_t_seconds(2451545.0) - 63.83).abs() < 0.02);
    // Monotone growth into the future, finite in the API range.
    let a = delta_t_seconds(julday(2100, 1, 1, 0.0, Calendar::Gregorian));
    let b = delta_t_seconds(julday(2400, 1, 1, 0.0, Calendar::Gregorian));
    assert!(a > 69.0 && b > a && b < 2000.0);
    assert!(DELTA_T_MODEL.starts_with("FE-DT/2"));
}

/// Delta T vs Swiss table for the report:
/// `cargo test --release --test time delta_t_report -- --ignored --nocapture`.
#[test]
#[ignore]
fn delta_t_report() {
    require_oracle!("time/time.json");
    let f = oracle();
    println!("| year | swe 2.10 (s) | FE default (s) | default − swe (s) | Moon Δλ (\") | SWISS_COMPAT − swe (s) |");
    println!("|---:|---:|---:|---:|---:|---:|");
    for c in f["swe_deltat"].as_array().unwrap() {
        let (y, jd, s) = (c["year"].as_f64().unwrap(), c["jd"].as_f64().unwrap(), c["swe"].as_f64().unwrap());
        let d = delta_t_seconds(jd);
        let k = delta_t_seconds_with(jd, &DeltaTModel::SWISS_COMPAT);
        println!(
            "| {y:.0} | {s:.3} | {d:.3} | {:+.3} | {:+.2} | {:+.3} |",
            d - s,
            (d - s) * 0.000152 * 3600.0,
            k - s
        );
    }
}
