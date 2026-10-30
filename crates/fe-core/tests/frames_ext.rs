// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Frames at extreme epochs (−13200 … +17191) vs ERFA: Vondrák 2011 precession is
//! valid over ±200 kyr; IAU 2000A nutation and the ERA are evaluated as published.
//! Fixtures: tests/fixtures/frames_ext/gen_frames_ext.py.

use fe_core::frames::{era, ltp, ltpb, nut00a};

#[test]
fn ltp_nutation_era_match_erfa_over_the_extended_range() {
    let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/frames_ext/frames_ext.json");
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    let (mut dm, mut dn, mut de) = (0.0f64, 0.0f64, 0.0f64);
    for r in v["rows"].as_array().unwrap() {
        let jd = f64::from_bits(r["jd_bits"].as_u64().unwrap());
        let a = ltpb(jd);
        let (b, _) = ltp(jd);
        for i in 0..3 {
            for j in 0..3 {
                dm = dm.max((a[i][j] - r["ltpb"][i][j].as_f64().unwrap()).abs());
                dm = dm.max((b[i][j] - r["ltp"][i][j].as_f64().unwrap()).abs());
            }
        }
        let (dp, dq) = nut00a(jd);
        dn = dn
            .max((dp - r["nut00a"][0].as_f64().unwrap()).abs())
            .max((dq - r["nut00a"][1].as_f64().unwrap()).abs());
        let d = (era(jd) - r["era00"].as_f64().unwrap()).rem_euclid(std::f64::consts::TAU);
        de = de.max(d.min(std::f64::consts::TAU - d));
    }
    eprintln!("max |Δ| ltp/ltpb {dm:.2e}, nut00a {dn:.2e} rad, era {de:.2e} rad");
    assert!(dm < 1e-14, "precession {dm}");
    assert!(dn < 1e-13, "nutation {dn}");
    assert!(de < 1e-10, "era {de}");
}
