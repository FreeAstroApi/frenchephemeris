// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Frames vs ERFA (pyerfa) and vs Swiss 2.10 frames recovered black-box.
//! Fixtures: tests/fixtures/frames/gen_frames.py.

#[macro_use]
mod common;

use fe_core::frames::*;
use fe_core::math::{Mat3, ARCSEC};
use serde_json::Value;

const UAS: f64 = ARCSEC * 1e-6;

fn fixture() -> Value {
    let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/frames/frames.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// Full fixture including the Swiss-oracle part (private corpus; see tests/common).
fn oracle() -> Value {
    serde_json::from_str(&common::read_oracle("frames/frames.json")).unwrap()
}
fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}
fn bits(v: &Value) -> f64 {
    f64::from_bits(v.as_u64().unwrap())
}
fn mat(v: &Value) -> Mat3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = f(&v[i][j]);
        }
    }
    m
}
fn dmat(a: &Mat3, b: &Mat3) -> f64 {
    let mut d: f64 = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            d = d.max((a[i][j] - b[i][j]).abs());
        }
    }
    d
}
fn dang(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(std::f64::consts::TAU);
    d.min(std::f64::consts::TAU - d)
}

struct Max(Vec<(&'static str, f64, f64)>);
impl Max {
    fn upd(&mut self, k: &'static str, v: f64, tol: f64) {
        if let Some(e) = self.0.iter_mut().find(|e| e.0 == k) {
            e.1 = e.1.max(v);
        } else {
            self.0.push((k, v, tol));
        }
    }
    fn check(&self) {
        let mut ok = true;
        for (k, v, tol) in &self.0 {
            eprintln!("{k:<28} max err {v:.3e} rad ({:.4} uas)  tol {tol:.1e}", v / UAS);
            ok &= v <= tol;
        }
        assert!(ok, "tolerance exceeded");
    }
}

#[test]
fn erfa_reference_values() {
    let fx = fixture();
    let mut m = Max(vec![]);
    for c in fx["erfa"].as_array().unwrap() {
        let jd = bits(&c["jd_tt_bits"]);
        let ut = bits(&c["jd_ut_bits"]);
        m.upd("obl06", (obl06(jd) - f(&c["obl06"])).abs(), 1e-15);
        let (g, p, s) = pfw06(jd);
        m.upd(
            "pfw06",
            (g - f(&c["pfw06"][0]))
                .abs()
                .max((p - f(&c["pfw06"][1])).abs())
                .max((s - f(&c["pfw06"][2])).abs()),
            1e-15,
        );
        let n = nut06a(jd);
        m.upd("nut06a", (n.0 - f(&c["nut06a"][0])).abs().max((n.1 - f(&c["nut06a"][1])).abs()), 1e-13);
        let n = nut00a(jd);
        m.upd("nut00a", (n.0 - f(&c["nut00a"][0])).abs().max((n.1 - f(&c["nut00a"][1])).abs()), 1e-13);
        let n = nut00b(jd);
        m.upd("nut00b", (n.0 - f(&c["nut00b"][0])).abs().max((n.1 - f(&c["nut00b"][1])).abs()), 1e-13);
        assert_eq!(nutation(jd), nut06a(jd));

        let rb = frame_bias();
        m.upd("bp06 rb (frame bias)", dmat(&rb, &mat(&c["bp06_rb"])), 1e-15);
        let fs = FrameSet::with_model(jd, FrameModel::IAU2006_2000A);
        m.upd("bp06 rbp", dmat(&fs.bp, &mat(&c["bp06_rbp"])), 1e-13);
        m.upd("pnm06a", dmat(&fs.bpn, &mat(&c["pnm06a"])), 1e-12);
        m.upd("gst06a (FrameSet)", dang(fs.gast(ut), f(&c["gst06a"])), 1e-12);
        m.upd("gst06a (free fn)", dang(gast(ut, jd), f(&c["gst06a"])), 1e-12);
        m.upd("gmst06 (FrameSet)", dang(fs.gmst(ut), f(&c["gmst06"])), 1e-12);
        m.upd("ee06a", (fs.eqeq - f(&c["ee06a"])).abs(), 1e-12);
        m.upd("ee06a (free fn)", (equation_of_equinoxes(jd) - f(&c["ee06a"])).abs(), 1e-12);
        let (x, y) = (fs.bpn[2][0], fs.bpn[2][1]);
        m.upd("s06", (s06(jd, x, y) - f(&c["s06"])).abs(), 1e-14);
        m.upd("eect00", (eect00(jd) - f(&c["eect00"])).abs(), 1e-15);
        m.upd("era00", dang(era(ut), f(&c["era00"])), 1e-12);
        m.upd("gmst06 (free fn)", dang(gmst(ut, jd), f(&c["gmst06"])), 1e-12);

        let fb = FrameSet::with_model(jd, FrameModel::IAU2006_2000B);
        m.upd("fw2m(pfw06+nut00b)", dmat(&fb.bpn, &mat(&c["pn06_00b"])), 1e-12);
        m.upd("gst06(pn 06/00b)", dang(fb.gast(ut), f(&c["gst06_00b"])), 1e-12);

        // Vondrak 2011
        let pec = ltp_ecliptic_pole(jd);
        let peq = ltp_equator_pole(jd);
        for k in 0..3 {
            m.upd("ltpecl", (pec[k] - f(&c["ltpecl"][k])).abs(), 1e-14);
            m.upd("ltpequ", (peq[k] - f(&c["ltpequ"][k])).abs(), 1e-14);
        }
        m.upd("ltp", dmat(&ltp(jd).0, &mat(&c["ltp"])), 1e-14);
        m.upd("ltpb", dmat(&ltpb(jd), &mat(&c["ltpb"])), 1e-14);
        let fv = FrameSet::new(jd);
        m.upd("DEFAULT bp == ltpb", dmat(&fv.bp, &mat(&c["ltpb"])), 1e-14);
    }
    m.check();
}

#[test]
fn ecliptic_j2000() {
    let e = gcrs_to_ecliptic_j2000();
    let b = frame_bias();
    // The ecliptic pole (row 2) is the bias-transformed equatorial pole rotated by 84381.406".
    let (s, c) = (84381.406 * ARCSEC).sin_cos();
    for j in 0..3 {
        assert!((e[2][j] - (-s * b[1][j] + c * b[2][j])).abs() < 1e-16);
    }
}

#[test]
fn swiss_compat_frames() {
    require_oracle!("frames/frames.json");
    let fx = oracle();
    let mut m = Max(vec![]);
    for c in fx["swiss"].as_array().unwrap() {
        let jd = bits(&c["jd_tt_bits"]);
        let sw = FrameSet::with_model(jd, FrameModel::SWISS_COMPAT);
        let df = FrameSet::new(jd);
        m.upd("swiss bias vs frame_bias", dmat(&frame_bias(), &mat(&c["bias"])), 2e-11);
        m.upd("swiss bp vs COMPAT", dmat(&sw.bp, &mat(&c["bp"])), 2e-11);
        m.upd("swiss bpn vs COMPAT", dmat(&sw.bpn, &mat(&c["bpn"])), 2e-11);
        // SWISS_COMPAT uses the Vondrak 2011 published eps_A series, which is Swiss's mean obliquity.
        m.upd("swiss ecl_true vs COMPAT", dmat(&sw.ecl_true, &mat(&c["ecl_true"])), 2e-11);
        m.upd("swiss ecl_j2000", dmat(&gcrs_to_ecliptic_j2000(), &mat(&c["ecl_j2000"])), 2e-11);
        // Default model: identical precession, nutation within a few mas.
        m.upd("swiss bp vs DEFAULT", dmat(&df.bp, &mat(&c["bp"])), 2e-11);
        m.upd("swiss bpn vs DEFAULT", dmat(&df.bpn, &mat(&c["bpn"])), 1.5e-8);
        // DEFAULT uses the Vondrak pole-angle obliquity (0.2" from Swiss's at year 0).
        m.upd("swiss ecl_true vs DEFAULT", dmat(&df.ecl_true, &mat(&c["ecl_true"])), 1.2e-6);
    }
    m.check();
}

/// Release-mode timing: `cargo test --release --test frames -- --ignored --nocapture`.
#[test]
#[ignore]
fn timing_frameset() {
    use std::hint::black_box;
    use std::time::Instant;
    for (name, model) in [
        ("DEFAULT (Vondrak+2000A)", FrameModel::DEFAULT),
        ("VONDRAK_2000B", FrameModel::VONDRAK_2000B),
        ("IAU2006_2000A", FrameModel::IAU2006_2000A),
        ("IAU2006_2000B", FrameModel::IAU2006_2000B),
        ("SWISS_COMPAT", FrameModel::SWISS_COMPAT),
    ] {
        let n = 20_000;
        let t0 = Instant::now();
        let mut acc = 0.0;
        for i in 0..n {
            let fs = FrameSet::with_model(black_box(2451545.0 + i as f64 * 0.37), model);
            acc += fs.bpn[0][1];
        }
        let dt = t0.elapsed().as_nanos() as f64 / n as f64;
        black_box(acc);
        let t1 = Instant::now();
        for i in 0..n {
            acc += bpn_rate(black_box(2451545.0 + i as f64 * 0.37), model)[0][1];
        }
        let dr = t1.elapsed().as_nanos() as f64 / n as f64;
        black_box(acc);
        eprintln!("{name:<26} FrameSet::with_model {dt:8.0} ns   bpn_rate {dr:8.0} ns");
    }
    let n = 20_000;
    let t0 = Instant::now();
    let mut acc = 0.0;
    for i in 0..n {
        acc += fe_core::time::tdb_minus_tt_seconds(black_box(2451545.0 + i as f64 * 0.37));
    }
    eprintln!("tdb_minus_tt_seconds (87 terms) {:.0} ns", t0.elapsed().as_nanos() as f64 / n as f64);
    let t0 = Instant::now();
    for i in 0..n {
        acc += fe_core::time::delta_t_seconds(black_box(2451545.0 + i as f64 * 0.37 - 3e5));
    }
    eprintln!("delta_t_seconds {:.0} ns", t0.elapsed().as_nanos() as f64 / n as f64);
    black_box(acc);
}
