// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Pack reader vs direct JPL evaluation (fixtures from jplephem, see
//! tests/fixtures/pack/gen_pack_fixtures.py).
//!
//! Packs are gitignored build products. The tests look for them in `$FE_PACK`
//! (a single file), `$FE_DATA_DIR/packs`, `<repo>/data/packs` and the main
//! checkout's data dir, and skip with a message when none is present.

use fe_core::pack::Pack;
use fe_core::Error;
use std::path::PathBuf;

const PACK_NAMES: [&str; 2] = ["fepk-de441-de442-v1.bin", "fepk-de442-1900-2100-v1.bin"];

fn pack_paths() -> Vec<PathBuf> {
    if let Ok(p) = std::env::var("FE_PACK") {
        return vec![PathBuf::from(p)].into_iter().filter(|p| p.exists()).collect();
    }
    let mut dirs = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        dirs.push(PathBuf::from(d).join("packs"));
    }
    dirs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    dirs.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    let mut out = Vec::new();
    for name in PACK_NAMES {
        if let Some(p) = dirs.iter().map(|d| d.join(name)).find(|p| p.exists()) {
            out.push(p);
        }
    }
    out
}

struct Row {
    target: i32,
    center: i32,
    jd: f64,
    p: [f64; 3],
    v: [f64; 3],
    region: String,
}

fn fixtures() -> Vec<Row> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/pack/pack_states.json");
    let text = std::fs::read_to_string(path).expect("fixture file");
    let json: serde_json::Value = serde_json::from_str(&text).unwrap();
    json["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| {
            let f = |i: usize| r[i].as_f64().unwrap();
            Row {
                target: r[0].as_i64().unwrap() as i32,
                center: r[1].as_i64().unwrap() as i32,
                // exact epoch from its bit pattern (decimal parse may be 1 ulp off)
                jd: f64::from_bits(r[10].as_u64().unwrap()),
                p: [f(3), f(4), f(5)],
                v: [f(6), f(7), f(8)],
                region: r[9].as_str().unwrap().to_string(),
            }
        })
        .collect()
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[test]
fn pack_matches_jpl_fixtures() {
    let paths = pack_paths();
    if paths.is_empty() {
        eprintln!("SKIP pack_matches_jpl_fixtures: no pack found (build with tools/build_pack.py)");
        return;
    }
    let rows = fixtures();
    for path in paths {
        let pack = Pack::open(&path).unwrap();
        let (mut n, mut max_dp, mut max_dv, mut max_rel) = (0usize, 0.0f64, 0.0f64, 0.0f64);
        let mut worst = String::new();
        for r in &rows {
            let Some((lo, hi)) = pack.chain_coverage(r.target, r.center) else {
                panic!("pair {}/{} not chainable", r.target, r.center);
            };
            if r.jd < lo || r.jd > hi {
                continue;
            }
            let (p, v) = pack.state(r.target, r.center, r.jd).unwrap();
            let p_only = pack.position(r.target, r.center, r.jd).unwrap();
            assert_eq!(p, p_only, "position() and state() disagree");
            let dp = norm([p[0] - r.p[0], p[1] - r.p[1], p[2] - r.p[2]]);
            let dv = norm([v[0] - r.v[0], v[1] - r.v[1], v[2] - r.v[2]]);
            let rel = dp / norm(r.p);
            // Tolerances: pure floating-point roundoff. Absolute 1e-6 km, or
            // 1e-15 relative for ~1e9 km barycentric outer-planet vectors.
            let tol_p = 1e-6_f64.max(1e-15 * norm(r.p));
            let tol_v = 1e-6_f64.max(1e-13 * norm(r.v));
            assert!(
                dp <= tol_p && dv <= tol_v,
                "{}: {}/{} jd {} ({}): dp {dp:e} km dv {dv:e} km/d",
                path.display(),
                r.target,
                r.center,
                r.jd,
                r.region
            );
            if dp > max_dp {
                max_dp = dp;
                worst = format!("{}/{} jd {} ({})", r.target, r.center, r.jd, r.region);
            }
            max_dv = max_dv.max(dv);
            max_rel = max_rel.max(rel);
            n += 1;
        }
        assert!(n > 0, "no fixture inside the coverage of {}", path.display());
        eprintln!(
            "{}: {n} states checked; max |dp| {max_dp:.3e} km (at {worst}), max |dv| {max_dv:.3e} km/day, \
             max rel {max_rel:.3e}",
            path.file_name().unwrap().to_string_lossy()
        );
    }
}

#[test]
fn pack_coverage_and_errors() {
    let paths = pack_paths();
    if paths.is_empty() {
        eprintln!("SKIP pack_coverage_and_errors: no pack found");
        return;
    }
    for path in paths {
        let pack = Pack::open(&path).unwrap();
        let m: serde_json::Value = serde_json::from_str(pack.manifest()).unwrap();
        assert_eq!(m["format_version"], 1);
        for (t, c) in [
            (1, 0),
            (2, 0),
            (3, 0),
            (4, 0),
            (5, 0),
            (6, 0),
            (7, 0),
            (8, 0),
            (9, 0),
            (10, 0),
            (301, 3),
            (399, 3),
        ] {
            let (lo, hi) = pack.coverage(t, c).expect("stored pair");
            assert!(pack.state(t, c, lo).is_ok());
            assert!(pack.state(t, c, hi).is_ok());
            let below = lo - 1e-3;
            assert_eq!(pack.state(t, c, below), Err(Error::OutOfRange { body: t, jd: below }));
            assert!(matches!(pack.position(t, c, hi + 1e-3), Err(Error::OutOfRange { .. })));
        }
        let cov = m["coverage_jd_tdb"].as_array().unwrap();
        let (lo, hi) = (cov[0].as_f64().unwrap(), cov[1].as_f64().unwrap());
        // every chainable pair works over the advertised common coverage
        let bodies = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399];
        for &t in &bodies {
            for &c in &bodies {
                for jd in [lo, 0.5 * (lo + hi), hi] {
                    pack.state(t, c, jd).unwrap();
                }
            }
        }
        // Mercury/Venus aliases equal their barycenters
        assert_eq!(pack.position(199, 0, lo).unwrap(), pack.position(1, 0, lo).unwrap());
        assert!(matches!(pack.position(500, 0, lo), Err(Error::UnknownBody(500))));
        // Send + Sync: share across threads
        let pack = std::sync::Arc::new(pack);
        let hs: Vec<_> = (0..4)
            .map(|k| {
                let p = pack.clone();
                std::thread::spawn(move || p.position(301, 399, lo + 10.0 * k as f64 + 1.0).unwrap())
            })
            .collect();
        for h in hs {
            h.join().unwrap();
        }
    }
}

/// Timing: `cargo test --release -p fe-core --test pack -- --ignored --nocapture`
#[test]
#[ignore]
fn pack_eval_speed() {
    let Some(path) = pack_paths().into_iter().next() else {
        eprintln!("SKIP pack_eval_speed: no pack found");
        return;
    };
    let pack = Pack::open(&path).unwrap();
    let (lo, hi) = pack.chain_coverage(301, 399).unwrap();
    const N: usize = 1_000_000;
    // pseudo-random times spread over the coverage (LCG, no deps)
    let mut x: u64 = 0x9E3779B97F4A7C15;
    let jds: Vec<f64> = (0..N)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            lo + (hi - lo) * ((x >> 11) as f64 / (1u64 << 53) as f64)
        })
        .collect();
    let seq: Vec<f64> = (0..N).map(|i| 2451545.0 + i as f64 * 0.01).filter(|t| *t <= hi).collect();
    eprintln!("pack {}", path.file_name().unwrap().to_string_lossy());
    for (label, t, c) in [
        ("Moon wrt EMB (stored)", 301, 3),
        ("Moon geocentric 301/399", 301, 399),
        ("Mars bary wrt SSB (stored)", 4, 0),
        ("Mars geocentric 4/399", 4, 399),
    ] {
        for (tl, times) in [("random", &jds), ("sequential", &seq)] {
            let t0 = std::time::Instant::now();
            let mut acc = 0.0;
            for &jd in times.iter() {
                acc += std::hint::black_box(pack.position(t, c, std::hint::black_box(jd)).unwrap())[0];
            }
            let ns_p = t0.elapsed().as_nanos() as f64 / times.len() as f64;
            let t0 = std::time::Instant::now();
            for &jd in times.iter() {
                acc += std::hint::black_box(pack.state(t, c, std::hint::black_box(jd)).unwrap()).1[0];
            }
            let ns_s = t0.elapsed().as_nanos() as f64 / times.len() as f64;
            eprintln!(
                "{label:<28} {tl:<10} position {ns_p:6.1} ns/call   state {ns_s:6.1} ns/call  ({acc:.0})"
            );
        }
    }
}
