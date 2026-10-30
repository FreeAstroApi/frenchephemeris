// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Extended-range packs (DE441 refits, FEPK v2) chained with the main pack.
//! Skips with a message when the packs are absent (they are gitignored build
//! products of tools/build_ext_pack.py).

use fe_core::pack::Pack;
use fe_core::Error;
use std::path::PathBuf;

const MAIN: &str = "fepk-de441-de442-v1.bin";
const PAST: &str = "fepk-de441-ext-past-v1.bin";
const FUTURE: &str = "fepk-de441-ext-future-v1.bin";
const EMRAT: f64 = 81.300_568_221_497_22;

fn dir() -> Option<PathBuf> {
    let mut c = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        c.push(PathBuf::from(d).join("packs"));
    }
    c.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    c.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    c.into_iter().find(|d| d.join(MAIN).exists() && d.join(PAST).exists() && d.join(FUTURE).exists())
}

fn open() -> Option<Pack> {
    let Some(d) = dir() else {
        eprintln!("SKIP: extended packs not found (build with tools/build_ext_pack.py)");
        return None;
    };
    Some(Pack::open_many(&[d.join(MAIN), d.join(PAST), d.join(FUTURE)]).unwrap())
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Fit tolerance (km) per stored pair (tools/build_ext_pack.py `BODIES`).
fn tol(t: i32) -> f64 {
    match t {
        301 => 0.1,
        3 => 0.04,
        1 => 0.3,
        2 => 0.15,
        4 => 0.2,
        10 => 0.6,
        _ => 0.9,
    }
}

#[test]
fn coverage_chaining_and_out_of_range() {
    let Some(p) = open() else { return };
    let (lo, hi) = p.chain_coverage(301, 399).unwrap();
    assert!(lo < -3_099_000.0 && hi > 7_990_000.0, "{lo} {hi}");
    for t in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] {
        let (a, b) = p.coverage(t, 0).unwrap();
        assert!(a < -3_099_000.0 && b > 7_990_000.0, "{t}: {a} {b}");
        assert!(matches!(p.position(t, 0, a - 1.0), Err(Error::OutOfRange { .. })));
        assert!(matches!(p.position(t, 0, b + 1.0), Err(Error::OutOfRange { .. })));
    }
    let (a, b) = p.coverage(399, 3).unwrap();
    assert!(matches!(p.state(399, 3, a - 0.5), Err(Error::OutOfRange { body: 399, .. })));
    assert!(p.state(399, 3, b).is_ok());
}

#[test]
fn earth_alias_is_minus_moon_over_emrat() {
    let Some(p) = open() else { return };
    for jd in [-3_000_000.25, -1_000_000.5, 1_000_000.0, 3_000_000.0, 7_000_000.75] {
        let (m, mv) = p.state(301, 3, jd).unwrap();
        let (e, ev) = p.state(399, 3, jd).unwrap();
        for k in 0..3 {
            assert!((e[k] + m[k] / EMRAT).abs() < 1e-9 && (ev[k] + mv[k] / EMRAT).abs() < 1e-9);
        }
    }
}

#[test]
fn seams_with_main_pack_are_continuous() {
    let Some(p) = open() else { return };
    let d = dir().unwrap();
    let main = Pack::open(&d.join(MAIN)).unwrap();
    for t in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 301] {
        let c = if t == 301 { 3 } else { 0 };
        let (lo, hi) = main.coverage(t, c).unwrap();
        // Past seam: DE441 (main, bit-identical) vs DE441 refit: same source.
        // (ε = 1e-7 d; the motion over ε is removed with the velocity)
        const E: f64 = 1e-7;
        let a = p.position(t, c, lo).unwrap(); // main pack (seam epoch -> later series)
        let (b, v) = p.state(t, c, lo - E).unwrap(); // refit, at its last record end
        let b = [b[0] + v[0] * E, b[1] + v[1] * E, b[2] + v[2] * E];
        let jump = norm(sub(a, b));
        assert!(jump < tol(t), "past seam {t}: {jump} km");
        // Future seam: DE442 (main) -> DE441 refit: the DE442−DE441 difference (Uranus ~1e4 km).
        let (a, v) = p.state(t, c, hi - E).unwrap();
        let a = [a[0] + v[0] * E, a[1] + v[1] * E, a[2] + v[2] * E];
        let b = p.position(t, c, hi).unwrap();
        let jump = norm(sub(a, b));
        eprintln!("future seam {t}/{c} at JD {hi}: DE442 -> DE441 jump {jump:.3} km");
        assert!(jump < 2.0e5, "future seam {t}: {jump}");
    }
}

#[test]
fn matches_horizons_de441_vectors() {
    let Some(p) = open() else { return };
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/pack_ext/horizons_vectors.json");
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut worst = 0.0f64;
    for r in v["rows"].as_array().unwrap() {
        let t = r[0].as_i64().unwrap() as i32;
        let jd = f64::from_bits(r[6].as_u64().unwrap());
        let h = [r[3].as_f64().unwrap(), r[4].as_f64().unwrap(), r[5].as_f64().unwrap()];
        let ours = p.position(t, 399, jd).unwrap();
        let dk = norm(sub(ours, h));
        // geocentric error = target error + EMB error + Earth (Moon/EMRAT) error, plus
        // Horizons' 16-digit print (1e-16 relative)
        let budget = if t == 301 {
            0.1
        } else if t == 3 {
            0.1 / EMRAT + 0.01
        } else {
            tol(t) + 0.04 + 0.002
        };
        let lim = budget + 2e-15 * norm(h);
        assert!(dk <= lim, "{t} at {jd}: {dk} km > {lim}");
        worst = worst.max(dk / lim);
    }
    eprintln!("Horizons vectors: worst error / budget = {worst:.3}");
}

/// Timing: `cargo test --release -p fe-core --test pack_ext -- --ignored --nocapture`
#[test]
#[ignore]
fn ext_eval_speed() {
    let Some(p) = open() else { return };
    const N: usize = 1_000_000;
    let mut x: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = |lo: f64, hi: f64| {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        lo + (hi - lo) * ((x >> 11) as f64 / (1u64 << 53) as f64)
    };
    let past: Vec<f64> = (0..N).map(|_| rnd(-3_000_000.0, 1_700_000.0)).collect();
    let main: Vec<f64> = (0..N).map(|_| rnd(1_721_000.0, 2_598_000.0)).collect();
    let seq: Vec<f64> = (0..N).map(|i| -1_000_000.0 + i as f64 * 0.01).collect();
    for (label, t, c) in
        [("Moon geocentric", 301, 399), ("Mars geocentric", 4, 399), ("Sun geocentric", 10, 399)]
    {
        for (tl, times) in [("ext random", &past), ("main random", &main), ("ext sequential", &seq)] {
            let t0 = std::time::Instant::now();
            let mut acc = 0.0;
            for &jd in times.iter() {
                acc += std::hint::black_box(p.position(t, c, std::hint::black_box(jd)).unwrap())[0];
            }
            let ns = t0.elapsed().as_nanos() as f64 / times.len() as f64;
            eprintln!("{label:<16} {tl:<15} position {ns:6.1} ns/call ({acc:.0})");
        }
    }
}
