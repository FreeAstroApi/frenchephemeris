// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Eclipses vs the NASA Five Millennium Canons (Espenak & Meeus) for years
//! 0–2400 and vs the Swiss oracle (black-box), all compared in TT.
//! Fixtures: tests/fixtures/eclipses/gen_eclipses.py. Run with
//! `cargo test -p fe-core --test eclipses -- --nocapture` for the statistics.
//! Skips (passes with a message) when no pack is available.

#[macro_use]
mod common;

use fe_core::context::Engine;
use fe_core::eclipses::*;
use fe_core::ephem::Observer;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::OnceLock;

fn pack_path() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FE_PACK") {
        return Some(PathBuf::from(p)).filter(|p| p.exists());
    }
    let mut dirs = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        dirs.push(PathBuf::from(d).join("packs"));
    }
    dirs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    dirs.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    dirs.iter().map(|d| d.join("fepk-de441-de442-v1.bin")).find(|p| p.exists())
}

fn engine() -> Option<&'static Engine> {
    static E: OnceLock<Option<Engine>> = OnceLock::new();
    E.get_or_init(|| pack_path().map(|p| Engine::open(&p).expect("open pack"))).as_ref()
}

fn fixture() -> &'static Value {
    static F: OnceLock<Value> = OnceLock::new();
    F.get_or_init(|| {
        let p = common::oracle_path("eclipses/eclipses.json");
        serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
    })
}

macro_rules! need_engine {
    () => {
        match engine() {
            Some(e) => e,
            None => {
                eprintln!("SKIP: no FEPK pack (set FE_PACK)");
                return;
            }
        }
    };
}

/// Percentiles of |x|: (p50, p95, max, n).
fn stats(v: &[f64]) -> (f64, f64, f64, usize) {
    let mut a: Vec<f64> = v.iter().map(|x| x.abs()).collect();
    a.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if a.is_empty() {
        return (0.0, 0.0, 0.0, 0);
    }
    let q = |p: f64| a[((a.len() - 1) as f64 * p).round() as usize];
    (q(0.5), q(0.95), a[a.len() - 1], a.len())
}

fn show(name: &str, v: &[f64]) -> (f64, f64, f64, usize) {
    let s = stats(v);
    println!("  {name:<34} p50 {:>10.4}  p95 {:>10.4}  max {:>10.4}  (n={})", s.0, s.1, s.2, s.3);
    s
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

const T0: f64 = 1_721_057.5; // 0000-01-01 (proleptic Julian)
const T1: f64 = 2_598_007.5; // ≈ 2401-01-01

fn solar_letter(k: SolarKind) -> char {
    match k {
        SolarKind::Partial => 'P',
        SolarKind::Annular => 'A',
        SolarKind::Total => 'T',
        SolarKind::Hybrid => 'H',
    }
}

fn lunar_letter(k: LunarKind) -> char {
    match k {
        LunarKind::Penumbral => 'N',
        LunarKind::Partial => 'P',
        LunarKind::Total => 'T',
    }
}

/// Swiss solar type letter from its flags.
fn swiss_solar_letter(fl: i64) -> char {
    if fl & 32 != 0 {
        'H'
    } else if fl & 4 != 0 {
        'T'
    } else if fl & 8 != 0 {
        'A'
    } else {
        'P'
    }
}

fn swiss_lunar_letter(fl: i64) -> char {
    if fl & 4 != 0 {
        'T'
    } else if fl & 16 != 0 {
        'P'
    } else {
        'N'
    }
}

#[test]
fn solar_canon_0_2400() {
    require_oracle!("eclipses/eclipses.json");
    let e = need_engine!();
    let rows = fixture()["solar"].as_array().unwrap();
    let t_start = std::time::Instant::now();
    let mut ours = Vec::new();
    let mut t = T0;
    loop {
        let Ok(ecl) = solar_eclipse_search(e, t, false) else { break };
        if ecl.max > T1 {
            break;
        }
        t = ecl.max;
        ours.push(ecl);
    }
    let elapsed = t_start.elapsed();
    println!(
        "solar: {} found 0–2400 in {:.2?} ({:.1} µs/eclipse)",
        ours.len(),
        elapsed,
        elapsed.as_secs_f64() * 1e6 / ours.len() as f64
    );

    let (mut missed, mut type_nasa, mut type_swiss) = (0, Vec::new(), Vec::new());
    let (mut dt_nasa, mut dg, mut dmag, mut dt_swiss, mut dmag_sw, mut dobs_sw, mut dpos_sw) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    let mut dcont: Vec<Vec<f64>> = vec![vec![]; 6];
    let mut dpos_nasa = vec![];
    let mut used = vec![false; ours.len()];
    let mut j = 0;
    for r in rows {
        let jd = f(&r[0]);
        while j < ours.len() && ours[j].max < jd - 1.0 {
            j += 1;
        }
        let Some(k) = (j..ours.len().min(j + 2)).find(|&k| (ours[k].max - jd).abs() < 0.5) else {
            missed += 1;
            println!("  MISSED NASA solar {} {}", jd, r[2]);
            continue;
        };
        used[k] = true;
        let o = &ours[k];
        let nasa_type = r[2].as_str().unwrap().chars().next().unwrap();
        if solar_letter(o.kind) != nasa_type {
            type_nasa.push(format!(
                "{jd} NASA {} ours {:?} gamma {:.4} mag {:.4}",
                r[2], o.kind, o.gamma, o.magnitude
            ));
        }
        let sw_type = swiss_solar_letter(r[7].as_i64().unwrap());
        if solar_letter(o.kind) != sw_type {
            type_swiss.push(format!("{jd} Swiss {sw_type} ours {:?} NASA {}", o.kind, r[2]));
        }
        dt_nasa.push((o.max - jd) * 86400.0);
        dg.push(o.gamma - f(&r[3]));
        dmag.push(o.magnitude - f(&r[4]));
        let sw_max = jd + f(&r[8]) / 86400.0;
        dt_swiss.push((o.max - sw_max) * 86400.0);
        let c = r[9].as_array().unwrap();
        let mine = [Some(o.p1), Some(o.p4), o.u1, o.u4, o.c1, o.c2];
        for (i, m) in mine.iter().enumerate() {
            if let (Some(m), Some(s)) = (m, c[i].as_f64()) {
                if s != 0.0 || i < 2 {
                    dcont[i].push((m - sw_max) * 86400.0 - s);
                }
            }
        }
        // Earth rotation from Swiss's own Delta T, so only geometry is compared.
        let w = solar_where_ut(e, sw_max, sw_max - f(&r[14]) / 86400.0).unwrap();
        dmag_sw.push(w.attrs.nasa_magnitude - f(&r[12]));
        dobs_sw.push(w.attrs.obscuration.clamp(0.0, 1.0) - f(&r[13]).clamp(0.0, 1.0));
        let dlon = ((w.lon_deg - f(&r[10]) + 540.0) % 360.0 - 180.0) * (w.lat_deg * fe_core::math::DEG).cos();
        dpos_sw.push(dlon.hypot(w.lat_deg - f(&r[11])) * 111.2);
        // NASA greatest-eclipse point (0.1° rounding), with NASA's Delta T.
        let wn = solar_where_ut(e, o.max, o.max - f(&r[1]) / 86400.0).unwrap();
        let dlon =
            ((wn.lon_deg - f(&r[6]) + 540.0) % 360.0 - 180.0) * (wn.lat_deg * fe_core::math::DEG).cos();
        dpos_nasa.push(dlon.hypot(wn.lat_deg - f(&r[5])) * 111.2);
    }
    let spurious: Vec<_> = ours.iter().zip(&used).filter(|(_, u)| !**u).map(|(o, _)| o.max).collect();
    println!(
        "solar vs NASA: {} catalogue, missed {missed}, spurious {} {:?}",
        rows.len(),
        spurious.len(),
        spurious
    );
    println!("  type mismatches vs NASA: {}", type_nasa.len());
    for s in &type_nasa {
        println!("    {s}");
    }
    println!("  type mismatches vs Swiss: {}", type_swiss.len());
    let mut cats = std::collections::BTreeMap::new();
    for s in &type_swiss {
        let w: Vec<&str> = s.split_whitespace().collect();
        *cats.entry(format!("Swiss {} / ours {} / NASA {}", w[2], w[4], w[6])).or_insert(0) += 1;
    }
    for (k, n) in &cats {
        println!("    {n:>3} × {k}");
    }
    let tn = show("max − NASA (s, TT)", &dt_nasa);
    let gs = show("gamma − NASA", &dg);
    let ms = show("magnitude − NASA", &dmag);
    let ts = show("max − Swiss (s, TT)", &dt_swiss);
    for (i, n) in ["P1", "P4", "U1 (tret4)", "U4 (tret5)", "C1 (tret6)", "C2 (tret7)"].iter().enumerate() {
        show(&format!("{n} − Swiss (s)"), &dcont[i]);
    }
    let mw = show("where magnitude − Swiss", &dmag_sw);
    show("where obscuration − Swiss", &dobs_sw);
    let pw = show("where position − Swiss (km)", &dpos_sw);
    show("greatest-eclipse point − NASA (km)", &dpos_nasa);
    assert_eq!(missed, 0);
    assert!(spurious.is_empty());
    assert!(type_nasa.len() <= 3, "{type_nasa:?}");
    assert!(tn.1 < 60.0 && tn.2 < 120.0);
    assert!(gs.2 < 2e-3 && ms.2 < 2e-3);
    assert!(ts.1 < 60.0);
    assert!(mw.1 < 2e-3);
    assert!(pw.1 < 30.0);
}

#[test]
fn lunar_canon_0_2400() {
    require_oracle!("eclipses/eclipses.json");
    let e = need_engine!();
    let rows = fixture()["lunar"].as_array().unwrap();
    let t_start = std::time::Instant::now();
    let mut ours = Vec::new();
    let mut t = T0;
    loop {
        // The pack ends shortly after 2400: the search past the last eclipse may run out.
        let Ok(ecl) = lunar_eclipse_search(e, t, false) else { break };
        if ecl.max > T1 {
            break;
        }
        t = ecl.max;
        ours.push(ecl);
    }
    let elapsed = t_start.elapsed();
    println!(
        "lunar: {} found 0–2400 in {:.2?} ({:.1} µs/eclipse)",
        ours.len(),
        elapsed,
        elapsed.as_secs_f64() * 1e6 / ours.len() as f64
    );

    let (mut missed, mut type_nasa, mut type_swiss) = (0, Vec::new(), Vec::new());
    let (mut dt_nasa, mut dg, mut dum, mut dpen, mut dt_swiss, mut dum_sw, mut dpen_sw) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    let mut dcont: Vec<Vec<f64>> = vec![vec![]; 6];
    let mut used = vec![false; ours.len()];
    let mut j = 0;
    for r in rows {
        let jd = f(&r[0]);
        while j < ours.len() && ours[j].max < jd - 1.0 {
            j += 1;
        }
        let Some(k) = (j..ours.len().min(j + 2)).find(|&k| (ours[k].max - jd).abs() < 0.5) else {
            missed += 1;
            println!("  MISSED NASA lunar {} {} gamma {} pen {}", jd, r[2], r[3], r[4]);
            continue;
        };
        used[k] = true;
        let o = &ours[k];
        let nasa_type = r[2].as_str().unwrap().chars().next().unwrap();
        if lunar_letter(o.kind) != nasa_type {
            type_nasa
                .push(format!("{jd} NASA {} ours {:?} um {:.5} (NASA {})", r[2], o.kind, o.umbral_mag, r[5]));
        }
        let sw_type = swiss_lunar_letter(r[6].as_i64().unwrap());
        if lunar_letter(o.kind) != sw_type {
            type_swiss.push(format!(
                "{jd} Swiss {sw_type} ours {:?} um {:.5} (Swiss {})",
                o.kind, o.umbral_mag, r[9]
            ));
        }
        dt_nasa.push((o.max - jd) * 86400.0);
        dg.push(o.gamma - f(&r[3]));
        dpen.push(o.penumbral_mag - f(&r[4]));
        dum.push(o.umbral_mag - f(&r[5]));
        let sw_max = jd + f(&r[7]) / 86400.0;
        dt_swiss.push((o.max - sw_max) * 86400.0);
        if f(&r[9]) > 0.0 {
            dum_sw.push(o.umbral_mag - f(&r[9]));
        }
        dpen_sw.push(o.penumbral_mag - f(&r[10]));
        // Swiss tret[2..7] = U1, U4, U2, U3, P1, P4.
        let c = r[8].as_array().unwrap();
        let mine = [o.u1, o.u4, o.u2, o.u3, Some(o.p1), Some(o.p4)];
        for (i, m) in mine.iter().enumerate() {
            if let (Some(m), Some(s)) = (m, c[i].as_f64()) {
                dcont[i].push((m - sw_max) * 86400.0 - s);
            }
        }
    }
    let spurious: Vec<_> =
        ours.iter().zip(&used).filter(|(_, u)| !**u).map(|(o, _)| (o.max, o.penumbral_mag)).collect();
    println!(
        "lunar vs NASA: {} catalogue, missed {missed}, spurious {} {:?}",
        rows.len(),
        spurious.len(),
        spurious
    );
    println!("  type mismatches vs NASA: {}", type_nasa.len());
    for s in &type_nasa {
        println!("    {s}");
    }
    println!("  type mismatches vs Swiss: {}", type_swiss.len());
    for s in type_swiss.iter().take(20) {
        println!("    {s}");
    }
    let tn = show("max − NASA (s, TT)", &dt_nasa);
    let gs = show("gamma − NASA", &dg);
    let us = show("umbral mag − NASA", &dum);
    let ps = show("penumbral mag − NASA", &dpen);
    let ts = show("max − Swiss (s, TT)", &dt_swiss);
    let usw = show("umbral mag − Swiss", &dum_sw);
    show("penumbral mag − Swiss", &dpen_sw);
    for (i, n) in ["U1", "U4", "U2", "U3", "P1", "P4"].iter().enumerate() {
        show(&format!("{n} − Swiss (s)"), &dcont[i]);
    }
    assert_eq!(missed, 0);
    assert!(spurious.is_empty());
    assert!(type_nasa.len() <= 3, "{type_nasa:?}");
    assert!(tn.1 < 60.0 && tn.2 < 120.0);
    assert!(gs.2 < 2e-3 && us.2 < 5e-3 && ps.2 < 5e-3);
    assert!(ts.1 < 60.0 && usw.1 < 5e-3);
}

#[test]
fn solar_local_vs_swiss() {
    require_oracle!("eclipses/eclipses.json");
    let e = need_engine!();
    let rows = fixture()["local"].as_array().unwrap();
    let (mut dmax, mut dc1, mut dc4, mut dc23, mut dmag, mut dobs, mut dalt) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    let mut dsun = vec![];
    let mut mismatched = Vec::new();
    for r in rows {
        let obs = Observer { lon_deg: f(&r["lon"]), lat_deg: f(&r["lat"]), height_m: f(&r["h"]) };
        let start = tt(f(&r["start"]));
        let (loc, _glob) = solar_local_search(e, start, &obs, false).unwrap();
        let sw_max = f(&r["max"]);
        if (loc.max - sw_max).abs() > 0.5 {
            mismatched.push(format!("{} start {} swiss {} ours {}", r["city"], r["start"], sw_max, loc.max));
            continue;
        }
        dmax.push((loc.max - sw_max) * 86400.0);
        let c = r["c"].as_array().unwrap();
        dc1.push((loc.c1 - sw_max) * 86400.0 - f(&c[0]));
        dc4.push((loc.c4 - sw_max) * 86400.0 - f(&c[3]));
        for (m, s) in [(loc.c2, &c[1]), (loc.c3, &c[2])] {
            if let (Some(m), Some(s)) = (m, s.as_f64()) {
                dc23.push((m - sw_max) * 86400.0 - s);
            }
        }
        for (m, s) in [(loc.sunrise, &c[4]), (loc.sunset, &c[5])] {
            if let (Some(m), Some(s)) = (m, s.as_f64()) {
                dsun.push((m - sw_max) * 86400.0 - s);
            }
        }
        let a = r["attr"].as_array().unwrap();
        dmag.push(loc.attrs.nasa_magnitude - f(&a[5]));
        dobs.push(loc.attrs.obscuration - f(&a[2]).clamp(0.0, 1.0));
        dalt.push(loc.attrs.true_alt_deg - f(&a[3]));
        let fl = r["flags"].as_i64().unwrap();
        let visible_bits = [(256, loc.max_visible), (512, loc.c1_visible), (4096, loc.c4_visible)];
        for (bit, mine) in visible_bits {
            if (fl & bit != 0) != mine {
                mismatched.push(format!(
                    "{} {} flag {bit}: swiss {} ours {mine}",
                    r["city"],
                    sw_max,
                    fl & bit != 0
                ));
            }
        }
    }
    println!("local solar: {} Swiss cases, {} mismatches", rows.len(), mismatched.len());
    for m in &mismatched {
        println!("    {m}");
    }
    let m = show("local max − Swiss (s)", &dmax);
    let c1 = show("C1 − Swiss (s)", &dc1);
    let c4 = show("C4 − Swiss (s)", &dc4);
    show("C2/C3 − Swiss (s)", &dc23);
    show("sunrise/sunset − Swiss (s)", &dsun);
    let mg = show("magnitude − Swiss", &dmag);
    let ob = show("obscuration − Swiss", &dobs);
    show("Sun true altitude − Swiss (deg)", &dalt);
    assert!(mismatched.len() <= rows.len() / 50, "{mismatched:?}");
    assert!(m.1 < 30.0 && c1.1 < 30.0 && c4.1 < 30.0);
    assert!(mg.1 < 3e-3 && ob.1 < 3e-3);
}

#[test]
fn lunar_local_visibility_vs_swiss() {
    require_oracle!("eclipses/eclipses.json");
    let e = need_engine!();
    let rows = fixture()["lunar_local"].as_array().unwrap();
    let (mut dalt, mut dapp, mut wrong) = (vec![], vec![], 0);
    for r in rows {
        let obs = Observer { lon_deg: f(&r["lon"]), lat_deg: f(&r["lat"]), height_m: 0.0 };
        let mut vis = false;
        for p in r["pts"].as_array().unwrap() {
            let l = lunar_local(e, f(&p[0]), &obs).unwrap();
            dalt.push(l.true_alt_deg - f(&p[1]));
            dapp.push(l.apparent_alt_deg - f(&p[2]));
            vis |= l.apparent_alt_deg > -0.833;
        }
        if vis != r["visible"].as_bool().unwrap() {
            wrong += 1;
            println!("    visibility differs: {} {}", r["city"], r["max"]);
        }
    }
    println!("lunar local: {} city×eclipse cases, {wrong} visibility differences", rows.len());
    let a = show("Moon true altitude − Swiss (deg)", &dalt);
    show("Moon apparent altitude − Swiss (deg)", &dapp);
    assert_eq!(wrong, 0);
    assert!(a.2 < 0.02);
}

#[test]
fn next_eclipse_timing() {
    let e = need_engine!();
    let n = 2000;
    for (name, solar) in [("lunar", false), ("solar", true)] {
        let t0 = tt(2_451_545.0);
        let start = std::time::Instant::now();
        for i in 0..n {
            let t = t0 + i as f64 * 7.3;
            if solar {
                std::hint::black_box(solar_eclipse_search(e, t, false).unwrap());
            } else {
                std::hint::black_box(lunar_eclipse_search(e, t, false).unwrap());
            }
        }
        println!(
            "next {name} eclipse: {:.1} µs/search (this build profile)",
            start.elapsed().as_secs_f64() * 1e6 / n as f64
        );
    }
    let obs = Observer { lon_deg: 2.3522, lat_deg: 48.8566, height_m: 35.0 };
    let start = std::time::Instant::now();
    for i in 0..200 {
        std::hint::black_box(solar_local_search(e, tt(2_451_545.0) + i as f64 * 30.0, &obs, false).unwrap());
    }
    println!(
        "next local solar eclipse (Paris): {:.1} µs/search",
        start.elapsed().as_secs_f64() * 1e6 / 200.0
    );
}

#[test]
fn backward_search_inverts_forward() {
    let e = need_engine!();
    let mut t = tt(2_460_000.0);
    for _ in 0..20 {
        let a = lunar_eclipse_search(e, t, false).unwrap();
        let b = lunar_eclipse_search(e, a.max + 1.0, true).unwrap();
        assert!((a.max - b.max).abs() < 1e-9);
        let c = solar_eclipse_search(e, t, false).unwrap();
        let d = solar_eclipse_search(e, c.max + 1.0, true).unwrap();
        assert!((c.max - d.max).abs() < 1e-9);
        t = a.max.min(c.max) + 1.0;
    }
}

#[test]
fn besselian_vs_nasa_elements() {
    require_oracle!("eclipses/eclipses.json");
    let e = need_engine!();
    let rows = fixture()["besselian"].as_array().unwrap();
    let (mut dx, mut dd, mut dl1, mut dl2, mut dmu, mut dtf) =
        (vec![], vec![], vec![], vec![], vec![], vec![]);
    let mut dge = vec![];
    let (mut c_nasa, mut c_swiss_nasa, mut c_ours_swiss) =
        (vec![vec![]; 3], vec![vec![]; 3], vec![vec![]; 3]);
    for r in rows {
        let t0 = f(&r["t0"]);
        let dt = f(&r["dt"]) / 86400.0;
        let ev = |k: &str, h: f64| -> f64 {
            r["co"][k].as_array().unwrap().iter().enumerate().map(|(i, c)| f(c) * h.powi(i as i32)).sum()
        };
        for h in [-2.0, -1.0, 0.0, 1.0, 2.0] {
            let t = t0 + h / 24.0;
            // NASA tabulates the *ephemeris* hour angle μ (sidereal time at TT).
            let b = besselian_elements(e, t, t).unwrap();
            dx.push((b.x - ev("x", h)).hypot(b.y - ev("y", h)));
            dd.push((b.d_deg - ev("d", h)) * 3600.0);
            dl1.push(b.l1 - ev("l1", h));
            dl2.push(b.l2 - ev("l2", h));
            dmu.push(((b.mu_deg - ev("mu", h) + 540.0) % 360.0 - 180.0) * 3600.0);
        }
        let b = besselian_elements(e, t0, t0).unwrap();
        dtf.push(b.tan_f1 - f(&r["tan_f1"]));
        dtf.push(b.tan_f2 - f(&r["tan_f2"]));
        let o = solar_eclipse_search(e, f(&r["max"]) - 1.0, false).unwrap();
        if let (Some(la), Some(lo)) = (r["ge_lat"].as_f64(), r["ge_lon"].as_f64()) {
            // NASA greatest-eclipse point (0.1′), at NASA's instant and Delta T.
            let m = f(&r["max"]);
            let w = solar_where_ut(e, m, m - dt).unwrap();
            let dlon = ((w.lon_deg - lo + 540.0) % 360.0 - 180.0) * (la * fe_core::math::DEG).cos();
            dge.push(dlon.hypot(w.lat_deg - la) * 111.2);
        }
        let sw = r["swiss"].as_array().unwrap();
        let mine = [[Some(o.p1), Some(o.p4)], [o.u1, o.u4], [o.c1, o.c2]];
        for (i, key) in ["p", "u", "c"].iter().enumerate() {
            for j in 0..2 {
                let (Some(n), Some(m)) = (r[*key][j].as_f64(), mine[i][j]) else { continue };
                c_nasa[i].push((m - n) * 86400.0);
                if let Some(s) = sw[2 + 2 * i + j].as_f64() {
                    c_swiss_nasa[i].push((s - n) * 86400.0);
                    c_ours_swiss[i].push((m - s) * 86400.0);
                }
            }
        }
    }
    println!("Besselian elements vs NASA (VSOP87/ELP2000-85), {} eclipses × 5 epochs:", rows.len());
    let x = show("|(x, y)| (Earth radii)", &dx);
    show("d (arcsec)", &dd);
    let l1 = show("l1", &dl1);
    show("l2", &dl2);
    let mu = show("mu (arcsec, ephemeris meridian)", &dmu);
    show("tan f1, tan f2", &dtf);
    let ge = show("greatest-eclipse point (km)", &dge);
    for (i, n) in ["P1/P4", "U1/U4", "C1/C2"].iter().enumerate() {
        show(&format!("{n} ours − NASA-elements (s)"), &c_nasa[i]);
        show(&format!("{n} Swiss − NASA-elements (s)"), &c_swiss_nasa[i]);
    }
    assert!(x.2 < 2e-4 && l1.2 < 1e-4 && mu.2 < 30.0);
    assert!(c_nasa.iter().flatten().all(|d| d.abs() < 5.0));
    assert!(ge.2 < 3.0);
}
