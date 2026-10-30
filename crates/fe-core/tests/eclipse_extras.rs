// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Eclipse extras: saros series/member numbers (vs NASA's canon and Swiss),
//! "eclipse at local apparent noon" (vs Swiss tret[1]) and the umbra/penumbra
//! limits of `sol_eclipse_where` (vs NASA path tables and Swiss geopos[2..9]).
//! Fixtures: tests/fixtures/eclipses/gen_eclipse_extras.py.

#[macro_use]
mod common;

use fe_core::context::Engine;
use fe_core::eclipses::*;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::OnceLock;

fn engine() -> Option<&'static Engine> {
    static E: OnceLock<Option<Engine>> = OnceLock::new();
    E.get_or_init(|| {
        let mut c = vec![];
        if let Ok(p) = std::env::var("FE_PACK") {
            c.push(PathBuf::from(p));
        }
        c.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs/fepk-de441-de442-v1.bin"));
        c.push(
            PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join("fepk-de441-de442-v1.bin"),
        );
        c.into_iter().find(|p| p.exists()).map(|p| Engine::open(&p).expect("open pack"))
    })
    .as_ref()
}

fn fixture() -> &'static Value {
    static F: OnceLock<Value> = OnceLock::new();
    F.get_or_init(|| {
        let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/eclipses/extras.json");
        serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
    })
}

/// Full fixture including the Swiss-oracle part (private corpus; see tests/common).
fn oracle() -> &'static Value {
    static O: OnceLock<Value> = OnceLock::new();
    O.get_or_init(|| serde_json::from_str(&common::read_oracle("eclipses/extras.json")).unwrap())
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

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
    println!("  {name:<44} p50 {:>9.3}  p95 {:>9.3}  max {:>9.3}  (n={})", s.0, s.1, s.2, s.3);
    s
}

fn km(lon: f64, lat: f64, slon: f64, slat: f64) -> f64 {
    let dl = ((lon - slon + 540.0) % 360.0 - 180.0) * (lat * fe_core::math::DEG).cos();
    dl.hypot(lat - slat) * 111.2
}

#[test]
fn saros_vs_nasa_canon() {
    let fx = fixture();
    for (key, solar) in [("nasa_solar_saros", true), ("nasa_lunar_saros", false)] {
        let (mut ok, mut bad, mut none) = (0, 0, 0);
        for r in fx[key].as_array().unwrap() {
            match saros(f(&r[0]), solar) {
                Some((s, _)) if s as i64 == r[1].as_i64().unwrap() => ok += 1,
                Some(_) => bad += 1,
                None => none += 1,
            }
        }
        println!("{key}: {ok} equal, {bad} different, {none} not in a catalogued series");
        assert_eq!(bad + none, 0);
    }
}

#[test]
fn saros_members_vs_swiss() {
    require_oracle!("eclipses/extras.json");
    let fx = oracle();
    let (mut n, mut same) = (0, 0);
    for (key, solar) in [("swiss_solar", true), ("swiss_lunar", false)] {
        for r in fx[key].as_array().unwrap() {
            let sw = if solar { &r["where"][0]["saros"] } else { &r["saros"] };
            let ours = saros(f(&r["max"]), solar).map(|(s, m)| (s as f64, m as f64));
            n += 1;
            if ours == Some((f(&sw[0]), f(&sw[1]))) {
                same += 1;
            } else {
                println!("    {key} {} ours {ours:?} swiss {sw}", r["max"]);
            }
        }
    }
    println!("saros series and member vs Swiss attr[9..10]: {same}/{n} identical");
    assert!(same * 100 >= n * 99);
}

#[test]
fn local_noon_and_limits() {
    require_oracle!("eclipses/extras.json");
    let Some(e) = engine() else {
        eprintln!("SKIP: no pack");
        return;
    };
    let fx = oracle();
    let (mut dnoon, mut dpos, mut dlim) = (vec![], vec![], [vec![], vec![], vec![], vec![]]);
    let (mut noon_presence, mut n, mut lim_presence) = (0, 0, 0);
    for r in fx["swiss_solar"].as_array().unwrap() {
        let ev = solar_eclipse_search(e, f(&r["max"]) - 1.0, false).unwrap();
        let noon = solar_local_noon(e, &ev).unwrap();
        let sw = f(&r["noon"]);
        n += 1;
        if noon.is_some() == (sw != 0.0) {
            noon_presence += 1;
        } else {
            println!("    noon presence differs at {}: ours {noon:?} swiss {sw}", r["max"]);
        }
        if let (Some(t), true) = (noon, sw != 0.0) {
            dnoon.push((t - sw) * 86400.0);
        }
        for w in r["where"].as_array().unwrap() {
            let ut = f(&w["ut"]);
            let tt = ut + f(&w["dt"]) / 86400.0;
            let ours = solar_where_ut(e, tt, ut).unwrap();
            let gp: Vec<f64> = w["geopos"].as_array().unwrap().iter().map(f).collect();
            if w["flags"].as_i64().unwrap() != 0 {
                dpos.push(km(ours.lon_deg, ours.lat_deg, gp[0], gp[1]));
            }
            for i in 0..4 {
                let sw_has = gp[2 + 2 * i] != 0.0 || gp[3 + 2 * i] != 0.0;
                if ours.limits[i].is_some() == sw_has {
                    lim_presence += 1;
                }
                if let (Some((lo, la)), true) = (ours.limits[i], sw_has) {
                    dlim[i].push(km(lo, la, gp[2 + 2 * i], gp[3 + 2 * i]));
                }
            }
        }
    }
    println!("local apparent noon: presence equal in {noon_presence}/{n}");
    let nn = show("tret[1] noon − Swiss (s)", &dnoon);
    show("where point − Swiss (km, Swiss ΔT)", &dpos);
    println!("limits present/absent like Swiss: {lim_presence}/{}", 3 * n * 4);
    for (i, name) in ["umbra N", "umbra S", "penumbra N", "penumbra S"].iter().enumerate() {
        show(&format!("{name} limit − Swiss (km)"), &dlim[i]);
    }
    // NASA path tables: umbral limits at fixed UT (NASA ΔT from its Besselian pages).
    let nasa_dt = |file: &str| match &file[2..12] {
        "2017Aug21T" => 70.0,
        "2023Oct14A" => 70.6,
        _ => 70.6,
    };
    let (mut dn, mut ds) = (vec![], vec![]);
    for r in fx["nasa_paths"].as_array().unwrap() {
        let ut = f(&r["jd_ut"]);
        let tt = ut + nasa_dt(r["file"].as_str().unwrap()) / 86400.0;
        let w = solar_where_ut(e, tt, ut).unwrap();
        let (n0, s0) = (&r["north"], &r["south"]);
        if let Some((lo, la)) = w.limits[0] {
            dn.push(km(lo, la, f(&n0[0]), f(&n0[1])));
        }
        if let Some((lo, la)) = w.limits[1] {
            ds.push(km(lo, la, f(&s0[0]), f(&s0[1])));
        }
    }
    let a = show("umbra N limit − NASA path table (km)", &dn);
    let b = show("umbra S limit − NASA path table (km)", &ds);
    assert!(nn.1 < 120.0);
    assert!(a.1 < 5.0 && b.1 < 5.0);
}
