// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Lunar occultations vs the Swiss oracle (black box) and vs IOTA's published
//! graze-limit predictions. Fixtures: tests/fixtures/occultations/gen_occultations.py.
//! `cargo test -p fe-core --test occultations -- --nocapture` prints the statistics.

#[macro_use]
mod common;

use fe_core::context::{Body, Engine};
use fe_core::eclipses::{tt, SolarKind};
use fe_core::ephem::Observer;
use fe_core::occultations::*;
use fe_core::stars;
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

/// Full fixture including the Swiss-oracle part (private corpus; see tests/common).
fn oracle() -> &'static Value {
    static O: OnceLock<Value> = OnceLock::new();
    O.get_or_init(|| serde_json::from_str(&common::read_oracle("occultations/occultations.json")).unwrap())
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
    println!("  {name:<40} p50 {:>9.3}  p95 {:>9.3}  max {:>9.3}  (n={})", s.0, s.1, s.2, s.3);
    s
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn obj(v: &Value) -> Occulted {
    match v {
        Value::String(s) => Occulted::Star(stars::lookup(s).unwrap()),
        n => Occulted::Body(Body::from_id(n.as_i64().unwrap() as i32).unwrap()),
    }
}

fn kind_flags(k: SolarKind, central: bool) -> i64 {
    (match k {
        SolarKind::Total => 4,
        SolarKind::Annular => 8,
        SolarKind::Partial => 16,
        SolarKind::Hybrid => 32,
    }) | if central { 1 } else { 2 }
}

fn found<T>(s: Search<T>) -> Option<T> {
    match s {
        Search::Found(x) => Some(x),
        Search::NoneAt(_) => None,
    }
}

#[test]
fn global_vs_swiss() {
    require_oracle!("occultations/occultations.json");
    let e = need_engine!();
    let rows = oracle()["glob"].as_array().unwrap();
    // Group the oracle's chained events per body and chain ours over the same span.
    let mut bodies: Vec<String> = rows.iter().map(|r| r["body"].to_string()).collect();
    bodies.dedup();
    let (mut missed, mut spurious, mut flag_diff) = (0, 0, vec![]);
    let (mut dmax, mut dp, mut du, mut dc, mut dpos, mut dlim, mut dmag) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    let start = std::time::Instant::now();
    let mut n_ours = 0;
    for b in &bodies {
        let rs: Vec<&Value> = rows.iter().filter(|r| &r["body"].to_string() == b).collect();
        let o = obj(&rs[0]["body"]);
        let end = f(&rs[rs.len() - 1]["tret"][0]) + 1.0;
        let mut ours = vec![];
        let mut t = tt(f(&rs[0]["start"]));
        loop {
            let ev = found(occultation_search(e, &o, t, false, false).unwrap()).unwrap();
            if ev.max > end {
                break;
            }
            t = ev.max;
            ours.push(ev);
        }
        n_ours += ours.len();
        let mut used = vec![false; ours.len()];
        for r in &rs {
            let tr: Vec<f64> = r["tret"].as_array().unwrap().iter().map(f).collect();
            let Some(k) = ours.iter().position(|x| (x.max - tr[0]).abs() < 0.3) else {
                missed += 1;
                println!("    MISSED {b} {}", tr[0]);
                continue;
            };
            used[k] = true;
            let ev = &ours[k];
            let fl = r["flags"].as_i64().unwrap();
            if kind_flags(ev.kind, ev.central) != fl {
                flag_diff.push(format!("{b} {} swiss {fl} ours {}", tr[0], kind_flags(ev.kind, ev.central)));
            }
            dmax.push((ev.max - tr[0]) * 86400.0);
            for (i, (m, s)) in [
                (Some(ev.p1), tr[2]),
                (Some(ev.p4), tr[3]),
                (ev.u1, tr[4]),
                (ev.u4, tr[5]),
                (ev.c1, tr[6]),
                (ev.c2, tr[7]),
            ]
            .into_iter()
            .enumerate()
            {
                if let (Some(m), true) = (m, s != 0.0) {
                    [&mut dp, &mut du, &mut dc][i / 2].push((m - s) * 86400.0);
                }
            }
            let w = occultation_where_ut(e, &o, tr[0], tr[0] - f(&r["dt"]) / 86400.0).unwrap();
            let gp: Vec<f64> = r["geopos"].as_array().unwrap().iter().map(f).collect();
            let km = |lon: f64, lat: f64, slon: f64, slat: f64| {
                let dl = ((lon - slon + 540.0) % 360.0 - 180.0) * (lat * fe_core::math::DEG).cos();
                dl.hypot(lat - slat) * 111.2
            };
            dpos.push(km(w.lon_deg, w.lat_deg, gp[0], gp[1]));
            for (i, l) in w.limits.iter().enumerate() {
                if let (Some((lo, la)), true) = (l, gp[2 + 2 * i] != 0.0 || gp[3 + 2 * i] != 0.0) {
                    dlim.push(km(*lo, *la, gp[2 + 2 * i], gp[3 + 2 * i]));
                }
            }
            let at: Vec<f64> = r["attr"].as_array().unwrap().iter().map(f).collect();
            if !o.is_star() {
                dmag.push(w.attrs.ratio / at[1] - 1.0);
            }
        }
        for (k, u) in used.iter().enumerate() {
            if !u {
                spurious += 1;
                println!(
                    "    SPURIOUS {b} {} kind {:?} gamma {:.4}",
                    ours[k].max, ours[k].kind, ours[k].gamma
                );
            }
        }
    }
    let el = start.elapsed();
    println!("global occultations: {} oracle events, {n_ours} ours ({:.1} µs/search); missed {missed}, spurious {spurious}, flag differences {}",
        rows.len(), el.as_secs_f64() * 1e6 / n_ours as f64, flag_diff.len());
    for s in flag_diff.iter().take(15) {
        println!("    {s}");
    }
    let m = show("max − Swiss (s, TT)", &dmax);
    show("P1/P4 (tret 2,3) − Swiss (s)", &dp);
    show("U1/U4 (tret 4,5) − Swiss (s)", &du);
    show("C1/C2 (tret 6,7) − Swiss (s)", &dc);
    show("where point − Swiss (km, Swiss ΔT)", &dpos);
    show("where limits − Swiss (km)", &dlim);
    show("planet Moon/planet ratio / Swiss − 1", &dmag);
    assert!(missed <= rows.len() / 100 && spurious <= rows.len() / 100);
    assert!(m.1 < 30.0);
}

#[test]
fn local_vs_swiss() {
    require_oracle!("occultations/occultations.json");
    let e = need_engine!();
    let rows = oracle()["local"].as_array().unwrap();
    let (mut dmax, mut dc14, mut dc23, mut drs) = (vec![], vec![], vec![], vec![]);
    let (mut mism, mut flagd) = (vec![], 0);
    let mut bits = std::collections::BTreeMap::new();
    let t_start = std::time::Instant::now();
    for r in rows {
        let o = obj(&r["body"]);
        let obs = Observer { lon_deg: f(&r["lon"]), lat_deg: f(&r["lat"]), height_m: f(&r["h"]) };
        let tr: Vec<f64> = r["tret"].as_array().unwrap().iter().map(f).collect();
        let Some((l, _)) =
            found(occultation_local_search(e, &o, tt(f(&r["start"])), &obs, false, false).unwrap())
        else {
            mism.push(format!("{} {} none", r["city"], r["body"]));
            continue;
        };
        if (l.max - tr[0]).abs() > 0.3 {
            mism.push(format!("{} {} swiss {} ours {}", r["city"], r["body"], tr[0], l.max));
            continue;
        }
        dmax.push((l.max - tr[0]) * 86400.0);
        dc14.push((l.c1 - tr[1]) * 86400.0);
        dc14.push((l.c4 - tr[4]) * 86400.0);
        for (m, s) in [(l.c2, tr[2]), (l.c3, tr[3])] {
            if let (Some(m), true) = (m, s != 0.0) {
                dc23.push((m - s) * 86400.0);
            }
        }
        for (m, s) in [(l.rise, tr[5]), (l.set, tr[6])] {
            if let (Some(m), true) = (m, s != 0.0) {
                drs.push((m - s) * 86400.0);
            }
        }
        let fl = r["flags"].as_i64().unwrap();
        let ours = [
            (4, l.kind == SolarKind::Total),
            (16, l.kind == SolarKind::Partial),
            (128, l.visible),
            (256, l.max_visible),
            (512, l.c1_visible),
            (1024, l.c2_visible),
            (2048, l.c3_visible),
            (4096, l.c4_visible),
            (8192, l.c1_daylight),
            (16384, l.c4_daylight),
        ];
        let mut bad = false;
        for (bit, mine) in ours {
            if (fl & bit != 0) != mine {
                *bits.entry(bit).or_insert(0) += 1;
                bad = true;
            }
        }
        flagd += bad as usize;
    }
    println!("local occultations: {} oracle events, {} not matched, {flagd} with flag differences {:?} ({:.0} µs/search)",
        rows.len(), mism.len(), bits, t_start.elapsed().as_secs_f64() * 1e6 / rows.len() as f64);
    for m in mism.iter().take(15) {
        println!("    {m}");
    }
    let m = show("local max − Swiss (s)", &dmax);
    let c = show("C1/C4 − Swiss (s)", &dc14);
    show("C2/C3 − Swiss (s)", &dc23);
    show("body rise/set in [C1,C4] − Swiss (s)", &drs);
    assert!(mism.len() <= rows.len() / 50, "{mism:?}");
    assert!(m.1 < 10.0 && c.1 < 15.0);
}

#[test]
fn one_try_and_backward_vs_swiss() {
    require_oracle!("occultations/occultations.json");
    let e = need_engine!();
    let (mut same, mut n, mut dt, mut dconj) = (0, 0, vec![], vec![]);
    for r in oracle()["one_try"].as_array().unwrap() {
        let o = obj(&r["body"]);
        let swiss_found = r["flags"].as_i64().unwrap() != 0;
        let res = occultation_search(e, &o, tt(f(&r["start"])), false, true).unwrap();
        let (mine, t0) = match res {
            Search::Found(ev) => (true, ev.max),
            Search::NoneAt(c) => (false, c),
        };
        n += 1;
        if mine == swiss_found {
            same += 1;
            // Found: maximum; not found: the tested conjunction.
            if mine { &mut dt } else { &mut dconj }.push((t0 - f(&r["t0"])) * 86400.0);
        } else {
            println!(
                "    ONE_TRY differs: {} start {} swiss {} ours {mine}",
                r["body"], r["start"], swiss_found
            );
        }
    }
    println!("ECL_ONE_TRY: {same}/{n} same found/not-found outcome");
    show("ONE_TRY found: max − Swiss (s)", &dt);
    show("ONE_TRY none: conjunction − Swiss tret[0] (s)", &dconj);
    let (mut bsame, mut bn) = (0, 0);
    for r in oracle()["backward"].as_array().unwrap() {
        let o = obj(&r["body"]);
        let ev = found(occultation_search(e, &o, tt(f(&r["start"])), true, false).unwrap()).unwrap();
        bn += 1;
        if (ev.max - f(&r["t0"])).abs() < 0.01 {
            bsame += 1;
        } else {
            println!(
                "    backward differs: {} start {} swiss {} ours {}",
                r["body"], r["start"], r["t0"], ev.max
            );
        }
    }
    println!("backward search: {bsame}/{bn} same event");
    assert!(same * 100 >= n * 95 && bsame * 100 >= bn * 95);
}

#[test]
fn graze_limits_vs_iota() {
    require_oracle!("occultations/occultations.json");
    let e = need_engine!();
    let mut by: std::collections::BTreeMap<String, (Vec<f64>, Vec<f64>)> = Default::default();
    for r in oracle()["iota"].as_array().unwrap() {
        let o = obj(&r["body"]);
        let obs = Observer { lon_deg: f(&r["lon"]), lat_deg: f(&r["lat"]), height_m: 0.0 };
        let t = tt(f(&r["jd_ut"]));
        let (tmin, sep, sm, ss) = closest_approach(e, &o, t, &obs).unwrap();
        let edge = if r["limit"] == "inner" { sm - ss } else { sm + ss };
        let key = format!(
            "{} {} {}",
            r["file"].as_str().unwrap(),
            r["side"].as_str().unwrap(),
            r["limit"].as_str().unwrap()
        );
        let ent = by.entry(key).or_default();
        // Gap on the sky (arcsec) between the limb and the body edge at closest approach.
        ent.0.push((sep - edge) * 3600.0);
        ent.1.push((tmin - t) * 86400.0);
    }
    println!("graze limits vs IOTA (site on the predicted limit, smooth Moon, sea level):");
    let mut worst: f64 = 0.0;
    for (k, (gap, dt)) in &by {
        let g = show(&format!("{k}: limb gap (arcsec)"), gap);
        show(&format!("{k}: graze time − IOTA (s)"), dt);
        worst = worst.max(g.1);
    }
    // 1″ of limb gap ≈ 1.9 km across the limit line.
    assert!(worst < 1.0);
}
