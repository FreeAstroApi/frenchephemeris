// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Eclipse visibility at a place: the classic `*_eclipse_how` flags (nothing is
//! reported while the body is below the apparent horizon), local lunar eclipse
//! searches at polar sites, and the NONCENTRAL bit of local solar eclipses.
//!
//! Public tests check physical consistency; the oracle comparison
//! (`eclipses/visibility.json`, gen_eclipse_visibility.py) needs the private corpus.

#[macro_use]
mod common;

use fe_core::context::{Body, Engine};
use fe_core::eclipses::*;
use fe_core::ephem::Observer;
use fe_core::events::{DiscPoint, RiseSetOptions};
use fe_core::lunar_local::lunar_eclipse_local_search;
use fe_core::time::{julday, Calendar};
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

fn obs(lon: f64, lat: f64, h: f64) -> Observer {
    Observer { lon_deg: lon, lat_deg: lat, height_m: h }
}

fn jd(y: i32, m: u32, d: u32) -> f64 {
    julday(y, m, d, 0.0, Calendar::Gregorian)
}

#[test]
fn lunar_how_needs_the_moon_above_the_horizon() {
    let e = need_engine!();
    let ecl = lunar_eclipse_search(e, e.ut_to_tt(jd(2018, 7, 1)), false).unwrap();
    assert_eq!(ecl.kind, LunarKind::Total);
    // New York: the Moon is ~42° below the horizon at greatest eclipse.
    let ny = lunar_local(e, ecl.max, &obs(-74.006, 40.7128, 10.0)).unwrap();
    assert!(ny.apparent_alt_deg < -30.0 && ny.shadow.umbral_mag > 1.0);
    assert_eq!(lunar_seen_kind(&ny), None);
    // Nairobi: Moon high in the sky.
    let nb = lunar_local(e, ecl.max, &obs(36.82, -1.29, 1700.0)).unwrap();
    assert!(nb.apparent_alt_deg > 30.0);
    assert_eq!(lunar_seen_kind(&nb), Some(LunarKind::Total));
    // No eclipse three days later: nothing seen, whatever the altitude.
    let none = lunar_local(e, ecl.max + 3.0, &obs(36.82, -1.29, 1700.0)).unwrap();
    assert!(none.shadow.penumbral_mag <= 0.0);
    assert_eq!(lunar_seen_kind(&none), None);
}

#[test]
fn solar_how_needs_the_sun_above_the_horizon() {
    let e = need_engine!();
    let ecl = solar_eclipse_search(e, e.ut_to_tt(jd(2024, 1, 1)), false).unwrap();
    let w = solar_where(e, ecl.max).unwrap();
    let on_path = solar_how(e, ecl.max, &obs(w.lon_deg, w.lat_deg, 0.0)).unwrap();
    assert_eq!(solar_seen_kind(&on_path), Some(SolarKind::Total));
    assert!(on_path.axis_on_earth);
    // Delhi is in the night: no eclipse seen.
    let night = solar_how(e, ecl.max, &obs(77.209, 28.6139, 216.0)).unwrap();
    assert!(night.apparent_alt_deg < 0.0);
    assert_eq!(solar_seen_kind(&night), None);
    // Sydney in daylight five days later: no eclipse.
    let none = solar_how(e, ecl.max + 5.0, &obs(-104.0, 25.0, 0.0)).unwrap();
    assert!(none.fraction_diameter <= 0.0);
    assert_eq!(solar_seen_kind(&none), None);
}

/// Lower-limb altitude above the rise/set horizon (deg), the visibility rule
/// of the local lunar eclipse search.
fn moon_up(e: &Engine, jd_ut: f64, o: &Observer) -> f64 {
    let opt = RiseSetOptions { disc: DiscPoint::LowerLimb, ..Default::default() };
    e.rise_set_altitude(Body::Moon, jd_ut, o, &opt).unwrap().to_degrees()
}

#[test]
fn polar_local_lunar_eclipses_are_really_visible() {
    let e = need_engine!();
    let sites = [
        ("North Pole 89.99N", obs(0.0, 89.99, 0.0)),
        ("South Pole 89.99S", obs(0.0, -89.99, 2835.0)),
        ("Longyearbyen", obs(15.6267, 78.2232, 10.0)),
        ("Alert", obs(-62.3481, 82.5018, 30.0)),
        ("Utqiagvik", obs(-156.7887, 71.2906, 10.0)),
        ("Murmansk", obs(33.0827, 68.9585, 50.0)),
        ("McMurdo", obs(166.676, -77.846, 10.0)),
    ];
    let mut n = 0;
    for (name, o) in sites {
        let mut t = jd(2000, 1, 1);
        for _ in 0..12 {
            let l = lunar_eclipse_local_search(e, t, &o, false).unwrap();
            let (p1, p4) = (e.tt_to_ut(l.global.p1), e.tt_to_ut(l.global.p4));
            // Some instant of the eclipse has the Moon up: the maximum when it is
            // flagged visible, else a moonrise/moonset inside the eclipse.
            if l.max_visible {
                assert!(moon_up(e, l.max, &o) > -1e-6, "{name} {}: max flagged visible, Moon below", l.max);
            } else {
                let edge = l.moonrise.or(l.moonset);
                assert!(edge.is_some_and(|x| x >= p1 && x <= p4), "{name} {}: no visible part", l.max);
            }
            for c in [l.p1, l.u1, l.u2, l.u3, l.u4, l.p4].into_iter().flatten() {
                assert!(moon_up(e, c, &o) > -1e-4, "{name}: contact {c} flagged visible below the horizon");
            }
            n += 1;
            t = l.global.max + 1.0;
        }
    }
    assert_eq!(n, 84);
    // The 2026-03-03 total eclipse is visible from 89.99° N (Moon ~6° up).
    let l = lunar_eclipse_local_search(e, jd(2026, 1, 1), &obs(0.0, 89.99, 0.0), false).unwrap();
    assert!((l.max - jd(2026, 3, 3)).abs() < 1.0, "found {}", l.max);
    assert_eq!(l.kind, LunarKind::Total);
}

#[test]
fn local_partial_solar_eclipse_is_noncentral_when_the_axis_misses_the_earth() {
    let e = need_engine!();
    // 2025-03-29: a partial eclipse, the axis never meets the Earth.
    let g = solar_eclipse_search(e, e.ut_to_tt(jd(2025, 3, 1)), false).unwrap();
    assert!(g.kind == SolarKind::Partial && !g.central);
    let l = solar_local(e, &g, &obs(2.3522, 48.8566, 35.0)).unwrap().unwrap();
    assert!(!l.attrs.axis_on_earth);
    // 2024-04-08 from Dallas: inside the umbra, the axis is on the Earth.
    let g = solar_eclipse_search(e, e.ut_to_tt(jd(2024, 4, 1)), false).unwrap();
    let l = solar_local(e, &g, &obs(-96.797, 32.7767, 140.0)).unwrap().unwrap();
    assert!(l.kind == SolarKind::Total && l.attrs.axis_on_earth);
}

// ---------------------------------------------------------------------------------------------
// Oracle comparison (private corpus)

const FIX: &str = "eclipses/visibility.json";

fn lunar_bits(k: Option<LunarKind>) -> i64 {
    k.map_or(0, |k| match k {
        LunarKind::Total => 4,
        LunarKind::Partial => 16,
        LunarKind::Penumbral => 64,
    })
}

fn solar_bits(a: &SolarAttrs) -> i64 {
    match solar_seen_kind(a) {
        None => 0,
        Some(k) => {
            (match k {
                SolarKind::Total => 4,
                SolarKind::Annular => 8,
                SolarKind::Partial => 16,
                SolarKind::Hybrid => 32,
            }) | 128
                | if a.axis_on_earth { 1 } else { 2 }
        }
    }
}

#[test]
fn how_flags_vs_oracle() {
    require_oracle!(FIX);
    let e = need_engine!();
    let fx = common::oracle_json(FIX).unwrap();
    let (mut n, mut bad) = (0, vec![]);
    for r in fx["how"].as_array().unwrap() {
        let f = |k: &str| r[k].as_f64().unwrap();
        let o = obs(f("lon"), f("lat"), f("h"));
        // The oracle's instant (UT) with this engine's Delta T.
        let t = e.ut_to_tt(f("ut"));
        let mine = if r["kind"] == "lunar" {
            lunar_bits(lunar_seen_kind(&lunar_local(e, t, &o).unwrap()))
        } else {
            solar_bits(&solar_how(e, t, &o).unwrap())
        };
        n += 1;
        if mine != r["flags"].as_i64().unwrap() {
            bad.push(format!(
                "{} {} {}: oracle {} ours {mine} (app alt {:.3})",
                r["kind"],
                r["place"],
                f("ut"),
                r["flags"],
                f("app_alt")
            ));
        }
    }
    println!("eclipse_how flags: {} of {n} differ", bad.len());
    for b in &bad {
        println!("    {b}");
    }
    assert!(bad.len() * 100 <= n, "{bad:?}");
}

#[test]
fn local_searches_vs_oracle() {
    require_oracle!(FIX);
    let e = need_engine!();
    let fx = common::oracle_json(FIX).unwrap();
    let (mut n, mut other_event, mut nc_diff) = (0, vec![], 0);
    for r in fx["local"].as_array().unwrap() {
        let f = |k: &str| r[k].as_f64().unwrap();
        let o = obs(f("lon"), f("lat"), f("h"));
        n += 1;
        if r["kind"] == "lunar" {
            let l = lunar_eclipse_local_search(e, f("start"), &o, false).unwrap();
            if (l.max - f("max")).abs() > 0.3 {
                other_event.push(format!(
                    "lunar {} start {}: oracle {} ours {}",
                    r["place"],
                    f("start"),
                    f("max"),
                    l.max
                ));
            }
        } else {
            let (l, _) = solar_local_search(e, e.ut_to_tt(f("start")), &o, false).unwrap();
            if (e.tt_to_ut(l.max) - f("max")).abs() > 0.3 {
                other_event.push(format!(
                    "solar {} start {}: oracle {} ours {}",
                    r["place"],
                    f("start"),
                    f("max"),
                    l.max
                ));
            } else if (r["flags"].as_i64().unwrap() & 2 != 0) == l.attrs.axis_on_earth {
                nc_diff += 1;
            }
        }
    }
    println!("local searches: {n}, different event {}, NONCENTRAL bit differs {nc_diff}", other_event.len());
    for b in &other_event {
        println!("    {b}");
    }
    assert!(other_event.len() * 50 <= n && nc_diff * 50 <= n);
}
