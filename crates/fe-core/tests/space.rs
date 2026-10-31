// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Observers away from Earth (feature `space`) against JPL Horizons
//! (fixtures/space/gen_space.py), plus geometric checks of houses on every world.
#![cfg(feature = "space")]

use fe_core::context::{Body, Engine};
use fe_core::ephem::{apparent, Corrections, LightTime};
use fe_core::math::{cross, dot, mat_t_vec, mat_vec, norm, to_spherical, DEG};
use fe_core::pack::Pack;
use fe_core::space::{Site, SpaceObserver, World, Zodiac};
use serde_json::Value;
use std::path::PathBuf;

fn engine() -> Option<Engine> {
    let name = "fepk-de441-de442-v1.bin";
    let path = [
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs").join(name),
        PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join(name),
    ]
    .into_iter()
    .find(|p| p.exists());
    match path {
        Some(p) => Some(Engine::new(Pack::open(&p).expect("open pack"))),
        None => {
            eprintln!("SKIP: no FEPK pack found (set FE_PACK_DIR)");
            None
        }
    }
}

fn fixture() -> Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/space/space.json");
    serde_json::from_str(&std::fs::read_to_string(p).expect("space fixture")).unwrap()
}

fn body(name: &str) -> Body {
    match name {
        "sun" => Body::Sun,
        "earth" => Body::Earth,
        "moon" => Body::Moon,
        "venus" => Body::Venus,
        "mars" => Body::Mars,
        "jupiter" => Body::Jupiter,
        "saturn" => Body::Saturn,
        other => panic!("{other}"),
    }
}

/// Angular separation, arcseconds, of two (RA, Dec) pairs in degrees.
fn sep_arcsec(ra1: f64, de1: f64, ra2: f64, de2: f64) -> f64 {
    let v = |ra: f64, de: f64| {
        [(de * DEG).cos() * (ra * DEG).cos(), (de * DEG).cos() * (ra * DEG).sin(), (de * DEG).sin()]
    };
    let (a, b) = (v(ra1, de1), v(ra2, de2));
    norm(cross(a, b)).atan2(dot(a, b)) / DEG * 3600.0
}

#[test]
fn positions_and_horizon_match_horizons() {
    let Some(eng) = engine() else { return };
    let fx = fixture();
    // world -> (worst astrometric gap, worst horizontal gap), arcseconds
    let mut worst: std::collections::BTreeMap<String, (f64, f64)> = Default::default();
    for case in fx["cases"].as_array().unwrap() {
        let world = World::from_name(case["world"].as_str().unwrap()).unwrap();
        let site = Site {
            world,
            lon_deg: case["lon"].as_f64().unwrap(),
            lat_deg: case["lat"].as_f64().unwrap(),
            height_km: case["height_km"].as_f64().unwrap(),
        };
        let target = body(case["target"].as_str().unwrap());
        let w = worst.entry(case["world"].as_str().unwrap().to_string()).or_default();
        for row in case["rows"].as_array().unwrap() {
            let jd_tt = row["jd_tt"].as_f64().unwrap();
            let jd_ut = eng.tt_to_ut(jd_tt);
            let ep = eng.space_epoch(jd_ut);
            let obs = eng.space_observer_state(&ep, &SpaceObserver::Surface(site)).unwrap();
            let naif = target.naif().unwrap();
            let v = apparent(eng.pack(), &ep, &obs, naif, Corrections::ASTROMETRIC, LightTime::Iterate)
                .unwrap()
                .0;
            let (ra, de, _) = to_spherical(v);
            let d = sep_arcsec(ra / DEG, de / DEG, row["ra"].as_f64().unwrap(), row["dec"].as_f64().unwrap());
            w.0 = w.0.max(d);

            let h = eng.space_horizontal(target, jd_ut, &site).unwrap();
            let (az, el) = (row["az"].as_f64().unwrap(), row["el"].as_f64().unwrap());
            // Compare as a direction on the sky: azimuth differences shrink near the zenith.
            let dh = sep_arcsec(h.azimuth, h.elevation, az, el);
            w.1 = w.1.max(dh);
        }
    }
    for (world, (pos, hor)) in &worst {
        eprintln!("{world:>8}: astrometric {pos:.4}\", horizontal {hor:.4}\"");
        if world == "earth" {
            // Control: Horizons rotates the Earth with measured UT1, the engine with its
            // Delta T model; seconds of UT become tens of arcseconds on the horizon.
            assert!(*pos < 1.0 && *hor < 120.0, "earth {pos} {hor}");
        } else if world == "moon" {
            // Horizons orients the Moon with the JPL libration angles; the IAU 2015 series
            // approximates them to ~100 m on the surface (~12" on the horizon).
            assert!(*pos < 0.2 && *hor < 20.0, "moon {pos} {hor}");
        } else {
            // The pack is DE442, Horizons serves DE441: the two place Uranus 0.35" apart
            // (Neptune 0.1"), and an observer standing there inherits the gap.
            assert!(*pos < 0.4, "{world} astrometric {pos}");
            assert!(*hor < 0.4, "{world} horizontal {hor}");
        }
    }
}

/// The Ascendant must be the zodiac point on the eastern horizon, rising; the MC on the
/// local meridian, on the observer's side of the spin axis.
#[test]
fn angles_lie_on_horizon_and_meridian_for_every_world() {
    let Some(eng) = engine() else { return };
    let jd_ut = 2460600.3;
    for world in World::ALL {
        for &(lon, lat) in &[(0.0, 0.0), (77.0, 18.0), (200.0, -45.0), (310.0, 60.0)] {
            for zodiac in [Zodiac::Earth, Zodiac::Native] {
                let site = Site::new(world, lon, lat);
                // Regiomontanus is defined at every latitude.
                let h = eng.space_houses(jd_ut, &site, 'R', zodiac).unwrap();
                let ep = eng.space_epoch(jd_ut);
                let f = eng.world_frame(&ep, world);
                let m = f.to_body;
                let up = {
                    let (a, c) = world.radii_km();
                    // ellipsoid normal at planetodetic latitude = the latitude itself
                    let _ = (a, c);
                    let (sl, cl) = (lat * DEG).sin_cos();
                    let (so, co) = (lon * DEG).sin_cos();
                    mat_t_vec(&m, [cl * co, cl * so, sl])
                };
                let a = f.spin_axis();
                let east = {
                    let e = cross(a, up);
                    let n = norm(e);
                    [e[0] / n, e[1] / n, e[2] / n]
                };
                let dir = |lon: f64| eng.space_zodiac_direction(&ep, Some(world), zodiac, lon).unwrap();
                let asc = dir(h.houses.asc);
                let mc = dir(h.houses.mc);
                let tag = format!("{world:?} {lon} {lat} {zodiac:?}");
                assert!(dot(asc, up).abs() < 1e-9, "Asc off the horizon: {tag}");
                // Inside the polar circle the house convention may take the setting point
                // and swap MC and IC (see the houses module docs).
                let polar = fe_core::houses::is_polar(h.spin_latitude, h.obliquity.min(180.0 - h.obliquity));
                if !polar {
                    assert!(dot(asc, east) > 0.0, "Asc not in the east: {tag}");
                }
                assert!(dot(mc, east).abs() < 1e-9, "MC off the meridian: {tag}");
                // MC on the meridian half that contains the observer's zenith projection.
                let proj = |v: [f64; 3]| {
                    let k = dot(v, a);
                    [v[0] - k * a[0], v[1] - k * a[1], v[2] - k * a[2]]
                };
                if !polar {
                    assert!(dot(proj(mc), proj(up)) > 0.0, "MC is the IC: {tag}");
                }
                let _ = mat_vec(&m, up);
            }
        }
    }
}

/// On Earth, with Earth's zodiac, the extension must agree with the ordinary chart up to the
/// IAU Earth rotation model, which ignores nutation and the exact sidereal time (< 0.02°).
#[test]
fn earth_site_matches_the_ordinary_chart() {
    let Some(eng) = engine() else { return };
    let jd_ut = 2448029.104167;
    let site = Site::new(World::Earth, 2.35, 48.85);
    let s = eng.space_houses(jd_ut, &site, 'P', Zodiac::Earth).unwrap();
    let e = eng.houses_ut(jd_ut, 48.85, 2.35, 'P', None).unwrap();
    for (a, b) in s.houses.cusps.iter().zip(e.cusps.iter()) {
        let d = fe_core::math::diff_deg(*a, *b).abs();
        assert!(d < 0.02, "cusp {a} vs {b}");
    }
    let p = eng
        .space_position(
            Body::Moon,
            jd_ut,
            &SpaceObserver::Center(World::Earth),
            Zodiac::Earth,
            Corrections::APPARENT,
        )
        .unwrap();
    let q = eng.calc_ut(Body::Moon, jd_ut, &Default::default()).unwrap();
    assert!(fe_core::math::diff_deg(p.lon, q.lon).abs() * 3600.0 < 0.01, "{} vs {}", p.lon, q.lon);
    assert!((p.lon_speed - q.lon_speed).abs() < 1e-4);
}

/// Earth's native zodiac is the ecliptic of its own orbit, so it must sit within a few
/// arcseconds of the classic tropical zodiac (osculating EMB orbit against the ecliptic).
#[test]
fn earth_native_zodiac_is_the_tropical_zodiac() {
    let Some(eng) = engine() else { return };
    let jd_ut = 2460000.5;
    let obs = SpaceObserver::Center(World::Earth);
    for b in [Body::Sun, Body::Mars, Body::Jupiter] {
        let n = eng.space_position(b, jd_ut, &obs, Zodiac::Native, Corrections::APPARENT).unwrap();
        let e = eng.space_position(b, jd_ut, &obs, Zodiac::Earth, Corrections::APPARENT).unwrap();
        let d = fe_core::math::diff_deg(n.lon, e.lon).abs() * 3600.0;
        assert!(d < 30.0, "{b:?}: native {} vs tropical {} ({d}\")", n.lon, e.lon);
    }
}

#[test]
fn observer_cannot_target_its_own_world() {
    let Some(eng) = engine() else { return };
    let site = Site::new(World::Mars, 0.0, 0.0);
    assert!(eng
        .space_position(
            Body::Mars,
            2460000.5,
            &SpaceObserver::Surface(site),
            Zodiac::Earth,
            Corrections::APPARENT
        )
        .is_err());
    assert!(eng.space_houses(2460000.5, &site, 'G', Zodiac::Earth).is_err());
    assert!(eng
        .space_position(
            Body::Sun,
            2460000.5,
            &SpaceObserver::Point { center: 0, offset_km: [1e8, 0.0, 0.0], velocity_km_day: [0.0; 3] },
            Zodiac::Native,
            Corrections::APPARENT
        )
        .is_err());
}
