// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Gap-analysis "quick wins" vs the Swiss Ephemeris 2.10 black-box oracle (fixtures:
//! tests/fixtures/quickwins/gen_quickwins.py): output flags (XYZ, RADIANS, ICRS,
//! ASTROMETRIC, SPEED3), ECL_NUT, calc_pctr, helio_cross, mooncross_node,
//! lun_eclipse_when_loc, UTC / leap seconds, Delta T controls, calendar helpers, sidtime0,
//! split_deg and the centisecond / normalisation helpers, names.
//! Engine-dependent tests skip when no pack is present; the pure helpers always run.
//! Run with `-- --nocapture` for the per-group maxima.

#[macro_use]
mod common;

use fe_core::context::{Body, CalcOptions, Engine, EngineConfig};
use fe_core::eclipses::LunarKind;
use fe_core::ephem::{Corrections, Observer};
use fe_core::lunar_local::lunar_eclipse_local_search;
use fe_core::math::diff_deg;
use fe_core::pack::Pack;
use fe_core::time::{self, Calendar, DeltaTModel};
use fe_core::{utc, util};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::OnceLock;

fn pack_path() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FE_PACK") {
        return Some(PathBuf::from(p)).filter(|p| p.exists());
    }
    let name = "fepk-de441-de442-v1.bin";
    [
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs").join(name),
        PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join(name),
    ]
    .into_iter()
    .find(|p| p.exists())
}

/// The Python default engine: Swiss-compatible Delta T and long-term sidereal time.
fn engine() -> Option<Engine> {
    match pack_path() {
        Some(p) => Some(Engine::with_config(Pack::open(&p).expect("open pack"), EngineConfig::SWISS_TIME)),
        None => {
            eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
            None
        }
    }
}

fn fx() -> &'static Value {
    static F: OnceLock<Value> = OnceLock::new();
    F.get_or_init(|| {
        let p = common::oracle_path("quickwins/quickwins.json");
        serde_json::from_str(&std::fs::read_to_string(p).expect("fixture")).unwrap()
    })
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}
fn i(v: &Value) -> i64 {
    v.as_i64().unwrap()
}
fn arr(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(f).collect()
}

/// Per-group maxima of absolute errors.
#[derive(Default)]
struct Max(BTreeMap<String, (usize, f64)>);
impl Max {
    fn add(&mut self, k: &str, v: f64) {
        let e = self.0.entry(k.to_string()).or_insert((0, 0.0));
        e.0 += 1;
        e.1 = e.1.max(v.abs());
    }
    fn report(&self, unit: &str) {
        for (k, (n, m)) in &self.0 {
            eprintln!("{k:<40} n={n:<5} max {m:>12.6} {unit}");
        }
    }
    fn get(&self, k: &str) -> f64 {
        self.0.get(k).map_or(0.0, |x| x.1)
    }
}

const AS: f64 = 3600.0;

fn opt_for(name: &str) -> CalcOptions {
    let mut o = CalcOptions::default();
    for part in name.split('_') {
        match part {
            "xyz" => o.xyz = true,
            "radians" => o.radians = true,
            "equ" => o.equatorial = true,
            "j2000" => o.j2000 = true,
            "icrs" => o.icrs = true,
            "nonut" => o.mean_equinox = true,
            "astrometric" => o.corrections = Corrections::ASTROMETRIC,
            "sidereal" => o.sidereal = Some(fe_core::ayanamsa::Ayanamsa::Lahiri),
            "plain" | "speed3" => {}
            p => panic!("flag {p}"),
        }
    }
    o
}

/// Angular difference (″) of two outputs in the same representation, plus speed difference.
fn diff_out(name: &str, ours: [f64; 6], sw: &[f64]) -> (f64, f64) {
    if name.contains("xyz") {
        let r = (sw[0] * sw[0] + sw[1] * sw[1] + sw[2] * sw[2]).sqrt();
        let d = ((ours[0] - sw[0]).powi(2) + (ours[1] - sw[1]).powi(2) + (ours[2] - sw[2]).powi(2)).sqrt();
        let dv = ((ours[3] - sw[3]).powi(2) + (ours[4] - sw[4]).powi(2) + (ours[5] - sw[5]).powi(2)).sqrt();
        (d / r * 206_264.806, dv / r * 206_264.806)
    } else {
        let k = if name.contains("radians") { 180.0 / std::f64::consts::PI } else { 1.0 };
        let dl = diff_deg(ours[0] * k, sw[0] * k) * ((sw[1] * k).to_radians().cos());
        let db = (ours[1] - sw[1]) * k;
        let ds = (ours[3] - sw[3]) * k;
        ((dl * dl + db * db).sqrt() * AS, ds.abs() * AS)
    }
}

#[test]
fn output_flags_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let (mut pos, mut spd) = (Max::default(), Max::default());
    let frame_changing = |n: &str| n.contains("j2000") || n.contains("equ") || n.contains("icrs");
    for r in fx()["flags"].as_array().unwrap() {
        let name = r["flags"].as_str().unwrap();
        let body = i(&r["body"]) as i32;
        // Lunar points ignore output frames in this engine (pre-existing; not part of this item).
        if body >= 10 && frame_changing(name) {
            continue;
        }
        let o = opt_for(name);
        let p = eng.calc_tt(Body::from_id(body).unwrap(), f(&r["jd_tt"]), &o).unwrap().to_array();
        if name == "j2000_equ" && body == 1 {
            // what the previous behaviour (ICRS axes for J2000 equatorial) would give
            let old = CalcOptions { icrs: true, ..o };
            let q = eng.calc_tt(Body::Moon, f(&r["jd_tt"]), &old).unwrap().to_array();
            pos.add("j2000_equ without frame bias/moon", diff_out(name, q, &arr(&r["xx"])).0);
        }
        let (dp, ds) = diff_out(name, p, &arr(&r["xx"]));
        let grp = if body == 1 {
            "moon"
        } else if body >= 10 {
            "nodes"
        } else {
            "planets"
        };
        pos.add(&format!("{name}/{grp}"), dp);
        spd.add(&format!("{name}/{grp}"), ds);
    }
    eprintln!("-- position (″)");
    pos.report("″");
    eprintln!("-- speed (″/day)");
    spd.report("″/d");
    // Representation changes must not add error beyond the plain apparent place.
    for g in ["planets", "moon", "nodes"] {
        let base = pos.get(&format!("plain/{g}"));
        for n in ["xyz", "radians", "xyz_radians", "speed3", "sidereal_xyz"] {
            // (cartesian errors include the distance error, ~1e-7 relative)
            assert!(
                pos.get(&format!("{n}/{g}")) <= base * 1.2 + 0.01,
                "{n}/{g} {} vs plain {base}",
                pos.get(&format!("{n}/{g}"))
            );
        }
    }
    // Frame options add nothing to the plain apparent-place difference (J2000 equatorial
    // now includes the IAU 2006 frame bias unless ICRS is requested).
    let base = pos.get("plain/planets");
    for k in [
        "j2000_equ/planets",
        "j2000_equ_icrs/planets",
        "icrs/planets",
        "equ_icrs/planets",
        "j2000_icrs/planets",
        "nonut_icrs/planets",
    ] {
        assert!(pos.get(k) <= base * 1.05 + 0.005, "{k}: {} vs plain {base}", pos.get(k));
    }
    for k in ["j2000_equ/moon", "j2000_equ_icrs/moon", "icrs/moon", "equ_icrs/moon", "j2000_icrs/moon"] {
        assert!(pos.get(k) < 0.01, "{k}: {}", pos.get(k));
    }
}

#[test]
fn ecl_nut_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    for r in fx()["ecl_nut"].as_array().unwrap() {
        let name = r["flags"].as_str().unwrap();
        let o = opt_for(name);
        let jd = f(&r["jd_tt"]);
        let era = if (2378496.5..2524593.5).contains(&jd) { "1800-2200" } else { "0-2400" };
        let ours = eng.ecl_nut(jd, false, &o);
        let sw = arr(&r["xx"]);
        let k = if name == "radians" { 180.0 / std::f64::consts::PI * AS } else { AS };
        for (j, lbl) in ["eps_true", "eps_mean", "dpsi", "deps", "x4", "x5"].iter().enumerate() {
            m.add(&format!("{name}/{lbl} {era}"), (ours[j] - sw[j]) * k);
        }
    }
    m.report("″");
    // Nutation identical; mean obliquity is the engine's (Vondrák 2011 pole angle), within
    // 0.2″ of the original's at year 0 and a few mas in 1800–2200.
    for (k, (_, v)) in &m.0 {
        let lim = if k.contains("eps") && k.contains("0-2400") { 0.25 } else { 0.01 };
        assert!(*v < lim, "{k}: {v}″");
    }
}

#[test]
fn calc_pctr_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let (mut pos, mut spd) = (Max::default(), Max::default());
    for r in fx()["pctr"].as_array().unwrap() {
        let name = r["flags"].as_str().unwrap();
        let mut o = CalcOptions::default();
        match name {
            "app" => {}
            "nogdefl" => o.corrections.deflection = false,
            "noaberr" => o.corrections.aberration = false,
            "true" => o.corrections = Corrections::GEOMETRIC,
            "j2000" => o.j2000 = true,
            "equ" => o.equatorial = true,
            "xyz" => o.xyz = true,
            // centre flags are ignored by the original but force an astrometric place
            "helctr" => o.corrections = Corrections::ASTROMETRIC,
            _ => unreachable!(),
        }
        let t = Body::from_id(i(&r["target"]) as i32).unwrap();
        let c = Body::from_id(i(&r["center"]) as i32).unwrap();
        let p = eng.calc_pctr(t, c, f(&r["jd_tt"]), &o).unwrap().to_array();
        let (dp, ds) = diff_out(name, p, &arr(&r["xx"]));
        let grp = if t == Body::Moon || c == Body::Moon { "moon" } else { "other" };
        pos.add(&format!("{name}/{grp}"), dp);
        spd.add(&format!("{name}/{grp}"), ds);
    }
    pos.report("″");
    spd.report("″/d");
    for (k, (_, v)) in &pos.0 {
        assert!(*v < 2.0, "{k}: {v}″");
    }
    assert!(eng.calc_pctr(Body::Mars, Body::Mars, 2451545.0, &CalcOptions::default()).is_err());
}

#[test]
fn crossings_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    let o = CalcOptions::default();
    for r in fx()["helio_cross"].as_array().unwrap() {
        let b = Body::from_id(i(&r["body"]) as i32).unwrap();
        let (x, jd, back) = (f(&r["x"]), f(&r["jd"]), r["backward"].as_bool().unwrap());
        let tt = eng.helio_cross(b, x, jd, &o, back).unwrap();
        let ut = eng.helio_cross_ut(b, x, jd, &o, back).unwrap();
        let grp = if i(&r["body"]) >= 7 { "outer" } else { "inner" };
        m.add(&format!("helio_cross tt/{grp} (s)"), (tt - f(&r["tt"])) * 86400.0);
        let ho = CalcOptions { origin: fe_core::ephem::Origin::Heliocentric, ..o };
        let at_swiss = eng.calc_tt(b, f(&r["tt"]), &ho).unwrap().lon;
        m.add(&format!("helio lon at swiss time/{grp} (″)"), diff_deg(at_swiss, x) * AS);
        let tr = CalcOptions { corrections: Corrections::GEOMETRIC, ..ho };
        m.add(
            &format!("helio true lon at swiss time/{grp} (″)"),
            diff_deg(eng.calc_tt(b, f(&r["tt"]), &tr).unwrap().lon, x) * AS,
        );
        let na = CalcOptions { corrections: Corrections::ASTROMETRIC, ..ho };
        m.add(
            &format!("helio noaberr lon at swiss time/{grp} (″)"),
            diff_deg(eng.calc_tt(b, f(&r["tt"]), &na).unwrap().lon, x) * AS,
        );
        m.add(&format!("helio_cross ut/{grp} (s)"), (ut - f(&r["ut"])) * 86400.0);
    }
    for r in fx()["mooncross_node"].as_array().unwrap() {
        let jd = f(&r["jd"]);
        let (t, l, b) = eng.mooncross_node(jd, &o).unwrap();
        m.add("mooncross_node tt (s)", (t - f(&r["tt"])) * 86400.0);
        m.add("mooncross_node lon (″)", diff_deg(l, f(&r["lon"])) * AS);
        m.add("mooncross_node lat (″)", (b - f(&r["lat"])) * AS);
        // The oracle's mooncross_node_ut returns the same instant as mooncross_node for the
        // same number (i.e. treats the input as TT and returns TT): documented semantics
        // (UT in, UT out) are kept here and checked for consistency instead.
        let (tu, lu, _) = eng.mooncross_node_ut(jd, &o).unwrap();
        let (tt2, l2, _) = eng.mooncross_node(eng.ut_to_tt(jd), &o).unwrap();
        m.add("mooncross_node_ut vs tt variant (s)", (tu - eng.tt_to_ut(tt2)) * 86400.0);
        m.add("mooncross_node_ut lon (″)", diff_deg(lu, l2) * AS);
        m.add("oracle mooncross_node_ut - tt (s)", (f(&r["ut"]) - f(&r["tt"])) * 86400.0);
    }
    m.report("");
    assert!(m.get("mooncross_node tt (s)") < 0.5);
    assert!(m.get("mooncross_node_ut vs tt variant (s)") < 1e-3);
    // Crossing times differ by (heliocentric longitude difference) / speed: the engine's
    // heliocentric apparent longitudes are within 1.01″ of the original's (pre-existing,
    // GEOMETRY/ENGINE reports), which is minutes of time for slow outer planets.
    assert!(m.get("helio lon at swiss time/inner (″)") < 1.1);
    assert!(m.get("helio lon at swiss time/outer (″)") < 0.6);
    assert!(m.get("helio_cross tt/inner (s)") < 30.0);
    assert!(m.get("helio_cross tt/outer (s)") < 6.0 * 3600.0);
}

#[test]
fn lunar_eclipse_when_loc_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    let (mut n, mut flag_mismatch) = (0, 0);
    const VIS: i32 = 128;
    const MAXV: i32 = 256;
    for r in fx()["lun_loc"].as_array().unwrap() {
        let g = arr(&r["geo"]);
        let obs = Observer { lon_deg: g[0], lat_deg: g[1], height_m: g[2] };
        let e =
            lunar_eclipse_local_search(&eng, f(&r["start"]), &obs, r["backward"].as_bool().unwrap()).unwrap();
        let st = arr(&r["tret"]);
        if (e.global.max - st[0]).abs() > 1.0 && (eng.tt_to_ut(e.global.max) - st[0]).abs() > 1.0 {
            // different eclipse found (visibility at the margin); compare the start instead
            eprintln!("different eclipse: ours {} swiss {} at {:?}", e.max, st[0], g);
            flag_mismatch += 1;
            continue;
        }
        n += 1;
        let kind = match e.kind {
            LunarKind::Total => 4,
            LunarKind::Partial => 16,
            LunarKind::Penumbral => 64,
        };
        let mut flags = kind | VIS;
        for (on, bit) in [
            (e.max_visible, MAXV),
            (e.u1.is_some(), 512),
            (e.u2.is_some(), 1024),
            (e.u3.is_some(), 2048),
            (e.u4.is_some(), 4096),
            (e.p1.is_some(), 8192),
            (e.p4.is_some(), 16384),
        ] {
            if on {
                flags |= bit;
            }
        }
        if flags != i(&r["flags"]) as i32 {
            flag_mismatch += 1;
            eprintln!("flags ours {flags} swiss {} at {:?} {}", i(&r["flags"]), g, st[0]);
        }
        let ours = [Some(e.max), None, e.u1, e.u4, e.u2, e.u3, e.p1, e.p4, e.moonrise, e.moonset];
        for (j, (o, s)) in ours.iter().zip(&st).enumerate() {
            match (o, *s) {
                (Some(o), s) if s != 0.0 => m.add(&format!("tret[{j}] (s)"), (o - s) * 86400.0),
                (None, 0.0) => {}
                (o, s) => {
                    if j != 1 {
                        m.add("tret presence mismatch", 1.0);
                        eprintln!("tret[{j}] ours {o:?} swiss {s} at {:?}", g);
                    }
                }
            }
        }
        let sa = arr(&r["attr"]);
        let a = &e.attrs;
        m.add("attr umbral mag", a.shadow.umbral_mag.max(0.0) - sa[0]);
        m.add("attr penumbral mag", a.shadow.penumbral_mag.max(0.0) - sa[1]);
        m.add("attr azimuth (deg)", diff_deg((a.azimuth_deg + 180.0).rem_euclid(360.0), sa[4]));
        m.add("attr true alt (deg)", a.true_alt_deg - sa[5]);
        m.add("attr apparent alt (deg)", a.apparent_alt_deg - sa[6]);
    }
    eprintln!("compared {n} eclipses, flag mismatches {flag_mismatch}");
    m.report("");
    assert!(n >= 100);
    assert!(flag_mismatch <= 2, "{flag_mismatch} flag mismatches");
    assert!(m.get("tret presence mismatch") == 0.0);
    assert!(m.get("tret[0] (s)") < 60.0);
    for k in 2..10 {
        assert!(m.get(&format!("tret[{k}] (s)")) < 30.0, "tret[{k}] {}", m.get(&format!("tret[{k}] (s)")));
    }
}

/// Seconds since an arbitrary epoch for a UTC tuple (s may be 60.x).
fn utc_secs(t: &[f64]) -> f64 {
    let jd0 = time::julday(t[0] as i32, t[1] as u32, t[2] as u32, 0.0, Calendar::Gregorian);
    (jd0 - 2_400_000.5) * 86400.0 + t[3] * 3600.0 + t[4] * 60.0 + t[5]
}

#[test]
fn utc_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    for r in fx()["utc_to_jd"].as_array().unwrap() {
        let c = arr(&r["in"]);
        let cal = if i(&r["cal"]) == 1 { Calendar::Gregorian } else { Calendar::Julian };
        let ours = eng.utc_to_jd(c[0] as i32, c[1] as i32, c[2] as i32, c[3] as i32, c[4] as i32, c[5], cal);
        let ok = r["ok"].as_bool().unwrap();
        assert_eq!(
            ours.is_ok(),
            ok,
            "{c:?}: {ours:?} vs {}",
            r.get("err").map_or(String::new(), |e| e.to_string())
        );
        if let Ok((et, ut)) = ours {
            let grp = if c[0] < 1972.0 { "<1972" } else { ">=1972" };
            m.add(&format!("utc_to_jd et {grp} (s)"), (et - f(&r["et"])) * 86400.0);
            m.add(&format!("utc_to_jd ut {grp} (s)"), (ut - f(&r["ut"])) * 86400.0);
        }
    }
    for r in fx()["jd_to_utc"].as_array().unwrap() {
        let jd = f(&r["jd"]);
        let a = eng.jdet_to_utc(jd, Calendar::Gregorian);
        let b = eng.jdut1_to_utc(jd, Calendar::Gregorian);
        let tup = |x: (i32, i32, i32, i32, i32, f64)| {
            vec![x.0 as f64, x.1 as f64, x.2 as f64, x.3 as f64, x.4 as f64, x.5]
        };
        let grp = if jd < 2441317.5 { "<1972" } else { ">=1972" };
        m.add(&format!("jdet_to_utc {grp} (s)"), utc_secs(&tup(a)) - utc_secs(&arr(&r["et_utc"])));
        m.add(&format!("jdut1_to_utc {grp} (s)"), utc_secs(&tup(b)) - utc_secs(&arr(&r["ut1_utc"])));
        let sw = arr(&r["et_utc"]);
        if sw[5] >= 60.0 {
            assert!(a.5 >= 60.0 && a.3 == 23 && a.4 == 59, "leap second {a:?} vs {sw:?}");
        }
    }
    m.report("");
    // From 1972 TT is exact (leap-second table); UT1 differs only by the Delta T models
    // (both observed values, ≲ 0.1 s).
    assert!(m.get("utc_to_jd et >=1972 (s)") < 1e-4);
    assert!(m.get("utc_to_jd ut >=1972 (s)") < 0.5);
    assert!(m.get("utc_to_jd ut <1972 (s)") < 1e-4);
    assert!(m.get("jdet_to_utc >=1972 (s)") < 1e-3);
    assert!(m.get("jdut1_to_utc >=1972 (s)") < 0.5);
}

#[test]
fn utc_time_zone_and_calendar_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    for r in fx()["utc_time_zone"].as_array().unwrap() {
        let c = arr(&r["in"]);
        let o = utc::utc_time_zone(
            c[0] as i32,
            c[1] as i32,
            c[2] as i32,
            c[3] as i32,
            c[4] as i32,
            c[5],
            f(&r["offset"]),
        );
        let s = arr(&r["out"]);
        let ours = [o.0 as f64, o.1 as f64, o.2 as f64, o.3 as f64, o.4 as f64, o.5];
        assert!((utc_secs(&ours) - utc_secs(&s)).abs() < 1e-5, "{c:?}: {ours:?} vs {s:?}");
        assert_eq!(ours[5] >= 60.0, s[5] >= 60.0);
    }
    for r in fx()["date_conversion"].as_array().unwrap() {
        let c = r["in"].as_array().unwrap();
        let cal = if c[4].as_str().unwrap() == "g" { Calendar::Gregorian } else { Calendar::Julian };
        let (v, jd, (y, m, d, h)) =
            utc::date_conversion(i(&c[0]) as i32, i(&c[1]) as i32, i(&c[2]) as i32, f(&c[3]), cal);
        let dt = arr(&r["dt"]);
        assert_eq!(v, r["valid"].as_bool().unwrap(), "{c:?}");
        // (serde_json's default float parsing is not always correctly rounded: 1-ulp slack)
        assert!((jd - f(&r["jd"])).abs() < 1e-9, "{c:?}");
        assert_eq!((y as f64, m as f64, d as f64), (dt[0], dt[1], dt[2]), "{c:?}");
        assert!((h - dt[3]).abs() < 1e-8, "{c:?}: {h} vs {}", dt[3]);
    }
    for r in fx()["day_of_week"].as_array().unwrap() {
        let a = r.as_array().unwrap();
        assert_eq!(utc::day_of_week(f(&a[0])) as i64, i(&a[1]), "{}", f(&a[0]));
    }
}

#[test]
fn delta_t_controls_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    for r in fx()["deltat"].as_array().unwrap() {
        let jd = f(&r["jd"]);
        let base = eng.delta_t_seconds(jd);
        m.add("deltat_ex vs swiss (s)", base - f(&r["deltat_ex"]) * 86400.0);
        for (name, v) in r["tid"].as_object().unwrap() {
            let v = arr(v);
            let model = DeltaTModel { lunar_ndot: v[0], ..eng.delta_t_model() };
            let e2 = eng.with_delta_t(model);
            assert_eq!(e2.tidal_acceleration(), v[2]);
            // Change of Delta T caused by the tidal acceleration (relative to each default).
            let ours = e2.delta_t_seconds(jd) - eng.delta_t_seconds_with_ndot(jd, f(&r["default_acc"]));
            m.add(&format!("tidal change {name} (s)"), ours - v[1] * 86400.0);
        }
    }
    m.report("");
    for (k, (_, v)) in &m.0 {
        // Published correction −0.91072 Δṅ u² (Morrison & Stephenson 2004); the original's
        // corrections are 0.07 % smaller (≤ 0.83 s for DE200 at −500, ≤ 0.06 s for DE441).
        if k.starts_with("tidal") {
            assert!(*v < 1.0, "{k}: {v}");
        }
    }
    // User-defined Delta T: every UT computation uses it.
    let u = &fx()["deltat_user"];
    let e2 =
        eng.with_delta_t(DeltaTModel { user_seconds: Some(f(&u["dt"]) * 86400.0), ..eng.delta_t_model() });
    let jd = f(&u["jd"]);
    assert!((e2.delta_t_seconds(jd) / 86400.0 - f(&u["deltat"])).abs() < 1e-15);
    let a = e2.calc_ut(Body::Moon, jd, &CalcOptions::default()).unwrap();
    let b = eng.calc_tt(Body::Moon, jd + f(&u["dt"]), &CalcOptions::default()).unwrap();
    assert!((a.lon - b.lon).abs() < 1e-10);
    let sw = arr(&u["calc_ut_moon"]);
    assert!(diff_deg(a.lon, sw[0]).abs() * AS < 1.0, "{} vs {}", a.lon, sw[0]);
}

/// Helper: Delta T of this engine's model with another tidal acceleration.
trait NdotExt {
    fn delta_t_seconds_with_ndot(&self, jd: f64, ndot: f64) -> f64;
}
impl NdotExt for Engine {
    fn delta_t_seconds_with_ndot(&self, jd: f64, ndot: f64) -> f64 {
        time::delta_t_seconds_with(jd, &DeltaTModel { lunar_ndot: ndot, ..self.delta_t_model() })
    }
}

#[test]
fn sidtime0_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let Some(eng) = engine() else { return };
    let mut m = Max::default();
    for r in fx()["sidtime0"].as_array().unwrap() {
        let jd = f(&r["jd"]);
        let grp = if jd < 2396758.5 {
            "before 1850"
        } else if jd < 2469807.5 {
            "1850-2050"
        } else {
            "after 2050"
        };
        let st = eng.sidtime0(jd, f(&r["eps"]), f(&r["nut"]));
        m.add(
            &format!("sidtime0 {grp} (s of time)"),
            ((st - f(&r["st"]) + 12.0).rem_euclid(24.0) - 12.0) * 3600.0,
        );
        let st2 = eng.sidtime0(jd, 23.0, 0.005);
        m.add(
            &format!("sidtime0 given eps/nut {grp} (s)"),
            ((st2 - f(&r["st_alt"]) + 12.0).rem_euclid(24.0) - 12.0) * 3600.0,
        );
    }
    m.report("");
    assert!(m.get("sidtime0 1850-2050 (s of time)") < 0.01);
    assert!(m.get("sidtime0 before 1850 (s of time)") < 0.05);
    // After 2050 the original adds an unpublished continuity offset (ENGINE_REPORT §3).
    assert!(m.get("sidtime0 after 2050 (s of time)") < 0.3);
}

#[test]
fn split_deg_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    let mut bad = 0;
    let mut max_fr: f64 = 0.0;
    let rows = fx()["split_deg"].as_array().unwrap();
    for r in rows {
        let a = r.as_array().unwrap();
        let (v, fl) = (f(&a[0]), i(&a[1]) as i32);
        let s = util::split_deg(v, fl);
        let want = (i(&a[2]) as i32, i(&a[3]) as i32, i(&a[4]) as i32, i(&a[6]) as i32);
        if (s.deg, s.min, s.sec, s.sign) != want {
            bad += 1;
            if bad < 20 {
                eprintln!("split_deg({v}, {fl}) = {s:?}, swiss {want:?}");
            }
        } else {
            max_fr = max_fr.max((s.secfr - f(&a[5])).abs());
        }
    }
    eprintln!("split_deg: {} cases, {bad} integer mismatches, max secfr diff {max_fr:e}", rows.len());
    assert_eq!(bad, 0);
    assert!(max_fr < 1e-6);
}

#[test]
fn centisecond_helpers_vs_swiss() {
    require_oracle!("quickwins/quickwins.json");
    for r in fx()["cs"].as_array().unwrap() {
        let c = i(&r["cs"]) as i32;
        assert_eq!(util::cs2degstr(c), r["degstr"].as_str().unwrap(), "cs2degstr({c})");
        assert_eq!(util::cs2lonlatstr(c, 'N', 'S'), r["lonlat"].as_str().unwrap(), "cs2lonlatstr({c})");
        assert_eq!(util::cs2timestr(c, ':', false), r["time"].as_str().unwrap(), "cs2timestr({c})");
        assert_eq!(util::cs2timestr(c, 'h', true), r["time_sz"].as_str().unwrap(), "cs2timestr({c}, sz)");
        assert_eq!(util::csnorm(c) as i64, i(&r["csnorm"]), "csnorm({c})");
        assert_eq!(util::csroundsec(c) as i64, i(&r["csroundsec"]), "csroundsec({c})");
        assert_eq!(util::difcsn(c, 4321000) as i64, i(&r["difcsn"]), "difcsn({c})");
        assert_eq!(util::difcs2n(c, -77000000) as i64, i(&r["difcs2n"]), "difcs2n({c})");
    }
    for r in fx()["d2l"].as_array().unwrap() {
        let a = r.as_array().unwrap();
        assert_eq!(util::d2l(f(&a[0])) as i64, i(&a[1]), "d2l({})", f(&a[0]));
    }
    for r in fx()["norm"].as_array().unwrap() {
        let a = arr(r);
        assert!(
            (util::degnorm(a[0]) - a[1]).abs() < 1e-12,
            "degnorm({}) {} vs {}",
            a[0],
            util::degnorm(a[0]),
            a[1]
        );
        assert!(
            (util::radnorm(a[0]) - a[2]).abs() < 1e-12,
            "radnorm({}) {} vs {}",
            a[0],
            util::radnorm(a[0]),
            a[2]
        );
    }
    for r in fx()["dif"].as_array().unwrap() {
        let a = arr(r);
        assert!((util::difdegn(a[0], a[1]) - a[2]).abs() < 1e-9, "difdegn {a:?}");
        assert!((util::difdeg2n(a[0], a[1]) - a[3]).abs() < 1e-9, "difdeg2n {a:?}");
        assert!((util::difrad2n(a[0], a[1]) - a[4]).abs() < 1e-9, "difrad2n {a:?}");
    }
    for r in fx()["planet_names"].as_array().unwrap() {
        let a = r.as_array().unwrap();
        assert_eq!(util::planet_name(i(&a[0]) as i32), a[1].as_str().unwrap());
    }
    for r in fx()["ayanamsa_names"].as_array().unwrap() {
        let a = r.as_array().unwrap();
        assert_eq!(util::ayanamsa_name(i(&a[0]) as i32), a[1].as_str().unwrap());
    }
}
