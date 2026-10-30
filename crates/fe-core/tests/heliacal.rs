// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Heliacal events / limiting magnitude vs the Swiss Ephemeris 2.10 black-box
//! oracle (fixtures: tests/fixtures/heliacal/gen_heliacal.py). Skips without a pack.
//! Reported, not just asserted: the visibility models differ (see docs/HELIACAL_REPORT.md).

#[macro_use]
mod common;

use fe_core::context::{Engine, EngineConfig};
use fe_core::ephem::Observer;
use fe_core::heliacal::*;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn engine() -> Option<Engine> {
    let name = "fepk-de441-de442-v1.bin";
    let p = std::env::var("FE_PACK").map(PathBuf::from).ok().or_else(|| {
        [
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs").join(name),
            PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()).join(name),
        ]
        .into_iter()
        .find(|p| p.exists())
    });
    match p.filter(|p| p.exists()) {
        Some(p) => {
            Some(Engine::with_config(fe_core::pack::Pack::open(&p).unwrap(), EngineConfig::SWISS_TIME))
        }
        None => {
            eprintln!("SKIP: no FEPK pack found (set FE_PACK)");
            None
        }
    }
}

fn fixture() -> Value {
    let p = common::oracle_path("heliacal/heliacal.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

#[derive(Default)]
struct Stats(BTreeMap<String, Vec<f64>>);
impl Stats {
    fn add(&mut self, k: &str, v: f64) {
        self.0.entry(k.to_string()).or_default().push(v);
    }
    /// Signed values: report median, mean |x|, p90 |x|, max |x|.
    fn report(&self, unit: &str) -> BTreeMap<String, (usize, f64, f64, f64)> {
        let mut out = BTreeMap::new();
        for (k, v) in &self.0 {
            let mut s = v.clone();
            s.sort_by(f64::total_cmp);
            let n = s.len();
            let med = s[n / 2];
            let mut a: Vec<f64> = v.iter().map(|x| x.abs()).collect();
            a.sort_by(f64::total_cmp);
            let mean = a.iter().sum::<f64>() / n as f64;
            let p90 = a[((n as f64 - 1.0) * 0.9).round() as usize];
            eprintln!("{k:<40} n={n:<4} median {med:>+8.3} mean|.| {mean:>7.3} p90|.| {p90:>7.3} max|.| {:>7.3} {unit}", a[n - 1]);
            out.insert(k.clone(), (n, mean, p90, a[n - 1]));
        }
        out
    }
}

fn atm(v: &Value) -> HelAtmosphere {
    let a: Vec<f64> = v.as_array().unwrap().iter().map(f).collect();
    HelAtmosphere { pressure_hpa: a[0], temperature_c: a[1], humidity_pct: a[2], vr_or_k: a[3] }
}

fn flags(v: i64) -> HelFlags {
    HelFlags { optical: v & 512 != 0, dark: v & 4096 != 0, no_moon: v & 8192 != 0, photopic: v & 16384 != 0 }
}

fn kind(obj: &str) -> &'static str {
    // (the oracle lower-cases the object-name buffer it is given, so compare case-insensitively)
    match obj.to_ascii_lowercase().as_str() {
        "venus" | "mercury" => "inner",
        "mars" | "jupiter" | "saturn" => "outer",
        "moon" => "moon",
        _ => "star",
    }
}

#[allow(dead_code)]
fn kind_exact(obj: &str) -> &'static str {
    match obj {
        "Venus" | "Mercury" => "inner",
        "Mars" | "Jupiter" | "Saturn" => "outer",
        "Moon" => "moon",
        _ => "star",
    }
}

#[test]
fn heliacal_vs_swiss() {
    require_oracle!("heliacal/heliacal.json");
    let Some(e) = engine() else { return };
    let fx = fixture();
    let mut days = Stats::default();
    let mut hours = Stats::default();
    let mut av = Stats::default();
    let (mut both_err, mut class_diff, mut n) = (0, Vec::new(), 0);
    let mut apparition: Vec<String> = Vec::new();
    let t0 = std::time::Instant::now();
    for c in fx["heliacal"].as_array().unwrap() {
        n += 1;
        let g: Vec<f64> = c["geo"].as_array().unwrap().iter().map(f).collect();
        let obs = Observer { lon_deg: g[0], lat_deg: g[1], height_m: g[2] };
        let obj = c["obj"].as_str().unwrap();
        let ev = c["ev"].as_i64().unwrap() as i32;
        let target = HelTarget::parse(obj).unwrap();
        let fl = flags(c["flags"].as_i64().unwrap());
        let ours = HelEvent::from_id(ev).and_then(|evt| {
            e.heliacal_ut(f(&c["jd"]), &obs, &atm(&c["atmo"]), &HelObserver::default(), target, evt, fl)
                .map(|r| (r, evt))
        });
        match (&ours, c.get("t")) {
            (Ok((r, evt)), Some(t)) => {
                let ts = f(&t[0]);
                let label = format!("{} ev{ev}", kind(obj));
                // Event day (local) difference, and start-time difference.
                let dd = (r[0] - ts + 0.5).floor();
                if (r[0] - ts).abs() > 20.0 {
                    // The two models disagree on whether one apparition is visible at all
                    // (marginal: low declination, faint, high latitude): the next event is
                    // a whole apparition away.
                    apparition.push(format!("{} {} {obj} ev{ev}: {:+.0} d", c["site"], c["year"], r[0] - ts));
                } else {
                    days.add(&label, (r[0] - ts).round());
                }
                if (r[0] - ts).abs() > 20.0 && std::env::var("FE_DEBUG").is_ok() {
                    eprintln!(
                        "BIG {} {} {obj} ev{ev}: fe {:.2} swiss {:.2} ({:+.1} d)",
                        c["site"],
                        c["year"],
                        r[0],
                        ts,
                        r[0] - ts
                    );
                }
                if dd == 0.0 {
                    hours.add(&label, (r[0] - ts) * 24.0);
                }
                // Arcus visionis (topocentric, object − Sun altitude) at the start of visibility.
                let p = e
                    .heliacal_pheno_ut(
                        r[0],
                        &obs,
                        &atm(&c["atmo"]),
                        &HelObserver::default(),
                        target,
                        *evt,
                        fl,
                    )
                    .unwrap();
                let ps = &c["pheno_start"];
                av.add(&format!("{label} TAV at first visibility"), p[6] - f(&ps[6]));
                if obj.eq_ignore_ascii_case("sirius") && ev == 1 {
                    av.add("Sirius ev1 TAV ours (abs value)", p[6]);
                    av.add("Sirius ev1 TAV swiss (abs value)", f(&ps[6]));
                    days.add(&format!("Sirius ev1 {}", c["site"].as_str().unwrap()), (r[0] - ts).round());
                }
            }
            (Err(_), None) => both_err += 1,
            _ => class_diff.push(format!(
                "{} {} {obj} ev{ev}: fe {:?} swiss {}",
                c["site"],
                c["year"],
                ours.as_ref().map(|x| x.0[0]),
                c.get("t").map(|t| t[0].to_string()).unwrap_or_else(|| c["err"].to_string())
            )),
        }
    }
    eprintln!(
        "heliacal: {n} cases, {:.0} ms/case (debug build)",
        t0.elapsed().as_secs_f64() * 1e3 / n as f64
    );
    eprintln!("-- event day difference (fe − swiss), days");
    let rd = days.report("d");
    eprintln!("-- start of visibility on the same day, hours");
    hours.report("h");
    eprintln!("-- arcus visionis, deg");
    av.report("deg");
    for m in &class_diff {
        eprintln!("CLASS {m}");
    }
    for m in &apparition {
        eprintln!("APPARITION {m}");
    }
    eprintln!(
        "both no-event: {both_err}; event/no-event differences: {}; different apparition (|Δ| > 20 d): {} / {n}",
        class_diff.len(),
        apparition.len()
    );
    // Model-level agreement bounds (the visibility models differ; not Swiss-exact).
    for (k, (_, mean, _, max)) in &rd {
        assert!(*mean < 4.0 && *max <= 20.0, "{k}: mean {mean} d, max {max} d");
    }
    assert!(class_diff.is_empty());
    assert!(apparition.len() * 25 <= n, "too many apparition mismatches");
}

#[test]
fn vis_limit_mag_vs_swiss() {
    require_oracle!("heliacal/heliacal.json");
    let Some(e) = engine() else { return };
    let fx = fixture();
    let mut st = Stats::default();
    let (mut res_diff, mut n) = (0, 0);
    for c in fx["vis"].as_array().unwrap() {
        if c.get("err").is_some() {
            continue;
        }
        let g: Vec<f64> = c["geo"].as_array().unwrap().iter().map(f).collect();
        let obs = Observer { lon_deg: g[0], lat_deg: g[1], height_m: g[2] };
        let target = HelTarget::parse(c["obj"].as_str().unwrap()).unwrap();
        let fl = flags(c["flags"].as_i64().unwrap());
        let r = e
            .vis_limit_mag(f(&c["jd"]), &obs, &atm(&c["atmo"]), &HelObserver::default(), target, fl)
            .unwrap();
        let d = c["dret"].as_array().unwrap();
        let sres = f(&c["res"]) as i32;
        n += 1;
        if std::env::var("FE_DEBUG").is_ok() {
            eprintln!("ROW {} res {}/{} lim {:.2}/{:.2} alt {:.2}/{:.2} sun {:.2}/{:.2} moon {:.1}/{:.1} mag {:.2}/{:.2} fl {:?}", c["obj"], r.res, sres, r.dret[0], f(&d[0]), r.dret[1], f(&d[1]), r.dret[3], f(&d[3]), r.dret[5], f(&d[5]), r.dret[7], f(&d[7]), c["flags"]);
        }
        if sres >= 0 {
            st.add("object altitude (deg)", r.dret[1] - f(&d[1]));
            st.add(&format!("object magnitude {}", kind(c["obj"].as_str().unwrap())), r.dret[7] - f(&d[7]));
            if !fl.dark {
                st.add("sun altitude (deg)", r.dret[3] - f(&d[3]));
            }
        }
        if (r.res == -2) != (sres == -2) {
            res_diff += 1;
            continue;
        }
        if sres < 0 {
            continue;
        }
        let sun = f(&d[3]);
        let regime = if fl.dark || sun < -18.0 {
            "dark"
        } else if sun < 0.0 {
            "twilight"
        } else {
            "day"
        };
        st.add(&format!("limiting mag, {regime}"), r.dret[0] - f(&d[0]));
        if (r.res == 0) != (sres == 0) {
            st.add("photopic/scotopic flag differs", 1.0);
        }
    }
    st.report("");
    eprintln!("above/below horizon differences: {res_diff} / {n}");
    assert!(res_diff * 50 <= n);
}
