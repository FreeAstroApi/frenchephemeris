// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Fixed stars vs ERFA (space motion + apparent chain, no pack needed) and vs
//! the Swiss Ephemeris oracle (apparent ecliptic places 0..2400 CE, needs a
//! pack; skipped when none is found). Fixtures: tests/fixtures/stars/gen_*.py.

#[macro_use]
mod common;

use fe_core::context::{CalcOptions, Engine};
use fe_core::ephem::{Corrections, ObsState};
use fe_core::math::*;
use fe_core::stars::{self, Astrometry};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

const UAS: f64 = ARCSEC * 1e-6;

fn fixture(name: &str) -> String {
    if name.starts_with("swiss") {
        return common::read_oracle(&format!("stars/{name}"));
    }
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/stars").join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn pack_path() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FE_PACK") {
        return Some(PathBuf::from(p)).filter(|p| p.exists());
    }
    let name = "fepk-de441-de442-v1.bin";
    let mut dirs = Vec::new();
    if let Ok(d) = std::env::var("FE_DATA_DIR") {
        dirs.push(PathBuf::from(d).join("packs"));
    }
    dirs.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/packs"));
    dirs.push(PathBuf::from(std::env::var("FE_PACK_DIR").unwrap_or_default()));
    dirs.into_iter().map(|d| d.join(name)).find(|p| p.exists())
}

fn angle(a: Vec3, b: Vec3) -> f64 {
    norm(cross(a, b)).atan2(dot(a, b))
}

#[derive(Deserialize)]
struct ErfaRow {
    /// Index into `stars::catalog()`.
    idx: usize,
    tdb: f64,
    pm: [f64; 6],
    obs: [f64; 9],
    astrometric: [f64; 3],
    apparent: [f64; 3],
    dist_au: f64,
}

#[derive(Deserialize)]
struct ErfaFixture {
    catalog_epoch: f64,
    rows: Vec<ErfaRow>,
}

#[test]
fn space_motion_and_apparent_chain_match_erfa() {
    let fx: ErfaFixture = serde_json::from_str(&fixture("erfa.json")).unwrap();
    assert_eq!(fx.catalog_epoch, stars::CATALOG_EPOCH_JD);
    let mut max_pos_rv0 = 0.0f64;
    let (mut max_pos, mut max_pm, mut max_px, mut max_rv, mut max_ast, mut max_app, mut max_dist) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for r in &fx.rows {
        let rec = &stars::catalog()[r.idx];
        // (a) pmsafe.
        let a: Astrometry = stars::pmsafe(&rec.astrometry(), stars::CATALOG_EPOCH_JD, r.tdb).unwrap();
        let [ra, dec, pmr, pmd, px, rv] = r.pm;
        let dpos =
            ((a.ra - ra).rem_euclid(TAU).min((ra - a.ra).rem_euclid(TAU)) * dec.cos()).hypot(a.dec - dec);
        if rec.rv_km_s == 0.0 {
            // SOFA's starpv divides by the (rounding-level) radial speed when RV = 0, so
            // libm differences give tens of µas after 2000 years; both are equally valid.
            max_pos_rv0 = max_pos_rv0.max(dpos);
            continue;
        }
        max_pos = max_pos.max(dpos);
        max_pm = max_pm.max(((a.pmr - pmr) * dec.cos()).hypot(a.pmd - pmd));
        max_px = max_px.max((a.px - px).abs());
        max_rv = max_rv.max((a.rv - rv).abs());

        // (b) parallax + ldsun + ab, from the same observer state.
        let km = |i: usize| [r.obs[i] * AU_KM, r.obs[i + 1] * AU_KM, r.obs[i + 2] * AU_KM];
        let obs = ObsState { pos: km(0), vel: km(3), sun: km(6), earthbound: true, barycentric: false };
        let p = stars::barycentric_au(rec, r.tdb);
        let (u0, d0) = stars::apparent_gcrs(p, &obs, Corrections::GEOMETRIC);
        let (u1, _) = stars::apparent_gcrs(p, &obs, Corrections::APPARENT);
        max_ast = max_ast.max(angle(u0, r.astrometric));
        max_app = max_app.max(angle(u1, r.apparent));
        max_dist = max_dist.max(((d0 - r.dist_au) / r.dist_au).abs());
    }
    println!(
        "ERFA ({} rows): pos {:.3e} µas, pm {:.3e} µas/yr, px {:.3e} µas, rv {:.3e} km/s, \
         astrometric dir {:.3e} µas, apparent dir {:.3e} µas, dist rel {:.1e}",
        fx.rows.len(),
        max_pos / UAS,
        max_pm / UAS,
        max_px * 1e6,
        max_rv,
        max_ast / UAS,
        max_app / UAS,
        max_dist
    );
    println!(
        "  stars with RV = 0 (ill-conditioned in SOFA starpv): max position difference {:.1} µas",
        max_pos_rv0 / UAS
    );
    assert!(max_pos < 1.0 * UAS, "pmsafe position {} µas", max_pos / UAS);
    assert!(max_pos_rv0 < 100.0 * UAS);
    assert!(max_pm < 1e-3 * UAS);
    assert!(max_px < 1e-9 && max_rv < 1e-9);
    assert!(max_ast < 1.0 * UAS && max_app < 1.0 * UAS, "{} {}", max_ast / UAS, max_app / UAS);
    assert!(max_dist < 1e-12);
}

#[derive(Deserialize)]
struct SwissRow {
    name: String,
    year: i32,
    jd_ut: f64,
    ecl: [f64; 6],
    equ: [f64; 2],
    #[serde(default)]
    fd: Option<[f64; 4]>,
}

#[derive(Deserialize)]
struct SwissFixture {
    swe_version: String,
    rows: Vec<SwissRow>,
}

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if v.is_empty() {
        return f64::NAN;
    }
    v[((v.len() - 1) as f64 * q).round() as usize]
}

/// Names whose classic star differs from the IAU one, and the alpha Cen names
/// (compatibility-table test).
const COMPAT: [&str; 9] = [
    "Algedi",
    "Deneb Kaitos",
    "Kurhah",
    "Menkar",
    "Sadalbari",
    "Ukdah",
    "Rigil Kentaurus",
    "Toliman",
    "Bungula",
];

/// FreeAstroAPI's `DEFAULT_ASTROLOGY_FIXED_STARS`.
const API_DEFAULT: [&str; 25] = [
    "Algol",
    "Alcyone",
    "Aldebaran",
    "Rigel",
    "Betelgeuse",
    "Sirius",
    "Canopus",
    "Procyon",
    "Castor",
    "Pollux",
    "Regulus",
    "Denebola",
    "Alkaid",
    "Spica",
    "Arcturus",
    "Agena",
    "Toliman",
    "Unukalhai",
    "Antares",
    "Vega",
    "Altair",
    "Deneb",
    "Fomalhaut",
    "Scheat",
    "Markab",
];

struct Stats {
    by_year: BTreeMap<i32, Vec<f64>>,
    /// Same, API default stars only, eras 0-1500 / 1600-2400.
    default_early: Vec<f64>,
    default_late: Vec<f64>,
    /// Max separation per COMPAT name over all eras.
    compat: BTreeMap<String, f64>,
    all: Vec<f64>,
    worst: Vec<(f64, String, i32)>,
    max_speed: f64,
    /// |our speed − Swiss's own positions differentiated over ±0.05 d|, and the same for Swiss's reported speed ("/day).
    speed_fd: Vec<f64>,
    swiss_speed_fd: Vec<f64>,
    dist_rel: Vec<f64>,
    max_equ: f64,
}

fn compare(eng: &Engine, fx: &SwissFixture) -> Stats {
    let opt = CalcOptions::default();
    let equ = CalcOptions { equatorial: true, speed: false, ..CalcOptions::default() };
    let mut st = Stats {
        by_year: BTreeMap::new(),
        default_early: vec![],
        default_late: vec![],
        compat: BTreeMap::new(),
        all: vec![],
        worst: vec![],
        max_speed: 0.0,
        speed_fd: vec![],
        swiss_speed_fd: vec![],
        dist_rel: vec![],
        max_equ: 0.0,
    };
    for r in &fx.rows {
        let (_, p) = eng.fixstar_ut(&r.name, r.jd_ut, &opt).unwrap();
        let dl = diff_deg(p.lon, r.ecl[0]) * 3600.0;
        let db = (p.lat - r.ecl[1]) * 3600.0;
        let sep = (dl * (r.ecl[1] * DEG).cos()).hypot(db);
        st.by_year.entry(r.year).or_default().push(sep);
        st.all.push(sep);
        if API_DEFAULT.contains(&r.name.as_str()) {
            if r.year <= 1500 {
                st.default_early.push(sep)
            } else {
                st.default_late.push(sep)
            }
        }
        st.worst.push((sep, r.name.clone(), r.year));
        if COMPAT.contains(&r.name.as_str()) {
            let m = st.compat.entry(r.name.clone()).or_insert(0.0);
            *m = m.max(sep);
        }
        st.max_speed = st.max_speed.max((p.lon_speed - r.ecl[3]).abs() * 3600.0);
        st.dist_rel.push(((p.dist - r.ecl[2]) / r.ecl[2]).abs());
        if let Some(fd) = r.fd {
            let cosb = (r.ecl[1] * DEG).cos();
            let fd_lon = diff_deg(fd[2], fd[0]) / 0.1;
            let fd_lat = (fd[3] - fd[1]) / 0.1;
            let e = |dl: f64, db: f64| (dl * cosb).hypot(db) * 3600.0;
            st.speed_fd.push(e(p.lon_speed - fd_lon, p.lat_speed - fd_lat));
            st.swiss_speed_fd.push(e(r.ecl[3] - fd_lon, r.ecl[4] - fd_lat));
        }
        let (_, q) = eng.fixstar_ut(&r.name, r.jd_ut, &equ).unwrap();
        let sq =
            (diff_deg(q.lon, r.equ[0]) * 3600.0 * (r.equ[1] * DEG).cos()).hypot((q.lat - r.equ[1]) * 3600.0);
        st.max_equ = st.max_equ.max(sq);
    }
    st
}

fn report(label: &str, st: &mut Stats) -> (f64, f64, f64) {
    println!("{label}: {} rows (ecliptic-of-date separation, arcsec)", st.all.len());
    println!("  year    median       p95       max");
    for (y, v) in st.by_year.iter_mut() {
        println!("  {y:>4}  {:>8.4}  {:>8.4}  {:>8.3}", pct(v, 0.5), pct(v, 0.95), pct(v, 1.0));
    }
    let (med, p95, max) = (pct(&mut st.all, 0.5), pct(&mut st.all, 0.95), pct(&mut st.all, 1.0));
    println!("  all   {med:>8.4}  {p95:>8.4}  {max:>8.3}");
    for (lbl, v) in [("0-1500", &mut st.default_early), ("1700-2400", &mut st.default_late)] {
        println!(
            "  API default 25, {lbl}: median {:.4} p95 {:.4} max {:.3}",
            pct(v, 0.5),
            pct(v, 0.95),
            pct(v, 1.0)
        );
    }
    let c: Vec<String> = st.compat.iter().map(|(n, m)| format!("{n} {m:.2}″")).collect();
    println!("  API-compat names, max over 0-2400: {}", c.join(", "));
    let (dmed, dp95, dmax) =
        (pct(&mut st.dist_rel, 0.5), pct(&mut st.dist_rel, 0.95), pct(&mut st.dist_rel, 1.0));
    let q = |v: &mut Vec<f64>| {
        format!("median {:.1e} p95 {:.1e} max {:.1e}", pct(v, 0.5), pct(v, 0.95), pct(v, 1.0))
    };
    println!("  speed: max |Δ lon speed| vs Swiss-reported {:.2e} ″/day", st.max_speed);
    println!(
        "  on-sky rate vs Swiss's own positions differentiated over ±0.05 d (1900/2000/2100), ″/day: ours {}; Swiss-reported {}",
        q(&mut st.speed_fd),
        q(&mut st.swiss_speed_fd)
    );
    println!("  |Δdist|/dist (parallax sources): median {dmed:.1e} p95 {dp95:.1e} max {dmax:.2}; max RA/Dec separation {:.3}″", st.max_equ);
    st.worst.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut seen = std::collections::HashSet::new();
    let w: Vec<String> = st
        .worst
        .iter()
        .filter(|w| seen.insert(w.1.clone()))
        .take(10)
        .map(|w| format!("{} {:.2}″@{}", w.1, w.0, w.2))
        .collect();
    println!("  worst: {}", w.join(", "));
    (med, p95, max)
}

#[test]
fn apparent_places_vs_swiss_oracle() {
    require_oracle!("stars/swiss.json");
    let Some(path) = pack_path() else {
        eprintln!("SKIP: no full pack (fepk-de441-de442-v1.bin) found");
        return;
    };
    let fx: SwissFixture = serde_json::from_str(&fixture("swiss.json")).unwrap();
    println!("Swiss {} oracle (IAU names)", fx.swe_version);
    let eng = Engine::open(&path).unwrap();

    // Same object: every IAU name in the fixture is < 30″ from Swiss's at J2000
    // (a different star would be degrees away).
    for r in fx.rows.iter().filter(|r| r.year == 2000) {
        let (full, p) = eng.fixstar_ut(&r.name, r.jd_ut, &CalcOptions::default()).unwrap();
        let sep = (diff_deg(p.lon, r.ecl[0]) * (r.ecl[1] * DEG).cos()).hypot(p.lat - r.ecl[1]) * 3600.0;
        assert!(sep < 30.0, "{}: ours {full} is {sep:.0}″ from Swiss", r.name);
    }
    println!(
        "  {} IAU names, all the same object as in Swiss",
        fx.rows.iter().filter(|r| r.year == 2000).count()
    );

    let (med, p95, _) = report("API frame model", &mut compare(&eng, &fx));
    let compat = Engine::with_model(
        fe_core::pack::Pack::open(&path).unwrap(),
        fe_core::frames::FrameModel::SWISS_COMPAT,
    );
    let mut st = compare(&compat, &fx);
    report("SWISS_COMPAT frame model", &mut st);
    // (The max is a star 0.3° from the Sun, where deflection is too curved for a ±0.05 d difference.)
    assert!(pct(&mut st.speed_fd, 0.95) < 1e-4, "speed vs differentiated Swiss positions");

    // Catalogue/model differences, not pipeline errors: see docs/STARS_REPORT.md.
    let mut at2000 = st.by_year[&2000].clone();
    assert!(pct(&mut at2000, 0.5) < 0.01, "median at 2000");
    assert!(med < 0.1 && p95 < 5.0, "median {med} p95 {p95}");
}

#[test]
fn batch_matches_single_and_timing() {
    let Some(path) = pack_path() else {
        eprintln!("SKIP: no pack");
        return;
    };
    let eng = Engine::open(&path).unwrap();
    let names = [
        "Algol",
        "Alcyone",
        "Aldebaran",
        "Rigel",
        "Betelgeuse",
        "Sirius",
        "Canopus",
        "Procyon",
        "Castor",
        "Pollux",
        "Regulus",
        "Denebola",
        "Alkaid",
        "Spica",
        "Arcturus",
        "Agena",
        "Toliman",
        "Unukalhai",
        "Antares",
        "Vega",
        "Altair",
        "Deneb",
        "Fomalhaut",
        "Scheat",
        "Markab",
    ];
    let opt = CalcOptions::default();
    let jd = 2460000.5;
    assert_eq!(names, API_DEFAULT);
    let batch = eng.fixstars_ut(&names, jd, &opt);
    for (n, b) in names.iter().zip(&batch) {
        let (s, p) = eng.fixstar_ut(n, jd, &opt).unwrap();
        let (bs, bp) = b.as_ref().unwrap();
        assert_eq!(&s, bs);
        assert_eq!(p, *bp);
    }
    assert!(eng.fixstars_ut(&["nosuchstar"], jd, &opt)[0].is_err());

    // Sidereal: tropical longitude minus the engine's ayanamsa (true equinox).
    let sid = CalcOptions { sidereal: Some(fe_core::ayanamsa::Ayanamsa::Lahiri), ..opt };
    let (_, t) = eng.fixstar_ut("Spica", jd, &opt).unwrap();
    let (_, sp) = eng.fixstar_ut("Spica", jd, &sid).unwrap();
    let ep = eng.prepare_ut(jd, false).epoch;
    let aya = eng.ayanamsa_deg(fe_core::ayanamsa::Ayanamsa::Lahiri, &ep, true).unwrap();
    assert!(diff_deg(sp.lon, t.lon - aya).abs() < 1e-12, "{} {} {aya}", sp.lon, t.lon);
    let aya_at = |d: f64| {
        let e = eng.prepare_ut(jd + d, false).epoch;
        eng.ayanamsa_deg(fe_core::ayanamsa::Ayanamsa::Lahiri, &e, true).unwrap()
    };
    let aya_rate = (aya_at(1e-3) - aya_at(-1e-3)) / 2e-3;
    assert!(
        ((t.lon_speed - sp.lon_speed) - aya_rate).abs() < 1e-8,
        "{} {aya_rate}",
        t.lon_speed - sp.lon_speed
    );

    let n = 400;
    let t0 = Instant::now();
    for i in 0..n {
        std::hint::black_box(eng.fixstar_ut("Spica", jd + i as f64 * 0.37, &opt).unwrap());
    }
    let single = t0.elapsed().as_secs_f64() / n as f64;
    let t0 = Instant::now();
    for i in 0..n / 10 {
        std::hint::black_box(eng.fixstars_ut(&names, jd + i as f64 * 0.37, &opt));
    }
    let batch = t0.elapsed().as_secs_f64() / (n / 10) as f64;
    let nospeed = CalcOptions { speed: false, ..opt };
    let t0 = Instant::now();
    for i in 0..n {
        std::hint::black_box(eng.fixstar_ut("Spica", jd + i as f64 * 0.37, &nospeed).unwrap());
    }
    let single_ns = t0.elapsed().as_secs_f64() / n as f64;
    println!(
        "timing ({} build): fixstar_ut with speed {:.1} µs, without speed {:.1} µs; 25-star batch with speed {:.1} µs ({:.2} µs/star)",
        if cfg!(debug_assertions) { "debug" } else { "release" },
        single * 1e6,
        single_ns * 1e6,
        batch * 1e6,
        batch * 1e6 / 25.0
    );
}

/// Oracle entries with no published source we could reproduce (reported in
/// docs/STARS_REPORT.md; counted here, not named).
const N_NO_SOURCE: usize = 2;
/// Catalogue records (by `StarRecord::id`) whose position legitimately differs from
/// the oracle's by more than 1″ at 2000 (reasons in docs/STARS_REPORT.md).
const POSITION_EXCEPTIONS: [&str; 14] = [
    "HIP 25281",       // eta Ori: the oracle's RA for this entry is 18 min off
    "SIMBAD NGC 2632", // cluster / nebula centres differ between catalogues
    "SIMBAD NGC 6475",
    "SIMBAD NGC 1039",
    "SIMBAD NGC 6405",
    "SIMBAD NGC 869",
    "SIMBAD M 42",
    "SIMBAD NAME Ori Trapezium",
    "SYS alf Cen AB", // AB barycentre vs the oracle's system place (2.4″)
    // Close visual binaries (pairs of 1-7″): the oracle uses another component/photocentre.
    "HIP 66249",
    "HIP 110960", // zeta Aqr pair (the oracle: zeta1 Aqr component)
    "HIP 40167",
    "SIMBAD * gam01 Leo",
    "HIP 78727",
];

#[derive(Deserialize)]
struct NamesFixture {
    swe_version: String,
    jd_ut: [f64; 3],
    rows: Vec<(String, String, Option<f64>, f64, f64, f64, f64, f64, f64)>,
}

/// The private compatibility table (`private/stars/star_compat_aliases.tsv`, or
/// `FRENCHEPHEMERIS_STAR_COMPAT`) and its oracle fixture (`$FE_ORACLE_DIR`, see tests/common).
fn private_compat() -> Option<(stars::StarCompat, NamesFixture)> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../private/stars");
    let table = std::env::var("FRENCHEPHEMERIS_STAR_COMPAT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join("star_compat_aliases.tsv"));
    let fx = common::oracle_path("stars/swiss_names.json");
    if !table.exists() || !fx.exists() {
        return None;
    }
    let c = stars::StarCompat::load(&table).unwrap();
    Some((c, serde_json::from_str(&std::fs::read_to_string(fx).unwrap()).unwrap()))
}

#[test]
fn every_swiss_accepted_name_resolves_identically_with_compat_table() {
    let Some((compat, fx)) = private_compat() else {
        eprintln!(
            "SKIP: private star compatibility table (private/stars/) or oracle fixture stars/swiss_names.json ({}) not found",
            common::oracle_dir().display()
        );
        return;
    };
    // Names (no pack needed): identical returned string for every accepted query.
    let mut missing = Vec::new();
    for (q, ret, ..) in &fx.rows {
        match stars::lookup_compat(q, Some(&compat)) {
            Ok((_, full)) => assert_eq!(&full, ret, "query {q:?}"),
            Err(_) => missing.push(ret.clone()),
        }
    }
    for n in API_DEFAULT.iter().chain(COMPAT.iter()) {
        assert!(fx.rows.iter().any(|r| &r.0 == n), "{n} in the oracle fixture");
    }
    missing.sort();
    missing.dedup();
    println!("Swiss {}: {} accepted queries, {} unresolved", fx.swe_version, fx.rows.len(), missing.len());
    assert_eq!(missing.len(), N_NO_SOURCE, "unresolved oracle entries");

    // Magnitudes: Johnson V (SIMBAD) vs Swiss's file magnitudes.
    let mut dm = Vec::new();
    let mut worst_m = Vec::new();
    for (q, ret, m, ..) in &fx.rows {
        if let (Some(m), Ok((s, _))) = (m, stars::lookup_compat(q, Some(&compat))) {
            let ours = s.magnitude();
            if *m < 900.0 && ours < 900.0 {
                dm.push((ours - m).abs());
                worst_m.push(((ours - m).abs(), ret.clone(), ours, *m));
            }
        }
    }
    worst_m.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    worst_m.dedup_by(|a, b| a.1 == b.1);
    println!(
        "  |Δmag|: n {} median {:.3} p95 {:.3} max {:.2}; > 0.1: {}; worst {:?}",
        dm.len(),
        pct(&mut dm.clone(), 0.5),
        pct(&mut dm.clone(), 0.95),
        pct(&mut dm.clone(), 1.0),
        dm.iter().filter(|d| **d > 0.1).count(),
        &worst_m[..8.min(worst_m.len())]
    );

    // Positions (needs the pack).
    let Some(path) = pack_path() else {
        eprintln!("SKIP positions: no full pack");
        return;
    };
    let eng = Engine::open(&path).unwrap().with_star_compat_table(Some(std::sync::Arc::new(compat)));
    let opt = CalcOptions { speed: false, ..CalcOptions::default() };
    let mut seps: [Vec<f64>; 3] = [vec![], vec![], vec![]];
    let mut over: Vec<(f64, String)> = Vec::new();
    let mut far: Vec<(f64, String, usize)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for (q, ret, _m, l0, b0, l1, b1, l2, b2) in &fx.rows {
        let Ok((star, _)) = eng.fixstar_lookup(q) else { continue };
        if !seen.insert(ret.clone()) {
            continue;
        }
        for (k, (l, b)) in [(*l0, *b0), (*l1, *b1), (*l2, *b2)].into_iter().enumerate() {
            let (_, p) = eng.fixstar_ut(q, fx.jd_ut[k], &opt).unwrap();
            let sep = (diff_deg(p.lon, l) * (b * DEG).cos()).hypot(p.lat - b) * 3600.0;
            if POSITION_EXCEPTIONS.contains(&star.record.id) {
                if k == 1 {
                    over.push((sep, format!("{ret} (exception)")));
                }
                continue;
            }
            seps[k].push(sep);
            far.push((sep, ret.clone(), k));
            if k == 1 && sep > 1.0 {
                over.push((sep, ret.clone()));
            }
        }
    }
    for (k, y) in ["0", "2000", "2400"].iter().enumerate() {
        let v = &mut seps[k];
        println!(
            "  year {y}: {} objects, median {:.4}″ p95 {:.4}″ max {:.3}″",
            v.len(),
            pct(v, 0.5),
            pct(v, 0.95),
            pct(v, 1.0)
        );
    }
    far.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    println!("  worst (sep, name, 0=year0/2=year2400): {:?}", &far[..15.min(far.len())]);
    over.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    println!("  > 1″ at 2000 and exceptions: {over:?}");
    let unexplained: Vec<_> = over.iter().filter(|o| !o.1.ends_with("(exception)")).collect();
    assert!(unexplained.is_empty(), "unexplained > 1″ at 2000: {unexplained:?}");
}

/// No catalogue value may come from Gaia (CC BY-SA 3.0 IGO): the provenance of every
/// record (`rv_ref`) is checked against the Gaia releases (VizieR I/337, I/345-I/358) and
/// the Gaia-derived publications blocked by tools/build_stars.py.
#[test]
fn no_gaia_derived_values_in_catalogue() {
    const GAIA_BIBCODES: [&str; 25] = [
        "2016A&A...595A...1G",
        "2016A&A...595A...2G",
        "2018A&A...616A...1G",
        "2021A&A...649A...1G",
        "2023A&A...674A...1G",
        "2023A&A...674A...5K",
        "2019A&A...622A.205K",
        "2018A&A...616A...5S",
        "2018A&A...616A...6K",
        "2018A&A...616A...7S",
        "2018A&A...616A..10G",
        "2018A&A...618A..93C",
        "2020A&A...640A...1C",
        "2021A&A...647A..19T",
        "2019MNRAS.482.5138B",
        "2019MNRAS.483.5026L",
        "2019ApJ...870...32K",
        "2021MNRAS.504..356D",
        "2021MNRAS.505.5978V",
        "2023A&A...674A..37G",
        "2023A&A...673A.114H",
        "2021A&A...649A...6G",
        "2021A&A...649A...7G",
        "2021MNRAS.506..150B",
        "2018MNRAS.473..849D",
    ];
    let vizier_gaia = |r: &str| {
        r.match_indices("yCat.13").any(|(i, _)| {
            let n: String = r[i + 7..].chars().take(2).collect();
            matches!(
                n.as_str(),
                "37" | "45"
                    | "46"
                    | "47"
                    | "48"
                    | "49"
                    | "50"
                    | "51"
                    | "52"
                    | "53"
                    | "54"
                    | "55"
                    | "56"
                    | "57"
                    | "58"
            )
        })
    };
    for r in stars::catalog() {
        assert!(!r.rv_ref.to_lowercase().contains("gaia"), "{}: {}", r.id, r.rv_ref);
        assert!(!vizier_gaia(r.rv_ref), "{}: {}", r.id, r.rv_ref);
        for b in GAIA_BIBCODES {
            assert!(!r.rv_ref.contains(b), "{}: Gaia-derived source {b}", r.id);
        }
    }
}
