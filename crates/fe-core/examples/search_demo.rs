// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Prints (as JSON) equinoxes, new moons and Mercury stations 2000–2030 with timings.
use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::search;
use std::time::Instant;

fn main() {
    let eng = Engine::open("data/packs/fepk-de441-de442-v1.bin".as_ref()).unwrap();
    let opt = CalcOptions::default();
    let (t0, t1) = (2451544.5, 2462502.5); // 2000-01-01 .. 2030-01-01 (TT)
    let s = Instant::now();
    let eq = search::longitude_crossings(&eng, Body::Sun, 0.0, t0, t1, &opt).unwrap();
    let d_eq = s.elapsed().as_secs_f64();
    let s = Instant::now();
    let nm = search::separation_crossings(&eng, Body::Moon, Body::Sun, 0.0, t0, t1, &opt).unwrap();
    let d_nm = s.elapsed().as_secs_f64();
    let s = Instant::now();
    let st = search::stations(&eng, Body::Mercury, t0, t1, &opt).unwrap();
    let d_st = s.elapsed().as_secs_f64();
    eprintln!(
        "30 years: {} equinoxes in {:.1} ms, {} new moons in {:.1} ms, {} Mercury stations in {:.1} ms",
        eq.len(),
        d_eq * 1e3,
        nm.len(),
        d_nm * 1e3,
        st.len(),
        d_st * 1e3
    );
    let st: Vec<f64> = st.iter().map(|x| x.0).collect();
    println!("{{\"equinox_tt\":{eq:?},\"newmoon_tt\":{nm:?},\"mercury_station_tt\":{st:?}}}");
}
