// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Evaluate a pack for many epochs (used by harness/check_pack_vs_jpl.py).
//!
//! `pack_eval <pack[,pack2,...]> <target> <center> <times.f64le> <out.f64le>`
//! (several comma-separated packs are opened together with `Pack::open_many`)
//! Input: raw little-endian f64 JD TDB values. Output: 6 f64 per epoch
//! (x y z km, vx vy vz km/day); NaN when the epoch is out of range.

use fe_core::pack::Pack;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 6 {
        eprintln!("usage: pack_eval <pack> <target> <center> <times.f64le> <out.f64le>");
        std::process::exit(2);
    }
    let paths: Vec<&str> = a[1].split(',').collect();
    let pack = Pack::open_many(&paths).expect("open pack");
    let (t, c): (i32, i32) = (a[2].parse().unwrap(), a[3].parse().unwrap());
    let raw = std::fs::read(&a[4]).expect("read times");
    let mut out = Vec::with_capacity(raw.len() * 6);
    for ch in raw.as_chunks::<8>().0 {
        let jd = f64::from_le_bytes(*ch);
        let (p, v) = pack.state(t, c, jd).unwrap_or(([f64::NAN; 3], [f64::NAN; 3]));
        for x in p.iter().chain(v.iter()) {
            out.extend_from_slice(&x.to_le_bytes());
        }
    }
    std::fs::write(&a[5], out).expect("write output");
}
