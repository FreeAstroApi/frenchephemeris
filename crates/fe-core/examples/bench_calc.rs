// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! cargo run --release -p fe-core --example bench_calc [pack]
use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::frames::FrameModel;
use fe_core::pack::Pack;
use std::time::Instant;

const BODIES: [Body; 10] = [
    Body::Sun,
    Body::Moon,
    Body::Mercury,
    Body::Venus,
    Body::Mars,
    Body::Jupiter,
    Body::Saturn,
    Body::Uranus,
    Body::Neptune,
    Body::Pluto,
];

fn main() {
    let path = std::env::args().nth(1).unwrap_or("data/packs/fepk-de441-de442-v1.bin".into());
    for (name, model) in [
        ("DEFAULT (Vondrak+2000A)", FrameModel::DEFAULT),
        ("API_FRAME_MODEL", fe_core::context::API_FRAME_MODEL),
    ] {
        let eng = Engine::with_model(Pack::open(path.as_ref()).unwrap(), model);
        let n = 20_000;
        let times: Vec<f64> = (0..n).map(|i| 2440000.5 + i as f64 * 1.37).collect();
        for (label, speed) in [("pos+speed", true), ("pos only", false)] {
            let opt = CalcOptions { speed, ..Default::default() };
            let t0 = Instant::now();
            let mut acc = 0.0;
            for r in eng.calc_batch_ut(&BODIES, &times, &opt) {
                acc += r.unwrap().lon;
            }
            let dt = t0.elapsed().as_secs_f64();
            println!(
                "{name:24} {label:9}: {:7.3} us per 10-body chart, {:6.3} us per body  (chk {acc:.1})",
                dt / n as f64 * 1e6,
                dt / (n * 10) as f64 * 1e6
            );
        }
    }
}
