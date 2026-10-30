// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! cargo run --release -p fe-core --example bench_events [pack]
use fe_core::context::{Body, CalcOptions, Engine};
use fe_core::ephem::Observer;
use fe_core::events::*;
use std::time::Instant;

fn bench(label: &str, n: usize, mut f: impl FnMut(usize) -> f64) {
    let mut acc = 0.0;
    for i in 0..n.min(200) {
        acc += f(i); // warm the pack pages
    }
    let t0 = Instant::now();
    for i in 0..n {
        acc += f(i);
    }
    let dt = t0.elapsed().as_secs_f64();
    println!("{label:34} {:8.2} us/call  (chk {acc:.3})", dt / n as f64 * 1e6);
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or("data/packs/fepk-de441-de442-v1.bin".into());
    let e = Engine::open(path.as_ref()).unwrap();
    let obs = Observer { lon_deg: 77.21, lat_deg: 28.61, height_m: 0.0 };
    let t0 = 2460000.5;
    let n = 20_000;
    let d = RiseSetOptions::default();
    let h = RiseSetOptions::hindu();
    let rs = |b, k, o: &RiseSetOptions, i: usize| {
        e.rise_set(b, t0 + i as f64 * 0.731, k, &obs, o).unwrap().jd().unwrap_or(0.0)
    };
    bench("rise Sun (default)", n, |i| rs(Body::Sun, RiseSetKind::Rise, &d, i));
    bench("set Sun (default)", n, |i| rs(Body::Sun, RiseSetKind::Set, &d, i));
    bench("rise Sun (hindu)", n, |i| rs(Body::Sun, RiseSetKind::Rise, &h, i));
    bench("rise Moon (default)", n, |i| rs(Body::Moon, RiseSetKind::Rise, &d, i));
    bench("rise Jupiter (default)", n, |i| rs(Body::Jupiter, RiseSetKind::Rise, &d, i));
    bench("upper transit Sun", n, |i| rs(Body::Sun, RiseSetKind::UpperTransit, &d, i));
    let polar = Observer { lon_deg: 18.95, lat_deg: 69.65, height_m: 0.0 };
    bench("rise Sun 69.65N (incl. polar day)", n, |i| {
        e.rise_set(Body::Sun, t0 + i as f64 * 0.731, RiseSetKind::Rise, &polar, &d)
            .unwrap()
            .jd()
            .unwrap_or(0.0)
    });
    let opt = CalcOptions::default();
    bench("pheno_ut Moon", n, |i| e.pheno_ut(Body::Moon, t0 + i as f64 * 0.731, &opt).unwrap().phase);
    bench("pheno_ut Saturn", n, |i| e.pheno_ut(Body::Saturn, t0 + i as f64 * 0.731, &opt).unwrap().magnitude);
    bench("solcross_ut", n, |i| {
        e.solcross_ut((i as f64 * 37.0) % 360.0, t0 + i as f64 * 0.731, &opt).unwrap()
    });
    bench("mooncross_ut", n, |i| {
        e.mooncross_ut((i as f64 * 37.0) % 360.0, t0 + i as f64 * 0.731, &opt).unwrap()
    });
    bench("time_equ", n, |i| e.time_equ(t0 + i as f64 * 0.731).unwrap());
    bench("cotrans", n, |i| cotrans([i as f64, 1.0, 1.0], -23.44)[0]);
    let m = e.frame_model();
    bench("  Epoch::from_ut_with (one instant)", n, |i| {
        fe_core::ephem::Epoch::from_ut_with(t0 + i as f64 * 0.731, m).jd_tt
    });
    bench("  FrameSet::with_model", n, |i| {
        fe_core::frames::FrameSet::with_model(t0 + i as f64 * 0.731, m).eo
    });
    bench("  delta_t", n, |i| fe_core::time::delta_t_seconds(t0 + i as f64 * 0.731));
    let ei = EpochInterp::new(t0, m);
    bench("  EpochInterp::at", n, |i| ei.at(t0 + (i % 100) as f64 * 0.012).jd_tt);
    bench("  topocentric place (interp epoch)", n, |i| {
        let ep = ei.at(t0 + (i % 100) as f64 * 0.012);
        e.rise_set_place_epoch(Body::Sun, &ep, &obs, &d).unwrap().0
    });
}
