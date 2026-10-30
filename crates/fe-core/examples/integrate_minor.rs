// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Test-particle integrator for minor bodies (used by tools/integrate_minor.py).
//!
//!     cargo run -p fe-core --example integrate_minor -- config.json
//!
//! Integrates one or more massless particles (barycentric ICRF, km, km/day,
//! JD TDB) from `t0` to `t1` with a fixed step `h` (sign taken from t1 - t0),
//! perturbed by point masses whose barycentric positions come from FEPK packs
//! (the main DE441/DE442 pack, optionally the minor pack's big asteroids), plus
//! the 1PN (Schwarzschild, PPN beta = gamma = 1) term of the Sun.
//!
//! Method: 6-stage Gauss-Legendre collocation (implicit Runge-Kutta, order 12,
//! symplectic), stage equations solved by fixed-point iteration to round-off,
//! compensated (Kahan) summation of the state update. Coefficients are
//! computed at start-up from the Legendre roots.
//!
//! config.json:
//! {"packs": [...], "t0": f, "t1": f, "h": f,
//!  "states": [[x,y,z,vx,vy,vz], ...],          // clone 0 = nominal
//!  "perturbers": [[naif, gm_km3_day2], ...],   // includes the Sun (10)
//!  "gm_sun": f, "c_km_day": f, "gr": bool,
//!  "out": path, "every": k}                    // clones >= 1 every k steps
//! Output (little-endian f64): clone 0 at every step: t x y z vx vy vz ax ay az;
//! then each further clone every `every` steps (and the last): t x y z vx vy vz.

use fe_core::pack::Pack;
use std::io::Write;

const S: usize = 6;

fn legendre_roots(n: usize) -> Vec<(f64, f64)> {
    // (root, weight) on [-1, 1], Newton from Chebyshev-like guesses
    let mut out = Vec::new();
    for i in 1..=n {
        let mut x = (std::f64::consts::PI * (i as f64 - 0.25) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-17 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

/// Gauss-Legendre collocation tableau (c, A, b) on [0, 1].
fn tableau() -> ([f64; S], [[f64; S]; S], [f64; S]) {
    let mut c = [0.0; S];
    let mut r = legendre_roots(S);
    r.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for i in 0..S {
        c[i] = 0.5 * (r[i].0 + 1.0);
    }
    let lag = |j: usize, t: f64| -> f64 {
        let mut v = 1.0;
        for m in 0..S {
            if m != j {
                v *= (t - c[m]) / (c[j] - c[m]);
            }
        }
        v
    };
    let q = legendre_roots(8); // exact for degree <= 15
    let integ = |j: usize, upper: f64| -> f64 {
        q.iter().map(|&(x, w)| 0.5 * upper * w * lag(j, 0.5 * upper * (x + 1.0))).sum()
    };
    let mut a = [[0.0; S]; S];
    let mut b = [0.0; S];
    for i in 0..S {
        for j in 0..S {
            a[i][j] = integ(j, c[i]);
        }
    }
    for j in 0..S {
        b[j] = integ(j, 1.0);
    }
    (c, a, b)
}

#[derive(serde::Deserialize)]
struct Config {
    packs: Vec<String>,
    t0: f64,
    t1: f64,
    h: f64,
    states: Vec<[f64; 6]>,
    perturbers: Vec<(i32, f64)>,
    gm_sun: f64,
    c_km_day: f64,
    gr: bool,
    out: String,
    every: usize,
    /// Optional constant offsets (km) added to perturber positions, for
    /// sensitivity experiments: [[naif, dx, dy, dz], ...].
    #[serde(default)]
    offsets: Vec<(i32, f64, f64, f64)>,
}

struct Model {
    pack: Pack,
    pert: Vec<(i32, f64)>,
    offsets: Vec<(i32, f64, f64, f64)>,
    gm_sun: f64,
    c2: f64,
    gr: bool,
}

impl Model {
    /// Accelerations for all clones at time t.
    fn accel(&self, t: f64, ys: &[[f64; 6]], out: &mut [[f64; 3]]) {
        let mut bodies = Vec::with_capacity(self.pert.len());
        for &(id, gm) in &self.pert {
            let mut p = self.pack.position(id, 0, t).unwrap_or_else(|e| panic!("perturber {id} at {t}: {e}"));
            for &(oid, dx, dy, dz) in &self.offsets {
                if oid == id {
                    p = [p[0] + dx, p[1] + dy, p[2] + dz];
                }
            }
            bodies.push((p, gm));
        }
        let sun = if self.gr { Some(self.pack.state(10, 0, t).unwrap()) } else { None };
        for (y, a) in ys.iter().zip(out.iter_mut()) {
            let mut acc = [0.0; 3];
            for &(p, gm) in &bodies {
                let d = [y[0] - p[0], y[1] - p[1], y[2] - p[2]];
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                let f = -gm / (r2 * r2.sqrt());
                for k in 0..3 {
                    acc[k] += f * d[k];
                }
            }
            if let Some((sp, sv)) = sun {
                let r = [y[0] - sp[0], y[1] - sp[1], y[2] - sp[2]];
                let v = [y[3] - sv[0], y[4] - sv[1], y[5] - sv[2]];
                let rn = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                let rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
                let k0 = self.gm_sun / (self.c2 * rn * rn * rn);
                let a1 = 4.0 * self.gm_sun / rn - v2;
                for k in 0..3 {
                    acc[k] += k0 * (a1 * r[k] + 4.0 * rv * v[k]);
                }
            }
            *a = acc;
        }
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("config.json");
    let cfg: Config = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let pack = Pack::open_many(&cfg.packs).unwrap();
    let m = Model {
        pack,
        pert: cfg.perturbers.clone(),
        offsets: cfg.offsets.clone(),
        gm_sun: cfg.gm_sun,
        c2: cfg.c_km_day * cfg.c_km_day,
        gr: cfg.gr,
    };
    let (c, a, b) = tableau();
    let nsteps = ((cfg.t1 - cfg.t0) / cfg.h).abs().round() as usize;
    let h = (cfg.t1 - cfg.t0) / nsteps as f64;
    let nc = cfg.states.len();
    let mut y: Vec<[f64; 6]> = cfg.states.clone();
    let mut comp: Vec<[f64; 6]> = vec![[0.0; 6]; nc]; // Kahan compensation
                                                      // stage derivatives K[i][clone] (6-vectors)
    let mut kst: Vec<Vec<[f64; 6]>> = vec![vec![[0.0; 6]; nc]; S];
    let mut acc0 = vec![[0.0; 3]; nc];
    let mut f = std::io::BufWriter::new(std::fs::File::create(&cfg.out).unwrap());
    let mut rows_other: Vec<Vec<f64>> = vec![Vec::new(); nc];
    let mut t = cfg.t0;
    let mut max_iter = 0;
    // initial guess for stages: F(t0, y0)
    m.accel(t, &y, &mut acc0);
    for i in 0..S {
        for j in 0..nc {
            kst[i][j] = [y[j][3], y[j][4], y[j][5], acc0[j][0], acc0[j][1], acc0[j][2]];
        }
    }
    let mut ystage = vec![vec![[0.0; 6]; nc]; S];
    let mut astage = vec![[0.0; 3]; nc];
    for step in 0..=nsteps {
        // output at t (state y, acceleration acc0)
        m.accel(t, &y, &mut acc0);
        let r = &y[0];
        for v in [t, r[0], r[1], r[2], r[3], r[4], r[5], acc0[0][0], acc0[0][1], acc0[0][2]] {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
        if step % cfg.every == 0 || step == nsteps {
            for j in 1..nc {
                rows_other[j].push(t);
                rows_other[j].extend_from_slice(&y[j]);
            }
        }
        if step == nsteps {
            break;
        }
        // solve stage equations by fixed-point iteration
        let mut it = 0;
        let mut prev = f64::INFINITY;
        loop {
            it += 1;
            for i in 0..S {
                for j in 0..nc {
                    let mut s = y[j];
                    for l in 0..S {
                        for k in 0..6 {
                            s[k] += h * a[i][l] * kst[l][j][k];
                        }
                    }
                    ystage[i][j] = s;
                }
            }
            let mut delta: f64 = 0.0;
            for i in 0..S {
                m.accel(t + c[i] * h, &ystage[i], &mut astage);
                for j in 0..nc {
                    let ys = &ystage[i][j];
                    let new = [ys[3], ys[4], ys[5], astage[j][0], astage[j][1], astage[j][2]];
                    for k in 0..6 {
                        let scale = if k < 3 {
                            ys[3].abs() + ys[4].abs() + ys[5].abs()
                        } else {
                            astage[j][0].abs() + astage[j][1].abs() + astage[j][2].abs()
                        };
                        delta = delta.max((new[k] - kst[i][j][k]).abs() / scale.max(1e-300));
                    }
                    kst[i][j] = new;
                }
            }
            // converged to round-off, or no longer contracting (round-off floor)
            if delta < 1e-15 || (it > 3 && delta >= 0.5 * prev) || it >= 30 {
                break;
            }
            prev = delta;
        }
        max_iter = max_iter.max(it);
        for j in 0..nc {
            for k in 0..6 {
                let mut inc = 0.0;
                for i in 0..S {
                    inc += h * b[i] * kst[i][j][k];
                }
                // Kahan
                let yk = inc - comp[j][k];
                let tt = y[j][k] + yk;
                comp[j][k] = (tt - y[j][k]) - yk;
                y[j][k] = tt;
            }
        }
        t = cfg.t0 + (step + 1) as f64 * h;
    }
    for j in 1..nc {
        for v in &rows_other[j] {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
    }
    f.flush().unwrap();
    eprintln!("integrated {nc} particle(s), {nsteps} steps of {h} d, max {max_iter} stage iterations");
}
