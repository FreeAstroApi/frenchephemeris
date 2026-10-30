// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Event searches: bracketed root finding on smooth functions of time.
//!
//! Everything reduces to "find t where g(t) = 0" with g continuous between
//! samples: longitude crossings (ingresses, returns, transits to a natal
//! point), separation crossings (aspects, lunations) and speed zeros
//! (stations). Brackets come from a uniform scan; each root is refined with
//! Brent's method to `TOL_DAYS`.

use crate::context::{Body, CalcOptions, Engine};

use crate::math::diff_deg;
use crate::Result;

/// Root tolerance, days (≈ 0.9 ms).
pub const TOL_DAYS: f64 = 1e-8;

/// Brent's method on [a, b] with g(a)·g(b) ≤ 0.
pub fn brent(
    mut g: impl FnMut(f64) -> Result<f64>,
    mut a: f64,
    mut b: f64,
    mut fa: f64,
    mut fb: f64,
) -> Result<f64> {
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    for _ in 0..100 {
        if fb.signum() == fc.signum() {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * f64::EPSILON * b.abs() + 0.5 * TOL_DAYS;
        let m = 0.5 * (c - b);
        if m.abs() <= tol || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                let qa = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * m * qa * (qa - r) - (b - a) * (r - 1.0));
                q = (qa - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            if 2.0 * p < (3.0 * m * q - (tol * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = d;
            }
        } else {
            d = m;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol { d } else { tol.copysign(m) };
        fb = g(b)?;
    }
    Ok(b)
}

/// All roots of `g` in [t0, t1] found by scanning with `step` days.
///
/// `max_jump` rejects sign changes where |g| is large on both sides — the
/// wrap-around discontinuity of an angular difference (±180°), not a root.
pub fn scan_roots(
    mut g: impl FnMut(f64) -> Result<f64>,
    t0: f64,
    t1: f64,
    step: f64,
    max_jump: f64,
) -> Result<Vec<f64>> {
    let mut roots = Vec::new();
    let mut a = t0;
    let mut fa = g(a)?;
    while a < t1 {
        let b = (a + step).min(t1);
        let fb = g(b)?;
        if fa == 0.0 {
            roots.push(a);
        } else if fa.signum() != fb.signum() && fb != 0.0 && (fa - fb).abs() < max_jump {
            roots.push(brent(&mut g, a, b, fa, fb)?);
        }
        a = b;
        fa = fb;
    }
    if fa == 0.0 {
        roots.push(t1);
    }
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-6);
    Ok(roots)
}

/// Roots of an angular function by speed-bounded stepping.
///
/// If |dg/dt| ≤ `vmax` everywhere, g cannot reach zero sooner than |g|/vmax
/// after the current sample, so the scan jumps that far (minus nothing — the
/// bound is strict) and only takes fixed `min_step` steps near a root. No root
/// is skipped as long as two roots are never closer than `min_step`.
pub fn scan_roots_bounded(
    mut g: impl FnMut(f64) -> Result<f64>,
    t0: f64,
    t1: f64,
    vmax: f64,
    min_step: f64,
    max_jump: f64,
) -> Result<Vec<f64>> {
    let mut roots = Vec::new();
    let mut a = t0;
    let mut fa = g(a)?;
    while a < t1 {
        let step = (fa.abs() / vmax).max(min_step);
        let b = (a + step).min(t1);
        let fb = g(b)?;
        if fa == 0.0 {
            roots.push(a);
        } else if fa.signum() != fb.signum() && fb != 0.0 && (fa - fb).abs() < max_jump {
            roots.push(brent(&mut g, a, b, fa, fb)?);
        }
        a = b;
        fa = fb;
    }
    if fa == 0.0 {
        roots.push(t1);
    }
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-6);
    Ok(roots)
}

/// Upper bound of |d lon / dt| (deg/day, geocentric apparent) with margin.
pub fn max_speed(body: Body) -> f64 {
    match body {
        Body::Moon => 16.5,
        Body::Mercury => 2.5,
        Body::Venus => 1.5,
        Body::Sun => 1.1,
        Body::Mars => 0.9,
        Body::Jupiter => 0.3,
        Body::Saturn => 0.15,
        Body::Uranus | Body::Neptune | Body::Pluto => 0.1,
        Body::MeanNode | Body::MeanApogee => 0.2,
        // Osculating points oscillate quickly.
        Body::TrueNode => 2.0,
        Body::OscuApogee => 30.0,
        Body::Earth => 1.1,
        // Minor bodies: max |dλ/dt| over their JPL ranges is checked in
        // harness/check_minor.py: measured 0.69 (Pallas), 0.95 (Eros), 0.15
        // (Chiron), 0.04°/d (Eris); bounds keep >= 30 % margin.
        Body::Ceres | Body::Pallas | Body::Juno | Body::Vesta => 1.0,
        Body::Eros => 3.5,
        Body::Chiron | Body::Pholus => 0.2,
        Body::Eris => 0.1,
        // Unknown orbit class: bound for near-Earth objects at ~0.05 AU.
        Body::Asteroid(_) => 20.0,
        // natural apsides: measured |dλ/dt| < 0.2 (apogee), < 1.0°/d (perigee)
        Body::IntpApogee | Body::IntpPerigee => 2.0,
        // hypothetical bodies: outer ones move like Neptune or slower; the
        // geocentric Waltemath moon ~3.0, Selena ~0.15°/d.
        Body::Waldemath => 4.0,
        Body::WhiteMoon => 0.3,
        Body::Cupido
        | Body::Hades
        | Body::Zeus
        | Body::Kronos
        | Body::Apollon
        | Body::Admetos
        | Body::Vulkanus
        | Body::Poseidon
        | Body::Harrington
        | Body::NeptuneLeverrier
        | Body::NeptuneAdams
        | Body::PlutoLowell => 0.1,
    }
}

/// Default scan step (days) for a body: small enough that no two crossings
/// of the same longitude fit in one step.
pub fn default_step(body: Body) -> f64 {
    match body {
        Body::Moon => 0.5,
        Body::Mercury | Body::Venus | Body::Sun | Body::TrueNode | Body::OscuApogee => 2.0,
        Body::Mars | Body::Eros => 4.0,
        _ => 8.0,
    }
}

fn lon_tt(engine: &Engine, body: Body, jd_tt: f64, opt: &CalcOptions) -> Result<(f64, f64)> {
    let p = engine.calc_tt(body, jd_tt, opt)?;
    Ok((p.lon, p.lon_speed))
}

/// Times (JD TT) when `body` reaches ecliptic longitude `target_deg`:
/// sign ingresses, returns, transits to a fixed natal point.
pub fn longitude_crossings(
    engine: &Engine,
    body: Body,
    target_deg: f64,
    t0_tt: f64,
    t1_tt: f64,
    opt: &CalcOptions,
) -> Result<Vec<f64>> {
    let opt = CalcOptions { speed: false, xyz: false, radians: false, ..*opt };
    let min_step = default_step(body) / 4.0;
    scan_roots_bounded(
        |t| Ok(diff_deg(lon_tt(engine, body, t, &opt)?.0, target_deg)),
        t0_tt,
        t1_tt,
        max_speed(body),
        min_step,
        180.0,
    )
}

/// Times when the longitude difference `lon(a) − lon(b)` equals `angle_deg`
/// (exact aspects between moving bodies; lunations with Sun/Moon).
pub fn separation_crossings(
    engine: &Engine,
    a: Body,
    b: Body,
    angle_deg: f64,
    t0_tt: f64,
    t1_tt: f64,
    opt: &CalcOptions,
) -> Result<Vec<f64>> {
    let opt = CalcOptions { speed: false, xyz: false, radians: false, ..*opt };
    let min_step = default_step(a).min(default_step(b)) / 4.0;
    scan_roots_bounded(
        |t| {
            let ep = engine.prepare_tt(t, false);
            let la = engine.calc_prepared(a, &ep, &opt)?.lon;
            let lb = engine.calc_prepared(b, &ep, &opt)?.lon;
            Ok(diff_deg(la - lb, angle_deg))
        },
        t0_tt,
        t1_tt,
        max_speed(a) + max_speed(b),
        min_step,
        180.0,
    )
}

/// Stations: times when the longitude speed of `body` changes sign.
/// Returns (jd_tt, turning_retrograde).
pub fn stations(
    engine: &Engine,
    body: Body,
    t0_tt: f64,
    t1_tt: f64,
    opt: &CalcOptions,
) -> Result<Vec<(f64, bool)>> {
    let opt = CalcOptions { speed: true, xyz: false, radians: false, ..*opt };
    let roots = scan_roots(|t| Ok(lon_tt(engine, body, t, &opt)?.1), t0_tt, t1_tt, 1.0, f64::INFINITY)?;
    roots.into_iter().map(|t| Ok((t, lon_tt(engine, body, t + 0.01, &opt)?.1 < 0.0))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brent_finds_cosine_zero() {
        let r = scan_roots(|t| Ok(t.cos()), 0.0, 10.0, 0.7, f64::INFINITY).unwrap();
        let expect = [
            std::f64::consts::FRAC_PI_2,
            3.0 * std::f64::consts::FRAC_PI_2,
            5.0 * std::f64::consts::FRAC_PI_2,
        ];
        assert_eq!(r.len(), 3);
        for (a, b) in r.iter().zip(expect) {
            assert!((a - b).abs() < 1e-7, "{a} vs {b}");
        }
    }

    #[test]
    fn wrap_is_not_a_root() {
        // Sawtooth angle difference jumping from +179 to −179: no root.
        let r = scan_roots(|t| Ok(diff_deg(t * 50.0, 180.0)), 3.0, 4.0, 0.1, 180.0).unwrap();
        assert!(r.iter().all(|t| (t * 50.0 - 180.0).abs() < 1e-5), "{r:?}");
    }
}
