// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Natural ("interpolated") lunar apogee and perigee — Swiss ids 21/22
//! (`SE_INTP_APOG`, natural Lilith; `SE_INTP_PERG`, Priapus).
//!
//! Definition (Swiss Ephemeris documentation §2.2.4; D. Koch, "Was ist
//! Lilith und welche Ephemeride ist richtig", Meridian 1/95; Koch & Rindgen,
//! "Lilith und Priapus", Frankfurt 2000; J. Revilla, "The Astronomical
//! Variants of the Lunar Apogee — Black Moon"): the natural apogee (perigee)
//! is the place of the Moon at its actual passages through apogee (perigee);
//! between passages it follows a smooth curve with a continuous derivative,
//! which Koch "derived from the analytical lunar theory" rather than by a
//! numerical interpolation of the passages (§2.2.4, "Note on implementation").
//!
//! Construction: for an instant t, keep the arguments of the lunar theory
//! that do not belong to the Moon's own motion (Sun, planets, and the
//! secular parts of perigee and node) at t, and move only the Moon along
//! its orbit (W1 → W1 + Δ: D, F, l and ζ shift by Δ) to the apsis. Δ solves
//! dr/dt = Σ ∂r/∂argᵢ · ṡᵢ = 0 near l + Δ = 0 (perigee) or π (apogee), with
//! the true rates ṡᵢ of all arguments — the apsis the Moon would pass if it
//! were there now. Its longitude, latitude and distance come from ELP 2000-82B
//! (main problem and perturbations, [`super::elp_main`]). When the Moon is
//! actually at an apsis Δ = 0 and the point is the Moon itself: at the
//! JPL DE441 apsis passages the curve is within ~0.001–0.004° of the JPL
//! Moon (docs/FICTITIOUS_BODIES.md).
//!
//! Coordinates: mean ecliptic and equinox of date (geometric); the engine
//! adds Δψ for the true equinox like for the other lunar points.

use super::elp_main as E;
use crate::math::*;

#[inline]
fn poly(p: &(f64, [f64; 4]), t: f64) -> f64 {
    p.0 + t * (p.1[0] + t * (p.1[1] + t * (p.1[2] + t * p.1[3])))
}

#[inline]
fn poly_rate(p: &(f64, [f64; 4]), t: f64) -> f64 {
    p.1[0] + t * (2.0 * p.1[1] + t * (3.0 * p.1[2] + t * 4.0 * p.1[3]))
}

/// Arguments shifted by moving the Moon along its orbit: D, F, l, ζ.
const SHIFT: [f64; 13] = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];

#[inline]
fn dot13(m: &[i8; 13], a: &[f64; 13]) -> f64 {
    m.iter().zip(a).map(|(&k, &x)| k as f64 * x).sum()
}

#[inline]
fn tpow(t: f64, k: u8) -> f64 {
    match k {
        0 => 1.0,
        1 => t,
        2 => t * t,
        _ => t * t * t,
    }
}

/// The virtual apsis at T (Julian centuries TDB from J2000): (longitude
/// minus the ELP mean perigee W2 (+π for the apogee), rad; latitude rad;
/// distance km).
pub fn theory(t: f64, apogee: bool) -> (f64, f64, f64) {
    // full polynomial arguments (main problem) and their linear part
    // (perturbation series), as in the authors' routine ELP82B_2
    let mut a0 = [0.0; 13];
    let mut rt = [0.0; 13];
    let mut b0 = [0.0; 13];
    let mut bt = [0.0; 13];
    for i in 0..13 {
        a0[i] = poly(&E::ARGS[i], t);
        rt[i] = poly_rate(&E::ARGS[i], t);
        b0[i] = E::ARGS[i].0 + E::ARGS[i].1[0] * t;
        bt[i] = E::ARGS[i].1[0];
    }
    let target = if apogee { PI } else { 0.0 };
    let mut dl = norm_pi(target - a0[2]);
    let at = |dl: f64| {
        let (mut a, mut b) = (a0, b0);
        for i in 0..13 {
            a[i] += dl * SHIFT[i];
            b[i] += dl * SHIFT[i];
        }
        (a, b)
    };
    for _ in 0..12 {
        let (a, b) = at(dl);
        let (mut d1, mut d2) = (0.0, 0.0);
        for (m, amp, ph, k, main) in E::DIST.iter() {
            let (s, c) = (dot13(m, if *main { &a } else { &b }) + ph).sin_cos();
            let rate = dot13(m, if *main { &rt } else { &bt });
            let shift = (m[0] + m[1] + m[2] + m[12]) as f64;
            let w = amp * tpow(t, *k);
            d1 += w * c * rate;
            d2 -= w * s * rate * shift;
        }
        let step = d1 / d2;
        dl -= step;
        if step.abs() < 1e-12 {
            break;
        }
    }
    let (a, b) = at(dl);
    let sum = |tab: &[([i8; 13], f64, f64, u8, bool)]| -> f64 {
        tab.iter()
            .map(|(m, amp, ph, k, main)| {
                amp * tpow(t, *k) * (dot13(m, if *main { &a } else { &b }) + ph).sin()
            })
            .sum()
    };
    let lon = sum(&E::LON);
    let lat = sum(&E::LAT);
    let r = sum(&E::DIST);
    // V = W1 + Δ + ΣL ; relative to W2 (+π)
    let rel = norm_pi(poly(&E::W1, t) - poly(&E::W2, t) + dl + lon - target);
    (rel, lat, r)
}

#[inline]
fn norm_pi(a: f64) -> f64 {
    let x = norm_rad(a);
    if x > PI {
        x - TAU
    } else {
        x
    }
}

/// Natural apogee (`apogee = true`) or perigee: (lon deg, lat deg, dist km),
/// mean ecliptic & equinox of date.
pub fn natural_apsis(jd_tt: f64, apogee: bool) -> (f64, f64, f64) {
    let t = (crate::time::tt_to_tdb(jd_tt) - J2000) / DAYS_PER_CENTURY;
    let (rel, lat, r) = theory(t, apogee);
    // mean apsis line of date (published ELP/MPP02 W2 + p_A, as `points`)
    let mean = super::mean_perigee_orbital_deg(jd_tt) + if apogee { 180.0 } else { 0.0 };
    (norm_deg(mean + rel / DEG), lat / DEG, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn theory_amplitudes() {
        // apogee within ~±6° of the mean, perigee within ~±30° (Koch: 5° / 25°)
        let (mut ma, mut mp) = (0.0f64, 0.0f64);
        for i in 0..2000 {
            let t = -0.3 + i as f64 * 0.0003;
            ma = ma.max(theory(t, true).0.abs());
            mp = mp.max(theory(t, false).0.abs());
        }
        assert!(ma > 4.0 * DEG && ma < 7.0 * DEG, "{}", ma / DEG);
        assert!(mp > 20.0 * DEG && mp < 32.0 * DEG, "{}", mp / DEG);
    }
}
