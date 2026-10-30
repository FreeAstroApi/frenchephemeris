// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! IAU 1976/1980 frames and the JPL Horizons convention.
//!
//! JPL Horizons refers "apparent of date" positions to the IERS Conventions 1996
//! frame (IERS Conventions 1996, p. 22): IAU 1976 precession (Lieske et al. 1977,
//! A&A 58, 1; Lieske 1979, A&A 73, 282), IAU 1980 nutation (Wahr 1981; Seidelmann
//! 1982) plus the observed celestial pole offsets δψ, δε published by the IERS
//! relative to that model, with no ICRS frame bias (the offsets absorb the pole
//! part of it). Outside the span of the observed offsets Horizons holds the first
//! (before 1962) or last value (documented in the Swiss Ephemeris manual, §2,
//! from J. Giorgini). Two variants:
//!
//! * `PoleOffsets::EopC04` — the pinned IERS EOP 14 C04 (IAU 1980) daily series
//!   1962-01-01 … 2026-01-05 (`eop_c04_iau1980.bin`, see tools/gen_iau1980.py;
//!   source sha256 40347f18a10eb72b15e4f84e647b5f1b7da4709534d935e2e1289697bdfa79df).
//! * `PoleOffsets::Model` — no data file: δψ, δε are the offsets that move the
//!   IAU 1976/1980 pole onto the IAU 2006/2000A celestial intermediate pole, which
//!   is what the observed offsets measure apart from the free core nutation
//!   (≲ 0.3 mas). Held constant outside the same 1962–2026 span as Horizons.
//!
//! Formulas follow the cited publications; the series coefficients are the
//! published IAU 1980 table (nut80_table.rs). Validated against ERFA
//! (`nut80`, `pmat76`, `obl80`, `gmst82`, `eqeq94`).

use super::{era, FrameModel, FrameSet};
use crate::math::*;

mod nut80_table;

/// Julian centuries TT from J2000.
fn cent(jd_tt: f64) -> f64 {
    (jd_tt - J2000) / DAYS_PER_CENTURY
}

fn anpm(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r > PI {
        r - TAU
    } else {
        r
    }
}

/// IAU 1980 nutation (Δψ, Δε), radians.
pub fn nut80(jd_tt: f64) -> (f64, f64) {
    let t = cent(jd_tt);
    let arg = |c: [f64; 4], rev: f64| {
        anpm((c[0] + (c[1] + (c[2] + c[3] * t) * t) * t) * ARCSEC + (rev * t).rem_euclid(1.0) * TAU)
    };
    let el = arg([485866.733, 715922.633, 31.310, 0.064], 1325.0);
    let elp = arg([1287099.804, 1292581.224, -0.577, -0.012], 99.0);
    let f = arg([335778.877, 295263.137, -13.257, 0.011], 1342.0);
    let d = arg([1072261.307, 1105601.328, -6.891, 0.019], 1236.0);
    let om = arg([450160.280, -482890.539, 7.455, 0.008], -5.0);
    let (mut dp, mut de) = (0.0, 0.0);
    for (m, c) in nut80_table::NUT80.iter().rev() {
        let a = m[0] as f64 * el + m[1] as f64 * elp + m[2] as f64 * f + m[3] as f64 * d + m[4] as f64 * om;
        dp += (c[0] + c[1] * t) * a.sin();
        de += (c[2] + c[3] * t) * a.cos();
    }
    let u = ARCSEC / 1e4;
    (dp * u, de * u)
}

/// IAU 1980 mean obliquity, radians.
pub fn obl80(jd_tt: f64) -> f64 {
    let t = cent(jd_tt);
    (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t) * ARCSEC
}

/// IAU 1976 precession matrix, mean J2000 → mean of date (Lieske 1979, eqs. 6–7).
pub fn pmat76(jd_tt: f64) -> Mat3 {
    let t = cent(jd_tt);
    let ta = t * ARCSEC;
    let w = 2306.2181;
    let zeta = (w + (0.30188 + 0.017998 * t) * t) * ta;
    let z = (w + (1.09468 + 0.018203 * t) * t) * ta;
    let theta = (2004.3109 + (-0.42665 - 0.041833 * t) * t) * ta;
    mat_mul(&rot_z(-z), &mat_mul(&rot_y(theta), &rot_z(-zeta)))
}

/// Greenwich mean sidereal time, IAU 1982 (Aoki et al. 1982), radians, JD UT1.
pub fn gmst82(jd_ut: f64) -> f64 {
    let t = (jd_ut - J2000) / DAYS_PER_CENTURY;
    let f = 86400.0 * jd_ut.rem_euclid(1.0);
    let s = 24110.54841 - 43200.0 + (8640184.812866 + (0.093104 - 6.2e-6 * t) * t) * t + f;
    norm_rad(s * TAU / 86400.0)
}

/// Equation of the equinoxes, IAU 1994 (Capitaine & Gontier 1993), radians.
pub fn eqeq94(jd_tt: f64, dpsi: f64, eps0: f64) -> f64 {
    let t = cent(jd_tt);
    let om = anpm(
        (450160.280 + (-482890.539 + (7.455 + 0.008 * t) * t) * t) * ARCSEC
            + (-5.0 * t).rem_euclid(1.0) * TAU,
    );
    dpsi * eps0.cos() + ARCSEC * (0.00264 * om.sin() + 0.000063 * (2.0 * om).sin())
}

/// Source of the celestial pole offsets added to IAU 1980 nutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PoleOffsets {
    /// None: plain IAU 1976/1980.
    None,
    /// IERS EOP 14 C04 (IAU 1980) observed daily values (`SEFLG_JPLHOR`).
    EopC04,
    /// IAU 2006/2000A pole expressed as offsets (`SEFLG_JPLHOR_APPROX`).
    Model,
}

static EOP: &[u8] = include_bytes!("eop_c04_iau1980.bin");
/// MJD (0h UTC) of the first daily value.
pub const EOP_FIRST_MJD: f64 = 37665.0;

fn eop_len() -> usize {
    EOP.len() / 8
}

fn eop_at(i: usize) -> (f64, f64) {
    let r = |o: usize| i32::from_le_bytes(EOP[o..o + 4].try_into().unwrap()) as f64 * 1e-6 * ARCSEC;
    (r(8 * i), r(8 * i + 4))
}

/// Span covered by the observed offsets, MJD.
pub fn eop_span_mjd() -> (f64, f64) {
    (EOP_FIRST_MJD, EOP_FIRST_MJD + (eop_len() - 1) as f64)
}

/// Observed (δψ, δε), radians, linearly interpolated in time (the table is
/// UTC-dated; TT is used, a ≤ 70 s shift), held at the end values outside.
pub fn eop_offsets(jd_tt: f64) -> (f64, f64) {
    let x = (jd_tt - 2_400_000.5 - EOP_FIRST_MJD).clamp(0.0, (eop_len() - 1) as f64);
    let i = (x.floor() as usize).min(eop_len() - 2);
    let w = x - i as f64;
    let (a, b) = (eop_at(i), eop_at(i + 1));
    (a.0 + w * (b.0 - a.0), a.1 + w * (b.1 - a.1))
}

/// (δψ, δε) that rotate the IAU 1976/1980 pole onto the IAU 2006/2000A CIP at
/// `jd_tt` (held constant outside the observed span, like Horizons).
pub fn model_offsets(jd_tt: f64) -> (f64, f64) {
    let (a, b) = eop_span_mjd();
    let jd = (jd_tt - 2_400_000.5).clamp(a, b) + 2_400_000.5;
    let target = FrameSet::with_model(jd, FrameModel::IAU2006_2000A).bpn[2];
    let p = pmat76(jd);
    let eps = obl80(jd);
    let (dp, de) = nut80(jd);
    let pole = |x: f64, y: f64| mat_mul(&super::nutation_matrix(eps, dp + x, de + y), &p)[2];
    // The pole in the mean-of-date frame moves by (δψ·sin ε, −δε) for small
    // offsets: two Gauss–Newton steps on the exact matrices.
    let (mut x, mut y) = (0.0, 0.0);
    for _ in 0..3 {
        let h = 1e-9;
        let c = pole(x, y);
        let gx = scale(sub(pole(x + h, y), c), 1.0 / h);
        let gy = scale(sub(pole(x, y + h), c), 1.0 / h);
        let r = sub(target, c);
        // Least squares in 3-D for the 2 unknowns.
        let (a11, a12, a22) = (dot(gx, gx), dot(gx, gy), dot(gy, gy));
        let (b1, b2) = (dot(gx, r), dot(gy, r));
        let det = a11 * a22 - a12 * a12;
        x += (b1 * a22 - b2 * a12) / det;
        y += (a11 * b2 - a12 * b1) / det;
    }
    (x, y)
}

/// Frames in the Horizons convention at `jd_tt`; `jd_ut` sets the Earth
/// rotation (GAST = GMST 1982 + equation of the equinoxes 1994, which only
/// matters for topocentric positions). `model` is a placeholder for code that
/// reads it; the matrices are not those of any `FrameModel`.
pub fn horizons_frames(jd_tt: f64, jd_ut: f64, offsets: PoleOffsets) -> FrameSet {
    let bp = pmat76(jd_tt);
    let eps_mean = obl80(jd_tt);
    let (dp80, de80) = nut80(jd_tt);
    let (ddp, dde) = match offsets {
        PoleOffsets::None => (0.0, 0.0),
        PoleOffsets::EopC04 => eop_offsets(jd_tt),
        PoleOffsets::Model => model_offsets(jd_tt),
    };
    let (dpsi, deps) = (dp80 + ddp, de80 + dde);
    let bpn = mat_mul(&super::nutation_matrix(eps_mean, dpsi, deps), &bp);
    let eps_true = eps_mean + deps;
    let eqeq = eqeq94(jd_tt, dpsi, eps_mean);
    let gast = gmst82(jd_ut) + eqeq;
    FrameSet {
        jd_tt,
        eps_mean,
        dpsi,
        deps,
        eps_true,
        bp,
        bpn,
        model: FrameModel::IAU2006_2000A,
        ecl_mean: mat_mul(&rot_x(eps_mean), &bp),
        ecl_true: mat_mul(&rot_x(eps_true), &bpn),
        eo: norm_rad(era(jd_ut) - gast),
        eqeq,
    }
}
