// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Long-term precession of the ecliptic and equator, Vondrák, Capitaine &
//! Wallace (2011), "New precession expressions, valid for long time
//! intervals", A&A 534, A22 (with the 2012 erratum). Valid ±200 000 years.
//!
//! Vectors are expressed in the J2000.0 mean equator & equinox frame
//! (no frame bias). Implemented here, locally, so that the ayanamsa module
//! has no dependency on `frames` (which may choose a different model).

use crate::math::{cross, dot, norm, scale, Vec3, ARCSEC, TAU};

/// Obliquity at J2000.0 used by Vondrák et al. (IAU 2006 value), radians.
pub const EPS0: f64 = 84381.406 * ARCSEC;

// Ecliptic pole: polynomial (arcsec) for P_A, Q_A.
const PQPOL: [[f64; 4]; 2] = [
    [5851.607687, -0.1189000, -0.00028913, 0.000000101],
    [-1600.886300, 1.1689818, -0.00000020, -0.000000437],
];
// Ecliptic pole periodic terms: period (centuries), C_P, C_Q, S_P, S_Q (arcsec).
const PQPER: [[f64; 5]; 8] = [
    [708.15, -5486.751211, -684.661560, 667.666730, -5523.863691],
    [2309.00, -17.127623, 2446.283880, -2354.886252, -549.747450],
    [1620.00, -617.517403, 399.671049, -428.152441, -310.998056],
    [492.20, 413.442940, -356.652376, 376.202861, 421.535876],
    [1183.00, 78.614193, -186.387003, 184.778874, -36.776172],
    [622.00, -180.732815, -316.800070, 335.321713, -145.278396],
    [882.00, -87.676083, 198.296701, -185.138669, -34.744450],
    [547.00, 46.140315, 101.135679, -120.972830, 22.885731],
];
// Equator pole: polynomial (arcsec) for X_A, Y_A.
const XYPOL: [[f64; 4]; 2] = [
    [5453.282155, 0.4252841, -0.00037173, -0.000000152],
    [-73750.930350, -0.7675452, -0.00018725, 0.000000231],
];
// Equator pole periodic terms: period (centuries), C_X, C_Y, S_X, S_Y (arcsec).
const XYPER: [[f64; 5]; 14] = [
    [256.75, -819.940624, 75004.344875, 81491.287984, 1558.515853],
    [708.15, -8444.676815, 624.033993, 787.163481, 7774.939698],
    [274.20, 2600.009459, 1251.136893, 1251.296102, -2219.534038],
    [241.45, 2755.175630, -1102.212834, -1257.950837, -2523.969396],
    [2309.00, -167.659835, -2660.664980, -2966.799730, 247.850422],
    [492.20, 871.855056, 699.291817, 639.744522, -846.485643],
    [396.10, 44.769698, 153.167220, 131.600209, -1393.124055],
    [288.90, -512.313065, -950.865637, -445.040117, 368.526116],
    [231.10, -819.415595, 499.754645, 584.522874, 749.045012],
    [1610.00, -538.071099, -145.188210, -89.756563, 444.704518],
    [620.00, -189.793622, 558.116553, 524.429630, 235.934465],
    [157.87, -402.922932, -23.923029, -13.549067, 374.049623],
    [220.30, 179.516345, -165.405086, -210.157124, -171.330180],
    [1200.00, -9.814756, 9.344131, -44.919798, -22.899655],
];

/// Julian centuries of 36525 days (TT) since J2000.0.
#[inline]
fn cent(jd_tt: f64) -> f64 {
    (jd_tt - crate::math::J2000) / crate::math::DAYS_PER_CENTURY
}

fn series(t: f64, pol: &[[f64; 4]; 2], per: &[[f64; 5]]) -> (f64, f64) {
    let (mut a, mut b) = (0.0, 0.0);
    for p in per {
        let w = TAU * t / p[0];
        let (s, c) = w.sin_cos();
        a += c * p[1] + s * p[3];
        b += c * p[2] + s * p[4];
    }
    let mut w = 1.0;
    for i in 0..4 {
        a += pol[0][i] * w;
        b += pol[1][i] * w;
        w *= t;
    }
    (a * ARCSEC, b * ARCSEC)
}

/// Unit vector of the ecliptic pole of date (J2000 mean equatorial frame).
pub fn ecliptic_pole(jd_tt: f64) -> Vec3 {
    let (p, q) = series(cent(jd_tt), &PQPOL, &PQPER);
    let w = (1.0 - p * p - q * q).max(0.0).sqrt();
    let (s, c) = EPS0.sin_cos();
    [p, -q * c - w * s, -q * s + w * c]
}

/// Unit vector of the mean equator pole of date (J2000 mean equatorial frame).
pub fn equator_pole(jd_tt: f64) -> Vec3 {
    let (x, y) = series(cent(jd_tt), &XYPOL, &XYPER);
    let w = (1.0 - x * x - y * y).max(0.0).sqrt();
    [x, y, w]
}

/// Rows of the rotation J2000 mean equatorial -> mean ecliptic & equinox of date.
pub fn ecliptic_of_date_matrix(jd_tt: f64) -> crate::math::Mat3 {
    let z = ecliptic_pole(jd_tt);
    let e = cross(equator_pole(jd_tt), z);
    let x = scale(e, 1.0 / norm(e));
    let y = cross(z, x);
    [x, y, z]
}

/// Ecliptic longitude (radians, (-π, π]) of the mean equinox of `jd_tt`
/// measured on the mean ecliptic of `jd0_tt` from the equinox of `jd0_tt`.
pub fn equinox_lon_on_ecliptic_of(jd_tt: f64, jd0_tt: f64) -> f64 {
    let m0 = ecliptic_of_date_matrix(jd0_tt);
    let z = ecliptic_pole(jd_tt);
    let e = cross(equator_pole(jd_tt), z);
    dot(e, m0[1]).atan2(dot(e, m0[0]))
}
