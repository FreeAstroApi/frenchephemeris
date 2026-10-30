// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Historical precession models, used only by the `SIDBIT_PREC_ORIG` option
//! ("compute the ayanamsa with the precession model its author used").
//!
//! - IAU 1976: Lieske, Lederle, Fricke & Morando (1977), A&A 58, 1 — the
//!   equatorial angles ζ_A, z_A, θ_A for an arbitrary start epoch, and the
//!   IAU 1976 mean obliquity ε = 84381.448″ − 46.8150″T − 0.00059″T² + 0.001813″T³.
//! - Newcomb: Kinoshita (1975, "Formulas for precession", SAO Special Report
//!   364) formulation of the Bessel–Newcomb precession — the formulation the
//!   Swiss documentation names — with the coefficients as published in
//!   SLALIB/PAL `palPrebn` (T from B1850.0, t in tropical centuries);
//!   obliquity = 23°27′08.26″ − 46.845″T − 0.0059″T² + 0.00181″T³ (T Julian
//!   centuries from 1900 Jan 0.5).

use crate::math::{cross, mat_mul, mat_t_vec, mat_vec, rot_x, rot_y, rot_z, Mat3, ARCSEC};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrigModel {
    Iau1976,
    Newcomb,
}

fn angles(m: OrigModel, jd0: f64, jd: f64) -> (f64, f64, f64) {
    match m {
        OrigModel::Iau1976 => {
            let tt = (jd0 - 2451545.0) / 36525.0;
            let t = (jd - jd0) / 36525.0;
            let w = 2306.2181 + 1.39656 * tt - 0.000139 * tt * tt;
            let zeta = (w + ((0.30188 - 0.000344 * tt) + 0.017998 * t) * t) * t;
            let z = (w + ((1.09468 + 0.000066 * tt) + 0.018203 * t) * t) * t;
            let th = 2004.3109 - 0.85330 * tt - 0.000217 * tt * tt;
            let theta = (th + (-(0.42665 + 0.000217 * tt) - 0.041833 * t) * t) * t;
            (zeta * ARCSEC, z * ARCSEC, theta * ARCSEC)
        }
        OrigModel::Newcomb => {
            // Kinoshita (1975) formulation of the Bessel–Newcomb precession,
            // as published in SLALIB/PAL `prebn` (Wallace): epochs Besselian,
            // T from B1850.0, t in tropical centuries.
            let bep = |jd: f64| 1900.0 + (jd - 2415020.31352) / 365.242198781;
            let bigt = (bep(jd0) - 1850.0) / 100.0;
            let t = (bep(jd) - bep(jd0)) / 100.0;
            let w = 2303.5548 + (1.39720 + 0.000059 * bigt) * bigt;
            let zeta = (w + (0.30242 - 0.000269 * bigt + 0.017996 * t) * t) * t;
            let z = (w + (1.09478 + 0.000387 * bigt + 0.018324 * t) * t) * t;
            let theta = (2005.1125
                + (-0.85294 - 0.000365 * bigt) * bigt
                + (-0.42647 - 0.000365 * bigt - 0.041802 * t) * t)
                * t;
            (zeta * ARCSEC, z * ARCSEC, theta * ARCSEC)
        }
    }
}

fn obliquity(m: OrigModel, jd: f64) -> f64 {
    match m {
        OrigModel::Iau1976 => {
            let t = (jd - 2451545.0) / 36525.0;
            (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t) * ARCSEC
        }
        OrigModel::Newcomb => {
            let t = (jd - 2415020.0) / 36525.0;
            (84428.26 + (-46.845 + (-0.0059 + 0.00181 * t) * t) * t) * ARCSEC
        }
    }
}

/// Precession matrix mean equator of `jd0` -> mean equator of `jd`.
fn prec(m: OrigModel, jd0: f64, jd: f64) -> Mat3 {
    let (zeta, z, theta) = angles(m, jd0, jd);
    mat_mul(&rot_z(-z), &mat_mul(&rot_y(theta), &rot_z(-zeta)))
}

/// Longitude (rad) of the mean equinox of `jd`, on the mean ecliptic of `jd0`,
/// counted from the equinox of `jd0`, in the given model.
pub fn equinox_lon_on_ecliptic_of(m: OrigModel, jd: f64, jd0: f64) -> f64 {
    let p = prec(m, jd0, jd);
    let eps = obliquity(m, jd);
    let eq_pole = mat_t_vec(&p, [0.0, 0.0, 1.0]);
    let ecl_pole = mat_t_vec(&p, [0.0, -eps.sin(), eps.cos()]);
    let e = mat_vec(&rot_x(obliquity(m, jd0)), cross(eq_pole, ecl_pole));
    e[1].atan2(e[0])
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn general_precession_rate() {
        // p_A over one century: IAU 1976 5029.0966" + 1.11113" (t² term); Newcomb ~5026"
        let a = -equinox_lon_on_ecliptic_of(OrigModel::Iau1976, 2451545.0 + 36525.0, 2451545.0) / ARCSEC;
        assert!((a - 5030.2077).abs() < 0.01, "{a}");
        let b = -equinox_lon_on_ecliptic_of(OrigModel::Newcomb, 2415020.0 + 36525.0, 2415020.0) / ARCSEC;
        assert!((b - 5026.75).abs() < 1.5, "{b}");
    }
}

/// Precession-model reconversion offset (arcsec) of a definition stated at
/// `t0` in model `m`: the difference Vondrák − `m` of the precession from t0
/// to J2000.0, i.e. the amount by which A0 must be changed so that sidereal
/// positions at J2000 are the same in both models. This reproduces the
/// offsets published in the Swiss documentation (IAU 1976 at 1956.2:
/// 0.13036″; Newcomb at B1950: 0.4135″ vs 0.41256″ published, at J1900:
/// 0.8303″ vs 0.828″) and gives the undocumented ones (Lahiri ICRC, 1956).
pub fn reconversion_offset_arcsec(m: OrigModel, t0: f64) -> f64 {
    let n = equinox_lon_on_ecliptic_of(m, 2451545.0, t0);
    let v = super::vondrak::equinox_lon_on_ecliptic_of(2451545.0, t0);
    (n - v) / ARCSEC
}

#[cfg(test)]
mod offset_tests {
    use super::*;
    #[test]
    fn documented_offsets() {
        let l = reconversion_offset_arcsec(OrigModel::Iau1976, 2435553.5);
        assert!((l + 0.13036).abs() < 2e-5, "{l}");
        let fb = reconversion_offset_arcsec(OrigModel::Newcomb, 2433282.42345905);
        assert!((fb - 0.41256).abs() < 0.0015, "{fb}");
        let j1900 = reconversion_offset_arcsec(OrigModel::Newcomb, 2415020.0);
        assert!((j1900 - 0.828).abs() < 0.0025, "{j1900}");
    }
}
