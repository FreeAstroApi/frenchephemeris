// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Lunar nodes and Black Moon Lilith.
//!
//! CONTRACT (owned by the astrology-geometry workstream):
//! - Degrees out. Mean points are analytic in JD TT (ecliptic & mean equinox of date).
//! - Osculating ("true") points take the geocentric Moon state in the
//!   ecliptic-of-date frame (km, km/day) so this module has no ephemeris dependency.
//!
//! Mean elements (published theory only — nothing fitted to any oracle).
//! The mean node Ω and mean perigee ϖ are the ELP/MPP02 secular polynomials
//! W3, W2 of Chapront & Francou (2003, A&A 404, 735; solution of Chapront,
//! Chapront-Touzé & Francou 2002, A&A 387, 700), which refer to the inertial
//! J2000 ecliptic, converted to the mean equinox of date by adding the IAU 2006
//! general precession in longitude p_A (Capitaine, Wallace & Chapront 2003,
//! A&A 412, 567, eq. 39). The small LLR/DE405 fit corrections (ΔW) of
//! ELP/MPP02 are not applied.
//! The mean apogee lies in the mean lunar orbit plane, inclined
//! [`MEAN_INCLINATION_DEG`] to the ecliptic of date, at orbital longitude
//! ϖ + 180°; it is projected onto the ecliptic (so mean Lilith has a
//! latitude of up to ±5.1°).
//!
//! Osculating points. Computed from the *geometric* geocentric Moon state at
//! the light-retarded time `t - |r|/c` (position AND velocity taken at that
//! time), in an inertial frame, then rotated to the mean ecliptic of date
//! (precession only; add Δψ for apparent output). Pass `P·r`, `P·v` here
//! (P = J2000 -> mean ecliptic of date), NOT a velocity that includes the
//! rotation rate of the date frame. GM is a parameter (use the ephemeris'
//! GM_Earth + GM_Moon).

use crate::math::{cross, dot, norm, norm_deg, scale, sub, Vec3, DEG, J2000};

const AS: f64 = 1.0 / 3600.0;

/// Natural (interpolated) apogee and perigee, Swiss ids 21/22.
pub mod intp;
#[rustfmt::skip] // generated
#[allow(dead_code)]
mod elp_main;

/// ELP/MPP02 W3 (mean node, J2000 inertial ecliptic): constant and T¹..T⁴ (arcsec).
const NODE0_DEG: f64 = 125.0 + 2.0 / 60.0 + 40.39816 * AS;
const NODE_W3: [f64; 4] = [-6967919.5383, 6.3590, 0.007625, -0.00003586];
/// ELP/MPP02 W2 (mean perigee, J2000 inertial ecliptic): constant and T¹..T⁴ (arcsec).
const PERI0_DEG: f64 = 83.0 + 21.0 / 60.0 + 11.67475 * AS;
const PERI_W2: [f64; 4] = [14643420.3171, -38.2631, -0.045047, 0.00021301];
/// IAU 2006 general precession in longitude p_A, T¹..T⁵ (arcsec).
const PREC_PA: [f64; 5] = [5028.796195, 1.1054348, 0.00007964, -0.000023857, -0.0000000383];

/// Mean inclination of the lunar orbit to the ecliptic, 5°08′43.427″, as
/// printed in the Astronomical Almanac (Moon section, p. D2, editions
/// before 2005; the USNO 2005 errata changed the printed value to
/// 5.156 6898°, which is the ELP "γ" mean inclination — see the report).
pub const MEAN_INCLINATION_DEG: f64 = 5.1453964;

/// Conventional distances reported for the mean points (not physical
/// positions): the mean Earth–Moon distance 384 400 km (Astronomical
/// Almanac lunar constants) for the node, and 384 400 km × (1 + e) with the
/// mean eccentricity e = 0.054 900 489 (Astronomical Almanac) for the apogee.
pub const MEAN_MOON_DIST_KM: f64 = 384_400.0;
pub const MEAN_ECCENTRICITY: f64 = 0.054900489;
pub const MEAN_NODE_DIST_AU: f64 = MEAN_MOON_DIST_KM / crate::math::AU_KM;
pub const MEAN_APOGEE_DIST_AU: f64 = MEAN_MOON_DIST_KM * (1.0 + MEAN_ECCENTRICITY) / crate::math::AU_KM;

#[inline]
fn cent(jd_tt: f64) -> f64 {
    (jd_tt - J2000) / crate::math::DAYS_PER_CENTURY
}

#[inline]
fn poly_as(c: &[f64], t: f64) -> f64 {
    // Σ c_i T^(i+1), arcsec
    c.iter().rev().fold(0.0, |acc, &ci| (acc + ci) * t)
}

#[inline]
fn poly_rate_as(c: &[f64], t: f64) -> f64 {
    // d/dT Σ c_i T^(i+1)
    c.iter().enumerate().rev().fold(0.0, |acc, (i, &ci)| acc * t + (i as f64 + 1.0) * ci)
}

#[inline]
fn of_date_deg(c0_deg: f64, w: &[f64; 4], t: f64) -> f64 {
    let v = poly_as(w, t) + poly_as(&PREC_PA, t);
    norm_deg(c0_deg + (v % 1_296_000.0) * AS)
}

/// Mean ascending lunar node longitude, degrees (mean equinox of date).
pub fn mean_node_deg(jd_tt: f64) -> f64 {
    of_date_deg(NODE0_DEG, &NODE_W3, cent(jd_tt))
}

/// Mean node daily motion, degrees/day.
pub fn mean_node_speed_deg_per_day(jd_tt: f64) -> f64 {
    let t = cent(jd_tt);
    (poly_rate_as(&NODE_W3, t) + poly_rate_as(&PREC_PA, t)) * AS / crate::math::DAYS_PER_CENTURY
}

/// Mean lunar perigee, orbital longitude (node + argument), degrees.
pub fn mean_perigee_orbital_deg(jd_tt: f64) -> f64 {
    of_date_deg(PERI0_DEG, &PERI_W2, cent(jd_tt))
}

/// Mean lunar apogee (mean Black Moon Lilith): (lon, lat) degrees,
/// ecliptic & mean equinox of date.
pub fn mean_apogee_deg(jd_tt: f64) -> (f64, f64) {
    project_apogee(mean_node_deg(jd_tt), mean_perigee_orbital_deg(jd_tt))
}

fn project_apogee(node: f64, perigee_orbital: f64) -> (f64, f64) {
    let u = (perigee_orbital + 180.0 - node) * DEG;
    let (si, ci) = (MEAN_INCLINATION_DEG * DEG).sin_cos();
    let (su, cu) = u.sin_cos();
    let lon = norm_deg(node + (su * ci).atan2(cu) / DEG);
    let lat = (su * si).asin() / DEG;
    (lon, lat)
}

/// Mean apogee (lon, lat) daily motion, degrees/day (central difference).
pub fn mean_apogee_speed_deg_per_day(jd_tt: f64) -> (f64, f64) {
    let h = 0.01;
    let (l0, b0) = mean_apogee_deg(jd_tt - h);
    let (l1, b1) = mean_apogee_deg(jd_tt + h);
    (crate::math::diff_deg(l1, l0) / (2.0 * h), (b1 - b0) / (2.0 * h))
}

/// Osculating (true) ascending node longitude, degrees.
pub fn osculating_node_deg(moon_pos: Vec3, moon_vel: Vec3) -> f64 {
    let h = cross(moon_pos, moon_vel);
    norm_deg(h[0].atan2(-h[1]) / DEG)
}

/// Osculating node: (lon deg, distance of the node point on the osculating
/// ellipse in km — the value Swiss reports as the node "distance").
pub fn osculating_node_full(moon_pos: Vec3, moon_vel: Vec3, gm_km3_day2: f64) -> (f64, f64) {
    let h = cross(moon_pos, moon_vel);
    let lon = norm_deg(h[0].atan2(-h[1]) / DEG);
    let ev = ecc_vector(moon_pos, moon_vel, h, gm_km3_day2);
    let nvec = [-h[1], h[0], 0.0];
    let nn = norm(nvec);
    let cosnu = if nn > 0.0 { dot(ev, nvec) / nn } else { 0.0 };
    let p = dot(h, h) / gm_km3_day2;
    (lon, p / (1.0 + cosnu))
}

#[inline]
fn ecc_vector(r: Vec3, v: Vec3, h: Vec3, gm: f64) -> Vec3 {
    sub(scale(cross(v, h), 1.0 / gm), scale(r, 1.0 / norm(r)))
}

/// Osculating apogee (true Lilith): (lon, lat, distance km).
/// `gm_km3_day2` is GM(Earth+Moon) in km^3/day^2.
pub fn osculating_apogee_deg(moon_pos: Vec3, moon_vel: Vec3, gm_km3_day2: f64) -> (f64, f64, f64) {
    let h = cross(moon_pos, moon_vel);
    let ev = ecc_vector(moon_pos, moon_vel, h, gm_km3_day2);
    let e = norm(ev);
    let p = dot(h, h) / gm_km3_day2;
    let a = p / (1.0 - e * e);
    let apo = scale(ev, -a * (1.0 + e) / e);
    let (lon, lat, r) = crate::math::to_spherical(apo);
    (lon / DEG, lat / DEG, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mean_node_j2000() {
        assert!((mean_node_deg(J2000) - NODE0_DEG).abs() < 1e-12);
        // daily motion about -0.0529539 deg/day
        assert!((mean_node_speed_deg_per_day(J2000) + 0.0529539).abs() < 1e-6);
        // analytic speed vs finite difference
        let jd = J2000 + 123456.7;
        let fd = crate::math::diff_deg(mean_node_deg(jd + 0.5), mean_node_deg(jd - 0.5));
        assert!((fd - mean_node_speed_deg_per_day(jd)).abs() < 1e-10);
    }
}
