// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Self-contained apparent place of δ Cancri (Asellus Australis, HIP 42911),
//! the yogatara of Pushya, for the True-Pushya ayanamsa.
//!
//! Catalogue: Hipparcos, the New Reduction (van Leeuwen 2007, A&A 474, 653;
//! VizieR I/311, HIP 42911), ICRS at epoch J1991.25:
//! α = 2.2893709279 rad, δ = 0.3168621473 rad, μα* = −17.67 mas/yr,
//! μδ = −229.26 mas/yr, ϖ = 24.98 mas (CDS catalogue I/311); radial velocity
//! +16.39 km/s (Massarotti et al. 2008, AJ 135, 209). Linear space motion.
//!
//! Earth's heliocentric state comes from the low-precision solar theory of
//! Meeus, *Astronomical Algorithms* (2nd ed.), ch. 25 — adequate here
//! because it only feeds annual parallax (25 mas × error), aberration
//! (20.5" × relative error) and light deflection. Then: frame bias
//! (IAU 2006), Vondrák 2011 precession to the mean ecliptic of date.
//! Measured accuracy of the resulting longitude vs the oracle: see
//! `docs/GEOMETRY_REPORT.md` (≲ 1e-5° over years 0–2400).

use super::vondrak::ecliptic_of_date_matrix;
use crate::math::{
    add, dot, mat_t_vec, mat_vec, norm, scale, sub, Vec3, ARCSEC, AU_KM, C_KM_DAY, DEG, J2000,
};

const MAS: f64 = ARCSEC / 1000.0;
const EPOCH_JD: f64 = 2448349.0625; // J1991.25
const RA0_DEG: f64 = 2.2893709279 / DEG;
const DEC0_DEG: f64 = 0.3168621473 / DEG;
const PM_RA_MAS: f64 = -17.67;
const PM_DEC_MAS: f64 = -229.26;
const PLX_MAS: f64 = 24.98;
const RV_KMS: f64 = 16.39;

/// Barycentric position of the star (AU, ICRS) at `jd_tt` (linear space motion).
fn star_bary(jd_tt: f64) -> Vec3 {
    let (sa, ca) = (RA0_DEG * DEG).sin_cos();
    let (sd, cd) = (DEC0_DEG * DEG).sin_cos();
    let r = 1.0 / (PLX_MAS * MAS);
    let u = [cd * ca, cd * sa, sd];
    let ea = [-sa, ca, 0.0];
    let ed = [-sd * ca, -sd * sa, cd];
    let v = add(
        scale(add(scale(ea, PM_RA_MAS * MAS), scale(ed, PM_DEC_MAS * MAS)), r / 365.25),
        scale(u, RV_KMS * 86400.0 / AU_KM),
    );
    add(scale(u, r), scale(v, jd_tt - EPOCH_JD))
}

/// IAU 2006 frame bias, ICRS -> J2000 mean equator & equinox (first order).
pub(crate) fn bias(v: Vec3) -> Vec3 {
    let dx = -0.016617 * ARCSEC;
    let de = -0.0068192 * ARCSEC;
    let dr = -0.0146 * ARCSEC;
    [v[0] + dr * v[1] - dx * v[2], -dr * v[0] + v[1] - de * v[2], dx * v[0] + de * v[1] + v[2]]
}

/// Heliocentric Earth, AU, J2000 mean equatorial frame (Meeus ch. 25).
fn earth_helio(jd_tt: f64) -> Vec3 {
    let t = (jd_tt - J2000) / 36525.0;
    let l0 = 280.46646 + 36000.76983 * t + 0.0003032 * t * t;
    let m = (357.52911 + 35999.05029 * t - 0.0001537 * t * t) * DEG;
    let e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;
    let c = (1.914602 - 0.004817 * t - 0.000014 * t * t) * m.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m).sin()
        + 0.000289 * (3.0 * m).sin();
    let lon = (l0 + c) * DEG;
    let nu = m + c * DEG;
    let r = 1.000001018 * (1.0 - e * e) / (1.0 + e * nu.cos());
    let ecl = [-r * lon.cos(), -r * lon.sin(), 0.0];
    mat_t_vec(&ecliptic_of_date_matrix(jd_tt), ecl)
}

/// Core apparent-place computation; all vectors in the J2000 mean equatorial
/// frame (AU, AU/day): Earth barycentric position/velocity, Sun barycentric
/// position.
fn apparent_core(jd_tt: f64, e: Vec3, ve: Vec3, sun: Vec3) -> (f64, f64) {
    apparent_of(jd_tt, bias(star_bary(jd_tt)), e, ve, sun, true)
}

/// Apparent mean-of-date ecliptic (lon, lat) in degrees of a star at
/// barycentric position `s` (AU, J2000 mean equatorial, bias applied), seen
/// from Earth at `e` with velocity `ve`, Sun at `sun` (same frame).
fn apparent_of(jd_tt: f64, s: Vec3, e: Vec3, ve: Vec3, sun: Vec3, deflect: bool) -> (f64, f64) {
    let p = sub(s, e);
    let mut u = scale(p, 1.0 / norm(p));
    // gravitational light deflection by the Sun (no limb/disc treatment)
    let g_on = if deflect { 1.0 } else { 0.0 };
    let es = sub(e, sun);
    let ed = norm(es);
    let eh = scale(es, 1.0 / ed);
    let qs = sub(s, sun);
    let q = scale(qs, 1.0 / norm(qs));
    // 2 GM_sun / c^2 in AU
    let g = g_on * 2.0 * 1.327_124_400_41e20 / (299_792_458.0f64 * 299_792_458.0) / 1.495_978_707e11;
    let ue = dot(u, eh);
    u = add(u, scale(sub(eh, scale(u, ue)), g / ed / (1.0 + dot(q, eh))));
    u = scale(u, 1.0 / norm(u));
    // annual aberration (first order in v/c)
    let b = scale(ve, AU_KM / C_KM_DAY);
    u = add(u, sub(b, scale(u, dot(u, b))));
    u = scale(u, 1.0 / norm(u));
    let w = mat_vec(&ecliptic_of_date_matrix(jd_tt), u);
    let lon = crate::math::norm_deg(w[1].atan2(w[0]) / DEG);
    let lat = w[2].clamp(-1.0, 1.0).asin() / DEG;
    (lon, lat)
}

/// Apparent (annual parallax, light deflection by the Sun, annual
/// aberration) geocentric ecliptic longitude and latitude of δ Cnc referred
/// to the MEAN ecliptic and equinox of date (no nutation), degrees, using
/// the built-in low-precision Earth (Sun at the origin).
pub fn delta_cnc_apparent_mean(jd_tt: f64) -> (f64, f64) {
    let h = 0.01;
    let e = earth_helio(jd_tt);
    let ve = scale(sub(earth_helio(jd_tt + h), earth_helio(jd_tt - h)), 0.5 / h);
    apparent_core(jd_tt, e, ve, [0.0; 3])
}

/// Apparent mean-of-date ecliptic (lon, lat), degrees, of any star whose
/// barycentric ICRS position (AU) at `jd_tt` is `star_icrs_au`, using the
/// built-in low-precision Earth (same model as [`delta_cnc_apparent_mean`]).
pub fn apparent_mean_builtin(jd_tt: f64, star_icrs_au: Vec3) -> (f64, f64) {
    let h = 0.01;
    let e = earth_helio(jd_tt);
    let ve = scale(sub(earth_helio(jd_tt + h), earth_helio(jd_tt - h)), 0.5 / h);
    apparent_of(jd_tt, bias(star_icrs_au), e, ve, [0.0; 3], true)
}

/// Geometric (no aberration, no deflection) geocentric mean-of-date
/// ecliptic (lon, lat), degrees, of a star at barycentric ICRS `star_icrs_au`,
/// with the built-in low-precision Earth (annual parallax only).
pub fn geometric_mean_builtin(jd_tt: f64, star_icrs_au: Vec3) -> (f64, f64) {
    let p = sub(bias(star_icrs_au), earth_helio(jd_tt));
    let w = mat_vec(&ecliptic_of_date_matrix(jd_tt), p);
    (crate::math::norm_deg(w[1].atan2(w[0]) / DEG), (w[2] / norm(w)).clamp(-1.0, 1.0).asin() / DEG)
}

/// Same as [`delta_cnc_apparent_mean`] but from exact ephemeris states:
/// barycentric Earth position (AU) and velocity (AU/day) and barycentric Sun
/// position (AU), all ICRS. With the oracle's own states this reproduces the
/// oracle to ~1e-9°, except while the star is behind the solar disc
/// (separation < 0.27°), where Swiss attenuates the deflection (not modelled).
pub fn delta_cnc_apparent_mean_from_states(
    jd_tt: f64,
    earth_bary_au: Vec3,
    earth_vel_au_day: Vec3,
    sun_bary_au: Vec3,
) -> (f64, f64) {
    apparent_core(jd_tt, bias(earth_bary_au), bias(earth_vel_au_day), bias(sun_bary_au))
}

/// [`delta_cnc_apparent_mean_from_states`] with the solar light deflection and the
/// annual aberration each optional, to match a body computed without them.
pub fn delta_cnc_mean_from_states_with(
    jd_tt: f64,
    earth_bary_au: Vec3,
    earth_vel_au_day: Vec3,
    sun_bary_au: Vec3,
    deflection: bool,
    aberration: bool,
) -> (f64, f64) {
    let ve = if aberration { bias(earth_vel_au_day) } else { [0.0; 3] };
    apparent_of(jd_tt, bias(star_bary(jd_tt)), bias(earth_bary_au), ve, bias(sun_bary_au), deflection)
}
