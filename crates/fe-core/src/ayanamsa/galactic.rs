// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Galactic reference directions for the galactic-centre and
//! galactic-equator ayanamsas (internal definitions; if the star catalogue
//! later gains a Galactic Center record, these should agree with it).
//!
//! - Galactic centre = Sgr A*: the star catalogue's "Galactic Center" record
//!   (α = 17h45m40.03599s, δ = −29°00′28.1699″ ICRS, Petrov et al. 2011;
//!   proper motion Reid & Brunthaler 2020; distance GRAVITY Collaboration 2019).
//! - North galactic pole, IAU 1958 system (Blaauw et al. 1960, defined as
//!   α = 12h49m, δ = +27.4° B1950), in the ICRS as derived rigorously by
//!   Liu, Zhu & Zhang (2011), A&A 526, A16, eq. (26):
//!   α = 12h51m26.27469s, δ = +27°07′41.7087″.
//! - North galactic pole of the "new" (Sgr A*-based) galactic system proposed
//!   by the same paper, eq. (29)/(31): α = 12h51m36.7151981s,
//!   δ = +27°06′11.193172″ — the "modern" pole the Swiss documentation cites
//!   for the true galactic-equator ayanamsas (Liu/Zhu/Zhang 2010, "galactic
//!   node corrected by 3′11″").

use super::star::bias;
use super::vondrak::{ecliptic_of_date_matrix, equator_pole};
use crate::math::{cross, dot, mat_vec, norm, scale, Vec3, DEG};

fn radec(h: f64, m: f64, s: f64, d: f64, dm: f64, ds: f64) -> (f64, f64) {
    let a = (h + m / 60.0 + s / 3600.0) * 15.0 * DEG;
    let sg = if d < 0.0 { -1.0 } else { 1.0 };
    (a, sg * (d.abs() + dm / 60.0 + ds / 3600.0) * DEG)
}

fn unit(a: f64, d: f64) -> Vec3 {
    [d.cos() * a.cos(), d.cos() * a.sin(), d.sin()]
}

/// Barycentric ICRS position of Sgr A* (AU) at JD TDB, from the star
/// catalogue's "Galactic Center" record (single source of truth; see
/// `stars::catalog`: SIMBAD place Petrov et al. 2011, AJ 142, 35; proper
/// motion Reid & Brunthaler 2020, ApJ 892, 39; distance GRAVITY 2019).
pub fn sgr_a_bary_au(jd_tdb: f64) -> Vec3 {
    let rec = crate::stars::lookup("Galactic Center").expect("Galactic Center in the star catalogue").record;
    crate::stars::barycentric_au(rec, jd_tdb)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GalacticPole {
    /// IAU 1958 pole, ICRS (Liu, Zhu & Zhang 2011, eq. 26).
    Iau1958,
    /// Sgr A*-based pole proposed by Liu, Zhu & Zhang (2011, eq. 31).
    Modern,
}

pub fn pole_icrs(p: GalacticPole) -> Vec3 {
    let (a, d) = match p {
        GalacticPole::Iau1958 => radec(12.0, 51.0, 26.27469, 27.0, 7.0, 41.7087),
        GalacticPole::Modern => radec(12.0, 51.0, 36.7151981, 27.0, 6.0, 11.193172),
    };
    unit(a, d)
}

/// Longitude (deg) on the mean ecliptic of date of the galactic node near
/// 0° Capricorn (the descending node of the galactic equator), from the
/// "true" (unaberrated, catalogue) direction of the galactic pole.
pub fn galactic_node_lon_deg(p: GalacticPole, jd_tt: f64) -> f64 {
    let g = mat_vec(&ecliptic_of_date_matrix(jd_tt), bias(pole_icrs(p)));
    let n = cross([0.0, 0.0, 1.0], g);
    crate::math::norm_deg(n[1].atan2(n[0]) / DEG)
}

/// Polar ("dhruva") projection onto the ecliptic of a direction given in the
/// mean ecliptic of date: the ecliptic point on the great circle through the
/// celestial pole and the object (same right ascension). `pole_ecl` is the
/// mean celestial pole in the same frame.
pub fn polar_longitude_deg(obj_ecl: Vec3, pole_ecl: Vec3) -> f64 {
    let n = cross(pole_ecl, obj_ecl);
    let mut d = cross(n, [0.0, 0.0, 1.0]);
    if dot(d, obj_ecl) < 0.0 {
        d = scale(d, -1.0);
    }
    let d = scale(d, 1.0 / norm(d));
    crate::math::norm_deg(d[1].atan2(d[0]) / DEG)
}

/// Mean celestial pole of date in the mean ecliptic of date (Vondrák 2011).
pub fn mean_pole_in_ecliptic(jd_tt: f64) -> Vec3 {
    mat_vec(&ecliptic_of_date_matrix(jd_tt), equator_pole(jd_tt))
}
