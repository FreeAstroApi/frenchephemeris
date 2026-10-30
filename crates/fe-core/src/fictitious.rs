// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Hypothetical ("fictitious") bodies, Swiss ids `SE_FICT_OFFSET` (40) + n.
//!
//! Each body is an unperturbed Kepler orbit about the Sun (or the Earth)
//! given by published orbital elements; provenance of every element set is
//! in docs/FICTITIOUS_BODIES.md. Conventions (Swiss Ephemeris programmer's
//! documentation, "Fictitious planets"): elements are referred to the mean
//! ecliptic and equinox of a stated epoch (or of date); the mean motion
//! follows from the semi-major axis and the Sun's gravity (Gaussian
//! constant, massless body) unless the source states a different one.
//!
//! The elements give a heliocentric (geocentric) vector in the ecliptic of
//! their equinox; it is rotated to GCRS with the engine's precession model
//! and then goes through the normal apparent-place pipeline (light time,
//! deflection, aberration, frames) like a planet.
//!
//! Not provided (no published element set could be obtained — see the
//! report): 48 Isis-Transpluto, 49 Nibiru, 54 Pluto (Pickering), 55 Vulcan,
//! 57 Proserpina.

use crate::frames::{self, FrameModel};
use crate::math::*;
use crate::pack::Pack;
use crate::Result;

/// Swiss `SE_FICT_OFFSET`.
pub const FICT_OFFSET: i32 = 40;

/// Gaussian gravitational constant k (rad/day, AU, solar masses).
pub const GAUSS_K: f64 = 0.017_202_098_95;
/// GM of the Earth, km³/s² (IERS 2010 / DE440 value 398 600.435 436).
pub const GM_EARTH_KM3_S2: f64 = 398_600.435_436;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Equinox {
    /// Mean ecliptic and equinox of the given JD (TT).
    Fixed(f64),
    /// Mean ecliptic and equinox of date.
    OfDate,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Center {
    Sun,
    Earth,
}

/// Keplerian elements (angles in degrees, distances in AU).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements {
    pub id: i32,
    pub name: &'static str,
    pub center: Center,
    /// Epoch of the mean anomaly, JD TT.
    pub epoch: f64,
    pub equinox: Equinox,
    /// Mean anomaly at epoch.
    pub m0: f64,
    /// Mean motion, deg/day; `None` = from `a` and the central body's GM.
    pub n: Option<f64>,
    pub a: f64,
    pub e: f64,
    /// Argument of perihelion (from the node).
    pub w: f64,
    /// Longitude of the ascending node.
    pub node: f64,
    pub i: f64,
}

const J1900: f64 = 2_415_020.0;
const J2000_: f64 = 2_451_545.0;
const fn dms(d: f64, m: f64, s: f64) -> f64 {
    d + m / 60.0 + s / 3600.0
}

/// Uranian (Hamburg) planets: J. Neely, "Orbital Elements for the
/// Transneptunian Planets", Matrix Journal VII (1980), Table I, epoch and
/// equinox 1900.0 (T from JD 2415020). Neely's node term +1.3960T is the
/// general precession, which a fixed equinox omits; mean motion from `a`.
const fn neely(id: i32, name: &'static str, a: f64, e: f64, i: f64, w: f64, m0: f64, node: f64) -> Elements {
    Elements {
        id,
        name,
        center: Center::Sun,
        epoch: J1900,
        equinox: Equinox::Fixed(J1900),
        m0,
        n: None,
        a,
        e,
        w,
        node,
        i,
    }
}

pub const ELEMENTS: [Elements; 14] = [
    neely(40, "Cupido", 40.99837, 0.00460, 1.0833, 171.4333, 163.7409, 129.8325),
    neely(41, "Hades", 50.66744, 0.00245, 1.0500, 148.1796, 27.6496, 161.3339),
    neely(42, "Zeus", 59.21436, 0.00120, 0.0, 299.0440, 165.1232, 0.0),
    neely(43, "Kronos", 64.81690, 0.00305, 0.0, 208.8801, 169.0193, 0.0),
    neely(44, "Apollon", 70.29949, 0.0, 0.0, 138.0533, 0.0, 0.0),
    neely(45, "Admetos", 73.62765, 0.0, 0.0, 351.3350, 0.0, 0.0),
    neely(46, "Vulkanus", 77.25568, 0.0, 0.0, 55.8983, 0.0, 0.0),
    neely(47, "Poseidon", 83.66907, 0.0, 0.0, 165.5163, 0.0, 0.0),
    // R. S. Harrington, "The Location of Planet X", AJ 96 (1988) 1476, p. 1478
    // ("typical good case"): perihelion 1789 Aug 6 0h (JD 2374696.5). The
    // paper states no equinox; J2000 is the reading Swiss follows (B1950
    // would put the body 0.70° further east; docs/FICTITIOUS_BODIES.md).
    Elements {
        id: 50,
        name: "Harrington",
        center: Center::Sun,
        epoch: 2_374_696.5,
        equinox: Equinox::Fixed(J2000_),
        m0: 0.0,
        n: None,
        a: 101.2,
        e: 0.411,
        w: 208.5,
        node: 275.4,
        i: 32.4,
    },
    // U.-J. Le Verrier, Comptes Rendus 23 (1846) 428–438, séance 31 Aug 1846,
    // p. 432: a 36,154; e 0,10761; ϖ 284°45′; mean longitude 318°47′ on
    // 1847 Jan 1, equinox of 1847 Jan 1 (JD 2395662.5); no inclination.
    Elements {
        id: 51,
        name: "Neptune (Leverrier)",
        center: Center::Sun,
        epoch: 2_395_662.5,
        equinox: Equinox::Fixed(2_395_662.5),
        m0: dms(318.0, 47.0, 0.0) - dms(284.0, 45.0, 0.0),
        n: None,
        a: 36.154,
        e: 0.10761,
        w: dms(284.0, 45.0, 0.0),
        node: 0.0,
        i: 0.0,
    },
    // J. C. Adams, MNRAS 7 (1846) 149, p. 151, hypothesis II (sent to Airy
    // Sept 1846): mean distance = a(Uranus)/0.515 with a(Uranus) = 19.18239;
    // mean longitude 323°2′ on 1846 Oct 6 (JD 2395575.5); ϖ 299°11′;
    // e 0.120615; no inclination; equinox of the epoch.
    Elements {
        id: 52,
        name: "Neptune (Adams)",
        center: Center::Sun,
        epoch: 2_395_575.5,
        equinox: Equinox::Fixed(2_395_575.5),
        m0: dms(323.0, 2.0, 0.0) - dms(299.0, 11.0, 0.0),
        n: None,
        a: 19.18239 / 0.515,
        e: 0.120615,
        w: dms(299.0, 11.0, 0.0),
        node: 0.0,
        i: 0.0,
    },
    // P. Lowell, Memoir on a Trans-Neptunian Planet (1915), solution X1, as
    // quoted by E. C. Bower, Lick Obs. Bull. 15 (1931) No. 437, p. 175:
    // a 43.0, e 0.202, heliocentric longitude 84.0 on 1914 Jul 0
    // (JD 2420313.5), ϖ 203.8; "i 10°±" without a node is not usable, so the
    // orbit is taken in the ecliptic; equinox of the epoch. The longitude
    // is the true heliocentric longitude (M0 computed below).
    Elements {
        id: 53,
        name: "Pluto (Lowell)",
        center: Center::Sun,
        epoch: 2_420_313.5,
        equinox: Equinox::Fixed(2_420_313.5),
        m0: f64::NAN, // filled by `lowell_m0`
        n: None,
        a: 43.0,
        e: 0.202,
        w: 203.8,
        node: 0.0,
        i: 0.0,
    },
    // White Moon / Selena (P. Globa): geocentric circle, period exactly 7
    // (Julian) years in tropical longitude, 88°03′ on 1990 Jan 1 0h UT
    // (JD 2447892.5); distance from Kepler's law with GM(Earth).
    Elements {
        id: 56,
        name: "Selena/White Moon",
        center: Center::Earth,
        epoch: 2_447_892.5,
        equinox: Equinox::OfDate,
        m0: dms(88.0, 3.0, 0.0),
        n: Some(360.0 / (7.0 * 365.25)),
        a: f64::NAN, // Kepler, see `geo_a`
        e: 0.0,
        w: 0.0,
        node: 0.0,
        i: 0.0,
    },
    // Waltemath's moon, Monthly Weather Review 26 (1898) 19–20: synodic
    // period ~177 d, predicted transit of the Sun on 1898 Feb 3. Circular
    // geocentric orbit; sidereal period from 1/P = 1/177 + 1/365.242189;
    // at 1898 Feb 3 12h UT (JD 2414324.0) it is at the Sun's geocentric
    // longitude (DE441: 314.7368°, mean ecliptic of date); equinox of epoch.
    Elements {
        id: 58,
        name: "Waldemath",
        center: Center::Earth,
        epoch: 2_414_324.0,
        equinox: Equinox::Fixed(2_414_324.0),
        m0: 314.7368,
        n: Some(360.0 * (1.0 / 177.0 + 1.0 / 365.242189)),
        a: f64::NAN,
        e: 0.0,
        w: 0.0,
        node: 0.0,
        i: 0.0,
    },
];

/// Elements of a Swiss fictitious-body id, if provided.
pub fn elements(id: i32) -> Option<&'static Elements> {
    ELEMENTS.iter().find(|e| e.id == id)
}

/// Mean anomaly (deg) at epoch reproducing a true anomaly `nu` (deg).
fn mean_from_true(nu_deg: f64, e: f64) -> f64 {
    let nu = nu_deg * DEG;
    let ea = 2.0 * (((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan()).atan();
    (ea - e * ea.sin()) / DEG
}

fn m0(el: &Elements) -> f64 {
    if el.id == 53 {
        // Lowell X1: true heliocentric longitude 84.0, ϖ 203.8
        mean_from_true(84.0 - 203.8, el.e)
    } else {
        el.m0
    }
}

/// Semi-major axis in AU (geocentric circles: Kepler's third law with GM(Earth)).
fn semi_axis(el: &Elements) -> f64 {
    if el.a.is_nan() {
        let gm = GM_EARTH_KM3_S2 * 86400.0 * 86400.0; // km³/day²
        let n = el.n.expect("geocentric body has a stated mean motion") * DEG;
        (gm / (n * n)).cbrt() / AU_KM
    } else {
        el.a
    }
}

/// Mean motion, rad/day.
fn mean_motion(el: &Elements, a: f64) -> f64 {
    match el.n {
        Some(n) => n * DEG,
        None => {
            let gm = match el.center {
                Center::Sun => GAUSS_K * GAUSS_K,
                Center::Earth => GM_EARTH_KM3_S2 * 86400.0 * 86400.0 / AU_KM.powi(3),
            };
            (gm / (a * a * a)).sqrt()
        }
    }
}

/// GCRS → mean ecliptic & equinox of `jd_tt`.
fn ecliptic_matrix(jd_tt: f64, model: FrameModel) -> Mat3 {
    let bp = frames::bias_precession_matrix(jd_tt, model.precession);
    let eps = frames::obliquity_mean_with(jd_tt, model.precession, model.obliquity);
    mat_mul(&rot_x(eps), &bp)
}

/// Position and velocity relative to the central body (km, km/day), GCRS axes,
/// at JD `jd` (TT ≈ TDB for these orbits).
pub fn relative_state(el: &Elements, jd: f64, model: FrameModel) -> (Vec3, Vec3) {
    let a = semi_axis(el);
    let n = mean_motion(el, a);
    let e = el.e;
    let m = m0(el) * DEG + n * (jd - el.epoch);
    let m = norm_rad(m);
    let mut ea = if e > 0.8 { PI } else { m };
    for _ in 0..30 {
        let d = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
        ea -= d;
        if d.abs() < 1e-14 {
            break;
        }
    }
    let (se, ce) = ea.sin_cos();
    let q = (1.0 - e * e).sqrt();
    let akm = a * AU_KM;
    let pos = [akm * (ce - e), akm * q * se, 0.0];
    let edot = n / (1.0 - e * ce);
    let vel = [-akm * se * edot, akm * q * ce * edot, 0.0];
    // perifocal -> ecliptic of the equinox: Rz(-Ω) Rx(-i) Rz(-ω)
    let r = mat_mul(&rot_z(-el.node * DEG), &mat_mul(&rot_x(-el.i * DEG), &rot_z(-el.w * DEG)));
    let (p, v) = (mat_vec(&r, pos), mat_vec(&r, vel));
    let eq = match el.equinox {
        Equinox::Fixed(j) => j,
        Equinox::OfDate => jd,
    };
    let m = ecliptic_matrix(eq, model);
    (mat_t_vec(&m, p), mat_t_vec(&m, v))
}

/// Barycentric GCRS state (km, km/day) at JD TDB.
pub fn barycentric_state(pack: &Pack, el: &Elements, jd_tdb: f64, model: FrameModel) -> Result<(Vec3, Vec3)> {
    let center = match el.center {
        Center::Sun => crate::ephem::SUN,
        Center::Earth => crate::ephem::EARTH,
    };
    let (c, cv) = pack.state(center, crate::ephem::SSB, jd_tdb)?;
    let (p, v) = relative_state(el, jd_tdb, model);
    Ok((add(c, p), add(cv, v)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geocentric_circles_match_published_sizes() {
        // Waltemath: 2.67 lunar distances, ~3°/day (Swiss docs §2.7.7);
        // Selena: > 20 lunar distances (§2.7.6).
        let w = semi_axis(elements(58).unwrap()) * AU_KM / 384_400.0;
        assert!((w - 2.66).abs() < 0.02, "{w}");
        assert!((elements(58).unwrap().n.unwrap() - 3.02).abs() < 0.01);
        let s = semi_axis(elements(56).unwrap()) * AU_KM / 384_400.0;
        assert!(s > 20.0 && s < 21.0, "{s}");
    }

    #[test]
    fn neely_mean_motion_is_keplerian() {
        // Neely's published M rates (deg/century) agree with Kepler's law to < 0.01 %.
        for (id, rate) in [(40, 137.13640), (41, 99.81803), (44, 61.0765), (47, 47.03868)] {
            let el = elements(id).unwrap();
            let n = mean_motion(el, el.a) / DEG * 36525.0;
            assert!((n / rate - 1.0).abs() < 1e-4, "{id}: {n} vs {rate}");
        }
    }

    #[test]
    fn lowell_true_longitude_at_epoch() {
        let el = elements(53).unwrap();
        let (p, _) = relative_state(el, el.epoch, FrameModel::DEFAULT);
        let m = ecliptic_matrix(el.epoch, FrameModel::DEFAULT);
        let (l, _, _) = to_spherical(mat_vec(&m, p));
        assert!((norm_deg(l / DEG) - 84.0).abs() < 1e-9);
    }
}
