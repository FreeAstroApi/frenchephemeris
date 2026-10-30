// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Engine-level house functions: houses with speeds for an instant and place
//! (Sunshine's Sun declination computed from the ephemeris), and Gauquelin
//! sector positions (`gauquelin_sector`, six methods).

use super::speeds::{house_cusps_speeds, house_cusps_speeds_swiss_compat, HouseSpeeds};
use super::systems::{house_cusps, house_cusps_swiss_compat};
use crate::ayanamsa::{self, Ayanamsa, Projection, SidBits};
use crate::context::{Body, CalcOptions, Engine};
use crate::ephem::Epoch;
use crate::ephem::Observer;
use crate::events::{Atmosphere, DiscPoint, RiseSet, RiseSetKind, RiseSetOptions};
use crate::houses;
use crate::math::{cross, dot, mat_mul, mat_t_vec, mat_vec, norm, rot_x, rot_z, scale};
use crate::math::{norm_deg, DEG};
use crate::{Error, Result};

/// Gauquelin "plus zones" (sector ranges found significant by M. Gauquelin,
/// as listed in the Swiss Ephemeris documentation §6.2.12): 36–3, 9–12,
/// 19–21, 28–30 (inclusive, clockwise).
pub const GAUQUELIN_PLUS_ZONES: [(u8, u8); 4] = [(36, 3), (9, 12), (19, 21), (28, 30)];

/// Method of [`Engine::gauquelin_sector`] (the classic `imeth` 0–5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GauquelinMethod {
    /// 0: Placidus-type sector of the ecliptic position with its latitude.
    WithLatitude,
    /// 1: same, ecliptic latitude set to 0.
    WithoutLatitude,
    /// 2: from rise and set times of the disc centre, no refraction.
    RiseSetCenter,
    /// 3: disc centre, with refraction.
    RiseSetCenterRefracted,
    /// 4: upper limb, no refraction.
    RiseSetLimb,
    /// 5: upper limb, with refraction.
    RiseSetLimbRefracted,
}

impl GauquelinMethod {
    pub fn from_code(code: i32) -> Result<GauquelinMethod> {
        use GauquelinMethod::*;
        Ok(match code {
            0 => WithLatitude,
            1 => WithoutLatitude,
            2 => RiseSetCenter,
            3 => RiseSetCenterRefracted,
            4 => RiseSetLimb,
            5 => RiseSetLimbRefracted,
            _ => return Err(Error::InvalidArgument(format!("Gauquelin method must be 0..5, got {code}"))),
        })
    }
}

/// Classic house-system code for Sunshine needs the Sun's declination.
fn is_sunshine(system: char) -> bool {
    system == 'I' || system == 'i'
}

impl Engine {
    /// ARMC (deg; apparent sidereal time + east longitude) and true obliquity
    /// (deg) at a UT instant, as `houses_ut` uses them.
    pub fn house_frame_ut(&self, jd_ut: f64, lon_deg: f64) -> (f64, f64) {
        let ep = self.prepare_ut(jd_ut, false).epoch;
        (norm_deg(ep.gast() / DEG + lon_deg), ep.frames.eps_true / DEG)
    }

    /// Apparent geocentric declination of the Sun (true equator of date),
    /// deg — the input of the Sunshine house systems (oracle semantics).
    pub fn sun_declination_ut(&self, jd_ut: f64) -> Result<f64> {
        let opt = CalcOptions { equatorial: true, speed: false, ..Default::default() };
        Ok(self.calc_ut(Body::Sun, jd_ut, &opt)?.lat)
    }

    /// The Sun's declination if `system` is a Sunshine system, else `None`.
    pub fn sunshine_declination(&self, system: char, jd_ut: f64) -> Result<Option<f64>> {
        if is_sunshine(system) {
            Ok(Some(self.sun_declination_ut(jd_ut)?))
        } else {
            Ok(None)
        }
    }

    /// Any house system (incl. 'G' and Sunshine) with cusp and angle speeds
    /// (deg/day) at an instant and place: the classic `houses_ex2`.
    /// `swiss_fallback` reproduces the oracle's Porphyry fallback (second
    /// value `true`) instead of returning `PolarHouses`. Sidereal output: see
    /// [`Engine::house_cusps_ut_bits`].
    pub fn house_cusps_ut(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        system: char,
        sidereal: Option<Ayanamsa>,
        swiss_fallback: bool,
    ) -> Result<(HouseSpeeds, bool)> {
        self.house_cusps_ut_bits(
            jd_ut,
            lat_deg,
            lon_deg,
            system,
            sidereal.map(|m| (m, SidBits(0))),
            swiss_fallback,
            true,
        )
    }

    /// [`Engine::house_cusps_ut`] with sidereal option bits; `speeds = false`
    /// skips the speed differences (speeds returned as 0).
    ///
    /// Sidereal, traditional modes: the true ayanamsa is subtracted from cusps
    /// and angles except ARMC. Projection modes (`SIDBIT_ECL_T0`,
    /// `SIDBIT_SSY_PLANE`, and the modes that project by default): the house
    /// circles are intersected with the plane of the projection instead of
    /// the ecliptic of date — houses are computed for that plane's inclination
    /// to the true equator of date and ARMC counted from its ascending node Ω
    /// on that equator, then the node's sidereal longitude is added;
    /// `ascmc[2]` is ARMC − Ω (oracle semantics, measured). Speeds stay
    /// tropical (the oracle does not shift them either).
    #[allow(clippy::too_many_arguments)]
    pub fn house_cusps_ut_bits(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        system: char,
        sidereal: Option<(Ayanamsa, SidBits)>,
        swiss_fallback: bool,
        speeds: bool,
    ) -> Result<(HouseSpeeds, bool)> {
        let (mut armc, mut eps) = self.house_frame_ut(jd_ut, lon_deg);
        let sd = self.sunshine_declination(system, jd_ut)?;
        let ep = self.prepare_ut(jd_ut, false).epoch;
        // (shift added to the longitudes, ARMC override)
        let mut shift = None;
        if let Some((mode, bits)) = sidereal {
            let proj = ayanamsa::projection(mode, bits);
            let proj_np = ayanamsa::projection(mode, SidBits(bits.0 | SidBits::NO_PREC_OFFSET));
            if proj == Projection::Traditional {
                shift = Some((-self.ayanamsa_deg_bits(mode, bits, &ep, true)?, None));
            } else {
                let (e2, om, lnode) = projection_frame(&ep, &proj_np);
                armc = norm_deg(armc - om);
                eps = e2;
                shift = Some((lnode, Some(armc)));
            }
        }
        let (mut h, fb) = match (speeds, swiss_fallback) {
            (true, true) => house_cusps_speeds_swiss_compat(system, armc, lat_deg, eps, sd)?,
            (true, false) => (house_cusps_speeds(system, armc, lat_deg, eps, sd)?, false),
            (false, fallback) => {
                let (c, fb) = if fallback {
                    house_cusps_swiss_compat(system, armc, lat_deg, eps, sd)?
                } else {
                    (house_cusps(system, armc, lat_deg, eps, sd)?, false)
                };
                let n = c.cusps.len();
                (HouseSpeeds { houses: c, cusp_speeds: vec![0.0; n], ascmc_speeds: [0.0; 8] }, fb)
            }
        };
        if let Some((d, armc_out)) = shift {
            let n = h.houses.cusps.len();
            for (i, c) in h.houses.cusps.iter_mut().enumerate() {
                // the G fallback's zero padding stays zero
                if !(fb && n == 36 && i >= 12) {
                    *c = norm_deg(*c + d);
                }
            }
            for (i, x) in h.houses.ascmc.iter_mut().enumerate() {
                if i != 2 {
                    *x = norm_deg(*x + d);
                }
            }
            if let Some(a) = armc_out {
                h.houses.ascmc[2] = a;
            }
            // Zodiac-anchored systems (oracle semantics, measured): 'N' starts
            // at sidereal 0° Aries in every mode; 'W' takes the sidereal sign of
            // the Asc in traditional modes (in projection modes the oracle keeps
            // the signs of the projection plane's own construction, as above).
            match houses::canonical_system(system)? {
                'N' => {
                    for (i, c) in h.houses.cusps.iter_mut().enumerate() {
                        *c = 30.0 * i as f64;
                    }
                }
                'W' if armc_out.is_none() => {
                    let base = (h.houses.ascmc[0] / 30.0).floor() * 30.0;
                    for (i, c) in h.houses.cusps.iter_mut().enumerate() {
                        *c = norm_deg(base + 30.0 * i as f64);
                    }
                }
                _ => {}
            }
        }
        Ok((h, fb))
    }

    /// Gauquelin sector position (1 ≤ s < 37, clockwise from the Asc) of a
    /// body at `jd_ut`, the classic `gauquelin_sector`:
    /// - methods 0/1: `house_position('G')` of the body's ecliptic place of
    ///   date (computed with `opt`'s origin and corrections; method 1 drops the
    ///   ecliptic latitude) for ARMC = LAST and the true obliquity;
    /// - methods 2–5: from the previous and next rise and set (disc centre or
    ///   upper limb, with or without refraction for `atm`): above the horizon
    ///   s = 1 + 18 (t − rise)/(set − rise), below s = 19 + 18 (t − set)/(rise − set).
    ///   A body that does not both rise and set within the search window
    ///   (circumpolar, polar day/night) is an error, as in the oracle.
    pub fn gauquelin_sector(
        &self,
        jd_ut: f64,
        body: Body,
        method: GauquelinMethod,
        obs: &Observer,
        atm: Atmosphere,
        opt: &CalcOptions,
    ) -> Result<f64> {
        use GauquelinMethod::*;
        match method {
            WithLatitude | WithoutLatitude => {
                let o = CalcOptions {
                    equatorial: false,
                    j2000: false,
                    mean_equinox: false,
                    sidereal: None,
                    speed: false,
                    ..*opt
                };
                let p = self.calc_ut(body, jd_ut, &o)?;
                let (armc, eps) = self.house_frame_ut(jd_ut, obs.lon_deg);
                let blat = if method == WithLatitude { p.lat } else { 0.0 };
                super::house_position_ex('G', armc, obs.lat_deg, eps, p.lon, blat, None)
            }
            _ => {
                let ro = RiseSetOptions {
                    disc: if matches!(method, RiseSetCenter | RiseSetCenterRefracted) {
                        DiscPoint::Center
                    } else {
                        DiscPoint::UpperLimb
                    },
                    refraction: if matches!(method, RiseSetCenterRefracted | RiseSetLimbRefracted) {
                        Some(atm)
                    } else {
                        None
                    },
                    ..Default::default()
                };
                let circumpolar = || {
                    Error::InvalidArgument(format!("{body:?} does not rise and set at this place and date"))
                };
                let next = |t: f64, kind: RiseSetKind| -> Result<Option<f64>> {
                    Ok(match self.rise_set(body, t, kind, obs, &ro)? {
                        RiseSet::Event(e) => Some(e),
                        _ => None,
                    })
                };
                // latest event of `kind` strictly before t
                let last_before = |t: f64, kind: RiseSetKind| -> Result<f64> {
                    let mut e = next(t - 1.2, kind)?.ok_or_else(circumpolar)?;
                    if e >= t {
                        // nothing in [t − 1.2, t): look further back
                        e = next(t - 2.4, kind)?.ok_or_else(circumpolar)?;
                        if e >= t {
                            return Err(circumpolar());
                        }
                    }
                    while let Some(e2) = next(e + 1e-4, kind)? {
                        if e2 >= t {
                            break;
                        }
                        e = e2;
                    }
                    Ok(e)
                };
                let rise = next(jd_ut, RiseSetKind::Rise)?.ok_or_else(circumpolar)?;
                let set = next(jd_ut, RiseSetKind::Set)?.ok_or_else(circumpolar)?;
                Ok(if set < rise {
                    let r0 = last_before(jd_ut, RiseSetKind::Rise)?;
                    1.0 + 18.0 * (jd_ut - r0) / (set - r0)
                } else {
                    let s0 = last_before(jd_ut, RiseSetKind::Set)?;
                    19.0 + 18.0 * (jd_ut - s0) / (rise - s0)
                })
            }
        }
    }
}

/// Plane of a sidereal projection seen from the true equator of date:
/// (inclination ε′, RA Ω of its ascending node, sidereal longitude of that
/// node), degrees. The plane's pole in the J2000 frame follows from the
/// projection itself (sin of the projected latitude of the basis vectors).
fn projection_frame(ep: &Epoch, proj: &Projection) -> (f64, f64, f64) {
    let lat_of = |v: [f64; 3]| proj.apply(v).expect("projection").1 * DEG;
    let pole = [lat_of([1.0, 0.0, 0.0]).sin(), lat_of([0.0, 1.0, 0.0]).sin(), lat_of([0.0, 0.0, 1.0]).sin()];
    // J2000 mean equatorial (Vondrák frame) -> true equator and equinox of date
    let m = mat_mul(
        &mat_mul(&rot_x(-ep.frames.eps_true), &rot_z(-ep.frames.dpsi)),
        &ayanamsa::vondrak::ecliptic_of_date_matrix(ep.jd_tt),
    );
    let n = mat_vec(&m, pole);
    let n = scale(n, 1.0 / norm(n));
    let inc = n[2].clamp(-1.0, 1.0).acos() / DEG;
    let node = cross([0.0, 0.0, 1.0], n);
    let node = scale(node, 1.0 / norm(node));
    let om = node[1].atan2(node[0]) / DEG;
    let lnode = proj.apply(mat_t_vec(&m, node)).expect("projection").0;
    debug_assert!(dot(node, n).abs() < 1e-12);
    (inc, om, lnode)
}
