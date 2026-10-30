// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Cusp and angle speeds (`houses_armc_ex2` / `houses_ex2`), deg/day.
//!
//! Semantics (measured on the oracle): the speeds are the derivatives of the
//! cusps and angles with respect to ARMC times the rate of sidereal time,
//! with the obliquity (and, for Sunshine, the Sun's declination) held fixed;
//! the ARMC speed itself is that rate. Computed here by central differences
//! (±1e-3° of ARMC); across a discontinuity (Whole-sign or Equal/Aries
//! cusps, Asc hemisphere flip) the one-sided difference on the continuous side
//! is used, so constant cusps get speed 0.

use super::systems::{house_cusps, house_cusps_swiss_compat, HouseCusps};
use crate::math::{diff_deg, norm_deg};
use crate::Result;

/// Rate of sidereal time, deg/day of UT: 360 × 1.002737909350795, the
/// ratio of mean sidereal to UT1 time at J2000 (Aoki et al. 1982, A&A 105,
/// 359; Explanatory Supplement 1992 eq. 2.24-1). Used as dARMC/dt.
pub const ARMC_RATE_DEG_PER_DAY: f64 = 360.0 * 1.002_737_909_350_795;

/// ARMC half-step (deg) of the central differences.
const STEP_DEG: f64 = 1e-3;

/// Cusps, angles and their speeds (deg/day), Swiss `houses_armc_ex2` layout.
#[derive(Debug, Clone, PartialEq)]
pub struct HouseSpeeds {
    pub houses: HouseCusps,
    /// Same length as `houses.cusps` (12, or 36 for 'G').
    pub cusp_speeds: Vec<f64>,
    pub ascmc_speeds: [f64; 8],
}

fn derivative(m: f64, c: f64, p: f64) -> f64 {
    let lo = diff_deg(c, m);
    let hi = diff_deg(p, c);
    let d = if (hi - lo).abs() > 1.0 {
        // a jump on one side: take the continuous (smaller) one-sided slope
        if hi.abs() < lo.abs() {
            hi / STEP_DEG
        } else {
            lo / STEP_DEG
        }
    } else {
        diff_deg(p, m) / (2.0 * STEP_DEG)
    };
    d * ARMC_RATE_DEG_PER_DAY
}

fn speeds_from(c: HouseCusps, m: &HouseCusps, p: &HouseCusps) -> HouseSpeeds {
    let cusp_speeds = (0..c.cusps.len()).map(|i| derivative(m.cusps[i], c.cusps[i], p.cusps[i])).collect();
    let mut ascmc_speeds: [f64; 8] = std::array::from_fn(|i| derivative(m.ascmc[i], c.ascmc[i], p.ascmc[i]));
    ascmc_speeds[2] = ARMC_RATE_DEG_PER_DAY;
    HouseSpeeds { houses: c, cusp_speeds, ascmc_speeds }
}

/// Cusps, angles and speeds of any system (see [`house_cusps`]).
pub fn house_cusps_speeds(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<HouseSpeeds> {
    let c = house_cusps(system, armc_deg, lat_deg, eps_deg, sun_dec_deg)?;
    let m = house_cusps(system, norm_deg(armc_deg - STEP_DEG), lat_deg, eps_deg, sun_dec_deg)?;
    let p = house_cusps(system, norm_deg(armc_deg + STEP_DEG), lat_deg, eps_deg, sun_dec_deg)?;
    Ok(speeds_from(c, &m, &p))
}

/// Swiss-compatible [`house_cusps_speeds`] (Porphyry fallback, see
/// [`house_cusps_swiss_compat`]); returns `(speeds, fell_back)`.
pub fn house_cusps_speeds_swiss_compat(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<(HouseSpeeds, bool)> {
    let (c, fb) = house_cusps_swiss_compat(system, armc_deg, lat_deg, eps_deg, sun_dec_deg)?;
    let (m, _) =
        house_cusps_swiss_compat(system, norm_deg(armc_deg - STEP_DEG), lat_deg, eps_deg, sun_dec_deg)?;
    let (p, _) =
        house_cusps_swiss_compat(system, norm_deg(armc_deg + STEP_DEG), lat_deg, eps_deg, sun_dec_deg)?;
    Ok((speeds_from(c, &m, &p), fb))
}
