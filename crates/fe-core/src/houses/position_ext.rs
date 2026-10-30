// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! `house_position` for the systems of `systems.rs`.
//!
//! Geometric positions ("house degree" mapped linearly between the house
//! circles of the system's own construction) for X F M B U H G Y I i;
//! linear in longitude between the cusps for A V D N L Q; Sripati = Porphyry
//! position + ½ house. Which form each system uses was chosen to agree with
//! the oracle's `house_pos` (measured; see docs/GEOMETRY_REPORT.md §1).

use super::systems::{
    ecl_dec, ecl_vec, eq_vec, equator_point_ra, horizon_frame, lon_to_ra, position_circle_k, semi_arc,
};
use super::{angles, clamp1, clamp_lat, ecl_to_eq, hdeg_to_pos, houses_ex, linear_pos, placidus_hdeg};
use crate::math::{diff_deg, dot, norm_deg, DEG};
use crate::{Error, Result};

pub(super) fn house_position_extra(
    sys: char,
    armc_deg: f64,
    lat_deg: f64,
    eps: f64,
    lon: f64,
    blat: f64,
    sun_dec: Option<f64>,
) -> Result<f64> {
    let lat = clamp_lat(lat_deg);
    let armc = norm_deg(armc_deg);
    let lon = norm_deg(lon);
    if matches!(sys, 'A' | 'V' | 'D' | 'N' | 'L' | 'Q') {
        let h = houses_ex(sys, armc, lat, eps, None)?;
        return Ok(linear_pos(&h.cusps, lon));
    }
    if sys == 'J' {
        // linear in the Campanus position between the cusps' Campanus
        // positions (oracle form, measured; the cusps lie on Campanus-type
        // circles, so this is the geometric position on that circle family)
        let h = houses_ex('J', armc, lat, eps, None)?;
        let cpos = |l: f64, b: f64| super::house_position('C', armc, lat, eps, l, b);
        let mut m = [0.0; 12];
        for (mi, &c) in m.iter_mut().zip(h.cusps.iter()) {
            *mi = (cpos(c, 0.0)? - 1.0) * 30.0;
        }
        let p = (cpos(lon, blat)? - 1.0) * 30.0;
        return Ok(linear_pos(&m, p));
    }
    if sys == 'S' {
        let p = super::house_position('O', armc, lat, eps, lon, blat)? + 0.5;
        return Ok(if p >= 13.0 { p - 12.0 } else { p });
    }
    let (ra, dec) = ecl_to_eq(lon, blat, eps);
    let a = angles(armc, lat, eps);
    let hdeg = match sys {
        // Meridian: house circles are hour circles at RA = ARMC + 30k.
        'X' => ra - armc,
        // Carter: hour circles at RA(Asc) + 30k.
        'F' => ra - lon_to_ra(a.asc, eps) + 90.0,
        // Morinus: circles through the ecliptic poles (latitude ignored).
        'M' => equator_point_ra(lon, eps) - armc,
        'B' => alcabitius_hdeg(armc, lat, eps, a.asc, ra),
        'U' => krusinski_hdeg(armc, lat, eps, a.asc, ra, dec),
        'H' => {
            let (_, s, e) = horizon_frame(armc, lat);
            let p = eq_vec(ra, dec);
            dot(p, e).atan2(dot(p, s)) / DEG
        }
        'G' => {
            let h = placidus_hdeg(armc, lat, ra, dec);
            // sectors clockwise: Asc (hdeg 90) = 1, MC (0) = 10, Desc = 19, IC = 28
            return Ok(norm_deg(90.0 - h) / 10.0 + 1.0);
        }
        'Y' => parallel_hdeg(armc, lat, ra, dec, ecl_dec(a.asc, eps)),
        'I' | 'i' => {
            let sd = sun_dec.ok_or_else(|| {
                Error::InvalidArgument("Sunshine house position needs the Sun's declination".into())
            })?;
            let t = (lat * DEG).tan() * (sd * DEG).tan();
            if t >= 1.0 {
                180.0 // no nocturnal arc: every point maps to the IC (oracle behaviour)
            } else if t <= -1.0 {
                0.0 // no diurnal arc: every point maps to the MC
            } else {
                parallel_hdeg(armc, lat, ra, dec, sd)
            }
        }
        _ => return Err(Error::InvalidArgument(format!("unknown house system {sys:?}"))),
    };
    Ok(hdeg_to_pos(hdeg))
}

/// Alcabitius: linear in RA between the cusp hour circles (the Asc's
/// semi-arcs trisected on the equator; west of the meridian the Desc's).
fn alcabitius_hdeg(armc: f64, lat: f64, eps: f64, asc: f64, ra: f64) -> f64 {
    let dsa = semi_arc(lat, ecl_dec(asc, eps));
    let nsa = 180.0 - dsa;
    let d = norm_deg(ra - armc);
    let seg = |x: f64, w: f64| if w > 0.0 { 90.0 * x / w } else { 0.0 };
    if d < dsa {
        seg(d, dsa)
    } else if d < 180.0 {
        90.0 + seg(d - dsa, nsa)
    } else if d < 180.0 + dsa {
        180.0 + seg(d - 180.0, dsa)
    } else {
        270.0 + seg(d - 180.0 - dsa, nsa)
    }
}

/// Krusinski: the point's hour circle meets the Asc–zenith circle at angle t
/// from the Asc (towards the zenith); house degree 90 − t.
fn krusinski_hdeg(armc: f64, lat: f64, eps: f64, asc: f64, ra: f64, dec: f64) -> f64 {
    let av = ecl_vec(asc, eps);
    let (z, _, _) = horizon_frame(armc, lat);
    let p = eq_vec(ra, dec);
    let n = [-p[1], p[0], 0.0]; // pole × P: normal of the hour-circle plane
    let mut t = (-dot(n, av)).atan2(dot(n, z)) / DEG;
    let (st, ct) = (t * DEG).sin_cos();
    let x = [av[0] * ct + z[0] * st, av[1] * ct + z[1] * st];
    if x[0] * p[0] + x[1] * p[1] < 0.0 {
        t += 180.0;
    }
    90.0 - t
}

/// APC / Sunshine: the point's circle through the N/S horizon points cuts the
/// reference parallel (declination `pdec`) at hour-angle offset d; house
/// degree = d as a fraction of that parallel's semi-arcs.
fn parallel_hdeg(armc: f64, lat: f64, ra: f64, dec: f64, pdec: f64) -> f64 {
    let ta = (lat * DEG).tan() * (pdec * DEG).tan();
    let dsa = 90.0 + clamp1(ta).asin() / DEG;
    let nsa = 180.0 - dsa;
    let k = position_circle_k(ra - armc, dec, lat);
    let d = diff_deg(k + clamp1(ta * (k * DEG).sin()).asin() / DEG, 0.0);
    if d.abs() <= dsa {
        if dsa > 0.0 {
            90.0 * d / dsa
        } else {
            0.0
        }
    } else if d > 0.0 {
        90.0 + 90.0 * (d - dsa) / nsa
    } else {
        -90.0 - 90.0 * (-d - dsa) / nsa
    }
}
