// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! House systems. Pure geometry from (ARMC, latitude, obliquity).
//!
//! The eight original systems live in this file; the other Swiss-documented
//! systems (A V D N B X F M H U Y L Q S I i G J), cusp speeds, Gauquelin sectors
//! and house names live in the `houses/` submodules (see docs/GEOMETRY_REPORT.md).
//!
//! CONTRACT (owned by the astrology-geometry workstream):
//! - Inputs and outputs in DEGREES here (matches the astrology domain and the
//!   oracle), longitudes normalized to [0, 360).
//! - System codes follow the widely used single-letter convention:
//!   'P' Placidus, 'W' Whole sign, 'E' Equal (from Asc), 'K' Koch,
//!   'O' Porphyry, 'R' Regiomontanus, 'C' Campanus, 'T' Topocentric (Polich/Page).
//! - Undefined geometry (Placidus/Koch above polar circle) returns
//!   `Error::PolarHouses`; the adapter decides fallback semantics
//!   (see [`houses_swiss_compat`] for the Swiss-compatible Porphyry fallback).
//!
//! Definitions and the edge rules measured against the Swiss oracle are
//! documented in `docs/GEOMETRY_REPORT.md`. Summary of the edge rules:
//! - |lat| is clamped to `90 - 1e-10` degrees (Swiss does the same, so the
//!   poles are not singular in the output).
//! - Placidus and Koch are undefined when `|lat| >= 90 - eps`.
//! - The Ascendant is always the point of the ecliptic on the eastern horizon
//!   that lies in the half-circle `(MC, MC + 180)`; inside the polar circle the
//!   raw rising point can fall in the western half and is then replaced by
//!   its opposite. For Regiomontanus, Campanus and Topocentric, when that flip
//!   happens all cusps except 1 and 7 are also turned by 180 (so cusp 10 and
//!   `mc` become the IC). Porphyry/Equal/Whole use the flipped Asc only.
//! - The Vertex is always in the western half `(MC - 180, MC)`.
//! - When the ecliptic coincides exactly with the horizon (only possible for
//!   `|lat| == 90 - eps` and ARMC = 90/270) the Asc is taken as MC + 90.

use crate::math::{diff_deg, norm_deg, DEG};
use crate::{Error, Result};

mod gauquelin;
mod names;
mod position_ext;
mod speeds;
mod systems;

pub use gauquelin::{GauquelinMethod, GAUQUELIN_PLUS_ZONES};
pub use names::house_name;
pub use speeds::{house_cusps_speeds, house_cusps_speeds_swiss_compat, HouseSpeeds, ARMC_RATE_DEG_PER_DAY};
pub use systems::{gauquelin_cusps, house_cusps, house_cusps_swiss_compat, HouseCusps};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Houses {
    /// cusps[0] = house 1 ... cusps[11] = house 12, degrees.
    pub cusps: [f64; 12],
    pub asc: f64,
    pub mc: f64,
    pub armc: f64,
    pub vertex: f64,
    pub equatorial_asc: f64,
    pub co_asc_koch: f64,
    pub co_asc_munkasey: f64,
    pub polar_asc: f64,
}

impl Houses {
    /// The 8 angles in the order of Swiss' `ascmc` array:
    /// Asc, MC, ARMC, Vertex, equatorial Asc, co-Asc (Koch), co-Asc (Munkasey), polar Asc.
    pub fn ascmc(&self) -> [f64; 8] {
        [
            self.asc,
            self.mc,
            self.armc,
            self.vertex,
            self.equatorial_asc,
            self.co_asc_koch,
            self.co_asc_munkasey,
            self.polar_asc,
        ]
    }
}

/// Numerical guard: latitudes are clamped to this magnitude (degrees) so the
/// poles are not singular (tan φ finite). Behaviour-level choice; the same
/// guard magnitude is observed in the oracle.
pub const MAX_LAT_DEG: f64 = 90.0 - 1e-10;

/// All system codes implemented (Swiss letters; case-insensitive except
/// 'I' Sunshine (variant I) vs 'i' Sunshine/Makransky). 'G' (Gauquelin, 36
/// sectors) is served by [`house_cusps`] / [`gauquelin_cusps`], not by
/// [`houses`]; 'I'/'i' need the Sun's declination ([`houses_ex`]).
pub const SYSTEMS: [char; 26] = [
    'P', 'K', 'O', 'R', 'C', 'T', 'E', 'W', 'A', 'V', 'D', 'N', 'B', 'X', 'F', 'M', 'H', 'U', 'Y', 'L', 'Q',
    'S', 'I', 'i', 'G', 'J',
];

/// Canonical system code: uppercase, except 'i' (Makransky Sunshine).
/// Errors for unknown codes.
pub fn canonical_system(system: char) -> Result<char> {
    let s = if system == 'i' { 'i' } else { system.to_ascii_uppercase() };
    if !SYSTEMS.contains(&s) {
        return Err(Error::InvalidArgument(format!("unknown house system {system:?}")));
    }
    Ok(s)
}

#[inline]
fn clamp_lat(lat: f64) -> f64 {
    lat.clamp(-MAX_LAT_DEG, MAX_LAT_DEG)
}

#[inline]
fn clamp1(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// (sin, cos) of an angle in degrees, exact at multiples of 90°.
#[inline]
fn sin_cos_deg(x: f64) -> (f64, f64) {
    let r = x.rem_euclid(360.0);
    if r == 0.0 {
        (0.0, 1.0)
    } else if r == 90.0 {
        (1.0, 0.0)
    } else if r == 180.0 {
        (0.0, -1.0)
    } else if r == 270.0 {
        (-1.0, 0.0)
    } else {
        (r * DEG).sin_cos()
    }
}

/// atan2 in degrees normalized to [0, 360), exact on the axes.
#[inline]
fn atan2_deg(y: f64, x: f64) -> f64 {
    if y == 0.0 && x > 0.0 {
        0.0
    } else if y == 0.0 && x < 0.0 {
        180.0
    } else if x == 0.0 && y > 0.0 {
        90.0
    } else if x == 0.0 && y < 0.0 {
        270.0
    } else {
        norm_deg(y.atan2(x) / DEG)
    }
}

const HALF_PI: f64 = std::f64::consts::FRAC_PI_2;
const PI: f64 = std::f64::consts::PI;

/// `true` if Placidus/Koch are undefined at this latitude/obliquity
/// (Swiss rule, measured by bisection: `|lat| >= 90 - eps`).
#[inline]
pub fn is_polar(lat_deg: f64, eps_deg: f64) -> bool {
    lat_deg.abs() >= 90.0 - eps_deg
}

/// Ecliptic longitude of the MC (degrees) for a given ARMC.
#[inline]
pub fn mc_deg(armc_deg: f64, eps_deg: f64) -> f64 {
    let (s, c) = sin_cos_deg(armc_deg);
    atan2_deg(s, c * (eps_deg * DEG).cos())
}

/// Raw "ascendant" of an oblique horizon: the ecliptic point that rises on a
/// horizon of pole height `pole_deg` when the meridian has right ascension
/// `ramc_deg`. `None` in the exactly degenerate case where the ecliptic
/// coincides with that horizon.
#[inline]
fn asc_raw(ramc_deg: f64, pole_deg: f64, eps_deg: f64) -> Option<f64> {
    let (st, ct) = sin_cos_deg(ramc_deg);
    let (se, ce) = (eps_deg * DEG).sin_cos();
    let tp = if pole_deg == 0.0 { 0.0 } else { (pole_deg * DEG).tan() };
    let y = ct;
    let x = -(st * ce + tp * se);
    if y.abs() < 1e-12 && x.abs() < 1e-12 {
        return None;
    }
    Some(atan2_deg(y, x))
}

#[inline]
fn asc_or(ramc_deg: f64, pole_deg: f64, eps_deg: f64, fallback: f64) -> f64 {
    asc_raw(ramc_deg, pole_deg, eps_deg).unwrap_or(fallback)
}

/// Ecliptic longitude of the ecliptic point with right ascension `ra_deg`.
#[inline]
fn ra_to_lon(ra_deg: f64, eps_deg: f64) -> f64 {
    let (s, c) = sin_cos_deg(ra_deg);
    atan2_deg(s, c * (eps_deg * DEG).cos())
}

struct Angles {
    asc: f64,
    mc: f64,
    flipped: bool,
    vertex: f64,
    equatorial_asc: f64,
    co_asc_koch: f64,
    co_asc_munkasey: f64,
    polar_asc: f64,
}

fn angles(armc: f64, lat: f64, eps: f64) -> Angles {
    let mc = mc_deg(armc, eps);
    let east = norm_deg(mc + 90.0);
    let raw = asc_or(armc, lat, eps, east);
    let flipped = diff_deg(raw, mc) < 0.0;
    let asc = if flipped { norm_deg(raw + 180.0) } else { raw };
    let mut vertex = asc_or(armc + 180.0, 90.0 - lat, eps, norm_deg(mc - 90.0));
    let dv = diff_deg(vertex, mc);
    if dv > 0.0 && dv < 180.0 {
        vertex = norm_deg(vertex + 180.0);
    }
    let equatorial_asc = asc_or(armc, 0.0, eps, east);
    let co_asc_koch = asc_or(armc, -lat, eps, east);
    let co_asc_munkasey = asc_or(armc, 90.0 - lat, eps, east);
    let polar_asc = norm_deg(co_asc_koch + 180.0);
    Angles { asc, mc, flipped, vertex, equatorial_asc, co_asc_koch, co_asc_munkasey, polar_asc }
}

/// Placidus cusp (11, 12, 2 or 3): solves the semi-arc condition
/// `RA = RAMC + f·DSA(δ)` (above the horizon) or `RA = RAMC + 180 - f·NSA(δ)`
/// (below), δ being the declination of the ecliptic point with that RA.
/// Newton iterations on RA with an analytic derivative; fixed-point fallback.
fn placidus_cusp(armc: f64, lat: f64, eps: f64, house: u8) -> f64 {
    let (frac, upper) = match house {
        11 => (1.0 / 3.0, true),
        12 => (2.0 / 3.0, true),
        2 => (2.0 / 3.0, false),
        _ => (1.0 / 3.0, false),
    };
    placidus_frac(armc, lat, eps, frac, upper)
}

/// Ecliptic point that has completed the fraction `frac` of its semi-arc,
/// counted from the upper meridian (`upper`: diurnal arc, east of the MC) or
/// from the lower meridian (nocturnal arc, east of the IC). Placidus uses
/// 1/3, 2/3; Gauquelin sectors use k/9.
fn placidus_frac(armc: f64, lat: f64, eps: f64, frac: f64, upper: bool) -> f64 {
    let start = if upper { 90.0 * frac } else { 180.0 - 90.0 * frac };
    let tf = (lat * DEG).tan();
    let (se, ce) = (eps * DEG).sin_cos();
    let th = armc * DEG;
    // (target RA, d target / d RA), radians
    let target = |ra: f64| -> (f64, f64) {
        let (sr, cr) = ra.sin_cos();
        // ecliptic point with this RA: sin λ ∝ sin RA, cos λ ∝ cos RA · cos ε
        let n2 = sr * sr + cr * cr * ce * ce;
        let n = n2.sqrt();
        let (sl, cl) = (sr / n, cr * ce / n);
        let dlam = ce / n2; // dλ/dRA
        let sd = se * sl;
        let cd = (1.0 - sd * sd).sqrt();
        let s = clamp1(tf * sd / cd);
        let ad = s.asin();
        let root = (1.0 - s * s).sqrt();
        let dad = if root > 1e-300 { tf / (cd * cd) / root * (se * cl / cd) * dlam } else { 0.0 };
        if upper {
            (th + frac * (HALF_PI + ad), frac * dad)
        } else {
            (th + PI - frac * (HALF_PI - ad), frac * dad)
        }
    };
    let mut ra = th + start * DEG;
    let mut converged = false;
    for _ in 0..60 {
        let (t, dt) = target(ra);
        let g = t - ra;
        if g.abs() < 1e-15 {
            converged = true;
            break;
        }
        let dg = dt - 1.0;
        let mut step = if dg.abs() > 1e-9 { -g / dg } else { g };
        if !step.is_finite() {
            break;
        }
        // damped Newton: the target is smooth but can be steep near the polar circle
        step = step.clamp(-0.2, 0.2);
        ra += step;
        if step.abs() < 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        ra = th + start * DEG;
        for _ in 0..5000 {
            let (t, _) = target(ra);
            let d = t - ra;
            ra = t;
            if d.abs() < 1e-15 {
                break;
            }
        }
    }
    ra_to_lon(ra / DEG, eps)
}

/// Diurnal semi-arc (degrees) of the MC degree.
#[inline]
fn mc_semi_arc(armc: f64, lat: f64, eps: f64) -> f64 {
    let mc = mc_deg(armc, eps);
    let dec_mc = ((eps * DEG).sin() * (mc * DEG).sin()).asin();
    90.0 + clamp1((lat * DEG).tan() * dec_mc.tan()).asin() / DEG
}

/// Koch: the cusps are the ascendants at the times that trisect the MC
/// degree's diurnal semi-arc D: RAMC - 2D/3, RAMC - D/3, RAMC + D/3, RAMC + 2D/3.
fn koch_cusps(armc: f64, lat: f64, eps: f64, fb: f64) -> [f64; 4] {
    let d = mc_semi_arc(armc, lat, eps);
    [
        asc_or(armc - 2.0 * d / 3.0, lat, eps, fb),
        asc_or(armc - d / 3.0, lat, eps, fb),
        asc_or(armc + d / 3.0, lat, eps, fb),
        asc_or(armc + 2.0 * d / 3.0, lat, eps, fb),
    ]
}

/// (cusp index, equator offset from the upper meridian) for houses 11, 12, 2, 3.
const QUAD: [(usize, f64); 4] = [(10, 30.0), (11, 60.0), (1, 120.0), (2, 150.0)];

/// Regiomontanus: equator trisected from the meridian; pole tan p = tan φ · sin H.
fn regio_cusps(armc: f64, lat: f64, eps: f64, fb: f64) -> [f64; 4] {
    let tf = (lat * DEG).tan();
    let mut r = [0.0; 4];
    for (i, &(_, k)) in QUAD.iter().enumerate() {
        let pole = (tf * (k * DEG).sin()).atan() / DEG;
        r[i] = asc_or(armc + k - 90.0, pole, eps, fb);
    }
    r
}

/// Topocentric (Polich/Page): equator offsets 30/60 with pole tan p = tan φ · n/3.
fn topo_cusps(armc: f64, lat: f64, eps: f64, fb: f64) -> [f64; 4] {
    let tf = (lat * DEG).tan();
    let frac = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0];
    let mut r = [0.0; 4];
    for (i, &(_, k)) in QUAD.iter().enumerate() {
        let pole = (tf * frac[i]).atan() / DEG;
        r[i] = asc_or(armc + k - 90.0, pole, eps, fb);
    }
    r
}

/// Campanus: prime vertical trisected; the house circle through prime-vertical
/// azimuth k has pole sin p = sin φ · sin k and meets the equator at
/// tan H = cos φ · tan k.
fn campanus_cusps(armc: f64, lat: f64, eps: f64, fb: f64) -> [f64; 4] {
    let (sf, cf) = (lat * DEG).sin_cos();
    let mut r = [0.0; 4];
    for (i, &(_, k)) in QUAD.iter().enumerate() {
        let (sk, ck) = (k * DEG).sin_cos();
        let pole = clamp1(sf * sk).asin() / DEG;
        let h = (sk * cf).atan2(ck) / DEG;
        r[i] = asc_or(armc + h - 90.0, pole, eps, fb);
    }
    r
}

fn porphyry(c: &mut [f64; 12], asc: f64, mc: f64) {
    let a = norm_deg(asc - mc); // arc MC -> Asc
    let b = 180.0 - a;
    c[9] = mc;
    c[10] = norm_deg(mc + a / 3.0);
    c[11] = norm_deg(mc + 2.0 * a / 3.0);
    c[0] = asc;
    c[1] = norm_deg(asc + b / 3.0);
    c[2] = norm_deg(asc + 2.0 * b / 3.0);
    fill_opposites(c);
}

#[inline]
fn fill_opposites(c: &mut [f64; 12]) {
    for i in [9usize, 10, 11] {
        c[i - 6] = norm_deg(c[i] + 180.0);
    }
    for i in [0usize, 1, 2] {
        c[i + 6] = norm_deg(c[i] + 180.0);
    }
}

fn simple_cusps(sys: char, asc: f64, mc: f64) -> [f64; 12] {
    let mut c = [0.0f64; 12];
    match sys {
        'E' => {
            for (i, ci) in c.iter_mut().enumerate() {
                *ci = norm_deg(asc + 30.0 * i as f64);
            }
        }
        'W' => {
            let base = (asc / 30.0).floor() * 30.0;
            for (i, ci) in c.iter_mut().enumerate() {
                *ci = norm_deg(base + 30.0 * i as f64);
            }
        }
        _ => porphyry(&mut c, asc, mc),
    }
    c
}

/// Compute houses (12-cusp systems). See the module docs for the edge
/// semantics. Sunshine ('I', 'i') needs the Sun's declination: use
/// [`houses_ex`]; Gauquelin ('G') has 36 sectors: use [`house_cusps`].
pub fn houses(system: char, armc_deg: f64, lat_deg: f64, eps_deg: f64) -> Result<Houses> {
    houses_ex(system, armc_deg, lat_deg, eps_deg, None)
}

/// [`houses`] with the Sun's declination (degrees, true equator of date)
/// for the Sunshine systems; ignored by the other systems.
pub fn houses_ex(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<Houses> {
    let sys = canonical_system(system)?;
    if !matches!(sys, 'P' | 'K' | 'O' | 'R' | 'C' | 'T' | 'E' | 'W') {
        return systems::houses_extra(sys, armc_deg, lat_deg, eps_deg, sun_dec_deg);
    }
    let lat = clamp_lat(lat_deg);
    let armc = norm_deg(armc_deg);
    if (sys == 'P' || sys == 'K') && is_polar(lat, eps_deg) {
        return Err(Error::PolarHouses { system: sys, lat_deg });
    }
    let a = angles(armc, lat, eps_deg);
    let mut mc = a.mc;
    let asc = a.asc;
    let c = match sys {
        'E' | 'W' | 'O' => simple_cusps(sys, asc, mc),
        _ => {
            let fb = norm_deg(mc + 90.0);
            let q = match sys {
                'P' => [
                    placidus_cusp(armc, lat, eps_deg, 11),
                    placidus_cusp(armc, lat, eps_deg, 12),
                    placidus_cusp(armc, lat, eps_deg, 2),
                    placidus_cusp(armc, lat, eps_deg, 3),
                ],
                'K' => koch_cusps(armc, lat, eps_deg, fb),
                'R' => regio_cusps(armc, lat, eps_deg, fb),
                'T' => topo_cusps(armc, lat, eps_deg, fb),
                _ => campanus_cusps(armc, lat, eps_deg, fb),
            };
            let flip = a.flipped && matches!(sys, 'R' | 'T' | 'C');
            if flip {
                mc = norm_deg(mc + 180.0);
            }
            let turn = if flip { 180.0 } else { 0.0 };
            let mut c = [0.0f64; 12];
            for (i, &(idx, _)) in QUAD.iter().enumerate() {
                c[idx] = norm_deg(q[i] + turn);
            }
            c[0] = asc;
            c[9] = mc;
            fill_opposites(&mut c);
            c
        }
    };
    Ok(Houses {
        cusps: c,
        asc,
        mc,
        armc,
        vertex: a.vertex,
        equatorial_asc: a.equatorial_asc,
        co_asc_koch: a.co_asc_koch,
        co_asc_munkasey: a.co_asc_munkasey,
        polar_asc: a.polar_asc,
    })
}

/// Swiss-compatible wrapper: Placidus/Koch inside the polar circle fall back
/// to Porphyry. Returns `(houses, fell_back)`.
///
/// This mirrors the Swiss C library, which returns ERR together with
/// Porphyry cusps *and* Porphyry `ascmc` ("within polar circle, switched to
/// Porphyry"; verified identical through the C API). pyswisseph turns that
/// ERR into a raised `swe.Error`, so FreeAstroAPI today *fails* for P/K
/// above the polar circle — the adapter chooses which behaviour to expose.
pub fn houses_swiss_compat(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
) -> Result<(Houses, bool)> {
    houses_ex_swiss_compat(system, armc_deg, lat_deg, eps_deg, None)
}

/// [`houses_swiss_compat`] with the Sun's declination (Sunshine systems;
/// 'i' falls back to Porphyry when the Sun is circumpolar).
pub fn houses_ex_swiss_compat(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<(Houses, bool)> {
    match houses_ex(system, armc_deg, lat_deg, eps_deg, sun_dec_deg) {
        Ok(h) => Ok((h, false)),
        Err(Error::PolarHouses { .. }) => Ok((houses('O', armc_deg, lat_deg, eps_deg)?, true)),
        Err(e) => Err(e),
    }
}

// ---------------------------------------------------------------------------
// House position
// ---------------------------------------------------------------------------

/// Equatorial (RA, Dec) in degrees from ecliptic (lon, lat).
#[inline]
fn ecl_to_eq(lon: f64, lat: f64, eps: f64) -> (f64, f64) {
    let (sl, cl) = (lon * DEG).sin_cos();
    let (sb, cb) = (lat * DEG).sin_cos();
    let (se, ce) = (eps * DEG).sin_cos();
    let x = cb * cl;
    let y = cb * sl * ce - sb * se;
    let z = cb * sl * se + sb * ce;
    (norm_deg(y.atan2(x) / DEG), clamp1(z).asin() / DEG)
}

/// Angle (degrees, (-180, 180]) of the point around the north–south axis of
/// the horizon, measured from the zenith towards the east.
fn psi(armc: f64, lat: f64, ra: f64, dec: f64) -> f64 {
    let h = (ra - armc) * DEG; // east positive
    let (sd, cd) = (dec * DEG).sin_cos();
    let (sh, ch) = h.sin_cos();
    let (sf, cf) = (lat * DEG).sin_cos();
    let zen = cf * cd * ch + sf * sd;
    (cd * sh).atan2(zen) / DEG
}

/// "House degree" (house 10 at 0, house 1 at 90, house 4 at 180, house 7 at
/// 270) to house number in [1, 13).
#[inline]
fn hdeg_to_pos(hdeg: f64) -> f64 {
    norm_deg(hdeg + 270.0) / 30.0 + 1.0
}

/// Placidus "house degree" of a point (0 upper meridian, 90 rising,
/// 180 lower meridian, 270 setting): fraction of its own semi-arcs;
/// circumpolar points use clamped semi-arcs (Otto Ludwig's rule).
fn placidus_hdeg(armc: f64, lat: f64, ra: f64, dec: f64) -> f64 {
    let s = (lat * DEG).tan() * (dec * DEG).tan();
    let ad = clamp1(s).asin() / DEG;
    let dsa = 90.0 + ad;
    let nsa = 90.0 - ad;
    let md = diff_deg(ra, armc); // > 0: east of the upper meridian
    if md >= -dsa && md < dsa {
        90.0 * md / dsa
    } else if md >= dsa {
        90.0 + 90.0 * (md - dsa) / nsa
    } else {
        -90.0 - 90.0 * (-md - dsa) / nsa
    }
}

fn linear_pos(c: &[f64; 12], lon: f64) -> f64 {
    for i in 0..12 {
        let w = norm_deg(c[(i + 1) % 12] - c[i]);
        let x = norm_deg(lon - c[i]);
        if x < w {
            return 1.0 + i as f64 + x / w;
        }
    }
    1.0 // all cusps coincide
}

/// Topocentric house degree H in [0, 180] for a point east of the meridian:
/// solves RA = RAMC + H + asin(tan p(H) · tan δ), tan p(H) = tan φ · min(H, 180-H)/90.
fn topo_hdeg(armc: f64, lat: f64, ra: f64, dec: f64) -> f64 {
    let tf = (lat * DEG).tan();
    let td = (dec * DEG).tan();
    let g = |h: f64| -> f64 {
        let hq = if h <= 90.0 { h } else { 180.0 - h };
        diff_deg(armc + h + clamp1(tf * hq / 90.0 * td).asin() / DEG, ra)
    };
    let (mut lo, mut hi) = (0.0f64, 180.0f64);
    // East of the meridian g(0) lies in [-180, 0]; a point on the IC gives ±180 with
    // the sign left to rounding, which would read as a crossing at h = 0 (the MC).
    let mut glo = -g(lo).abs();
    for _ in 0..64 {
        let m = 0.5 * (lo + hi);
        let gm = g(m);
        if (gm > 0.0) == (glo > 0.0) {
            lo = m;
            glo = gm;
        } else {
            hi = m;
        }
    }
    lo
}

/// House position of an ecliptic point as a continuous value in [1, 13).
///
/// Reproduces `swe.house_pos(armc, lat, eps, (lon, lat_ecl), hsys)`:
/// - P: fraction of the point's own diurnal / nocturnal semi-arc
///   (circumpolar points use clamped semi-arcs); works inside the polar circle.
/// - K: time of rising (east) / setting (west) relative to the trisected
///   semi-arc of the MC degree; returns `Ok(0.0)`, as Swiss does, where the
///   Koch position is undefined (points near the meridian with latitude).
/// - R/C: angle around the north–south horizon axis (C: uniform on the prime
///   vertical; R: mapped to the equator, tan H = cos φ · tan ψ).
/// - T: Polich/Page pole law, solved numerically (Swiss differs by ≲ 1e-5°,
///   and its circumpolar branch is not reproduced — see the report).
/// - E/W/O: linear in longitude between the cusps of [`houses`].
///
/// A point on a cusp, or less than 1 mas before it, gets that cusp's number (Swiss
/// likewise shifts results by 1 mas), so a cusp or angle evaluated at its own
/// longitude returns an integer despite rounding, and `floor` gives its house.
pub fn house_position(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    lon_deg: f64,
    lat_ecl_deg: f64,
) -> Result<f64> {
    house_position_ex(system, armc_deg, lat_deg, eps_deg, lon_deg, lat_ecl_deg, None)
}

/// [`house_position`] for every system, with the Sun's declination for the
/// Sunshine systems ('I', 'i'). 'G' returns a sector in [1, 37), counted
/// clockwise. Definitions per system: docs/GEOMETRY_REPORT.md §1.
pub fn house_position_ex(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    lon_deg: f64,
    lat_ecl_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<f64> {
    let sectors = if canonical_system(system)? == 'G' { 36.0 } else { 12.0 };
    house_position_raw(system, armc_deg, lat_deg, eps_deg, lon_deg, lat_ecl_deg, sun_dec_deg)
        .map(|pos| snap_to_cusp(pos, sectors))
}

/// Moves a position lying less than 1 mas before a cusp onto that cusp, in [1, sectors + 1).
fn snap_to_cusp(pos: f64, sectors: f64) -> f64 {
    let tol = (1.0 / 3_600_000.0) / (360.0 / sectors);
    let next = pos.floor() + 1.0;
    if pos >= 1.0 && next - pos < tol {
        if next >= sectors + 1.0 {
            1.0
        } else {
            next
        }
    } else {
        pos
    }
}

fn house_position_raw(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    lon_deg: f64,
    lat_ecl_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<f64> {
    let sys = canonical_system(system)?;
    if !matches!(sys, 'P' | 'K' | 'O' | 'R' | 'C' | 'T' | 'E' | 'W') {
        return position_ext::house_position_extra(
            sys,
            armc_deg,
            lat_deg,
            eps_deg,
            lon_deg,
            lat_ecl_deg,
            sun_dec_deg,
        );
    }
    let lat = clamp_lat(lat_deg);
    let armc = norm_deg(armc_deg);
    let eps = eps_deg;
    if matches!(sys, 'E' | 'W' | 'O') {
        let a = angles(armc, lat, eps);
        let c = simple_cusps(sys, a.asc, a.mc);
        return Ok(linear_pos(&c, norm_deg(lon_deg)));
    }
    let (ra, dec) = ecl_to_eq(lon_deg, lat_ecl_deg, eps);
    let pos = match sys {
        'C' => hdeg_to_pos(psi(armc, lat, ra, dec)),
        'R' => {
            let p = psi(armc, lat, ra, dec) * DEG;
            let h = (p.sin() * (lat * DEG).cos()).atan2(p.cos()) / DEG;
            hdeg_to_pos(h)
        }
        'P' => hdeg_to_pos(placidus_hdeg(armc, lat, ra, dec)),
        'K' => {
            let d = mc_semi_arc(armc, lat, eps);
            let ad = clamp1((lat * DEG).tan() * (dec * DEG).tan()).asin() / DEG;
            let hdeg = if diff_deg(ra, armc) > 0.0 {
                let x = diff_deg(ra - 90.0 - ad, armc);
                if x.abs() > d + 1e-9 {
                    return Ok(0.0);
                }
                90.0 + 90.0 * x / d
            } else {
                let y = diff_deg(ra + 90.0 + ad, armc);
                if y.abs() > d + 1e-9 {
                    return Ok(0.0);
                }
                -90.0 + 90.0 * y / d
            };
            hdeg_to_pos(hdeg)
        }
        _ => {
            if diff_deg(ra, armc) >= 0.0 {
                hdeg_to_pos(topo_hdeg(armc, lat, ra, dec))
            } else {
                let mirrored = norm_deg(2.0 * armc - ra);
                hdeg_to_pos(-topo_hdeg(armc, lat, mirrored, dec))
            }
        }
    };
    Ok(pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cusps_evaluate_to_their_own_number() {
        for sys in ['P', 'K', 'R', 'C', 'E', 'W', 'O', 'B', 'M', 'T'] {
            for (armc, lat) in [(10.0, 48.85), (200.5, -34.6), (333.3, 51.5), (77.7, 28.6)] {
                let eps = 23.4392911;
                let Ok(h) = houses(sys, armc, lat, eps) else { continue };
                for (i, &cusp) in h.cusps.iter().enumerate() {
                    let pos = house_position(sys, armc, lat, eps, cusp, 0.0).unwrap();
                    assert_eq!(
                        pos.floor() as usize,
                        i + 1,
                        "{sys} armc {armc} lat {lat} cusp {} pos {pos}",
                        i + 1
                    );
                }
            }
        }
    }

    #[test]
    fn known_case_placidus() {
        // swe.houses_armc(100, 45, 23.44, b'P')
        let h = houses('P', 100.0, 45.0, 23.44).unwrap();
        let exp = [
            187.6006086435921,
            213.55031668475442,
            244.43574964100873,
            279.1894652314951,
            313.391891753572,
            343.12758558827136,
        ];
        for i in 0..6 {
            assert!(diff_deg(h.cusps[i], exp[i]).abs() < 1e-7, "{i}");
        }
        assert!(diff_deg(h.vertex, 18.949795709273314).abs() < 1e-9);
    }

    #[test]
    fn polar_rule() {
        assert!(houses('P', 10.0, 66.55, 23.44).is_ok());
        assert!(matches!(houses('K', 10.0, 66.56, 23.44), Err(Error::PolarHouses { .. })));
        let (h, fb) = houses_swiss_compat('P', 10.0, 70.0, 23.44).unwrap();
        assert!(fb);
        assert_eq!(h, houses('O', 10.0, 70.0, 23.44).unwrap());
    }
}
