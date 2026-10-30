// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! The house systems beyond the original eight: A V D N (equal family),
//! B Alcabitius, X Meridian (axial rotation), F Carter poli-equatorial,
//! M Morinus, H Horizon (azimuthal), U Krusinski-Pisa-Goelzer, Y APC,
//! L Pullen SD, Q Pullen SR, S Sripati, I/i Sunshine (variant I / Makransky),
//! G Gauquelin (36 sectors).
//!
//! Each is built from its published definition (sources and the measured
//! Swiss edge rules: docs/GEOMETRY_REPORT.md §1). Edge rules shared with the
//! original systems: latitude clamp, Asc hemisphere rule (the Asc is always in
//! `(MC, MC + 180)`), and for Y and I the Regiomontanus-style MC/IC exchange
//! when that rule flips the Asc.

use super::{
    angles, asc_or, atan2_deg, canonical_system, clamp1, clamp_lat, fill_opposites, houses_ex, is_polar,
    mc_deg, placidus_frac, porphyry, ra_to_lon, sin_cos_deg, Angles, Houses,
};
use crate::math::{cross, norm_deg, Vec3, DEG};
use crate::{Error, Result};

/// Cusps of any system (12, or 36 for 'G') and the 8 Swiss `ascmc` angles.
#[derive(Debug, Clone, PartialEq)]
pub struct HouseCusps {
    pub system: char,
    /// cusps[0] = house (sector) 1 ...; 12 values, or 36 for 'G'
    /// (Gauquelin sectors are numbered clockwise).
    pub cusps: Vec<f64>,
    /// Asc, MC, ARMC, Vertex, equatorial Asc, co-Asc (Koch), co-Asc (Munkasey), polar Asc.
    pub ascmc: [f64; 8],
}

/// Any system, including 'G' (36 sectors). `sun_dec_deg` is required for
/// 'I'/'i' and ignored otherwise.
pub fn house_cusps(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<HouseCusps> {
    let sys = canonical_system(system)?;
    if sys == 'G' {
        let (c, ascmc) = gauquelin_cusps(armc_deg, lat_deg, eps_deg)?;
        return Ok(HouseCusps { system: sys, cusps: c.to_vec(), ascmc });
    }
    let h = houses_ex(sys, armc_deg, lat_deg, eps_deg, sun_dec_deg)?;
    Ok(HouseCusps { system: sys, cusps: h.cusps.to_vec(), ascmc: h.ascmc() })
}

/// Swiss-compatible [`house_cusps`]: where the oracle returns ERR and
/// Porphyry (P, K, G inside the polar circle; 'i' with a circumpolar Sun),
/// returns Porphyry and `true`. For 'G' the fallback fills sectors 1–12 with
/// the 12 Porphyry cusps and sectors 13–36 with 0.0, as the oracle does.
pub fn house_cusps_swiss_compat(
    system: char,
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    sun_dec_deg: Option<f64>,
) -> Result<(HouseCusps, bool)> {
    match house_cusps(system, armc_deg, lat_deg, eps_deg, sun_dec_deg) {
        Ok(h) => Ok((h, false)),
        Err(Error::PolarHouses { system: s, .. }) => {
            let o = super::houses('O', armc_deg, lat_deg, eps_deg)?;
            let mut cusps = o.cusps.to_vec();
            if s == 'G' {
                cusps.resize(36, 0.0);
            }
            Ok((HouseCusps { system: s, cusps, ascmc: o.ascmc() }, true))
        }
        Err(e) => Err(e),
    }
}

/// Gauquelin sectors: 36 cusps, clockwise from the Asc (sector 1), MC = 10,
/// Desc = 19, IC = 28. Each quadrant's semi-arc divided in 9 (Placidus
/// construction with fractions k/9). Undefined (`PolarHouses`) where
/// Placidus is.
pub fn gauquelin_cusps(armc_deg: f64, lat_deg: f64, eps_deg: f64) -> Result<([f64; 36], [f64; 8])> {
    let lat = clamp_lat(lat_deg);
    let armc = norm_deg(armc_deg);
    if is_polar(lat, eps_deg) {
        return Err(Error::PolarHouses { system: 'G', lat_deg });
    }
    let a = angles(armc, lat, eps_deg);
    let mut c = [0.0f64; 36];
    c[0] = a.asc;
    c[9] = a.mc;
    for j in 1..9 {
        c[j] = placidus_frac(armc, lat, eps_deg, (9 - j) as f64 / 9.0, true);
        c[27 + j] = placidus_frac(armc, lat, eps_deg, j as f64 / 9.0, false);
    }
    c[18] = norm_deg(a.asc + 180.0);
    c[27] = norm_deg(a.mc + 180.0);
    for j in 1..9 {
        c[9 + j] = norm_deg(c[27 + j] + 180.0);
        c[18 + j] = norm_deg(c[j] + 180.0);
    }
    Ok((c, ascmc_of(&a, a.mc, armc)))
}

fn ascmc_of(a: &Angles, mc: f64, armc: f64) -> [f64; 8] {
    [a.asc, mc, armc, a.vertex, a.equatorial_asc, a.co_asc_koch, a.co_asc_munkasey, a.polar_asc]
}

/// Dispatch for the 12-cusp systems not handled in the parent module.
pub(super) fn houses_extra(
    sys: char,
    armc_deg: f64,
    lat_deg: f64,
    eps: f64,
    sun_dec: Option<f64>,
) -> Result<Houses> {
    let lat = clamp_lat(lat_deg);
    let armc = norm_deg(armc_deg);
    let a = angles(armc, lat, eps);
    let mut mc = a.mc;
    let c = match sys {
        'A' => equal_from(a.asc),
        'V' => equal_from(norm_deg(a.asc - 15.0)),
        'D' => equal_from(norm_deg(a.mc + 90.0)),
        'N' => equal_from(0.0),
        'X' => std::array::from_fn(|i| ra_to_lon(armc + 90.0 + 30.0 * i as f64, eps)),
        'F' => {
            let ra_asc = lon_to_ra(a.asc, eps);
            std::array::from_fn(|i| ra_to_lon(ra_asc + 30.0 * i as f64, eps))
        }
        'M' => std::array::from_fn(|i| equator_point_lon(armc + 90.0 + 30.0 * i as f64, eps)),
        'B' => alcabitius(armc, lat, eps, &a),
        'H' => horizon(armc, lat, eps),
        'U' => krusinski(armc, lat, eps, a.asc),
        'J' => {
            let mut c = savard_a(armc, lat, eps, a.mc);
            if a.flipped {
                // Campanus rule: the whole system turns with the flipped Asc
                for x in c.iter_mut() {
                    *x = norm_deg(*x + 180.0);
                }
                mc = norm_deg(mc + 180.0);
            }
            c[0] = a.asc;
            c[6] = norm_deg(a.asc + 180.0);
            c
        }
        'Y' => {
            // Built on the raw rising point; when the hemisphere rule flipped
            // the Asc, the raw construction's house 7 is the Asc.
            let raw_asc = if a.flipped { norm_deg(a.asc + 180.0) } else { a.asc };
            let mut c = apc(armc, lat, eps, raw_asc, a.mc);
            if a.flipped {
                // Same handling as the R/T/C flip: the quadrants are
                // re-labelled so that house 10 starts at the IC.
                c.rotate_left(6);
                mc = norm_deg(mc + 180.0);
            }
            c
        }
        'L' | 'Q' | 'S' => {
            let mut p = [0.0; 12];
            porphyry(&mut p, a.asc, a.mc);
            match sys {
                'S' => std::array::from_fn(|i| {
                    let prev = p[(i + 11) % 12];
                    norm_deg(prev + norm_deg(p[i] - prev) / 2.0)
                }),
                _ => pullen(sys, a.asc, a.mc),
            }
        }
        'I' | 'i' => {
            let dec = sun_dec.ok_or_else(|| {
                Error::InvalidArgument("Sunshine houses need the Sun's declination (houses_ex)".into())
            })?;
            let (c, flip_mc) = sunshine(sys, armc, lat, eps, dec, &a, lat_deg)?;
            if flip_mc {
                mc = norm_deg(mc + 180.0);
            }
            c
        }
        _ => {
            return Err(Error::InvalidArgument(format!(
                "house system {sys:?} has no 12-cusp form (use house_cusps)"
            )))
        }
    };
    let ascmc = ascmc_of(&a, mc, armc);
    Ok(Houses {
        cusps: c,
        asc: ascmc[0],
        mc,
        armc,
        vertex: a.vertex,
        equatorial_asc: a.equatorial_asc,
        co_asc_koch: a.co_asc_koch,
        co_asc_munkasey: a.co_asc_munkasey,
        polar_asc: a.polar_asc,
    })
}

fn equal_from(start: f64) -> [f64; 12] {
    std::array::from_fn(|i| norm_deg(start + 30.0 * i as f64))
}

/// Right ascension (deg) of the ecliptic point of longitude `lon`.
#[inline]
pub(super) fn lon_to_ra(lon: f64, eps: f64) -> f64 {
    let (s, c) = sin_cos_deg(lon);
    atan2_deg(s * (eps * DEG).cos(), c)
}

/// Ecliptic longitude of the equator point with right ascension `ra`
/// (projection along a circle through the ecliptic poles: Morinus).
#[inline]
pub(super) fn equator_point_lon(ra: f64, eps: f64) -> f64 {
    let (s, c) = sin_cos_deg(ra);
    atan2_deg(s * (eps * DEG).cos(), c)
}

/// Inverse of [`equator_point_lon`]: RA of the equator point whose
/// ecliptic projection has longitude `lon`.
#[inline]
pub(super) fn equator_point_ra(lon: f64, eps: f64) -> f64 {
    let (s, c) = sin_cos_deg(lon);
    atan2_deg(s, c * (eps * DEG).cos())
}

/// Diurnal semi-arc (deg) of a declination at a latitude, clamped to
/// [0, 180] for circumpolar declinations.
#[inline]
pub(super) fn semi_arc(lat: f64, dec: f64) -> f64 {
    90.0 + clamp1((lat * DEG).tan() * (dec * DEG).tan()).asin() / DEG
}

/// Declination (deg) of the ecliptic point `lon`.
#[inline]
pub(super) fn ecl_dec(lon: f64, eps: f64) -> f64 {
    clamp1((eps * DEG).sin() * (lon * DEG).sin()).asin() / DEG
}

/// Alcabitius: the Asc's diurnal and nocturnal semi-arcs trisected on the
/// equator from the meridian; cusps are the ecliptic points with those RAs
/// (house circles through the celestial poles).
fn alcabitius(armc: f64, lat: f64, eps: f64, a: &Angles) -> [f64; 12] {
    let dsa = semi_arc(lat, ecl_dec(a.asc, eps));
    let nsa = 180.0 - dsa;
    let mut c = [0.0; 12];
    c[0] = a.asc;
    c[9] = a.mc;
    c[10] = ra_to_lon(armc + dsa / 3.0, eps);
    c[11] = ra_to_lon(armc + 2.0 * dsa / 3.0, eps);
    c[1] = ra_to_lon(armc + dsa + nsa / 3.0, eps);
    c[2] = ra_to_lon(armc + dsa + 2.0 * nsa / 3.0, eps);
    fill_opposites(&mut c);
    c
}

#[inline]
pub(super) fn eq_vec(ra: f64, dec: f64) -> Vec3 {
    let (sa, ca) = sin_cos_deg(ra);
    let (sd, cd) = sin_cos_deg(dec);
    [cd * ca, cd * sa, sd]
}

/// Ecliptic longitude (deg) of an equatorial direction.
#[inline]
pub(super) fn vec_lon(v: Vec3, eps: f64) -> f64 {
    let (se, ce) = (eps * DEG).sin_cos();
    atan2_deg(v[1] * ce + v[2] * se, v[0])
}

/// Equatorial unit vector of the ecliptic point `lon`.
#[inline]
pub(super) fn ecl_vec(lon: f64, eps: f64) -> Vec3 {
    let (sl, cl) = sin_cos_deg(lon);
    let (se, ce) = (eps * DEG).sin_cos();
    [cl, sl * ce, sl * se]
}

#[inline]
pub(super) fn vec_ra(v: Vec3) -> f64 {
    atan2_deg(v[1], v[0])
}

/// Zenith, south point and east point of the horizon (equatorial vectors).
pub(super) fn horizon_frame(armc: f64, lat: f64) -> (Vec3, Vec3, Vec3) {
    let (sf, cf) = sin_cos_deg(lat);
    let (st, ct) = sin_cos_deg(armc);
    let zenith = [cf * ct, cf * st, sf];
    let south = [sf * ct, sf * st, -cf];
    let east = [-st, ct, 0.0];
    (zenith, south, east)
}

/// Horizon (azimuthal) houses: the horizon divided into 12 from the south
/// point (cusp 10) through the east point (cusp 1); the house circles are the
/// verticals through those azimuths. The intersection with the ecliptic is
/// taken along `(Z × d) × n_ecl` (d: horizontal direction of the azimuth).
/// For φ ≤ 0 the azimuths are counted from the north point and the result
/// turned by 180°. This orientation is the oracle's (measured); its
/// `house_pos` numbers plain azimuths from the south point in both
/// hemispheres, which is not always consistent with these cusps (report).
fn horizon(armc: f64, lat: f64, eps: f64) -> [f64; 12] {
    let (z, s, e) = horizon_frame(armc, lat);
    let (se, ce) = (eps * DEG).sin_cos();
    let n_ecl = [0.0, -se, ce];
    // southern hemisphere (and the equator itself): azimuths counted from the
    // north point, result turned by 180°
    let south = lat <= 0.0;
    let s = if south { [-s[0], -s[1], -s[2]] } else { s };
    let turn = if south { 180.0 } else { 0.0 };
    let mut c = [0.0; 12];
    for j in 0..12 {
        let (sa, ca) = sin_cos_deg(30.0 * j as f64);
        let d = [ca * s[0] + sa * e[0], ca * s[1] + sa * e[1], ca * s[2] + sa * e[2]];
        let u = cross(cross(z, d), n_ecl);
        c[(9 + j) % 12] = norm_deg(vec_lon(u, eps) + turn);
    }
    c
}

/// Prime-vertical angles (deg from the zenith, eastwards) of the Savard-A
/// cusps 11, 12, 2, 3: the prime-vertical points whose declinations are 2/3,
/// 1/3, −1/3, −2/3 of the zenith's (= φ); a point at altitude h on the prime
/// vertical has sin δ = sin φ sin h. On the equator (φ = 0) the limit
/// sin h = 2/3, 1/3 is used.
pub(super) fn savard_k(lat: f64) -> [f64; 4] {
    let sf = (lat * DEG).sin();
    let alt = |f: f64| {
        if sf == 0.0 {
            f.asin() / DEG
        } else {
            clamp1((f * lat * DEG).sin() / sf).asin() / DEG
        }
    };
    let (a11, a12) = (alt(2.0 / 3.0), alt(1.0 / 3.0));
    [90.0 - a11, 90.0 - a12, 90.0 + a12, 90.0 + a11]
}

/// Savard-A (J. Savard's "Albategnus" houses): parallels of declination
/// divide the prime vertical between the east/west points and the zenith
/// (nadir) in thirds of the zenith's declination; house circles are the great
/// circles through the north and south points of the horizon and those
/// division points (Campanus construction with these prime-vertical points).
fn savard_a(armc: f64, lat: f64, eps: f64, mc: f64) -> [f64; 12] {
    let (sf, cf) = (lat * DEG).sin_cos();
    let fb = norm_deg(mc + 90.0);
    let mut c = [0.0; 12];
    for (&i, k) in [10usize, 11, 1, 2].iter().zip(savard_k(lat)) {
        let (sk, ck) = (k * DEG).sin_cos();
        let pole = clamp1(sf * sk).asin() / DEG;
        let h = (sk * cf).atan2(ck) / DEG;
        c[i] = asc_or(armc + h - 90.0, pole, eps, fb);
    }
    c[9] = mc;
    c[3] = norm_deg(mc + 180.0);
    fill_opposites(&mut c);
    c
}

/// Krusinski-Pisa-Goelzer: the great circle through the Asc and the zenith
/// divided into 12 from the Asc (cusp 1; zenith = cusp 10); each division
/// point projected to the ecliptic along its hour circle (same RA).
fn krusinski(armc: f64, lat: f64, eps: f64, asc: f64) -> [f64; 12] {
    let av = ecl_vec(asc, eps);
    let (z, _, _) = horizon_frame(armc, lat);
    let mut c = [0.0; 12];
    for j in 0..12 {
        let (s, co) = sin_cos_deg(30.0 * j as f64);
        let p = [av[0] * co + z[0] * s, av[1] * co + z[1] * s, av[2] * co + z[2] * s];
        c[(12 - j) % 12] = ra_to_lon(vec_ra(p), eps);
    }
    c
}

/// Equator crossing `k` (deg from the upper meridian, eastwards) of the
/// circle through the north and south points of the horizon that passes
/// through the point at hour-angle offset `d` (RA − ARMC) and declination
/// `dec`: tan k = sin d / (cos d + tan φ tan δ).
#[inline]
pub(super) fn position_circle_k(d: f64, dec: f64, lat: f64) -> f64 {
    let t = (lat * DEG).tan() * (dec * DEG).tan();
    let (sd, cd) = sin_cos_deg(d);
    sd.atan2(cd + t) / DEG
}

/// Ecliptic cusp on the circle through the N/S horizon points crossing the
/// equator at `k` (pole tan p = tan φ sin k, Regiomontanus construction).
#[inline]
fn circle_cusp(armc: f64, lat: f64, eps: f64, k: f64) -> f64 {
    let pole = ((lat * DEG).tan() * (k * DEG).sin()).atan() / DEG;
    asc_or(armc + k - 90.0, pole, eps, norm_deg(mc_deg(armc, eps) + 90.0))
}

/// Hour-angle offsets (RA − ARMC) of the trisection points of the semi-arcs
/// of a parallel with diurnal semi-arc `dsa`, for houses 11, 12, 2, 3, 9, 8, 6, 5.
fn trisection_offsets(dsa: f64) -> [(usize, f64); 8] {
    let nsa = 180.0 - dsa;
    [
        (10, dsa / 3.0),
        (11, 2.0 * dsa / 3.0),
        (1, dsa + nsa / 3.0),
        (2, dsa + 2.0 * nsa / 3.0),
        (8, -dsa / 3.0),
        (7, -2.0 * dsa / 3.0),
        (5, -dsa - nsa / 3.0),
        (4, -dsa - 2.0 * nsa / 3.0),
    ]
}

/// APC (Knegt): the Asc's parallel of declination divided in six above and
/// six below the horizon; cusps on the circles through the N/S horizon
/// points and those division points.
fn apc(armc: f64, lat: f64, eps: f64, asc: f64, mc: f64) -> [f64; 12] {
    let dec = ecl_dec(asc, eps);
    let dsa = semi_arc(lat, dec);
    let mut c = [0.0; 12];
    c[0] = asc;
    c[6] = norm_deg(asc + 180.0);
    c[9] = mc;
    c[3] = norm_deg(mc + 180.0);
    for (i, d) in trisection_offsets(dsa) {
        c[i] = circle_cusp(armc, lat, eps, position_circle_k(d, dec, lat));
    }
    c
}

/// Pullen sinusoidal delta ('L') and sinusoidal ratio ('Q') from the
/// ecliptic quadrant sizes.
fn pullen(sys: char, asc: f64, mc: f64) -> [f64; 12] {
    let q = norm_deg(asc - mc); // MC -> Asc: houses 10, 11, 12
    let (w1, w2) = if sys == 'L' {
        (pullen_sd(q), pullen_sd(180.0 - q))
    } else {
        let small = q.min(180.0 - q);
        let (sm, lg) = pullen_sr(small);
        if q <= 90.0 {
            (sm, lg)
        } else {
            (lg, sm)
        }
    };
    let mut c = [0.0; 12];
    c[9] = mc;
    c[10] = norm_deg(mc + w1[0]);
    c[11] = norm_deg(c[10] + w1[1]);
    c[0] = asc;
    c[1] = norm_deg(asc + w2[0]);
    c[2] = norm_deg(c[1] + w2[1]);
    fill_opposites(&mut c);
    c
}

/// Sinusoidal delta: middle house 30 + (q − 90)/2, outer houses
/// 30 + (q − 90)/4; a negative middle house is set to 0 (agreement recorded
/// in the Swiss documentation, Jan 2016) and the outer ones share the quadrant.
fn pullen_sd(q: f64) -> [f64; 3] {
    let dev = q - 90.0;
    let mid = 30.0 + dev / 2.0;
    if mid < 0.0 {
        [q / 2.0, 0.0, q / 2.0]
    } else {
        let side = 30.0 + dev / 4.0;
        [side, mid, side]
    }
}

/// Sinusoidal ratio (Pullen, 24 Jan 2016): small quadrant `s` = rx + x + rx,
/// large one 180 − s = r³x + r⁴x + r³x. Returns (small, large) house widths.
fn pullen_sr(s: f64) -> ([f64; 3], [f64; 3]) {
    if s <= 0.0 {
        return ([0.0; 3], [54.0, 72.0, 54.0]);
    }
    // f(r) = s r³ (2 + r) − (180 − s)(1 + 2r), increasing for r ≥ 1, f(1) ≤ 0.
    let f = |r: f64| s * r * r * r * (2.0 + r) - (180.0 - s) * (1.0 + 2.0 * r);
    let (mut lo, mut hi) = (1.0f64, 2.0f64);
    while f(hi) < 0.0 {
        hi *= 2.0;
    }
    for _ in 0..200 {
        let m = 0.5 * (lo + hi);
        if f(m) < 0.0 {
            lo = m;
        } else {
            hi = m;
        }
        if hi - lo <= 1e-15 * hi {
            break;
        }
    }
    let r = 0.5 * (lo + hi);
    let x = s / (1.0 + 2.0 * r);
    let r3 = r * r * r;
    ([r * x, x, r * x], [r3 * x, r3 * r * x, r3 * x])
}

/// Sunshine (Makransky): the Sun's diurnal and nocturnal semi-arcs
/// trisected; cusps on the circles through the N/S horizon points and the
/// trisection points. 'I' (variant I): a circumpolar Sun has semi-arcs
/// 180/0 and the cusps of the empty half coincide with the meridian; an MC
/// below the horizon is exchanged with the IC (Regiomontanus rule). 'i'
/// (Makransky's solution as the oracle implements it): undefined for a
/// circumpolar Sun (`PolarHouses`, oracle: Porphyry), and no MC/IC exchange.
/// Returns (cusps, mc exchanged).
fn sunshine(
    sys: char,
    armc: f64,
    lat: f64,
    eps: f64,
    dec: f64,
    a: &Angles,
    lat_in: f64,
) -> Result<([f64; 12], bool)> {
    let t = (lat * DEG).tan() * (dec * DEG).tan();
    if sys == 'i' && t.abs() >= 1.0 {
        return Err(Error::PolarHouses { system: 'i', lat_deg: lat_in });
    }
    let dsa = if t >= 1.0 {
        180.0
    } else if t <= -1.0 {
        0.0
    } else {
        90.0 + t.asin() / DEG
    };
    let mut c = [0.0; 12];
    c[9] = a.mc;
    c[3] = norm_deg(a.mc + 180.0);
    for (i, d) in trisection_offsets(dsa) {
        c[i] = circle_cusp(armc, lat, eps, position_circle_k(d, dec, lat));
    }
    if dsa == 0.0 {
        for i in [7, 8, 10, 11] {
            c[i] = a.mc;
        }
    }
    if dsa == 180.0 {
        for i in [1, 2, 4, 5] {
            c[i] = norm_deg(a.mc + 180.0);
        }
    }
    let flip = sys == 'I' && a.flipped;
    if flip {
        for x in c.iter_mut() {
            *x = norm_deg(*x + 180.0);
        }
    }
    c[0] = a.asc;
    c[6] = norm_deg(a.asc + 180.0);
    Ok((c, flip))
}
