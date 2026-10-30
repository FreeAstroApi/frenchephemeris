// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Precession, nutation, obliquity, Earth rotation.
//!
//! CONTRACT (owned by the fundamental-astronomy workstream):
//! - Angles in radians, times JD TT unless named `jd_ut`.
//! - Matrices rotate vectors: `v_of_date = M · v_gcrs` (math::mat_vec).
//! - `FrameSet::new(jd_tt)` computes everything once so a batch of bodies at
//!   one instant shares the work.
//!
//! # Models
//!
//! * Precession: IAU 2006 (Capitaine et al. 2003, Fukushima–Williams angles, Hilton et al.
//!   2006) or Vondrák, Capitaine & Wallace (2011, A&A 534 A22; corrigendum A&A 541 C1)
//!   long-term precession. Frame bias: IAU 2006 (ERFA `bp06` / `ltpb` conventions).
//! * Nutation: IAU 2000A (MHB2000, 1365 terms) with the IAU 2006 adjustments (`nut06a`),
//!   IAU 2000B (77 terms + fixed planetary offsets, `nut00b`), or the Swiss Ephemeris 2.10
//!   variant identified black-box (the 77 IAU 2000B luni-solar terms evaluated with the full
//!   IERS 2003 Delaunay polynomials, no planetary offsets, no IAU 2006 adjustment).
//! * Mean obliquity: IAU 2006 (`obl06`), the angle between the Vondrák equator and
//!   ecliptic poles (consistent with the Vondrák precession matrix), or the separately
//!   published Vondrák 2011 ε_A series (`vondrak_obliquity`), which is Swiss 2.10's value.
//! * Earth rotation: ERA (IAU 2000), GMST IAU 2006, GAST = ERA − EO with s06 (IAU 2006/2000A,
//!   equivalent to ERFA `gst06a`, which contains the equation-of-equinoxes complementary terms);
//!   optionally (`SiderealTime::SwissLongTerm`) the long-term method documented for Swiss
//!   Ephemeris outside 1850–2050 (Simon et al. 1994 mean Sun on Vondrák 2011 precession).
//!
//! `FrameModel::DEFAULT` = Vondrák 2011 precession + IAU 2006 bias + IAU 2000A/2006 nutation
//! and Vondrák-consistent obliquity: sound over ±several millennia and identical to Swiss
//! 2.10 in precession and bias (≤ 0.001 mas), within ~2 mas of Swiss in nutation. See
//! `docs/TIME_FRAMES_REPORT.md`.

#[rustfmt::skip] // generated
mod tables;
pub mod iau1980; // IAU 1976/1980 frames, JPL Horizons convention (SEFLG_JPLHOR)

use crate::math::{mat_mul, norm_rad, rot_x, rot_z, transpose, Mat3, ARCSEC, DAYS_PER_CENTURY, J2000, TAU};
use tables::{FaTerm, LsTerm};

/// Arcseconds in a full circle.
const TURNAS: f64 = 1_296_000.0;
/// 0.1 microarcsecond in radians.
const U2R: f64 = ARCSEC / 1e7;
/// J2000 mean obliquity (IAU 2006), rad.
pub const EPS0_J2000: f64 = 84381.406 * ARCSEC;

#[inline]
fn centuries(jd_tt: f64) -> f64 {
    (jd_tt - J2000) / DAYS_PER_CENTURY
}

// ---------------------------------------------------------------------------------------------
// Model selection
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Precession {
    /// IAU 2006 (P03) via Fukushima–Williams angles; intended for ± a few centuries.
    Iau2006,
    /// Vondrák, Capitaine & Wallace 2011 long-term model (±200 kyr), with IAU 2006 frame bias.
    Vondrak2011,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nutation {
    /// IAU 2000A with the IAU 2006 adjustments (ERFA `nut06a`).
    Iau2000A,
    /// IAU 2000B (ERFA `nut00b`, linear Delaunay arguments, fixed planetary offsets).
    Iau2000B,
    /// Swiss Ephemeris 2.10 behaviour: the 77 IAU 2000B luni-solar terms with the full
    /// IERS 2003 (quartic) Delaunay arguments, no planetary offset, no IAU 2006 adjustment.
    SwissCompat,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Obliquity {
    /// IAU 2006 polynomial (ERFA `obl06`).
    Iau2006,
    /// Consistent with the precession model: `obl06` for IAU 2006, the angle between the
    /// equator and ecliptic poles for Vondrák 2011.
    Consistent,
    /// Vondrák, Capitaine & Wallace (2011) published long-term series for ε_A (their Eq. 10
    /// with the periodic terms of Table 3), `vondrak_obliquity`. This is the mean obliquity
    /// Swiss Ephemeris 2.10 reports (identical to < 0.05 µas over −13000 … +17000, measured
    /// black-box with `swe.calc(jd, SE_ECL_NUT)`).
    Vondrak2011,
}

/// Sidereal-time (Earth rotation) model used by `FrameSet::gast` / `FrameSet::gmst`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SiderealTime {
    /// IAU 2006/2000A: GAST = ERA(UT1) − EO (ERFA `gst06a`) at all dates.
    Iau2006,
    /// The method documented for Swiss Ephemeris 2.x (swisseph.htm §6.1): IAU 2006 inside
    /// 1850-01-01 … 2050-01-01 UT; outside, the mean Sun of Simon et al. (1994) carried to the
    /// ecliptic of date by Vondrák et al. (2011) precession, plus nutation, with constants
    /// that make the function continuous at the two boundaries. See
    /// `gmst_long_term` and `docs/TIME_FRAMES_REPORT.md` §7.
    SwissLongTerm,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FrameModel {
    pub precession: Precession,
    pub nutation: Nutation,
    pub obliquity: Obliquity,
    pub sidereal: SiderealTime,
}

impl FrameModel {
    /// Default: Vondrák 2011 + IAU 2000A (2006-adjusted) + consistent obliquity.
    pub const DEFAULT: FrameModel = FrameModel {
        precession: Precession::Vondrak2011,
        nutation: Nutation::Iau2000A,
        obliquity: Obliquity::Consistent,
        sidereal: SiderealTime::Iau2006,
    };
    /// Fast variant of the default (IAU 2000B nutation).
    pub const VONDRAK_2000B: FrameModel = FrameModel {
        precession: Precession::Vondrak2011,
        nutation: Nutation::Iau2000B,
        obliquity: Obliquity::Consistent,
        sidereal: SiderealTime::Iau2006,
    };
    /// IAU 2006/2000A, exactly ERFA `pnm06a` / `gst06a`.
    pub const IAU2006_2000A: FrameModel = FrameModel {
        precession: Precession::Iau2006,
        nutation: Nutation::Iau2000A,
        obliquity: Obliquity::Consistent,
        sidereal: SiderealTime::Iau2006,
    };
    /// IAU 2006 precession with IAU 2000B nutation (ERFA `fw2m(pfw06 + nut00b)`).
    pub const IAU2006_2000B: FrameModel = FrameModel {
        precession: Precession::Iau2006,
        nutation: Nutation::Iau2000B,
        obliquity: Obliquity::Consistent,
        sidereal: SiderealTime::Iau2006,
    };
    /// Reproduces Swiss Ephemeris 2.10 frames (measured black-box): Vondrák 2011 + bias,
    /// IAU 2000B luni-solar nutation with full arguments, the Vondrák 2011 published ε_A series, and the
    /// documented long-term sidereal time outside 1850–2050.
    pub const SWISS_COMPAT: FrameModel = FrameModel {
        precession: Precession::Vondrak2011,
        nutation: Nutation::SwissCompat,
        obliquity: Obliquity::Vondrak2011,
        sidereal: SiderealTime::SwissLongTerm,
    };
}

impl Default for FrameModel {
    fn default() -> Self {
        FrameModel::DEFAULT
    }
}

// ---------------------------------------------------------------------------------------------
// Fundamental arguments (IERS Conventions 2003), t in Julian centuries TT
// ---------------------------------------------------------------------------------------------

/// Mean anomaly of the Moon.
pub fn fal03(t: f64) -> f64 {
    ((485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470))))) % TURNAS)
        * ARCSEC
}
/// Mean anomaly of the Sun.
pub fn falp03(t: f64) -> f64 {
    ((1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149))))) % TURNAS)
        * ARCSEC
}
/// Mean argument of latitude of the Moon (L − Ω).
pub fn faf03(t: f64) -> f64 {
    ((335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417))))) % TURNAS)
        * ARCSEC
}
/// Mean elongation of the Moon from the Sun.
pub fn fad03(t: f64) -> f64 {
    ((1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169))))) % TURNAS)
        * ARCSEC
}
/// Mean longitude of the Moon's ascending node.
pub fn faom03(t: f64) -> f64 {
    ((450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939))))) % TURNAS)
        * ARCSEC
}
fn fame03(t: f64) -> f64 {
    (4.402608842 + 2608.7903141574 * t) % TAU
}
fn fave03(t: f64) -> f64 {
    (3.176146697 + 1021.3285546211 * t) % TAU
}
fn fae03(t: f64) -> f64 {
    (1.753470314 + 628.3075849991 * t) % TAU
}
fn fama03(t: f64) -> f64 {
    (6.203480913 + 334.0612426700 * t) % TAU
}
fn faju03(t: f64) -> f64 {
    (0.599546497 + 52.9690962641 * t) % TAU
}
fn fasa03(t: f64) -> f64 {
    (0.874016757 + 21.3299104960 * t) % TAU
}
fn faur03(t: f64) -> f64 {
    (5.481293872 + 7.4781598567 * t) % TAU
}
fn fapa03(t: f64) -> f64 {
    (0.024381750 + 0.00000538691 * t) * t
}

// ---------------------------------------------------------------------------------------------
// Nutation
// ---------------------------------------------------------------------------------------------

#[inline]
fn sum_luni_solar(terms: &[LsTerm], fa: &[f64; 5], t: f64) -> (f64, f64) {
    let (mut dp, mut de) = (0.0, 0.0);
    for &(n, sp, spt, cp, ce, cet, se) in terms.iter().rev() {
        // Arguments are already reduced to [0, 2π); the integer combination stays below
        // ~100 rad, so the per-term fmod of ERFA is unnecessary (effect < 1e-17 rad).
        let arg = n[0] as f64 * fa[0]
            + n[1] as f64 * fa[1]
            + n[2] as f64 * fa[2]
            + n[3] as f64 * fa[3]
            + n[4] as f64 * fa[4];
        let (s, c) = arg.sin_cos();
        dp += (sp + spt * t) * s + cp * c;
        de += (ce + cet * t) * c + se * s;
    }
    (dp * U2R, de * U2R)
}

/// Table of e^{i k a} for k = 0..=K (unit complex numbers as (cos, sin)).
struct Powers<const K: usize>([(f64, f64); K]);

impl<const K: usize> Powers<K> {
    #[inline]
    fn new(a: f64) -> Self {
        let mut p = [(1.0, 0.0); K];
        let (s, c) = a.sin_cos();
        if K > 1 {
            p[1] = (c, s);
        }
        for k in 2..K {
            // e^{ika} = e^{i(k-1)a} e^{ia}, re-seeded from sin_cos every 8 steps to bound
            // rounding growth.
            p[k] = if k % 8 == 0 {
                let (s, c) = (k as f64 * a).sin_cos();
                (c, s)
            } else {
                let (c0, s0) = p[k - 1];
                (c0 * c - s0 * s, s0 * c + c0 * s)
            };
        }
        Powers(p)
    }
    /// e^{i n a} for signed n.
    #[inline]
    fn get(&self, n: i8) -> (f64, f64) {
        if n >= 0 {
            self.0[n as usize]
        } else {
            let (c, s) = self.0[(-n) as usize];
            (c, -s)
        }
    }
}

#[inline]
fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

/// IAU 2000A nutation (MHB2000), no IAU 2006 adjustment (ERFA `nut00a`).
///
/// The 1365 term arguments are integer combinations of 5 (luni-solar) or 13 (planetary)
/// fundamental arguments; their sines/cosines are built from tables of e^{ika} by complex
/// multiplication instead of one `sin_cos` per term (identical to ERFA to ~1e-16 rad).
pub fn nut00a(jd_tt: f64) -> (f64, f64) {
    let t = centuries(jd_tt);
    // Luni-solar: IERS 2003 l, F, Om; MHB2000 l', D.
    let fa = [
        fal03(t),
        ((1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))))
            % TURNAS)
            * ARCSEC,
        faf03(t),
        ((1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))))
            % TURNAS)
            * ARCSEC,
        faom03(t),
    ];
    let pw: [Powers<7>; 5] = std::array::from_fn(|k| Powers::new(fa[k]));
    let (mut dp, mut de) = (0.0, 0.0);
    for &(n, sp, spt, cp, ce, cet, se) in tables::NUT00A_LS.iter().rev() {
        let mut z = (1.0, 0.0);
        for k in 0..5 {
            if n[k] != 0 {
                z = cmul(z, pw[k].get(n[k]));
            }
        }
        let (c, s) = z;
        dp += (sp + spt * t) * s + cp * c;
        de += (ce + cet * t) * c + se * s;
    }
    let (dpls, dels) = (dp * U2R, de * U2R);

    // Planetary: MHB2000 linear Delaunay arguments, IERS 2003 planetary longitudes.
    let a = [
        (2.35555598 + 8328.6914269554 * t) % TAU,
        (1.627905234 + 8433.466158131 * t) % TAU,
        (5.198466741 + 7771.3771468121 * t) % TAU,
        (2.18243920 - 33.757045 * t) % TAU,
        fame03(t),
        fave03(t),
        fae03(t),
        fama03(t),
        faju03(t),
        fasa03(t),
        faur03(t),
        (5.321159000 + 3.8127774000 * t) % TAU,
        fapa03(t),
    ];
    let pw: [Powers<22>; 13] = std::array::from_fn(|k| Powers::new(a[k]));
    let (mut dp, mut de) = (0.0, 0.0);
    for &(n, sp, cp, se, ce) in tables::NUT00A_PL.iter().rev() {
        let mut z = (1.0, 0.0);
        for k in 0..13 {
            if n[k] != 0 {
                z = cmul(z, pw[k].get(n[k]));
            }
        }
        let (c, s) = z;
        dp += sp * s + cp * c;
        de += se * s + ce * c;
    }
    (dpls + dp * U2R, dels + de * U2R)
}

/// IAU 2000A nutation with the IAU 2006 adjustments (ERFA `nut06a`).
pub fn nut06a(jd_tt: f64) -> (f64, f64) {
    let t = centuries(jd_tt);
    let fj2 = -2.7774e-6 * t;
    let (dp, de) = nut00a(jd_tt);
    (dp + dp * (0.4697e-6 + fj2), de + de * fj2)
}

/// IAU 2000B nutation (ERFA `nut00b`).
pub fn nut00b(jd_tt: f64) -> (f64, f64) {
    let t = centuries(jd_tt);
    let fa = [
        ((485868.249036 + 1717915923.2178 * t) % TURNAS) * ARCSEC,
        ((1287104.79305 + 129596581.0481 * t) % TURNAS) * ARCSEC,
        ((335779.526232 + 1739527262.8478 * t) % TURNAS) * ARCSEC,
        ((1072260.70369 + 1602961601.2090 * t) % TURNAS) * ARCSEC,
        ((450160.398036 - 6962890.5431 * t) % TURNAS) * ARCSEC,
    ];
    let (dp, de) = sum_luni_solar(&tables::NUT00B_LS, &fa, t);
    const DPPLAN: f64 = -0.135e-3 * ARCSEC;
    const DEPLAN: f64 = 0.388e-3 * ARCSEC;
    (dp + DPPLAN, de + DEPLAN)
}

/// Swiss Ephemeris 2.10 nutation (identified black-box to 3e-8 arcsec): IAU 2000B
/// luni-solar terms with full IERS 2003 Delaunay arguments, without planetary offsets.
pub fn nut_swiss_compat(jd_tt: f64) -> (f64, f64) {
    let t = centuries(jd_tt);
    let fa = [fal03(t), falp03(t), faf03(t), fad03(t), faom03(t)];
    sum_luni_solar(&tables::NUT00B_LS, &fa, t)
}

/// Nutation (dpsi, deps) in rad with an explicit model.
pub fn nutation_with(jd_tt: f64, model: Nutation) -> (f64, f64) {
    match model {
        Nutation::Iau2000A => nut06a(jd_tt),
        Nutation::Iau2000B => nut00b(jd_tt),
        Nutation::SwissCompat => nut_swiss_compat(jd_tt),
    }
}

/// Nutation (dpsi, deps), rad — default model (IAU 2000A with IAU 2006 adjustments).
pub fn nutation(jd_tt: f64) -> (f64, f64) {
    nutation_with(jd_tt, FrameModel::DEFAULT.nutation)
}

/// Nutation matrix N = R1(−(ε+Δε)) · R3(−Δψ) · R1(ε): mean → true equator & equinox of date.
pub fn nutation_matrix(eps_mean: f64, dpsi: f64, deps: f64) -> Mat3 {
    mat_mul(&rot_x(-(eps_mean + deps)), &mat_mul(&rot_z(-dpsi), &rot_x(eps_mean)))
}

// ---------------------------------------------------------------------------------------------
// Precession, bias, obliquity
// ---------------------------------------------------------------------------------------------

/// IAU 2006 mean obliquity (ERFA `obl06`), rad.
pub fn obl06(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    (84381.406
        + t * (-46.836769 + t * (-0.0001831 + t * (0.00200340 + t * (-0.000000576 + t * (-0.0000000434))))))
        * ARCSEC
}

/// IAU 2006 Fukushima–Williams bias-precession angles (γ̄, φ̄, ψ̄), rad (ERFA `pfw06` without ε_A).
pub fn pfw06(jd_tt: f64) -> (f64, f64, f64) {
    let t = centuries(jd_tt);
    let gamb = (-0.052928
        + t * (10.556378 + t * (0.4932044 + t * (-0.00031238 + t * (-0.000002788 + t * 0.0000000260)))))
        * ARCSEC;
    let phib = (84381.412819
        + t * (-46.811016 + t * (0.0511268 + t * (0.00053289 + t * (-0.000000440 + t * (-0.0000000176))))))
        * ARCSEC;
    let psib = (-0.041775
        + t * (5038.481484 + t * (1.5584175 + t * (-0.00018522 + t * (-0.000026452 + t * (-0.0000000148))))))
        * ARCSEC;
    (gamb, phib, psib)
}

/// Fukushima–Williams angles to matrix: R1(−ε) · R3(−ψ) · R1(φ) · R3(γ) (ERFA `fw2m`).
pub fn fw2m(gamb: f64, phib: f64, psi: f64, eps: f64) -> Mat3 {
    mat_mul(&rot_x(-eps), &mat_mul(&rot_z(-psi), &mat_mul(&rot_x(phib), &rot_z(gamb))))
}

/// IAU 2006 frame bias matrix GCRS → mean J2000 (ERFA `bp06` rb).
pub fn frame_bias() -> Mat3 {
    let (g, p, s) = pfw06(J2000);
    fw2m(g, p, s, obl06(J2000))
}

/// Julian epoch (TT) for the Vondrák series.
#[inline]
fn julian_epoch(jd_tt: f64) -> f64 {
    2000.0 + (jd_tt - J2000) / 365.25
}

/// Vondrák 2011 long-term ecliptic pole (unit vector in the J2000 mean equatorial frame).
pub fn ltp_ecliptic_pole(jd_tt: f64) -> [f64; 3] {
    const POL: [[f64; 4]; 2] = [
        [5851.607687, -0.1189000, -0.00028913, 0.000000101],
        [-1600.886300, 1.1689818, -0.00000020, -0.000000437],
    ];
    const PER: [[f64; 5]; 8] = [
        [708.15, -5486.751211, -684.661560, 667.666730, -5523.863691],
        [2309.00, -17.127623, 2446.283880, -2354.886252, -549.747450],
        [1620.00, -617.517403, 399.671049, -428.152441, -310.998056],
        [492.20, 413.442940, -356.652376, 376.202861, 421.535876],
        [1183.00, 78.614193, -186.387003, 184.778874, -36.776172],
        [622.00, -180.732815, -316.800070, 335.321713, -145.278396],
        [882.00, -87.676083, 198.296701, -185.138669, -34.744450],
        [547.00, 46.140315, 101.135679, -120.972830, 22.885731],
    ];
    let t = (julian_epoch(jd_tt) - 2000.0) / 100.0;
    let (mut p, mut q) = (0.0, 0.0);
    let w = TAU * t;
    for r in PER.iter() {
        let (s, c) = (w / r[0]).sin_cos();
        p += c * r[1] + s * r[3];
        q += c * r[2] + s * r[4];
    }
    let mut w = 1.0;
    for i in 0..4 {
        p += POL[0][i] * w;
        q += POL[1][i] * w;
        w *= t;
    }
    p *= ARCSEC;
    q *= ARCSEC;
    let z = (1.0 - p * p - q * q).max(0.0).sqrt();
    let (s, c) = EPS0_J2000.sin_cos();
    [p, -q * c - z * s, -q * s + z * c]
}

/// Vondrák 2011 long-term equator pole (unit vector in the J2000 mean equatorial frame).
pub fn ltp_equator_pole(jd_tt: f64) -> [f64; 3] {
    const POL: [[f64; 4]; 2] = [
        [5453.282155, 0.4252841, -0.00037173, -0.000000152],
        [-73750.930350, -0.7675452, -0.00018725, 0.000000231],
    ];
    const PER: [[f64; 5]; 14] = [
        [256.75, -819.940624, 75004.344875, 81491.287984, 1558.515853],
        [708.15, -8444.676815, 624.033993, 787.163481, 7774.939698],
        [274.20, 2600.009459, 1251.136893, 1251.296102, -2219.534038],
        [241.45, 2755.175630, -1102.212834, -1257.950837, -2523.969396],
        [2309.00, -167.659835, -2660.664980, -2966.799730, 247.850422],
        [492.20, 871.855056, 699.291817, 639.744522, -846.485643],
        [396.10, 44.769698, 153.167220, 131.600209, -1393.124055],
        [288.90, -512.313065, -950.865637, -445.040117, 368.526116],
        [231.10, -819.415595, 499.754645, 584.522874, 749.045012],
        [1610.00, -538.071099, -145.188210, -89.756563, 444.704518],
        [620.00, -189.793622, 558.116553, 524.429630, 235.934465],
        [157.87, -402.922932, -23.923029, -13.549067, 374.049623],
        [220.30, 179.516345, -165.405086, -210.157124, -171.330180],
        [1200.00, -9.814756, 9.344131, -44.919798, -22.899655],
    ];
    let t = (julian_epoch(jd_tt) - 2000.0) / 100.0;
    let (mut x, mut y) = (0.0, 0.0);
    let w = TAU * t;
    for r in PER.iter() {
        let (s, c) = (w / r[0]).sin_cos();
        x += c * r[1] + s * r[3];
        y += c * r[2] + s * r[4];
    }
    let mut w = 1.0;
    for i in 0..4 {
        x += POL[0][i] * w;
        y += POL[1][i] * w;
        w *= t;
    }
    x *= ARCSEC;
    y *= ARCSEC;
    [x, y, (1.0 - x * x - y * y).max(0.0).sqrt()]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    crate::math::cross(a, b)
}

/// Vondrák 2011 precession matrix, mean J2000 → mean of date (ERFA `ltp`), and the
/// consistent mean obliquity of date (angle between the poles).
pub fn ltp(jd_tt: f64) -> (Mat3, f64) {
    let peq = ltp_equator_pole(jd_tt);
    let pec = ltp_ecliptic_pole(jd_tt);
    let v = cross(peq, pec);
    let w = crate::math::norm(v);
    let eqx = [v[0] / w, v[1] / w, v[2] / w];
    let v = cross(peq, eqx);
    let eps = w.atan2(crate::math::dot(peq, pec));
    ([eqx, v, peq], eps)
}

/// Vondrák 2011 precession including the IAU 2006 frame bias (first order, as ERFA `ltpb`):
/// GCRS → mean equator & equinox of date.
pub fn ltpb(jd_tt: f64) -> Mat3 {
    const DX: f64 = -0.016617 * ARCSEC;
    const DE: f64 = -0.0068192 * ARCSEC;
    const DR: f64 = -0.0146 * ARCSEC;
    let (rp, _) = ltp(jd_tt);
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        r[i][0] = rp[i][0] - rp[i][1] * DR + rp[i][2] * DX;
        r[i][1] = rp[i][0] * DR + rp[i][1] + rp[i][2] * DE;
        r[i][2] = -rp[i][0] * DX - rp[i][1] * DE + rp[i][2];
    }
    r
}

/// Vondrák, Capitaine & Wallace (2011, A&A 534 A22) long-term mean obliquity of date ε_A (rad),
/// Eq. (10): ε_A = 84028.206305″ + 0.3624445″T − 0.00004039″T² − 110e-9″T³ + Σ_ε, with the
/// ten cosine/sine terms of their Table 3 (C cos 2πT/P + S sin 2πT/P, P in Julian centuries),
/// T = Julian centuries TT from J2000. Valid ±200 000 years; equals 84381.406″ at J2000 and
/// stays within 0.02″ of IAU 2006 `obl06` over 0 … 3000 (0.29″ at −1000). Unaffected by the
/// 2012 corrigendum (A&A 541 C1), which concerns Table 1 only.
/// Source text: https://hal.science/hal-03785916 (open-access copy of the A&A article).
pub fn vondrak_obliquity(jd_tt: f64) -> f64 {
    // (P [cy], C_eps ["], S_eps ["]) — Table 3, ε_A columns.
    const PER: [[f64; 3]; 10] = [
        [409.90, 753.872780, -1704.720302],
        [396.15, -247.805823, -862.308358],
        [537.22, 379.471484, 447.832178],
        [402.90, -53.880558, -889.571909],
        [417.15, -90.109153, 190.402846],
        [288.92, -353.600190, -56.564991],
        [4043.00, -63.115353, -296.222622],
        [306.00, -28.248187, -75.859952],
        [277.00, 17.703387, 67.473503],
        [203.00, 38.911307, 3.014055],
    ];
    let t = centuries(jd_tt);
    let mut e = 84028.206305 + t * (0.3624445 + t * (-0.00004039 + t * (-110e-9)));
    let w = TAU * t;
    for r in PER.iter() {
        let (s, c) = (w / r[0]).sin_cos();
        e += r[1] * c + r[2] * s;
    }
    e * ARCSEC
}

/// Mean obliquity of date for a given model (rad).
pub fn obliquity_mean_with(jd_tt: f64, precession: Precession, obliquity: Obliquity) -> f64 {
    match (precession, obliquity) {
        (_, Obliquity::Vondrak2011) => vondrak_obliquity(jd_tt),
        (_, Obliquity::Iau2006) | (Precession::Iau2006, _) => obl06(jd_tt),
        (Precession::Vondrak2011, Obliquity::Consistent) => ltp(jd_tt).1,
    }
}

/// Mean obliquity of the ecliptic of date (rad) — default model (Vondrák 2011, from the poles).
pub fn obliquity_mean(jd_tt: f64) -> f64 {
    obliquity_mean_with(jd_tt, FrameModel::DEFAULT.precession, FrameModel::DEFAULT.obliquity)
}

/// GCRS → mean equator & equinox of date (bias × precession) for a precession model.
pub fn bias_precession_matrix(jd_tt: f64, precession: Precession) -> Mat3 {
    match precession {
        Precession::Iau2006 => {
            let (g, p, s) = pfw06(jd_tt);
            fw2m(g, p, s, obl06(jd_tt))
        }
        Precession::Vondrak2011 => ltpb(jd_tt),
    }
}

/// GCRS -> mean ecliptic & equinox of J2000 (frame bias + rotation by J2000 obliquity).
pub fn gcrs_to_ecliptic_j2000() -> Mat3 {
    mat_mul(&rot_x(EPS0_J2000), &frame_bias())
}

// ---------------------------------------------------------------------------------------------
// Earth rotation, CIO locator, sidereal time
// ---------------------------------------------------------------------------------------------

/// Earth rotation angle (IAU 2000), rad, for JD(UT1).
pub fn era(jd_ut: f64) -> f64 {
    let t = jd_ut - J2000;
    let f = jd_ut % 1.0;
    norm_rad(TAU * (f + 0.7790572732640 + 0.00273781191135448 * t))
}

/// IAU 2006 GMST − ERA polynomial (rad).
fn gmst06_poly(t: f64) -> f64 {
    (0.014506
        + t * (4612.156534 + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * (-0.0000000368))))))
        * ARCSEC
}

#[inline]
fn fa8(t: f64) -> [f64; 8] {
    [fal03(t), falp03(t), faf03(t), fad03(t), faom03(t), fave03(t), fae03(t), fapa03(t)]
}

#[inline]
fn sum_fa(terms: &[FaTerm], fa: &[f64; 8]) -> f64 {
    let mut w = 0.0;
    for &(n, s, c) in terms.iter().rev() {
        let mut a = 0.0;
        for k in 0..8 {
            a += n[k] as f64 * fa[k];
        }
        let (sa, ca) = a.sin_cos();
        w += s * sa + c * ca;
    }
    w
}

/// CIO locator s (IAU 2006/2000A, ERFA `s06`) given the CIP X, Y.
pub fn s06(jd_tt: f64, x: f64, y: f64) -> f64 {
    let t = centuries(jd_tt);
    let fa = fa8(t);
    let p = &tables::S06_POLY;
    let w0 = p[0] + sum_fa(&tables::S06_T0, &fa);
    let w1 = p[1] + sum_fa(&tables::S06_T1, &fa);
    let w2 = p[2] + sum_fa(&tables::S06_T2, &fa);
    let w3 = p[3] + sum_fa(&tables::S06_T3, &fa);
    let w4 = p[4] + sum_fa(&tables::S06_T4, &fa);
    let w5 = p[5];
    (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) * ARCSEC - x * y / 2.0
}

/// Equation of the origins (ERFA `eors`) from the NPB matrix and s.
pub fn eors(rnpb: &Mat3, s: f64) -> f64 {
    let x = rnpb[2][0];
    let ax = x / (1.0 + rnpb[2][2]);
    let xs = 1.0 - ax * x;
    let ys = -ax * rnpb[2][1];
    let zs = -x;
    let p = rnpb[0][0] * xs + rnpb[0][1] * ys + rnpb[0][2] * zs;
    let q = rnpb[1][0] * xs + rnpb[1][1] * ys + rnpb[1][2] * zs;
    if p != 0.0 || q != 0.0 {
        s - q.atan2(p)
    } else {
        s
    }
}

/// Equation-of-the-equinoxes complementary terms (IERS 2003, ERFA `eect00`), rad.
pub fn eect00(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    let fa = fa8(t);
    (sum_fa(&tables::EECT_T0, &fa) + sum_fa(&tables::EECT_T1, &fa) * t) * ARCSEC
}

/// Greenwich mean sidereal time (rad), IAU 2006 (ERFA `gmst06`).
pub fn gmst(jd_ut: f64, jd_tt: f64) -> f64 {
    norm_rad(era(jd_ut) + gmst06_poly(centuries(jd_tt)))
}

/// Greenwich apparent sidereal time (rad), IAU 2006/2000A (ERFA `gst06a`: ERA − EO, which
/// includes the equation of the equinoxes with its complementary terms).
pub fn gast(jd_ut: f64, jd_tt: f64) -> f64 {
    FrameSet::with_model(jd_tt, FrameModel::IAU2006_2000A).gast(jd_ut)
}

/// Equation of the equinoxes (rad), IAU 2006/2000A (ERFA `ee06a`).
pub fn equation_of_equinoxes(jd_tt: f64) -> f64 {
    FrameSet::with_model(jd_tt, FrameModel::IAU2006_2000A).eqeq
}

// ---------------------------------------------------------------------------------------------
// Long-term sidereal time (the method documented for Swiss Ephemeris 2.x)
// ---------------------------------------------------------------------------------------------

/// 1850-01-01 0h and 2050-01-01 0h (JD UT): the documented validity limits of the IAU 2006
/// sidereal time in Swiss Ephemeris 2.x (swisseph.htm §6.1).
pub const LONG_TERM_ST_LIMITS: (f64, f64) = (2_396_758.5, 2_469_807.5);

/// 0 before 1850, 1 from 2050 on, `None` inside the IAU interval.
#[inline]
fn long_term_side(jd_ut: f64) -> Option<usize> {
    if jd_ut < LONG_TERM_ST_LIMITS.0 {
        Some(0)
    } else if jd_ut >= LONG_TERM_ST_LIMITS.1 {
        Some(1)
    } else {
        None
    }
}

/// Mean longitude of the Earth–Moon barycentre referred to the fixed mean ecliptic and
/// equinox J2000 (rad), Simon et al. (1994, A&A 282, 663) §3.5, with t in Julian millennia TDB
/// (TT here; the difference is < 2 ms): 100°27′59.22059″ + 1295977422.83429″t − 2.04411″t² −
/// 0.00523″t³. Its linear part is the IERS Conventions (2010) Table 5.2e argument λ_E.
pub fn emb_mean_longitude_j2000(jd_tt: f64) -> f64 {
    let t = (jd_tt - J2000) / 365_250.0;
    (361_679.220_59 + t * (1_295_977_422.834_29 + t * (-2.04411 + t * -0.00523))) * ARCSEC
}

/// Uncorrected long-term Greenwich mean sidereal time (rad): the longitude of the mean Sun
/// (Simon 1994 EMB mean longitude + 180°, a point on the J2000 ecliptic) on the mean ecliptic
/// of date — rotated to the mean equator of date with the Vondrák 2011 precession matrix
/// (`ltp`) and to the ecliptic of date with the Vondrák 2011 published mean obliquity
/// (`vondrak_obliquity`; with `obl06` instead the residual vs Swiss drifts to 23″ at ±13 kyr,
/// with ε_A it is constant to 0.02″ over −13000 … +17000) — plus the Earth's
/// rotation in mean solar days, 2π·frac(UT − J2000). This is "the long term precession model
/// Vondrak 2011 … and the mean motion of the Earth according to Simon & alii 1994" of the
/// Swiss Ephemeris documentation (swisseph.htm §6.1). The mean Sun is taken at TT, so the value
/// moves with ΔT at the solar mean motion (0.0411″ per second of ΔT).
/// `FrameSet::gmst` adds the continuity constant (`long_term_offsets`).
pub fn gmst_long_term(jd_ut: f64, jd_tt: f64) -> f64 {
    let l = emb_mean_longitude_j2000(jd_tt) + std::f64::consts::PI;
    let (sl, cl) = l.sin_cos();
    // J2000 ecliptic -> J2000 mean equator.
    let (se, ce) = EPS0_J2000.sin_cos();
    let v = [cl, sl * ce, sl * se];
    let (p, _) = ltp(jd_tt);
    let w = crate::math::mat_vec(&p, v);
    // Mean equator of date -> mean ecliptic of date.
    let (so, co) = vondrak_obliquity(jd_tt).sin_cos();
    let y = w[1] * co + w[2] * so;
    let lon = y.atan2(w[0]);
    norm_rad(lon + TAU * (jd_ut - J2000).rem_euclid(1.0))
}

/// Constants (rad) added to `gmst_long_term` before 1850 and from 2050 on so that the long-term
/// GMST equals the IAU 2006 GMST (`gmst`) at 1850-01-01 0h and 2050-01-01 0h UT. Computed here
/// from the published formulas, with TT − UT at the boundaries from
/// `DeltaTModel::SWISS_COMPAT` (published SMH 2016 / USNO data). The Swiss constants themselves
/// are not published and are not used.
/// Delta T used inside the long-term sidereal time: SMH 2016 corrected to the DE431 tidal
/// acceleration (−25.80″/cy², published with DE431). Measured black-box, Swiss's sidereal
/// time and houses use this Delta T at every date, even where its planet positions use the
/// DE441 value (−25.936); following it keeps house cusps as Swiss users see them.
pub const LONG_TERM_ST_DELTA_T: crate::time::DeltaTModel = crate::time::DeltaTModel {
    spline: crate::time::DeltaTSpline::Smh2016,
    lunar_ndot: -25.80,
    user_seconds: None,
};

pub fn long_term_offsets() -> [f64; 2] {
    static OFF: std::sync::OnceLock<[f64; 2]> = std::sync::OnceLock::new();
    *OFF.get_or_init(|| {
        let at = |jd_ut: f64| {
            let jd_tt = crate::time::ut_to_tt_with(jd_ut, &LONG_TERM_ST_DELTA_T);
            let d = gmst(jd_ut, jd_tt) - gmst_long_term(jd_ut, jd_tt);
            (d + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI
        };
        [at(LONG_TERM_ST_LIMITS.0), at(LONG_TERM_ST_LIMITS.1)]
    })
}

// ---------------------------------------------------------------------------------------------
// FrameSet
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct FrameSet {
    pub jd_tt: f64,
    /// Mean obliquity of date.
    pub eps_mean: f64,
    /// Nutation in longitude and obliquity.
    pub dpsi: f64,
    pub deps: f64,
    /// True obliquity = eps_mean + deps.
    pub eps_true: f64,
    /// GCRS -> mean equator & equinox of date (frame bias × precession).
    pub bp: Mat3,
    /// GCRS -> true equator & equinox of date (bias × precession × nutation).
    pub bpn: Mat3,
    /// Models used.
    pub model: FrameModel,
    /// GCRS -> mean ecliptic & equinox of date: R1(eps_mean) · bp.
    pub ecl_mean: Mat3,
    /// GCRS -> true ecliptic & equinox of date: R1(eps_true) · bpn.
    pub ecl_true: Mat3,
    /// Equation of the origins (rad): GAST = ERA − eo.
    pub eo: f64,
    /// Equation of the equinoxes (rad): GAST − GMST.
    pub eqeq: f64,
}

impl FrameSet {
    /// Frames at `jd_tt` with the default model (`FrameModel::DEFAULT`).
    pub fn new(jd_tt: f64) -> FrameSet {
        FrameSet::with_model(jd_tt, FrameModel::DEFAULT)
    }

    /// Frames at `jd_tt` with an explicit model.
    pub fn with_model(jd_tt: f64, model: FrameModel) -> FrameSet {
        let t = centuries(jd_tt);
        let (dpsi, deps) = nutation_with(jd_tt, model.nutation);
        let (bp, bpn, eps_mean) = match model.precession {
            Precession::Iau2006 => {
                let (g, p, s) = pfw06(jd_tt);
                let eps = obl06(jd_tt);
                (fw2m(g, p, s, eps), fw2m(g, p, s + dpsi, eps + deps), eps)
            }
            Precession::Vondrak2011 => {
                let bp = ltpb(jd_tt);
                let eps = obliquity_mean_with(jd_tt, model.precession, model.obliquity);
                (bp, mat_mul(&nutation_matrix(eps, dpsi, deps), &bp), eps)
            }
        };
        let eps_true = eps_mean + deps;
        let eo = eors(&bpn, s06(jd_tt, bpn[2][0], bpn[2][1]));
        let eqeq = match model.precession {
            // gst06a − gmst06, independent of UT.
            Precession::Iau2006 => {
                let e = -eo - gmst06_poly(t);
                (e + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI
            }
            Precession::Vondrak2011 => dpsi * eps_mean.cos() + eect00(jd_tt),
        };
        FrameSet {
            jd_tt,
            eps_mean,
            dpsi,
            deps,
            eps_true,
            bp,
            bpn,
            model,
            ecl_mean: mat_mul(&rot_x(eps_mean), &bp),
            ecl_true: mat_mul(&rot_x(eps_true), &bpn),
            eo,
            eqeq,
        }
    }

    /// Nutation matrix (mean of date → true of date).
    pub fn nutation_matrix(&self) -> Mat3 {
        mat_mul(&self.bpn, &transpose(&self.bp))
    }

    /// Greenwich apparent sidereal time (rad): consistent with this frame's true equinox
    /// (`SiderealTime::Iau2006`), or the Swiss-documented long-term model outside 1850–2050.
    pub fn gast(&self, jd_ut: f64) -> f64 {
        if self.model.sidereal == SiderealTime::SwissLongTerm && long_term_side(jd_ut).is_some() {
            return norm_rad(self.gmst(jd_ut) + self.dpsi * self.eps_true.cos());
        }
        norm_rad(era(jd_ut) - self.eo)
    }

    /// Greenwich mean sidereal time (rad), same model choice as `gast`.
    pub fn gmst(&self, jd_ut: f64) -> f64 {
        if self.model.sidereal == SiderealTime::SwissLongTerm {
            if let Some(side) = long_term_side(jd_ut) {
                let jd_tt = crate::time::ut_to_tt_with(jd_ut, &LONG_TERM_ST_DELTA_T);
                return norm_rad(gmst_long_term(jd_ut, jd_tt) + long_term_offsets()[side]);
            }
        }
        norm_rad(era(jd_ut) - self.eo - self.eqeq)
    }

    /// Time derivative of `bpn` (per day), for transforming velocities:
    /// `v_date = bpn·v + bpn_rate·r`. Precession is differenced exactly; the nutation part
    /// uses the 77-term luni-solar series (rate error < 1e-9 rad/day vs IAU 2000A).
    pub fn bpn_rate(&self) -> Mat3 {
        bpn_rate(self.jd_tt, self.model)
    }

    /// Time derivative of `ecl_true` (per day).
    pub fn ecl_true_rate(&self) -> Mat3 {
        let fast = FrameModel { nutation: Nutation::SwissCompat, ..self.model };
        mat_rate(self.jd_tt, |j| FrameSet::cheap(j, fast).1)
    }
}

impl FrameSet {
    /// (bpn, ecl_true) with no sidereal-time work; used for numerical rates.
    fn cheap(jd_tt: f64, model: FrameModel) -> (Mat3, Mat3) {
        let (dpsi, deps) = nutation_with(jd_tt, model.nutation);
        let (bpn, eps) = match model.precession {
            Precession::Iau2006 => {
                let (g, p, s) = pfw06(jd_tt);
                let eps = obl06(jd_tt);
                (fw2m(g, p, s + dpsi, eps + deps), eps)
            }
            Precession::Vondrak2011 => {
                let eps = obliquity_mean_with(jd_tt, model.precession, model.obliquity);
                (mat_mul(&nutation_matrix(eps, dpsi, deps), &ltpb(jd_tt)), eps)
            }
        };
        (bpn, mat_mul(&rot_x(eps + deps), &bpn))
    }
}

fn mat_rate(jd_tt: f64, f: impl Fn(f64) -> Mat3) -> Mat3 {
    const H: f64 = 0.05; // days
    let a = f(jd_tt - H);
    let b = f(jd_tt + H);
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (b[i][j] - a[i][j]) / (2.0 * H);
        }
    }
    r
}

/// d(bpn)/dt (per day) for a model; see `FrameSet::bpn_rate`.
pub fn bpn_rate(jd_tt: f64, model: FrameModel) -> Mat3 {
    let fast = FrameModel { nutation: Nutation::SwissCompat, ..model };
    mat_rate(jd_tt, |j| FrameSet::cheap(j, fast).0)
}

/// Identifies the precession/nutation models in diagnostics.
pub const FRAME_MODEL: &str = "FE-FR/1: precession Vondrak-Capitaine-Wallace 2011 + IAU2006 frame bias; \
     nutation IAU2000A (MHB2000) with IAU2006 adjustments; obliquity from Vondrak poles; \
     GMST/GAST IAU2006 (ERA-based, s06)";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_matches_full_difference() {
        let jd = 2460000.3;
        let r = bpn_rate(jd, FrameModel::DEFAULT);
        let a = FrameSet::new(jd - 0.01).bpn;
        let b = FrameSet::new(jd + 0.01).bpn;
        for i in 0..3 {
            for j in 0..3 {
                let d = (b[i][j] - a[i][j]) / 0.02 - r[i][j];
                assert!(d.abs() < 1e-9, "{i}{j} {d:e}");
            }
        }
    }
}
