// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Calendars, time scales and Delta T.
//!
//! CONTRACT (owned by the fundamental-astronomy workstream):
//! - All functions are pure.
//! - Julian dates are f64 days. Astronomical year numbering (year 0 = 1 BC).
//!
//! # Delta T model (`DELTA_T_MODEL`)
//!
//! `delta_t_seconds` = TT − UT1 from, by date (UT):
//!
//! | range | source |
//! |---|---|
//! | before −720 | long-term parabola of Stephenson, Morrison & Hohenkerk (2016), −320 + 32.5 u², u = (Y−1825)/100, offset to join the spline continuously at −720 |
//! | −720 … 1950.0 | cubic spline of Morrison, Stephenson, Hohenkerk & Zawilski (2021), "Addendum 2020" (Table S15.2020), doi:10.1098/rspa.2020.0776 |
//! | 1950.0 … 1955.0 | the same spline plus a linear ramp of the (observed − spline) offset at 1955.0, so the curve is continuous |
//! | 1955.0 … 1973.0 | USNO semi-annual observed values (`historic_deltat.data`) |
//! | 1973-02-01 … last month | USNO monthly observed values (`deltat.data`), linearly interpolated |
//! | after the last value | ΔT_last + v·τ + 32.5·τ² s (τ in Julian centuries after the last value, v = mean rate over the last 12 months): today's state propagated with the long-term acceleration of the SMH parabola (LOD +1.78 ms/cy) |
//!
//! Y is the decimal Gregorian year `2000 + (jd − 2451545)/365.2425`. The spline refers to a
//! lunar tidal acceleration ṅ = −25.85″/cy² (LE430, per HMNAO); `DeltaTModel::lunar_ndot`
//! applies the standard correction −0.91072·(ṅ − (−25.85))·u², u = (Y − 1955)/100, to the
//! pre-1955 (spline) part. See `docs/TIME_FRAMES_REPORT.md` for sources and validation.

#[rustfmt::skip] // generated
mod tables;

use crate::math::J2000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Calendar {
    Gregorian,
    Julian,
}

/// Julian day for a calendar date; `hour` is decimal hours. Same semantics as
/// Swiss `julday` (no automatic switch: caller picks the calendar).
pub fn julday(year: i32, month: u32, day: u32, hour: f64, cal: Calendar) -> f64 {
    let (mut y, mut m) = (year as i64, month as i64);
    if m <= 2 {
        y -= 1;
        m += 12;
    }
    let b = match cal {
        Calendar::Gregorian => {
            let a = y.div_euclid(100);
            2 - a + a.div_euclid(4)
        }
        Calendar::Julian => 0,
    };
    // floor(365.25 (y + 4716)) + floor(30.6001 (m + 1)) + day + b − 1524.5, in exact integers.
    let n = (1461 * (y + 4716)).div_euclid(4) + (306_001 * (m + 1)).div_euclid(10_000) + day as i64 + b;
    (n as f64 - 1524.5) + hour / 24.0
}

/// Inverse of `julday`: (year, month, day, decimal hour).
pub fn revjul(jd: f64, cal: Calendar) -> (i32, u32, u32, f64) {
    let z = jd + 0.5;
    let zf = z.floor();
    let hour = (z - zf) * 24.0;
    let zi = zf as i64;
    let a = match cal {
        Calendar::Gregorian => {
            // alpha = floor((Z − 1867216.25) / 36524.25)
            let alpha = (4 * zi - 7_468_865).div_euclid(146_097);
            zi + 1 + alpha - alpha.div_euclid(4)
        }
        Calendar::Julian => zi,
    };
    let b = a + 1524;
    let c = (20 * b - 2442).div_euclid(7305); // floor((B − 122.1)/365.25)
    let d = (1461 * c).div_euclid(4); // floor(365.25 C)
    let e = ((b - d) * 10_000).div_euclid(306_001); // floor((B − D)/30.6001)
    let day = b - d - (306_001 * e).div_euclid(10_000);
    let month = if e < 14 { e - 1 } else { e - 13 };
    let year = if month > 2 { c - 4716 } else { c - 4715 };
    (year as i32, month as u32, day as u32, hour)
}

// ---------------------------------------------------------------------------------------------
// Delta T
// ---------------------------------------------------------------------------------------------

/// Lunar tidal acceleration (″/cy²) the SMH spline refers to (LE430, per HMNAO data files).
pub const SMH_LUNAR_NDOT: f64 = -25.85;

/// Which published Stephenson–Morrison–Hohenkerk spline to use before 1955.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DeltaTSpline {
    /// Stephenson, Morrison & Hohenkerk 2016, Table S15 (−720 … 2016).
    Smh2016,
    /// Morrison, Stephenson, Hohenkerk & Zawilski 2021, Addendum 2020, Table S15.2020 (−720 … 2019).
    Smh2020,
}

/// Delta T model selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeltaTModel {
    pub spline: DeltaTSpline,
    /// Lunar tidal acceleration (″/cy²) of the ephemeris the result is used with. The
    /// pre-1955 part is corrected from `SMH_LUNAR_NDOT` to this value.
    pub lunar_ndot: f64,
    /// Fixed user-defined Delta T (seconds) returned for every date instead of the model
    /// (the classic `set_delta_t_userdef`, as a value instead of global state); `None` = model.
    pub user_seconds: Option<f64>,
}

/// Lunar tidal acceleration of JPL DE440/DE441 (″/cy²), the lunar ephemeris our packs
/// contain: −25.936 per JPL Horizons (J. Giorgini), as quoted in the Swiss Ephemeris
/// programmer's documentation (release notes, "Correct handling of new JPL Ephemeris DE441").
pub const DE441_LUNAR_NDOT: f64 = -25.936;

impl DeltaTModel {
    /// Physically best model: SMH 2020 addendum spline, corrected to the tidal
    /// acceleration of the DE441 Moon actually used (Delta T derived from lunar
    /// observations must be consistent with the lunar ephemeris it is applied to).
    pub const DEFAULT: DeltaTModel =
        DeltaTModel { spline: DeltaTSpline::Smh2020, lunar_ndot: DE441_LUNAR_NDOT, user_seconds: None };
    /// Matches what Swiss Ephemeris 2.10 charts use: SMH 2016 spline corrected to the
    /// DE441 tidal acceleration. Swiss applies the ṅ of the data file covering the date —
    /// −25.936 (DE441) before 1800, −25.80 (DE431) after — so this agrees with Swiss's
    /// in-chart Delta T to ~0.01 s before 1800 and ≤ 0.3 s in 1800–1955 (Moon ≤ 0.15″).
    /// From 1955 on both use observed values; Swiss's post-2026 predictions are not
    /// reproduced (they match no published formula).
    pub const SWISS_COMPAT: DeltaTModel =
        DeltaTModel { spline: DeltaTSpline::Smh2016, lunar_ndot: DE441_LUNAR_NDOT, user_seconds: None };
}

impl Default for DeltaTModel {
    fn default() -> Self {
        DeltaTModel::DEFAULT
    }
}

/// Decimal Gregorian year used as the argument of the Delta T splines.
#[inline]
pub fn decimal_year(jd: f64) -> f64 {
    2000.0 + (jd - J2000) / 365.2425
}

/// Start of the spline → observed blend (decimal year); it ends at the first observed value (1955.0).
const BLEND_START_YEAR: f64 = 1950.0;

fn spline_table(s: DeltaTSpline) -> &'static [tables::SplineRow] {
    match s {
        DeltaTSpline::Smh2016 => &tables::SMH2016,
        DeltaTSpline::Smh2020 => &tables::SMH2020,
    }
}

fn eval_spline(tab: &[tables::SplineRow], y: f64) -> f64 {
    // Rows are contiguous; take the last row with K_i <= y (clamped to the table).
    let i = tab.partition_point(|r| r.0 <= y).saturating_sub(1).min(tab.len() - 1);
    let (k0, k1, a0, a1, a2, a3) = tab[i];
    let t = (y - k0) / (k1 - k0);
    a0 + t * (a1 + t * (a2 + t * a3))
}

fn smh_parabola(y: f64) -> f64 {
    let u = (y - 1825.0) / 100.0;
    -320.0 + 32.5 * u * u
}

/// Correction (s) to a lunar-derived Delta T when the lunar tidal acceleration changes from
/// `ndot_from` to `ndot_to` (″/cy²): −0.91072 (ṅ_to − ṅ_from) u², u = (Y − 1955)/100
/// (Morrison & Stephenson 2004; Explanatory Supplement 3rd ed. §15.3).
pub fn tidal_correction_seconds(jd_ut: f64, ndot_from: f64, ndot_to: f64) -> f64 {
    let u = (decimal_year(jd_ut) - 1955.0) / 100.0;
    -0.91072 * (ndot_to - ndot_from) * u * u
}

/// Historical (spline / long-term) part at decimal year `y`, ṅ-corrected.
fn historical(y: f64, m: &DeltaTModel) -> f64 {
    let tab = spline_table(m.spline);
    let k_first = tab[0].0;
    let base = if y >= k_first {
        eval_spline(tab, y)
    } else {
        smh_parabola(y) + (eval_spline(tab, k_first) - smh_parabola(k_first))
    };
    let u = (y - 1955.0) / 100.0;
    base - 0.91072 * (m.lunar_ndot - SMH_LUNAR_NDOT) * u * u
}

/// Linear interpolation in the observed table (caller guarantees range).
fn observed(jd: f64) -> f64 {
    let obs = &tables::OBSERVED;
    let i = obs.partition_point(|r| r.0 <= jd).clamp(1, obs.len() - 1);
    let (j0, v0) = obs[i - 1];
    let (j1, v1) = obs[i];
    v0 + (v1 - v0) * (jd - j0) / (j1 - j0)
}

/// Long-term quadratic coefficient for the future extrapolation (s/cy², the SMH parabola's 32.5 u²).
const FUTURE_QUADRATIC: f64 = 32.5;

/// Baseline for the extrapolation's linear rate (days): the mean rate over the last 10 years of
/// the table. A one-year rate swings with Earth's decade fluctuations (0.04 to 0.20 s/yr in
/// 2006–2027); 5- and 10-year baselines agree to 0.6 s at 2100.
const FUTURE_RATE_BASELINE_DAYS: f64 = 10.0 * 365.25;

/// JD (UT) range covered by observed values: (first, last).
pub fn delta_t_observed_range() -> (f64, f64) {
    (tables::OBSERVED[0].0, tables::OBSERVED[tables::OBSERVED.len() - 1].0)
}

/// Delta T = TT − UT1 in seconds for a JD in UT with an explicit model.
pub fn delta_t_seconds_with(jd_ut: f64, m: &DeltaTModel) -> f64 {
    if let Some(s) = m.user_seconds {
        return s;
    }
    let obs = &tables::OBSERVED;
    let (j_first, v_first) = obs[0];
    let (j_last, v_last) = obs[obs.len() - 1];
    if jd_ut >= j_first && jd_ut <= j_last {
        return observed(jd_ut);
    }
    if jd_ut > j_last {
        let v = (v_last - observed(j_last - FUTURE_RATE_BASELINE_DAYS)) * 36525.0 / FUTURE_RATE_BASELINE_DAYS; // s/cy
        let tau = (jd_ut - j_last) / 36525.0;
        return v_last + v * tau + FUTURE_QUADRATIC * tau * tau;
    }
    let y = decimal_year(jd_ut);
    let h = historical(y, m);
    if y > BLEND_START_YEAR {
        let y1 = decimal_year(j_first);
        let w = (y - BLEND_START_YEAR) / (y1 - BLEND_START_YEAR);
        h + w * (v_first - historical(y1, m))
    } else {
        h
    }
}

/// Delta T = TT − UT1 in seconds, for a JD in UT, using the default versioned model.
pub fn delta_t_seconds(jd_ut: f64) -> f64 {
    delta_t_seconds_with(jd_ut, &DeltaTModel::DEFAULT)
}

/// JD(TT) from JD(UT1).
pub fn ut_to_tt(jd_ut: f64) -> f64 {
    jd_ut + delta_t_seconds(jd_ut) / 86400.0
}

/// JD(TT) from JD(UT1) with an explicit Delta T model.
pub fn ut_to_tt_with(jd_ut: f64, m: &DeltaTModel) -> f64 {
    jd_ut + delta_t_seconds_with(jd_ut, m) / 86400.0
}

/// JD(UT1) from JD(TT) with an explicit Delta T model (fixed-point iteration; |dΔT/dt| < 1e-6,
/// so it converges to < 1e-12 day in a few steps).
pub fn tt_to_ut_with(jd_tt: f64, m: &DeltaTModel) -> f64 {
    let mut ut = jd_tt - delta_t_seconds_with(jd_tt, m) / 86400.0;
    for _ in 0..5 {
        let next = jd_tt - delta_t_seconds_with(ut, m) / 86400.0;
        if (next - ut).abs() < 1e-12 {
            return next;
        }
        ut = next;
    }
    ut
}

/// JD(UT1) from JD(TT) (iterative inverse of `ut_to_tt`).
pub fn tt_to_ut(jd_tt: f64) -> f64 {
    tt_to_ut_with(jd_tt, &DeltaTModel::DEFAULT)
}

// ---------------------------------------------------------------------------------------------
// TDB − TT
// ---------------------------------------------------------------------------------------------

fn fairhead(terms: &[(f64, f64, f64, u8)], jd_tt: f64) -> f64 {
    // Julian millennia from J2000 (TT stands in for TDB: the difference is negligible here).
    let t = (jd_tt - J2000) / 365_250.0;
    let mut w = [0.0f64; 5];
    for &(a, f, p, n) in terms.iter().rev() {
        w[n as usize] += a * (f * t + p).sin();
    }
    let wf = t * (t * (t * (t * w[4] + w[3]) + w[2]) + w[1]) + w[0];
    // Adjustments to use JPL planetary masses instead of IAU (IERS Conventions / SOFA dtdb).
    let wj = 0.00065e-6 * (6069.776754 * t + 4.021194).sin() + 0.00033e-6 * (213.299095 * t + 5.543132).sin()
        - 0.00196e-6 * (6208.294251 * t + 5.696701).sin()
        - 0.00173e-6 * (74.781599 * t + 2.435900).sin()
        + 0.03638e-6 * t * t;
    wf + wj
}

/// TDB − TT in seconds (periodic terms; |value| < 2 ms), geocentric.
///
/// Fairhead & Bretagnon (1990) series truncated to 87 terms (see `tables::FB_TRUNC`);
/// max error vs the full 787-term series 0.42 µs over years 0–2400.
pub fn tdb_minus_tt_seconds(jd_tt: f64) -> f64 {
    fairhead(&tables::FB_TRUNC, jd_tt)
}

/// TDB − TT in seconds from the full 787-term Fairhead & Bretagnon series (geocentric;
/// equals SOFA/ERFA `dtdb(jd, 0, 0, 0, 0, 0)`).
pub fn tdb_minus_tt_seconds_full(jd_tt: f64) -> f64 {
    fairhead(&tables::FB_FULL, jd_tt)
}

/// JD(TDB) from JD(TT).
pub fn tt_to_tdb(jd_tt: f64) -> f64 {
    jd_tt + tdb_minus_tt_seconds(jd_tt) / 86400.0
}

/// Identifies the Delta T model in diagnostics.
pub const DELTA_T_MODEL: &str =
    "FE-DT/2: SMH2020 spline (-720..1955, ndot -25.85) + USNO observed (1955.0..2026-10-01) + \
     IERS Bulletin A predictions (to 2027-11-01) + extrapolation (last + v*tau + 32.5*tau^2 s, \
     v = mean rate over the last 10 years, tau in Julian centuries)";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn julday_known_dates() {
        assert_eq!(julday(2000, 1, 1, 12.0, Calendar::Gregorian), 2451545.0);
        assert_eq!(julday(-4712, 1, 1, 12.0, Calendar::Julian), 0.0);
        assert_eq!(julday(1582, 10, 15, 0.0, Calendar::Gregorian), 2299160.5);
        assert_eq!(julday(1582, 10, 4, 0.0, Calendar::Julian), 2299159.5);
    }

    #[test]
    fn revjul_roundtrip() {
        for &cal in &[Calendar::Gregorian, Calendar::Julian] {
            let mut jd = 1_000_000.25;
            while jd < 2_700_000.0 {
                let (y, m, d, h) = revjul(jd, cal);
                assert_eq!(julday(y, m, d, h, cal), jd);
                jd += 3337.375;
            }
        }
    }

    #[test]
    fn delta_t_is_continuous() {
        let (j0, j1) = delta_t_observed_range();
        for &j in &[j0, j1, julday(1950, 1, 1, 0.0, Calendar::Gregorian)] {
            let a = delta_t_seconds(j - 1e-6);
            let b = delta_t_seconds(j + 1e-6);
            assert!((a - b).abs() < 1e-5, "jump at {j}: {a} {b}");
        }
        let jm720 = J2000 + (-720.0 - 2000.0) * 365.2425;
        assert!((delta_t_seconds(jm720 - 1e-6) - delta_t_seconds(jm720 + 1e-6)).abs() < 1e-3);
    }

    #[test]
    fn tt_ut_inverse() {
        for &y in &[-1.0, 0.0, 800.0, 1620.0, 1999.0, 2026.9, 2100.0, 2400.0] {
            let ut = J2000 + (y - 2000.0) * 365.25 + 0.3;
            assert!((tt_to_ut(ut_to_tt(ut)) - ut).abs() < 1e-9);
        }
    }
}
