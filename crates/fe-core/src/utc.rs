// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! UTC with leap seconds, time zones and calendar-date checks (the classic `utc_to_jd`,
//! `jdet_to_utc`, `jdut1_to_utc`, `utc_time_zone`, `date_conversion`, `day_of_week`).
//!
//! # Leap seconds
//!
//! `LEAP_SECONDS` is TAI − UTC from IERS Bulletin C, file `Leap_Second.dat`
//! (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat), "Updated through IERS Bulletin 72
//! issued in July 2026", "File expires on 28 June 2027", SHA-256
//! 6cb6f5d4b819f2e568e25db4b0b26d89dedf031fdffb18bc94d40f4e94e268d7. The last leap second
//! was inserted at the end of 2016-12-31 (TAI − UTC = 37 s since 2017-01-01).
//!
//! # Semantics (Swiss Ephemeris programmer's documentation §9.2–9.3)
//!
//! * Before 1972-01-01 the input/output time is treated as UT1 (UTC with rubber seconds,
//!   1961–1971, stayed within 0.1 s of UT1).
//! * From 1972: TT = UTC + (TAI − UTC) + 32.184 s; UT1 = TT − ΔT.
//! * If ΔT − (TAI − UTC) − 32.184 s > 1 s (a leap second the table does not know yet, i.e.
//!   the future), the time is treated as UT1 so the error stays below 1 s.
//! * A leap second is written 23:59:60.x on the day it is inserted; 60 ≤ s < 61 on any
//!   other minute is an error.
//!
//! Delta T comes from the caller's `DeltaTModel` (the engine's), so the UT1 values follow
//! the same model as every other UT computation.

use crate::time::{self, Calendar, DeltaTModel};
use crate::{Error, Result};

/// (JD of 0h UTC from which the value applies, TAI − UTC in seconds), IERS Bulletin C 72.
pub const LEAP_SECONDS: [(f64, f64); 28] = [
    (2441317.5, 10.0), // 1972-01-01
    (2441499.5, 11.0), // 1972-07-01
    (2441683.5, 12.0), // 1973-01-01
    (2442048.5, 13.0), // 1974-01-01
    (2442413.5, 14.0), // 1975-01-01
    (2442778.5, 15.0), // 1976-01-01
    (2443144.5, 16.0), // 1977-01-01
    (2443509.5, 17.0), // 1978-01-01
    (2443874.5, 18.0), // 1979-01-01
    (2444239.5, 19.0), // 1980-01-01
    (2444786.5, 20.0), // 1981-07-01
    (2445151.5, 21.0), // 1982-07-01
    (2445516.5, 22.0), // 1983-07-01
    (2446247.5, 23.0), // 1985-07-01
    (2447161.5, 24.0), // 1988-01-01
    (2447892.5, 25.0), // 1990-01-01
    (2448257.5, 26.0), // 1991-01-01
    (2448804.5, 27.0), // 1992-07-01
    (2449169.5, 28.0), // 1993-07-01
    (2449534.5, 29.0), // 1994-07-01
    (2450083.5, 30.0), // 1996-01-01
    (2450630.5, 31.0), // 1997-07-01
    (2451179.5, 32.0), // 1999-01-01
    (2453736.5, 33.0), // 2006-01-01
    (2454832.5, 34.0), // 2009-01-01
    (2456109.5, 35.0), // 2012-07-01
    (2457204.5, 36.0), // 2015-07-01
    (2457754.5, 37.0), // 2017-01-01
];

/// Identifies the leap-second table in diagnostics.
pub const LEAP_SECONDS_SOURCE: &str = "IERS Bulletin C 72 (Leap_Second.dat, July 2026, expires 2027-06-28)";

/// TT − TAI, seconds.
pub const TT_MINUS_TAI: f64 = 32.184;

const JD_1972: f64 = 2441317.5;

/// TAI − UTC (s) for a UTC Julian date (0 before 1972).
pub fn tai_minus_utc(jd_utc: f64) -> f64 {
    let i = LEAP_SECONDS.partition_point(|r| r.0 <= jd_utc);
    if i == 0 {
        0.0
    } else {
        LEAP_SECONDS[i - 1].1
    }
}

/// Whether a leap second is inserted at the end of the UTC day starting at `jd_day0` (0h).
fn leap_second_after(jd_day0: f64) -> bool {
    LEAP_SECONDS[1..].iter().any(|r| r.0 == jd_day0 + 1.0)
}

/// Date check: `(valid, jd, (y, m, d, hour))` — the classic `date_conversion`. `jd` is
/// `julday` of the input; if the date does not exist in the calendar or the hour is outside
/// [0, 24), `valid` is false and the tuple is the normalised date of `jd`.
pub fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    cal: Calendar,
) -> (bool, f64, (i32, i32, i32, f64)) {
    let jd = julday_any(year, month, day, hour, cal);
    let (y, m, d, _) = time::revjul(julday_any(year, month, day, 0.0, cal), cal);
    let valid = (y, m as i32, d as i32) == (year, month, day) && (0.0..24.0).contains(&hour);
    if valid {
        (true, jd, (year, month, day, hour))
    } else {
        let (y, m, d, h) = time::revjul(jd, cal);
        (false, jd, (y, m as i32, d as i32, h))
    }
}

/// `julday` accepting out-of-range months/days (normalised the way the formula does).
fn julday_any(year: i32, month: i32, day: i32, hour: f64, cal: Calendar) -> f64 {
    // Months outside 1..12 roll the year; days outside the month are linear offsets.
    let m0 = month - 1;
    let (y, m) = (year + m0.div_euclid(12), m0.rem_euclid(12) + 1);
    if day >= 1 {
        time::julday(y, m as u32, day as u32, hour, cal)
    } else {
        time::julday(y, m as u32, 1, hour, cal) + (day - 1) as f64
    }
}

/// Day of week, Monday = 0 … Sunday = 6, of the civil day containing `jd`.
pub fn day_of_week(jd: f64) -> i32 {
    ((jd + 0.5).floor() as i64).rem_euclid(7) as i32
}

/// Split a decimal hour into (h, min, sec).
fn split_hour(hour: f64) -> (i32, i32, f64) {
    let h = hour.floor();
    let mf = (hour - h) * 60.0;
    let m = mf.floor();
    let s = (mf - m) * 60.0;
    (h as i32, m as i32, s)
}

/// Local time ↔ UTC: subtract `timezone_hours` (east positive) from the date and time. Use
/// +offset for local → UTC, −offset for UTC → local. A leap second (s ≥ 60) is carried
/// through unchanged. Gregorian calendar.
pub fn utc_time_zone(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    min: i32,
    sec: f64,
    timezone_hours: f64,
) -> (i32, i32, i32, i32, i32, f64) {
    let leap = sec >= 60.0;
    let s = if leap { sec - 1.0 } else { sec };
    // Shift the time of day in hours and carry whole days through the calendar, so the
    // seconds keep their precision (no round trip through a JD).
    let h = hour as f64 + min as f64 / 60.0 + s / 3600.0 - timezone_hours;
    let days = (h / 24.0).floor();
    let h = h - days * 24.0;
    let jd0 = julday_any(year, month, day, 0.0, Calendar::Gregorian) + days;
    let (y, m, d, _) = time::revjul(jd0 + 0.5, Calendar::Gregorian);
    let (h, mi, mut s) = split_hour(h);
    if leap {
        s += 1.0;
    }
    (y, m as i32, d as i32, h, mi, s)
}

/// UTC date and time → (JD TT, JD UT1). Errors for impossible dates/times, and for second
/// 60 outside a leap second.
#[allow(clippy::too_many_arguments)]
pub fn utc_to_jd(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    min: i32,
    sec: f64,
    cal: Calendar,
    dt: &DeltaTModel,
) -> Result<(f64, f64)> {
    let (ok, jd_day0, _) = date_conversion(year, month, day, 0.0, cal);
    if !ok {
        return Err(Error::InvalidArgument(format!(
            "invalid date: year = {year}, month = {month}, day = {day}"
        )));
    }
    if !(0..24).contains(&hour) || !(0..60).contains(&min) || !(0.0..61.0).contains(&sec) {
        return Err(Error::InvalidArgument(format!("invalid time: {hour}:{min}:{sec:.2}")));
    }
    let jd = jd_day0 + (hour as f64 + min as f64 / 60.0 + sec / 3600.0) / 24.0;
    let as_ut1 = |ut: f64| (ut + time::delta_t_seconds_with(ut, dt) / 86400.0, ut);
    if sec >= 60.0 && (jd_day0 < JD_1972 || hour != 23 || min != 59 || !leap_second_after(jd_day0)) {
        return Err(Error::InvalidArgument(format!("invalid time (no leap second!): {hour}:{min}:{sec:.2}")));
    }
    if jd_day0 < JD_1972 {
        return Ok(as_ut1(jd));
    }
    // TAI − UTC of the day (a leap second still counts with the old value).
    let tai_utc = tai_minus_utc(jd_day0);
    let tt = jd + (tai_utc + TT_MINUS_TAI) / 86400.0;
    let ut = time::tt_to_ut_with(tt, dt);
    if (tt - ut) * 86400.0 - tai_utc - TT_MINUS_TAI > 1.0 {
        return Ok(as_ut1(jd));
    }
    Ok((tt, ut))
}

/// JD TT → UTC (year, month, day, hour, min, sec); sec may be 60.x during a leap second.
pub fn jdet_to_utc(jd_tt: f64, cal: Calendar, dt: &DeltaTModel) -> (i32, i32, i32, i32, i32, f64) {
    let ut = time::tt_to_ut_with(jd_tt, dt);
    let as_ut1 = || {
        let (y, m, d, h) = time::revjul(ut, cal);
        let (h, mi, s) = split_hour(h);
        (y, m as i32, d as i32, h, mi, s)
    };
    let first_tt = JD_1972 + (LEAP_SECONDS[0].1 + TT_MINUS_TAI) / 86400.0;
    if jd_tt < first_tt {
        return as_ut1();
    }
    // Last table entry whose insertion instant (in TT, with its new value) is <= jd_tt.
    let i = LEAP_SECONDS.partition_point(|r| r.0 + (r.1 + TT_MINUS_TAI) / 86400.0 <= jd_tt) - 1;
    let (_, val) = LEAP_SECONDS[i];
    if (jd_tt - ut) * 86400.0 - val - TT_MINUS_TAI > 1.0 {
        return as_ut1();
    }
    let mut leap = false;
    let mut utc = jd_tt - (val + TT_MINUS_TAI) / 86400.0;
    if let Some(&(next_jd, next_val)) = LEAP_SECONDS.get(i + 1) {
        // Within the inserted second: UTC (old count) falls in [next_jd, next_jd + 1 s).
        let step = next_val - val;
        if utc >= next_jd && step > 0.0 {
            leap = true;
            utc -= step / 86400.0;
        }
    }
    let (y, m, d, h) = time::revjul(utc, cal);
    let (h, mi, mut s) = split_hour(h);
    if leap {
        s += 1.0;
    }
    (y, m as i32, d as i32, h, mi, s)
}

/// JD UT1 → UTC (via TT with the given Delta T model).
pub fn jdut1_to_utc(jd_ut: f64, cal: Calendar, dt: &DeltaTModel) -> (i32, i32, i32, i32, i32, f64) {
    jdet_to_utc(jd_ut + time::delta_t_seconds_with(jd_ut, dt) / 86400.0, cal, dt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leap_second_round_trip() {
        let m = DeltaTModel::SWISS_COMPAT;
        let (tt, _) = utc_to_jd(2016, 12, 31, 23, 59, 60.5, Calendar::Gregorian, &m).unwrap();
        let (y, mo, d, h, mi, s) = jdet_to_utc(tt, Calendar::Gregorian, &m);
        assert_eq!((y, mo, d, h, mi), (2016, 12, 31, 23, 59));
        assert!((s - 60.5).abs() < 1e-4, "{s}");
        let (tt2, _) = utc_to_jd(2017, 1, 1, 0, 0, 0.0, Calendar::Gregorian, &m).unwrap();
        assert!(((tt2 - tt) * 86400.0 - 0.5).abs() < 1e-4);
        assert!(utc_to_jd(2015, 6, 29, 23, 59, 60.0, Calendar::Gregorian, &m).is_err());
    }

    #[test]
    fn dow() {
        assert_eq!(day_of_week(2451545.0), 5); // Saturday 2000-01-01
        assert_eq!(day_of_week(-1.0), 6);
    }
}
