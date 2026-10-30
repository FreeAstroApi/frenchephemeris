// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Small formatting / normalisation helpers and name tables with the semantics of the classic
//! astrology-engine API (`split_deg`, centisecond helpers, `degnorm`, `difdeg2n`, body and
//! ayanamsa names). Behaviour was established black-box against the Swiss Ephemeris 2.10
//! oracle (pyswisseph 2.10.3.2 and its exported C symbols, public API only); see
//! `tests/fixtures/quickwins/gen_quickwins.py`. Integer helpers use 32-bit wrapping
//! arithmetic and C (truncating) division like the original, so out-of-range inputs give
//! the same results.

use crate::math::{PI, TAU};

/// Centiseconds of arc: 1/100 arcsecond.
pub type Centisec = i32;

pub const DEG30_CS: i32 = 30 * 360_000;
pub const DEG180_CS: i32 = 180 * 360_000;
pub const DEG360_CS: i32 = 360 * 360_000;

pub const SPLIT_DEG_ROUND_SEC: i32 = 1;
pub const SPLIT_DEG_ROUND_MIN: i32 = 2;
pub const SPLIT_DEG_ROUND_DEG: i32 = 4;
pub const SPLIT_DEG_ZODIACAL: i32 = 8;
pub const SPLIT_DEG_KEEP_SIGN: i32 = 16;
pub const SPLIT_DEG_KEEP_DEG: i32 = 32;
pub const SPLIT_DEG_NAKSHATRA: i32 = 1024;

/// Result of `split_deg`: degrees, minutes, seconds, fraction of second (or the rounded
/// seconds as a float when a rounding flag is set) and the sign: zodiac sign 0–11 /
/// nakshatra 0–26 with `ZODIACAL` / `NAKSHATRA`, else ±1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SplitDeg {
    pub deg: i32,
    pub min: i32,
    pub sec: i32,
    pub secfr: f64,
    pub sign: i32,
}

/// Split decimal degrees into (sign/nakshatra,) degrees, minutes, seconds.
///
/// * The absolute value is split; without a zodiac/nakshatra flag `sign` is −1 for negative
///   input, else 1. Nakshatra splitting applies to non-negative input only (negative input
///   falls back to the zodiac split if `ZODIACAL` is also set, else to ±1).
/// * Rounding adds half a unit (`ROUND_DEG` > `ROUND_MIN` > `ROUND_SEC`) before the value
///   is split into a sign (30°) or nakshatra (13°20′) index and the remainder. The lower
///   fields are those of the shifted value (`ROUND_MIN` leaves the seconds of value + 30″),
///   and `secfr` is the rounded seconds as a float. `KEEP_DEG` skips the rounding if it
///   would carry the remainder into its next degree, `KEEP_SIGN` if it would carry into the
///   next sign / nakshatra (30° steps when splitting by neither).
/// * A sign index of exactly 12 (zodiac) / 27 (nakshatra) — 360° — wraps to 0; larger
///   indices (input ≥ 390°) are returned as they are.
/// * The nakshatra width is 13.33333333333333° (40/3 written to 16 significant digits), as
///   the original's results at exact boundaries show (360° → nakshatra 0).
pub fn split_deg(ddeg: f64, roundflag: i32) -> SplitDeg {
    const NAK: f64 = 13.333_333_333_333_33;
    let x = ddeg.abs();
    let nakshatra = roundflag & SPLIT_DEG_NAKSHATRA != 0 && ddeg >= 0.0;
    let zodiacal = roundflag & SPLIT_DEG_ZODIACAL != 0 && !nakshatra;
    let rounding = roundflag & (SPLIT_DEG_ROUND_SEC | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_DEG) != 0;
    let unit = if nakshatra { NAK } else { 30.0 };
    let split = |x: f64| -> (i64, f64) {
        if nakshatra || zodiacal {
            let n = (x / unit) as i64;
            (n, x - n as f64 * unit)
        } else {
            (0, x)
        }
    };
    let mut x = x;
    if rounding {
        let r = if roundflag & SPLIT_DEG_ROUND_DEG != 0 {
            0.5
        } else if roundflag & SPLIT_DEG_ROUND_MIN != 0 {
            0.5 / 60.0
        } else {
            0.5 / 3600.0
        };
        let (_, rel) = split(x);
        let crosses_deg = ((rel + r) as i64) > (rel as i64);
        let crosses_sign = if nakshatra || zodiacal {
            rel + r >= unit
        } else {
            (((x + r) / 30.0) as i64) > ((x / 30.0) as i64)
        };
        let keep = (roundflag & SPLIT_DEG_KEEP_DEG != 0 && crosses_deg)
            || (roundflag & SPLIT_DEG_KEEP_SIGN != 0 && crosses_sign);
        if !keep {
            x += r;
        }
    }
    let (n, rel) = split(x);
    let sign = if nakshatra || zodiacal {
        let full = if nakshatra { 27 } else { 12 };
        if n == full {
            0
        } else {
            n as i32
        }
    } else if ddeg < 0.0 {
        -1
    } else {
        1
    };
    let deg = rel as i32;
    let mut f = (rel - deg as f64) * 60.0;
    let min = f as i32;
    f = (f - min as f64) * 60.0;
    let sec = f as i32;
    let secfr = if rounding { sec as f64 } else { f - sec as f64 };
    SplitDeg { deg, min, sec, secfr, sign }
}

/// Normalise degrees to [0, 360). Values within 1e-13 of a multiple of 360 give 0.
pub fn degnorm(x: f64) -> f64 {
    let mut y = x % 360.0;
    if y.abs() < 1e-13 {
        y = 0.0;
    }
    if y < 0.0 {
        y += 360.0;
    }
    y
}

/// Normalise radians to [0, 2π).
pub fn radnorm(x: f64) -> f64 {
    let mut y = x % TAU;
    if y.abs() < 1e-13 {
        y = 0.0;
    }
    if y < 0.0 {
        y += TAU;
    }
    y
}

/// p1 − p2 normalised to [0, 360).
pub fn difdegn(p1: f64, p2: f64) -> f64 {
    degnorm(p1 - p2)
}

/// p1 − p2 normalised to [−180, 180).
pub fn difdeg2n(p1: f64, p2: f64) -> f64 {
    let d = degnorm(p1 - p2);
    if d >= 180.0 {
        d - 360.0
    } else {
        d
    }
}

/// p1 − p2 (radians) normalised to [−π, π).
pub fn difrad2n(p1: f64, p2: f64) -> f64 {
    let d = radnorm(p1 - p2);
    if d >= PI {
        d - TAU
    } else {
        d
    }
}

/// Centiseconds normalised to [0, 360°).
pub fn csnorm(p: Centisec) -> Centisec {
    p.rem_euclid(DEG360_CS)
}

/// p1 − p2 in centiseconds, normalised to [0, 360°).
pub fn difcsn(p1: Centisec, p2: Centisec) -> Centisec {
    csnorm(p1.wrapping_sub(p2))
}

/// p1 − p2 in centiseconds, normalised to [−180°, 180°).
pub fn difcs2n(p1: Centisec, p2: Centisec) -> Centisec {
    let d = csnorm(p1.wrapping_sub(p2));
    if d >= DEG180_CS {
        d - DEG360_CS
    } else {
        d
    }
}

/// Round centiseconds to whole seconds of arc, except that 29°59′59.5″ (and the like) is
/// rounded down rather than into the next 30° sign.
pub fn csroundsec(x: Centisec) -> Centisec {
    let t = (x.wrapping_add(50) / 100).wrapping_mul(100);
    if t > x && t % DEG30_CS == 0 {
        (x / 100) * 100
    } else {
        t
    }
}

/// Double to 32-bit integer, rounding half away from zero (saturating; NaN → 0).
pub fn d2l(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        -((-x + 0.5) as i32)
    }
}

/// C-style digit character for a (possibly negative) single-digit value, like `'0' + d`.
fn digit(d: i32) -> char {
    char::from_u32((b'0' as i32 + d) as u32).unwrap_or('?')
}

/// Time string "hh:mm:ss" (`sep` between fields) from centiseconds of time, rounded to the
/// second; hours modulo 24. `suppress_zero` drops ":ss" when the seconds are 0.
pub fn cs2timestr(t: Centisec, sep: char, suppress_zero: bool) -> String {
    let t = t.wrapping_add(50) / 100;
    let (h, m, s) = ((t / 3600) % 24, (t / 60) % 60, t % 60);
    let mut out = String::new();
    out.push(digit(h / 10));
    out.push(digit(h % 10));
    out.push(sep);
    out.push(digit(m / 10));
    out.push(digit(m % 10));
    if s != 0 || !suppress_zero {
        out.push(sep);
        out.push(digit(s / 10));
        out.push(digit(s % 10));
    }
    out
}

/// Longitude/latitude string "dddPmm'ss" from centiseconds of arc: `plus` / `minus` is the
/// hemisphere letter, degrees are printed without padding (modulo 1000), seconds rounded and
/// omitted when 0.
pub fn cs2lonlatstr(t: Centisec, plus: char, minus: char) -> String {
    let (letter, t) = if t < 0 { (minus, t.wrapping_neg()) } else { (plus, t) };
    let t = t.wrapping_add(50) / 100;
    let (d, m, s) = ((t / 3600) % 1000, (t / 60) % 60, t % 60);
    let mut out = String::new();
    if d >= 100 {
        out.push(digit(d / 100));
    }
    if d >= 10 {
        out.push(digit((d / 10) % 10));
    }
    out.push(digit(d % 10));
    out.push(letter);
    out.push(digit(m / 10));
    out.push(digit(m % 10));
    if s != 0 {
        out.push('\'');
        out.push(digit(s / 10));
        out.push(digit(s % 10));
    }
    out
}

/// Degrees-within-sign string "dd°mm'ss" (C `%2d°%02d'%02d`, truncated to the second).
pub fn cs2degstr(t: Centisec) -> String {
    let t = t / 100;
    format!("{:2}°{:02}'{:02}", (t / 3600) % 30, (t / 60) % 60, t % 60)
}

/// Body name in the classic API spelling, for ids FrenchEphemeris knows (−1 = obliquity /
/// nutation pseudo-body gives "", as the original). Numbered asteroids without a known name
/// give their MPC number. Unknown ids give "".
pub fn planet_name(id: i32) -> String {
    const NAMES: [&str; 23] = [
        "Sun",
        "Moon",
        "Mercury",
        "Venus",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
        "mean Node",
        "true Node",
        "mean Apogee",
        "osc. Apogee",
        "Earth",
        "Chiron",
        "Pholus",
        "Ceres",
        "Pallas",
        "Juno",
        "Vesta",
        "intp. Apogee",
        "intp. Perigee",
    ];
    match id {
        0..=22 => NAMES[id as usize].to_string(),
        10001 => "Ceres".into(),
        10002 => "Pallas".into(),
        10003 => "Juno".into(),
        10004 => "Vesta".into(),
        10433 => "Eros".into(),
        12060 => "Chiron".into(),
        15145 => "Pholus".into(),
        146199 => "Eris".into(),
        n if n > crate::context::AST_OFFSET => (n - crate::context::AST_OFFSET).to_string(),
        _ => String::new(),
    }
}

/// Sidereal-mode name for the classic `SE_SIDM_*` number (the original's
/// labels, swephprg.htm §12.1) for every mode 0–46 (option bits ignored);
/// "" for other numbers and SE_SIDM_USER.
pub fn ayanamsa_name(sid_mode: i32) -> &'static str {
    crate::ayanamsa::Ayanamsa::from_swiss(sid_mode & 0xFF).map_or("", |m| m.name())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_basic() {
        let s = split_deg(12.3456789, 0);
        assert_eq!((s.deg, s.min, s.sec, s.sign), (12, 20, 44, 1));
        let s = split_deg(29.9999998, SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_SEC);
        assert_eq!((s.deg, s.min, s.sec, s.sign), (0, 0, 0, 1));
        let s = split_deg(29.9999998, SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_SEC | SPLIT_DEG_KEEP_SIGN);
        assert_eq!((s.deg, s.min, s.sec, s.sign), (29, 59, 59, 0));
    }

    #[test]
    fn strings() {
        assert_eq!(cs2degstr(123456789), "12°56'07");
        assert_eq!(cs2lonlatstr(123456789, 'E', 'W'), "342E56'08");
        assert_eq!(cs2timestr(123456789, ':', false), "06:56:08");
        assert_eq!(cs2timestr(-123456789, ':', false), "0*:+*:0)");
        assert_eq!(csroundsec(10799950), 10799900);
    }
}
