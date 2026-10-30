// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Sidereal zodiac offsets (ayanamsas): every Swiss-numbered sidereal mode
//! (0–46), a user-defined mode, and the sidereal option bits.
//!
//! CONTRACT (owned by the astrology-geometry workstream):
//! - Degrees. `ayanamsa_mean_deg` excludes nutation; the adapter adds
//!   nutation in longitude when the requested positions are apparent
//!   (true ayanamsa = mean + Δψ for every mode, as in Swiss).
//!
//! Kinds of definition (see `docs/GEOMETRY_REPORT.md` §2 for the provenance
//! of every constant):
//!
//! 1. **Fixed epoch** (t0, A0): the sidereal zero point is fixed on the mean
//!    ecliptic of t0, A0 degrees behind the mean equinox of t0, and
//!    `ayanamsa(t) = A0 − λ_{ecl(t0)}(mean equinox of t)` with Vondrák et al.
//!    (2011) precession. t0 is TT, or a calendar instant in UT converted with
//!    `DeltaTModel::SWISS_COMPAT` (the published definitions give dates; the
//!    UT reading was chosen black-box among {UT, TT}).
//! 2. **Star-anchored** ("true" ayanamsas): a catalogue star's apparent
//!    longitude on the mean ecliptic of date minus its defined sidereal
//!    longitude (True Citra, Revati, Pushya, Mula; Galactic Centre modes use
//!    Sgr A*, [`galactic`]).
//! 3. **Galactic equator**: the node of the galactic equator on the mean
//!    ecliptic of date minus its defined sidereal longitude.
//! 4. **Projection** (J2000, J1900, B1950, Skydram/Mardyks by default, and any
//!    fixed-epoch mode with [`SidBits::ECL_T0`] or [`SidBits::SSY_PLANE`]):
//!    positions are projected onto the mean ecliptic of t0 (or the invariable
//!    plane) instead of subtracting an ayanamsa; the ayanamsa value itself is
//!    still the kind-1 value.

pub mod galactic;
pub mod orig_precession;
pub mod star;
pub mod vondrak;

use crate::math::{norm_deg, DEG};
use crate::time::{julday, Calendar, DeltaTModel};

/// Sidereal mode. Numbering (`swiss_number`) and display names
/// (`name`) follow the Swiss Ephemeris public documentation (swephprg.htm
/// §12.1, `swe_get_ayanamsa_name`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Ayanamsa {
    FaganBradley,
    Lahiri,
    DeLuce,
    Raman,
    UshaShashi,
    Krishnamurti,
    DjwhalKhul,
    Yukteshwar,
    JnBhasin,
    BabylKugler1,
    BabylKugler2,
    BabylKugler3,
    BabylHuber,
    BabylEtaPiscium,
    BabylAldebaran15Tau,
    Hipparchos,
    Sassanian,
    GalCent0Sag,
    J2000,
    J1900,
    B1950,
    Suryasiddhanta,
    SuryasiddhantaMeanSun,
    Aryabhata,
    AryabhataMeanSun,
    SsRevati,
    SsCitra,
    TrueCitra,
    TrueRevati,
    TruePushya,
    GalCentRgBrand,
    GalEquIau1958,
    GalEquTrue,
    GalEquMula,
    GalAlignMardyks,
    TrueMula,
    GalCentMulaWilhelm,
    Aryabhata522,
    BabylBritton,
    TrueSheoran,
    GalCentCochrane,
    GalEquFiorenza,
    ValensMoon,
    Lahiri1940,
    LahiriVp285,
    KrishnamurtiVp291,
    LahiriIcrc,
    /// FreeAstroAPI maps `"pushya"` to Swiss `SIDM_TRUE_PUSHYA` as well, so
    /// this is the same mode as `TruePushya`.
    Pushya,
    /// User-defined (Swiss `SE_SIDM_USER`): ayanamsa `ayan_t0` (deg) at `t0`
    /// (JD TT, or UT when `t0_is_ut`, i.e. `SIDBIT_USER_UT`).
    User {
        t0: f64,
        ayan_t0: f64,
        t0_is_ut: bool,
    },
}

/// Swiss sidereal option bits, passed per call (no global state).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SidBits(pub u32);

impl SidBits {
    /// Project positions onto the mean ecliptic of t0.
    pub const ECL_T0: u32 = 256;
    /// Project positions onto the solar-system invariable plane.
    pub const SSY_PLANE: u32 = 512;
    /// User mode: t0 is UT.
    pub const USER_UT: u32 = 1024;
    /// Ayanamsa measured on the ecliptic of date ("consistent" method 2).
    pub const ECL_DATE: u32 = 2048;
    /// Do not apply the precession-model reconversion offset.
    pub const NO_PREC_OFFSET: u32 = 4096;
    /// Use the definition's original precession model (test option).
    pub const PREC_ORIG: u32 = 8192;
    pub const ALL: u32 = 256 | 512 | 1024 | 2048 | 4096 | 8192;
    pub fn has(self, b: u32) -> bool {
        self.0 & b != 0
    }
}

/// Swiss `SE_SIDM_USER`.
pub const SIDM_USER: i32 = 255;

const ALL_MODES: [Ayanamsa; 47] = {
    use Ayanamsa::*;
    [
        FaganBradley,
        Lahiri,
        DeLuce,
        Raman,
        UshaShashi,
        Krishnamurti,
        DjwhalKhul,
        Yukteshwar,
        JnBhasin,
        BabylKugler1,
        BabylKugler2,
        BabylKugler3,
        BabylHuber,
        BabylEtaPiscium,
        BabylAldebaran15Tau,
        Hipparchos,
        Sassanian,
        GalCent0Sag,
        J2000,
        J1900,
        B1950,
        Suryasiddhanta,
        SuryasiddhantaMeanSun,
        Aryabhata,
        AryabhataMeanSun,
        SsRevati,
        SsCitra,
        TrueCitra,
        TrueRevati,
        TruePushya,
        GalCentRgBrand,
        GalEquIau1958,
        GalEquTrue,
        GalEquMula,
        GalAlignMardyks,
        TrueMula,
        GalCentMulaWilhelm,
        Aryabhata522,
        BabylBritton,
        TrueSheoran,
        GalCentCochrane,
        GalEquFiorenza,
        ValensMoon,
        Lahiri1940,
        LahiriVp285,
        KrishnamurtiVp291,
        LahiriIcrc,
    ]
};

/// Display names (swephprg.htm §12.1) and snake_case identifiers.
const NAMES: [(&str, &str); 47] = [
    ("Fagan/Bradley", "fagan_bradley"),
    ("Lahiri", "lahiri"),
    ("De Luce", "deluce"),
    ("Raman", "raman"),
    ("Usha/Shashi", "ushashashi"),
    ("Krishnamurti", "krishnamurti"),
    ("Djwhal Khul", "djwhal_khul"),
    ("Yukteshwar", "yukteshwar"),
    ("J.N. Bhasin", "jn_bhasin"),
    ("Babylonian/Kugler 1", "babyl_kugler1"),
    ("Babylonian/Kugler 2", "babyl_kugler2"),
    ("Babylonian/Kugler 3", "babyl_kugler3"),
    ("Babylonian/Huber", "babyl_huber"),
    ("Babylonian/Eta Piscium", "babyl_etpsc"),
    ("Babylonian/Aldebaran = 15 Tau", "aldebaran_15tau"),
    ("Hipparchos", "hipparchos"),
    ("Sassanian", "sassanian"),
    ("Galact. Center = 0 Sag", "galcent_0sag"),
    ("J2000", "j2000"),
    ("J1900", "j1900"),
    ("B1950", "b1950"),
    ("Suryasiddhanta", "suryasiddhanta"),
    ("Suryasiddhanta, mean Sun", "suryasiddhanta_msun"),
    ("Aryabhata", "aryabhata"),
    ("Aryabhata, mean Sun", "aryabhata_msun"),
    ("SS Revati", "ss_revati"),
    ("SS Citra", "ss_citra"),
    ("True Citra", "true_citra"),
    ("True Revati", "true_revati"),
    ("True Pushya (PVRN Rao)", "true_pushya"),
    ("Galactic Center (Gil Brand)", "galcent_rgilbrand"),
    ("Galactic Equator (IAU1958)", "galequ_iau1958"),
    ("Galactic Equator", "galequ_true"),
    ("Galactic Equator mid-Mula", "galequ_mula"),
    ("Skydram (Mardyks)", "galalign_mardyks"),
    ("True Mula (Chandra Hari)", "true_mula"),
    ("Dhruva/Gal.Center/Mula (Wilhelm)", "galcent_mula_wilhelm"),
    ("Aryabhata 522", "aryabhata_522"),
    ("Babylonian/Britton", "babyl_britton"),
    ("\"Vedic\"/Sheoran", "true_sheoran"),
    ("Cochrane (Gal.Center = 0 Cap)", "galcent_cochrane"),
    ("Galactic Equator (Fiorenza)", "galequ_fiorenza"),
    ("Vettius Valens", "valens_moon"),
    ("Lahiri 1940", "lahiri_1940"),
    ("Lahiri VP285", "lahiri_vp285"),
    ("Krishnamurti-Senthilathiban", "krishnamurti_vp291"),
    ("Lahiri ICRC", "lahiri_icrc"),
];

impl Ayanamsa {
    /// All predefined modes, in Swiss order (index = Swiss number).
    pub fn all() -> &'static [Ayanamsa; 47] {
        &ALL_MODES
    }

    /// Mode for a Swiss `SE_SIDM_*` number (0–46). `SE_SIDM_USER` needs
    /// parameters: build `Ayanamsa::User` directly.
    pub fn from_swiss(n: i32) -> Option<Ayanamsa> {
        usize::try_from(n).ok().and_then(|i| ALL_MODES.get(i).copied())
    }

    /// Swiss `SE_SIDM_*` number (255 for `User`).
    pub fn swiss_number(self) -> i32 {
        match self {
            Ayanamsa::Pushya => 29,
            Ayanamsa::User { .. } => SIDM_USER,
            m => ALL_MODES.iter().position(|&x| x == m).unwrap() as i32,
        }
    }

    /// Display name (Swiss `get_ayanamsa_name`); `User` gives "".
    pub fn name(self) -> &'static str {
        match self {
            Ayanamsa::User { .. } => "",
            m => NAMES[m.swiss_number() as usize].0,
        }
    }

    /// snake_case identifier (the Swiss constant name without `SIDM_`).
    pub fn ident(self) -> &'static str {
        match self {
            Ayanamsa::User { .. } => "user",
            m => NAMES[m.swiss_number() as usize].1,
        }
    }

    /// Parse a Swiss number ("1"), a constant name ("SIDM_LAHIRI"), an
    /// identifier ("lahiri", "true-citra"), a display name ("Lahiri ICRC") or
    /// a common alias ("kp", "fagan", "pushya").
    pub fn parse(s: &str) -> Option<Ayanamsa> {
        let t = s.trim();
        if let Ok(n) = t.parse::<i32>() {
            return Ayanamsa::from_swiss(n);
        }
        let squash = |x: &str| {
            x.chars().filter(|c| c.is_ascii_alphanumeric()).collect::<String>().to_ascii_lowercase()
        };
        let low = t.to_ascii_lowercase();
        let low = low.strip_prefix("se_").unwrap_or(&low);
        let ks = squash(low.strip_prefix("sidm_").unwrap_or(low));
        let alias = match ks.as_str() {
            "kp" => Some(Ayanamsa::Krishnamurti),
            "fagan" => Some(Ayanamsa::FaganBradley),
            "pushya" => Some(Ayanamsa::TruePushya),
            "galcentrgbrand" => Some(Ayanamsa::GalCentRgBrand),
            "chitrapaksha" | "truechitra" | "truechitrapaksha" => Some(Ayanamsa::TrueCitra),
            "usha" => Some(Ayanamsa::UshaShashi),
            "skydram" | "mardyks" => Some(Ayanamsa::GalAlignMardyks),
            "vedic" | "sheoran" => Some(Ayanamsa::TrueSheoran),
            "kpsenthilathiban" => Some(Ayanamsa::KrishnamurtiVp291),
            _ => None,
        };
        alias.or_else(|| {
            NAMES.iter().position(|(n, id)| squash(id) == ks || squash(n) == ks).map(|i| ALL_MODES[i])
        })
    }

    /// The Swiss modes whose sidereal positions are, by default, projected
    /// onto the ecliptic of t0 (swephprg.htm §3.3 note on equatorial
    /// coordinates; confirmed black-box).
    pub fn projects_by_default(self) -> bool {
        matches!(self, Ayanamsa::J2000 | Ayanamsa::J1900 | Ayanamsa::B1950 | Ayanamsa::GalAlignMardyks)
    }
}

/// What a mode is anchored to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Kind {
    /// Fixed epoch: (t0 JD TT, mean A0 deg).
    Fixed(f64, f64),
    /// Catalogue star (name for `stars::lookup`) at sidereal longitude (deg).
    Star(&'static str, f64),
    /// δ Cancri, built-in Hipparcos solution (True Pushya), at this longitude.
    DeltaCnc(f64),
    /// Galactic centre (Sgr A*) at this sidereal longitude.
    GalCenter(f64),
    /// Galactic centre, polar projection, at this sidereal longitude.
    GalCenterPolar(f64),
    /// Galactic-equator node at this sidereal longitude.
    GalNode(galactic::GalacticPole, f64),
    /// A0 is the geometric mean-of-date longitude of a star at t0 minus
    /// this longitude (Aldebaran at 15° Taurus): (t0 TT, star, lon).
    StarAtEpoch(f64, &'static str, f64),
    /// No published definition available to us (see the report).
    Unsupported,
}

const AS: f64 = 1.0 / 3600.0;
/// Golden section between 0° Scorpio and 0° Aquarius (R. Gil Brand):
/// 300° − 90°/φ = 210° + 90°/φ².
#[allow(clippy::approx_constant)] // φ written out: f64::sqrt is not const
const GIL_BRAND_GC_LON: f64 = 300.0 - 90.0 / 1.618_033_988_749_895;

/// Δψ at the Lahiri epoch, degrees, with the IAU 2000B nutation model
/// (McCarthy & Luzum 2003; `frames::nut00b`) — the model used by the oracle
/// to separate true from mean ayanamsa. Checked by a unit test.
const DPSI_LAHIRI_T0_DEG: f64 = 0.004660197838803416;

/// Precession-model reconversion offsets (arcsec) published in the Swiss
/// Ephemeris documentation (swisseph.htm §2.8.11: the value stated by each
/// author refers to the precession model of their time — Newcomb for
/// Fagan–Bradley, Raman, Krishnamurti and Lahiri 1940, IAU 1976 for the IAE
/// Lahiri value — and is re-expressed for Vondrák 2011 precession by these
/// documented amounts). Taken from that documentation, not from oracle output.
const OFFSET_LAHIRI_IAU1976_AS: f64 = -0.13036;
const OFFSET_NEWCOMB_FB_AS: f64 = 0.41256;
const OFFSET_NEWCOMB_1900_AS: f64 = 0.82800;

pub const J1900: f64 = 2415020.0;
pub const B1950: f64 = 2433282.42345905;
pub const J2000: f64 = 2451545.0;

/// Apparent autumnal equinox 1998 (Sun's apparent longitude 180°), JD TT,
/// computed with this engine (JPL DE442; re-checked by
/// `tests/ayanamsa_all.rs` when the pack is present). Mardyks: ayanamsa 30°
/// at that instant.
pub const MARDYKS_T0_TT: f64 = 2451079.734892;

fn dms(d: f64, m: f64, s: f64) -> f64 {
    let sg = if d < 0.0 { -1.0 } else { 1.0 };
    sg * (d.abs() + m / 60.0 + s / 3600.0)
}

fn ut(jd_ut: f64) -> f64 {
    crate::time::ut_to_tt_with(jd_ut, &DeltaTModel::SWISS_COMPAT)
}

fn jul(y: i32, m: u32, d: u32, h: f64, mi: f64, s: f64) -> f64 {
    julday(y, m, d, h + mi / 60.0 + s / 3600.0, Calendar::Julian)
}

/// Epoch of the Suryasiddhanta ayanamshas (Kali 3600, Suryasiddhanta year):
/// 499-03-21 07:30:31.57 UT (Julian calendar).
fn t0_ss() -> f64 {
    ut(jul(499, 3, 21, 7.0, 30.0, 31.57))
}

fn t0_m100() -> f64 {
    ut(jul(-100, 1, 1, 0.0, 0.0, 0.0))
}

/// Published-definition values, with the documented interpretation
/// (true vs mean, precession-model offset). `bits` may remove the offset
/// (NO_PREC_OFFSET / PREC_ORIG).
pub fn kind(mode: Ayanamsa, bits: SidBits) -> Kind {
    use Ayanamsa::*;
    let off =
        |a: f64| if bits.has(SidBits::NO_PREC_OFFSET) || bits.has(SidBits::PREC_ORIG) { 0.0 } else { a * AS };
    match mode {
        // Lahiri (IAE 1985): 23°15'00.658" TRUE at 1956-03-21 0h TDT, IAU 1976 offset.
        Lahiri => Kind::Fixed(
            2435553.5,
            dms(23.0, 15.0, 0.658) - off(OFFSET_LAHIRI_IAU1976_AS) - DPSI_LAHIRI_T0_DEG,
        ),
        // Lahiri (ICRC 1955, pre-1985 IAE/IENA): 23°15'00" TRUE, same epoch.
        // The ICRC value refers to Newcomb precession (IENA before 1985, per the
        // Swiss documentation, App. E); its reconversion offset is not
        // tabulated, so it is computed from the published models
        // (`orig_precession::reconversion_offset_arcsec`, Kinoshita 1975).
        LahiriIcrc => Kind::Fixed(
            2435553.5,
            dms(23.0, 15.0, 0.0)
                - off(orig_precession::reconversion_offset_arcsec(
                    orig_precession::OrigModel::Newcomb,
                    2435553.5,
                ))
                - DPSI_LAHIRI_T0_DEG,
        ),
        Raman => Kind::Fixed(self::J1900, 21.01444 - off(OFFSET_NEWCOMB_1900_AS)),
        Krishnamurti => Kind::Fixed(self::J1900, 22.363889 - off(OFFSET_NEWCOMB_1900_AS)),
        FaganBradley => Kind::Fixed(self::B1950, dms(24.0, 2.0, 31.36) - off(OFFSET_NEWCOMB_FB_AS)),
        Yukteshwar => Kind::Fixed(self::J1900, 21.082222),
        // Lahiri, Panchanga Darpan (1940): mean 22°26'45.50" at J1900.0, Newcomb offset.
        Lahiri1940 => Kind::Fixed(self::J1900, dms(22.0, 26.0, 45.5) - off(OFFSET_NEWCOMB_1900_AS)),
        UshaShashi => Kind::Fixed(self::J1900, dms(18.0, 39.0, 39.46)),
        // De Luce: zero at the birth of Jesus, 1 Jan 1 BCE (year 0) 0h UT.
        DeLuce => Kind::Fixed(ut(jul(0, 1, 1, 0.0, 0.0, 0.0)), 0.0),
        // Britton (2010): −3.2° at 1 Jan year 0, 0h UT.
        BabylBritton => Kind::Fixed(ut(jul(0, 1, 1, 0.0, 0.0, 0.0)), -3.2),
        // Djwhal Khul: Age of Aquarius (ayanamsa 30°) on 1 July 2117 (G. Dawson).
        DjwhalKhul => Kind::Fixed(ut(julday(2117, 7, 1, 0.0, Calendar::Gregorian)), 30.0),
        // Kugler 1/2/3, Huber: values at the year −100 (1 Jan −100, 0h UT).
        BabylKugler1 => Kind::Fixed(t0_m100(), dms(-5.0, 40.0, 0.0)),
        BabylKugler2 => Kind::Fixed(t0_m100(), dms(-4.0, 16.0, 0.0)),
        BabylKugler3 => Kind::Fixed(t0_m100(), dms(-3.0, 25.0, 0.0)),
        BabylHuber => Kind::Fixed(t0_m100(), dms(-4.0, 28.0, 0.0)),
        // Mercier: eta Piscium culminating with the zero point, −5°04'46" at −129.
        BabylEtaPiscium => Kind::Fixed(ut(jul(-129, 1, 1, 0.0, 0.0, 0.0)), dms(-5.0, 4.0, 46.0)),
        // Aldebaran at 15°00'00" Taurus around −100 (geometric position).
        BabylAldebaran15Tau => Kind::StarAtEpoch(t0_m100(), "Aldebaran", 45.0),
        // Mercier: −9°20' on 27 June −128 (JD 1674484).
        Hipparchos => Kind::Fixed(1674484.0, dms(-9.0, 20.0, 0.0)),
        // Sassanian: zero on 18 Mar 564 (JD 1927135.8747793, read as UT).
        Sassanian => Kind::Fixed(ut(1927135.8747793), 0.0),
        J2000 => Kind::Fixed(self::J2000, 0.0),
        J1900 => Kind::Fixed(self::J1900, 0.0),
        B1950 => Kind::Fixed(self::B1950, 0.0),
        Suryasiddhanta => Kind::Fixed(t0_ss(), 0.0),
        SuryasiddhantaMeanSun => Kind::Fixed(t0_ss(), -0.21463395),
        SsRevati => Kind::Fixed(t0_ss(), -0.79167046),
        SsCitra => Kind::Fixed(t0_ss(), 2.11070444),
        // Aryabhata, Kali 3600: 499-03-21 06:56:55.57 UT.
        Aryabhata => Kind::Fixed(ut(jul(499, 3, 21, 6.0, 56.0, 55.57)), 0.0),
        AryabhataMeanSun => Kind::Fixed(ut(jul(499, 3, 21, 6.0, 56.0, 55.57)), -0.23763238),
        // Govindasvamin: zero on 21 March 522, 05:46:44 UT.
        Aryabhata522 => Kind::Fixed(ut(jul(522, 3, 21, 5.0, 46.0, 44.0)), 0.0),
        // Lahiri, Tables of the Sun (1972): mean equinox 22 March 285, 17:54:02 TT.
        LahiriVp285 => Kind::Fixed(jul(285, 3, 22, 17.0, 54.0, 2.0), 0.0),
        // Krishnamurti/Senthilathiban: equinox 21 March 291, 06:10:29 TT.
        KrishnamurtiVp291 => Kind::Fixed(jul(291, 3, 21, 6.0, 10.0, 29.0), 0.0),
        // Fiorenza: exactly 25° on 1 Jan 2000 (0h UT).
        GalEquFiorenza => Kind::Fixed(ut(julday(2000, 1, 1, 0.0, Calendar::Gregorian)), 25.0),
        // Mardyks: 30° on the autumnal equinox 1998.
        GalAlignMardyks => Kind::Fixed(MARDYKS_T0_TT, 30.0),
        TrueCitra => Kind::Star("Spica", 180.0),
        TrueRevati => Kind::Star("Revati", dms(359.0, 50.0, 0.0)),
        TrueMula => Kind::Star("Shaula", 240.0),
        TruePushya | Pushya => Kind::DeltaCnc(TRUE_PUSHYA_STAR_LON_DEG),
        GalCent0Sag => Kind::GalCenter(240.0),
        GalCentCochrane => Kind::GalCenter(270.0),
        GalCentRgBrand => Kind::GalCenter(GIL_BRAND_GC_LON),
        GalCentMulaWilhelm => Kind::GalCenterPolar(dms(246.0, 40.0, 0.0)),
        GalEquIau1958 => Kind::GalNode(galactic::GalacticPole::Iau1958, 240.0),
        GalEquTrue => Kind::GalNode(galactic::GalacticPole::Modern, 240.0),
        GalEquMula => Kind::GalNode(galactic::GalacticPole::Modern, dms(246.0, 40.0, 0.0)),
        User { t0, ayan_t0, t0_is_ut } => {
            Kind::Fixed(if t0_is_ut || bits.has(SidBits::USER_UT) { ut(t0) } else { t0 }, ayan_t0)
        }
        JnBhasin | TrueSheoran | ValensMoon => Kind::Unsupported,
    }
}

/// Original precession model of a definition (for `SIDBIT_PREC_ORIG`).
pub fn original_model(mode: Ayanamsa) -> Option<orig_precession::OrigModel> {
    use orig_precession::OrigModel::*;
    match mode {
        Ayanamsa::Lahiri => Some(Iau1976),
        Ayanamsa::FaganBradley
        | Ayanamsa::Raman
        | Ayanamsa::Krishnamurti
        | Ayanamsa::Lahiri1940
        | Ayanamsa::LahiriIcrc => Some(Newcomb),
        _ => None,
    }
}

/// Reference epoch (JD TT) and MEAN ayanamsa at that epoch (degrees) of a
/// fixed-epoch mode, with default bits.
pub fn reference(mode: Ayanamsa) -> Option<(f64, f64)> {
    match kind(mode, SidBits::default()) {
        Kind::Fixed(t0, a0) => Some((t0, a0)),
        _ => None,
    }
}

/// Sidereal longitude of δ Cancri in the True-Pushya definition
/// (P. V. R. Narasimha Rao): 16°00' Cancer = 106°.
pub const TRUE_PUSHYA_STAR_LON_DEG: f64 = 106.0;

/// Mean fixed-epoch ayanamsa at `jd_tt` with option bits.
pub fn fixed_mean_deg(mode: Ayanamsa, t0: f64, a0: f64, jd_tt: f64, bits: SidBits) -> f64 {
    if bits.has(SidBits::PREC_ORIG) {
        if let Some(m) = original_model(mode) {
            return norm_deg(a0 - orig_precession::equinox_lon_on_ecliptic_of(m, jd_tt, t0) / DEG);
        }
    }
    if bits.has(SidBits::ECL_DATE) && !mode.projects_by_default() {
        // Zero point (longitude A0 on the ecliptic of t0) re-measured on the
        // mean ecliptic of date.
        let m0 = vondrak::ecliptic_of_date_matrix(t0);
        let (s, c) = (a0 * DEG).sin_cos();
        let z = crate::math::mat_t_vec(&m0, [c, s, 0.0]);
        let w = crate::math::mat_vec(&vondrak::ecliptic_of_date_matrix(jd_tt), z);
        return norm_deg(w[1].atan2(w[0]) / DEG);
    }
    norm_deg(a0 - vondrak::equinox_lon_on_ecliptic_of(jd_tt, t0) / DEG)
}

/// Mean ayanamsa at JD TT, degrees, in [0, 360), using only built-in models
/// (no ephemeris pack). Star-anchored modes use the built-in low-precision
/// Earth for aberration/parallax (≲ 5e-6°); the engine
/// (`Engine::ayanamsa_deg`) uses JPL states instead. Returns NaN for modes
/// without an available definition (`Kind::Unsupported`).
pub fn ayanamsa_mean_deg(mode: Ayanamsa, jd_tt: f64) -> f64 {
    ayanamsa_mean_deg_bits(mode, jd_tt, SidBits::default())
}

/// [`ayanamsa_mean_deg`] with explicit sidereal option bits.
pub fn ayanamsa_mean_deg_bits(mode: Ayanamsa, jd_tt: f64, bits: SidBits) -> f64 {
    match kind(mode, bits) {
        Kind::Fixed(t0, a0) => fixed_mean_deg(mode, t0, a0, jd_tt, bits),
        Kind::StarAtEpoch(t0, name, lon) => {
            fixed_mean_deg(mode, t0, star_at_epoch_a0(t0, name, lon), jd_tt, bits)
        }
        Kind::DeltaCnc(lon) => norm_deg(star::delta_cnc_apparent_mean(jd_tt).0 - lon),
        Kind::Star(name, lon) => {
            let rec = crate::stars::lookup(name).expect("catalogue star").record;
            norm_deg(star::apparent_mean_builtin(jd_tt, crate::stars::barycentric_au(rec, jd_tt)).0 - lon)
        }
        Kind::GalCenter(lon) => {
            norm_deg(star::apparent_mean_builtin(jd_tt, galactic::sgr_a_bary_au(jd_tt)).0 - lon)
        }
        Kind::GalCenterPolar(lon) => {
            let (l, b) = star::apparent_mean_builtin(jd_tt, galactic::sgr_a_bary_au(jd_tt));
            let (sl, cl) = (l * DEG).sin_cos();
            let (sb, cb) = (b * DEG).sin_cos();
            let p = galactic::mean_pole_in_ecliptic(jd_tt);
            norm_deg(galactic::polar_longitude_deg([cb * cl, cb * sl, sb], p) - lon)
        }
        Kind::GalNode(p, lon) => norm_deg(galactic::galactic_node_lon_deg(p, jd_tt) - lon),
        Kind::Unsupported => f64::NAN,
    }
}

/// A0 (deg) of a star-at-epoch definition: geometric geocentric
/// mean-of-date longitude of the star at t0 minus `lon`.
pub fn star_at_epoch_a0(t0: f64, name: &str, lon: f64) -> f64 {
    let rec = crate::stars::lookup(name).expect("catalogue star").record;
    star::geometric_mean_builtin(t0, crate::stars::barycentric_au(rec, t0)).0 - lon
}

/// Rate of the mean ayanamsa, degrees/day (central difference, h = 1 day).
pub fn ayanamsa_mean_rate_deg_per_day(mode: Ayanamsa, jd_tt: f64) -> f64 {
    let a = ayanamsa_mean_deg(mode, jd_tt - 1.0);
    let b = ayanamsa_mean_deg(mode, jd_tt + 1.0);
    crate::math::diff_deg(b, a) / 2.0
}

/// True Pushya ayanamsa from the ecliptic longitude of δ Cancri, degrees.
///
/// Swiss behaviour (measured): pass the star's *apparent* longitude
/// (aberration + light deflection) referred to the mean equinox of date to
/// get the mean ayanamsa, or referred to the true equinox (with nutation) to
/// get the true ayanamsa; the relation is exactly `lon - 106°`.
pub fn true_pushya_deg(delta_cnc_lon_deg: f64) -> f64 {
    norm_deg(delta_cnc_lon_deg - TRUE_PUSHYA_STAR_LON_DEG)
}

/// Invariable plane of the solar system relative to the mean ecliptic and
/// equinox J2000: inclination 1.578694°, ascending node 107.582222°
/// (Souami & Souchay 2012, A&A 543, A133).
pub const SSY_INCLINATION_DEG: f64 = 1.578694;
pub const SSY_NODE_DEG: f64 = 107.582222;

/// Rotation J2000 mean equator (Vondrák frame, no bias) -> invariable-plane
/// frame with x along the plane's ascending node on the J2000 ecliptic.
pub fn ssy_matrix() -> crate::math::Mat3 {
    use crate::math::{mat_mul, rot_x, rot_z};
    let ecl = vondrak::ecliptic_of_date_matrix(J2000);
    mat_mul(&rot_x(SSY_INCLINATION_DEG * DEG), &mat_mul(&rot_z(SSY_NODE_DEG * DEG), &ecl))
}

/// How sidereal positions are formed for a mode and option bits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Projection {
    /// sidereal = tropical (same flags) − ayanamsa.
    Traditional,
    /// Project onto the mean ecliptic of t0: lon = λ_{t0} − A0, lat = β_{t0}.
    EclT0 { t0: f64, a0: f64 },
    /// Project onto the invariable plane; lon = λ_ssy − λ_ssy(equinox t0) − A0.
    Ssy { t0: f64, a0: f64 },
}

/// The projection rule (fixed-epoch modes only; star-anchored and galactic
/// modes ignore the projection bits, as the oracle does).
pub fn projection(mode: Ayanamsa, bits: SidBits) -> Projection {
    let (t0, a0) = match kind(mode, bits) {
        Kind::Fixed(t0, a0) => (t0, a0),
        Kind::StarAtEpoch(t0, n, l) => (t0, star_at_epoch_a0(t0, n, l)),
        _ => return Projection::Traditional,
    };
    if bits.has(SidBits::SSY_PLANE) && !mode.projects_by_default() {
        Projection::Ssy { t0, a0 }
    } else if bits.has(SidBits::ECL_T0) || mode.projects_by_default() {
        Projection::EclT0 { t0, a0 }
    } else {
        Projection::Traditional
    }
}

impl Projection {
    /// Sidereal (lon, lat) in degrees of a direction `v` given in the J2000
    /// mean equatorial frame (Vondrák frame, no bias). `None` = traditional.
    pub fn apply(&self, v: crate::math::Vec3) -> Option<(f64, f64)> {
        use crate::math::{mat_vec, norm, Vec3};
        let sph =
            |w: Vec3| (norm_deg(w[1].atan2(w[0]) / DEG), (w[2] / norm(w)).clamp(-1.0, 1.0).asin() / DEG);
        match *self {
            Projection::Traditional => None,
            Projection::EclT0 { t0, a0 } => {
                let (l, b) = sph(mat_vec(&vondrak::ecliptic_of_date_matrix(t0), v));
                Some((norm_deg(l - a0), b))
            }
            Projection::Ssy { t0, a0 } => {
                let m = ssy_matrix();
                let (l, b) = sph(mat_vec(&m, v));
                let (z, _) = sph(mat_vec(&m, vondrak::ecliptic_of_date_matrix(t0)[0]));
                Some((norm_deg(l - z - a0), b))
            }
        }
    }

    /// Mean equator of t0 (for RA/Dec output with `EclT0`), rows in the J2000
    /// mean equatorial frame.
    pub fn equator_of_t0(&self) -> Option<crate::math::Mat3> {
        use crate::math::{cross, norm, scale};
        match *self {
            Projection::EclT0 { t0, .. } => {
                let z = vondrak::equator_pole(t0);
                let e = cross(z, vondrak::ecliptic_pole(t0));
                let x = scale(e, 1.0 / norm(e));
                Some([x, cross(z, x), z])
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lahiri_at_epoch() {
        let v = ayanamsa_mean_deg(Ayanamsa::Lahiri, 2435553.5);
        let dpsi = crate::frames::nut00b(2435553.5).0 / DEG;
        assert!((dpsi - DPSI_LAHIRI_T0_DEG).abs() < 1e-12);
        let exp = 23.0 + 15.0 / 60.0 + (0.658 + 0.13036) / 3600.0 - dpsi;
        assert!((v - exp).abs() < 1e-12);
    }

    #[test]
    fn numbering_and_names() {
        for (i, m) in Ayanamsa::all().iter().enumerate() {
            assert_eq!(m.swiss_number(), i as i32);
            assert_eq!(Ayanamsa::from_swiss(i as i32), Some(*m));
            assert_eq!(Ayanamsa::parse(m.ident()), Some(*m), "{}", m.ident());
            assert_eq!(Ayanamsa::parse(m.name()), Some(*m), "{}", m.name());
            assert_eq!(Ayanamsa::parse(&i.to_string()), Some(*m));
        }
        assert_eq!(Ayanamsa::parse("KP"), Some(Ayanamsa::Krishnamurti));
        assert_eq!(Ayanamsa::parse("SIDM_LAHIRI_ICRC"), Some(Ayanamsa::LahiriIcrc));
        assert_eq!(Ayanamsa::parse("SE_SIDM_J2000"), Some(Ayanamsa::J2000));
        assert_eq!(Ayanamsa::parse("true-pushya"), Some(Ayanamsa::TruePushya));
        assert_eq!(Ayanamsa::parse("fagan"), Some(Ayanamsa::FaganBradley));
        assert_eq!(Ayanamsa::Pushya.swiss_number(), 29);
        assert_eq!(Ayanamsa::from_swiss(47), None);
        assert_eq!(Ayanamsa::Lahiri.name(), "Lahiri");
    }

    #[test]
    fn eclt0_equinox_of_t0() {
        let p = Projection::EclT0 { t0: J1900, a0: 0.0 };
        let x = vondrak::ecliptic_of_date_matrix(J1900)[0];
        let (l, b) = p.apply(x).unwrap();
        assert!(crate::math::diff_deg(l, 0.0).abs() < 1e-12 && b.abs() < 1e-12);
    }
}
