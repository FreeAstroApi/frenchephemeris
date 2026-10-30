// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Events and observational utilities: rise / set / meridian transits,
//! planetary phenomena (phase, elongation, diameter, magnitude), Sun/Moon
//! longitude crossings, equation of time, ecliptic↔equatorial rotation,
//! horizontal coordinates and refraction.
//!
//! Semantics follow the classic astrology-engine API that FreeAstroAPI calls
//! (`rise_trans`, `pheno_ut`, `solcross_ut`, `time_equ`, `cotrans`), measured
//! black-box against pyswisseph 2.10 (public API only) and implemented from
//! published algorithms. See docs/EVENTS_REPORT.md for the measurements.
//!
//! Rise/set model (semantics observed black-box, see the report):
//! - the body position is the *topocentric apparent* place of date (true
//!   equinox) for the observer, whatever ephemeris flags the caller passes;
//! - the event is the instant the chosen disc point (upper limb by default)
//!   reaches the *true* altitude `h0 = horizon − R(horizon)`, where `R` is a
//!   published refraction formula (default Sinclair 1982 / Explanatory
//!   Supplement 1992 §3.28, see
//!   `HorizonRefraction`) for (pressure, temperature); pressure 0 means
//!   "standard atmosphere at the observer height"; the dip of the horizon is
//!   not applied. No numerical value in this module is fitted to the oracle:
//!   see "Provenance" in docs/EVENTS_REPORT.md;
//! - `geocentric_no_lat` ("Hindu rising"): geocentric apparent place with the
//!   ecliptic latitude set to zero;
//! - the next event after `jd_ut` is searched inside a fixed window of 26 h;
//!   if there is none (circumpolar, polar day/night, or the Moon skipping a
//!   rise/set) the result is `AlwaysAbove` / `AlwaysBelow` (the classic API's
//!   return code −2);
//! - meridian transits (hour angle 0 / 180°) always use the topocentric
//!   apparent place; disc, refraction and "Hindu" options do not change them.
//!   (The classic API's Hindu-transit branch reads the process-global
//!   `set_topo` location instead of the given one — a thread-safety bug we do
//!   not reproduce.)
//!
//! Altitude is a smooth function of time whose extrema are (to first order)
//! the meridian transits, so the search finds the transits in the window by
//! secant iteration on the hour angle, evaluates the altitude there, and
//! refines the bracketing interval with Brent. Epochs inside the window come
//! from `EpochInterp` (three exact epochs + quadratic interpolation).

use crate::context::{Body, CalcOptions, Engine};
use crate::ephem::{self, Corrections, Epoch, Frame, LightTime, ObsState, Observer, Origin};
use crate::frames::{FrameModel, FrameSet};
use crate::math::*;
use crate::search::brent;
use crate::time::{self, DeltaTModel};
use crate::{Error, Result};

/// Search window of `rise_set`, days. A semantic choice (what "no event"
/// means), set to 26 h so the circumpolar / no-event classification is the
/// classic API's (observed: events up to +1.0795 d reported, none from
/// +1.0840 d); not a fitted numerical constant.
pub const RISE_SET_WINDOW_DAYS: f64 = 26.0 / 24.0;

/// Mean radii (km): IAU WGCCRE 2015 report (Archinal et al. 2018) for the
/// planets, Pluto and Moon; Sun 696 000 km (Allen, Astrophysical Quantities);
/// Earth WGS84.
pub fn body_radius_km(body: Body) -> f64 {
    use Body::*;
    match body {
        Sun => 696_000.0,
        Moon => 1737.4,
        Mercury => 2439.4,
        Venus => 6051.8,
        Mars => 3389.5,
        Jupiter => 69_911.0,
        Saturn => 58_232.0,
        Uranus => 25_362.0,
        Neptune => 24_622.0,
        Pluto => 1188.3,
        Earth => 6378.137,
        _ => 0.0,
    }
}

/// Mean Earth–Moon distance (km, conventional value, e.g. Explanatory
/// Supplement / AA) used for the fixed-disc-size option.
const MOON_MEAN_DIST_KM: f64 = 384_400.0;
const EARTH_EQ_RADIUS_KM: f64 = 6378.137;

// ---------------------------------------------------------------------------------------------
// Refraction
// ---------------------------------------------------------------------------------------------

/// Standard pressure (hPa) at height `h_m` above sea level: barometric formula
/// of the ICAO / US Standard Atmosphere 1976 troposphere (P0 = 1013.25 hPa,
/// T0 = 288 K, lapse 6.5 K/km, exponent g·M/(R·L) = 5.255).
pub fn standard_pressure_hpa(height_m: f64) -> f64 {
    1013.25 * (1.0 - 0.0065 * height_m / 288.0).max(0.0).powf(5.255)
}

/// Atmosphere for refraction. `pressure_hpa == 0` means "standard pressure
/// at the observer height".
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atmosphere {
    pub pressure_hpa: f64,
    pub temperature_c: f64,
}

impl Default for Atmosphere {
    /// The classic API defaults (pressure 0 → estimated, temperature 0 °C).
    fn default() -> Self {
        Atmosphere { pressure_hpa: 0.0, temperature_c: 0.0 }
    }
}

impl Atmosphere {
    pub fn pressure_at(&self, height_m: f64) -> f64 {
        if self.pressure_hpa == 0.0 {
            standard_pressure_hpa(height_m)
        } else {
            self.pressure_hpa
        }
    }
}

/// Published refraction models. All take the *apparent* altitude where noted
/// and return the refraction R (deg) = apparent − true altitude.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RefractionModel {
    /// Default. Sinclair (1982) as printed in the Explanatory Supplement to the
    /// Astronomical Almanac (1992), §3.28 (eq. 3.283-1/-2; HM Nautical Almanac
    /// Office low-altitude formula), h = apparent altitude in degrees:
    ///   R0 = (34.46 + 4.23 h + 0.004 h²) / (1 + 0.505 h + 0.0845 h²)  arcmin,
    /// above the altitude where it meets the high-altitude asymptote
    ///   R0 = 0.97′ · cot h  (58.2″ cot h, standard conditions 1010 hPa / 10 °C),
    /// both scaled to the actual atmosphere by
    ///   R = R0 · (P − 80)/930 / (1 + 8·10⁻⁵ (R0 + 39)(T − 10)).
    /// Identified (coordinator, then verified here black-box) as the model of
    /// the classic API's extended refraction (written with V. Reijs, whose
    /// ARCHAEOCOSMO pages document "the inverse of Sinclair's formula");
    /// it reproduces it to < 0.001″ at 0–90° for every P/T tested.
    #[default]
    Sinclair1982,
    /// Sæmundsson (1986) true→apparent formula, P/T scaled (Meeus ch. 16),
    /// inverted numerically for an apparent altitude.
    Saemundsson,
    /// Bennett (1982) apparent→true formula, P/T scaled (Meeus eq. 16.3).
    Bennett,
    /// Explanatory Supplement (1992) eq. 3.283-2 simplified Sinclair form
    /// P(0.1594 + 0.0196a + 0.00002a²)/((273+T)(1 + 0.505a + 0.0845a²)) (a < 15°).
    SinclairSimple,
    /// Astronomical Almanac rise/set convention: a fixed 34′ at the horizon.
    Almanac34,
}

/// Former name of `RefractionModel` (rise/set horizon refraction).
pub type HorizonRefraction = RefractionModel;

/// Bennett (1982) refraction (deg) for an *apparent* altitude `h_app` (deg).
pub fn refraction_bennett(h_app_deg: f64, pressure_hpa: f64, temperature_c: f64) -> f64 {
    let h = h_app_deg;
    let r_arcmin = 1.0 / ((h + 7.31 / (h + 4.4)) * DEG).tan();
    r_arcmin / 60.0 * (pressure_hpa / 1010.0) * (283.0 / (273.0 + temperature_c))
}

/// Sæmundsson (1986) refraction (deg) for a *true* altitude `h_true` (deg)
/// (Meeus eq. 16.5, with P/T scaling).
pub fn refraction_saemundsson(h_true_deg: f64, pressure_hpa: f64, temperature_c: f64) -> f64 {
    let h = h_true_deg;
    let r_arcmin = 1.02 / ((h + 10.3 / (h + 5.11)) * DEG).tan();
    r_arcmin / 60.0 * (pressure_hpa / 1010.0) * (283.0 / (273.0 + temperature_c))
}

/// Simplified Sinclair form of the Explanatory Supplement (1992) eq. 3.283-2,
/// apparent altitude below 15°.
pub fn refraction_sinclair_simple(h_app_deg: f64, pressure_hpa: f64, temperature_c: f64) -> f64 {
    let a = h_app_deg;
    pressure_hpa * (0.1594 + 0.0196 * a + 0.00002 * a * a)
        / ((273.0 + temperature_c) * (1.0 + 0.505 * a + 0.0845 * a * a))
}

fn sinclair_low_arcmin(h: f64) -> f64 {
    (34.46 + 4.23 * h + 0.004 * h * h) / (1.0 + 0.505 * h + 0.0845 * h * h)
}

fn sinclair_high_arcmin(h: f64) -> f64 {
    0.97 / (h * DEG).tan()
}

/// Altitude (deg) where Sinclair's low-altitude formula meets the 0.97′ cot h
/// asymptote (≈ 17.904°); the model switches there, so it is continuous.
pub fn sinclair_switch_altitude() -> f64 {
    static H: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *H.get_or_init(|| {
        let (mut lo, mut hi) = (15.0, 20.0);
        for _ in 0..80 {
            let m = 0.5 * (lo + hi);
            if sinclair_low_arcmin(m) > sinclair_high_arcmin(m) {
                hi = m;
            } else {
                lo = m;
            }
        }
        0.5 * (lo + hi)
    })
}

/// Sinclair (1982) refraction (deg) for an *apparent* altitude `h_app` (deg),
/// Explanatory Supplement (1992) §3.28 form with the (P − 80)/930 pressure
/// and 1/(1 + 8e-5 (R0 + 39)(T − 10)) temperature factors.
pub fn refraction_sinclair(h_app_deg: f64, pressure_hpa: f64, temperature_c: f64) -> f64 {
    let h = h_app_deg;
    let r0 = if h < sinclair_switch_altitude() { sinclair_low_arcmin(h) } else { sinclair_high_arcmin(h) };
    let r = r0 * (pressure_hpa - 80.0) / 930.0 / (1.0 + 8e-5 * (r0 + 39.0) * (temperature_c - 10.0));
    r / 60.0
}

/// Refraction R (deg) at the *apparent* altitude `h_app_deg` for a model.
pub fn refraction_apparent(
    model: RefractionModel,
    h_app_deg: f64,
    pressure_hpa: f64,
    temperature_c: f64,
) -> f64 {
    let (p, t) = (pressure_hpa, temperature_c);
    match model {
        RefractionModel::Sinclair1982 => refraction_sinclair(h_app_deg, p, t),
        RefractionModel::Bennett => refraction_bennett(h_app_deg, p, t),
        RefractionModel::SinclairSimple => refraction_sinclair_simple(h_app_deg, p, t),
        RefractionModel::Almanac34 => 34.0 / 60.0,
        RefractionModel::Saemundsson => {
            // Solve h + R_S(h) = h_app by fixed-point iteration (|dR/dh| ≈ 0.2).
            let mut h = h_app_deg - refraction_saemundsson(h_app_deg, p, t);
            for _ in 0..40 {
                let hn = h_app_deg - refraction_saemundsson(h, p, t);
                if (hn - h).abs() < 1e-12 {
                    h = hn;
                    break;
                }
                h = hn;
            }
            h_app_deg - h
        }
    }
}

/// Refraction R (deg) such that a true altitude `horizon − R` is seen at the
/// apparent altitude `horizon_deg`. A horizon below apparent altitude 0
/// gets none (observed semantics at sea level; no dip model).
pub fn horizon_refraction(model: RefractionModel, horizon_deg: f64, atm: &Atmosphere, height_m: f64) -> f64 {
    if horizon_deg < 0.0 {
        return 0.0;
    }
    refraction_apparent(model, horizon_deg, atm.pressure_at(height_m), atm.temperature_c)
}

/// Direction of a `refrac` conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefracDir {
    TrueToApparent,
    ApparentToTrue,
}

/// Refraction conversion (deg) with the default model (`Sinclair1982`).
pub fn refrac(alt_deg: f64, pressure_hpa: f64, temperature_c: f64, dir: RefracDir) -> f64 {
    refrac_with(RefractionModel::default(), alt_deg, pressure_hpa, temperature_c, dir)
}

/// Refraction conversion for a model. The models are functions of the
/// apparent altitude h_a ≥ 0 (no refraction is applied below the apparent
/// horizon, observed semantics for a sea-level observer):
/// apparent→true is h_a − R(h_a); true→apparent solves h_a − R(h_a) = h_t
/// numerically (no solution with h_a ≥ 0 when h_t < −R(0): returned unchanged).
pub fn refrac_with(
    model: RefractionModel,
    alt_deg: f64,
    pressure_hpa: f64,
    temperature_c: f64,
    dir: RefracDir,
) -> f64 {
    let r = |h: f64| refraction_apparent(model, h.max(1e-9), pressure_hpa, temperature_c);
    match dir {
        RefracDir::ApparentToTrue => {
            if alt_deg < 0.0 {
                alt_deg
            } else {
                alt_deg - r(alt_deg)
            }
        }
        RefracDir::TrueToApparent => {
            if alt_deg < -r(0.0) {
                return alt_deg;
            }
            // Fixed point h_a = h_t + R(h_a); |dR/dh_a| ≤ ~0.2, so it contracts.
            let mut ha = (alt_deg + r(alt_deg.max(0.0))).max(0.0);
            for _ in 0..60 {
                let hn = (alt_deg + r(ha)).max(0.0);
                if (hn - ha).abs() < 1e-13 {
                    ha = hn;
                    break;
                }
                ha = hn;
            }
            ha
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Coordinate helpers
// ---------------------------------------------------------------------------------------------

fn polar_to_cart(lon_deg: f64, lat_deg: f64, r: f64) -> Vec3 {
    let (sl, cl) = (lon_deg * DEG).sin_cos();
    let (sb, cb) = (lat_deg * DEG).sin_cos();
    [r * cb * cl, r * cb * sl, r * sb]
}

/// Rotate spherical coordinates about the x axis by `eps_deg` (the classic
/// `cotrans`): `eps < 0` converts ecliptic → equatorial, `eps > 0`
/// equatorial → ecliptic. Input/output (lon, lat, dist) in degrees.
pub fn cotrans(coord: [f64; 3], eps_deg: f64) -> [f64; 3] {
    let v = mat_vec(&rot_x(eps_deg * DEG), polar_to_cart(coord[0], coord[1], 1.0));
    let (l, b, _) = to_spherical(v);
    [l / DEG, b / DEG, coord[2]]
}

/// `cotrans` for position and daily motion (lon, lat, dist, dlon, dlat, ddist).
pub fn cotrans_sp(coord: [f64; 6], eps_deg: f64) -> [f64; 6] {
    let m = rot_x(eps_deg * DEG);
    let (l, b) = (coord[0] * DEG, coord[1] * DEG);
    let (dl, db) = (coord[3] * DEG, coord[4] * DEG);
    let (sl, cl) = l.sin_cos();
    let (sb, cb) = b.sin_cos();
    let p = [cb * cl, cb * sl, sb];
    let v = [-sb * cl * db - cb * sl * dl, -sb * sl * db + cb * cl * dl, cb * db];
    let (p2, v2) = (mat_vec(&m, p), mat_vec(&m, v));
    let (l2, b2, _) = to_spherical(p2);
    let (dl2, db2, _) = spherical_rates(p2, v2);
    [l2 / DEG, b2 / DEG, coord[2], dl2 / DEG, db2 / DEG, coord[5]]
}

/// True altitude (rad) and azimuth (rad, from south, positive westwards) of an
/// equatorial-of-date direction at local apparent sidereal time `last` (rad)
/// and geodetic latitude `lat` (rad).
fn horizontal(ra: f64, dec: f64, last: f64, lat: f64) -> (f64, f64) {
    let h = last - ra;
    let (sh, ch) = h.sin_cos();
    let (sd, cd) = dec.sin_cos();
    let (sp, cp) = lat.sin_cos();
    let alt = (sp * sd + cp * cd * ch).clamp(-1.0, 1.0).asin();
    let az = norm_rad((cd * sh).atan2(ch * cd * sp - sd * cp));
    (alt, az)
}

/// Input frame of `azalt`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AzaltInput {
    /// Ecliptic of date (true equinox) longitude/latitude.
    Ecliptic,
    /// Equator of date (true equinox) RA/Dec.
    Equatorial,
}

/// Rise/set/transit selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RiseSetKind {
    Rise,
    Set,
    /// Upper meridian transit (culmination).
    UpperTransit,
    /// Lower meridian transit.
    LowerTransit,
}

/// Point of the disc that defines rise/set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscPoint {
    UpperLimb,
    Center,
    LowerLimb,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RiseSetOptions {
    pub disc: DiscPoint,
    /// `None` = no refraction (geometric horizon).
    pub refraction: Option<Atmosphere>,
    pub refraction_model: HorizonRefraction,
    /// Apparent altitude of the local horizon, degrees.
    pub horizon_deg: f64,
    /// Geocentric place with ecliptic latitude zeroed ("Hindu rising").
    pub geocentric_no_lat: bool,
    /// Semidiameter at mean distance (Sun 1 AU, Moon 384 400 km).
    pub fixed_disc_size: bool,
    pub window_days: f64,
}

impl Default for RiseSetOptions {
    fn default() -> Self {
        RiseSetOptions {
            disc: DiscPoint::UpperLimb,
            refraction: Some(Atmosphere::default()),
            refraction_model: RefractionModel::Sinclair1982,
            horizon_deg: 0.0,
            geocentric_no_lat: false,
            fixed_disc_size: false,
            window_days: RISE_SET_WINDOW_DAYS,
        }
    }
}

impl RiseSetOptions {
    /// The classic "Hindu rising" combination: disc centre, no refraction,
    /// geocentric without ecliptic latitude.
    pub fn hindu() -> Self {
        RiseSetOptions {
            disc: DiscPoint::Center,
            refraction: None,
            geocentric_no_lat: true,
            ..Default::default()
        }
    }
}

/// Outcome of a rise/set/transit search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RiseSet {
    /// Event at this JD (UT).
    Event(f64),
    /// No rise/set inside the window: the disc point stays above (circumpolar,
    /// polar day) ...
    AlwaysAbove,
    /// ... or below (never rises, polar night).
    AlwaysBelow,
}

impl RiseSet {
    pub fn jd(self) -> Option<f64> {
        match self {
            RiseSet::Event(t) => Some(t),
            _ => None,
        }
    }
}

/// Apparent place of a target for rise/set: (RA, Dec) of date (true equinox)
/// in rad and semidiameter in rad, at an epoch (use `ep.frames` to precess a
/// catalogue place). Implemented by solar-system bodies here and by fixed
/// stars elsewhere (semidiameter 0).
pub type TargetFn<'a> = dyn Fn(&Epoch) -> Result<(f64, f64, f64)> + 'a;

/// An `Epoch` over a ~1.5-day search window, from three exact epochs and
/// quadratic (Lagrange) interpolation of the time-scale offsets and every
/// frame matrix/angle. The cubic remainder of the fastest nutation term
/// (0.23″, 13.66 d) over the 0.75-d node spacing is < 0.001″, while one exact
/// epoch (Delta T, precession, nutation, CIO locator) costs more than the
/// whole topocentric place, so root finding evaluates ~5× faster.
pub struct EpochInterp {
    t: [f64; 3],
    e: [Epoch; 3],
    lo: f64,
    hi: f64,
    model: FrameModel,
    delta_t: DeltaTModel,
}

/// Exact epoch at a UT instant under explicit frame and Delta T models.
fn epoch_ut(jd_ut: f64, model: FrameModel, delta_t: &DeltaTModel) -> Epoch {
    Epoch::with_scales(jd_ut, time::ut_to_tt_with(jd_ut, delta_t), model)
}

impl EpochInterp {
    /// Interpolator with the default Delta T model.
    pub fn new(t0_ut: f64, model: FrameModel) -> EpochInterp {
        EpochInterp::new_with(t0_ut, model, DeltaTModel::DEFAULT)
    }

    /// Interpolator with an explicit Delta T model (an engine's).
    pub fn new_with(t0_ut: f64, model: FrameModel, delta_t: DeltaTModel) -> EpochInterp {
        let t = [t0_ut - 0.1, t0_ut + 0.65, t0_ut + 1.4];
        EpochInterp {
            t,
            e: t.map(|x| epoch_ut(x, model, &delta_t)),
            lo: t0_ut - 0.2,
            hi: t0_ut + 1.5,
            model,
            delta_t,
        }
    }

    pub fn at(&self, jd_ut: f64) -> Epoch {
        if jd_ut < self.lo || jd_ut > self.hi {
            return epoch_ut(jd_ut, self.model, &self.delta_t);
        }
        let [t0, t1, t2] = self.t;
        let w = [
            (jd_ut - t1) * (jd_ut - t2) / ((t0 - t1) * (t0 - t2)),
            (jd_ut - t0) * (jd_ut - t2) / ((t1 - t0) * (t1 - t2)),
            (jd_ut - t0) * (jd_ut - t1) / ((t2 - t0) * (t2 - t1)),
        ];
        let [a, b, c] = &self.e;
        let mix = |fa: f64, fb: f64, fc: f64| w[0] * fa + w[1] * fb + w[2] * fc;
        let mixm = |ma: &Mat3, mb: &Mat3, mc: &Mat3| {
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] = mix(ma[i][j], mb[i][j], mc[i][j]);
                }
            }
            m
        };
        // Angles are unwrapped around the middle node.
        let mixa = |fa: f64, fb: f64, fc: f64| fb + mix(wrap_pi(fa - fb), 0.0, wrap_pi(fc - fb));
        let (fa, fb, fc) = (&a.frames, &b.frames, &c.frames);
        let jd_tt = jd_ut + mix(a.jd_tt - a.jd_ut, b.jd_tt - b.jd_ut, c.jd_tt - c.jd_ut);
        let jd_tdb = jd_tt + mix(a.jd_tdb - a.jd_tt, b.jd_tdb - b.jd_tt, c.jd_tdb - c.jd_tt);
        let frames = FrameSet {
            jd_tt,
            eps_mean: mix(fa.eps_mean, fb.eps_mean, fc.eps_mean),
            dpsi: mix(fa.dpsi, fb.dpsi, fc.dpsi),
            deps: mix(fa.deps, fb.deps, fc.deps),
            eps_true: mix(fa.eps_true, fb.eps_true, fc.eps_true),
            bp: mixm(&fa.bp, &fb.bp, &fc.bp),
            bpn: mixm(&fa.bpn, &fb.bpn, &fc.bpn),
            model: self.model,
            ecl_mean: mixm(&fa.ecl_mean, &fb.ecl_mean, &fc.ecl_mean),
            ecl_true: mixm(&fa.ecl_true, &fb.ecl_true, &fc.ecl_true),
            eo: mixa(fa.eo, fb.eo, fc.eo),
            eqeq: mixa(fa.eqeq, fb.eqeq, fc.eqeq),
        };
        Epoch { jd_ut, jd_tt, jd_tdb, frames }
    }
}

/// Planetary phenomena, the classic `pheno` layout.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Pheno {
    /// Phase angle Sun–body–observer, deg.
    pub phase_angle: f64,
    /// Illuminated fraction of the disc, (1 + cos i)/2.
    pub phase: f64,
    /// Elongation from the Sun, deg.
    pub elongation: f64,
    /// Apparent diameter, deg.
    pub diameter: f64,
    /// Apparent visual magnitude (see `magnitude`).
    pub magnitude: f64,
    /// Horizontal parallax (Moon only), deg.
    pub horizontal_parallax: f64,
}

impl Pheno {
    pub fn to_array(self) -> [f64; 20] {
        let mut a = [0.0; 20];
        a[..6].copy_from_slice(&[
            self.phase_angle,
            self.phase,
            self.elongation,
            self.diameter,
            self.magnitude,
            self.horizontal_parallax,
        ]);
        a
    }
}

/// Unit vector of an ICRF pole (RA, Dec in degrees).
fn pole(ra_deg: f64, dec_deg: f64) -> Vec3 {
    polar_to_cart(ra_deg, dec_deg, 1.0)
}

/// Apparent V magnitude.
///
/// Planets Mercury…Neptune: Mallama & Hilton (2018), "Computing apparent
/// planetary magnitudes for The Astronomical Almanac", Astronomy and Computing
/// 25, 10–24 (the AA 2018+ model). Saturn includes the rings (sub-observer
/// ring latitude β from the IAU 2015 pole) where the paper's ring formula is
/// valid (α ≤ 6.5°, β ≤ 27°), the globe-only formula otherwise; Uranus includes
/// the sub-observer latitude term; Mars omits the optional orientation terms.
/// Pluto: H = −1.01, no phase term (AA). Moon: Allen (1976) as in AA,
/// V = 0.21 + 5 log10(rΔ) + 0.026|α| + 4·10⁻⁹ α⁴. Sun: V = −26.74 + 5 log10 Δ.
/// `r`, `delta` in AU; `alpha` in deg; `u_obs` the unit vector body→observer
/// (ICRF); `year` for Neptune's secular brightening.
pub fn magnitude(body: Body, r: f64, delta: f64, alpha: f64, u_obs: Vec3, year: f64) -> f64 {
    use Body::*;
    let m5 = 5.0 * (r * delta).log10();
    let a = alpha;
    match body {
        Sun => -26.74 + 5.0 * delta.log10(),
        Moon => 0.21 + m5 + 0.026 * a.abs() + 4e-9 * a.powi(4),
        Mercury => {
            m5 - 0.613 + 6.3280e-02 * a - 1.6336e-03 * a.powi(2) + 3.3644e-05 * a.powi(3)
                - 3.4265e-07 * a.powi(4)
                + 1.6893e-09 * a.powi(5)
                - 3.0334e-12 * a.powi(6)
        }
        Venus => {
            if a <= 163.7 {
                m5 - 4.384 - 1.044e-03 * a + 3.687e-04 * a.powi(2) - 2.814e-06 * a.powi(3)
                    + 8.938e-09 * a.powi(4)
            } else {
                m5 + 236.05828 - 2.81914 * a + 8.39034e-03 * a.powi(2)
            }
        }
        Mars => {
            if a <= 50.0 {
                m5 - 1.601 + 2.267e-02 * a - 1.302e-04 * a.powi(2)
            } else {
                m5 - 0.367 - 0.02573 * a + 3.445e-04 * a.powi(2)
            }
        }
        Jupiter => {
            if a <= 12.0 {
                m5 - 9.395 - 3.7e-04 * a + 6.16e-04 * a.powi(2)
            } else {
                let x = a / 180.0;
                m5 - 9.428
                    - 2.5
                        * (1.0 - 1.507 * x - 0.363 * x.powi(2) - 0.062 * x.powi(3) + 2.809 * x.powi(4)
                            - 1.876 * x.powi(5))
                        .log10()
            }
        }
        Saturn => {
            let sb = dot(u_obs, pole(40.589, 83.537)).abs();
            let beta = sb.asin() / DEG;
            if a <= 6.5 && beta <= 27.0 {
                m5 - 8.914 - 1.825 * sb + 0.026 * a - 0.378 * sb * (-2.25 * a).exp()
            } else if a <= 6.0 {
                m5 - 8.95 - 3.7e-04 * a + 6.16e-04 * a.powi(2)
            } else {
                m5 - 8.94 + 2.446e-04 * a + 2.672e-04 * a.powi(2) - 1.505e-06 * a.powi(3)
                    + 4.767e-09 * a.powi(4)
            }
        }
        Uranus => {
            let lat = dot(u_obs, pole(257.311, -15.175)).clamp(-1.0, 1.0).asin().abs() / DEG;
            m5 - 7.110 - 8.4e-04 * lat + 6.587e-03 * a + 1.045e-04 * a.powi(2)
        }
        Neptune => {
            let h = if year > 2000.0 {
                -7.00
            } else if year < 1980.0 {
                -6.89
            } else {
                -6.89 - 0.0054 * (year - 1980.0)
            };
            m5 + h + 7.944e-03 * a + 9.617e-05 * a.powi(2)
        }
        Pluto => m5 - 1.01,
        _ => f64::NAN,
    }
}

fn wrap_pi(a: f64) -> f64 {
    let r = norm_rad(a + PI) - PI;
    if r <= -PI {
        r + TAU
    } else {
        r
    }
}

// ---------------------------------------------------------------------------------------------
// Engine API
// ---------------------------------------------------------------------------------------------

impl Engine {
    fn obs_opts(origin: Origin) -> CalcOptions {
        CalcOptions { origin, equatorial: true, speed: false, ..Default::default() }
    }

    /// Topocentric (or "Hindu" geocentric, latitude-free) apparent RA/Dec of
    /// date (rad), and semidiameter (rad), of a body at a UT instant.
    pub fn rise_set_place(
        &self,
        body: Body,
        jd_ut: f64,
        obs: &Observer,
        o: &RiseSetOptions,
    ) -> Result<(f64, f64, f64)> {
        self.rise_set_place_epoch(body, &epoch_ut(jd_ut, self.frame_model(), &self.delta_t_model()), obs, o)
    }

    /// `rise_set_place` at a prepared epoch.
    pub fn rise_set_place_epoch(
        &self,
        body: Body,
        ep: &Epoch,
        obs: &Observer,
        o: &RiseSetOptions,
    ) -> Result<(f64, f64, f64)> {
        let naif = body
            .naif()
            .ok_or(Error::InvalidArgument(format!("rise/set needs a physical body, got {body:?}")))?;
        let pack = self.pack();
        let corr = Corrections::APPARENT;
        let (ra, dec, dist_km) = if o.geocentric_no_lat {
            let os = ObsState::new(pack, ep, &Origin::Geocentric)?;
            let frame = Frame::EclipticOfDate { true_equinox: true };
            let ((l, _, r), ..) = ephem::place_obs(pack, ep, &os, naif, corr, frame, LightTime::Iterate)?;
            let eq = cotrans([l / DEG, 0.0, 1.0], -ep.frames.eps_true / DEG);
            (eq[0] * DEG, eq[1] * DEG, r)
        } else {
            let os = ObsState::new(pack, ep, &Origin::Topocentric(*obs))?;
            let frame = Frame::EquatorialOfDate { true_equinox: true };
            let ((ra, dec, r), ..) = ephem::place_obs(pack, ep, &os, naif, corr, frame, LightTime::Iterate)?;
            (ra, dec, r)
        };
        let dist_km = match (o.fixed_disc_size, body) {
            (true, Body::Sun) => AU_KM,
            (true, Body::Moon) => MOON_MEAN_DIST_KM,
            _ => dist_km,
        };
        let sd = (body_radius_km(body) / dist_km).min(1.0).asin();
        Ok((ra, dec, sd))
    }

    /// Local apparent sidereal time (rad) at a UT instant and east longitude.
    pub fn last(&self, jd_ut: f64, lon_deg: f64) -> f64 {
        norm_rad(self.prepare_ut(jd_ut, false).epoch.gast() + lon_deg * DEG)
    }

    /// Next rise / set / meridian transit of `body` after `jd_ut` (UT).
    pub fn rise_set(
        &self,
        body: Body,
        jd_ut: f64,
        kind: RiseSetKind,
        obs: &Observer,
        o: &RiseSetOptions,
    ) -> Result<RiseSet> {
        // Meridian transits always use the topocentric apparent place: the
        // disc / refraction / "Hindu" options do not affect them (measured).
        let plain = RiseSetOptions::default();
        let po =
            if matches!(kind, RiseSetKind::UpperTransit | RiseSetKind::LowerTransit) { &plain } else { o };
        let place = |ep: &Epoch| self.rise_set_place_epoch(body, ep, obs, po);
        self.rise_set_target(&place, jd_ut, kind, obs, o)
    }

    /// Altitude (rad) of the disc point chosen in `o` above the rise/set horizon
    /// of `o` at a UT instant: positive while the body counts as risen in the
    /// sense of `rise_set` (same place, refraction and horizon rules).
    pub fn rise_set_altitude(
        &self,
        body: Body,
        jd_ut: f64,
        obs: &Observer,
        o: &RiseSetOptions,
    ) -> Result<f64> {
        let lat = obs.lat_deg * DEG;
        let h0 = match o.refraction {
            None => o.horizon_deg,
            Some(atm) if o.horizon_deg < 0.0 => {
                o.horizon_deg.min(-horizon_refraction(o.refraction_model, 0.0, &atm, obs.height_m))
            }
            Some(atm) => {
                o.horizon_deg - horizon_refraction(o.refraction_model, o.horizon_deg, &atm, obs.height_m)
            }
        } * DEG;
        let limb = match o.disc {
            DiscPoint::UpperLimb => 1.0,
            DiscPoint::Center => 0.0,
            DiscPoint::LowerLimb => -1.0,
        };
        let ep = epoch_ut(jd_ut, self.frame_model(), &self.delta_t_model());
        let (ra, dec, sd) = self.rise_set_place_epoch(body, &ep, obs, o)?;
        let (alt, _) = horizontal(ra, dec, norm_rad(ep.gast() + obs.lon_deg * DEG), lat);
        Ok(alt + limb * sd - h0)
    }

    /// `rise_set` for any target given by its apparent place (hook for fixed stars).
    pub fn rise_set_target(
        &self,
        place: &TargetFn<'_>,
        jd_ut: f64,
        kind: RiseSetKind,
        obs: &Observer,
        o: &RiseSetOptions,
    ) -> Result<RiseSet> {
        let ei = EpochInterp::new_with(jd_ut, self.frame_model(), self.delta_t_model());
        match kind {
            RiseSetKind::UpperTransit => {
                return Ok(RiseSet::Event(self.next_transit(&ei, place, jd_ut, obs, true)?))
            }
            RiseSetKind::LowerTransit => {
                return Ok(RiseSet::Event(self.next_transit(&ei, place, jd_ut, obs, false)?))
            }
            _ => {}
        }
        let rising = kind == RiseSetKind::Rise;
        let lat = obs.lat_deg * DEG;
        let h0 = match o.refraction {
            None => o.horizon_deg,
            // A horizon below the apparent horizon: the classic API keeps the
            // true altitude of the refracted zero-altitude horizon unless the
            // given horizon is lower still (measured at sea level).
            Some(atm) if o.horizon_deg < 0.0 => {
                o.horizon_deg.min(-horizon_refraction(o.refraction_model, 0.0, &atm, obs.height_m))
            }
            Some(atm) => {
                o.horizon_deg - horizon_refraction(o.refraction_model, o.horizon_deg, &atm, obs.height_m)
            }
        } * DEG;
        let limb = match o.disc {
            DiscPoint::UpperLimb => 1.0,
            DiscPoint::Center => 0.0,
            DiscPoint::LowerLimb => -1.0,
        };
        // f(t) = true altitude of the disc point − h0.
        let f = |t: f64| -> Result<f64> {
            let ep = ei.at(t);
            let (ra, dec, sd) = place(&ep)?;
            let (alt, _) = horizontal(ra, dec, norm_rad(ep.gast() + obs.lon_deg * DEG), lat);
            Ok(alt + limb * sd - h0)
        };
        let t1 = jd_ut + o.window_days;
        // Altitude extrema ≈ meridian transits: bracket between consecutive ones.
        let mut knots = vec![jd_ut];
        let mut tr = self.transits(&ei, place, jd_ut, t1, obs, true)?;
        tr.extend(self.transits(&ei, place, jd_ut, t1, obs, false)?);
        tr.sort_by(f64::total_cmp);
        knots.extend(tr.into_iter().filter(|&t| t > jd_ut && t < t1));
        knots.push(t1);
        let mut fa = f(knots[0])?;
        let (mut any_above, mut any_below) = (fa >= 0.0, fa < 0.0);
        for w in knots.windows(2) {
            let fb = f(w[1])?;
            any_above |= fb >= 0.0;
            any_below |= fb < 0.0;
            let hit = if rising { fa < 0.0 && fb >= 0.0 } else { fa >= 0.0 && fb < 0.0 };
            if hit {
                let t = brent(&f, w[0], w[1], fa, fb)?;
                return Ok(RiseSet::Event(t));
            }
            fa = fb;
        }
        // No event: classify by the altitude at the knots.
        Ok(if any_above && !any_below {
            RiseSet::AlwaysAbove
        } else if any_below && !any_above {
            RiseSet::AlwaysBelow
        } else if rising {
            // e.g. it sets inside the window and does not rise again before its end.
            RiseSet::AlwaysBelow
        } else {
            RiseSet::AlwaysAbove
        })
    }

    /// Upper (`upper = true`) or lower meridian transits in (t0, t1].
    fn transits(
        &self,
        ei: &EpochInterp,
        place: &TargetFn<'_>,
        t0: f64,
        t1: f64,
        obs: &Observer,
        upper: bool,
    ) -> Result<Vec<f64>> {
        // Consecutive transits are ≥ 0.997 d apart (sidereal day; slower for
        // bodies moving east), so stop as soon as the next cannot fall in (t0, t1].
        const MIN_PERIOD: f64 = 0.99;
        let mut out = Vec::new();
        let mut t = t0;
        loop {
            let tr = self.next_transit(ei, place, t, obs, upper)?;
            if tr > t1 {
                return Ok(out);
            }
            out.push(tr);
            if tr + MIN_PERIOD > t1 {
                return Ok(out);
            }
            t = tr + 1e-4;
        }
    }

    /// First transit strictly after `t0`: secant iteration on the wrapped
    /// hour angle H − (0 | π), which increases monotonically at ≈ 2π/day.
    fn next_transit(
        &self,
        ei: &EpochInterp,
        place: &TargetFn<'_>,
        t0: f64,
        obs: &Observer,
        upper: bool,
    ) -> Result<f64> {
        let off = if upper { 0.0 } else { PI };
        let g = |t: f64| -> Result<f64> {
            let ep = ei.at(t);
            let (ra, _, _) = place(&ep)?;
            Ok(wrap_pi(ep.gast() + obs.lon_deg * DEG - ra - off))
        };
        // Sidereal rate; the target's own RA motion (Moon ≤ ~16°/day) slows it by < 5 %.
        let mut rate = TAU * 1.002_737_909_35;
        let g0 = g(t0)?;
        let need = if g0 < 0.0 { -g0 } else { TAU - g0 };
        let mut t = t0 + need / (rate * 0.98);
        let mut gt = g(t)?;
        for _ in 0..30 {
            let dt = -gt / rate;
            let tn = t + dt;
            if dt.abs() < 1e-9 {
                return Ok(tn);
            }
            let gn = g(tn)?;
            if (gn - gt).abs() > 1e-15 && (gn - gt).abs() < 1.0 {
                rate = (gn - gt) / (tn - t);
            }
            t = tn;
            gt = gn;
        }
        Ok(t)
    }

    /// Horizontal coordinates of a direction of date, the classic `azalt`:
    /// returns (azimuth from south westwards, true altitude, apparent altitude), deg.
    /// The apparent altitude uses `refrac` (pressure 0 → standard at height).
    pub fn azalt(
        &self,
        jd_ut: f64,
        input: AzaltInput,
        obs: &Observer,
        atm: &Atmosphere,
        coord: [f64; 3],
    ) -> [f64; 3] {
        let ep = self.prepare_ut(jd_ut, false).epoch;
        let (ra, dec) = match input {
            AzaltInput::Equatorial => (coord[0], coord[1]),
            AzaltInput::Ecliptic => {
                let e = cotrans(coord, -ep.frames.eps_true / DEG);
                (e[0], e[1])
            }
        };
        let last = norm_rad(ep.gast() + obs.lon_deg * DEG);
        let (alt, az) = horizontal(ra * DEG, dec * DEG, last, obs.lat_deg * DEG);
        let alt = alt / DEG;
        let app = refrac(alt, atm.pressure_at(obs.height_m), atm.temperature_c, RefracDir::TrueToApparent);
        [az / DEG, alt, app]
    }

    /// Planetary phenomena of `body` at `jd_ut` (UT) for the observer and
    /// corrections in `opt` (ecliptic/equatorial/sidereal settings are irrelevant).
    pub fn pheno_ut(&self, body: Body, jd_ut: f64, opt: &CalcOptions) -> Result<Pheno> {
        let ep = self.prepare_ut(jd_ut, false).epoch;
        self.pheno_epoch(body, &ep, opt)
    }

    pub fn pheno_tt(&self, body: Body, jd_tt: f64, opt: &CalcOptions) -> Result<Pheno> {
        let ep = self.prepare_tt(jd_tt, false).epoch;
        self.pheno_epoch(body, &ep, opt)
    }

    fn pheno_epoch(&self, body: Body, ep: &ephem::Epoch, opt: &CalcOptions) -> Result<Pheno> {
        let naif = body
            .naif()
            .ok_or(Error::InvalidArgument(format!("pheno needs a physical body, got {body:?}")))?;
        let pack = self.pack();
        let origin = match opt.origin {
            o @ (Origin::Geocentric | Origin::Topocentric(_)) => o,
            _ => Origin::Geocentric,
        };
        let obs = ObsState::new(pack, ep, &origin)?;
        let corr = opt.corrections;
        let (u_body, tau, _) = ephem::apparent(pack, ep, &obs, naif, corr, LightTime::Iterate)?;
        let delta_km = norm(u_body);
        let diameter = 2.0 * (body_radius_km(body) / delta_km).min(1.0).asin() / DEG;
        if body == Body::Sun {
            return Ok(Pheno {
                diameter,
                magnitude: magnitude(body, 1.0, delta_km / AU_KM, 0.0, [0.0; 3], 0.0),
                ..Default::default()
            });
        }
        let (u_sun, ..) = ephem::apparent(pack, ep, &obs, ephem::SUN, corr, LightTime::Iterate)?;
        let elongation = (dot(u_body, u_sun) / (norm(u_body) * norm(u_sun))).clamp(-1.0, 1.0).acos() / DEG;
        // Distances for the magnitude and the planets' phase angle: geometric,
        // body at the emission time t − τ, Sun at t − τ, observer at t.
        let t_emit = ep.jd_tdb - if corr.light_time { tau } else { 0.0 };
        let p = pack.position(naif, ephem::SSB, t_emit)?;
        let s = pack.position(ephem::SUN, ephem::SSB, t_emit)?;
        let to_obs = sub(obs.pos, p);
        let to_sun = sub(s, p);
        let alpha = if body == Body::Moon {
            // Moon: from the apparent geocentric and apparent heliocentric
            // vectors, the classic definition (matches it to < 1″; differs
            // from the geometric one by up to ~40″, i.e. ~1e-4 in the phase).
            let helio = ObsState::new(pack, ep, &Origin::Heliocentric)?;
            let (h_vec, ..) = ephem::apparent(pack, ep, &helio, naif, corr, LightTime::Iterate)?;
            (dot(u_body, h_vec) / (norm(u_body) * norm(h_vec))).clamp(-1.0, 1.0).acos() / DEG
        } else {
            (dot(to_obs, to_sun) / (norm(to_obs) * norm(to_sun))).clamp(-1.0, 1.0).acos() / DEG
        };
        let r_au = norm(to_sun) / AU_KM;
        let d_au = norm(to_obs) / AU_KM;
        let year = crate::time::decimal_year(ep.jd_tt);
        let mag = magnitude(body, r_au, d_au, alpha, scale(to_obs, 1.0 / norm(to_obs)), year);
        let hpar = if body == Body::Moon { (EARTH_EQ_RADIUS_KM / delta_km).asin() / DEG } else { 0.0 };
        Ok(Pheno {
            phase_angle: alpha,
            phase: (1.0 + (alpha * DEG).cos()) / 2.0,
            elongation,
            diameter,
            magnitude: mag,
            horizontal_parallax: hpar,
        })
    }

    /// Next time (UT) after `jd_ut` when `body`'s longitude (per `opt`:
    /// tropical/sidereal, apparent/true) equals `lon_deg`. Sun and Moon only
    /// (monotonic longitude), the classic `solcross_ut` / `mooncross_ut`.
    pub fn longitude_cross_ut(&self, body: Body, lon_deg: f64, jd_ut: f64, opt: &CalcOptions) -> Result<f64> {
        let mean_speed = match body {
            Body::Sun => 0.985_647_36,
            Body::Moon => 13.176_358,
            _ => return Err(Error::InvalidArgument("longitude crossing: Sun or Moon only".into())),
        };
        let opt = CalcOptions { equatorial: false, speed: true, xyz: false, radians: false, ..*opt };
        let pos = self.calc_ut(body, jd_ut, &opt)?;
        let mut t = jd_ut + norm_deg(lon_deg - pos.lon) / mean_speed;
        for _ in 0..50 {
            let p = self.calc_ut(body, t, &opt)?;
            let dt = diff_deg(lon_deg, p.lon) / p.lon_speed;
            t += dt;
            if dt.abs() < 1e-9 {
                break;
            }
        }
        if t <= jd_ut {
            // Started exactly on (or a hair past) the target: take the next one.
            return self.longitude_cross_ut(body, lon_deg, jd_ut + 1e-6, &opt).map(|x| x.max(t));
        }
        Ok(t)
    }

    pub fn solcross_ut(&self, lon_deg: f64, jd_ut: f64, opt: &CalcOptions) -> Result<f64> {
        self.longitude_cross_ut(Body::Sun, lon_deg, jd_ut, opt)
    }

    pub fn mooncross_ut(&self, lon_deg: f64, jd_ut: f64, opt: &CalcOptions) -> Result<f64> {
        self.longitude_cross_ut(Body::Moon, lon_deg, jd_ut, opt)
    }

    /// Equation of time (days) at `jd_ut`: local apparent time − local mean
    /// time, E = GAST − α☉ − (UT − 12h), with α☉ the apparent geocentric right
    /// ascension of date (true equinox) — the definition the classic
    /// `time_equ` follows (measured to < 1 µs of this formula).
    pub fn time_equ(&self, jd_ut: f64) -> Result<f64> {
        let p = self.prepare_ut(jd_ut, false);
        let sun = self.calc_prepared(Body::Sun, &p, &Self::obs_opts(Origin::Geocentric))?;
        let ut = (jd_ut - 0.5).rem_euclid(1.0) * TAU;
        Ok(wrap_pi(p.epoch.gast() - sun.lon * DEG + PI - ut) / TAU)
    }

    /// Local mean time → local apparent (true solar) time, JD in both.
    pub fn lmt_to_lat(&self, jd_lmt: f64, lon_deg: f64) -> Result<f64> {
        Ok(jd_lmt + self.time_equ(jd_lmt - lon_deg / 360.0)?)
    }

    /// Local apparent time → local mean time (inverse of `lmt_to_lat`).
    pub fn lat_to_lmt(&self, jd_lat: f64, lon_deg: f64) -> Result<f64> {
        let mut lmt = jd_lat;
        for _ in 0..4 {
            lmt = jd_lat - self.time_equ(lmt - lon_deg / 360.0)?;
        }
        Ok(lmt)
    }

    /// True (apparent) solar time at east longitude `lon_deg` for a UT instant,
    /// as a JD-like value: UT + lon/360 + E(UT). The fractional day − 0.5 is
    /// the sundial time (what the Chinese true-solar-time correction needs).
    pub fn true_solar_time(&self, jd_ut: f64, lon_deg: f64) -> Result<f64> {
        Ok(jd_ut + lon_deg / 360.0 + self.time_equ(jd_ut)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cotrans_roundtrip() {
        let e = 23.4392911;
        let eq = cotrans([123.4, 5.6, 1.0], -e);
        let ec = cotrans(eq, e);
        assert!((ec[0] - 123.4).abs() < 1e-12 && (ec[1] - 5.6).abs() < 1e-12);
        // Ecliptic longitude 90° → RA 90°, Dec = ε.
        let q = cotrans([90.0, 0.0, 1.0], -e);
        assert!((q[0] - 90.0).abs() < 1e-12 && (q[1] - e).abs() < 1e-12);
    }

    #[test]
    fn cotrans_sp_matches_finite_difference() {
        let c = [200.0, -3.0, 1.0, 13.0, 0.5, 0.0];
        let e = -23.44;
        let h = 1e-5;
        let a = cotrans([c[0] - c[3] * h, c[1] - c[4] * h, 1.0], e);
        let b = cotrans([c[0] + c[3] * h, c[1] + c[4] * h, 1.0], e);
        let r = cotrans_sp(c, e);
        assert!((diff_deg(b[0], a[0]) / (2.0 * h) - r[3]).abs() < 1e-6);
        assert!(((b[1] - a[1]) / (2.0 * h) - r[4]).abs() < 1e-6);
    }

    #[test]
    fn epoch_interp_matches_exact() {
        let m = crate::context::API_FRAME_MODEL;
        for t0 in [1721424.3, 2451545.2, 2460400.7, 2597000.9] {
            let ei = EpochInterp::new(t0, m);
            for k in 0..=30 {
                let t = t0 - 0.15 + k as f64 * 0.05;
                let (a, b) = (ei.at(t), Epoch::from_ut_with(t, m));
                let mut d: f64 = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        d = d.max((a.frames.bpn[i][j] - b.frames.bpn[i][j]).abs());
                        d = d.max((a.frames.ecl_true[i][j] - b.frames.ecl_true[i][j]).abs());
                    }
                }
                assert!(d < 1e-3 * ARCSEC, "matrix {} uas at {t}", d / ARCSEC * 1e6); // < 1 mas (seen: 0.25 mas)
                assert!(
                    (a.gast() - b.gast()).abs() < 1e-3 * ARCSEC,
                    "gast {} uas",
                    (a.gast() - b.gast()).abs() / ARCSEC * 1e6
                );
                assert!(((a.jd_tt - b.jd_tt) * 86400.0).abs() < 1e-4, "tt"); // JD f64 resolution ~40 µs
            }
        }
    }

    #[test]
    fn refraction_formulas() {
        // Meeus example 16.a-ish magnitudes: ~34.5' at the horizon, ~1' at 45°.
        assert!((refraction_bennett(0.0, 1010.0, 10.0) * 60.0 - 34.48).abs() < 0.05);
        assert!((refraction_saemundsson(45.0, 1010.0, 10.0) * 60.0 - 1.0).abs() < 0.05);
        // Sæmundsson inverted at the horizon: h + R_S(h) = 0.
        let atm = Atmosphere { pressure_hpa: 1013.25, temperature_c: 0.0 };
        let r = horizon_refraction(RefractionModel::Saemundsson, 0.0, &atm, 0.0);
        assert!((-r + refraction_saemundsson(-r, 1013.25, 0.0)).abs() < 1e-10);
        assert!((r * 60.0 - 36.10).abs() < 0.01, "{}", r * 60.0);
        // Sinclair at the horizon: 34.46′·(933.25/930)/(1 + 8e-5·73.46·(−10)).
        let r = horizon_refraction(RefractionModel::Sinclair1982, 0.0, &atm, 0.0) * 60.0;
        assert!((r - 34.46 * 933.25 / 930.0 / (1.0 - 8e-4 * 73.46)).abs() < 1e-9);
        assert!((sinclair_switch_altitude() - 17.904).abs() < 0.001);
        // Continuity at the switch and true↔apparent round trip.
        let hx = sinclair_switch_altitude();
        assert!(
            (refraction_sinclair(hx - 1e-9, 1013.25, 0.0) - refraction_sinclair(hx + 1e-9, 1013.25, 0.0))
                .abs()
                < 1e-9
        );
        for h in [0.0, 0.3, 2.0, 10.0, 17.9, 18.0, 45.0, 89.0] {
            let t = refrac(h, 1013.25, 15.0, RefracDir::ApparentToTrue);
            assert!((refrac(t, 1013.25, 15.0, RefracDir::TrueToApparent) - h).abs() < 1e-10, "{h}");
        }
        assert!((standard_pressure_hpa(1000.0) - 898.707).abs() < 0.01);
    }
}
