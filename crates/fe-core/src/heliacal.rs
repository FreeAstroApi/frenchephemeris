// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Heliacal events and visual limiting magnitude: the classic `heliacal_ut`,
//! `heliacal_pheno_ut` and `vis_limit_mag`.
//!
//! Visibility model: B. E. Schaefer, "To the visual limits", Sky & Telescope
//! 95(5), 57 (May 1998) — the published BASIC program VISLIMIT (extinction
//! from Rayleigh, aerosol, ozone and water-vapour terms; night-sky, moonlight,
//! twilight and daylight sky brightness; Blackwell–Hecht threshold). Background
//! in Schaefer, "Astronomy and the limits of vision", Vistas Astron. 36, 311
//! (1993) and "Telescopic limiting magnitudes", PASP 102, 212 (1990) (optical
//! factors). Only the V band enters the result (`B(V)` and `DM(V)`), so only it
//! is evaluated.
//!
//! Heliacal search ("VR" method of Reijs's documentation of Schaefer's model):
//! on each local morning (rising events) or evening (setting events), the
//! object is visible at an instant when its apparent magnitude is brighter
//! than the limiting magnitude at its (topocentric, unrefracted) altitude with
//! the Sun below the horizon. The event day is the first visible morning /
//! evening after a run of invisible ones (morning first, evening first) or the
//! last visible one before invisibility (evening last, morning last); the
//! result is (first visible instant, best instant, last visible instant) of
//! that day.
//!
//! Positions: our apparent topocentric places (JPL), magnitudes from `pheno`
//! (Mallama & Hilton 2018) and the star catalogue (SIMBAD V), refraction for
//! the apparent altitude from Sinclair (1982) / Explanatory Supplement (1992).

use crate::context::{Body, CalcOptions, Engine};
use crate::ephem::{Observer, Origin};
use crate::events::{
    self, Atmosphere, AzaltInput, DiscPoint, RefracDir, RiseSet, RiseSetKind, RiseSetOptions,
};
use crate::math::*;
use crate::stars::{self, Star};
use crate::{Error, Result};

/// Event types (classic numbering).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HelEvent {
    /// Heliacal rising: first morning visibility.
    MorningFirst = 1,
    /// Heliacal setting: last evening visibility.
    EveningLast = 2,
    /// Evening first (Mercury, Venus, Moon).
    EveningFirst = 3,
    /// Morning last (Mercury, Venus, Moon).
    MorningLast = 4,
}

impl HelEvent {
    pub fn from_id(id: i32) -> Result<HelEvent> {
        Ok(match id {
            1 => HelEvent::MorningFirst,
            2 => HelEvent::EveningLast,
            3 => HelEvent::EveningFirst,
            4 => HelEvent::MorningLast,
            5 | 6 => {
                return Err(Error::InvalidArgument(format!(
                    "acronychal rising/setting (event type {id}) is not provided"
                )))
            }
            _ => return Err(Error::InvalidArgument(format!("unknown heliacal event type {id}"))),
        })
    }
    fn morning(self) -> bool {
        matches!(self, HelEvent::MorningFirst | HelEvent::MorningLast)
    }
    fn first(self) -> bool {
        matches!(self, HelEvent::MorningFirst | HelEvent::EveningFirst)
    }
}

/// The object: a solar-system body or a catalogue star.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HelTarget {
    Body(Body),
    Star(Star),
}

impl HelTarget {
    /// Planet name (`"Venus"`, case-insensitive) or any star name/designation.
    pub fn parse(name: &str) -> Result<HelTarget> {
        let n = name.trim().to_ascii_lowercase();
        let body = match n.as_str() {
            "sun" => Some(Body::Sun),
            "moon" => Some(Body::Moon),
            "mercury" => Some(Body::Mercury),
            "venus" => Some(Body::Venus),
            "mars" => Some(Body::Mars),
            "jupiter" => Some(Body::Jupiter),
            "saturn" => Some(Body::Saturn),
            "uranus" => Some(Body::Uranus),
            "neptune" => Some(Body::Neptune),
            "pluto" => Some(Body::Pluto),
            _ => None,
        };
        match body {
            Some(b) => Ok(HelTarget::Body(b)),
            None => Ok(HelTarget::Star(stars::lookup(name)?)),
        }
    }

    fn label(&self) -> String {
        match self {
            HelTarget::Body(b) => format!("{b:?}"),
            HelTarget::Star(s) => s.full_name(),
        }
    }

    /// Which event types exist for this object (classic rules).
    pub fn check_event(&self, ev: HelEvent) -> Result<()> {
        let ok = match self {
            HelTarget::Body(Body::Sun) => false,
            HelTarget::Body(Body::Moon) => matches!(ev, HelEvent::EveningFirst | HelEvent::MorningLast),
            HelTarget::Body(Body::Mercury | Body::Venus) => true,
            _ => matches!(ev, HelEvent::MorningFirst | HelEvent::EveningLast),
        };
        if ok {
            Ok(())
        } else {
            Err(Error::InvalidArgument(format!(
                "event type {} does not exist for {}",
                ev as i32,
                self.label()
            )))
        }
    }
}

/// Atmosphere (`datm`): pressure hPa, temperature °C, relative humidity %, and
/// `vr_or_k`: ≥ 1 meteorological range (km), in (0, 1) the total V extinction
/// coefficient (mag/airmass), 0 → computed from the other parameters.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HelAtmosphere {
    pub pressure_hpa: f64,
    pub temperature_c: f64,
    pub humidity_pct: f64,
    pub vr_or_k: f64,
}

impl HelAtmosphere {
    /// All-zero input → defaults 1013.25 hPa, 15 °C, 40 %, adjusted to the
    /// observer height by the standard atmosphere (pressure, 6.5 K/km lapse).
    pub fn resolved(&self, height_m: f64) -> HelAtmosphere {
        let mut a = *self;
        if a.pressure_hpa == 0.0 && a.temperature_c == 0.0 && a.humidity_pct == 0.0 {
            a.pressure_hpa = events::standard_pressure_hpa(height_m);
            a.temperature_c = 15.0 - 0.0065 * height_m;
            a.humidity_pct = 40.0;
        }
        if a.pressure_hpa == 0.0 {
            a.pressure_hpa = events::standard_pressure_hpa(height_m);
        }
        a.humidity_pct = a.humidity_pct.clamp(1.0, 99.9);
        a
    }
}

/// Observer (`dobs`): age (years), Snellen ratio, and for optical aids
/// (`optical = true`): binocular flag, magnification, aperture (mm), transmission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelObserver {
    pub age: f64,
    pub snellen: f64,
    pub binocular: bool,
    pub magnification: f64,
    pub aperture_mm: f64,
    pub transmission: f64,
}

impl Default for HelObserver {
    fn default() -> Self {
        HelObserver {
            age: 36.0,
            snellen: 1.0,
            binocular: true,
            magnification: 1.0,
            aperture_mm: 0.0,
            transmission: 0.0,
        }
    }
}

impl HelObserver {
    /// Zero age / Snellen → defaults 36 / 1 (classic semantics).
    pub fn resolved(&self) -> HelObserver {
        let mut o = *self;
        if o.age <= 0.0 {
            o.age = 36.0;
        }
        if o.snellen <= 0.0 {
            o.snellen = 1.0;
        }
        o
    }
}

/// Flags (classic `HELFLAG_*` subset that changes the physics).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct HelFlags {
    /// `HELFLAG_OPTICAL_PARAMS`: use the binocular/magnification/aperture/transmission inputs.
    pub optical: bool,
    /// `HELFLAG_VISLIM_DARK`: behave as if the Sun (and, as observed, the Moon)
    /// were at the nadir.
    pub dark: bool,
    /// `HELFLAG_VISLIM_NOMOON`: behave as if the Moon were at the nadir.
    pub no_moon: bool,
    /// `HELFLAG_VISLIM_PHOTOPIC`: force photopic vision.
    pub photopic: bool,
}

/// Vision regime of the limiting magnitude.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Vision {
    Photopic = 0,
    Scotopic = 1,
}

// ---------------------------------------------------------------------------------------------
// Schaefer (1998) VISLIMIT, V band
// ---------------------------------------------------------------------------------------------

/// Inputs of the visibility model (degrees; zenith distances > 90 mean "below").
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkyGeometry {
    /// Zenith distance of the object.
    pub z_obj: f64,
    pub z_sun: f64,
    pub z_moon: f64,
    /// Angular distances object–Sun and object–Moon.
    pub sep_sun: f64,
    pub sep_moon: f64,
    /// Moon phase angle (0 = full).
    pub moon_phase_angle: f64,
    pub lat_deg: f64,
    pub height_m: f64,
    /// Month 1–12 (seasonal aerosol/ozone terms) and decimal year (solar cycle).
    pub month: f64,
    pub year: f64,
}

/// Result of the visibility model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VisModel {
    /// Limiting V magnitude for a point source at the object's position.
    pub limiting_mag: f64,
    /// Sky brightness, nanolamberts.
    pub sky_nl: f64,
    /// V extinction coefficient, mag/airmass.
    pub k: f64,
    /// Extinction along the line of sight, mag.
    pub dm: f64,
    pub vision: Vision,
}

const RD: f64 = PI / 180.0;
// V-band constants of VISLIMIT (Schaefer 1998).
const WA_V: f64 = 0.55;
const MO_V: f64 = -11.05;
const OZ_V: f64 = 0.031;
const WT_V: f64 = 0.031;
const BO_V: f64 = 1.0e-13;
const CM_V: f64 = 0.00;
const MS_V: f64 = -26.74;

/// V extinction coefficient k and the line-of-sight extinction ΔM (mag) at
/// zenith distance `z` (deg), Schaefer (1998) extinction subroutine; `k_total`
/// or a meteorological range from `atm.vr_or_k` replace the computed values.
pub fn extinction_v(z: f64, atm: &HelAtmosphere, lat_deg: f64, height_m: f64, month: f64) -> (f64, f64) {
    let lt = lat_deg * RD;
    let ra = (month - 3.0) * 30.0 * RD;
    let sl = if lat_deg >= 0.0 { 1.0 } else { -1.0 };
    let zz = z.min(90.0) * RD;
    let c = zz.cos();
    let xg = 1.0 / (c + 0.0286 * (-10.5 * c).exp());
    let xa = 1.0 / (c + 0.0123 * (-24.5 * c).exp());
    let xo = 1.0 / (1.0 - (zz.sin() / (1.0 + 20.0 / 6378.0)).powi(2)).sqrt();
    let rh = atm.humidity_pct;
    let kr = 0.1066 * (-height_m / 8200.0).exp() * (WA_V / 0.55).powf(-4.0);
    let mut ka = 0.1 * (WA_V / 0.55).powf(-1.3) * (-height_m / 1500.0).exp();
    ka *= (1.0 - 0.32 / (rh / 100.0).ln()).powf(1.33) * (1.0 + sl * ra.sin());
    let ko = OZ_V * (3.0 + 0.4 * (lt * ra.cos() - (3.0 * lt).cos())) / 3.0;
    let kw = WT_V * 0.94 * (rh / 100.0) * (atm.temperature_c / 15.0).exp() * (-height_m / 8200.0).exp();
    let v = atm.vr_or_k;
    if v > 0.0 && v < 1.0 {
        // Total coefficient given: single (gas) airmass.
        return (v, v * xg);
    }
    if v >= 1.0 {
        // Meteorological range (Koschmieder: β = 3.912/VR per km) as the aerosol
        // term, over Schaefer's 1.5 km aerosol scale height, in magnitudes.
        ka = 1.086 * 3.912 / v * 1.5 * (-height_m / 1500.0).exp();
    }
    (kr + ka + ko + kw, kr * xg + ka * xa + ko * xo + kw * xg)
}

fn airmass_sky(z_deg: f64) -> f64 {
    if z_deg > 90.0 {
        return 40.0;
    }
    let c = (z_deg * RD).cos();
    1.0 / (c + 0.025 * (-11.0 * c).exp())
}

/// Schaefer (1998) limiting magnitude for the given geometry, atmosphere and observer.
pub fn schaefer_limit(g: &SkyGeometry, atm: &HelAtmosphere, ob: &HelObserver, flags: &HelFlags) -> VisModel {
    let z = g.z_obj.min(90.0);
    let (k, dm) = extinction_v(z, atm, g.lat_deg, g.height_m, g.month);
    let zz = z * RD;
    let x = airmass_sky(z);
    let (z_sun, z_moon) = (
        if flags.dark { 180.0 } else { g.z_sun },
        if flags.no_moon || flags.dark { 180.0 } else { g.z_moon },
    );
    let xm = airmass_sky(z_moon);
    let xs = airmass_sky(z_sun);
    let fx = 1.0 - 10f64.powf(-0.4 * k * x);
    // Dark night sky (with the 11-year solar cycle).
    #[allow(clippy::approx_constant)] // Schaefer's published 6.283, kept as printed
    let mut bn = BO_V * (1.0 + 0.3 * (6.283 * (g.year - 1992.0) / 11.0).cos());
    bn *= 0.4 + 0.6 / (1.0 - 0.96 * zz.sin().powi(2)).sqrt();
    bn *= 10f64.powf(-0.4 * k * x);
    // Moonlight.
    let am = g.moon_phase_angle;
    let mm = -12.73 + 0.026 * am.abs() + 4e-9 * am.powi(4) + CM_V;
    let c3 = 10f64.powf(-0.4 * k * xm);
    let rm = g.sep_moon.max(0.1);
    let fm = 6.2e7 / (rm * rm)
        + 10f64.powf(6.15 - rm / 40.0)
        + 10f64.powf(5.36) * (1.06 + (rm * RD).cos().powi(2));
    let bm = 10f64.powf(-0.4 * (mm - MO_V + 43.27)) * fx * (fm * c3 + 440_000.0 * (1.0 - c3));
    // Twilight.
    let hs = 90.0 - z_sun;
    let rs = g.sep_sun.max(0.1);
    let bt = 10f64.powf(-0.4 * (MS_V - MO_V + 32.5 - hs - z / (360.0 * k))) * (100.0 / rs) * fx;
    // Daylight.
    let c4 = 10f64.powf(-0.4 * k * xs);
    let fs = 6.2e7 / (rs * rs)
        + 10f64.powf(6.15 - rs / 40.0)
        + 10f64.powf(5.36) * (1.06 + (rs * RD).cos().powi(2));
    let bd = 10f64.powf(-0.4 * (MS_V - MO_V + 43.27)) * fx * (fs * c4 + 440_000.0 * (1.0 - c4));
    let mut b = bn + if bd < bt { bd } else { bt };
    if z_moon < 90.0 {
        b += bm;
    }
    let mut bl = b * 1e12 / 1.11e-3; // picoergs → nanolamberts
                                     // Optical aid (Schaefer 1990): sky seen through the instrument, light gain.
    let mut gain_mag = 0.0;
    if flags.optical {
        let de = 7.0 * (-0.5 * (ob.age / 100.0).powi(2)).exp(); // pupil (mm), Schaefer 1993
        let t = if ob.transmission > 0.0 { ob.transmission } else { 1.0 };
        let mag = ob.magnification.max(1.0);
        if !ob.binocular {
            gain_mag -= 2.5 * 1.41f64.log10();
        }
        if ob.aperture_mm > 0.0 {
            let fa = (ob.aperture_mm / de).powi(2);
            gain_mag += 2.5 * (fa * t).log10();
            bl *= t * fa / (mag * mag);
        }
    }
    let scotopic = bl < 1500.0 && !flags.photopic;
    let (c1, c2) = if scotopic {
        (10f64.powf(-9.8), 10f64.powf(-1.9))
    } else {
        (10f64.powf(-8.350001), 10f64.powf(-5.9))
    };
    let th = c1 * (1.0 + (c2 * bl).sqrt()).powi(2);
    let mn = -16.57 - 2.5 * th.log10() - dm + 5.0 * ob.snellen.log10() + gain_mag;
    VisModel {
        limiting_mag: mn,
        sky_nl: bl,
        k,
        dm,
        vision: if scotopic { Vision::Scotopic } else { Vision::Photopic },
    }
}

// ---------------------------------------------------------------------------------------------
// Ephemeris side
// ---------------------------------------------------------------------------------------------

/// Everything the model needs at one instant (degrees; azimuth from south, westwards).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelSky {
    pub obj_alt: f64,
    pub obj_app_alt: f64,
    pub obj_geo_alt: f64,
    pub obj_az: f64,
    pub sun_alt: f64,
    pub sun_geo_alt: f64,
    pub sun_az: f64,
    pub moon_alt: f64,
    pub moon_az: f64,
    pub sep_sun: f64,
    pub sep_moon: f64,
    pub moon_phase_angle: f64,
    pub obj_mag: f64,
    /// Ecliptic longitude difference object − Sun (deg, signed).
    pub dlon: f64,
    /// Object horizontal parallax (deg).
    pub parallax: f64,
    pub geometry: SkyGeometry,
}

fn sep(ra1: f64, de1: f64, ra2: f64, de2: f64) -> f64 {
    let (a, b) = (polar(ra1, de1), polar(ra2, de2));
    dot(a, b).clamp(-1.0, 1.0).acos() / DEG
}

fn polar(ra: f64, de: f64) -> Vec3 {
    let (sr, cr) = (ra * DEG).sin_cos();
    let (sd, cd) = (de * DEG).sin_cos();
    [cd * cr, cd * sr, sd]
}

pub struct HelContext<'a> {
    pub engine: &'a Engine,
    pub target: HelTarget,
    pub obs: Observer,
    pub atm: HelAtmosphere,
    pub observer: HelObserver,
    pub flags: HelFlags,
}

impl HelContext<'_> {
    fn place(
        &self,
        target: HelTarget,
        jd_ut: f64,
        origin: Origin,
        ecliptic: bool,
    ) -> Result<(f64, f64, f64)> {
        let opt = CalcOptions { origin, equatorial: !ecliptic, speed: false, ..Default::default() };
        let pos = match target {
            HelTarget::Body(b) => self.engine.calc_ut(b, jd_ut, &opt)?,
            HelTarget::Star(s) => {
                let p = self.engine.prepare_ut(jd_ut, false);
                self.engine.fixstar_prepared(&s, &p, &opt)?
            }
        };
        Ok((pos.lon, pos.lat, pos.dist))
    }

    fn altaz(&self, jd_ut: f64, ra: f64, de: f64) -> (f64, f64) {
        let r = self.engine.azalt(
            jd_ut,
            AzaltInput::Equatorial,
            &self.obs,
            &Atmosphere::default(),
            [ra, de, 1.0],
        );
        (r[1], r[0])
    }

    fn magnitude(&self, jd_ut: f64) -> Result<f64> {
        Ok(match self.target {
            HelTarget::Star(s) => s.magnitude(),
            HelTarget::Body(b) => self.engine.pheno_ut(b, jd_ut, &CalcOptions::default())?.magnitude,
        })
    }

    /// Geometry and magnitudes at an instant; `mag` reuses a known magnitude.
    pub fn sky(&self, jd_ut: f64, mag: Option<f64>) -> Result<HelSky> {
        let topo = Origin::Topocentric(self.obs);
        let (ro, dobj, dist) = self.place(self.target, jd_ut, topo, false)?;
        let (rs, ds, _) = self.place(HelTarget::Body(Body::Sun), jd_ut, topo, false)?;
        let (rm, dm, _) = self.place(HelTarget::Body(Body::Moon), jd_ut, topo, false)?;
        let (rog, dog, _) = self.place(self.target, jd_ut, Origin::Geocentric, false)?;
        let (rsg, dsg, _) = self.place(HelTarget::Body(Body::Sun), jd_ut, Origin::Geocentric, false)?;
        let (lo, ..) = self.place(self.target, jd_ut, Origin::Geocentric, true)?;
        let (ls, ..) = self.place(HelTarget::Body(Body::Sun), jd_ut, Origin::Geocentric, true)?;
        let (obj_alt, obj_az) = self.altaz(jd_ut, ro, dobj);
        let (sun_alt, sun_az) = self.altaz(jd_ut, rs, ds);
        let (moon_alt, moon_az) = self.altaz(jd_ut, rm, dm);
        let (obj_geo_alt, _) = self.altaz(jd_ut, rog, dog);
        let (sun_geo_alt, _) = self.altaz(jd_ut, rsg, dsg);
        let moon_pa = self.engine.pheno_ut(Body::Moon, jd_ut, &CalcOptions::default())?.phase_angle;
        let obj_mag = match mag {
            Some(m) => m,
            None => self.magnitude(jd_ut)?,
        };
        let ymd = crate::time::revjul(jd_ut, crate::time::Calendar::Gregorian);
        let atm = &self.atm;
        let app = events::refrac(obj_alt, atm.pressure_hpa, atm.temperature_c, RefracDir::TrueToApparent);
        let parallax = match self.target {
            HelTarget::Body(_) if dist > 0.0 => (6378.137 / (dist * AU_KM)).min(1.0).asin() / DEG,
            _ => 0.0,
        };
        let (sep_sun, sep_moon) = (sep(ro, dobj, rs, ds), sep(ro, dobj, rm, dm));
        let geometry = SkyGeometry {
            z_obj: 90.0 - obj_alt,
            z_sun: 90.0 - sun_alt,
            z_moon: 90.0 - moon_alt,
            sep_sun,
            sep_moon,
            moon_phase_angle: moon_pa,
            lat_deg: self.obs.lat_deg,
            height_m: self.obs.height_m,
            month: ymd.1 as f64,
            year: crate::time::decimal_year(jd_ut),
        };
        Ok(HelSky {
            obj_alt,
            obj_app_alt: app,
            obj_geo_alt,
            obj_az,
            sun_alt,
            sun_geo_alt,
            sun_az,
            moon_alt,
            moon_az,
            sep_sun,
            sep_moon,
            moon_phase_angle: moon_pa,
            obj_mag,
            dlon: diff_deg(lo, ls),
            parallax,
            geometry,
        })
    }

    /// Limiting magnitude at an instant (None when the object is below the horizon).
    pub fn limit(&self, sky: &HelSky) -> Option<VisModel> {
        if sky.obj_alt <= 0.0 {
            return None;
        }
        let mut g = sky.geometry;
        if self.target == HelTarget::Body(Body::Moon) {
            g.z_moon = 180.0; // the Moon does not brighten the sky against itself
        }
        Some(schaefer_limit(&g, &self.atm, &self.observer, &self.flags))
    }

    /// Visibility margin (limiting − object magnitude); −∞ when not observable
    /// (object below the horizon or Sun above it). `fixed` = (object
    /// magnitude, Moon phase angle, month, decimal year), constant over a window.
    /// `margin` restricted to the object being east (`Some(true)`, rising side)
    /// or west (`Some(false)`, setting side) of the meridian: heliacal events
    /// are seen in the eastern morning sky or the western evening sky.
    fn margin_side(&self, jd_ut: f64, fixed: (f64, f64, f64, f64), east: Option<bool>) -> Result<f64> {
        let topo = Origin::Topocentric(self.obs);
        let (rs, ds, _) = self.place(HelTarget::Body(Body::Sun), jd_ut, topo, false)?;
        let (sun_alt, _) = self.altaz(jd_ut, rs, ds);
        if sun_alt >= 0.0 && !self.flags.dark {
            return Ok(f64::NEG_INFINITY);
        }
        let (ro, dobj, _) = self.place(self.target, jd_ut, topo, false)?;
        let (obj_alt, obj_az) = self.altaz(jd_ut, ro, dobj);
        if obj_alt <= 0.0 {
            return Ok(f64::NEG_INFINITY);
        }
        if let Some(e) = east {
            // azimuth from south, westwards: 180–360 is east of the meridian
            if (obj_az > 180.0) != e {
                return Ok(f64::NEG_INFINITY);
            }
        }
        let (rm, dm, _) = self.place(HelTarget::Body(Body::Moon), jd_ut, topo, false)?;
        let (moon_alt, _) = self.altaz(jd_ut, rm, dm);
        let (mag, moon_pa, month, year) = fixed;
        let mut g = SkyGeometry {
            z_obj: 90.0 - obj_alt,
            z_sun: 90.0 - sun_alt,
            z_moon: 90.0 - moon_alt,
            sep_sun: sep(ro, dobj, rs, ds),
            sep_moon: sep(ro, dobj, rm, dm),
            moon_phase_angle: moon_pa,
            lat_deg: self.obs.lat_deg,
            height_m: self.obs.height_m,
            month,
            year,
        };
        if self.target == HelTarget::Body(Body::Moon) {
            g.z_moon = 180.0;
        }
        Ok(schaefer_limit(&g, &self.atm, &self.observer, &self.flags).limiting_mag - mag)
    }

    fn fixed(&self, jd_ut: f64) -> Result<(f64, f64, f64, f64)> {
        let ymd = crate::time::revjul(jd_ut, crate::time::Calendar::Gregorian);
        Ok((
            self.magnitude(jd_ut)?,
            self.engine.pheno_ut(Body::Moon, jd_ut, &CalcOptions::default())?.phase_angle,
            ymd.1 as f64,
            crate::time::decimal_year(jd_ut),
        ))
    }

    /// Search window of a local morning or evening: (start, end) JD UT, from
    /// astronomical twilight (Sun −18°, or 4 h before sunrise/after sunset when
    /// it does not get that dark) to sunrise / from sunset.
    fn window(&self, day0: f64, morning: bool) -> Result<Option<(f64, f64)>> {
        let e = self.engine;
        let sun_ev = |kind, t0: f64, o: &RiseSetOptions| -> Result<Option<f64>> {
            Ok(match e.rise_set(Body::Sun, t0, kind, &self.obs, o)? {
                RiseSet::Event(t) => Some(t),
                _ => None,
            })
        };
        let plain = RiseSetOptions::default();
        let twl = RiseSetOptions {
            disc: DiscPoint::Center,
            refraction: None,
            horizon_deg: -18.0,
            ..Default::default()
        };
        if morning {
            // Sunrise of this local morning (search from local midnight).
            let Some(rise) = sun_ev(RiseSetKind::Rise, day0, &plain)? else { return Ok(None) };
            if rise > day0 + 0.75 {
                return Ok(None);
            }
            let start = match sun_ev(RiseSetKind::Rise, rise - 0.5, &twl)? {
                Some(t) if t < rise => t,
                _ => rise - 4.0 / 24.0,
            };
            Ok(Some((start.max(rise - 6.0 / 24.0), rise)))
        } else {
            let Some(set) = sun_ev(RiseSetKind::Set, day0 + 0.25, &plain)? else { return Ok(None) };
            if set > day0 + 1.0 {
                return Ok(None);
            }
            let end = match sun_ev(RiseSetKind::Set, set, &twl)? {
                Some(t) if t < set + 0.5 => t,
                _ => set + 4.0 / 24.0,
            };
            Ok(Some((set, end.min(set + 6.0 / 24.0))))
        }
    }

    /// Visibility on one local morning/evening: (first, best, last) JD UT.
    /// Quick yes/no visibility on one local morning/evening (12 samples).
    pub fn day_visible(&self, day0: f64, morning: bool) -> Result<bool> {
        let Some((t0, t1)) = self.window(day0, morning)? else { return Ok(false) };
        let fx = self.fixed(0.5 * (t0 + t1))?;
        const N: usize = 12;
        for i in 0..=N {
            if self.margin_side(t0 + (t1 - t0) * i as f64 / N as f64, fx, Some(morning))? > 0.0 {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub fn day_visibility(&self, day0: f64, morning: bool) -> Result<Option<[f64; 3]>> {
        let Some((t0, t1)) = self.window(day0, morning)? else { return Ok(None) };
        let mag = self.fixed(0.5 * (t0 + t1))?;
        const N: usize = 40;
        let ts: Vec<f64> = (0..=N).map(|i| t0 + (t1 - t0) * i as f64 / N as f64).collect();
        let mut m = Vec::with_capacity(N + 1);
        for &t in &ts {
            m.push(self.margin_side(t, mag, Some(morning))?);
        }
        let Some(ib) = (0..=N).filter(|&i| m[i] > 0.0).max_by(|&a, &b| m[a].total_cmp(&m[b])) else {
            return Ok(None);
        };
        let first = (0..=N).find(|&i| m[i] > 0.0).unwrap();
        let last = (0..=N).rev().find(|&i| m[i] > 0.0).unwrap();
        let refine = |a: f64, b: f64, rising: bool| -> Result<f64> {
            // Bisection on the sign of the margin between a (invisible) and b (visible) or reverse.
            let (mut lo, mut hi) = (a, b);
            for _ in 0..30 {
                let mid = 0.5 * (lo + hi);
                let v = self.margin_side(mid, mag, Some(morning))? > 0.0;
                if v == rising {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            Ok(0.5 * (lo + hi))
        };
        let tf = if first == 0 { ts[0] } else { refine(ts[first - 1], ts[first], true)? };
        let tl = if last == N { ts[N] } else { refine(ts[last], ts[last + 1], false)? };
        // Best instant: golden-section maximum of the margin around the best sample.
        let (mut a, mut b) = (ts[ib.saturating_sub(1)], ts[(ib + 1).min(N)]);
        let gr = 0.618_033_988_75;
        for _ in 0..30 {
            let c = b - gr * (b - a);
            let d = a + gr * (b - a);
            if self.margin_side(c, mag, Some(morning))? > self.margin_side(d, mag, Some(morning))? {
                b = d;
            } else {
                a = c;
            }
        }
        Ok(Some([tf, (0.5 * (a + b)).clamp(tf, tl), tl]))
    }
}

/// Classic `vis_limit_mag` result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VisLimit {
    /// −2 object below the horizon, 0 photopic, 1 scotopic.
    pub res: i32,
    /// [limiting mag, obj alt, obj az, sun alt, sun az, moon alt, moon az, obj mag, 0, 0]
    /// (azimuths from south, westwards, as the classic `azalt`).
    pub dret: [f64; 10],
}

/// Local midnight (JD UT) of the local civil day containing `jd_ut` (mean solar time).
fn local_midnight(jd_ut: f64, lon_deg: f64) -> f64 {
    let local = jd_ut + lon_deg / 360.0;
    (local - 0.5).floor() + 0.5 - lon_deg / 360.0
}

impl Engine {
    /// Limiting visual magnitude at the object's position (classic `vis_limit_mag`).
    pub fn vis_limit_mag(
        &self,
        jd_ut: f64,
        obs: &Observer,
        atm: &HelAtmosphere,
        observer: &HelObserver,
        target: HelTarget,
        flags: HelFlags,
    ) -> Result<VisLimit> {
        let cx = HelContext {
            engine: self,
            target,
            obs: *obs,
            atm: atm.resolved(obs.height_m),
            observer: observer.resolved(),
            flags,
        };
        let s = cx.sky(jd_ut, None)?;
        let mut dret = [0.0; 10];
        // Substituted bodies are reported at the nadir (classic semantics).
        let sun_alt = if flags.dark { -90.0 } else { s.sun_alt };
        let moon_alt = if flags.dark || flags.no_moon { -90.0 } else { s.moon_alt };
        dret[1..8].copy_from_slice(&[s.obj_alt, s.obj_az, sun_alt, s.sun_az, moon_alt, s.moon_az, s.obj_mag]);
        match cx.limit(&s) {
            None => Ok(VisLimit { res: -2, dret }),
            Some(v) => {
                dret[0] = v.limiting_mag;
                Ok(VisLimit { res: v.vision as i32, dret })
            }
        }
    }

    /// Next heliacal event after `jd_ut` (classic `heliacal_ut`): (first visible,
    /// best, last visible) instants of the event day, JD UT.
    pub fn heliacal_ut(
        &self,
        jd_ut: f64,
        obs: &Observer,
        atm: &HelAtmosphere,
        observer: &HelObserver,
        target: HelTarget,
        event: HelEvent,
        flags: HelFlags,
    ) -> Result<[f64; 3]> {
        target.check_event(event)?;
        let cx = HelContext {
            engine: self,
            target,
            obs: *obs,
            atm: atm.resolved(obs.height_m),
            observer: observer.resolved(),
            flags,
        };
        let morning = event.morning();
        let day0 = local_midnight(jd_ut, obs.lon_deg);
        // Coarse scan for the visibility transition, then day by day inside it.
        let step: i64 = match target {
            HelTarget::Body(Body::Moon) => 1,
            HelTarget::Body(Body::Mercury | Body::Venus) => 2,
            _ => 4,
        };
        const MAX_DAYS: i64 = 1200; // > 1.5 × the longest synodic period (Mars, 780 d)
        let vis = |d: i64| cx.day_visible(day0 + d as f64, morning);
        let mut prev = vis(0)?;
        let mut seen_other = if event.first() { !prev } else { prev };
        let mut k = 0;
        while k < MAX_DAYS {
            let cur = vis(k + step)?;
            let transition = if event.first() { !prev && cur } else { prev && !cur };
            if transition && seen_other {
                // Locate the exact day inside (k, k + step].
                let mut d = k + 1;
                let mut last_vis = if event.first() { None } else { Some(k) };
                while d <= k + step {
                    let v = if d == k + step { cur } else { vis(d)? };
                    if event.first() && v {
                        break;
                    }
                    if !event.first() {
                        if v {
                            last_vis = Some(d);
                        } else {
                            break;
                        }
                    }
                    d += 1;
                }
                let ev_day = if event.first() { d } else { last_vis.unwrap() };
                if let Some(r) = cx.day_visibility(day0 + ev_day as f64, morning)? {
                    if r[0] >= jd_ut {
                        return Ok(r);
                    }
                }
            }
            if event.first() && !cur || !event.first() && cur {
                seen_other = true;
            }
            prev = cur;
            k += step;
        }
        Err(Error::InvalidArgument(format!("heliacal event does not happen for {}", target.label())))
    }

    /// Heliacal phenomena at an instant (classic `heliacal_pheno_ut`, 50 values).
    /// Filled: 0 AltO, 1 AppAltO, 2 GeoAltO, 3 AziO, 4 AltS, 5 AziS, 6 TAVact,
    /// 7 ARCVact, 8 DAZact, 9 ARCLact, 10 kact, 11 minTAV, 12–14 first/best/last
    /// visible instants of that day (VR), 16 WMoon, 17 qYal, 18 qCrit (Moon),
    /// 19 ParO, 20 Magn, 21 RiseO, 22 RiseS, 23 Lag, 24 TvisVR, 25 LMoon,
    /// 27 Illum (%). Others 0.
    pub fn heliacal_pheno_ut(
        &self,
        jd_ut: f64,
        obs: &Observer,
        atm: &HelAtmosphere,
        observer: &HelObserver,
        target: HelTarget,
        event: HelEvent,
        flags: HelFlags,
    ) -> Result<[f64; 50]> {
        let cx = HelContext {
            engine: self,
            target,
            obs: *obs,
            atm: atm.resolved(obs.height_m),
            observer: observer.resolved(),
            flags,
        };
        let s = cx.sky(jd_ut, None)?;
        let mut d = [0.0; 50];
        d[0] = s.obj_alt;
        d[1] = s.obj_app_alt;
        d[2] = s.obj_geo_alt;
        d[3] = s.obj_az;
        d[4] = s.sun_alt;
        d[5] = s.sun_az;
        d[6] = s.obj_alt - s.sun_alt;
        d[7] = s.obj_geo_alt - s.sun_geo_alt;
        d[8] = diff_deg(s.obj_az, s.sun_az);
        d[9] = s.dlon;
        let (k, _) =
            extinction_v(s.geometry.z_obj.min(90.0), &cx.atm, obs.lat_deg, obs.height_m, s.geometry.month);
        d[10] = k;
        let day0 = local_midnight(jd_ut, obs.lon_deg);
        if let Some(v) = cx.day_visibility(day0, event.morning())? {
            d[12..15].copy_from_slice(&v);
            d[24] = v[2] - v[0];
            let sf = cx.sky(v[0], Some(s.obj_mag))?;
            d[11] = sf.obj_alt - sf.sun_alt;
        }
        d[19] = s.parallax;
        d[20] = s.obj_mag;
        // Rise (morning events) or set (evening) of object and Sun on that day.
        let kind = if event.morning() { RiseSetKind::Rise } else { RiseSetKind::Set };
        let t0 = if event.morning() { day0 } else { day0 + 0.25 };
        if let RiseSet::Event(ts) = self.rise_set(Body::Sun, t0, kind, obs, &RiseSetOptions::default())? {
            d[22] = ts;
            let to = match target {
                HelTarget::Body(b) => self.rise_set(b, ts - 0.5, kind, obs, &RiseSetOptions::default())?.jd(),
                HelTarget::Star(_) => None,
            };
            if let Some(to) = to {
                d[21] = to;
                d[23] = to - ts;
            }
        }
        if let HelTarget::Body(Body::Moon) = target {
            // Yallop (1997, NAO Technical Note 69) crescent width and q test.
            let p = self.pheno_ut(Body::Moon, jd_ut, &CalcOptions::default())?;
            let sd_arcmin = p.diameter / 2.0 * 60.0;
            let arcl = s.sep_sun;
            let w = sd_arcmin * (1.0 - (arcl * DEG).cos());
            let arcv = d[7];
            d[16] = w / 60.0;
            d[17] = (arcv - (11.8371 - 6.3226 * w + 0.7319 * w * w - 0.1018 * w * w * w)) / 10.0;
            d[18] = if d[17] > 0.216 {
                1.0
            } else if d[17] > -0.014 {
                2.0
            } else if d[17] > -0.160 {
                3.0
            } else if d[17] > -0.232 {
                4.0
            } else {
                5.0
            };
            d[25] = arcl;
            d[27] = p.phase * 100.0;
        } else if let HelTarget::Body(b) = target {
            d[27] = self.pheno_ut(b, jd_ut, &CalcOptions::default())?.phase * 100.0;
        }
        Ok(d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schaefer_dark_sky_zenith() {
        // Dark moonless sky, star at the zenith, sea level, 40 % RH, 15 °C:
        // Schaefer's program gives a limit of about 6.1–6.5.
        let atm =
            HelAtmosphere { pressure_hpa: 1013.25, temperature_c: 15.0, humidity_pct: 40.0, vr_or_k: 0.0 };
        let g = SkyGeometry {
            z_obj: 0.0,
            z_sun: 150.0,
            z_moon: 150.0,
            sep_sun: 150.0,
            sep_moon: 150.0,
            moon_phase_angle: 180.0,
            lat_deg: 30.0,
            height_m: 0.0,
            month: 3.0,
            year: 1992.0 + 5.5,
        };
        let v = schaefer_limit(&g, &atm, &HelObserver::default(), &HelFlags::default());
        assert!(v.limiting_mag > 5.8 && v.limiting_mag < 6.8, "{v:?}");
        assert_eq!(v.vision, Vision::Scotopic);
    }
}
