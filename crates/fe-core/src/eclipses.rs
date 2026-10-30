// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Solar and lunar eclipses computed directly from JPL Sun/Moon states.
//!
//! Method (see docs/ECLIPSES_REPORT.md for validation):
//!
//! * **Geometry.** Geocentric *astrometric* vectors (light time, no
//!   aberration) of the Sun and Moon in GCRS, from the pack (`Sky::at`).
//!   Annual aberration displaces both bodies by the same ~20″ as seen by one
//!   observer, so it cancels in every alignment condition; it is ignored for
//!   eclipse geometry (only horizon altitudes would see it, at the 0.006°
//!   level).
//! * **Candidates.** Mean lunations (Meeus, *Astronomical Algorithms* ch. 49)
//!   with the argument of latitude F; lunations with |sin F| > 0.5 (F more
//!   than 30° from a node, well outside the ≈ 21° eclipse limit even allowing
//!   for the mean-vs-true phase offset) are skipped without ephemeris work.
//! * **Lunar eclipses.** Shadow centre = anti-Sun direction. Greatest eclipse
//!   = minimum angular distance σ of the Moon from the shadow axis (Newton on
//!   the chord vector, analytic rates). Shadow radii use **Danjon's rule**
//!   (Earth radius enlarged by 1 %, as in the NASA Five Millennium Canon; the
//!   Swiss oracle measures to the same model, see the report):
//!   ρ_umbra = 1.01·π_☾ + π_☉ − s_☉, ρ_penumbra = 1.01·π_☾ + π_☉ + s_☉,
//!   with the Moon radius k = 0.2725076 Earth radii.
//!   Magnitudes: (ρ − σ + s_☾) / (2 s_☾). Contacts: σ = ρ ± s_☾ (Brent).
//! * **Solar eclipses.** Besselian elements (Explanatory Supplement 3rd ed.
//!   ch. 11; Meeus, *Elements of Solar Eclipses*) evaluated *directly* at each
//!   instant from the vectors instead of from fitted polynomials: shadow axis
//!   ẑ = (☉ − ☾)/|☉ − ☾|, fundamental plane (x̂ east, ŷ north w.r.t. the true
//!   pole of date), x, y, z of the Moon in Earth equatorial radii, cone angles
//!   f1, f2 and radii l1, l2. Penumbra uses k1 = 0.272488, umbra
//!   k2 = 0.2722810 (the constants printed on NASA's Besselian-element pages). Greatest eclipse = minimum of √(x²+y²)
//!   (axis closest to the geocentre, NASA's definition), gamma = that
//!   distance signed by y. Earth's silhouette on the fundamental plane is the
//!   ellipse of semi-axes 1 and ρ1 = √(1 − e² cos² d) (WGS84), and contacts
//!   use the exact point–ellipse distance. Central-line points intersect the
//!   axis with the WGS84 ellipsoid.
//! * **Local circumstances.** Topocentric (WGS84 site) separation and
//!   semidiameters; contacts are roots of sep − (s☉ ± s☾).
//!
//! Constants: Earth equatorial radius 6378.137 km and flattening
//! 1/298.257223563 (WGS84 = IERS/IAU 2015 nominal values), Sun radius
//! 696 000 km (959.63″ at 1 au, the NASA canon value; IAU 2015 nominal is
//! 695 700 km), Moon k-values as above (lunar eclipses: k = 0.2725076).
//!
//! All core functions take and return **TT** Julian dates; the `*_ut`
//! helpers convert with the engine's Delta T model.

use crate::context::Engine;
use crate::ephem::{Observer, EARTH, MOON, SSB, SUN};
use crate::frames::FrameSet;
use crate::math::*;
use crate::pack::Pack;
use crate::search::brent;
use crate::time;
use crate::{Error, Result};

mod saros;

/// Earth equatorial radius, km (WGS84; IAU 2015 nominal).
pub const EARTH_EQ_RADIUS_KM: f64 = 6378.137;
/// Earth flattening (WGS84).
pub const EARTH_FLATTENING: f64 = 1.0 / 298.257_223_563;
/// Solar radius, km (959.63″ at 1 au; NASA canon).
pub const SUN_RADIUS_KM: f64 = 696_000.0;
/// Moon radius in Earth radii for lunar eclipses (NASA lunar canon, mean radius).
pub const K_MOON_MEAN: f64 = 0.272_507_6;
/// Moon radius in Earth radii for the solar penumbra / partial phases (the
/// k1 printed on NASA's solar Besselian-element pages).
pub const K_MOON_PENUMBRA: f64 = 0.272_488;
/// Moon radius in Earth radii for umbral (total/annular) phases (NASA, IAU 1982).
pub const K_MOON_UMBRA: f64 = 0.272_281_0;
/// Danjon enlargement of the Earth's radius for the lunar-eclipse shadow.
pub const DANJON_FACTOR: f64 = 1.01;
/// "Next" eclipse means greatest eclipse more than this after the start (days);
/// measured on the Swiss oracle (`*_when*` skip an eclipse up to 1e-4 d past its maximum).
pub const SEARCH_MARGIN_DAYS: f64 = 1e-4;

pub(crate) const RE: f64 = EARTH_EQ_RADIUS_KM;
const E2: f64 = EARTH_FLATTENING * (2.0 - EARTH_FLATTENING);
const SYNODIC: f64 = 29.530_588_861;
/// Earth rotation rate (rad/day of UT1): dERA/dt.
const OMEGA_UT: f64 = TAU * 1.002_737_811_911_354_6;
/// Max lunations scanned by one search (≈ 8 years: there is at least one
/// lunar and one solar eclipse every year).
const MAX_LUNATIONS: usize = 100;
/// |sin F| above which a lunation cannot produce an eclipse.
const NODE_FILTER: f64 = 0.5;

// ---------------------------------------------------------------------------------------------
// Shared geometry
// ---------------------------------------------------------------------------------------------

/// Geocentric astrometric occulted body (`s`: the Sun for eclipses) and Moon
/// (GCRS, km, km/day).
#[derive(Debug, Clone, Copy)]
pub(crate) struct SunMoon {
    pub(crate) s: Vec3,
    pub(crate) vs: Vec3,
    pub(crate) m: Vec3,
    pub(crate) vm: Vec3,
}

/// Size of the occulted body and the Moon radii (Earth radii) used for the
/// outer (k1: first/last contact) and inner (k2: second/third) contacts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Shape {
    pub(crate) radius_km: f64,
    pub(crate) k1: f64,
    pub(crate) k2: f64,
}

pub(crate) const SUN_SHAPE: Shape = Shape { radius_km: SUN_RADIUS_KM, k1: K_MOON_PENUMBRA, k2: K_MOON_UMBRA };

/// Distance (km) at which a star's unit direction is placed so that the
/// eclipse geometry applies unchanged (parallax of a 6 378 km offset:
/// 0.013 mas).
pub(crate) const STAR_DISTANCE_KM: f64 = 1e14;

/// The body behind the Moon.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Target {
    Sun,
    /// Pack body (NAIF id), light-time corrected.
    Body(i32),
    /// Catalogue star (space motion and parallax applied).
    Star(&'static crate::stars::StarRecord),
}

/// Target/Moon evaluator for one event (TDB − TT frozen at the candidate:
/// its ≤ 1.7 ms variation moves the Moon by ≤ 2 m).
pub(crate) struct Sky<'a> {
    pack: &'a Pack,
    dtdb: f64,
    target: Target,
    pub(crate) shape: Shape,
}

impl<'a> Sky<'a> {
    pub(crate) fn new(pack: &'a Pack, jd_tt: f64) -> Sky<'a> {
        Sky::with_target(pack, jd_tt, Target::Sun, SUN_SHAPE)
    }

    pub(crate) fn with_target(pack: &'a Pack, jd_tt: f64, target: Target, shape: Shape) -> Sky<'a> {
        Sky { pack, dtdb: time::tdb_minus_tt_seconds(jd_tt) / 86400.0, target, shape }
    }

    /// Earth barycentric, target, Moon geocentric (three or four pack
    /// evaluations). Light time is applied linearly (Sun: 500 s × 0.01 km/s;
    /// Moon: 1.3 s; planets: one re-evaluation at t − τ), leaving < 1e-4 km.
    pub(crate) fn at(&self, jd_tt: f64) -> Result<SunMoon> {
        let tdb = jd_tt + self.dtdb;
        let (e, ve) = self.pack.state(EARTH, SSB, tdb)?;
        let (s, vs) = match self.target {
            Target::Sun => {
                let (sb, vsb) = self.pack.state(SUN, SSB, tdb)?;
                let tau_s = norm(sub(sb, e)) / C_KM_DAY;
                (sub(sub(sb, scale(vsb, tau_s)), e), sub(vsb, ve))
            }
            Target::Body(id) => {
                let p0 = self.pack.position(id, SSB, tdb)?;
                let tau0 = norm(sub(p0, e)) / C_KM_DAY;
                let (p1, v1) = self.pack.state(id, SSB, tdb - tau0)?;
                let tau1 = norm(sub(p1, e)) / C_KM_DAY;
                (sub(sub(p1, scale(v1, tau1 - tau0)), e), sub(v1, ve))
            }
            Target::Star(rec) => {
                let p = scale(crate::stars::barycentric_au(rec, tdb), AU_KM);
                let d = sub(p, e);
                (scale(d, STAR_DISTANCE_KM / norm(d)), [0.0; 3])
            }
        };
        let (m0, vm) = self.pack.state(MOON, EARTH, tdb)?;
        let tau_m = norm(m0) / C_KM_DAY;
        // Moon_bary(t−τ) − Earth(t) = (Moon − Earth)(t−τ) + Earth(t−τ) − Earth(t).
        let m = sub(m0, scale(add(vm, ve), tau_m));
        Ok(SunMoon { s, vs, m, vm })
    }

    /// Lunar-eclipse geometry for an observation (geocentric) at `jd_tt`.
    /// The Moon is seen at t = T − τ☾; the shadow it is in at t was cast by
    /// sunlight passing the Earth at t − τ☾, so (physically)
    /// `m` = Moon(t) − Earth(t − τ☾) and the shadow axis is
    /// Earth(t − τ☾) − Sun(t − τ☾ − τ☉); `s` is returned as minus that axis.
    /// To < 1″ this equals the classical "apparent Moon vs apparent Sun + 180°"
    /// convention of the NASA canon, whereas astrometric vectors (light time
    /// only) would offset the Moon by ≈ 42″ (≈ 80 s in time).
    /// `vs`/`vm` are geometric rates (only used as Newton derivatives).
    fn lunar(&self, jd_tt: f64) -> Result<SunMoon> {
        let tdb = jd_tt + self.dtdb;
        let (e, ve) = self.pack.state(EARTH, SSB, tdb)?;
        let (sb, vsb) = self.pack.state(SUN, SSB, tdb)?;
        let (m0, vm) = self.pack.state(MOON, EARTH, tdb)?;
        let tau_m = norm(m0) / C_KM_DAY;
        let tau_s = norm(sub(sb, e)) / C_KM_DAY;
        // MoonGeo(t) + Earth(t) − Earth(t − τ☾).
        let m = add(sub(m0, scale(vm, tau_m)), scale(ve, tau_m));
        let earth = sub(e, scale(ve, 2.0 * tau_m));
        let sun = sub(sb, scale(vsb, 2.0 * tau_m + tau_s));
        Ok(SunMoon { s: sub(sun, earth), vs: sub(vsb, ve), m, vm })
    }
}

#[inline]
pub(crate) fn unit_rate(r: Vec3, v: Vec3) -> (Vec3, Vec3) {
    let n = norm(r);
    let u = scale(r, 1.0 / n);
    (u, scale(sub(v, scale(u, dot(u, v))), 1.0 / n))
}

#[inline]
pub(crate) fn angle(a: Vec3, b: Vec3) -> f64 {
    // atan2(|a×b|, a·b): accurate at small angles.
    norm(cross(a, b)).atan2(dot(a, b))
}

/// One Gauss–Newton step towards the minimum of |p(t)|.
pub(crate) fn newton_min(mut eval: impl FnMut(f64) -> Result<(Vec3, Vec3)>, t0: f64) -> Result<f64> {
    let mut t = t0;
    for _ in 0..20 {
        let (p, dp) = eval(t)?;
        let dt = (-dot(p, dp) / dot(dp, dp)).clamp(-1.0, 1.0);
        t += dt;
        if dt.abs() < 1e-9 {
            break;
        }
    }
    Ok(t)
}

/// Root of `g` between `a` and `b` (either order), if it changes sign.
pub(crate) fn root(mut g: impl FnMut(f64) -> Result<f64>, a: f64, b: f64) -> Result<Option<f64>> {
    let (fa, fb) = (g(a)?, g(b)?);
    if fa.signum() == fb.signum() && fa != 0.0 && fb != 0.0 {
        return Ok(None);
    }
    let (lo, hi, flo, fhi) = if a < b { (a, b, fa, fb) } else { (b, a, fb, fa) };
    Ok(Some(brent(&mut g, lo, hi, flo, fhi)?))
}

/// Mean lunation `k` (integer: new moon, +0.5: full moon) → (JDE, F rad). Meeus ch. 49.
fn mean_phase(k: f64) -> (f64, f64) {
    let t = k / 1236.85;
    let jde =
        2_451_550.097_66 + SYNODIC * k + t * t * (0.000_154_37 + t * (-0.000_000_150 + t * 0.000_000_000_73));
    let f = 160.7108 + 390.670_502_84 * k + t * t * (-0.001_611_8 + t * (-0.000_002_27 + t * 0.000_000_011));
    (jde, f * DEG)
}

fn first_lunation(jd_tt: f64, offset: f64, backward: bool) -> f64 {
    let k = ((jd_tt - 2_451_550.097_66) / SYNODIC).floor() + offset;
    if backward {
        k + 2.0
    } else {
        k - 1.0
    }
}

pub(crate) fn frames_at(engine: &Engine, jd_tt: f64) -> FrameSet {
    FrameSet::with_model(jd_tt, engine.frame_model())
}

/// True celestial pole of date (GCRS unit vector): third row of NPB.
pub(crate) fn pole(fr: &FrameSet) -> Vec3 {
    fr.bpn[2]
}

/// GCRS position (km) → (geodetic lon, lat) in degrees.
pub(crate) fn geodetic_deg(fr: &FrameSet, jd_ut: f64, r: Vec3) -> (f64, f64) {
    let itrs = mat_vec(&rot_z(fr.gast(jd_ut)), mat_vec(&fr.bpn, r));
    let lon = itrs[1].atan2(itrs[0]);
    let rho = itrs[0].hypot(itrs[1]);
    let lat = itrs[2].atan2((1.0 - EARTH_FLATTENING).powi(2) * rho);
    (lon / DEG, lat / DEG)
}

// ---------------------------------------------------------------------------------------------
// Observer on the ellipsoid
// ---------------------------------------------------------------------------------------------

/// Observer prepared for one eclipse: fixed frames (NPB changes < 0.01″/h)
/// and Earth rotation from ERA (exact).
pub(crate) struct Site {
    itrs: Vec3,
    up: Vec3,
    east: Vec3,
    north: Vec3,
    pub(crate) height_m: f64,
}

impl Site {
    pub(crate) fn new(obs: &Observer) -> Site {
        let (sp, cp) = (obs.lat_deg * DEG).sin_cos();
        let (sl, cl) = (obs.lon_deg * DEG).sin_cos();
        Site {
            itrs: obs.itrs_km(),
            up: [cp * cl, cp * sl, sp],
            east: [-sl, cl, 0.0],
            north: [-sp * cl, -sp * sl, cp],
            height_m: obs.height_m,
        }
    }

    /// GCRS position and velocity (km, km/day) at JD UT.
    pub(crate) fn gcrs(&self, fr: &FrameSet, jd_ut: f64) -> (Vec3, Vec3) {
        let r_tod = mat_vec(&rot_z(-fr.gast(jd_ut)), self.itrs);
        let v_tod = [-OMEGA_UT * r_tod[1], OMEGA_UT * r_tod[0], 0.0];
        (mat_t_vec(&fr.bpn, r_tod), mat_t_vec(&fr.bpn, v_tod))
    }

    /// (true altitude, azimuth from north through east), radians, of a GCRS direction.
    pub(crate) fn alt_az(&self, fr: &FrameSet, jd_ut: f64, v: Vec3) -> (f64, f64) {
        let it = mat_vec(&rot_z(fr.gast(jd_ut)), mat_vec(&fr.bpn, v));
        let n = norm(it);
        let alt = (dot(it, self.up) / n).clamp(-1.0, 1.0).asin();
        let az = norm_rad(dot(it, self.east).atan2(dot(it, self.north)));
        (alt, az)
    }
}

/// Apparent altitude (deg) from true altitude (deg) at a site of height
/// `height_m`: Sæmundsson's refraction formula (Meeus, Astronomical Algorithms
/// eq. 16.4; 10 °C) scaled to the standard-atmosphere pressure at that height
/// (`events::standard_pressure_hpa`). No refraction is applied while the
/// refracted body would still be below the (dip-free) horizon, so below about
/// −0.6° (sea level) the apparent altitude equals the true one. This is the
/// altitude the classic `*_eclipse_how` attributes report and test.
pub fn apparent_altitude_deg(true_alt_deg: f64, height_m: f64) -> f64 {
    let h = true_alt_deg;
    if h < -2.0 {
        return h;
    }
    let p = crate::events::standard_pressure_hpa(height_m);
    let r_arcmin = 1.02 / ((h + 10.3 / (h + 5.11)) * DEG).tan() * (p / 1010.0);
    let app = h + r_arcmin / 60.0;
    if app < 0.0 && h < 0.0 {
        h
    } else {
        app
    }
}

// ---------------------------------------------------------------------------------------------
// Lunar eclipses
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LunarKind {
    Penumbral,
    Partial,
    Total,
}

/// Instantaneous lunar-eclipse geometry (angles in radians, magnitudes unclamped).
#[derive(Debug, Clone, Copy)]
pub struct LunarShadow {
    /// Angular distance of the Moon's centre from the shadow axis.
    pub sigma: f64,
    pub rho_umbra: f64,
    pub rho_penumbra: f64,
    pub moon_sd: f64,
    pub umbral_mag: f64,
    pub penumbral_mag: f64,
}

fn lunar_shadow(g: &SunMoon) -> LunarShadow {
    let dm = norm(g.m);
    let ds = norm(g.s);
    let pi_m = (RE / dm).asin();
    let pi_s = (RE / ds).asin();
    let s_s = (SUN_RADIUS_KM / ds).asin();
    let s_m = (K_MOON_MEAN * RE / dm).asin();
    let sigma = angle(g.m, scale(g.s, -1.0));
    let rho_umbra = DANJON_FACTOR * pi_m + pi_s - s_s;
    let rho_penumbra = DANJON_FACTOR * pi_m + pi_s + s_s;
    LunarShadow {
        sigma,
        rho_umbra,
        rho_penumbra,
        moon_sd: s_m,
        umbral_mag: (rho_umbra - sigma + s_m) / (2.0 * s_m),
        penumbral_mag: (rho_penumbra - sigma + s_m) / (2.0 * s_m),
    }
}

/// Lunar eclipse; all times JD TT.
#[derive(Debug, Clone, Copy)]
pub struct LunarEclipse {
    pub kind: LunarKind,
    pub max: f64,
    pub umbral_mag: f64,
    pub penumbral_mag: f64,
    /// Moon–axis distance at maximum, Earth equatorial radii, + when the Moon passes north.
    pub gamma: f64,
    pub p1: f64,
    pub u1: Option<f64>,
    pub u2: Option<f64>,
    pub u3: Option<f64>,
    pub u4: Option<f64>,
    pub p4: f64,
}

/// Lunar-eclipse geometry at an instant (JD TT).
pub fn lunar_shadow_at(engine: &Engine, jd_tt: f64) -> Result<LunarShadow> {
    Ok(lunar_shadow(&Sky::new(engine.pack(), jd_tt).lunar(jd_tt)?))
}

fn lunar_at_lunation(engine: &Engine, k: f64, accept: impl Fn(f64) -> bool) -> Result<Option<LunarEclipse>> {
    let (jde, f) = mean_phase(k);
    if f.sin().abs() > NODE_FILTER {
        return Ok(None);
    }
    let sky = Sky::new(engine.pack(), jde);
    let max = newton_min(
        |t| {
            let g = sky.lunar(t)?;
            let (mu, dmu) = unit_rate(g.m, g.vm);
            let (su, dsu) = unit_rate(g.s, g.vs);
            Ok((add(mu, su), add(dmu, dsu)))
        },
        jde,
    )?;
    let g = sky.lunar(max)?;
    let sh = lunar_shadow(&g);
    if sh.penumbral_mag <= 0.0 || !accept(max) {
        return Ok(None);
    }
    let fr = frames_at(engine, max);
    let chord = add(scale(g.m, 1.0 / norm(g.m)), scale(g.s, 1.0 / norm(g.s)));
    let gamma = norm(g.m) * sh.sigma.sin() / RE * dot(chord, pole(&fr)).signum();

    let contact = |sign: f64, which: fn(&LunarShadow) -> f64| -> Result<Option<f64>> {
        root(
            |t| {
                let s = lunar_shadow(&sky.lunar(t)?);
                Ok(s.sigma - which(&s))
            },
            max,
            max + sign * 0.3,
        )
    };
    let pen = |s: &LunarShadow| s.rho_penumbra + s.moon_sd;
    let um_out = |s: &LunarShadow| s.rho_umbra + s.moon_sd;
    let um_in = |s: &LunarShadow| s.rho_umbra - s.moon_sd;
    let p1 = contact(-1.0, pen)?.unwrap_or(max);
    let p4 = contact(1.0, pen)?.unwrap_or(max);
    let (mut u1, mut u2, mut u3, mut u4) = (None, None, None, None);
    let mut kind = LunarKind::Penumbral;
    if sh.umbral_mag > 0.0 {
        kind = LunarKind::Partial;
        u1 = contact(-1.0, um_out)?;
        u4 = contact(1.0, um_out)?;
        if sh.umbral_mag >= 1.0 {
            kind = LunarKind::Total;
            u2 = contact(-1.0, um_in)?;
            u3 = contact(1.0, um_in)?;
        }
    }
    Ok(Some(LunarEclipse {
        kind,
        max,
        umbral_mag: sh.umbral_mag,
        penumbral_mag: sh.penumbral_mag,
        gamma,
        p1,
        u1,
        u2,
        u3,
        u4,
        p4,
    }))
}

/// Next (or previous) lunar eclipse whose greatest eclipse is more than
/// `SEARCH_MARGIN_DAYS` after (before) `jd_tt`.
pub fn lunar_eclipse_search(engine: &Engine, jd_tt: f64, backward: bool) -> Result<LunarEclipse> {
    let mut k = first_lunation(jd_tt, 0.5, backward);
    let step = if backward { -1.0 } else { 1.0 };
    for _ in 0..MAX_LUNATIONS {
        let accept =
            |t: f64| if backward { t < jd_tt - SEARCH_MARGIN_DAYS } else { t > jd_tt + SEARCH_MARGIN_DAYS };
        if let Some(e) = lunar_at_lunation(engine, k, accept)? {
            return Ok(e);
        }
        k += step;
    }
    Err(Error::InvalidArgument("no lunar eclipse found".into()))
}

/// Local lunar-eclipse attributes at an instant: magnitudes plus the Moon's
/// topocentric horizontal coordinates.
#[derive(Debug, Clone, Copy)]
pub struct LunarLocal {
    pub shadow: LunarShadow,
    /// Azimuth from north through east, degrees.
    pub azimuth_deg: f64,
    pub true_alt_deg: f64,
    pub apparent_alt_deg: f64,
    /// 180° − elongation of the Moon (distance from opposition), degrees.
    pub opposition_dist_deg: f64,
}

/// Type of lunar eclipse seen at a place at one instant, as the classic
/// `lun_eclipse_how` flags report it: `None` when there is no eclipse or the
/// Moon's centre is not above the apparent horizon (refracted altitude ≤ 0).
pub fn lunar_seen_kind(l: &LunarLocal) -> Option<LunarKind> {
    let (um, pm) = (l.shadow.umbral_mag, l.shadow.penumbral_mag);
    if pm <= 0.0 || l.apparent_alt_deg <= 0.0 {
        None
    } else if um >= 1.0 {
        Some(LunarKind::Total)
    } else if um > 0.0 {
        Some(LunarKind::Partial)
    } else {
        Some(LunarKind::Penumbral)
    }
}

pub fn lunar_local(engine: &Engine, jd_tt: f64, obs: &Observer) -> Result<LunarLocal> {
    let sky = Sky::new(engine.pack(), jd_tt);
    let g = sky.lunar(jd_tt)?;
    let fr = frames_at(engine, jd_tt);
    let jd_ut = engine.tt_to_ut(jd_tt);
    let site = Site::new(obs);
    let (o, _) = site.gcrs(&fr, jd_ut);
    // Apparent direction ≈ Moon(t) − Earth(t) (light time + aberration).
    let tau_m = norm(g.m) / C_KM_DAY;
    let (_, ve) = engine.pack().state(EARTH, SSB, jd_tt + sky.dtdb)?;
    let m_app = sub(g.m, scale(ve, tau_m));
    let (alt, az) = site.alt_az(&fr, jd_ut, sub(m_app, o));
    let ta = alt / DEG;
    Ok(LunarLocal {
        shadow: lunar_shadow(&g),
        azimuth_deg: az / DEG,
        true_alt_deg: ta,
        apparent_alt_deg: apparent_altitude_deg(ta, obs.height_m),
        opposition_dist_deg: 180.0 - angle(g.m, g.s) / DEG,
    })
}

// ---------------------------------------------------------------------------------------------
// Solar eclipses: Besselian elements
// ---------------------------------------------------------------------------------------------

/// Besselian elements at one instant (lengths in Earth equatorial radii).
#[derive(Debug, Clone, Copy)]
pub struct Besselian {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Declination of the shadow axis (towards the Sun), rad.
    pub d: f64,
    pub l1: f64,
    pub l2: f64,
    pub tan_f1: f64,
    pub tan_f2: f64,
    /// Semi-minor axis of the Earth's silhouette on the fundamental plane.
    pub rho1: f64,
    pub(crate) zh: Vec3,
    pub(crate) xh: Vec3,
    pub(crate) yh: Vec3,
    /// Moon projected on the fundamental plane, GCRS km.
    pub(crate) p0: Vec3,
}

pub(crate) fn besselian(g: &SunMoon, k: Vec3, sh: &Shape) -> Besselian {
    let gv = sub(g.s, g.m);
    let gn = norm(gv);
    let zh = scale(gv, 1.0 / gn);
    let xc = cross(k, zh);
    let xh = scale(xc, 1.0 / norm(xc));
    let yh = cross(zh, xh);
    let z = dot(g.m, zh) / RE;
    let sin_f1 = (sh.radius_km + sh.k1 * RE) / gn;
    let sin_f2 = (sh.radius_km - sh.k2 * RE) / gn;
    let cos_f1 = (1.0 - sin_f1 * sin_f1).sqrt();
    let cos_f2 = (1.0 - sin_f2 * sin_f2).sqrt();
    let (tan_f1, tan_f2) = (sin_f1 / cos_f1, sin_f2 / cos_f2);
    let sin_d = dot(k, zh);
    let cos2_d = 1.0 - sin_d * sin_d;
    Besselian {
        x: dot(g.m, xh) / RE,
        y: dot(g.m, yh) / RE,
        z,
        d: sin_d.asin(),
        l1: z * tan_f1 + sh.k1 / cos_f1,
        l2: z * tan_f2 - sh.k2 / cos_f2,
        tan_f1,
        tan_f2,
        rho1: (1.0 - E2 * cos2_d).sqrt(),
        zh,
        xh,
        yh,
        p0: sub(g.m, scale(zh, z * RE)),
    }
}

/// Moon offset from the shadow axis through the geocentre and its rate.
pub(crate) fn axis_offset(g: &SunMoon) -> (Vec3, Vec3) {
    let (zh, dz) = unit_rate(sub(g.s, g.m), sub(g.vs, g.vm));
    let mz = dot(g.m, zh);
    let dmz = dot(g.vm, zh) + dot(g.m, dz);
    let p = sub(g.m, scale(zh, mz));
    let dp = sub(sub(g.vm, scale(zh, dmz)), scale(dz, mz));
    (p, dp)
}

/// Signed distance from (x, y) to the ellipse X² + (Y/b)² = 1 (negative inside).
fn ellipse_dist(x: f64, y: f64, b: f64) -> f64 {
    let inside = x * x + (y / b) * (y / b) < 1.0;
    let mut th = (y / b).atan2(x);
    let bb = b * b - 1.0;
    for _ in 0..6 {
        let (s, c) = th.sin_cos();
        let f = x * s - b * y * c + bb * s * c;
        let fp = x * c + b * y * s + bb * (c * c - s * s);
        let step = f / fp;
        th -= step;
        if step.abs() < 1e-13 {
            break;
        }
    }
    let (s, c) = th.sin_cos();
    let d = (x - c).hypot(y - b * s);
    if inside {
        -d
    } else {
        d
    }
}

impl Besselian {
    /// Distance of the axis outside Earth's silhouette (negative: axis hits Earth).
    fn edge(&self) -> f64 {
        ellipse_dist(self.x, self.y, self.rho1)
    }
}

/// Intersection of the shadow axis with the WGS84 ellipsoid on the sunward
/// side: (GCRS point km, ζ in Earth radii).
pub(crate) fn axis_surface(b: &Besselian, k: Vec3) -> Option<(Vec3, f64)> {
    line_surface(b.p0, b.zh, k)
}

/// Sunward intersection of the line p0 + s·ẑ (p0 on the fundamental plane,
/// GCRS km) with the WGS84 ellipsoid: (point, ζ in Earth radii).
fn line_surface(p0: Vec3, zh: Vec3, k: Vec3) -> Option<(Vec3, f64)> {
    let c = 1.0 / (1.0 - EARTH_FLATTENING).powi(2) - 1.0;
    let zk = dot(zh, k);
    let pk = dot(p0, k);
    let qa = 1.0 + c * zk * zk;
    let qb = 2.0 * (dot(p0, zh) + c * pk * zk);
    let qc = dot(p0, p0) + c * pk * pk - RE * RE;
    let disc = qb * qb - 4.0 * qa * qc;
    if disc < 0.0 {
        return None;
    }
    let s = (-qb + disc.sqrt()) / (2.0 * qa);
    Some((add(p0, scale(zh, s)), s / RE))
}

/// Point of the ellipsoid closest to the axis when the axis misses Earth.
pub(crate) fn limb_point(b: &Besselian, k: Vec3) -> Vec3 {
    let st = 1.0 / (1.0 - EARTH_FLATTENING) - 1.0;
    let stretch = |v: Vec3| add(v, scale(k, st * dot(v, k)));
    let q0 = stretch(b.p0);
    let w = stretch(b.zh);
    let w = scale(w, 1.0 / norm(w));
    let foot = sub(q0, scale(w, dot(q0, w)));
    let q = scale(foot, RE / norm(foot));
    add(q, scale(k, -EARTH_FLATTENING * dot(q, k)))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolarKind {
    Partial,
    Annular,
    Total,
    Hybrid,
}

/// Global solar eclipse; all times JD TT.
#[derive(Debug, Clone, Copy)]
pub struct SolarEclipse {
    pub kind: SolarKind,
    /// Shadow axis meets the Earth.
    pub central: bool,
    /// Umbra/antumbra touches Earth but the axis misses it.
    pub noncentral: bool,
    pub max: f64,
    pub gamma: f64,
    /// NASA-style magnitude at greatest eclipse (Moon/Sun diameter ratio on the
    /// central line; fraction of the diameter covered for partial eclipses).
    pub magnitude: f64,
    /// Penumbra first / last touches Earth.
    pub p1: f64,
    pub p4: f64,
    /// Umbra (or antumbra) first / last touches Earth.
    pub u1: Option<f64>,
    pub u4: Option<f64>,
    /// Central line begins / ends.
    pub c1: Option<f64>,
    pub c2: Option<f64>,
    /// Hybrid: central line turns total / turns annular again.
    pub h1: Option<f64>,
    pub h2: Option<f64>,
}

/// Classical Besselian elements at an instant (Explanatory Supplement
/// notation; lengths in Earth equatorial radii, angles in degrees).
#[derive(Debug, Clone, Copy)]
pub struct BesselianElements {
    pub x: f64,
    pub y: f64,
    pub d_deg: f64,
    pub l1: f64,
    pub l2: f64,
    /// Greenwich hour angle of the shadow axis (true equinox, from `jd_ut`).
    pub mu_deg: f64,
    pub tan_f1: f64,
    pub tan_f2: f64,
}

/// Besselian elements at `jd_tt`, with `jd_ut` for the Earth's rotation (μ).
pub fn besselian_elements(engine: &Engine, jd_tt: f64, jd_ut: f64) -> Result<BesselianElements> {
    let g = Sky::new(engine.pack(), jd_tt).at(jd_tt)?;
    let fr = frames_at(engine, jd_tt);
    let b = besselian(&g, pole(&fr), &SUN_SHAPE);
    let z = mat_vec(&fr.bpn, b.zh);
    let ra = z[1].atan2(z[0]);
    Ok(BesselianElements {
        x: b.x,
        y: b.y,
        d_deg: b.d / DEG,
        l1: b.l1,
        l2: b.l2,
        mu_deg: norm_rad(fr.gast(jd_ut) - ra) / DEG,
        tan_f1: b.tan_f1,
        tan_f2: b.tan_f2,
    })
}

fn solar_at_lunation(engine: &Engine, k: f64, accept: impl Fn(f64) -> bool) -> Result<Option<SolarEclipse>> {
    let (jde, f) = mean_phase(k);
    if f.sin().abs() > NODE_FILTER {
        return Ok(None);
    }
    shadow_event(engine, &Sky::new(engine.pack(), jde), jde, accept)
}

/// Global event (eclipse or occultation) near `t_guess`: greatest eclipse,
/// contacts, type. `None` if the penumbra misses the Earth or `accept(max)`
/// is false.
pub(crate) fn shadow_event(
    engine: &Engine,
    sky: &Sky,
    t_guess: f64,
    accept: impl Fn(f64) -> bool,
) -> Result<Option<SolarEclipse>> {
    let max = newton_min(|t| Ok(axis_offset(&sky.at(t)?)), t_guess)?;
    let fr = frames_at(engine, max);
    let kp = pole(&fr);
    let bes = |t: f64| -> Result<Besselian> { Ok(besselian(&sky.at(t)?, kp, &sky.shape)) };
    let b = bes(max)?;
    let edge = b.edge();
    if edge >= b.l1 || !accept(max) {
        return Ok(None);
    }
    let gamma = b.x.hypot(b.y) * b.y.signum();
    let span = 0.35;
    let side = |g: &dyn Fn(&Besselian) -> f64, dir: f64| -> Result<Option<f64>> {
        root(|t| Ok(g(&bes(t)?)), max, max + dir * span)
    };
    let pen = |b: &Besselian| b.edge() - b.l1;
    let p1 = side(&pen, -1.0)?.unwrap_or(max);
    let p4 = side(&pen, 1.0)?.unwrap_or(max);
    let (mut u1, mut u4, mut c1, mut c2, mut h1, mut h2) = (None, None, None, None, None, None);
    let central = edge < 0.0;
    let umbral = edge < b.l2.abs();
    let mut kind = SolarKind::Partial;
    if umbral {
        let umb = |b: &Besselian| b.edge() - b.l2.abs();
        u1 = side(&umb, -1.0)?;
        u4 = side(&umb, 1.0)?;
        kind = if b.l2 < 0.0 { SolarKind::Total } else { SolarKind::Annular };
    }
    if central {
        let cen = |b: &Besselian| b.edge();
        c1 = side(&cen, -1.0)?;
        c2 = side(&cen, 1.0)?;
        // L2 on the central line: negative = total there.
        let l2c = |t: f64| -> Result<f64> {
            let b = bes(t)?;
            let zeta = axis_surface(&b, kp).map_or(0.0, |(_, z)| z);
            Ok(b.l2 - zeta * b.tan_f2)
        };
        let at_max = l2c(max)?;
        if at_max < 0.0 {
            kind = SolarKind::Total;
            if let Some(a) = c1 {
                h1 = root(&l2c, a, max)?;
            }
            if let Some(e) = c2 {
                h2 = root(&l2c, max, e)?;
            }
            if h1.is_some() || h2.is_some() {
                kind = SolarKind::Hybrid;
            }
        } else {
            kind = SolarKind::Annular;
        }
    }
    let wh = solar_where_impl(sky, &fr, max, engine.tt_to_ut(max))?;
    // NASA convention: diameter ratio on the central line, else the fraction
    // of the diameter covered at the point nearest the axis (also for
    // non-central total/annular eclipses).
    let magnitude = if central { wh.attrs.ratio } else { wh.attrs.fraction_diameter };
    Ok(Some(SolarEclipse {
        kind,
        central,
        noncentral: umbral && !central,
        max,
        gamma,
        magnitude,
        p1,
        p4,
        u1,
        u4,
        c1,
        c2,
        h1,
        h2,
    }))
}

/// Next (or previous) global solar eclipse whose greatest eclipse is more than
/// `SEARCH_MARGIN_DAYS` after (before) `jd_tt`.
pub fn solar_eclipse_search(engine: &Engine, jd_tt: f64, backward: bool) -> Result<SolarEclipse> {
    let mut k = first_lunation(jd_tt, 0.0, backward);
    let step = if backward { -1.0 } else { 1.0 };
    for _ in 0..MAX_LUNATIONS {
        let accept =
            |t: f64| if backward { t < jd_tt - SEARCH_MARGIN_DAYS } else { t > jd_tt + SEARCH_MARGIN_DAYS };
        if let Some(e) = solar_at_lunation(engine, k, accept)? {
            return Ok(e);
        }
        k += step;
    }
    Err(Error::InvalidArgument("no solar eclipse found".into()))
}

/// Local solar-eclipse attributes at one instant and place.
#[derive(Debug, Clone, Copy, Default)]
pub struct SolarAttrs {
    /// Fraction of the solar diameter covered (can exceed 1).
    pub fraction_diameter: f64,
    /// Moon / Sun apparent diameter ratio.
    pub ratio: f64,
    /// Fraction of the solar disc area covered.
    pub obscuration: f64,
    /// Diameter of the umbra (negative) or antumbra (positive) at the site's
    /// height along the axis, km (sign convention of the Swiss `attr[3]`).
    pub core_shadow_km: f64,
    /// Sun azimuth from north through east, degrees.
    pub azimuth_deg: f64,
    pub true_alt_deg: f64,
    pub apparent_alt_deg: f64,
    /// Sun–Moon separation, degrees.
    pub separation_deg: f64,
    /// NASA magnitude: `fraction_diameter` for a partial phase, `ratio` inside the umbra/antumbra.
    pub nasa_magnitude: f64,
    /// The site is inside the umbra (ratio ≥ 1) or antumbra (ratio < 1).
    pub in_umbra: bool,
    /// The shadow axis meets the Earth at this instant (the eclipse is central
    /// somewhere then): the classic CENTRAL / NONCENTRAL bit of local results.
    pub axis_on_earth: bool,
}

/// Overlap area of two discs (radii r1 = Sun, r2 = Moon, centre distance d) / area of disc 1.
pub(crate) fn obscuration(r1: f64, r2: f64, d: f64) -> f64 {
    if d >= r1 + r2 {
        return 0.0;
    }
    if d <= (r1 - r2).abs() {
        return if r2 >= r1 { 1.0 } else { (r2 / r1).powi(2) };
    }
    if r1 == 0.0 {
        return 0.0;
    }
    let a1 = ((d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1)).clamp(-1.0, 1.0).acos();
    let a2 = ((d * d + r2 * r2 - r1 * r1) / (2.0 * d * r2)).clamp(-1.0, 1.0).acos();
    let area = r1 * r1 * (a1 - a1.sin() * a1.cos()) + r2 * r2 * (a2 - a2.sin() * a2.cos());
    area / (PI * r1 * r1)
}

/// Topocentric Sun/Moon angles: (separation, Sun SD, Moon SD (mean k), Moon SD (umbral k)).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Topo {
    pub(crate) sep: f64,
    pub(crate) ss: f64,
    pub(crate) sm: f64,
    pub(crate) sm_u: f64,
    pub(crate) sun: Vec3,
    pub(crate) w: Vec3,
    pub(crate) dw: Vec3,
}

pub(crate) fn topo(g: &SunMoon, o: Vec3, vo: Vec3, sh: &Shape) -> Topo {
    let sv = sub(g.s, o);
    let mv = sub(g.m, o);
    let (us, dus) = unit_rate(sv, sub(g.vs, vo));
    let (um, dum) = unit_rate(mv, sub(g.vm, vo));
    let ds = norm(sv);
    let dm = norm(mv);
    Topo {
        sep: angle(sv, mv),
        ss: (sh.radius_km / ds).asin(),
        sm: (sh.k1 * RE / dm).asin(),
        sm_u: (sh.k2 * RE / dm).asin(),
        sun: sv,
        w: sub(um, us),
        dw: sub(dum, dus),
    }
}

pub(crate) fn solar_attrs(
    tp: &Topo,
    b: &Besselian,
    site: &Site,
    o: Vec3,
    fr: &FrameSet,
    jd_ut: f64,
) -> SolarAttrs {
    let (alt, az) = site.alt_az(fr, jd_ut, tp.sun);
    let central = tp.sep < (tp.sm_u - tp.ss).abs();
    // A point source (star) is either fully covered or not.
    let (frac, ratio) = if tp.ss > 0.0 {
        ((tp.ss + tp.sm - tp.sep) / (2.0 * tp.ss), tp.sm_u / tp.ss)
    } else {
        (if tp.sep < tp.sm { 1.0 } else { 0.0 }, 0.0)
    };
    let zeta = dot(o, b.zh) / RE;
    let l2 = b.l2 - zeta * b.tan_f2;
    let ta = alt / DEG;
    SolarAttrs {
        fraction_diameter: frac,
        ratio,
        obscuration: obscuration(tp.ss, if central { tp.sm_u } else { tp.sm }, tp.sep),
        core_shadow_km: 2.0 * l2 * RE,
        azimuth_deg: az / DEG,
        true_alt_deg: ta,
        apparent_alt_deg: apparent_altitude_deg(ta, site.height_m),
        separation_deg: tp.sep / DEG,
        nasa_magnitude: if central && tp.ss > 0.0 { ratio } else { frac.max(0.0) },
        in_umbra: central,
        axis_on_earth: b.edge() < 0.0,
    }
}

/// Type of solar eclipse seen at a place at one instant (`sol_eclipse_how`
/// flags): `None` when the Moon does not cover the Sun there or the Sun's
/// centre is not above the apparent horizon.
pub fn solar_seen_kind(a: &SolarAttrs) -> Option<SolarKind> {
    if a.fraction_diameter <= 0.0 || a.apparent_alt_deg <= 0.0 {
        None
    } else if a.in_umbra && (a.ratio >= 1.0 || a.ratio == 0.0) {
        Some(SolarKind::Total)
    } else if a.in_umbra {
        Some(SolarKind::Annular)
    } else {
        Some(SolarKind::Partial)
    }
}

/// Where the eclipse is central (or greatest, if the axis misses Earth) at an instant.
#[derive(Debug, Clone, Copy)]
pub struct SolarWhere {
    pub lon_deg: f64,
    pub lat_deg: f64,
    pub central: bool,
    /// Total (`Some(true)`) / annular (`Some(false)`) at that point, `None` if partial there.
    pub umbral: Option<bool>,
    pub attrs: SolarAttrs,
    /// (lon, lat) degrees of the northern and southern limits of the umbra
    /// (or antumbra) and of the penumbra at this instant; `None` where that
    /// limit does not touch the Earth. Filled by the public `*_where*` calls.
    pub limits: [Option<(f64, f64)>; 4],
}

pub(crate) fn solar_where_impl(sky: &Sky, fr: &FrameSet, jd_tt: f64, jd_ut: f64) -> Result<SolarWhere> {
    let kp = pole(fr);
    let g = sky.at(jd_tt)?;
    let b = besselian(&g, kp, &sky.shape);
    let (r, central) = match axis_surface(&b, kp) {
        Some((r, _)) => (r, true),
        None => (limb_point(&b, kp), false),
    };
    let (lon_deg, lat_deg) = geodetic_deg(fr, jd_ut, r);
    let site = Site::new(&Observer { lon_deg, lat_deg, height_m: 0.0 });
    let (o, vo) = site.gcrs(fr, jd_ut);
    let tp = topo(&g, o, vo, &sky.shape);
    let attrs = solar_attrs(&tp, &b, &site, o, fr, jd_ut);
    let umbral = if tp.sep < (tp.sm_u - tp.ss).abs() { Some(tp.sm_u > tp.ss) } else { None };
    Ok(SolarWhere { lon_deg, lat_deg, central, umbral, attrs, limits: [None; 4] })
}

/// Northern/southern limits of the umbra and penumbra at an instant: the
/// points of the shadow outline on the Earth (radius L = l − ζ·tan f at the
/// observer's height ζ along the axis) where the outline is tangent to the
/// path, i.e. where the shadow's motion *relative to the rotating Earth* is
/// parallel to the outline (Explanatory Supplement 3rd ed. §11.3; Meeus,
/// Elements of Solar Eclipses, path limits).
pub(crate) fn shadow_limits(
    sky: &Sky,
    fr: &FrameSet,
    jd_tt: f64,
    jd_ut: f64,
) -> Result<[Option<(f64, f64)>; 4]> {
    let kp = pole(fr);
    let b = besselian(&sky.at(jd_tt)?, kp, &sky.shape);
    let h = 1.0 / 1440.0;
    let (b0, b1) =
        (besselian(&sky.at(jd_tt - h)?, kp, &sky.shape), besselian(&sky.at(jd_tt + h)?, kp, &sky.shape));
    // Shadow-axis velocity on the fundamental plane, Earth radii per day.
    let (vx, vy) = ((b1.x - b0.x) / (2.0 * h), (b1.y - b0.y) / (2.0 * h));
    let zeta0 = axis_surface(&b, kp).map_or(0.0, |(_, z)| z);
    let limit = |l: f64, tan_f: f64, north: bool| -> Option<(f64, f64)> {
        let (mut zeta, mut rel) = (zeta0, (vx, vy));
        let mut pt = None;
        for _ in 0..6 {
            let rn = rel.0.hypot(rel.1);
            // Unit normal to the relative motion, on the requested side.
            let (mut nx, mut ny) = (-rel.1 / rn, rel.0 / rn);
            if (ny < 0.0) == north {
                (nx, ny) = (-nx, -ny);
            }
            let radius = (l - zeta * tan_f).abs() * RE;
            let p = add(b.p0, add(scale(b.xh, nx * radius), scale(b.yh, ny * radius)));
            let (q, z) = line_surface(p, b.zh, kp)?;
            zeta = z;
            // Observer velocity from the Earth's rotation, projected (Earth radii/day).
            let vq = scale(cross(kp, q), OMEGA_UT / RE);
            rel = (vx - dot(vq, b.xh), vy - dot(vq, b.yh));
            pt = Some(q);
        }
        pt.map(|q| geodetic_deg(fr, jd_ut, q))
    };
    Ok([
        limit(b.l2, b.tan_f2, true),
        limit(b.l2, b.tan_f2, false),
        limit(b.l1, b.tan_f1, true),
        limit(b.l1, b.tan_f1, false),
    ])
}

/// Time (JD TT) between first and last contact when the body is on the
/// meridian (hour angle 0, upper transit) at the point of greatest eclipse of
/// that instant ("eclipse at local apparent noon"); `None` if it never is.
pub(crate) fn local_noon_impl(engine: &Engine, sky: &Sky, ev: &SolarEclipse) -> Result<Option<f64>> {
    let fr = frames_at(engine, ev.max);
    let dt = ev.max - engine.tt_to_ut(ev.max);
    let ha = |t: f64| -> Result<f64> {
        let w = solar_where_impl(sky, &fr, t, t - dt)?;
        let v = mat_vec(&fr.bpn, sky.at(t)?.s);
        let ra = v[1].atan2(v[0]);
        Ok(diff_deg((fr.gast(t - dt) - ra) / DEG + w.lon_deg, 0.0))
    };
    let n = (((ev.p4 - ev.p1) / (20.0 / 1440.0)).ceil() as usize).max(1);
    let (mut ta, mut fa) = (ev.p1, ha(ev.p1)?);
    let mut lower = None;
    for i in 1..=n {
        let tb = ev.p1 + (ev.p4 - ev.p1) * i as f64 / n as f64;
        let fb = ha(tb)?;
        // The point of greatest eclipse can move faster than the Earth turns,
        // so the hour angle may cross the meridian in either direction.
        if (fa > 0.0) != (fb > 0.0) && fa.abs() < 90.0 && fb.abs() < 90.0 {
            return Ok(Some(brent(&ha, ta, tb, fa, fb)?));
        }
        // Lower transit (hour angle through 180°), used when there is no
        // upper one: polar points of greatest eclipse in the midnight Sun.
        if lower.is_none() && (fa > 0.0) != (fb > 0.0) && fa.abs() > 90.0 && fb.abs() > 90.0 {
            let g = |t: f64| ha(t).map(|h| norm_deg(h) - 180.0);
            lower = Some(brent(&g, ta, tb, norm_deg(fa) - 180.0, norm_deg(fb) - 180.0)?);
        }
        ta = tb;
        fa = fb;
    }
    Ok(lower)
}

/// `sol_eclipse_tret[1]`: the eclipse at local apparent noon (JD TT), if any.
pub fn solar_local_noon(engine: &Engine, ev: &SolarEclipse) -> Result<Option<f64>> {
    local_noon_impl(engine, &Sky::new(engine.pack(), ev.max), ev)
}

/// `sol_eclipse_where` equivalent at JD TT.
pub fn solar_where(engine: &Engine, jd_tt: f64) -> Result<SolarWhere> {
    solar_where_ut(engine, jd_tt, engine.tt_to_ut(jd_tt))
}

/// `solar_where` with an explicit UT1 for the Earth's rotation (isolates
/// Delta T model differences when comparing with other ephemerides).
pub fn solar_where_ut(engine: &Engine, jd_tt: f64, jd_ut: f64) -> Result<SolarWhere> {
    let sky = Sky::new(engine.pack(), jd_tt);
    where_with_limits(&sky, &frames_at(engine, jd_tt), jd_tt, jd_ut)
}

pub(crate) fn where_with_limits(sky: &Sky, fr: &FrameSet, jd_tt: f64, jd_ut: f64) -> Result<SolarWhere> {
    let mut w = solar_where_impl(sky, fr, jd_tt, jd_ut)?;
    w.limits = shadow_limits(sky, fr, jd_tt, jd_ut)?;
    Ok(w)
}

/// Local attributes of the Sun/Moon configuration at an instant and place
/// (`sol_eclipse_how` equivalent).
pub fn solar_how(engine: &Engine, jd_tt: f64, obs: &Observer) -> Result<SolarAttrs> {
    let sky = Sky::new(engine.pack(), jd_tt);
    let fr = frames_at(engine, jd_tt);
    let g = sky.at(jd_tt)?;
    let b = besselian(&g, pole(&fr), &SUN_SHAPE);
    let jd_ut = engine.tt_to_ut(jd_tt);
    let site = Site::new(obs);
    let (o, vo) = site.gcrs(&fr, jd_ut);
    Ok(solar_attrs(&topo(&g, o, vo, &SUN_SHAPE), &b, &site, o, &fr, jd_ut))
}

/// Local circumstances of a solar eclipse; times JD TT. Contact times are
/// geometric (also when the Sun is below the horizon); the `*_visible` flags
/// say whether the Sun is up then.
#[derive(Debug, Clone, Copy)]
pub struct SolarLocal {
    /// Local type: Partial, Annular or Total.
    pub kind: SolarKind,
    pub max: f64,
    pub c1: f64,
    pub c2: Option<f64>,
    pub c3: Option<f64>,
    pub c4: f64,
    /// Sunrise / sunset between C1 and C4.
    pub sunrise: Option<f64>,
    pub sunset: Option<f64>,
    pub visible: bool,
    pub max_visible: bool,
    pub c1_visible: bool,
    pub c2_visible: bool,
    pub c3_visible: bool,
    pub c4_visible: bool,
    pub attrs: SolarAttrs,
}

/// True altitude of the Sun's centre (deg) at which the Swiss oracle's
/// `sol_eclipse_when_loc` puts sunrise/sunset and its visibility flags,
/// measured black-box (see docs/ECLIPSES_REPORT.md): s☉ + 0.053° − 0.665°·P/P₀
/// with the standard-atmosphere pressure ratio at the site height
/// (−0.349° at sea level in June). It corresponds to the apparent *lower*
/// limb on the horizon, i.e. stricter than the usual upper-limb sunrise.
pub fn sun_horizon_true_alt_deg(sun_sd_deg: f64, height_m: f64) -> f64 {
    let p_ratio = (1.0 - 0.0065 * height_m.max(0.0) / 288.0).max(0.0).powf(5.255);
    sun_sd_deg + 0.053 - 0.665 * p_ratio
}

fn solar_local_for(engine: &Engine, glob: &SolarEclipse, obs: &Observer) -> Result<Option<SolarLocal>> {
    let sky = Sky::new(engine.pack(), glob.max);
    let fr = frames_at(engine, glob.max);
    let dt = glob.max - engine.tt_to_ut(glob.max);
    let site = Site::new(obs);
    let tp_at = |t: f64| -> Result<(Topo, SunMoon, Vec3)> {
        let g = sky.at(t)?;
        let (o, vo) = site.gcrs(&fr, t - dt);
        Ok((topo(&g, o, vo, &SUN_SHAPE), g, o))
    };
    let max = newton_min(
        |t| {
            let (tp, ..) = tp_at(t)?;
            Ok((tp.w, tp.dw))
        },
        glob.max,
    )?;
    if (max - glob.max).abs() > 0.4 {
        return Ok(None);
    }
    let (tp, g, o) = tp_at(max)?;
    if tp.sep >= tp.ss + tp.sm {
        return Ok(None);
    }
    let outer = |t: f64| -> Result<f64> {
        let (tp, ..) = tp_at(t)?;
        Ok(tp.sep - (tp.ss + tp.sm))
    };
    let inner = |t: f64| -> Result<f64> {
        let (tp, ..) = tp_at(t)?;
        Ok(tp.sep - (tp.sm_u - tp.ss).abs())
    };
    let span = 0.25;
    let c1 = root(outer, max, max - span)?.unwrap_or(max);
    let c4 = root(outer, max, max + span)?.unwrap_or(max);
    let (mut c2, mut c3) = (None, None);
    let mut kind = SolarKind::Partial;
    if tp.sep < (tp.sm_u - tp.ss).abs() {
        kind = if tp.sm_u > tp.ss { SolarKind::Total } else { SolarKind::Annular };
        c2 = root(inner, max, c1)?;
        c3 = root(inner, max, c4)?;
    }
    let alt = |t: f64| -> Result<f64> {
        let g = sky.at(t)?;
        let (o, _) = site.gcrs(&fr, t - dt);
        let sv = sub(g.s, o);
        let sd = (SUN_RADIUS_KM / norm(sv)).asin() / DEG;
        Ok(site.alt_az(&fr, t - dt, sv).0 / DEG - sun_horizon_true_alt_deg(sd, site.height_m))
    };
    // Visibility of an instant (contacts, maximum): apparent altitude of the
    // centre above the horizon, i.e. true altitude above ≈ −0.58° (measured
    // on the oracle; it differs from the sunrise/sunset rule above).
    let up = |t: f64| -> Result<bool> {
        let g = sky.at(t)?;
        let (o, _) = site.gcrs(&fr, t - dt);
        let ta = site.alt_az(&fr, t - dt, sub(g.s, o)).0 / DEG;
        Ok(apparent_altitude_deg(ta, site.height_m) > 0.0)
    };
    // Sunrise / sunset inside [C1, C4]: scan at 5-minute steps (the Sun's
    // altitude is monotonic over far longer near the horizon except at polar sites).
    let (mut sunrise, mut sunset) = (None, None);
    let n = (((c4 - c1) / (5.0 / 1440.0)).ceil() as usize).max(1);
    let mut ta = c1;
    let mut fa = alt(ta)?;
    for i in 1..=n {
        let tb = c1 + (c4 - c1) * i as f64 / n as f64;
        let fb = alt(tb)?;
        if fa.signum() != fb.signum() {
            let r = brent(&alt, ta, tb, fa, fb)?;
            if fb > 0.0 {
                sunrise.get_or_insert(r);
            } else {
                sunset.get_or_insert(r);
            }
        }
        ta = tb;
        fa = fb;
    }
    let c1_visible = up(c1)?;
    let c4_visible = up(c4)?;
    let max_visible = up(max)?;
    // Whole-eclipse visibility (and the search's skip rule), as measured on
    // the oracle: the maximum is visible (per-instant rule), or the Sun is
    // above the sunrise/sunset horizon at some time in [C1, C4]. A contact
    // seen only in the gap between the two horizons does not count.
    let visible = max_visible || alt(c1)? > 0.0 || alt(c4)? > 0.0 || sunrise.is_some() || sunset.is_some();
    let c2_visible = match c2 {
        Some(t) => up(t)?,
        None => false,
    };
    let c3_visible = match c3 {
        Some(t) => up(t)?,
        None => false,
    };
    // Like the oracle, report the greatest *visible* phase: a maximum below
    // the horizon moves to the sunrise / sunset inside the eclipse.
    let mut max = max;
    let (mut tp, mut g, mut o) = (tp, g, o);
    if !max_visible {
        let shifted = match (sunrise, sunset) {
            (Some(r), _) if max < r => Some(r),
            (_, Some(s)) if max > s => Some(s),
            _ => None,
        };
        if let Some(t) = shifted {
            max = t;
            (tp, g, o) = tp_at(t)?;
            // The local type then describes the visible part (oracle behaviour).
            if tp.sep >= (tp.sm_u - tp.ss).abs() {
                kind = SolarKind::Partial;
            }
        }
    }
    let b = besselian(&g, pole(&fr), &SUN_SHAPE);
    let attrs = solar_attrs(&tp, &b, &site, o, &fr, max - dt);
    Ok(Some(SolarLocal {
        kind,
        max,
        c1,
        c2,
        c3,
        c4,
        sunrise,
        sunset,
        visible,
        max_visible,
        c1_visible,
        c2_visible,
        c3_visible,
        c4_visible,
        attrs,
    }))
}

/// Next (previous) solar eclipse *visible* at `obs` (Sun up at some time
/// between C1 and C4) whose local maximum is more than `SEARCH_MARGIN_DAYS`
/// after (before) `jd_tt`. Returns the global eclipse too.
pub fn solar_local_search(
    engine: &Engine,
    jd_tt: f64,
    obs: &Observer,
    backward: bool,
) -> Result<(SolarLocal, SolarEclipse)> {
    let mut t = if backward { jd_tt + 0.5 } else { jd_tt - 0.5 };
    for _ in 0..200 {
        let glob = solar_eclipse_search(engine, t, backward)?;
        if let Some(loc) = solar_local_for(engine, &glob, obs)? {
            let ok = if backward {
                loc.max < jd_tt - SEARCH_MARGIN_DAYS
            } else {
                loc.max > jd_tt + SEARCH_MARGIN_DAYS
            };
            if ok && loc.visible {
                return Ok((loc, glob));
            }
        }
        t = glob.max;
    }
    Err(Error::InvalidArgument("no local solar eclipse found".into()))
}

/// Local circumstances of a given global eclipse (no visibility filter).
pub fn solar_local(engine: &Engine, glob: &SolarEclipse, obs: &Observer) -> Result<Option<SolarLocal>> {
    solar_local_for(engine, glob, obs)
}

// ---------------------------------------------------------------------------------------------
// Saros
// ---------------------------------------------------------------------------------------------

/// Saros series and member number of an eclipse (greatest eclipse at
/// `max_tt`), from NASA's series catalogue: members of one series are 223
/// lunations apart, numbered from the series' first eclipse. `None` outside
/// the catalogued series (e.g. a marginal eclipse NASA does not list).
pub fn saros(max_tt: f64, solar: bool) -> Option<(i32, i32)> {
    let off = if solar { 0.0 } else { 0.5 };
    let lunation = |jd: f64| ((jd - 2_451_550.097_66) / SYNODIC - off).round() as i64;
    let k = lunation(max_tt);
    let table: &[(i32, i32, f64)] = if solar { &saros::SOLAR } else { &saros::LUNAR };
    table.iter().find_map(|&(series, members, first)| {
        // The catalogue gives the civil date (0h): use mid-day.
        let d = k - lunation(first + 0.5);
        (d >= 0 && d % 223 == 0 && d / 223 < members as i64).then(|| (series, (d / 223 + 1) as i32))
    })
}

// ---------------------------------------------------------------------------------------------
// UT helpers
// ---------------------------------------------------------------------------------------------

/// UT → TT with the default Delta T model. Engine-level code should use
/// `Engine::ut_to_tt`, which honours the engine's Delta T model.
pub fn tt(jd_ut: f64) -> f64 {
    time::ut_to_tt(jd_ut)
}

/// TT → UT with the default Delta T model (see `tt`).
pub fn ut(jd_tt: f64) -> f64 {
    time::tt_to_ut(jd_tt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ellipse_distance() {
        assert!((ellipse_dist(2.0, 0.0, 0.99) - 1.0).abs() < 1e-12);
        assert!((ellipse_dist(0.0, 2.0, 0.99) - 1.01).abs() < 1e-12);
        assert!((ellipse_dist(0.0, 0.5, 0.99) + 0.49).abs() < 1e-12);
        let (x, y) = (1.2f64 * 0.6, 1.2f64 * 0.8);
        let d = ellipse_dist(x, y, 1.0);
        assert!((d - 0.2).abs() < 1e-12, "{d}");
    }

    #[test]
    fn obscuration_limits() {
        assert_eq!(obscuration(1.0, 1.1, 0.05), 1.0);
        assert!((obscuration(1.0, 0.9, 0.05) - 0.81).abs() < 1e-12);
        assert_eq!(obscuration(1.0, 1.0, 2.0), 0.0);
        // Half-overlap of equal discs at d = r: 2(π/3 − √3/4) / π.
        let v = obscuration(1.0, 1.0, 1.0);
        assert!((v - 2.0 * (PI / 3.0 - 3f64.sqrt() / 4.0) / PI).abs() < 1e-12);
    }

    #[test]
    fn mean_phase_meeus_example() {
        // Meeus example 49.a: new moon 1977 Feb, k = −283, JDE 2443192.94102 (mean 2443192.94101).
        let (jde, _) = mean_phase(-283.0);
        assert!((jde - 2_443_192.941_01).abs() < 1e-4, "{jde}");
    }
}
