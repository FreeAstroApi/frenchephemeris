// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Position pipeline: barycentric JPL states -> apparent/true, geo/topocentric,
//! equatorial/ecliptic coordinates.
//!
//! Order of operations for an apparent geocentric position (IERS Conventions /
//! SOFA "astrometric -> apparent" chain):
//! 1. observer barycentric state at t (Earth, plus topocentric offset);
//! 2. iterative light time: target at t − τ;
//! 3. gravitational light deflection by the Sun;
//! 4. annual (relativistic) aberration from the observer velocity;
//! 5. GCRS -> equator/ecliptic of date (or J2000) rotation.
//!
//! Every step is controlled by `Corrections`, so Swiss-style flags such as
//! TRUEPOS / NOABERR / NOGDEFL map one-to-one.

use crate::frames::{self, FrameModel, FrameSet};
use crate::math::*;
use crate::pack::Pack;
use crate::time;
use crate::Result;

mod solar_mass;

/// Schwarzschild radius of the Sun in AU (2·GM☉/c²), as used by SOFA.
const SRS_AU: f64 = 1.974_125_743_36e-8;
/// Earth rotation rate, rad/s (IERS nominal).
const OMEGA_EARTH: f64 = 7.292_115_0e-5;
/// WGS84 ellipsoid.
const WGS84_A_KM: f64 = 6378.137;
const WGS84_F: f64 = 1.0 / 298.257_223_563;

pub const SSB: i32 = 0;
pub const EMB: i32 = 3;
pub const SUN: i32 = 10;
pub const MOON: i32 = 301;
pub const EARTH: i32 = 399;

/// Geodetic observer on the WGS84 ellipsoid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observer {
    pub lon_deg: f64,
    pub lat_deg: f64,
    pub height_m: f64,
}

impl Observer {
    /// Earth-fixed (ITRS, polar motion neglected) position in km.
    pub fn itrs_km(&self) -> Vec3 {
        let (sl, cl) = (self.lat_deg * DEG).sin_cos();
        let (so, co) = (self.lon_deg * DEG).sin_cos();
        let e2 = WGS84_F * (2.0 - WGS84_F);
        let n = WGS84_A_KM / (1.0 - e2 * sl * sl).sqrt();
        let h = self.height_m / 1000.0;
        [(n + h) * cl * co, (n + h) * cl * so, (n * (1.0 - e2) + h) * sl]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Origin {
    Geocentric,
    Topocentric(Observer),
    Heliocentric,
    Barycentric,
    /// Centre of any pack body (NAIF id): planet-centric places (the classic
    /// `calc_pctr`). Light time, solar deflection and aberration (from the
    /// centre's barycentric velocity) apply as for a geocentric observer.
    Body(i32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Corrections {
    pub light_time: bool,
    pub deflection: bool,
    pub aberration: bool,
}

impl Corrections {
    /// Apparent place (Swiss default).
    pub const APPARENT: Corrections = Corrections { light_time: true, deflection: true, aberration: true };
    /// Geometric place at the instant (Swiss TRUEPOS).
    pub const GEOMETRIC: Corrections =
        Corrections { light_time: false, deflection: false, aberration: false };
    /// Astrometric place (light time + deflection, no aberration).
    pub const ASTROMETRIC: Corrections =
        Corrections { light_time: true, deflection: true, aberration: false };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Frame {
    /// Ecliptic of date; `true_equinox` adds nutation (Swiss default), else mean (NONUT).
    EclipticOfDate { true_equinox: bool },
    /// Equator of date (RA/Dec).
    EquatorialOfDate { true_equinox: bool },
    /// Mean ecliptic & equinox J2000.
    EclipticJ2000,
    /// ICRF / GCRS axes.
    Icrf,
    /// Mean equator & equinox J2000 (ICRS rotated by the IAU 2006 frame bias).
    EquatorialJ2000,
}

/// Everything about one instant that does not depend on the target body.
/// Build once, reuse for every body of a chart.
#[derive(Debug, Clone, Copy)]
pub struct Epoch {
    pub jd_ut: f64,
    pub jd_tt: f64,
    pub jd_tdb: f64,
    pub frames: FrameSet,
}

impl Epoch {
    pub fn from_ut(jd_ut: f64) -> Epoch {
        Epoch::from_ut_with(jd_ut, FrameModel::DEFAULT)
    }
    pub fn from_tt(jd_tt: f64) -> Epoch {
        Epoch::with_scales(time::tt_to_ut(jd_tt), jd_tt, FrameModel::DEFAULT)
    }
    pub fn from_ut_with(jd_ut: f64, model: FrameModel) -> Epoch {
        Epoch::with_scales(jd_ut, time::ut_to_tt(jd_ut), model)
    }
    pub fn from_tt_with(jd_tt: f64, model: FrameModel) -> Epoch {
        Epoch::with_scales(time::tt_to_ut(jd_tt), jd_tt, model)
    }
    pub fn with_scales(jd_ut: f64, jd_tt: f64, model: FrameModel) -> Epoch {
        Epoch { jd_ut, jd_tt, jd_tdb: time::tt_to_tdb(jd_tt), frames: FrameSet::with_model(jd_tt, model) }
    }
    /// Greenwich apparent sidereal time, rad, consistent with this epoch's frames.
    pub fn gast(&self) -> f64 {
        self.frames.gast(self.jd_ut)
    }

    /// This epoch moved by `h` days, with frames extrapolated linearly from
    /// `rates`. For |h| ≲ 0.01 d the matrix error is O(h²·M̈) < 1e-12 rad, so
    /// numerical speeds need one full frame evaluation instead of three.
    pub fn shifted(&self, h: f64, rates: &FrameRates) -> Epoch {
        let f = &self.frames;
        let mv = |m: &Mat3, r: &Mat3| {
            let mut o = *m;
            for i in 0..3 {
                for j in 0..3 {
                    o[i][j] += h * r[i][j];
                }
            }
            o
        };
        let frames = FrameSet {
            jd_tt: f.jd_tt + h,
            dpsi: f.dpsi + h * rates.dpsi,
            bp: mv(&f.bp, &rates.bp),
            bpn: mv(&f.bpn, &rates.bpn),
            ecl_mean: mv(&f.ecl_mean, &rates.ecl_mean),
            ecl_true: mv(&f.ecl_true, &rates.ecl_true),
            ..*f
        };
        Epoch { jd_ut: self.jd_ut + h, jd_tt: self.jd_tt + h, jd_tdb: self.jd_tdb + h, frames }
    }
}

/// Time derivatives (per day) of the frame matrices of an epoch.
#[derive(Debug, Clone, Copy)]
pub struct FrameRates {
    pub bp: Mat3,
    pub bpn: Mat3,
    pub ecl_mean: Mat3,
    pub ecl_true: Mat3,
    pub dpsi: f64,
}

impl FrameRates {
    /// Central difference over ±0.05 d of the 77-term nutation variant of
    /// `model`: rates good to < 1e-9 rad/day (IAU 2000A differs only in
    /// small, slow terms), at ~2 cheap frame evaluations.
    pub fn new(jd_tt: f64, model: FrameModel) -> FrameRates {
        const H: f64 = 0.05;
        let fast = FrameModel { nutation: frames::Nutation::SwissCompat, ..model };
        let mats = |jd: f64| {
            let (dpsi, deps) = frames::nutation_with(jd, fast.nutation);
            let bp = frames::bias_precession_matrix(jd, fast.precession);
            let eps = frames::obliquity_mean_with(jd, fast.precession, fast.obliquity);
            let bpn = mat_mul(&frames::nutation_matrix(eps, dpsi, deps), &bp);
            (bp, bpn, mat_mul(&rot_x(eps), &bp), mat_mul(&rot_x(eps + deps), &bpn), dpsi)
        };
        let a = mats(jd_tt - H);
        let b = mats(jd_tt + H);
        let d = |x: &Mat3, y: &Mat3| {
            let mut o = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    o[i][j] = (y[i][j] - x[i][j]) / (2.0 * H);
                }
            }
            o
        };
        FrameRates {
            bp: d(&a.0, &b.0),
            bpn: d(&a.1, &b.1),
            ecl_mean: d(&a.2, &b.2),
            ecl_true: d(&a.3, &b.3),
            dpsi: (b.4 - a.4) / (2.0 * H),
        }
    }
}

/// Barycentric state of a pack body (target wrt SSB).
#[inline]
fn bary(pack: &Pack, target: i32, jd_tdb: f64) -> Result<(Vec3, Vec3)> {
    pack.state(target, SSB, jd_tdb)
}

/// Observer barycentric state (km, km/day) in GCRS axes.
pub fn observer_state(pack: &Pack, ep: &Epoch, origin: &Origin) -> Result<(Vec3, Vec3)> {
    match origin {
        Origin::Barycentric => Ok(([0.0; 3], [0.0; 3])),
        Origin::Heliocentric => bary(pack, SUN, ep.jd_tdb),
        Origin::Geocentric => bary(pack, EARTH, ep.jd_tdb),
        Origin::Body(naif) => bary(pack, *naif, ep.jd_tdb),
        Origin::Topocentric(obs) => {
            let (e, ev) = bary(pack, EARTH, ep.jd_tdb)?;
            let (r, v) = topocentric_offset_gcrs(ep, obs);
            Ok((add(e, r), add(ev, v)))
        }
    }
}

/// Observer-side state shared by every target at one epoch: barycentric
/// observer position/velocity and the barycentric Sun position.
#[derive(Debug, Clone, Copy)]
pub struct ObsState {
    pub pos: Vec3,
    pub vel: Vec3,
    pub sun: Vec3,
    /// Geocentric or topocentric: deflection and aberration apply.
    pub earthbound: bool,
    pub barycentric: bool,
}

impl ObsState {
    pub fn new(pack: &Pack, ep: &Epoch, origin: &Origin) -> Result<ObsState> {
        let (pos, vel) = observer_state(pack, ep, origin)?;
        let earthbound = match origin {
            Origin::Geocentric | Origin::Topocentric(_) => true,
            Origin::Body(naif) => *naif != SUN,
            Origin::Heliocentric | Origin::Barycentric => false,
        };
        let sun = if earthbound { pack.position(SUN, SSB, ep.jd_tdb)? } else { [0.0; 3] };
        Ok(ObsState { pos, vel, sun, earthbound, barycentric: matches!(origin, Origin::Barycentric) })
    }
}

/// Observer position/velocity relative to the geocenter, GCRS, km and km/day.
pub fn topocentric_offset_gcrs(ep: &Epoch, obs: &Observer) -> (Vec3, Vec3) {
    let r_itrs = obs.itrs_km();
    let gast = ep.gast();
    // ITRS -> true equator of date: undo the Earth rotation by GAST.
    let r_tod = mat_vec(&rot_z(-gast), r_itrs);
    let w = OMEGA_EARTH * 86400.0;
    let v_tod = [-w * r_tod[1], w * r_tod[0], 0.0];
    (mat_t_vec(&ep.frames.bpn, r_tod), mat_t_vec(&ep.frames.bpn, v_tod))
}

/// GCRS vector of `target` as seen from `origin` with the requested corrections.
/// Returns the direction-corrected vector scaled to the light-time distance.
pub fn gcrs_vector(pack: &Pack, ep: &Epoch, target: i32, origin: &Origin, corr: Corrections) -> Result<Vec3> {
    gcrs_vector_obs(pack, ep, &ObsState::new(pack, ep, origin)?, target, corr)
}

/// `gcrs_vector` with a precomputed observer state (reuse it across bodies).
pub fn gcrs_vector_obs(
    pack: &Pack,
    ep: &Epoch,
    obs: &ObsState,
    target: i32,
    corr: Corrections,
) -> Result<Vec3> {
    Ok(apparent(pack, ep, obs, target, corr, LightTime::Iterate)?.0)
}

/// How to obtain the light time τ for one evaluation.
#[derive(Debug, Clone, Copy)]
pub enum LightTime {
    /// Iterate from the geometric distance (two iterations, < 1e-7 s for Pluto).
    Iterate,
    /// Use a known τ: e.g. τ(t) + τ̇·h for a speed sample at t + h, whose
    /// error O(h²·τ̈) is ~1e-15 d for h = 1e-3 d.
    Known(f64),
}

/// Core of the pipeline. Returns (GCRS vector, τ days, dτ/dt).
pub fn apparent(
    pack: &Pack,
    ep: &Epoch,
    obs: &ObsState,
    target: i32,
    corr: Corrections,
    lt: LightTime,
) -> Result<(Vec3, f64, f64)> {
    let pos = |t: f64| pack.position(target, SSB, t);
    let state = |t: f64| pack.state(target, SSB, t);
    apparent_with(ep, obs, target == SUN, corr, lt, pos, state)
}

/// [`apparent`] for a target given by functions of JD TDB returning its
/// barycentric GCRS position / state (km, km/day) — e.g. a body computed
/// from orbital elements (`fictitious`). `is_sun` disables self-deflection.
pub fn apparent_with<P, S>(
    ep: &Epoch,
    obs: &ObsState,
    is_sun: bool,
    corr: Corrections,
    lt: LightTime,
    position: P,
    state: S,
) -> Result<(Vec3, f64, f64)>
where
    P: Fn(f64) -> Result<Vec3>,
    S: Fn(f64) -> Result<(Vec3, Vec3)>,
{
    let o = obs.pos;

    // 1-2. light time. Each iteration contracts the error by ~v/c (1e-4), so
    // two iterations after the geometric guess leave < 1e-7 s even for Pluto.
    let (mut tau, mut taudot) = (0.0, 0.0);
    let p = if !corr.light_time {
        position(ep.jd_tdb)?
    } else {
        match lt {
            LightTime::Known(t) => {
                tau = t;
                position(ep.jd_tdb - t)?
            }
            LightTime::Iterate => {
                let p0 = position(ep.jd_tdb)?;
                tau = norm(sub(p0, o)) / C_KM_DAY;
                let p1 = position(ep.jd_tdb - tau)?;
                tau = norm(sub(p1, o)) / C_KM_DAY;
                let (p2, v2) = state(ep.jd_tdb - tau)?;
                let rel = sub(p2, o);
                let r = norm(rel);
                // d|P(t−τ) − O(t)|/dt / c, to first order in v/c.
                taudot = dot(rel, sub(v2, obs.vel)) / (r * C_KM_DAY);
                p2
            }
        }
    };
    let rel = sub(p, o);
    let dist = norm(rel);
    let mut u = scale(rel, 1.0 / dist);
    let s = obs.sun;

    // 3. light deflection by the Sun (not for the Sun itself).
    if corr.deflection && obs.earthbound && !is_sun {
        let q_vec = sub(p, s);
        let e_vec = sub(o, s);
        let em_au = norm(e_vec) / AU_KM;
        let q = scale(q_vec, 1.0 / norm(q_vec));
        let e = scale(e_vec, 1.0 / norm(e_vec));
        let qpe = add(q, e);
        let qdqpe = dot(q, qpe);
        // A body seen behind the solar disc: deflect its light by the mass enclosed in the
        // line-of-sight cylinder of radius b (impact parameter), the Sun treated as
        // transparent. Continuous at the limb, zero at the centre, instead of diverging as
        // q·(q+e) -> 0 at exact superior conjunction (solar Model S, see ephem/solar_mass.rs).
        let em_km = norm(e_vec);
        let cos_sep = -dot(u, e); // angle between the body and the Sun's centre, seen by the observer
        let b_km = em_km * (1.0 - cos_sep * cos_sep).max(0.0).sqrt();
        let behind = cos_sep > 0.0 && dist > em_km && b_km < solar_mass::SUN_RADIUS_KM;
        let mass_fraction =
            if behind { solar_mass::cylinder_mass_fraction(b_km / solar_mass::SUN_RADIUS_KM) } else { 1.0 };
        let w = mass_fraction * SRS_AU / em_au / qdqpe.max(1e-15);
        let peq = cross(u, cross(e, q));
        u = add(u, scale(peq, w));
    }

    // 4. aberration (SOFA iauAb, relativistic).
    if corr.aberration && !obs.barycentric {
        let v = scale(obs.vel, 1.0 / C_KM_DAY);
        let bm1 = (1.0 - dot(v, v)).sqrt();
        let s_au = if obs.earthbound { norm(sub(o, s)) / AU_KM } else { f64::INFINITY };
        let pdv = dot(u, v);
        let w1 = 1.0 + pdv / (1.0 + bm1);
        let w2 = SRS_AU / s_au;
        let mut p1 = [0.0; 3];
        for i in 0..3 {
            p1[i] = u[i] * bm1 + w1 * v[i] + w2 * (v[i] - pdv * u[i]);
        }
        let n = norm(p1);
        u = scale(p1, 1.0 / n);
    }

    Ok((scale(u, dist), tau, taudot))
}

/// IAU 2006 frame bias matrix (ICRS -> mean J2000), computed once.
pub fn frame_bias() -> &'static Mat3 {
    static B: std::sync::OnceLock<Mat3> = std::sync::OnceLock::new();
    B.get_or_init(frames::frame_bias)
}

/// Rotate a GCRS vector into the requested output frame.
pub fn to_frame(ep: &Epoch, frame: Frame, v_gcrs: Vec3) -> Vec3 {
    let f = &ep.frames;
    match frame {
        Frame::Icrf => v_gcrs,
        Frame::EquatorialJ2000 => mat_vec(frame_bias(), v_gcrs),
        Frame::EclipticJ2000 => mat_vec(&frames::gcrs_to_ecliptic_j2000(), v_gcrs),
        Frame::EquatorialOfDate { true_equinox } => {
            mat_vec(if true_equinox { &f.bpn } else { &f.bp }, v_gcrs)
        }
        Frame::EclipticOfDate { true_equinox } => {
            mat_vec(if true_equinox { &f.ecl_true } else { &f.ecl_mean }, v_gcrs)
        }
    }
}

/// Spherical place: (lon/RA rad, lat/Dec rad, distance km).
pub fn place(
    pack: &Pack,
    ep: &Epoch,
    target: i32,
    origin: &Origin,
    corr: Corrections,
    frame: Frame,
) -> Result<(f64, f64, f64)> {
    let v = gcrs_vector(pack, ep, target, origin, corr)?;
    Ok(to_spherical(to_frame(ep, frame, v)))
}

/// `place` with a precomputed observer state and light-time mode; also
/// returns (τ, dτ/dt) for reuse by speed samples.
pub fn place_obs(
    pack: &Pack,
    ep: &Epoch,
    obs: &ObsState,
    target: i32,
    corr: Corrections,
    frame: Frame,
    lt: LightTime,
) -> Result<((f64, f64, f64), f64, f64)> {
    let (v, tau, taudot) = apparent(pack, ep, obs, target, corr, lt)?;
    Ok((to_spherical(to_frame(ep, frame, v)), tau, taudot))
}
