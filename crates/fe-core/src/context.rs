// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Immutable calculation context: the public, degree-based API of the kernel.
//!
//! `Engine` owns the coefficient pack and is shared (`&Engine` / `Arc<Engine>`)
//! across threads. `CalcOptions` carries every setting that Swiss keeps as
//! process-global state (sidereal mode, topocentric observer, flags), so two
//! concurrent requests can never observe each other's configuration.

use crate::ayanamsa::{self, Ayanamsa, Kind, Projection, SidBits};
use crate::ephem::{self, Corrections, Epoch, Frame, FrameRates, LightTime, ObsState, Origin};
use crate::frames::FrameModel;
use crate::time::{self, DeltaTModel, DeltaTSpline};

type MemoKey = (u64, bool, FrameModel, DeltaTSpline, u64, Option<u64>);
use crate::houses::{self, Houses};
use crate::math::*;
use crate::pack::Pack;
use crate::points;
use crate::{Error, Result};

/// Body / point ids. Numbering follows the long-standing astrology convention
/// (Sun = 0 ... Pluto = 9, nodes 10/11, apogees 12/13, Chiron 15 ... Vesta 20,
/// numbered asteroids 10000 + MPC number) so adapters map 1:1.
///
/// Minor bodies need the minor-body pack (`Pack::open_many`); each has an
/// explicit JPL-covered range (docs/MINOR_BODIES_REPORT.md) outside which
/// positions are `Error::OutOfRange`, and without the pack `UnknownBody`.
/// Any other numbered asteroid is `Asteroid(mpc_number)` (id 10000 + number,
/// NAIF 2000000 + number), served by the asteroid pack or an on-demand pack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum Body {
    Sun,
    Moon,
    Mercury,
    Venus,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    MeanNode,
    TrueNode,
    MeanApogee,
    OscuApogee,
    Earth,
    Chiron,
    Pholus,
    Ceres,
    Pallas,
    Juno,
    Vesta,
    /// Natural ("interpolated") lunar apogee, natural Lilith (Swiss 21).
    IntpApogee,
    /// Natural ("interpolated") lunar perigee, Priapus (Swiss 22).
    IntpPerigee,
    // Hypothetical bodies (Swiss SE_FICT_OFFSET + n; crate::fictitious).
    Cupido = 40,
    Hades = 41,
    Zeus = 42,
    Kronos = 43,
    Apollon = 44,
    Admetos = 45,
    Vulkanus = 46,
    Poseidon = 47,
    Harrington = 50,
    NeptuneLeverrier = 51,
    NeptuneAdams = 52,
    PlutoLowell = 53,
    WhiteMoon = 56,
    Waldemath = 58,
    /// 136199 Eris (10000 + MPC number).
    Eris = 146199,
    /// 433 Eros (10000 + MPC number).
    Eros = 10433,
    /// Any other numbered asteroid by MPC number (id `AST_OFFSET + n`, 1 <= n < 1e6).
    Asteroid(u32),
}

/// Offset of numbered-asteroid ids (`10000 + MPC number`).
pub const AST_OFFSET: i32 = 10000;

impl Body {
    pub fn from_id(id: i32) -> Result<Body> {
        use Body::*;
        Ok(match id {
            0 => Sun,
            1 => Moon,
            2 => Mercury,
            3 => Venus,
            4 => Mars,
            5 => Jupiter,
            6 => Saturn,
            7 => Uranus,
            8 => Neptune,
            9 => Pluto,
            10 => MeanNode,
            11 => TrueNode,
            12 => MeanApogee,
            13 => OscuApogee,
            14 => Earth,
            15 => Chiron,
            16 => Pholus,
            17 => Ceres,
            18 => Pallas,
            19 => Juno,
            20 => Vesta,
            21 => IntpApogee,
            22 => IntpPerigee,
            40 => Cupido,
            41 => Hades,
            42 => Zeus,
            43 => Kronos,
            44 => Apollon,
            45 => Admetos,
            46 => Vulkanus,
            47 => Poseidon,
            50 => Harrington,
            51 => NeptuneLeverrier,
            52 => NeptuneAdams,
            53 => PlutoLowell,
            56 => WhiteMoon,
            58 => Waldemath,
            146199 => Eris,
            10433 => Eros,
            // the same bodies by MPC number
            10001 => Ceres,
            10002 => Pallas,
            10003 => Juno,
            10004 => Vesta,
            12060 => Chiron,
            15145 => Pholus,
            id if id > AST_OFFSET && id < AST_OFFSET + 1_000_000 => Asteroid((id - AST_OFFSET) as u32),
            _ => return Err(Error::UnknownBody(id)),
        })
    }

    pub fn id(self) -> i32 {
        match self {
            Body::Asteroid(n) => AST_OFFSET + n as i32,
            // SAFETY: `#[repr(i32)]` enums start with their i32 discriminant
            // (RFC 2195), and every other variant carries its id as discriminant.
            other => unsafe { *(&other as *const Body as *const i32) },
        }
    }

    /// NAIF id of the pack body, for physical bodies.
    pub fn naif(self) -> Option<i32> {
        use Body::*;
        Some(match self {
            Sun => ephem::SUN,
            Moon => ephem::MOON,
            Mercury => 1,
            Venus => 2,
            Mars => 4,
            Jupiter => 5,
            Saturn => 6,
            Uranus => 7,
            Neptune => 8,
            Pluto => 9,
            Earth => ephem::EARTH,
            // NAIF small-body ids: 2000000 + MPC number
            Chiron => 2002060,
            Pholus => 2005145,
            Ceres => 2000001,
            Pallas => 2000002,
            Juno => 2000003,
            Vesta => 2000004,
            Eris => 2136199,
            Eros => 2000433,
            Asteroid(n) => 2_000_000 + n as i32,
            _ => return None,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalcOptions {
    pub origin: Origin,
    pub corrections: Corrections,
    /// Equatorial (RA/Dec) instead of ecliptic output.
    pub equatorial: bool,
    /// Mean equinox (no nutation) instead of true equinox.
    pub mean_equinox: bool,
    /// Output in the J2000 frame instead of "of date".
    pub j2000: bool,
    /// Sidereal zodiac (ecliptic output only, except projection modes, whose
    /// RA/Dec refer to the mean equator of t0).
    pub sidereal: Option<Ayanamsa>,
    /// Swiss `SE_SIDBIT_*` option bits for `sidereal` (per call).
    pub sid_bits: SidBits,
    pub speed: bool,
    /// Cartesian x, y, z (AU) and their daily rates instead of polar output (`SEFLG_XYZ`).
    pub xyz: bool,
    /// Angles and angular rates in radians instead of degrees (`SEFLG_RADIANS`; no effect
    /// with `xyz`).
    pub radians: bool,
    /// No IAU 2006 frame bias (`SEFLG_ICRS`): the ICRS vector is precessed as if it were
    /// mean J2000, and J2000 equatorial output is ICRS.
    pub icrs: bool,
}

impl Default for CalcOptions {
    fn default() -> Self {
        CalcOptions {
            origin: Origin::Geocentric,
            corrections: Corrections::APPARENT,
            equatorial: false,
            mean_equinox: false,
            j2000: false,
            sidereal: None,
            sid_bits: SidBits(0),
            speed: true,
            xyz: false,
            radians: false,
            icrs: false,
        }
    }
}

impl CalcOptions {
    fn frame(&self) -> Frame {
        match (self.j2000, self.equatorial) {
            (true, false) => Frame::EclipticJ2000,
            (true, true) => Frame::EquatorialJ2000,
            (false, false) => Frame::EclipticOfDate { true_equinox: !self.mean_equinox },
            (false, true) => Frame::EquatorialOfDate { true_equinox: !self.mean_equinox },
        }
    }
}

/// Swiss-compatible 6-vector: lon/RA, lat/Dec (deg), distance (AU), and their daily rates.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Position {
    pub lon: f64,
    pub lat: f64,
    pub dist: f64,
    pub lon_speed: f64,
    pub lat_speed: f64,
    pub dist_speed: f64,
}

impl Position {
    pub fn to_array(self) -> [f64; 6] {
        [self.lon, self.lat, self.dist, self.lon_speed, self.lat_speed, self.dist_speed]
    }

    /// Apply the `xyz` / `radians` output options to a polar-degree position: cartesian
    /// coordinates with velocities by the chain rule, or angles in radians. Fields keep their
    /// names (`lon` = x, ...).
    pub fn formatted(self, opt: &CalcOptions) -> Position {
        if opt.xyz {
            let (l, b, r) = (self.lon * DEG, self.lat * DEG, self.dist);
            let (dl, db, dr) = (self.lon_speed * DEG, self.lat_speed * DEG, self.dist_speed);
            let (sl, cl) = l.sin_cos();
            let (sb, cb) = b.sin_cos();
            Position {
                lon: r * cb * cl,
                lat: r * cb * sl,
                dist: r * sb,
                lon_speed: dr * cb * cl - r * sb * db * cl - r * cb * sl * dl,
                lat_speed: dr * cb * sl - r * sb * db * sl + r * cb * cl * dl,
                dist_speed: dr * sb + r * cb * db,
            }
        } else if opt.radians {
            Position {
                lon: self.lon * DEG,
                lat: self.lat * DEG,
                lon_speed: self.lon_speed * DEG,
                lat_speed: self.lat_speed * DEG,
                ..self
            }
        } else {
            self
        }
    }
}

/// Half-step for numerical speeds, days. Central differences at this step are
/// accurate to ~1e-7 deg/day for the Moon and far better for planets.
const SPEED_STEP: f64 = 1e-3;

/// GM(Earth + Moon), km^3/day^2 (DE440/DE442 value 403503.235648 km^3/s^2).
const GM_EMB_KM3_DAY2: f64 = 403_503.235_648 * 86400.0 * 86400.0;

/// Frame model of `Engine::new`: Vondrák 2011 precession and the 77
/// luni-solar IAU 2000B terms evaluated with the full IERS 2003 fundamental
/// arguments (no linearisation, so no drift over 0–2400 — ERFA's `nut00b`
/// drifts to ~0.2″ by year 0). Within ~1–2 mas of full IAU 2000A, thousands of
/// times below the API's 0.001° rounding, at a fifth of the cost; use
/// `Engine::with_model(pack, FrameModel::DEFAULT)` for full 2000A.
/// Mean obliquity (ecliptic of date, houses, ECL_NUT): the published Vondrák et al. (2011)
/// ε_A series, identical to Swiss 2.10 over −13000 … +17000; the precession matrix itself
/// stays the Vondrák pole construction.
pub const API_FRAME_MODEL: FrameModel = FrameModel {
    precession: crate::frames::Precession::Vondrak2011,
    nutation: crate::frames::Nutation::SwissCompat,
    obliquity: crate::frames::Obliquity::Vondrak2011,
    sidereal: crate::frames::SiderealTime::Iau2006,
};

/// Model choices fixed for the lifetime of an `Engine`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EngineConfig {
    pub frames: FrameModel,
    pub delta_t: DeltaTModel,
}

impl EngineConfig {
    /// Physically best Delta T, fast frames (`API_FRAME_MODEL`).
    pub const DEFAULT: EngineConfig = EngineConfig { frames: API_FRAME_MODEL, delta_t: DeltaTModel::DEFAULT };
    /// Delta T reproducing the Swiss oracle before 1950 (to ~0.006 s), for
    /// routes that must keep historical charts unchanged.
    pub const SWISS_DELTA_T: EngineConfig =
        EngineConfig { frames: API_FRAME_MODEL, delta_t: DeltaTModel::SWISS_COMPAT };
    /// `SWISS_DELTA_T` plus the long-term sidereal time documented for Swiss Ephemeris
    /// (`SiderealTime::SwissLongTerm`) for ARMC/houses, topocentric places and events outside
    /// 1850–2050.
    pub const SWISS_TIME: EngineConfig = EngineConfig {
        frames: FrameModel { sidereal: crate::frames::SiderealTime::SwissLongTerm, ..API_FRAME_MODEL },
        delta_t: DeltaTModel::SWISS_COMPAT,
    };
}

impl Default for EngineConfig {
    fn default() -> Self {
        EngineConfig::DEFAULT
    }
}

pub struct Engine {
    pack: std::sync::Arc<Pack>,
    model: FrameModel,
    delta_t: DeltaTModel,
    /// Optional fixed-star compatibility alias table (`stars::StarCompat`).
    pub(crate) star_compat: Option<std::sync::Arc<crate::stars::StarCompat>>,
}

/// One instant prepared for any number of bodies: time scales, frames, the
/// geocentric observer state, and (for speeds) the t ± h sample epochs.
#[derive(Debug, Clone, Copy)]
pub struct Prepared {
    pub epoch: Epoch,
    /// (t − h, t + h) with linearly extrapolated frames, when prepared for speeds.
    pub shifted: Option<[Epoch; 2]>,
    /// Geocentric observer state at [t, t − h, t + h]; `None` if not computed
    /// (or out of range, in which case the per-body call reports the error).
    geo: [Option<ObsState>; 3],
}

impl Prepared {
    fn new(epoch: Epoch) -> Prepared {
        Prepared { epoch, shifted: None, geo: [None; 3] }
    }
}

thread_local! {
    /// Memo of the last prepared instant per thread: (jd bits, is_ut, model).
    /// Pure memoization of a deterministic function, so it carries no
    /// configuration between calls; it only makes body-by-body calls at the
    /// same time as cheap as a batch.
    static LAST: std::cell::Cell<Option<(MemoKey, Prepared)>> = const { std::cell::Cell::new(None) };
}

impl Engine {
    /// Engine with `API_FRAME_MODEL`.
    pub fn new(pack: Pack) -> Engine {
        Engine::with_model(pack, API_FRAME_MODEL)
    }

    pub fn with_model(pack: Pack, model: FrameModel) -> Engine {
        Engine::with_config(pack, EngineConfig { frames: model, ..EngineConfig::DEFAULT })
    }

    pub fn with_config(pack: Pack, cfg: EngineConfig) -> Engine {
        Engine { pack: std::sync::Arc::new(pack), model: cfg.frames, delta_t: cfg.delta_t, star_compat: None }
    }

    /// The same engine (shared pack, same frames) with another Delta T model, e.g. a
    /// user-defined Delta T or tidal acceleration (the classic `set_delta_t_userdef` /
    /// `set_tid_acc`, as a value instead of global state). Cheap: the pack is shared.
    pub fn with_delta_t(&self, delta_t: DeltaTModel) -> Engine {
        Engine { pack: self.pack.clone(), model: self.model, delta_t, star_compat: self.star_compat.clone() }
    }

    pub fn config(&self) -> EngineConfig {
        EngineConfig { frames: self.model, delta_t: self.delta_t }
    }

    /// Delta T (s) under this engine's model.
    pub fn delta_t_model(&self) -> DeltaTModel {
        self.delta_t
    }

    pub fn delta_t_seconds(&self, jd_ut: f64) -> f64 {
        time::delta_t_seconds_with(jd_ut, &self.delta_t)
    }

    fn memo_key(&self, jd: f64, is_ut: bool) -> MemoKey {
        let d = &self.delta_t;
        (jd.to_bits(), is_ut, self.model, d.spline, d.lunar_ndot.to_bits(), d.user_seconds.map(f64::to_bits))
    }

    pub fn open(path: &std::path::Path) -> Result<Engine> {
        Ok(Engine::new(Pack::open(path)?))
    }

    pub fn pack(&self) -> &Pack {
        &self.pack
    }

    pub fn frame_model(&self) -> FrameModel {
        self.model
    }

    fn prepare(&self, jd: f64, is_ut: bool, speed: bool) -> Prepared {
        let key = self.memo_key(jd, is_ut);
        let mut p = match LAST.get() {
            Some((k, p)) if k == key => p,
            _ => {
                let (ut, tt) = if is_ut {
                    (jd, time::ut_to_tt_with(jd, &self.delta_t))
                } else {
                    (time::tt_to_ut_with(jd, &self.delta_t), jd)
                };
                let epoch = Epoch::with_scales(ut, tt, self.model);
                let mut p = Prepared::new(epoch);
                p.geo[0] = ObsState::new(&self.pack, &epoch, &Origin::Geocentric).ok();
                p
            }
        };
        if speed && p.shifted.is_none() {
            let rates = FrameRates::new(p.epoch.jd_tt, self.model);
            let s = [p.epoch.shifted(-SPEED_STEP, &rates), p.epoch.shifted(SPEED_STEP, &rates)];
            p.geo[1] = ObsState::new(&self.pack, &s[0], &Origin::Geocentric).ok();
            p.geo[2] = ObsState::new(&self.pack, &s[1], &Origin::Geocentric).ok();
            p.shifted = Some(s);
        }
        LAST.set(Some((key, p)));
        p
    }

    /// Prepare a UT instant; `speed` also computes frame rates.
    pub fn prepare_ut(&self, jd_ut: f64, speed: bool) -> Prepared {
        self.prepare(jd_ut, true, speed)
    }

    pub fn prepare_tt(&self, jd_tt: f64, speed: bool) -> Prepared {
        self.prepare(jd_tt, false, speed)
    }

    /// An instant prepared in the JPL Horizons frame convention (IAU 1976
    /// precession, IAU 1980 nutation plus celestial pole offsets, no frame
    /// bias; `frames::iau1980`): the classic `SEFLG_JPLHOR` (`PoleOffsets::EopC04`)
    /// and `SEFLG_JPLHOR_APPROX` (`PoleOffsets::Model`) modes. Use with
    /// `calc_prepared` / `fixstar_prepared`; only "of date" outputs change.
    pub fn prepare_horizons(
        &self,
        jd: f64,
        is_ut: bool,
        speed: bool,
        offsets: crate::frames::iau1980::PoleOffsets,
    ) -> Prepared {
        let base = self.prepare(jd, is_ut, false);
        let mk = |dt: f64| {
            let (ut, tt) = (base.epoch.jd_ut + dt, base.epoch.jd_tt + dt);
            Epoch {
                jd_ut: ut,
                jd_tt: tt,
                jd_tdb: base.epoch.jd_tdb + dt,
                frames: crate::frames::iau1980::horizons_frames(tt, ut, offsets),
            }
        };
        let mut p = Prepared::new(mk(0.0));
        p.geo[0] = base.geo[0];
        if speed {
            let s = [mk(-SPEED_STEP), mk(SPEED_STEP)];
            p.geo[1] = ObsState::new(&self.pack, &s[0], &Origin::Geocentric).ok();
            p.geo[2] = ObsState::new(&self.pack, &s[1], &Origin::Geocentric).ok();
            p.shifted = Some(s);
        }
        p
    }

    /// Ayanamsa (deg) for positions in the given equinox convention: the
    /// mean value, plus Δψ for the true equinox (Swiss `get_ayanamsa_ex`).
    /// True Pushya is the apparent longitude of δ Cancri − 106°, from JPL
    /// Earth/Sun states.
    pub fn ayanamsa_deg(&self, mode: Ayanamsa, ep: &Epoch, true_equinox: bool) -> Result<f64> {
        self.ayanamsa_deg_bits(mode, SidBits(0), ep, true_equinox)
    }

    /// [`Engine::ayanamsa_deg`] with Swiss sidereal option bits. Star-anchored
    /// and galactic-centre modes use the catalogue / Sgr A* apparent place
    /// (parallax, deflection, aberration) from JPL Earth and Sun states.
    pub fn ayanamsa_deg_bits(
        &self,
        mode: Ayanamsa,
        bits: SidBits,
        ep: &Epoch,
        true_equinox: bool,
    ) -> Result<f64> {
        self.ayanamsa_deg_for(mode, bits, ep, true_equinox, Corrections::APPARENT)
    }

    /// Ayanamsa for a body computed with the given corrections. As in Swiss
    /// Ephemeris, a star-anchored or galactic-centre mode places its reference point
    /// with the body's own aberration and light-deflection settings, so
    /// `FLG_TRUEPOS` / `FLG_NOABERR` / `FLG_NOGDEFL` sidereal positions stay
    /// star-relative.
    fn ayanamsa_deg_for(
        &self,
        mode: Ayanamsa,
        bits: SidBits,
        ep: &Epoch,
        true_equinox: bool,
        corr: Corrections,
    ) -> Result<f64> {
        let mean = match ayanamsa::kind(mode, bits) {
            Kind::DeltaCnc(lon) => {
                let (e, ev) = self.pack.state(ephem::EARTH, ephem::SSB, ep.jd_tdb)?;
                let s = self.pack.position(ephem::SUN, ephem::SSB, ep.jd_tdb)?;
                let au = |v: Vec3| scale(v, 1.0 / AU_KM);
                let (l, _) = ayanamsa::star::delta_cnc_mean_from_states_with(
                    ep.jd_tt,
                    au(e),
                    au(ev),
                    au(s),
                    corr.deflection,
                    corr.aberration,
                );
                norm_deg(l - lon)
            }
            Kind::Star(name, lon) => {
                let rec = crate::stars::lookup(name)?.record;
                let w =
                    self.apparent_mean_ecl_with(crate::stars::barycentric_au(rec, ep.jd_tdb), ep, corr)?;
                norm_deg(w[1].atan2(w[0]) / DEG - lon)
            }
            Kind::GalCenter(lon) => {
                let w =
                    self.apparent_mean_ecl_with(ayanamsa::galactic::sgr_a_bary_au(ep.jd_tdb), ep, corr)?;
                norm_deg(w[1].atan2(w[0]) / DEG - lon)
            }
            Kind::GalCenterPolar(lon) => {
                let w =
                    self.apparent_mean_ecl_with(ayanamsa::galactic::sgr_a_bary_au(ep.jd_tdb), ep, corr)?;
                let (s, c) = ep.frames.eps_mean.sin_cos();
                norm_deg(ayanamsa::galactic::polar_longitude_deg(w, [0.0, s, c]) - lon)
            }
            Kind::Unsupported => {
                return Err(Error::InvalidArgument(format!(
                    "ayanamsa {:?} has no published definition available in FrenchEphemeris",
                    mode.name()
                )))
            }
            _ => ayanamsa::ayanamsa_mean_deg_bits(mode, ep.jd_tt, bits),
        };
        Ok(norm_deg(if true_equinox { mean + ep.frames.dpsi / DEG } else { mean }))
    }

    /// Sidereal (lon, lat) or sidereal (RA, Dec), degrees, of a direction
    /// given in the GCRS (used by fixed stars), with the same rules as planets:
    /// projection modes onto the plane of t0 / invariable plane, traditional
    /// modes tropical − ayanamsa, sidereal RA/Dec on the mean equator.
    /// `None` when `opt.sidereal` is not set.
    pub fn sidereal_direction(
        &self,
        u_gcrs: Vec3,
        ep: &Epoch,
        opt: &CalcOptions,
    ) -> Result<Option<(f64, f64)>> {
        let Some(mode) = opt.sidereal else { return Ok(None) };
        let proj = ayanamsa::projection(mode, opt.sid_bits);
        if proj != Projection::Traditional
            && !opt.j2000
            && (!opt.equatorial || proj.equator_of_t0().is_some())
        {
            let w = ephem::to_frame(ep, Frame::EclipticOfDate { true_equinox: false }, u_gcrs);
            let v = mat_t_vec(&ayanamsa::vondrak::ecliptic_of_date_matrix(ep.jd_tt), w);
            return Ok(Some(match proj.equator_of_t0() {
                Some(m) if opt.equatorial => {
                    let (a, d, _) = to_spherical(mat_vec(&m, v));
                    (norm_deg(a / DEG), d / DEG)
                }
                _ => proj.apply(v).expect("projection"),
            }));
        }
        let frame = match (opt.j2000, opt.equatorial) {
            (true, false) => Frame::EclipticJ2000,
            (true, true) => Frame::Icrf,
            (false, false) => Frame::EclipticOfDate { true_equinox: !opt.mean_equinox },
            (false, true) => Frame::EquatorialOfDate { true_equinox: false },
        };
        let (l, b, _) = to_spherical(ephem::to_frame(ep, frame, u_gcrs));
        let mut l = l / DEG;
        if !opt.equatorial {
            l = norm_deg(
                l - self.ayanamsa_deg_for(mode, opt.sid_bits, ep, !opt.mean_equinox, opt.corrections)?,
            );
        }
        Ok(Some((l, b / DEG)))
    }

    /// Apparent geocentric direction (unit, mean ecliptic & equinox of date) of a star at
    /// barycentric ICRS position `p_au`, with the body's deflection and aberration settings.
    fn apparent_mean_ecl_with(&self, p_au: Vec3, ep: &Epoch, body: Corrections) -> Result<Vec3> {
        let obs = ObsState::new(&self.pack, ep, &Origin::Geocentric)?;
        let corr = Corrections { light_time: true, deflection: body.deflection, aberration: body.aberration };
        let (u, _) = crate::stars::apparent_gcrs(p_au, &obs, corr);
        Ok(ephem::to_frame(ep, Frame::EclipticOfDate { true_equinox: false }, u))
    }

    /// Raw (lon, lat, dist_au) in degrees at one epoch, before sidereal shift.
    /// Also returns (τ, dτ/dt) for physical bodies (zero otherwise).
    fn raw(
        &self,
        body: Body,
        ep: &Epoch,
        geo: Option<&ObsState>,
        lt: LightTime,
        opt: &CalcOptions,
    ) -> Result<((f64, f64, f64), f64, f64)> {
        match body {
            // Lunar points are defined on the mean ecliptic of date (see `lunar_point`).
            Body::MeanNode => Ok((
                self.lunar_point((points::mean_node_deg(ep.jd_tt), 0.0, points::MEAN_NODE_DIST_AU), ep, opt),
                0.0,
                0.0,
            )),
            Body::MeanApogee => {
                let (l, b) = points::mean_apogee_deg(ep.jd_tt);
                Ok((self.lunar_point((l, b, points::MEAN_APOGEE_DIST_AU), ep, opt), 0.0, 0.0))
            }
            Body::TrueNode | Body::OscuApogee => {
                // Osculating elements of the geometric geocentric Moon at the
                // light-retarded time t − |r|/c, from an inertial state rotated
                // (position and velocity alike, no frame-rotation term) to the
                // mean ecliptic of date.
                let r = norm(self.pack.position(ephem::MOON, ephem::EARTH, ep.jd_tdb)?);
                let (p, v) = self.pack.state(ephem::MOON, ephem::EARTH, ep.jd_tdb - r / C_KM_DAY)?;
                let m = &ep.frames.ecl_mean;
                let (p, v) = (mat_vec(m, p), mat_vec(m, v));
                if body == Body::TrueNode {
                    let (l, d) = points::osculating_node_full(p, v, GM_EMB_KM3_DAY2);
                    Ok((self.lunar_point((l, 0.0, d / AU_KM), ep, opt), 0.0, 0.0))
                } else {
                    let (l, b, d) = points::osculating_apogee_deg(p, v, GM_EMB_KM3_DAY2);
                    Ok((self.lunar_point((l, b, d / AU_KM), ep, opt), 0.0, 0.0))
                }
            }
            Body::IntpApogee | Body::IntpPerigee => {
                let (l, b, d) = points::intp::natural_apsis(ep.jd_tt, body == Body::IntpApogee);
                Ok((self.lunar_point((l, b, d / AU_KM), ep, opt), 0.0, 0.0))
            }
            b if b.id() >= crate::fictitious::FICT_OFFSET && b.id() < 100 => {
                let el = crate::fictitious::elements(b.id()).ok_or(Error::UnknownBody(b.id()))?;
                let owned;
                let obs = match (geo, opt.origin) {
                    (Some(obs), Origin::Geocentric) => obs,
                    _ => {
                        owned = ObsState::new(&self.pack, ep, &opt.origin)?;
                        &owned
                    }
                };
                let model = ep.frames.model;
                let st = |t: f64| crate::fictitious::barycentric_state(&self.pack, el, t, model);
                let (v, tau, taudot) =
                    ephem::apparent_with(ep, obs, false, opt.corrections, lt, |t| Ok(st(t)?.0), st)?;
                let (l, b, r) = to_spherical(ephem::to_frame(ep, opt.frame(), v));
                Ok(((l / DEG, b / DEG, r / AU_KM), tau, taudot))
            }
            _ => {
                let naif = body.naif().expect("physical body");
                let owned;
                let obs = match (geo, opt.origin) {
                    (Some(obs), Origin::Geocentric) => obs,
                    _ => {
                        owned = ObsState::new(&self.pack, ep, &opt.origin)?;
                        &owned
                    }
                };
                let (v, tau, taudot) = ephem::apparent(&self.pack, ep, obs, naif, opt.corrections, lt)?;
                // ICRS: undo the frame bias so the precession acts on the ICRS vector itself.
                let v = if opt.icrs { mat_t_vec(ephem::frame_bias(), v) } else { v };
                let (l, b, r) = to_spherical(ephem::to_frame(ep, opt.frame(), v));
                Ok(((l / DEG, b / DEG, r / AU_KM), tau, taudot))
            }
        }
    }

    /// A lunar point given on the mean ecliptic of date, in the requested frame.
    /// Ecliptic-of-date output adds Δψ to the longitude only; equatorial and J2000
    /// output rotate the direction exactly as for a planet.
    fn lunar_point(&self, (l, b, d): (f64, f64, f64), ep: &Epoch, opt: &CalcOptions) -> (f64, f64, f64) {
        match opt.frame() {
            Frame::EclipticOfDate { .. } => (self.nutate(l, ep, opt), b, d),
            frame => {
                let (sl, cl) = (l * DEG).sin_cos();
                let (sb, cb) = (b * DEG).sin_cos();
                let u = mat_t_vec(&ep.frames.ecl_mean, [cb * cl, cb * sl, sb]);
                let (a, dec, _) = to_spherical(ephem::to_frame(ep, frame, u));
                (norm_deg(a / DEG), dec / DEG, d)
            }
        }
    }

    #[inline]
    fn nutate(&self, mean_lon: f64, ep: &Epoch, opt: &CalcOptions) -> f64 {
        if opt.mean_equinox {
            mean_lon
        } else {
            norm_deg(mean_lon + ep.frames.dpsi / DEG)
        }
    }

    fn at_epoch(
        &self,
        body: Body,
        ep: &Epoch,
        geo: Option<&ObsState>,
        lt: LightTime,
        opt: &CalcOptions,
    ) -> Result<((f64, f64, f64), f64, f64)> {
        if let Some(mode) = opt.sidereal {
            let proj = ayanamsa::projection(mode, opt.sid_bits);
            if proj != Projection::Traditional
                && !opt.j2000
                && (!opt.equatorial || proj.equator_of_t0().is_some())
            {
                // Projection modes: the body's mean-of-date ecliptic direction,
                // back to the J2000 mean frame, onto the plane of t0 / the
                // invariable plane (nutation never enters).
                let mopt = CalcOptions { mean_equinox: true, equatorial: false, sidereal: None, ..*opt };
                let ((l, b, r), tau, taudot) = self.raw(body, ep, geo, lt, &mopt)?;
                let (sl, cl) = (l * DEG).sin_cos();
                let (sb, cb) = (b * DEG).sin_cos();
                let v =
                    mat_t_vec(&ayanamsa::vondrak::ecliptic_of_date_matrix(ep.jd_tt), [cb * cl, cb * sl, sb]);
                let (l, b) = match proj.equator_of_t0() {
                    Some(m) if opt.equatorial => {
                        let (a, d, _) = to_spherical(mat_vec(&m, v));
                        (norm_deg(a / DEG), d / DEG)
                    }
                    _ => proj.apply(v).expect("projection"),
                };
                return Ok(((l, b, r), tau, taudot));
            }
        }
        // Swiss: sidereal RA/Dec (traditional modes) refer to the mean equator
        // and equinox of date — nutation is dropped (measured black-box).
        let mopt;
        let ropt = if opt.sidereal.is_some() && opt.equatorial && !opt.j2000 {
            mopt = CalcOptions { mean_equinox: true, ..*opt };
            &mopt
        } else {
            opt
        };
        let ((mut l, b, r), tau, taudot) = self.raw(body, ep, geo, lt, ropt)?;
        if let Some(mode) = opt.sidereal {
            if !opt.equatorial {
                l = norm_deg(
                    l - self.ayanamsa_deg_for(mode, opt.sid_bits, ep, !opt.mean_equinox, opt.corrections)?,
                );
            }
        }
        Ok(((l, b, r), tau, taudot))
    }

    /// Position of `body` at a prepared instant.
    pub fn calc_prepared(&self, body: Body, p: &Prepared, opt: &CalcOptions) -> Result<Position> {
        let ep = &p.epoch;
        let ((l, b, r), tau, taudot) = self.at_epoch(body, ep, p.geo[0].as_ref(), LightTime::Iterate, opt)?;
        let mut pos = Position { lon: l, lat: b, dist: r, ..Default::default() };
        if opt.speed {
            let h = SPEED_STEP;
            let s = match p.shifted {
                Some(s) => s,
                None => {
                    let rates = FrameRates::new(ep.jd_tt, ep.frames.model);
                    [ep.shifted(-h, &rates), ep.shifted(h, &rates)]
                }
            };
            let ((l0, b0, r0), ..) =
                self.at_epoch(body, &s[0], p.geo[1].as_ref(), LightTime::Known(tau - taudot * h), opt)?;
            let ((l1, b1, r1), ..) =
                self.at_epoch(body, &s[1], p.geo[2].as_ref(), LightTime::Known(tau + taudot * h), opt)?;
            pos.lon_speed = diff_deg(l1, l0) / (2.0 * h);
            pos.lat_speed = (b1 - b0) / (2.0 * h);
            pos.dist_speed = (r1 - r0) / (2.0 * h);
        }
        Ok(pos.formatted(opt))
    }

    /// Position of `body` at an epoch built by the caller.
    pub fn calc_epoch(&self, body: Body, ep: &Epoch, opt: &CalcOptions) -> Result<Position> {
        self.calc_prepared(body, &Prepared::new(*ep), opt)
    }

    pub fn calc_ut(&self, body: Body, jd_ut: f64, opt: &CalcOptions) -> Result<Position> {
        self.calc_prepared(body, &self.prepare_ut(jd_ut, opt.speed), opt)
    }

    pub fn calc_tt(&self, body: Body, jd_tt: f64, opt: &CalcOptions) -> Result<Position> {
        self.calc_prepared(body, &self.prepare_tt(jd_tt, opt.speed), opt)
    }

    /// All bodies at all UT instants; result is row-major `[time][body]`.
    /// Instant work (Delta T, precession-nutation) is shared across bodies,
    /// and large batches are split across CPU cores (the engine is immutable,
    /// so workers share it without locks).
    pub fn calc_batch_ut(&self, bodies: &[Body], jd_ut: &[f64], opt: &CalcOptions) -> Vec<Result<Position>> {
        const MIN_TIMES_PER_THREAD: usize = 64;
        let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
        let threads = threads.min(jd_ut.len() / MIN_TIMES_PER_THREAD).max(1);
        if threads == 1 {
            return self.calc_batch_serial(bodies, jd_ut, opt);
        }
        let chunk = jd_ut.len().div_ceil(threads);
        std::thread::scope(|sc| {
            let handles: Vec<_> = jd_ut
                .chunks(chunk)
                .map(|ts| sc.spawn(move || self.calc_batch_serial(bodies, ts, opt)))
                .collect();
            handles.into_iter().flat_map(|h| h.join().expect("batch worker panicked")).collect()
        })
    }

    fn calc_batch_serial(&self, bodies: &[Body], jd_ut: &[f64], opt: &CalcOptions) -> Vec<Result<Position>> {
        let mut out = Vec::with_capacity(bodies.len() * jd_ut.len());
        for &t in jd_ut {
            let p = self.prepare_ut(t, opt.speed);
            for &b in bodies {
                out.push(self.calc_prepared(b, &p, opt));
            }
        }
        out
    }

    /// House cusps for an instant and place (12-cusp systems; Sunshine uses the
    /// Sun's apparent declination; for 'G' and speeds see `house_cusps_ut`).
    /// `sidereal` shifts cusps and angles.
    pub fn houses_ut(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        system: char,
        sidereal: Option<Ayanamsa>,
    ) -> Result<Houses> {
        self.houses_ut_bits(jd_ut, lat_deg, lon_deg, system, sidereal, SidBits(0))
    }

    /// [`Engine::houses_ut`] with Swiss sidereal option bits. Traditional
    /// modes subtract the true ayanamsa; projection modes intersect the house
    /// circles with the plane of t0 / the invariable plane, as the oracle does
    /// (see `Engine::house_cusps_ut_bits`).
    pub fn houses_ut_bits(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        system: char,
        sidereal: Option<Ayanamsa>,
        bits: SidBits,
    ) -> Result<Houses> {
        if houses::canonical_system(system)? == 'G' {
            return Err(Error::InvalidArgument("Gauquelin sectors have 36 cusps: use house_cusps_ut".into()));
        }
        let (h, _) = self.house_cusps_ut_bits(
            jd_ut,
            lat_deg,
            lon_deg,
            system,
            sidereal.map(|m| (m, bits)),
            false,
            false,
        )?;
        let (c, a) = (&h.houses.cusps, h.houses.ascmc);
        Ok(Houses {
            cusps: std::array::from_fn(|i| c[i]),
            asc: a[0],
            mc: a[1],
            armc: a[2],
            vertex: a[3],
            equatorial_asc: a[4],
            co_asc_koch: a[5],
            co_asc_munkasey: a[6],
            polar_asc: a[7],
        })
    }
}
