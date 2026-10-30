// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Python bindings. Flags use the bit values of the classic astrology-engine
//! API so the FreeAstroAPI adapter can pass its existing flag words through;
//! every other setting (sidereal mode, observer) is an explicit argument.

use fe_core::ayanamsa::{Ayanamsa, SidBits};
use fe_core::context::{Body, CalcOptions, Engine as CoreEngine, EngineConfig, Position};
use fe_core::ephem::{Corrections, Observer, Origin};
use fe_core::{time, Error};
use numpy::{PyArray3, PyArrayMethods, PyReadonlyArray1};
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
mod stars;

mod eclipses;
mod events;
mod extras;
mod heliacal;
mod houses;

pub const FLG_HELCTR: i32 = 8;
pub const FLG_TRUEPOS: i32 = 16;
pub const FLG_J2000: i32 = 32;
pub const FLG_NONUT: i32 = 64;
pub const FLG_SPEED: i32 = 256;
pub const FLG_NOGDEFL: i32 = 512;
pub const FLG_NOABERR: i32 = 1024;
pub const FLG_EQUATORIAL: i32 = 2048;
pub const FLG_BARYCTR: i32 = 16384;
pub const FLG_TOPOCTR: i32 = 32 * 1024;
pub const FLG_SIDEREAL: i32 = 64 * 1024;
pub const FLG_SPEED3: i32 = 128;
pub const FLG_XYZ: i32 = 4096;
pub const FLG_RADIANS: i32 = 8192;
pub const FLG_ICRS: i32 = 128 * 1024;
pub const FLG_ASTROMETRIC: i32 = FLG_NOABERR | FLG_NOGDEFL;
/// JPL Horizons frames: IAU 1976/1980 + IERS EOP C04 pole offsets (= FLG_DPSIDEPS_1980).
pub const FLG_JPLHOR: i32 = 256 * 1024;
/// JPL Horizons frames with modelled pole offsets (no data file).
pub const FLG_JPLHOR_APPROX: i32 = 512 * 1024;

/// Horizons frame mode requested by `flags`, if any.
fn horizons_mode(flags: i32) -> Option<fe_core::frames::iau1980::PoleOffsets> {
    use fe_core::frames::iau1980::PoleOffsets;
    if flags & FLG_JPLHOR_APPROX != 0 {
        Some(PoleOffsets::Model)
    } else if flags & FLG_JPLHOR != 0 {
        Some(PoleOffsets::EopC04)
    } else {
        None
    }
}

/// Position with the Horizons frames for physical bodies (nodes and apsides
/// keep the standard method, as documented for the classic flag).
fn calc_any(
    e: &CoreEngine,
    body: Body,
    jd: f64,
    is_ut: bool,
    flags: i32,
    opt: &CalcOptions,
) -> Result<Position, Error> {
    match horizons_mode(flags) {
        Some(m) if body.naif().is_some() => {
            e.calc_prepared(body, &e.prepare_horizons(jd, is_ut, opt.speed, m), opt)
        }
        _ if is_ut => e.calc_ut(body, jd, opt),
        _ => e.calc_tt(body, jd, opt),
    }
}

pyo3::create_exception!(_native, EphemerisError, pyo3::exceptions::PyException);
pyo3::create_exception!(_native, OutOfRangeError, EphemerisError);
pyo3::create_exception!(_native, PolarHousesError, EphemerisError);

fn to_py(e: Error) -> PyErr {
    match e {
        Error::OutOfRange { .. } => OutOfRangeError::new_err(e.to_string()),
        Error::PolarHouses { .. } => PolarHousesError::new_err(e.to_string()),
        Error::UnknownBody(_) | Error::InvalidArgument(_) => PyValueError::new_err(e.to_string()),
        _ => EphemerisError::new_err(e.to_string()),
    }
}

/// Sidereal mode spec: a Swiss number (`"1"`, or `"257"` = mode + SIDBIT_*
/// bits), a name or identifier (`"lahiri"`, `"True Citra"`, `"kp"`), or
/// `"user:<t0>:<ayan_t0>[:ut]"` (SE_SIDM_USER), optionally followed by
/// `|<bits>` (e.g. `"lahiri|256"`).
fn parse_sidereal(spec: &str) -> PyResult<(Ayanamsa, SidBits)> {
    let bad = || PyValueError::new_err(format!("unknown ayanamsa {spec:?}"));
    let mut parts = spec.split('|');
    let head = parts.next().unwrap_or("").trim();
    let mut bits = 0u32;
    for b in parts {
        bits |= b.trim().parse::<u32>().map_err(|_| bad())?;
    }
    let low = head.to_ascii_lowercase();
    let mode = if let Some(rest) = low.strip_prefix("user:") {
        let f: Vec<&str> = rest.split(':').collect();
        if f.len() < 2 {
            return Err(bad());
        }
        let t0 = f[0].trim().parse::<f64>().map_err(|_| bad())?;
        let ayan_t0 = f[1].trim().parse::<f64>().map_err(|_| bad())?;
        Ayanamsa::User { t0, ayan_t0, t0_is_ut: f.get(2).is_some_and(|x| x.trim() == "ut") }
    } else if let Ok(n) = head.parse::<i64>() {
        bits |= (n as u32) & SidBits::ALL;
        Ayanamsa::from_swiss((n & 0xFF) as i32).ok_or_else(bad)?
    } else {
        Ayanamsa::parse(head).ok_or_else(bad)?
    };
    if bits & !SidBits::ALL != 0 {
        return Err(PyValueError::new_err(format!("unknown sidereal bits in {spec:?}")));
    }
    Ok((mode, SidBits(bits)))
}

/// `sidereal=` argument: a mode name / spec string or a Swiss SIDM number.
#[derive(FromPyObject)]
enum SidArg {
    Int(i64),
    Str(String),
}

impl SidArg {
    fn spec(&self) -> String {
        match self {
            SidArg::Int(n) => n.to_string(),
            SidArg::Str(s) => s.clone(),
        }
    }
}

fn options(flags: i32, sidereal: Option<&str>, topo: Option<(f64, f64, f64)>) -> PyResult<CalcOptions> {
    let origin = if flags & FLG_TOPOCTR != 0 {
        let (lon_deg, lat_deg, height_m) =
            topo.ok_or_else(|| PyValueError::new_err("FLG_TOPOCTR requires topo=(lon, lat, height_m)"))?;
        Origin::Topocentric(Observer { lon_deg, lat_deg, height_m })
    } else if flags & FLG_HELCTR != 0 {
        Origin::Heliocentric
    } else if flags & FLG_BARYCTR != 0 {
        Origin::Barycentric
    } else {
        Origin::Geocentric
    };
    let corrections = if flags & FLG_TRUEPOS != 0 {
        Corrections::GEOMETRIC
    } else {
        Corrections {
            light_time: true,
            deflection: flags & FLG_NOGDEFL == 0,
            aberration: flags & FLG_NOABERR == 0,
        }
    };
    let (sidereal, sid_bits) = if flags & FLG_SIDEREAL != 0 {
        let (m, b) = parse_sidereal(
            sidereal.ok_or_else(|| PyValueError::new_err("FLG_SIDEREAL requires sidereal=<mode>"))?,
        )?;
        (Some(m), b)
    } else {
        (None, SidBits(0))
    };
    Ok(CalcOptions {
        origin,
        corrections,
        equatorial: flags & FLG_EQUATORIAL != 0,
        mean_equinox: flags & FLG_NONUT != 0,
        j2000: flags & FLG_J2000 != 0,
        sidereal,
        sid_bits,
        // SPEED3 (three-position speed) is accepted as SPEED: ours is always the derivative.
        speed: flags & (FLG_SPEED | FLG_SPEED3) != 0,
        xyz: flags & FLG_XYZ != 0,
        radians: flags & FLG_RADIANS != 0,
        icrs: flags & FLG_ICRS != 0,
    })
}

/// Thread-safe engine handle; share one per process.
#[pyclass(frozen)]
struct Engine {
    inner: CoreEngine,
}

#[pymethods]
impl Engine {
    /// Delta T: by default the published SMH 2016 spline with the DE431
    /// tidal correction, which keeps historical charts as Swiss-based tools
    /// show them; `swiss_delta_t=False` selects the newer SMH 2020 addendum
    /// (physically preferred, moves the Moon by up to ~0.03° before 1300).
    /// `extra_packs` adds packs such as `fepk-minor-v2.bin`
    /// (minor bodies: ids 15-20, 146199 Eris, 10433 Eros).
    /// `swiss_sidereal_time=True` selects the long-term sidereal time documented
    /// for Swiss Ephemeris (Simon 1994 mean Sun + Vondrák 2011 precession) outside
    /// 1850–2050 for ARMC/houses, topocentric places and rise/set; inside that
    /// interval both models are IAU 2006. `False`: IAU 2006 at all dates.
    #[new]
    /// `star_compat`: optional fixed-star compatibility alias table (TSV path).
    /// `fallback_packs`: packs used whole for dates or bodies the main packs do not cover.
    #[pyo3(signature = (path, swiss_delta_t=true, extra_packs=None, swiss_sidereal_time=true, star_compat=None, fallback_packs=None))]
    fn new(
        path: &str,
        swiss_delta_t: bool,
        extra_packs: Option<Vec<String>>,
        swiss_sidereal_time: bool,
        star_compat: Option<String>,
        fallback_packs: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let mut paths = vec![path.to_string()];
        paths.extend(extra_packs.unwrap_or_default());
        let mut pack = fe_core::pack::Pack::open_many(&paths).map_err(to_py)?;
        if let Some(fb) = fallback_packs.filter(|v| !v.is_empty()) {
            pack = pack.with_fallback(fe_core::pack::Pack::open_many(&fb).map_err(to_py)?);
        }
        let mut cfg = if swiss_delta_t { EngineConfig::SWISS_DELTA_T } else { EngineConfig::DEFAULT };
        if swiss_sidereal_time {
            cfg.frames.sidereal = fe_core::frames::SiderealTime::SwissLongTerm;
        }
        let inner = CoreEngine::with_config(pack, cfg);
        let inner = match star_compat {
            Some(p) => inner.with_star_compat(std::path::Path::new(&p)).map_err(to_py)?,
            None => inner,
        };
        Ok(Engine { inner })
    }

    /// Delta T in days under this engine's model.
    fn deltat(&self, jd_ut: f64) -> f64 {
        self.inner.delta_t_seconds(jd_ut) / 86400.0
    }

    /// (lon, lat, dist, lon_speed, lat_speed, dist_speed) for a body at JD UT.
    #[pyo3(signature = (jd_ut, body, flags=FLG_SPEED, sidereal=None, topo=None))]
    fn calc_ut(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        sidereal: Option<SidArg>,
        topo: Option<(f64, f64, f64)>,
    ) -> PyResult<[f64; 6]> {
        let opt = options(flags, sidereal.map(|s| s.spec()).as_deref(), topo)?;
        if body == fe_core::extras::ECL_NUT {
            return Ok(self.inner.ecl_nut(jd_ut, true, &opt));
        }
        let body = Body::from_id(body).map_err(to_py)?;
        // A single position takes a few microseconds: releasing the GIL for it costs more than it
        // saves and makes threads contend. The batch API below runs without the GIL.
        calc_any(&self.inner, body, jd_ut, true, flags, &opt).map(Position::to_array).map_err(to_py)
    }

    /// Same as `calc_ut` with JD in TT.
    #[pyo3(signature = (jd_tt, body, flags=FLG_SPEED, sidereal=None, topo=None))]
    fn calc_tt(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        sidereal: Option<SidArg>,
        topo: Option<(f64, f64, f64)>,
    ) -> PyResult<[f64; 6]> {
        let opt = options(flags, sidereal.map(|s| s.spec()).as_deref(), topo)?;
        if body == fe_core::extras::ECL_NUT {
            return Ok(self.inner.ecl_nut(jd_tt, false, &opt));
        }
        let body = Body::from_id(body).map_err(to_py)?;
        calc_any(&self.inner, body, jd_tt, false, flags, &opt).map(Position::to_array).map_err(to_py)
    }

    /// Batch: array of shape (len(jd_ut), len(bodies), 6). Runs without the GIL.
    #[pyo3(signature = (jd_ut, bodies, flags=FLG_SPEED, sidereal=None, topo=None))]
    fn calc_batch_ut<'py>(
        &self,
        py: Python<'py>,
        jd_ut: PyReadonlyArray1<'py, f64>,
        bodies: Vec<i32>,
        flags: i32,
        sidereal: Option<SidArg>,
        topo: Option<(f64, f64, f64)>,
    ) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let opt = options(flags, sidereal.map(|s| s.spec()).as_deref(), topo)?;
        let bodies: Vec<Body> =
            bodies.into_iter().map(Body::from_id).collect::<Result<_, _>>().map_err(to_py)?;
        let times = jd_ut.as_slice()?.to_vec();
        let nb = bodies.len();
        let flat = py.detach(|| -> Result<Vec<f64>, Error> {
            let mut flat = Vec::with_capacity(times.len() * nb * 6);
            for r in self.inner.calc_batch_ut(&bodies, &times, &opt) {
                flat.extend_from_slice(&r?.to_array());
            }
            Ok(flat)
        });
        let flat = flat.map_err(to_py)?;
        let arr = numpy::PyArray1::from_vec(py, flat);
        arr.reshape([times.len(), nb, 6]).map_err(|e| PyIndexError::new_err(e.to_string()))
    }

    /// (cusps[12], ascmc[8]) with ascmc = asc, mc, armc, vertex, equ_asc,
    /// co_asc_koch, co_asc_munkasey, polar_asc.
    #[pyo3(signature = (jd_ut, lat, lon, hsys="P", sidereal=None))]
    fn houses_ut(
        &self,
        jd_ut: f64,
        lat: f64,
        lon: f64,
        hsys: &str,
        sidereal: Option<SidArg>,
    ) -> PyResult<([f64; 12], [f64; 8])> {
        let system = hsys.chars().next().ok_or_else(|| PyValueError::new_err("empty hsys"))?;
        let sid = sidereal.map(|s| parse_sidereal(&s.spec())).transpose()?;
        let (mode, bits) = match sid {
            Some((m, b)) => (Some(m), b),
            None => (None, SidBits(0)),
        };
        let h = self.inner.houses_ut_bits(jd_ut, lat, lon, system, mode, bits).map_err(to_py)?;
        Ok((
            h.cusps,
            [h.asc, h.mc, h.armc, h.vertex, h.equatorial_asc, h.co_asc_koch, h.co_asc_munkasey, h.polar_asc],
        ))
    }

    /// Ayanamsa in degrees at JD UT (`true_equinox` adds nutation in longitude).
    /// `mode` is a name, a Swiss SIDM number, or a spec (see `parse_sidereal`).
    #[pyo3(signature = (jd_ut, mode, true_equinox=true))]
    fn ayanamsa_ut(&self, jd_ut: f64, mode: SidArg, true_equinox: bool) -> PyResult<f64> {
        let (m, b) = parse_sidereal(&mode.spec())?;
        let ep = self.inner.prepare_ut(jd_ut, false).epoch;
        self.inner.ayanamsa_deg_bits(m, b, &ep, true_equinox).map_err(to_py)
    }

    /// Same as `ayanamsa_ut` with JD in TT.
    #[pyo3(signature = (jd_tt, mode, true_equinox=true))]
    fn ayanamsa_tt(&self, jd_tt: f64, mode: SidArg, true_equinox: bool) -> PyResult<f64> {
        let (m, b) = parse_sidereal(&mode.spec())?;
        let ep = self.inner.prepare_tt(jd_tt, false).epoch;
        self.inner.ayanamsa_deg_bits(m, b, &ep, true_equinox).map_err(to_py)
    }

    fn manifest(&self) -> String {
        self.inner.pack().manifest().to_string()
    }
}

/// Delta T in days (TT − UT), matching the classic `deltat` return unit.
#[pyfunction]
fn deltat(jd_ut: f64) -> f64 {
    time::delta_t_seconds(jd_ut) / 86400.0
}

#[pyfunction]
#[pyo3(signature = (year, month, day, hour=12.0, gregorian=true))]
fn julday(year: i32, month: u32, day: u32, hour: f64, gregorian: bool) -> f64 {
    time::julday(
        year,
        month,
        day,
        hour,
        if gregorian { time::Calendar::Gregorian } else { time::Calendar::Julian },
    )
}

#[pyfunction]
#[pyo3(signature = (jd, gregorian=true))]
fn revjul(jd: f64, gregorian: bool) -> (i32, u32, u32, f64) {
    time::revjul(jd, if gregorian { time::Calendar::Gregorian } else { time::Calendar::Julian })
}

/// Display name of a Swiss sidereal mode number (`get_ayanamsa_name`);
/// "" for unknown numbers and SE_SIDM_USER, like the classic function.
#[pyfunction]
fn get_ayanamsa_name(sid_mode: i64) -> &'static str {
    Ayanamsa::from_swiss((sid_mode & 0xFF) as i32).map_or("", Ayanamsa::name)
}

/// Snake-case identifiers of all sidereal modes, indexed by Swiss number.
#[pyfunction]
fn ayanamsa_idents() -> Vec<&'static str> {
    Ayanamsa::all().iter().map(|m| m.ident()).collect()
}

/// Greenwich apparent sidereal time in hours: IAU 2006 in 1850–2050 and the
/// Swiss-documented long-term model outside (as the engine's houses and Swiss's
/// `sidtime`); the IAU 2006 GMST polynomial diverges by hours at ±10,000 years.
#[pyfunction]
fn sidtime(jd_ut: f64) -> f64 {
    use fe_core::frames::{FrameModel, FrameSet, SiderealTime};
    let jd_tt = time::ut_to_tt_with(jd_ut, &time::DeltaTModel::SWISS_COMPAT);
    let model = FrameModel { sidereal: SiderealTime::SwissLongTerm, ..FrameModel::DEFAULT };
    FrameSet::with_model(jd_tt, model).gast(jd_ut) / fe_core::math::TAU * 24.0
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Engine>()?;
    events::register(m)?;
    stars::register(m)?;
    m.add_function(wrap_pyfunction!(deltat, m)?)?;
    m.add_function(wrap_pyfunction!(julday, m)?)?;
    m.add_function(wrap_pyfunction!(revjul, m)?)?;
    m.add_function(wrap_pyfunction!(sidtime, m)?)?;
    houses::register(m)?;
    m.add_function(wrap_pyfunction!(get_ayanamsa_name, m)?)?;
    m.add_function(wrap_pyfunction!(ayanamsa_idents, m)?)?;
    for (k, v) in [
        ("SIDBIT_ECL_T0", SidBits::ECL_T0),
        ("SIDBIT_SSY_PLANE", SidBits::SSY_PLANE),
        ("SIDBIT_USER_UT", SidBits::USER_UT),
        ("SIDBIT_ECL_DATE", SidBits::ECL_DATE),
        ("SIDBIT_NO_PREC_OFFSET", SidBits::NO_PREC_OFFSET),
        ("SIDBIT_PREC_ORIG", SidBits::PREC_ORIG),
    ] {
        m.add(k, v)?;
    }
    eclipses::register(m)?;
    extras::register(m)?;
    heliacal::register(m)?;
    m.add("EphemerisError", m.py().get_type::<EphemerisError>())?;
    m.add("OutOfRangeError", m.py().get_type::<OutOfRangeError>())?;
    m.add("PolarHousesError", m.py().get_type::<PolarHousesError>())?;
    m.add("DELTA_T_MODEL", time::DELTA_T_MODEL)?;
    m.add("FRAME_MODEL", fe_core::frames::FRAME_MODEL)?;
    for (k, v) in [
        ("FLG_HELCTR", FLG_HELCTR),
        ("FLG_TRUEPOS", FLG_TRUEPOS),
        ("FLG_J2000", FLG_J2000),
        ("FLG_NONUT", FLG_NONUT),
        ("FLG_SPEED", FLG_SPEED),
        ("FLG_NOGDEFL", FLG_NOGDEFL),
        ("FLG_NOABERR", FLG_NOABERR),
        ("FLG_EQUATORIAL", FLG_EQUATORIAL),
        ("FLG_BARYCTR", FLG_BARYCTR),
        ("FLG_TOPOCTR", FLG_TOPOCTR),
        ("FLG_SIDEREAL", FLG_SIDEREAL),
        ("FLG_SPEED3", FLG_SPEED3),
        ("FLG_XYZ", FLG_XYZ),
        ("FLG_RADIANS", FLG_RADIANS),
        ("FLG_ICRS", FLG_ICRS),
        ("FLG_JPLHOR", FLG_JPLHOR),
        ("FLG_DPSIDEPS_1980", FLG_JPLHOR),
        ("FLG_JPLHOR_APPROX", FLG_JPLHOR_APPROX),
        ("FLG_ASTROMETRIC", FLG_ASTROMETRIC),
    ] {
        m.add(k, v)?;
    }
    Ok(())
}
