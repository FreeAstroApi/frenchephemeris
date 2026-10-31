// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! WebAssembly bindings (npm package `frenchephemeris`). Flat functions returning numbers and
//! `Float64Array`s; `js/index.js` builds the public JavaScript API (body names, objects) on top.
//!
//! The engine uses the same models as the Python package's default engine
//! (`frenchephemeris.open_default()`): Swiss-compatible Delta T and long-term sidereal time.
//! Packs are passed as bytes, since a browser has no file system.

use fe_core::ayanamsa::Ayanamsa;
use fe_core::context::{Body, CalcOptions, Engine as CoreEngine, EngineConfig};
use fe_core::ephem::{Corrections, Observer, Origin};
use fe_core::pack::Pack;
use fe_core::space::{self, Site, SpaceObserver, World, Zodiac};
use fe_core::time::{self, Calendar};
use wasm_bindgen::prelude::*;

/// Output and origin options of `calcUt` (bit set).
pub const EQUATORIAL: u32 = 1;
pub const HELIOCENTRIC: u32 = 2;
pub const TRUE_POSITION: u32 = 4;
pub const NO_NUTATION: u32 = 8;
pub const J2000: u32 = 16;

fn err(e: impl std::fmt::Display) -> JsError {
    JsError::new(&e.to_string())
}

fn sidereal(name: Option<String>) -> Result<Option<Ayanamsa>, JsError> {
    match name {
        None => Ok(None),
        Some(s) if s.is_empty() => Ok(None),
        Some(s) => {
            Ayanamsa::parse(&s).map(Some).ok_or_else(|| JsError::new(&format!("unknown ayanamsa {s:?}")))
        }
    }
}

fn one_char(s: &str) -> Result<char, JsError> {
    let mut it = s.chars();
    match (it.next(), it.next()) {
        (Some(c), None) => Ok(c),
        _ => Err(JsError::new("house system must be one character, such as \"P\"")),
    }
}

fn zodiac(name: &str) -> Result<Zodiac, JsError> {
    match name {
        "earth" | "" => Ok(Zodiac::Earth),
        "native" => Ok(Zodiac::Native),
        other => Err(JsError::new(&format!("unknown zodiac {other:?}: use \"earth\" or \"native\""))),
    }
}

/// An ephemeris engine over one or more data packs.
#[wasm_bindgen]
pub struct Engine {
    inner: CoreEngine,
}

#[wasm_bindgen]
impl Engine {
    /// `main`: the planetary pack; `minor`: the optional minor-body pack (Chiron, Ceres, ...).
    #[wasm_bindgen(constructor)]
    pub fn new(main: Vec<u8>, minor: Option<Vec<u8>>) -> Result<Engine, JsError> {
        let mut files = vec![main];
        files.extend(minor);
        let pack = Pack::from_bytes_many(files).map_err(err)?;
        let mut cfg = EngineConfig::SWISS_DELTA_T;
        cfg.frames.sidereal = fe_core::frames::SiderealTime::SwissLongTerm;
        Ok(Engine { inner: CoreEngine::with_config(pack, cfg) })
    }

    /// Position of `body` (classic id: 0 Sun, 1 Moon, ... 9 Pluto, 11 true node, 14 Earth,
    /// 15 Chiron, ...) at JD UT: [lon, lat, dist AU, lon speed, lat speed, dist speed].
    /// `topo`: [east longitude, latitude, height m] for a topocentric place.
    #[wasm_bindgen(js_name = calcUt)]
    pub fn calc_ut(
        &self,
        jd_ut: f64,
        body: i32,
        flags: u32,
        sidereal_mode: Option<String>,
        topo: Option<Vec<f64>>,
    ) -> Result<Vec<f64>, JsError> {
        let body = Body::from_id(body).map_err(err)?;
        let origin = match (&topo, flags & HELIOCENTRIC != 0) {
            (_, true) => Origin::Heliocentric,
            (Some(t), false) if t.len() >= 2 => Origin::Topocentric(Observer {
                lon_deg: t[0],
                lat_deg: t[1],
                height_m: t.get(2).copied().unwrap_or(0.0),
            }),
            _ => Origin::Geocentric,
        };
        let opt = CalcOptions {
            origin,
            corrections: if flags & TRUE_POSITION != 0 {
                Corrections::GEOMETRIC
            } else {
                Corrections::APPARENT
            },
            equatorial: flags & EQUATORIAL != 0,
            mean_equinox: flags & NO_NUTATION != 0,
            j2000: flags & J2000 != 0,
            sidereal: sidereal(sidereal_mode)?,
            ..CalcOptions::default()
        };
        Ok(self.inner.calc_ut(body, jd_ut, &opt).map_err(err)?.to_array().to_vec())
    }

    /// House cusps and angles: 12 cusps then Asc, MC, ARMC, Vertex, equatorial Asc,
    /// co-Asc (Koch), co-Asc (Munkasey), polar Asc.
    #[wasm_bindgen(js_name = housesUt)]
    pub fn houses_ut(
        &self,
        jd_ut: f64,
        lat: f64,
        lon: f64,
        system: &str,
        sidereal_mode: Option<String>,
    ) -> Result<Vec<f64>, JsError> {
        let h = self
            .inner
            .houses_ut(jd_ut, lat, lon, one_char(system)?, sidereal(sidereal_mode)?)
            .map_err(err)?;
        let mut out = h.cusps.to_vec();
        out.extend_from_slice(&h.ascmc());
        Ok(out)
    }

    /// Ayanamsa in degrees at JD UT (mean, without nutation, unless `true_equinox`).
    #[wasm_bindgen(js_name = ayanamsaUt)]
    pub fn ayanamsa_ut(&self, jd_ut: f64, mode: &str, true_equinox: bool) -> Result<f64, JsError> {
        let m = Ayanamsa::parse(mode).ok_or_else(|| JsError::new(&format!("unknown ayanamsa {mode:?}")))?;
        let ep = self.inner.prepare_ut(jd_ut, false).epoch;
        self.inner.ayanamsa_deg(m, &ep, true_equinox).map_err(err)
    }

    /// Fixed star position at JD UT (same layout as `calcUt`).
    #[wasm_bindgen(js_name = fixstarUt)]
    pub fn fixstar_ut(
        &self,
        name: &str,
        jd_ut: f64,
        flags: u32,
        sidereal_mode: Option<String>,
    ) -> Result<Vec<f64>, JsError> {
        let opt = CalcOptions {
            equatorial: flags & EQUATORIAL != 0,
            mean_equinox: flags & NO_NUTATION != 0,
            j2000: flags & J2000 != 0,
            sidereal: sidereal(sidereal_mode)?,
            ..CalcOptions::default()
        };
        Ok(self.inner.fixstar_ut(name, jd_ut, &opt).map_err(err)?.1.to_array().to_vec())
    }

    /// The `"Name,designation"` string a star name resolves to.
    #[wasm_bindgen(js_name = fixstarName)]
    pub fn fixstar_name(&self, name: &str) -> Result<String, JsError> {
        Ok(self.inner.fixstar_lookup(name).map_err(err)?.1)
    }

    /// V magnitude of a fixed star.
    #[wasm_bindgen(js_name = fixstarMag)]
    pub fn fixstar_mag(&self, name: &str) -> Result<f64, JsError> {
        Ok(self.inner.fixstar_mag(name).map_err(err)?.1)
    }

    /// Delta T (TT - UT) in seconds at JD UT.
    #[wasm_bindgen(js_name = deltaT)]
    pub fn delta_t(&self, jd_ut: f64) -> f64 {
        self.inner.delta_t_seconds(jd_ut)
    }

    /// JD TT for JD UT, with this engine's Delta T.
    #[wasm_bindgen(js_name = utToTt)]
    pub fn ut_to_tt(&self, jd_ut: f64) -> f64 {
        self.inner.ut_to_tt(jd_ut)
    }

    /// JD UT for JD TT, with this engine's Delta T.
    #[wasm_bindgen(js_name = ttToUt)]
    pub fn tt_to_ut(&self, jd_tt: f64) -> f64 {
        self.inner.tt_to_ut(jd_tt)
    }

    /// Position of `body` seen from a site on another world (`world` name, east longitude,
    /// latitude): [lon, lat, dist AU, lon speed, lat speed, dist speed] in `zodiac`
    /// ("earth" or "native").
    #[wasm_bindgen(js_name = spacePosition)]
    pub fn space_position(
        &self,
        jd_ut: f64,
        body: i32,
        world: &str,
        lon: f64,
        lat: f64,
        zodiac_name: &str,
    ) -> Result<Vec<f64>, JsError> {
        let site = Site::new(World::from_name(world).map_err(err)?, lon, lat);
        let p = self
            .inner
            .space_position(
                Body::from_id(body).map_err(err)?,
                jd_ut,
                &SpaceObserver::Surface(site),
                zodiac(zodiac_name)?,
                Corrections::APPARENT,
            )
            .map_err(err)?;
        Ok(vec![p.lon, p.lat, p.dist_au, p.lon_speed, p.lat_speed, p.dist_speed])
    }

    /// Houses on another world: 12 cusps, then Asc, MC, sidereal angle, spin latitude, obliquity.
    #[wasm_bindgen(js_name = spaceHouses)]
    pub fn space_houses(
        &self,
        jd_ut: f64,
        world: &str,
        lon: f64,
        lat: f64,
        system: &str,
        zodiac_name: &str,
    ) -> Result<Vec<f64>, JsError> {
        let site = Site::new(World::from_name(world).map_err(err)?, lon, lat);
        let h =
            self.inner.space_houses(jd_ut, &site, one_char(system)?, zodiac(zodiac_name)?).map_err(err)?;
        let mut out = h.houses.cusps.to_vec();
        out.extend([h.houses.asc, h.houses.mc, h.sidereal_angle, h.spin_latitude, h.obliquity]);
        Ok(out)
    }

    /// Manifest text of the loaded packs (sources, coverage, checksums).
    pub fn manifest(&self) -> String {
        self.inner.pack().manifests().join("\n")
    }
}

/// Julian day for a calendar date and decimal hour (Gregorian, or Julian when `gregorian` is false).
#[wasm_bindgen]
pub fn julday(year: i32, month: u32, day: u32, hour: f64, gregorian: bool) -> f64 {
    time::julday(year, month, day, hour, if gregorian { Calendar::Gregorian } else { Calendar::Julian })
}

/// Calendar date for a Julian day: [year, month, day, hour].
#[wasm_bindgen]
pub fn revjul(jd: f64, gregorian: bool) -> Vec<f64> {
    let (y, m, d, h) = time::revjul(jd, if gregorian { Calendar::Gregorian } else { Calendar::Julian });
    vec![y as f64, m as f64, d as f64, h]
}

/// Mars Sol Date for JD TT.
#[wasm_bindgen(js_name = marsSolDate)]
pub fn mars_sol_date(jd_tt: f64) -> f64 {
    space::mars_sol_date(jd_tt)
}

/// JD TT for a Mars Sol Date.
#[wasm_bindgen(js_name = jdTtFromMarsSolDate)]
pub fn jd_tt_from_mars_sol_date(msd: f64) -> f64 {
    space::jd_tt_from_mars_sol_date(msd)
}

/// Local mean solar time on Mars (hours) at east longitude `lon`.
#[wasm_bindgen(js_name = marsLocalMeanTime)]
pub fn mars_local_mean_time(jd_tt: f64, lon: f64) -> f64 {
    space::mars_local_mean_time(jd_tt, lon)
}
