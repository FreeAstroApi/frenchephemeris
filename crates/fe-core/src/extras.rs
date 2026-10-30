// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Smaller classic-API capabilities built on the engine: the obliquity/nutation pseudo-body
//! (`ECL_NUT`), planet-centric positions (`calc_pctr`), heliocentric longitude crossings
//! (`helio_cross`), Moon node crossings (`mooncross_node`), sidereal time from given
//! obliquity and nutation (`sidtime0`) and UTC conversions with the engine's Delta T.
//! Semantics follow the Swiss Ephemeris programmer's documentation; return shapes are
//! validated black-box against the Swiss oracle in `tests/quickwins.rs`.

use crate::context::{Body, CalcOptions, Engine};
use crate::ephem::Origin;
use crate::math::*;
use crate::search::brent;
use crate::time::Calendar;
use crate::utc;
use crate::{Error, Result};

/// Body id of the obliquity / nutation pseudo-body.
pub const ECL_NUT: i32 = -1;

impl Engine {
    /// The `ECL_NUT` pseudo-body at a JD (UT if `is_ut`, else TT): `[true obliquity, mean
    /// obliquity, nutation in longitude, nutation in obliquity, 0, 0]` in degrees (radians with
    /// `opt.radians`), from this engine's frame model. As in the original, cartesian or
    /// equatorial output is not defined for it and gives zeros; the other options are ignored.
    pub fn ecl_nut(&self, jd: f64, is_ut: bool, opt: &CalcOptions) -> [f64; 6] {
        if opt.xyz || opt.equatorial {
            return [0.0; 6];
        }
        let p = if is_ut { self.prepare_ut(jd, false) } else { self.prepare_tt(jd, false) };
        let f = &p.epoch.frames;
        let k = if opt.radians { 1.0 } else { 1.0 / DEG };
        [f.eps_true * k, f.eps_mean * k, f.dpsi * k, f.deps * k, 0.0, 0.0]
    }

    /// Position of `target` seen from the centre of `center` (JD TT) — the classic
    /// `calc_pctr`. Both must be physical bodies of the pack and differ. Light time,
    /// solar light deflection and aberration (centre's barycentric velocity) follow
    /// `opt.corrections`; frames, sidereal and output options as for `calc_tt`.
    /// `opt.origin` is ignored.
    pub fn calc_pctr(
        &self,
        target: Body,
        center: Body,
        jd_tt: f64,
        opt: &CalcOptions,
    ) -> Result<crate::context::Position> {
        if target == center {
            return Err(Error::InvalidArgument(format!(
                "calc_pctr: target and center must differ ({})",
                target.id()
            )));
        }
        let need = |b: Body| {
            b.naif().ok_or(Error::InvalidArgument(format!("calc_pctr needs physical bodies, got {b:?}")))
        };
        need(target)?;
        let c = need(center)?;
        let opt = CalcOptions { origin: Origin::Body(c), ..*opt };
        self.calc_tt(target, jd_tt, &opt)
    }

    /// Next (or, `backward`, previous) time (JD TT) after `jd_tt` when the heliocentric
    /// longitude of `body` (per `opt`: apparent/true, tropical/sidereal) equals `lon_deg` —
    /// the classic `helio_cross`. Bodies with a heliocentric orbit only (planets, Earth,
    /// minor bodies). Uses the monotonic heliocentric motion: the scan accumulates the
    /// longitude travelled in steps of ≤ 5° at the starting speed (≤ 70° anywhere on the
    /// most eccentric orbits served) and refines the bracketing step with Brent.
    pub fn helio_cross(
        &self,
        body: Body,
        lon_deg: f64,
        jd_tt: f64,
        opt: &CalcOptions,
        backward: bool,
    ) -> Result<f64> {
        if matches!(body, Body::Sun | Body::Moon) || body.naif().is_none() {
            return Err(Error::InvalidArgument(format!("helio_cross: {body:?} has no heliocentric orbit")));
        }
        let opt = CalcOptions {
            origin: Origin::Heliocentric,
            equatorial: false,
            speed: true,
            xyz: false,
            radians: false,
            ..*opt
        };
        let lon = |t: f64| self.calc_tt(body, t, &opt).map(|p| p.lon);
        let p0 = self.calc_tt(body, jd_tt, &opt)?;
        let dir = if backward { -1.0 } else { 1.0 };
        let need = if backward { norm_deg(p0.lon - lon_deg) } else { norm_deg(lon_deg - p0.lon) };
        let need = if need == 0.0 { 360.0 } else { need };
        let step = dir * 5.0 / p0.lon_speed.abs().max(1e-6);
        let g = |t: f64| lon(t).map(|l| diff_deg(l, lon_deg));
        let (mut t, mut prev, mut acc) = (jd_tt, p0.lon, 0.0);
        for _ in 0..100_000 {
            let t1 = t + step;
            let l = lon(t1)?;
            acc += dir * diff_deg(l, prev);
            if acc >= need {
                let (ga, gb) = (g(t)?, g(t1)?);
                let (a, b, fa, fb) = if t < t1 { (t, t1, ga, gb) } else { (t1, t, gb, ga) };
                return brent(g, a, b, fa, fb);
            }
            t = t1;
            prev = l;
        }
        Err(Error::InvalidArgument("helio_cross: no crossing found".into()))
    }

    /// `helio_cross` with UT input and result.
    pub fn helio_cross_ut(
        &self,
        body: Body,
        lon_deg: f64,
        jd_ut: f64,
        opt: &CalcOptions,
        backward: bool,
    ) -> Result<f64> {
        let tt = self.helio_cross(body, lon_deg, self.ut_to_tt(jd_ut), opt, backward)?;
        Ok(self.tt_to_ut(tt))
    }

    /// Next crossing of the Moon through the ecliptic (geocentric latitude 0, per `opt`) after
    /// `jd_tt`: (JD TT, longitude, latitude) in degrees — the classic `mooncross_node`.
    pub fn mooncross_node(&self, jd_tt: f64, opt: &CalcOptions) -> Result<(f64, f64, f64)> {
        let opt = CalcOptions { equatorial: false, speed: false, xyz: false, radians: false, ..*opt };
        let g = |t: f64| self.calc_tt(Body::Moon, t, &opt).map(|p| p.lat);
        // Latitude changes sign every ~13.6 days; 1-day steps cannot skip a crossing.
        let mut a = jd_tt + 1e-7;
        let mut fa = g(a)?;
        for _ in 0..40 {
            let b = a + 1.0;
            let fb = g(b)?;
            if fa == 0.0 || fa.signum() != fb.signum() {
                let t = brent(g, a, b, fa, fb)?;
                let p = self.calc_tt(Body::Moon, t, &opt)?;
                return Ok((t, p.lon, p.lat));
            }
            a = b;
            fa = fb;
        }
        Err(Error::InvalidArgument("mooncross_node: no crossing found".into()))
    }

    /// `mooncross_node` with UT input and result time.
    pub fn mooncross_node_ut(&self, jd_ut: f64, opt: &CalcOptions) -> Result<(f64, f64, f64)> {
        let (t, l, b) = self.mooncross_node(self.ut_to_tt(jd_ut), opt)?;
        Ok((self.tt_to_ut(t), l, b))
    }

    /// Sidereal time (hours) at JD UT from a given obliquity and nutation in longitude
    /// (degrees): GMST of this engine's model + Δψ cos ε — the classic `sidtime0`.
    pub fn sidtime0(&self, jd_ut: f64, eps_deg: f64, nut_deg: f64) -> f64 {
        let gmst = self.prepare_ut(jd_ut, false).epoch.frames.gmst(jd_ut);
        norm_rad(gmst + nut_deg * DEG * (eps_deg * DEG).cos()) / TAU * 24.0
    }

    /// JD TT from JD UT with this engine's Delta T.
    pub fn ut_to_tt(&self, jd_ut: f64) -> f64 {
        crate::time::ut_to_tt_with(jd_ut, &self.delta_t_model())
    }

    /// JD UT from JD TT with this engine's Delta T.
    pub fn tt_to_ut(&self, jd_tt: f64) -> f64 {
        crate::time::tt_to_ut_with(jd_tt, &self.delta_t_model())
    }

    /// Lunar tidal acceleration (″/cy²) of this engine's Delta T model (`get_tid_acc`).
    pub fn tidal_acceleration(&self) -> f64 {
        self.delta_t_model().lunar_ndot
    }

    /// UTC → (JD TT, JD UT1) with this engine's Delta T (`utc_to_jd`).
    #[allow(clippy::too_many_arguments)]
    pub fn utc_to_jd(
        &self,
        y: i32,
        m: i32,
        d: i32,
        h: i32,
        mi: i32,
        s: f64,
        cal: Calendar,
    ) -> Result<(f64, f64)> {
        utc::utc_to_jd(y, m, d, h, mi, s, cal, &self.delta_t_model())
    }

    /// JD TT → UTC (`jdet_to_utc`).
    pub fn jdet_to_utc(&self, jd_tt: f64, cal: Calendar) -> (i32, i32, i32, i32, i32, f64) {
        utc::jdet_to_utc(jd_tt, cal, &self.delta_t_model())
    }

    /// JD UT1 → UTC (`jdut1_to_utc`).
    pub fn jdut1_to_utc(&self, jd_ut: f64, cal: Calendar) -> (i32, i32, i32, i32, i32, f64) {
        utc::jdut1_to_utc(jd_ut, cal, &self.delta_t_model())
    }
}
