// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Lunar occultations of planets and fixed stars (`lun_occult_when_glob`,
//! `lun_occult_when_loc`, `lun_occult_where`).
//!
//! The solar-eclipse machinery (`eclipses.rs`) is reused unchanged with the
//! occulted body in place of the Sun: Besselian elements of the Moon's
//! "shadow" cast by the body, exact WGS84 Earth outline, topocentric contacts.
//! A star is a point source (radius 0; its direction from the catalogue with
//! space motion and parallax, placed at `STAR_DISTANCE_KM`); a planet is a disc
//! of its IAU 2015 equatorial radius, light-time corrected. For occultations
//! the Moon's limb is the mean radius k = 0.2725076 Earth radii for every
//! contact (the eclipse k1/k2 split only matters for the solar umbra).
//!
//! Candidates are the Moon–body conjunctions in ecliptic longitude (J2000
//! ecliptic of the GCRS vectors), filtered by the latitude difference, then
//! refined like a solar eclipse. All times are JD TT.

use crate::context::{Body, Engine};
use crate::eclipses::{
    self as ecl, apparent_altitude_deg, frames_at, newton_min, pole, root, shadow_event, solar_attrs, topo,
    Shape, Site, Sky, SolarAttrs, SolarEclipse, SolarKind, SolarWhere, Target, K_MOON_MEAN,
    SEARCH_MARGIN_DAYS, SUN_SHAPE,
};
use crate::ephem::Observer;
use crate::math::*;
use crate::search::brent;
use crate::stars::Star;
use crate::{Error, Result};

/// Mean (volumetric) radii, km: IAU WGCCRE 2015 (Archinal et al. 2018,
/// Celest. Mech. Dyn. Astron. 130:22, Table 1; Pluto from Table 4). The mean
/// radius stands for the planet's (slightly elliptical) apparent disc.
pub fn planet_radius_km(body: Body) -> Option<f64> {
    Some(match body {
        Body::Mercury => 2439.4,
        Body::Venus => 6051.8,
        Body::Mars => 3389.50,
        Body::Jupiter => 69911.0,
        Body::Saturn => 58232.0,
        Body::Uranus => 25362.0,
        Body::Neptune => 24622.0,
        Body::Pluto => 1188.3,
        _ => return None,
    })
}

/// True altitude of the Sun's centre (deg) above which an occultation counts
/// as in daylight: the standard sunrise/sunset altitude (upper limb with 34′
/// of horizontal refraction, Explanatory Supplement 3rd ed. §7.1).
pub const SUN_UP_TRUE_ALT_DEG: f64 = -0.8333;

/// The occulted object.
#[derive(Debug, Clone, Copy)]
pub enum Occulted {
    /// Sun (then the result is a solar eclipse) or a planet.
    Body(Body),
    Star(Star),
}

impl Occulted {
    fn target(&self) -> Result<(Target, Shape)> {
        let lunar = |r: f64| Shape { radius_km: r, k1: K_MOON_MEAN, k2: K_MOON_MEAN };
        match *self {
            Occulted::Body(Body::Sun) => Ok((Target::Sun, SUN_SHAPE)),
            Occulted::Body(b) => {
                let r = planet_radius_km(b).ok_or_else(|| {
                    Error::InvalidArgument(format!("occultations of body {} are not supported", b.id()))
                })?;
                Ok((Target::Body(b.naif().expect("planet")), lunar(r)))
            }
            Occulted::Star(s) => Ok((Target::Star(s.record), lunar(0.0))),
        }
    }

    pub fn is_star(&self) -> bool {
        matches!(self, Occulted::Star(_))
    }
}

/// Largest ecliptic-latitude difference (deg) at conjunction that can still
/// give an occultation somewhere: Moon parallax ≤ 1.0°, semidiameter ≤ 0.28°,
/// plus margin for the latitude change between conjunction and maximum.
const BETA_LIMIT_DEG: f64 = 1.75;
/// Stars beyond this |ecliptic latitude| are never occulted (Moon's
/// latitude ≤ 5.3°, plus parallax, semidiameter and the J2000-ecliptic tilt).
const STAR_LAT_LIMIT_DEG: f64 = 7.5;
/// Guard against endless loops for bodies that are never occulted.
const MAX_CONJUNCTIONS: usize = 20_000;

/// Longitude difference Moon − target (deg, wrapped to (−180, 180]), its rate
/// (deg/day), and the latitude difference (deg), J2000 ecliptic.
fn elong(sky: &Sky, t: f64) -> Result<(f64, f64, f64)> {
    let g = sky.at(t)?;
    let r = crate::frames::gcrs_to_ecliptic_j2000();
    let (m, vm, s) = (mat_vec(&r, g.m), mat_vec(&r, g.vm), mat_vec(&r, g.s));
    let vs = mat_vec(&r, g.vs);
    let (lm, bm, _) = to_spherical(m);
    let (ls, bs, _) = to_spherical(s);
    let (dlm, ..) = spherical_rates(m, vm);
    let (dls, ..) = spherical_rates(s, vs);
    Ok((diff_deg(lm / DEG, ls / DEG), (dlm - dls) / DEG, (bm - bs) / DEG))
}

/// Next (or previous) conjunction in longitude strictly after (before) `t0`.
fn conjunction(sky: &Sky, t0: f64, backward: bool) -> Result<f64> {
    let (d, rate, _) = elong(sky, t0)?;
    let dm = norm_deg(d);
    let mut t = if backward { t0 - dm / rate } else { t0 + (360.0 - dm) / rate };
    if !backward && dm == 0.0 {
        t = t0 + 360.0 / rate;
    }
    for _ in 0..12 {
        let (d, rate, _) = elong(sky, t)?;
        let step = d / rate;
        t -= step;
        if step.abs() < 1e-9 {
            break;
        }
    }
    Ok(t)
}

/// Result of a search that may stop after one conjunction (`ECL_ONE_TRY`).
#[derive(Debug, Clone, Copy)]
pub enum Search<T> {
    Found(T),
    /// No occultation at the tested conjunction (JD TT of that conjunction).
    NoneAt(f64),
}

fn star_check(obj: &Occulted, engine: &Engine, jd_tt: f64) -> Result<()> {
    if let Occulted::Star(s) = obj {
        let (target, shape) = obj.target()?;
        let sky = Sky::with_target(engine.pack(), jd_tt, target, shape);
        let g = sky.at(jd_tt)?;
        let (_, b, _) = to_spherical(mat_vec(&crate::frames::gcrs_to_ecliptic_j2000(), g.s));
        if (b / DEG).abs() > STAR_LAT_LIMIT_DEG {
            return Err(Error::InvalidArgument(format!(
                "{} (ecliptic latitude {:.1}°) is never occulted by the Moon",
                s.display(),
                b / DEG
            )));
        }
    }
    Ok(())
}

/// Next (previous) occultation of `obj` whose maximum is more than
/// `SEARCH_MARGIN_DAYS` after (before) `jd_tt`. With `one_try`, only the
/// first conjunction after (before) `jd_tt` is examined.
pub fn occultation_search(
    engine: &Engine,
    obj: &Occulted,
    jd_tt: f64,
    backward: bool,
    one_try: bool,
) -> Result<Search<SolarEclipse>> {
    star_check(obj, engine, jd_tt)?;
    let (target, shape) = obj.target()?;
    let sky = Sky::with_target(engine.pack(), jd_tt, target, shape);
    let accept =
        |t: f64| if backward { t < jd_tt - SEARCH_MARGIN_DAYS } else { t > jd_tt + SEARCH_MARGIN_DAYS };
    // Start one day early so that an event whose conjunction precedes the
    // start but whose maximum follows it is not missed.
    let mut t = if backward { jd_tt + 1.0 } else { jd_tt - 1.0 };
    for _ in 0..MAX_CONJUNCTIONS {
        let c = conjunction(&sky, t, backward)?;
        let (_, _, db) = elong(&sky, c)?;
        let ev = if db.abs() < BETA_LIMIT_DEG {
            let local_sky = Sky::with_target(engine.pack(), c, target, shape);
            shadow_event(engine, &local_sky, c, accept)?
        } else {
            None
        };
        if let Some(e) = ev {
            return Ok(Search::Found(e));
        }
        let past_start = if backward { c < jd_tt } else { c > jd_tt };
        if one_try && past_start {
            return Ok(Search::NoneAt(c));
        }
        // Conjunctions are ≥ 25 days apart: restart a day past this one.
        t = if backward { c - 1.0 } else { c + 1.0 };
    }
    Err(Error::InvalidArgument("no occultation found".into()))
}

/// `lun_occult_where` equivalent: where the body is centred on the lunar disc
/// (or closest to it) at `jd_tt`, with `jd_ut` for the Earth's rotation.
pub fn occultation_where_ut(engine: &Engine, obj: &Occulted, jd_tt: f64, jd_ut: f64) -> Result<SolarWhere> {
    let (target, shape) = obj.target()?;
    let sky = Sky::with_target(engine.pack(), jd_tt, target, shape);
    ecl::where_with_limits(&sky, &frames_at(engine, jd_tt), jd_tt, jd_ut)
}

/// Time (JD TT) when the event takes place at local apparent noon (body on
/// the meridian of the point of greatest eclipse), if it does between first
/// and last contact.
pub fn local_noon(engine: &Engine, obj: &Occulted, ev: &SolarEclipse) -> Result<Option<f64>> {
    let (target, shape) = obj.target()?;
    let sky = Sky::with_target(engine.pack(), ev.max, target, shape);
    ecl::local_noon_impl(engine, &sky, ev)
}

/// Local circumstances of an occultation; times JD TT. Contacts are
/// geometric; for a star C2 = C1 and C3 = C4.
#[derive(Debug, Clone, Copy)]
pub struct OccultLocal {
    /// `Total` (body fully hidden) or `Partial` (grazing a planet's disc).
    pub kind: SolarKind,
    pub max: f64,
    pub c1: f64,
    pub c2: Option<f64>,
    pub c3: Option<f64>,
    pub c4: f64,
    /// The occulted body rises / sets between C1 and C4.
    pub rise: Option<f64>,
    pub set: Option<f64>,
    pub visible: bool,
    pub max_visible: bool,
    pub c1_visible: bool,
    pub c2_visible: bool,
    pub c3_visible: bool,
    pub c4_visible: bool,
    /// The Sun is above the horizon at C1 / C4.
    pub c1_daylight: bool,
    pub c4_daylight: bool,
    /// Attributes at the maximum (azimuth/altitudes are the occulted body's).
    pub attrs: SolarAttrs,
}

fn local_for(
    engine: &Engine,
    obj: &Occulted,
    glob: &SolarEclipse,
    obs: &Observer,
) -> Result<Option<OccultLocal>> {
    let (target, shape) = obj.target()?;
    let sky = Sky::with_target(engine.pack(), glob.max, target, shape);
    let sun = Sky::new(engine.pack(), glob.max);
    let fr = frames_at(engine, glob.max);
    let kp = pole(&fr);
    let dt = glob.max - engine.tt_to_ut(glob.max);
    let site = Site::new(obs);
    let tp_at = |t: f64| {
        let g = sky.at(t)?;
        let (o, vo) = site.gcrs(&fr, t - dt);
        Ok::<_, Error>((topo(&g, o, vo, &shape), g, o))
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
    if tp.sep < tp.sm_u - tp.ss {
        kind = SolarKind::Total;
        if shape.radius_km == 0.0 {
            (c2, c3) = (Some(c1), Some(c4));
        } else {
            c2 = root(inner, max, c1)?;
            c3 = root(inner, max, c4)?;
        }
    }
    let alt = |sky: &Sky, t: f64| -> Result<f64> {
        let g = sky.at(t)?;
        let (o, _) = site.gcrs(&fr, t - dt);
        Ok(apparent_altitude_deg(site.alt_az(&fr, t - dt, sub(g.s, o)).0 / DEG, site.height_m))
    };
    let up = |t: f64| -> Result<bool> { Ok(alt(&sky, t)? > 0.0) };
    let sun_alt = |t: f64| -> Result<f64> {
        let g = sun.at(t)?;
        let (o, _) = site.gcrs(&fr, t - dt);
        Ok(site.alt_az(&fr, t - dt, sub(g.s, o)).0 / DEG)
    };
    // Rising / setting of the body inside [C1, C4] (5-minute scan + Brent).
    let (mut rise, mut set) = (None, None);
    let n = (((c4 - c1) / (5.0 / 1440.0)).ceil() as usize).max(1);
    let f = |t: f64| alt(&sky, t);
    let (mut ta, mut fa) = (c1, f(c1)?);
    for i in 1..=n {
        let tb = c1 + (c4 - c1) * i as f64 / n as f64;
        let fb = f(tb)?;
        if (fa > 0.0) != (fb > 0.0) {
            let r = brent(&f, ta, tb, fa, fb)?;
            if fb > 0.0 {
                rise.get_or_insert(r);
            } else {
                set.get_or_insert(r);
            }
        }
        ta = tb;
        fa = fb;
    }
    let (c1_visible, c4_visible, max_visible) = (up(c1)?, up(c4)?, up(max)?);
    let c2_visible = match c2 {
        Some(t) => up(t)?,
        None => false,
    };
    let c3_visible = match c3 {
        Some(t) => up(t)?,
        None => false,
    };
    let visible = c1_visible || c4_visible || max_visible || rise.is_some() || set.is_some();
    let b = ecl::besselian(&g, kp, &shape);
    let attrs = solar_attrs(&tp, &b, &site, o, &fr, max - dt);
    Ok(Some(OccultLocal {
        kind,
        max,
        c1,
        c2,
        c3,
        c4,
        rise,
        set,
        visible,
        max_visible,
        c1_visible,
        c2_visible,
        c3_visible,
        c4_visible,
        c1_daylight: sun_alt(c1)? > SUN_UP_TRUE_ALT_DEG,
        c4_daylight: sun_alt(c4)? > SUN_UP_TRUE_ALT_DEG,
        attrs,
    }))
}

/// Next (previous) occultation of `obj` visible at `obs` (the body above the
/// horizon at some time between C1 and C4) with its local maximum more than
/// `SEARCH_MARGIN_DAYS` after (before) `jd_tt`. With `one_try`, only the first
/// conjunction is examined.
pub fn occultation_local_search(
    engine: &Engine,
    obj: &Occulted,
    jd_tt: f64,
    obs: &Observer,
    backward: bool,
    one_try: bool,
) -> Result<Search<(OccultLocal, SolarEclipse)>> {
    let ok = |t: f64| if backward { t < jd_tt - SEARCH_MARGIN_DAYS } else { t > jd_tt + SEARCH_MARGIN_DAYS };
    // Global events are taken from slightly before the start: a local maximum
    // can follow the global one by up to a few hours.
    let mut t = if backward { jd_tt + 0.3 } else { jd_tt - 0.3 };
    for _ in 0..MAX_CONJUNCTIONS {
        let glob = match occultation_search(engine, obj, t, backward, one_try)? {
            Search::Found(g) => g,
            Search::NoneAt(c) => return Ok(Search::NoneAt(c)),
        };
        if let Some(loc) = local_for(engine, obj, &glob, obs)? {
            if ok(loc.max) && loc.visible {
                return Ok(Search::Found((loc, glob)));
            }
        }
        if one_try && (if backward { glob.max < jd_tt } else { glob.max > jd_tt }) {
            return Ok(Search::NoneAt(glob.max));
        }
        t = glob.max;
    }
    Err(Error::InvalidArgument("no local occultation found".into()))
}

/// Closest approach at a site to the Moon's centre near `t_guess` (JD TT):
/// (time, separation, Moon semidiameter, body semidiameter), angles in degrees.
/// A graze (limit line) is where separation = Moon SD ± body SD.
pub fn closest_approach(
    engine: &Engine,
    obj: &Occulted,
    t_guess: f64,
    obs: &Observer,
) -> Result<(f64, f64, f64, f64)> {
    let (target, shape) = obj.target()?;
    let sky = Sky::with_target(engine.pack(), t_guess, target, shape);
    let fr = frames_at(engine, t_guess);
    let dt = t_guess - engine.tt_to_ut(t_guess);
    let site = Site::new(obs);
    let tp_at = |t: f64| -> Result<ecl::Topo> {
        let g = sky.at(t)?;
        let (o, vo) = site.gcrs(&fr, t - dt);
        Ok(topo(&g, o, vo, &shape))
    };
    let t = newton_min(
        |t| {
            let tp = tp_at(t)?;
            Ok((tp.w, tp.dw))
        },
        t_guess,
    )?;
    let tp = tp_at(t)?;
    Ok((t, tp.sep / DEG, tp.sm / DEG, tp.ss / DEG))
}

/// Local circumstances of a given global occultation (no visibility filter).
pub fn occultation_local(
    engine: &Engine,
    obj: &Occulted,
    glob: &SolarEclipse,
    obs: &Observer,
) -> Result<Option<OccultLocal>> {
    local_for(engine, obj, glob, obs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radii_present_for_planets() {
        for b in [
            Body::Mercury,
            Body::Venus,
            Body::Mars,
            Body::Jupiter,
            Body::Saturn,
            Body::Uranus,
            Body::Neptune,
            Body::Pluto,
        ] {
            assert!(planet_radius_km(b).unwrap() > 1000.0);
        }
        assert!(planet_radius_km(Body::Moon).is_none());
    }
}
