// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Lunar eclipses visible from a place: the classic `lun_eclipse_when_loc`.
//!
//! Semantics (Swiss Ephemeris programmer's documentation §8.9 and release notes; return
//! layout and visibility rules confirmed black-box on the Swiss oracle):
//!
//! * The Moon is *visible* while its whole disc is above the horizon: from the rising to the
//!   setting of its lower limb, apparent (refracted, standard atmosphere at the site height),
//!   topocentric — `rise_set` with `DiscPoint::LowerLimb`.
//! * The search returns the next (previous) global lunar eclipse of which some phase is
//!   visible. Contacts (P1, U1, U2, U3, U4, P4) that occur while the Moon is not visible are
//!   reported as absent; moonrise / moonset inside [P1, P4] are reported.
//! * The reported maximum is the greatest eclipse when visible, else the visible instant
//!   nearest to it (the moonrise or moonset that bounds the visible part); the attributes are
//!   evaluated at that instant.
//! * The type is the deepest phase visible: total if the Moon is visible during totality,
//!   partial if during the umbral phase, else penumbral.
//!
//! All times are JD UT (engine Delta T).

use crate::context::{Body, Engine};
use crate::eclipses::{self, LunarEclipse, LunarKind, LunarLocal};
use crate::ephem::Observer;
use crate::events::{DiscPoint, RiseSetOptions};
use crate::search::brent;
use crate::Result;

/// A lunar eclipse as seen from one place (times JD UT).
#[derive(Debug, Clone, Copy)]
pub struct LunarLocalEclipse {
    /// The global eclipse (times JD TT).
    pub global: LunarEclipse,
    /// Deepest visible phase.
    pub kind: LunarKind,
    /// Greatest eclipse if visible, else the nearest visible instant.
    pub max: f64,
    pub max_visible: bool,
    /// Visible contacts, `None` when absent or not visible: P1, U1, U2, U3, U4, P4.
    pub p1: Option<f64>,
    pub u1: Option<f64>,
    pub u2: Option<f64>,
    pub u3: Option<f64>,
    pub u4: Option<f64>,
    pub p4: Option<f64>,
    /// Moonrise / moonset (lower limb) between P1 and P4.
    pub moonrise: Option<f64>,
    pub moonset: Option<f64>,
    /// Magnitudes and the Moon's horizontal coordinates at `max`.
    pub attrs: LunarLocal,
}

/// Visible intervals of the Moon (lower limb above the apparent horizon) within [t0, t1],
/// plus the moonrise / moonset inside the window.
///
/// The altitude is sampled every 5 minutes over the window (≤ 6 h) and each sign
/// change is refined; this needs no assumption about the next rising or setting,
/// which may lie days away (or not exist) at polar sites.
fn visibility(
    engine: &Engine,
    t0: f64,
    t1: f64,
    obs: &Observer,
) -> Result<(Vec<(f64, f64)>, Option<f64>, Option<f64>)> {
    let o = RiseSetOptions { disc: DiscPoint::LowerLimb, ..Default::default() };
    let f = |t: f64| engine.rise_set_altitude(Body::Moon, t, obs, &o);
    let n = (((t1 - t0) / (5.0 / 1440.0)).ceil() as usize).max(1);
    let (mut out, mut rise, mut set) = (Vec::new(), None, None);
    let (mut ta, mut fa) = (t0, f(t0)?);
    let mut start = (fa > 0.0).then_some(t0);
    for i in 1..=n {
        let tb = t0 + (t1 - t0) * i as f64 / n as f64;
        let fb = f(tb)?;
        if (fa > 0.0) != (fb > 0.0) {
            let t = brent(&f, ta, tb, fa, fb)?;
            if fb > 0.0 {
                rise.get_or_insert(t);
                start = Some(t);
            } else {
                set.get_or_insert(t);
                if let Some(a) = start.take() {
                    out.push((a, t));
                }
            }
        }
        ta = tb;
        fa = fb;
    }
    if let Some(a) = start {
        out.push((a, t1));
    }
    Ok((out, rise, set))
}

/// Next (previous) lunar eclipse with a visible phase at `obs`, after (before) `jd_ut`.
pub fn lunar_eclipse_local_search(
    engine: &Engine,
    jd_ut: f64,
    obs: &Observer,
    backward: bool,
) -> Result<LunarLocalEclipse> {
    let mut t = engine.ut_to_tt(jd_ut);
    loop {
        let e = eclipses::lunar_eclipse_search(engine, t, backward)?;
        t = e.max;
        if let Some(l) = local(engine, &e, obs)? {
            return Ok(l);
        }
    }
}

fn local(engine: &Engine, e: &LunarEclipse, obs: &Observer) -> Result<Option<LunarLocalEclipse>> {
    let ut = |t: f64| engine.tt_to_ut(t);
    let (p1, p4, max) = (ut(e.p1), ut(e.p4), ut(e.max));
    let (iv, moonrise, moonset) = visibility(engine, p1, p4, obs)?;
    if iv.is_empty() {
        return Ok(None);
    }
    let vis = |t: f64| iv.iter().any(|&(a, b)| t >= a && t <= b);
    let overlaps = |a: f64, b: f64| iv.iter().any(|&(c, d)| c.max(a) <= d.min(b));
    let seen = |t: Option<f64>| t.map(ut).filter(|&x| vis(x));
    let kind = match (e.kind, e.u1.zip(e.u4), e.u2.zip(e.u3)) {
        (LunarKind::Total, _, Some((a, b))) if overlaps(ut(a), ut(b)) => LunarKind::Total,
        (LunarKind::Total | LunarKind::Partial, Some((a, b)), _) if overlaps(ut(a), ut(b)) => {
            LunarKind::Partial
        }
        _ => LunarKind::Penumbral,
    };
    let max_visible = vis(max);
    let shown = if max_visible {
        max
    } else {
        iv.iter()
            .flat_map(|&(a, b)| [a, b])
            .min_by(|a, b| (a - max).abs().total_cmp(&(b - max).abs()))
            .unwrap_or(max)
    };
    let attrs = eclipses::lunar_local(engine, engine.ut_to_tt(shown), obs)?;
    Ok(Some(LunarLocalEclipse {
        global: *e,
        kind,
        max: shown,
        max_visible,
        p1: seen(Some(e.p1)),
        u1: seen(e.u1),
        u2: seen(e.u2),
        u3: seen(e.u3),
        u4: seen(e.u4),
        p4: seen(Some(e.p4)),
        moonrise,
        moonset,
        attrs,
    }))
}
