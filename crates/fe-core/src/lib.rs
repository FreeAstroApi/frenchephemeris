// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! FrenchEphemeris core: a JPL-derived ephemeris and astrology-geometry kernel.
//!
//! Design rules (see docs/ARCHITECTURE.md):
//! - No global mutable state. Every setting lives in an immutable context or
//!   is passed as an argument, so every function here is thread-safe.
//! - Angles are radians internally; degrees only at the public adapter edge
//!   (`context` and the Python bindings).
//! - Time arguments are Julian dates; the suffix says the scale
//!   (`jd_ut`, `jd_tt`, `jd_tdb`).

pub mod error;
pub mod math;

pub mod ayanamsa; // sidereal zodiac modes
pub mod context; // immutable calculation context + body ids
pub mod eclipses; // solar and lunar eclipses (Besselian elements, Danjon shadow)
pub mod ephem; // apparent / true / topocentric positions pipeline
pub mod events; // rise/set/transit, phenomena, crossings, equation of time
pub mod extras; // ECL_NUT, calc_pctr, helio_cross, mooncross_node, sidtime0, UTC via the engine
pub mod fictitious; // hypothetical bodies from published orbital elements
pub mod frames; // precession, nutation, obliquity, sidereal time
pub mod heliacal;
pub mod houses; // house systems
pub mod lunar_local; // lunar eclipses visible from a place
pub mod occultations; // lunar occultations of planets and stars
pub mod pack; // JPL-derived coefficient pack reader
pub mod points; // lunar nodes, Lilith
pub mod search; // ingress / aspect / station / return root finding
pub mod stars; // fixed stars: catalogue, space motion, apparent places
pub mod time; // calendar, UT/TT/TDB, Delta T
pub mod utc; // UTC with leap seconds, time zones, date checks
pub mod util; // split_deg, centisecond helpers, normalisation, names // heliacal events and visual limiting magnitude (Schaefer)

pub use error::{Error, Result};
