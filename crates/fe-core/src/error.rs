// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum Error {
    #[error("date out of range: jd {jd} not covered for body {body}")]
    OutOfRange { body: i32, jd: f64 },
    #[error("unknown body {0}")]
    UnknownBody(i32),
    #[error("invalid pack: {0}")]
    InvalidPack(String),
    #[error("io: {0}")]
    Io(String),
    /// House system undefined at this latitude (e.g. Placidus/Koch inside polar circle).
    #[error("house system {system} undefined at latitude {lat_deg}")]
    PolarHouses { system: char, lat_deg: f64 },
    #[error("invalid argument: {0}")]
    InvalidArgument(String),
}

pub type Result<T> = std::result::Result<T, Error>;
