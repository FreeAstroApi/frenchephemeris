// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Shared test helpers: locating the private Swiss-oracle corpus.
//!
//! Fixtures produced by calling the Swiss Ephemeris (pyswisseph) black box are not part of
//! the public tree. They live in a private corpus that mirrors the repository layout:
//! `$FE_ORACLE_DIR/crates/fe-core/tests/fixtures/<area>/<file>.json`, with `FE_ORACLE_DIR`
//! defaulting to `<repo>/private/oracle` (gitignored). Tests that compare against the oracle
//! call `require_oracle!("<area>/<file>.json")` first and are skipped, with a message, when
//! the corpus is absent, so the public test suite passes without it.
//! See tools/privatize_oracle_fixtures.py.
#![allow(dead_code)]

use std::path::PathBuf;

/// Root of the private oracle corpus (`FE_ORACLE_DIR`, default `<repo>/private/oracle`).
pub fn oracle_dir() -> PathBuf {
    match std::env::var_os("FE_ORACLE_DIR") {
        Some(d) if !d.is_empty() => PathBuf::from(d),
        _ => PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../private/oracle"),
    }
}

/// Path of an oracle fixture, `rel` relative to `crates/fe-core/tests/fixtures` (may not exist).
pub fn oracle_path(rel: &str) -> PathBuf {
    oracle_dir().join("crates/fe-core/tests/fixtures").join(rel)
}

/// Whether the oracle fixture `rel` is present.
pub fn has_oracle(rel: &str) -> bool {
    oracle_path(rel).is_file()
}

/// Read an oracle fixture (call `require_oracle!` first).
pub fn read_oracle(rel: &str) -> String {
    let p = oracle_path(rel);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("oracle fixture {}: {e}", p.display()))
}

/// Parse an oracle fixture as JSON, or `None` when the private corpus is absent.
pub fn oracle_json(rel: &str) -> Option<serde_json::Value> {
    has_oracle(rel).then(|| serde_json::from_str(&read_oracle(rel)).expect("oracle fixture JSON"))
}

/// Return early (test passes as skipped) when the oracle fixture `rel` is absent.
#[macro_export]
macro_rules! require_oracle {
    ($rel:expr) => {
        if !common::has_oracle($rel) {
            eprintln!(
                "SKIP {}: Swiss-oracle fixture {} not found in the private corpus ({}); set FE_ORACLE_DIR",
                module_path!(),
                $rel,
                common::oracle_dir().display()
            );
            return;
        }
    };
}
