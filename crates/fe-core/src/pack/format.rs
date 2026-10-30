// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! FEPK v1/v2 binary layout: header + series table parsing and validation.
//!
//! v2 (extended-range refit packs) uses the v1 reserved bytes of a series row:
//! `n_coeff_f64` (mixed f64/f32 coefficient storage) and `kind`/`alias_target`
//! (a series defined as a fixed multiple of another pair). See docs/FORMAT.md.

use crate::{Error, Result};

pub(super) const MAGIC: &[u8; 8] = b"FEPK\0\0\0\x01";
pub(super) const FORMAT_VERSION: u32 = 2; // highest version read; 1 is also accepted
pub(super) const HEADER_LEN: usize = 64;
pub(super) const SERIES_LEN: usize = 64;

/// Memory holding the pack. Both variants guarantee 8-byte alignment of the
/// start of the file, so 8-aligned coefficient offsets can be viewed as `&[f64]`.
pub(super) enum Backing {
    Mmap(memmap2::Mmap),
    Owned { words: Vec<u64>, len: usize },
}

impl Backing {
    pub(super) fn owned(bytes: &[u8]) -> Backing {
        let mut words = vec![0u64; bytes.len().div_ceil(8)];
        // SAFETY: `words` has at least `bytes.len()` bytes of storage.
        unsafe {
            std::ptr::copy_nonoverlapping(bytes.as_ptr(), words.as_mut_ptr() as *mut u8, bytes.len());
        }
        Backing::Owned { words, len: bytes.len() }
    }

    #[inline]
    pub(super) fn bytes(&self) -> &[u8] {
        match self {
            Backing::Mmap(m) => m,
            // SAFETY: u64 storage reinterpreted as bytes, within `len`.
            Backing::Owned { words, len } => unsafe {
                std::slice::from_raw_parts(words.as_ptr() as *const u8, *len)
            },
        }
    }

    /// Whole file viewed as f64 words (index = byte offset / 8).
    #[inline]
    pub(super) fn f64s(&self) -> &[f64] {
        let b = self.bytes();
        // SAFETY: the base pointer is 8-aligned (page-aligned mmap or Vec<u64>,
        // checked in parse()); every bit pattern is a valid f64; the host is
        // little-endian (checked in parse()), matching the file encoding.
        unsafe { std::slice::from_raw_parts(b.as_ptr() as *const f64, b.len() / 8) }
    }

    /// Whole file viewed as f32 words (index = byte offset / 4).
    #[inline]
    pub(super) fn f32s(&self) -> &[f32] {
        let b = self.bytes();
        // SAFETY: as for f64s (alignment 8 >= 4).
        unsafe { std::slice::from_raw_parts(b.as_ptr() as *const f32, b.len() / 4) }
    }
}

/// Series kinds (v2).
pub(super) const KIND_CHEBYSHEV: u32 = 0;
pub(super) const KIND_ALIAS: u32 = 1;
/// Quantized Chebyshev series (compact packs). The coefficient block starts with a header:
/// q (f64, km per integer step), rec_bits (u32), n (u32), then 3·n bit widths (u8, padded to
/// 8 bytes); then the records, each 3·n signed integers in that order (two's complement,
/// little-endian bit stream, fixed rec_bits per record), followed by 16 zero bytes.
pub(super) const KIND_QUANT: u32 = 2;

/// One series row, with derived values pre-computed for evaluation.
#[derive(Clone, Debug)]
pub(super) struct Series {
    pub target: i32,
    pub center: i32,
    pub n_coeff: usize,
    pub n_records: usize,
    pub start: f64,
    pub interval: f64,
    pub inv_interval: f64,
    pub end: f64,
    /// Index of the first coefficient in `Backing::f64s()`.
    pub coeff_index: usize,
    /// Leading coefficients per component stored as f64 (== n_coeff: all f64).
    pub n_f64: usize,
    /// Record size in 8-byte words.
    pub rec_words: usize,
    /// KIND_CHEBYSHEV or KIND_ALIAS.
    pub kind: u32,
    /// KIND_ALIAS: value = scale × state(alias_target wrt center).
    pub alias_target: i32,
    pub scale: f64,
    /// KIND_ALIAS: index of the aliased pair, resolved by `Pack::build`.
    pub alias_pair: usize,
    /// Which file of a multi-pack `Pack` holds the coefficients (0 = first).
    pub backing: usize,
    /// KIND_QUANT: integer step (km), bits per record, bit width per [component][degree],
    /// byte offset of the first record.
    pub quant_q: f64,
    pub quant_rec_bits: usize,
    pub quant_widths: Vec<u8>,
    pub quant_data: usize,
    #[allow(dead_code)]
    pub source_index: u32,
}

pub(super) struct Parsed {
    pub series: Vec<Series>,
    pub manifest: (usize, usize),
}

fn bad(msg: impl Into<String>) -> Error {
    Error::InvalidPack(msg.into())
}

#[inline]
fn u32_at(b: &[u8], o: usize) -> u32 {
    u32::from_le_bytes(b[o..o + 4].try_into().unwrap())
}
#[inline]
fn i32_at(b: &[u8], o: usize) -> i32 {
    i32::from_le_bytes(b[o..o + 4].try_into().unwrap())
}
#[inline]
fn u64_at(b: &[u8], o: usize) -> u64 {
    u64::from_le_bytes(b[o..o + 8].try_into().unwrap())
}
#[inline]
fn f64_at(b: &[u8], o: usize) -> f64 {
    f64::from_le_bytes(b[o..o + 8].try_into().unwrap())
}

fn checked_range(off: u64, len: u64, total: usize, what: &str) -> Result<(usize, usize)> {
    let end = off.checked_add(len).ok_or_else(|| bad(format!("{what}: overflow")))?;
    if end > total as u64 {
        return Err(bad(format!("{what}: [{off}, {end}) beyond file size {total}")));
    }
    Ok((off as usize, end as usize))
}

pub(super) fn parse(b: &[u8]) -> Result<Parsed> {
    if !cfg!(target_endian = "little") {
        return Err(bad("big-endian hosts are not supported"));
    }
    if !(b.as_ptr() as usize).is_multiple_of(8) {
        return Err(bad("pack memory not 8-byte aligned"));
    }
    if b.len() < HEADER_LEN {
        return Err(bad("file shorter than header"));
    }
    if &b[0..8] != MAGIC {
        return Err(bad("bad magic"));
    }
    let version = u32_at(b, 8);
    if version != 1 && version != FORMAT_VERSION {
        return Err(bad(format!("unsupported format_version {version}")));
    }
    let n_series = u32_at(b, 12) as u64;
    let table_off = u64_at(b, 16);
    let man_off = u64_at(b, 24);
    let man_len = u64_at(b, 32);
    let data_off = u64_at(b, 40);
    let (t0, _) = checked_range(table_off, n_series * SERIES_LEN as u64, b.len(), "series table")?;
    let manifest = checked_range(man_off, man_len, b.len(), "manifest")?;
    std::str::from_utf8(&b[manifest.0..manifest.1]).map_err(|_| bad("manifest is not UTF-8"))?;
    if !data_off.is_multiple_of(8) || data_off > b.len() as u64 {
        return Err(bad("bad data_offset"));
    }

    let mut series = Vec::with_capacity(n_series as usize);
    for k in 0..n_series as usize {
        let o = t0 + k * SERIES_LEN;
        let target = i32_at(b, o);
        let center = i32_at(b, o + 4);
        let n_coeff = u32_at(b, o + 8) as usize;
        let n_records = u32_at(b, o + 12) as usize;
        let start = f64_at(b, o + 16);
        let interval = f64_at(b, o + 24);
        let end = f64_at(b, o + 32);
        let coeff_off = u64_at(b, o + 40);
        let source_index = u32_at(b, o + 48);
        let n_f64_raw = u32_at(b, o + 52) as usize;
        let kind = u32_at(b, o + 56);
        let alias_target = i32_at(b, o + 60);
        let id = format!("series {k} ({target} wrt {center})");
        if version == 1 && (n_f64_raw != 0 || kind != 0 || alias_target != 0) {
            return Err(bad(format!("{id}: v1 reserved bytes not zero")));
        }
        if kind == KIND_ALIAS {
            if n_records != 1 || !(start.is_finite() && end.is_finite() && end > start) {
                return Err(bad(format!("{id}: alias series must be one record with start < end")));
            }
            if alias_target == target || !coeff_off.is_multiple_of(8) || coeff_off < data_off {
                return Err(bad(format!("{id}: bad alias")));
            }
            checked_range(coeff_off, 8, b.len(), &id)?;
            let scale = f64_at(b, coeff_off as usize);
            if !scale.is_finite() {
                return Err(bad(format!("{id}: alias scale not finite")));
            }
            series.push(Series {
                target,
                center,
                n_coeff: 0,
                n_records: 1,
                start,
                interval: end - start,
                inv_interval: 1.0 / (end - start),
                end,
                coeff_index: 0,
                n_f64: 0,
                rec_words: 0,
                kind,
                alias_target,
                scale,
                alias_pair: usize::MAX,
                backing: 0,
                quant_q: 0.0,
                quant_rec_bits: 0,
                quant_widths: Vec::new(),
                quant_data: 0,
                source_index,
            });
            continue;
        }
        if kind != KIND_CHEBYSHEV && kind != KIND_QUANT {
            return Err(bad(format!("{id}: unknown series kind {kind}")));
        }
        if !(2..=super::MAX_COEFF).contains(&n_coeff) {
            return Err(bad(format!("{id}: n_coeff {n_coeff} outside 2..={}", super::MAX_COEFF)));
        }
        if n_records == 0 {
            return Err(bad(format!("{id}: no records")));
        }
        if !(start.is_finite() && interval.is_finite() && interval > 0.0 && end.is_finite()) {
            return Err(bad(format!("{id}: bad time grid")));
        }
        let expect_end = start + n_records as f64 * interval;
        if (expect_end - end).abs() > 1e-9 * interval.max(1.0) {
            return Err(bad(format!("{id}: end {end} != start + n*interval = {expect_end}")));
        }
        if !coeff_off.is_multiple_of(8) || coeff_off < data_off {
            return Err(bad(format!("{id}: misaligned or misplaced coefficients")));
        }
        if kind == KIND_QUANT {
            let co = coeff_off as usize;
            checked_range(coeff_off, 16, b.len(), &id)?;
            let q = f64_at(b, co);
            let rec_bits = u32_at(b, co + 8) as usize;
            if u32_at(b, co + 12) as usize != n_coeff || !(q.is_finite() && q > 0.0) {
                return Err(bad(format!("{id}: bad quantized header")));
            }
            let wlen = (3 * n_coeff).div_ceil(8) * 8;
            let (w0, w1) = checked_range(coeff_off + 16, 3 * n_coeff as u64, b.len(), &id)?;
            let widths = b[w0..w1].to_vec();
            if widths.iter().any(|&w| w > 64) || widths.iter().map(|&w| w as usize).sum::<usize>() != rec_bits
            {
                return Err(bad(format!("{id}: bad quantized bit widths")));
            }
            let data = co + 16 + wlen;
            let nbytes = (n_records * rec_bits).div_ceil(8) + 16;
            checked_range(data as u64, nbytes as u64, b.len(), &id)?;
            series.push(Series {
                target,
                center,
                n_coeff,
                n_records,
                start,
                interval,
                inv_interval: 1.0 / interval,
                end,
                coeff_index: 0,
                n_f64: n_coeff,
                rec_words: 0,
                kind,
                alias_target: 0,
                scale: 1.0,
                alias_pair: usize::MAX,
                backing: 0,
                quant_q: q,
                quant_rec_bits: rec_bits,
                quant_widths: widths,
                quant_data: data,
                source_index,
            });
            continue;
        }
        let n_f64 = if n_f64_raw == 0 || n_f64_raw >= n_coeff { n_coeff } else { n_f64_raw };
        // f64 block [3][n_f64], then f32 block [3][n_coeff - n_f64] padded to 8 bytes
        let rec_words = 3 * n_f64 + (3 * (n_coeff - n_f64)).div_ceil(2);
        let nbytes = (n_records as u64) * (rec_words as u64) * 8;
        checked_range(coeff_off, nbytes, b.len(), &id)?;
        series.push(Series {
            target,
            center,
            n_coeff,
            n_records,
            start,
            interval,
            inv_interval: 1.0 / interval,
            end,
            coeff_index: (coeff_off / 8) as usize,
            n_f64,
            rec_words,
            kind,
            alias_target: 0,
            scale: 1.0,
            alias_pair: usize::MAX,
            backing: 0,
            quant_q: 0.0,
            quant_rec_bits: 0,
            quant_widths: Vec::new(),
            quant_data: 0,
            source_index,
        });
    }
    for w in series.windows(2) {
        let (a, c) = (&w[0], &w[1]);
        if (a.target, a.center) > (c.target, c.center)
            || ((a.target, a.center) == (c.target, c.center) && a.start >= c.start)
        {
            return Err(bad("series table not sorted by (target, center, start)"));
        }
    }
    Ok(Parsed { series, manifest })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Series description for the in-memory test writer.
    pub(crate) struct TestSeries {
        pub target: i32,
        pub center: i32,
        pub n_coeff: u32,
        pub start: f64,
        pub interval: f64,
        /// records x [x[n], y[n], z[n]]
        pub records: Vec<Vec<f64>>,
    }

    /// Minimal FEPK writer (mirrors tools/build_pack.py) used by unit tests.
    pub(crate) fn write_pack(series: &[TestSeries], manifest: &str) -> Vec<u8> {
        let table_off = HEADER_LEN;
        let man_off = table_off + series.len() * SERIES_LEN;
        let mut data_off = man_off + manifest.len();
        data_off = data_off.div_ceil(8) * 8;
        let mut out = vec![0u8; data_off];
        out[0..8].copy_from_slice(MAGIC);
        out[8..12].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
        out[12..16].copy_from_slice(&(series.len() as u32).to_le_bytes());
        out[16..24].copy_from_slice(&(table_off as u64).to_le_bytes());
        out[24..32].copy_from_slice(&(man_off as u64).to_le_bytes());
        out[32..40].copy_from_slice(&(manifest.len() as u64).to_le_bytes());
        out[40..48].copy_from_slice(&(data_off as u64).to_le_bytes());
        out[man_off..man_off + manifest.len()].copy_from_slice(manifest.as_bytes());
        let mut coeff_off = data_off;
        for (k, s) in series.iter().enumerate() {
            let o = table_off + k * SERIES_LEN;
            let n = s.records.len();
            out[o..o + 4].copy_from_slice(&s.target.to_le_bytes());
            out[o + 4..o + 8].copy_from_slice(&s.center.to_le_bytes());
            out[o + 8..o + 12].copy_from_slice(&s.n_coeff.to_le_bytes());
            out[o + 12..o + 16].copy_from_slice(&(n as u32).to_le_bytes());
            out[o + 16..o + 24].copy_from_slice(&s.start.to_le_bytes());
            out[o + 24..o + 32].copy_from_slice(&s.interval.to_le_bytes());
            out[o + 32..o + 40].copy_from_slice(&(s.start + n as f64 * s.interval).to_le_bytes());
            out[o + 40..o + 48].copy_from_slice(&(coeff_off as u64).to_le_bytes());
            coeff_off += n * 3 * s.n_coeff as usize * 8;
        }
        for s in series {
            for r in &s.records {
                assert_eq!(r.len(), 3 * s.n_coeff as usize);
                for v in r {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
        out
    }

    pub(crate) fn synthetic_pack() -> Vec<u8> {
        let rec = |x: &[f64], y: &[f64], z: &[f64]| [x, y, z].concat();
        write_pack(
            &[
                TestSeries {
                    target: 3,
                    center: 0,
                    n_coeff: 2,
                    start: 100.0,
                    interval: 10.0,
                    records: vec![rec(&[1., 2.], &[0., 0.], &[0., 0.]); 2],
                },
                TestSeries {
                    target: 3,
                    center: 0,
                    n_coeff: 3,
                    start: 120.0,
                    interval: 20.0,
                    records: vec![rec(&[7., 0., 0.], &[0., 0., 0.], &[0., 0., 0.])],
                },
                TestSeries {
                    target: 301,
                    center: 3,
                    n_coeff: 3,
                    start: 100.0,
                    interval: 10.0,
                    records: vec![rec(&[0., 0., 0.], &[0., 0., 5.], &[0., 0., 0.]); 2],
                },
                TestSeries {
                    target: 399,
                    center: 3,
                    n_coeff: 2,
                    start: 100.0,
                    interval: 10.0,
                    records: vec![rec(&[0., 0.], &[0., 0.], &[-1., 0.]); 2],
                },
            ],
            r#"{"name":"synthetic"}"#,
        )
    }

    /// One v2 row for `write_pack_v2`: payload is the raw coefficient bytes.
    pub(crate) struct RowV2 {
        pub target: i32,
        pub center: i32,
        pub n_coeff: u32,
        pub n_f64: u32,
        pub kind: u32,
        pub alias_target: i32,
        pub start: f64,
        pub interval: f64,
        pub n_records: u32,
        pub payload: Vec<u8>,
    }

    /// Minimal FEPK v2 writer (mirrors tools/build_ext_pack.py: data first).
    pub(crate) fn write_pack_v2(rows: &[RowV2], manifest: &str) -> Vec<u8> {
        let mut out = vec![0u8; HEADER_LEN];
        let mut offs = Vec::new();
        for r in rows {
            offs.push(out.len() as u64);
            out.extend_from_slice(&r.payload);
            while !out.len().is_multiple_of(8) {
                out.push(0);
            }
        }
        let table_off = out.len();
        for (r, off) in rows.iter().zip(&offs) {
            out.extend_from_slice(&r.target.to_le_bytes());
            out.extend_from_slice(&r.center.to_le_bytes());
            out.extend_from_slice(&r.n_coeff.to_le_bytes());
            out.extend_from_slice(&r.n_records.to_le_bytes());
            out.extend_from_slice(&r.start.to_le_bytes());
            out.extend_from_slice(&r.interval.to_le_bytes());
            out.extend_from_slice(&(r.start + r.n_records as f64 * r.interval).to_le_bytes());
            out.extend_from_slice(&off.to_le_bytes());
            out.extend_from_slice(&0u32.to_le_bytes());
            out.extend_from_slice(&r.n_f64.to_le_bytes());
            out.extend_from_slice(&r.kind.to_le_bytes());
            out.extend_from_slice(&r.alias_target.to_le_bytes());
        }
        let man_off = out.len();
        out.extend_from_slice(manifest.as_bytes());
        out[0..8].copy_from_slice(MAGIC);
        out[8..12].copy_from_slice(&2u32.to_le_bytes());
        out[12..16].copy_from_slice(&(rows.len() as u32).to_le_bytes());
        out[16..24].copy_from_slice(&(table_off as u64).to_le_bytes());
        out[24..32].copy_from_slice(&(man_off as u64).to_le_bytes());
        out[32..40].copy_from_slice(&(manifest.len() as u64).to_le_bytes());
        out[40..48].copy_from_slice(&(HEADER_LEN as u64).to_le_bytes());
        out
    }
}
