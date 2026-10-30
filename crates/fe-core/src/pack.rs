// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Reader for FEPK coefficient packs (docs/FORMAT.md).
//!
//! CONTRACT (owned by the data-pipeline workstream):
//! - `Pack::open(path)` memory-maps a pack; `Pack::from_bytes(Vec<u8>)` for tests.
//! - `Pack::state(target, center, jd_tdb)` returns (position km, velocity km/day)
//!   of `target` relative to `center` for any pair that can be chained through
//!   the stored series (e.g. 399 wrt 0 = (3 wrt 0) + (399 wrt 3); 10 wrt 399 ...).
//! - `Pack` is `Send + Sync` and immutable.
//!
//! Chainable bodies: 0 (SSB), 1..=10, 301, 399, plus the aliases 199 → 1 and
//! 299 → 2 (Mercury/Venus coincide with their system barycenters in DE44x),
//! and the minor bodies of `MINOR_BODIES` (NAIF 2000000 + MPC number, stored
//! wrt the Sun or the SSB in a separate minor-body pack).
//!
//! Several pack files can be combined into one `Pack` (`Pack::open_many`):
//! their series tables are merged, each series remembering which file holds
//! its coefficients. A pair may come from one file only (or continue
//! contiguously in another); overlapping series are rejected.
//!
//! Extended-range packs (FEPK v2, tools/build_ext_pack.py) are refits of
//! DE441 stored with mixed f64/f32 coefficients; their Earth wrt EMB series is
//! an *alias*: −1/EMRAT × (Moon wrt EMB), evaluated from the Moon's series.
//! `open_many([main, ext_past, ext_future])` chains them with the main pack
//! (same pair, contiguous series).
//!
//! Evaluation does no heap allocation: body → stored-pair lookup is a fixed
//! table, series selection is a binary search over the (one or two) series of
//! the pair, the record index is a direct division (records are uniform within
//! a series), and the Chebyshev polynomials live in stack arrays.

use crate::math::Vec3;
use crate::{Error, Result};

mod format;
use format::{parse, Backing, Series, KIND_ALIAS};

/// Maximum Chebyshev coefficients per component accepted by the reader.
pub const MAX_COEFF: usize = 64;

/// NAIF ids (2000000 + MPC number) of the chainable minor bodies:
/// Ceres, Pallas, Juno, Vesta, Eros, Chiron, Pholus, Eris.
pub const MINOR_BODIES: [i32; 8] = [2000001, 2000002, 2000003, 2000004, 2000433, 2002060, 2005145, 2136199];

/// Dense node index for every body we know how to chain.
const N_NODES: usize = 13 + MINOR_BODIES.len(); // 0..=10, 301, 399, minor bodies

/// Any other NAIF small-body id (2000000 + number, or 20000000 + number for
/// numbers >= 1e6) stored in a pack gets a node assigned when the pack is built.
#[inline]
fn is_small_body(id: i32) -> bool {
    id >= 2_000_000
}

#[inline]
fn node_of(body: i32) -> Option<usize> {
    match body {
        0..=10 => Some(body as usize),
        199 => Some(1),
        299 => Some(2),
        301 => Some(11),
        399 => Some(12),
        2000001 => Some(13),
        2000002 => Some(14),
        2000003 => Some(15),
        2000004 => Some(16),
        2000433 => Some(17),
        2002060 => Some(18),
        2005145 => Some(19),
        2136199 => Some(20),
        _ => None,
    }
}

/// Link of one node to its parent through a stored pair.
#[derive(Clone, Copy, Debug)]
struct Link {
    pair: u16,
    parent: u32,
}

/// All series of one stored (target, center) pair, contiguous in time.
#[derive(Clone, Debug)]
struct Pair {
    target: i32,
    center: i32,
    first: usize, // index range into Pack::series
    last: usize,  // exclusive
    start: f64,
    end: f64,
}

/// Up to 4 links from a body to its chain root.
type Chain = ([u16; 4], usize, usize); // (pair indices, length, root node)

pub struct Pack {
    backings: Vec<Backing>,
    series: Vec<Series>,
    pairs: Vec<Pair>,
    /// Indexed by node: the fixed nodes, then `extra` small bodies.
    links: Vec<Option<Link>>,
    /// Sorted ids of chainable small bodies beyond `MINOR_BODIES` (node N_NODES + i).
    extra: Vec<i32>,
    manifest_ranges: Vec<(usize, usize)>,
    /// Used whole when this pack cannot answer (date outside its coverage, or a body it does
    /// not hold): a lossy full-range pack behind exact packs of a shorter span.
    fallback: Option<Box<Pack>>,
}

// Backing owns its memory (Mmap or Vec<u64>); everything else is plain data.
const _: () = {
    const fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Pack>();
};

impl std::fmt::Debug for Pack {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Pack")
            .field("files", &self.backings.len())
            .field("bytes", &self.backings.iter().map(|b| b.bytes().len()).sum::<usize>())
            .field("n_series", &self.series.len())
            .finish()
    }
}

impl Pack {
    pub fn open(path: &std::path::Path) -> Result<Pack> {
        Self::build(vec![Self::map(path)?])
    }

    /// One `Pack` over several files, e.g. the main planetary pack plus the
    /// minor-body pack. The first file's manifest is `manifest()`.
    pub fn open_many<P: AsRef<std::path::Path>>(paths: &[P]) -> Result<Pack> {
        let b = paths.iter().map(|p| Self::map(p.as_ref())).collect::<Result<Vec<_>>>()?;
        if b.is_empty() {
            return Err(Error::InvalidArgument("open_many: no pack files".into()));
        }
        Self::build(b)
    }

    fn map(path: &std::path::Path) -> Result<Backing> {
        let file = std::fs::File::open(path).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
        // SAFETY: packs are immutable files; modifying one while it is mapped
        // is unsupported (as for any mmap-based reader).
        let mmap = unsafe { memmap2::Mmap::map(&file) }
            .map_err(|e| Error::Io(format!("mmap {}: {e}", path.display())))?;
        Ok(Backing::Mmap(mmap))
    }

    pub fn from_bytes(bytes: Vec<u8>) -> Result<Pack> {
        Self::build(vec![Backing::owned(&bytes)])
    }

    /// In-memory variant of `open_many`.
    pub fn from_bytes_many(files: Vec<Vec<u8>>) -> Result<Pack> {
        Self::build(files.iter().map(|b| Backing::owned(b)).collect())
    }

    fn build(backings: Vec<Backing>) -> Result<Pack> {
        let mut series = Vec::new();
        let mut manifest_ranges = Vec::new();
        for (k, b) in backings.iter().enumerate() {
            let parsed = parse(b.bytes())?;
            manifest_ranges.push(parsed.manifest);
            series.extend(parsed.series.into_iter().map(|mut s| {
                s.backing = k;
                s
            }));
        }
        // Each table is sorted; the merge keeps (target, center, start) order.
        series.sort_by(|a, b| {
            (a.target, a.center).cmp(&(b.target, b.center)).then(a.start.total_cmp(&b.start))
        });
        // Group consecutive series of the same pair (the table is sorted).
        let mut pairs: Vec<Pair> = Vec::new();
        for (i, s) in series.iter().enumerate() {
            match pairs.last_mut() {
                Some(p) if p.target == s.target && p.center == s.center => {
                    let prev = &series[i - 1];
                    if s.start != prev.end {
                        return Err(Error::InvalidPack(format!(
                            "series {} wrt {} not contiguous: {} then {}",
                            s.target, s.center, prev.end, s.start
                        )));
                    }
                    p.last = i + 1;
                    p.end = s.end;
                }
                _ => pairs.push(Pair {
                    target: s.target,
                    center: s.center,
                    first: i,
                    last: i + 1,
                    start: s.start,
                    end: s.end,
                }),
            }
        }
        // Resolve alias series to the aliased pair (same center, plain Chebyshev series).
        for i in 0..series.len() {
            if series[i].kind != KIND_ALIAS {
                continue;
            }
            let (at, c) = (series[i].alias_target, series[i].center);
            let Some(pi) = pairs.iter().position(|p| p.target == at && p.center == c) else {
                return Err(Error::InvalidPack(format!(
                    "series {} wrt {c} aliases missing pair {at} wrt {c}",
                    series[i].target
                )));
            };
            let ap = &pairs[pi];
            if series[ap.first..ap.last].iter().any(|s| s.kind == KIND_ALIAS) {
                return Err(Error::InvalidPack("alias of an alias".into()));
            }
            if series[i].start < ap.start || series[i].end > ap.end {
                return Err(Error::InvalidPack(format!(
                    "alias series {} wrt {c} outside the coverage of {at} wrt {c}",
                    series[i].target
                )));
            }
            series[i].alias_pair = pi;
        }
        if pairs.len() > u16::MAX as usize {
            return Err(Error::InvalidPack(format!("too many stored pairs ({})", pairs.len())));
        }
        let mut extra: Vec<i32> =
            pairs.iter().map(|p| p.target).filter(|&t| node_of(t).is_none() && is_small_body(t)).collect();
        extra.sort_unstable();
        extra.dedup();
        let node = |id: i32| node_of(id).or_else(|| extra.binary_search(&id).ok().map(|i| N_NODES + i));
        let n_nodes = N_NODES + extra.len();
        let mut links: Vec<Option<Link>> = vec![None; n_nodes];
        for (pi, p) in pairs.iter().enumerate() {
            if p.target == 199 || p.target == 299 {
                continue; // aliases never own a link
            }
            let (Some(t), Some(c)) = (node(p.target), node(p.center)) else {
                continue; // pair kept in the file but not chainable
            };
            if links[t].is_some() {
                return Err(Error::InvalidPack(format!("two stored parents for body {}", p.target)));
            }
            links[t] = Some(Link { pair: pi as u16, parent: c as u32 });
        }
        for n in 0..n_nodes {
            let mut cur = n;
            let mut depth = 0;
            while let Some(l) = links[cur] {
                cur = l.parent as usize;
                depth += 1;
                if depth > 4 {
                    return Err(Error::InvalidPack("body chain cyclic or deeper than 4".into()));
                }
            }
        }
        Ok(Pack { manifest_ranges, backings, series, pairs, links, extra, fallback: None })
    }

    #[inline]
    fn node(&self, body: i32) -> Option<usize> {
        node_of(body).or_else(|| {
            if is_small_body(body) {
                self.extra.binary_search(&body).ok().map(|i| N_NODES + i)
            } else {
                None
            }
        })
    }

    /// Small bodies beyond the fixed list that this pack can chain (NAIF ids).
    pub fn extra_small_bodies(&self) -> &[i32] {
        &self.extra
    }

    #[inline]
    fn chain(&self, body: i32) -> Result<Chain> {
        let mut node = self.node(body).ok_or(Error::UnknownBody(body))?;
        let mut out = [0u16; 4];
        let mut n = 0;
        while let Some(l) = self.links[node] {
            out[n] = l.pair;
            n += 1;
            node = l.parent as usize;
        }
        Ok((out, n, node))
    }

    /// Locate the series of `pair` containing `t` and evaluate it.
    #[inline]
    fn eval_pair<const VEL: bool>(&self, pair: &Pair, body: i32, t: f64) -> Result<(Vec3, Vec3)> {
        // Written so that NaN fails the test.
        if !(t >= pair.start && t <= pair.end) {
            return Err(Error::OutOfRange { body, jd: t });
        }
        let s = self.find_series(pair, t);
        if s.kind == KIND_ALIAS {
            // Validated in build(): the aliased pair covers this series and is plain.
            let s2 = self.find_series(&self.pairs[s.alias_pair], t);
            let (p, v) = eval_series::<VEL>(s2, &self.backings[s2.backing], t);
            let k = s.scale;
            return Ok(([p[0] * k, p[1] * k, p[2] * k], [v[0] * k, v[1] * k, v[2] * k]));
        }
        Ok(eval_series::<VEL>(s, &self.backings[s.backing], t))
    }

    #[inline]
    fn find_series(&self, pair: &Pair, t: f64) -> &Series {
        let ser = &self.series[pair.first..pair.last];
        if ser.len() == 1 {
            &ser[0]
        } else {
            let k = ser.partition_point(|s| s.start <= t); // last series with start <= t
            &ser[k.max(1) - 1]
        }
    }

    /// Use `fallback` for every request this pack cannot answer (out of range or unknown body).
    /// The whole chain is then evaluated in the fallback, so one result never mixes packs.
    pub fn with_fallback(mut self, fallback: Pack) -> Pack {
        self.fallback = Some(Box::new(fallback));
        self
    }

    fn state_impl<const VEL: bool>(&self, target: i32, center: i32, jd: f64) -> Result<(Vec3, Vec3)> {
        match self.state_own::<VEL>(target, center, jd) {
            Err(e @ (Error::OutOfRange { .. } | Error::UnknownBody(_))) => match &self.fallback {
                Some(f) => f.state_impl::<VEL>(target, center, jd).map_err(|_| e),
                None => Err(e),
            },
            r => r,
        }
    }

    fn state_own<const VEL: bool>(&self, target: i32, center: i32, jd: f64) -> Result<(Vec3, Vec3)> {
        let (ct, mut nt, rt) = self.chain(target)?;
        let (cc, mut nc, rc) = self.chain(center)?;
        if rt != rc {
            // One side does not reach the common root: some link is missing.
            return Err(Error::UnknownBody(if rt != 0 { target } else { center }));
        }
        // Drop the shared tail (e.g. Moon wrt Earth skips EMB→SSB): cheaper and
        // avoids cancelling large barycentric vectors.
        while nt > 0 && nc > 0 && ct[nt - 1] == cc[nc - 1] {
            nt -= 1;
            nc -= 1;
        }
        let mut p = [0.0; 3];
        let mut v = [0.0; 3];
        for &pi in &ct[..nt] {
            let (dp, dv) = self.eval_pair::<VEL>(&self.pairs[pi as usize], target, jd)?;
            for k in 0..3 {
                p[k] += dp[k];
                v[k] += dv[k];
            }
        }
        for &pi in &cc[..nc] {
            let (dp, dv) = self.eval_pair::<VEL>(&self.pairs[pi as usize], center, jd)?;
            for k in 0..3 {
                p[k] -= dp[k];
                v[k] -= dv[k];
            }
        }
        Ok((p, v))
    }

    /// Position (km) and velocity (km/day) of `target` relative to `center`, ICRF, at JD TDB.
    #[inline]
    pub fn state(&self, target: i32, center: i32, jd_tdb: f64) -> Result<(Vec3, Vec3)> {
        self.state_impl::<true>(target, center, jd_tdb)
    }

    /// Position only (cheaper: skips derivative recurrence).
    #[inline]
    pub fn position(&self, target: i32, center: i32, jd_tdb: f64) -> Result<Vec3> {
        self.state_impl::<false>(target, center, jd_tdb).map(|(p, _)| p)
    }

    /// Batch helper: states at many times into `out` (same length as `jds`).
    /// Stops at the first error.
    pub fn states_into(&self, target: i32, center: i32, jds: &[f64], out: &mut [(Vec3, Vec3)]) -> Result<()> {
        if jds.len() != out.len() {
            return Err(Error::InvalidArgument("states_into: length mismatch".into()));
        }
        for (o, &t) in out.iter_mut().zip(jds) {
            *o = self.state(target, center, t)?;
        }
        Ok(())
    }

    /// Batch helper: positions at many times into `out`.
    pub fn positions_into(&self, target: i32, center: i32, jds: &[f64], out: &mut [Vec3]) -> Result<()> {
        if jds.len() != out.len() {
            return Err(Error::InvalidArgument("positions_into: length mismatch".into()));
        }
        for (o, &t) in out.iter_mut().zip(jds) {
            *o = self.position(target, center, t)?;
        }
        Ok(())
    }

    /// Covered TDB range for a stored pair, if present.
    pub fn coverage(&self, target: i32, center: i32) -> Option<(f64, f64)> {
        self.pairs.iter().find(|p| p.target == target && p.center == center).map(|p| (p.start, p.end))
    }

    /// TDB range over which `state(target, center, ..)` succeeds (intersection over the chain).
    pub fn chain_coverage(&self, target: i32, center: i32) -> Option<(f64, f64)> {
        let own = self.chain_coverage_own(target, center);
        match (own, self.fallback.as_ref().and_then(|f| f.chain_coverage(target, center))) {
            (Some((a, b)), Some((c, d))) => Some((a.min(c), b.max(d))),
            (o, f) => o.or(f),
        }
    }

    fn chain_coverage_own(&self, target: i32, center: i32) -> Option<(f64, f64)> {
        let (ct, nt, rt) = self.chain(target).ok()?;
        let (cc, nc, rc) = self.chain(center).ok()?;
        if rt != rc {
            return None;
        }
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        for &pi in ct[..nt].iter().chain(&cc[..nc]) {
            let p = &self.pairs[pi as usize];
            lo = lo.max(p.start);
            hi = hi.min(p.end);
        }
        Some((lo, hi))
    }

    /// Stored (target, center, start, end) pairs.
    pub fn pairs(&self) -> Vec<(i32, i32, f64, f64)> {
        self.pairs.iter().map(|p| (p.target, p.center, p.start, p.end)).collect()
    }

    /// Number of series rows in the table.
    pub fn n_series(&self) -> usize {
        self.series.len()
    }

    /// Raw manifest JSON (of the first file).
    pub fn manifest(&self) -> &str {
        self.manifest_of(0)
    }

    /// Raw manifest JSON of every file, in `open_many` order.
    pub fn manifests(&self) -> Vec<&str> {
        let mut out: Vec<&str> = (0..self.backings.len()).map(|k| self.manifest_of(k)).collect();
        if let Some(f) = &self.fallback {
            out.extend(f.manifests());
        }
        out
    }

    fn manifest_of(&self, k: usize) -> &str {
        let (a, b) = self.manifest_ranges[k];
        // validated as UTF-8 in parse()
        std::str::from_utf8(&self.backings[k].bytes()[a..b]).unwrap_or("")
    }
}

/// Evaluate one Chebyshev series at `t` (caller guarantees start <= t <= end).
#[inline(always)]
fn eval_series<const VEL: bool>(s: &Series, b: &Backing, t: f64) -> (Vec3, Vec3) {
    let n = s.n_coeff;
    let dt = t - s.start;
    let mut i = (dt * s.inv_interval) as usize; // dt >= 0: truncation == floor
    if i >= s.n_records {
        i = s.n_records - 1; // t == end
    }
    let t0 = s.start + i as f64 * s.interval;
    let u = 2.0 * (t - t0) * s.inv_interval - 1.0;
    if s.kind == format::KIND_QUANT {
        return eval_quant::<VEL>(s, b.bytes(), i, u);
    }
    let base = s.coeff_index + i * s.rec_words;
    if s.n_f64 == n {
        let rec = &b.f64s()[base..base + 3 * n];
        let (cx, rest) = rec.split_at(n);
        let (cy, cz) = rest.split_at(n);
        cheb::<VEL>(cx, cy, cz, n, u, s.inv_interval)
    } else {
        // Mixed record: f64 [3][n64] then f32 [3][n - n64] (FEPK v2).
        let n64 = s.n_f64;
        let m = n - n64;
        let hi = &b.f64s()[base..base + 3 * n64];
        let f0 = 2 * (base + 3 * n64);
        let lo = &b.f32s()[f0..f0 + 3 * m];
        cheb_mixed::<VEL>(hi, lo, n64, m, u, s.inv_interval)
    }
}

/// Decode record `i` of a quantized series (`format::KIND_QUANT`) and evaluate it.
#[inline(always)]
fn eval_quant<const VEL: bool>(s: &Series, bytes: &[u8], i: usize, u: f64) -> (Vec3, Vec3) {
    let n = s.n_coeff;
    let mut c = [[0.0f64; MAX_COEFF]; 3];
    let mut pos = s.quant_data * 8 + i * s.quant_rec_bits;
    for (comp, row) in c.iter_mut().enumerate() {
        for (k, out) in row.iter_mut().take(n).enumerate() {
            let w = s.quant_widths[comp * n + k] as u32;
            *out = read_signed(bytes, pos, w) as f64 * s.quant_q;
            pos += w as usize;
        }
    }
    cheb::<VEL>(&c[0][..n], &c[1][..n], &c[2][..n], n, u, s.inv_interval)
}

/// Signed `w`-bit integer at bit `pos` of a little-endian bit stream (w <= 64; the stream is
/// followed by 16 zero bytes, so the 16-byte window never runs past the file).
#[inline(always)]
fn read_signed(bytes: &[u8], pos: usize, w: u32) -> i64 {
    if w == 0 {
        return 0;
    }
    let o = pos >> 3;
    let x = u128::from_le_bytes(bytes[o..o + 16].try_into().unwrap()) >> (pos & 7);
    let v = (x as i128) << (128 - w);
    (v >> (128 - w)) as i64
}

/// `cheb` for a mixed record, reading the f32 tail in place (no copy).
#[inline(always)]
fn cheb_mixed<const VEL: bool>(
    hi: &[f64],
    lo: &[f32],
    n64: usize,
    m: usize,
    u: f64,
    inv_interval: f64,
) -> (Vec3, Vec3) {
    let n = n64 + m;
    let mut tk = [0.0f64; MAX_COEFF];
    tk[0] = 1.0;
    tk[1] = u;
    let u2 = 2.0 * u;
    for k in 2..n {
        tk[k] = u2 * tk[k - 1] - tk[k - 2];
    }
    let mut p = [0.0f64; 3];
    for c in 0..3 {
        let (h, l) = (&hi[c * n64..(c + 1) * n64], &lo[c * m..(c + 1) * m]);
        let mut acc = 0.0;
        for k in (0..m).rev() {
            acc += l[k] as f64 * tk[n64 + k];
        }
        for k in (0..n64).rev() {
            acc += h[k] * tk[k];
        }
        p[c] = acc;
    }
    if !VEL {
        return (p, [0.0; 3]);
    }
    let mut dk = [0.0f64; MAX_COEFF];
    dk[1] = 1.0;
    for k in 2..n {
        dk[k] = u2 * dk[k - 1] + 2.0 * tk[k - 1] - dk[k - 2];
    }
    let mut v = [0.0f64; 3];
    let sc = 2.0 * inv_interval;
    for c in 0..3 {
        let (h, l) = (&hi[c * n64..(c + 1) * n64], &lo[c * m..(c + 1) * m]);
        let mut acc = 0.0;
        for k in (0..m).rev() {
            acc += l[k] as f64 * dk[n64 + k];
        }
        for k in (1..n64).rev() {
            acc += h[k] * dk[k];
        }
        v[c] = acc * sc;
    }
    (p, v)
}

#[inline(always)]
fn cheb<const VEL: bool>(
    cx: &[f64],
    cy: &[f64],
    cz: &[f64],
    n: usize,
    u: f64,
    inv_interval: f64,
) -> (Vec3, Vec3) {
    // T_k(u), forward recurrence (|u| <= 1: stable).
    let mut tk = [0.0f64; MAX_COEFF];
    tk[0] = 1.0;
    tk[1] = u; // n >= 2 is guaranteed by parse()
    let u2 = 2.0 * u;
    for k in 2..n {
        tk[k] = u2 * tk[k - 1] - tk[k - 2];
    }
    let mut p = [0.0f64; 3];
    for k in (0..n).rev() {
        p[0] += cx[k] * tk[k];
        p[1] += cy[k] * tk[k];
        p[2] += cz[k] * tk[k];
    }
    if !VEL {
        return (p, [0.0; 3]);
    }
    // T'_k(u): T'_0 = 0, T'_1 = 1, T'_k = 2 T_{k-1} + 2u T'_{k-1} - T'_{k-2}.
    let mut dk = [0.0f64; MAX_COEFF];
    dk[1] = 1.0;
    for k in 2..n {
        dk[k] = u2 * dk[k - 1] + 2.0 * tk[k - 1] - dk[k - 2];
    }
    let mut v = [0.0f64; 3];
    for k in (1..n).rev() {
        v[0] += cx[k] * dk[k];
        v[1] += cy[k] * dk[k];
        v[2] += cz[k] * dk[k];
    }
    let sc = 2.0 * inv_interval;
    (p, [v[0] * sc, v[1] * sc, v[2] * sc])
}

#[cfg(test)]
mod tests {
    #[test]
    fn fallback_answers_outside_the_primary_pack() {
        use super::format::tests::{write_pack, TestSeries};
        // fallback: 3 wrt 0 constant x = 9 over [0, 300]
        let rec = vec![9.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let fb = write_pack(
            &[TestSeries {
                target: 3,
                center: 0,
                n_coeff: 2,
                start: 0.0,
                interval: 300.0,
                records: vec![rec],
            }],
            r#"{"name":"fallback"}"#,
        );
        let pack = Pack::from_bytes(synthetic_pack()).unwrap().with_fallback(Pack::from_bytes(fb).unwrap());
        assert_eq!(
            pack.position(3, 0, 105.0).unwrap()[0],
            Pack::from_bytes(synthetic_pack()).unwrap().position(3, 0, 105.0).unwrap()[0]
        );
        assert_eq!(pack.position(3, 0, 250.0).unwrap(), [9.0, 0.0, 0.0]); // outside the primary: fallback
        assert!(matches!(pack.position(3, 0, 400.0), Err(Error::OutOfRange { body: 3, .. })));
        assert_eq!(pack.chain_coverage(3, 0), Some((0.0, 300.0)));
        assert_eq!(pack.manifests().len(), 2);
    }

    #[test]
    fn quantized_bit_stream_round_trip() {
        // values packed little-endian at mixed widths, as tools/build_ext_pack.py writes them
        let vals: [(i64, u32); 6] =
            [(-5, 4), (123_456_789_012, 41), (0, 0), (-1, 1), (7, 4), (-(1 << 40), 42)];
        let (mut acc, mut nbits) = (0u128, 0u32);
        for (v, w) in vals {
            if w > 0 {
                acc |= ((v as u128) & ((1u128 << w) - 1)) << nbits;
                nbits += w;
            }
        }
        let mut bytes = acc.to_le_bytes().to_vec();
        bytes.extend([0u8; 16]);
        let mut pos = 0;
        for (v, w) in vals {
            assert_eq!(super::read_signed(&bytes, pos, w), v);
            pos += w as usize;
        }
    }

    use super::format::tests::synthetic_pack;
    use super::*;

    #[test]
    fn synthetic_eval_and_chaining() {
        let pack = Pack::from_bytes(synthetic_pack()).unwrap();
        // 3 wrt 0: x = 1 + 2 T1(u); 301 wrt 3: y = 5 T2(u); 399 wrt 3: z = -T0.
        // Records of 10 days from JD 100; t = 107.5 -> u = 0.5.
        let (p, v) = pack.state(3, 0, 107.5).unwrap();
        assert!((p[0] - 2.0).abs() < 1e-15 && p[1] == 0.0 && p[2] == 0.0);
        assert!((v[0] - 2.0 * 2.0 / 10.0).abs() < 1e-15);
        let (p, v) = pack.state(301, 0, 107.5).unwrap();
        // T2(0.5) = -0.5 -> y = -2.5; T2' = 4u = 2 -> vy = 5*2*(2/10) = 2
        assert!((p[1] + 2.5).abs() < 1e-15 && (p[0] - 2.0).abs() < 1e-15);
        assert!((v[1] - 2.0).abs() < 1e-15);
        // Moon wrt Earth: skips the EMB link entirely
        let p = pack.position(301, 399, 107.5).unwrap();
        assert!((p[1] + 2.5).abs() < 1e-15 && p[0] == 0.0 && (p[2] - 1.0).abs() < 1e-15);
        assert_eq!(pack.position(399, 301, 107.5).unwrap(), [-p[0], -p[1], -p[2]]);
        // 2nd series (second "source") is used after the seam at 120
        let p = pack.position(3, 0, 125.0).unwrap();
        assert!((p[0] - 7.0).abs() < 1e-15);
        // end inclusive, outside is an error
        assert!(pack.position(3, 0, 140.0).is_ok());
        assert!(matches!(pack.position(3, 0, 140.0001), Err(Error::OutOfRange { body: 3, .. })));
        assert!(matches!(pack.position(3, 0, 99.9), Err(Error::OutOfRange { .. })));
        assert!(matches!(pack.position(3, 0, f64::NAN), Err(Error::OutOfRange { .. })));
        assert!(matches!(pack.position(5, 0, 105.0), Err(Error::UnknownBody(5))));
        assert!(matches!(pack.position(0, 5, 105.0), Err(Error::UnknownBody(5))));
        assert!(matches!(pack.position(42, 0, 105.0), Err(Error::UnknownBody(42))));
        assert_eq!(pack.coverage(3, 0), Some((100.0, 140.0)));
        assert_eq!(pack.coverage(301, 3), Some((100.0, 120.0)));
        assert_eq!(pack.chain_coverage(301, 0), Some((100.0, 120.0)));
        assert_eq!(pack.position(0, 0, 105.0).unwrap(), [0.0; 3]);
        assert!(pack.manifest().contains("synthetic"));
        let jds = [101.0, 102.0];
        let mut out = [[0.0; 3]; 2];
        pack.positions_into(3, 0, &jds, &mut out).unwrap();
        assert_eq!(out[1], pack.position(3, 0, 102.0).unwrap());
    }

    #[test]
    fn multi_file_pack_chains_minor_bodies() {
        use super::format::tests::{write_pack, TestSeries};
        let minor = write_pack(
            &[TestSeries {
                target: 2002060,
                center: 3,
                n_coeff: 2,
                start: 100.0,
                interval: 20.0,
                records: vec![vec![10., 0., 0., 0., 0., 0.]],
            }],
            r#"{"name":"minor"}"#,
        );
        let pack = Pack::from_bytes_many(vec![synthetic_pack(), minor.clone()]).unwrap();
        // 2002060 wrt 0 = (2002060 wrt 3) + (3 wrt 0); at 107.5: 10 + 2
        let p = pack.position(2002060, 0, 107.5).unwrap();
        assert!((p[0] - 12.0).abs() < 1e-15);
        let p = pack.position(2002060, 399, 107.5).unwrap();
        assert!((p[0] - 10.0).abs() < 1e-15 && (p[2] - 1.0).abs() < 1e-15);
        assert!(matches!(pack.position(2002060, 0, 120.5), Err(Error::OutOfRange { body: 2002060, .. })));
        assert_eq!(pack.chain_coverage(2002060, 0), Some((100.0, 120.0)));
        assert_eq!(pack.manifests().len(), 2);
        assert!(pack.manifests()[1].contains("minor"));
        // the same pair in two files overlaps: rejected
        assert!(Pack::from_bytes_many(vec![synthetic_pack(), minor.clone(), minor]).is_err());
        // without the minor file the body is unknown, not out of range
        let main = Pack::from_bytes(synthetic_pack()).unwrap();
        assert!(matches!(main.position(2002060, 0, 107.5), Err(Error::UnknownBody(2002060))));
    }

    #[test]
    fn v2_mixed_precision_alias_and_multi_file_chaining() {
        use super::format::tests::{write_pack_v2, RowV2};
        // "past" file: 3 wrt 0 on [80, 100] as one record of 4 coeffs, x = 1 + 2 T1 + 0.5 T3,
        // stored with n_f64 = 2 (T2, T3 as f32); 301 wrt 3 likewise; 399 = -0.5 x 301.
        let rec = |x: [f64; 4]| -> Vec<u8> {
            let mut b = Vec::new();
            for c in 0..3 {
                let v = if c == 0 { x } else { [0.0; 4] };
                b.extend(v[..2].iter().flat_map(|f| f.to_le_bytes()));
            }
            for c in 0..3 {
                let v = if c == 0 { x } else { [0.0; 4] };
                b.extend(v[2..].iter().flat_map(|f| (*f as f32).to_le_bytes()));
            }
            b
        };
        let row = |t, c, x, start| RowV2 {
            target: t,
            center: c,
            n_coeff: 4,
            n_f64: 2,
            kind: 0,
            alias_target: 0,
            start,
            interval: 20.0,
            n_records: 1,
            payload: rec(x),
        };
        let past = write_pack_v2(
            &[
                row(3, 0, [1.0, 2.0, 0.0, 0.5], 80.0),
                row(301, 3, [3.0, 0.0, 0.25, 0.0], 80.0),
                RowV2 {
                    target: 399,
                    center: 3,
                    n_coeff: 0,
                    n_f64: 0,
                    kind: 1,
                    alias_target: 301,
                    start: 80.0,
                    interval: 20.0,
                    n_records: 1,
                    payload: (-0.5f64).to_le_bytes().to_vec(),
                },
            ],
            r#"{"name":"past"}"#,
        );
        let pack = Pack::from_bytes_many(vec![synthetic_pack(), past]).unwrap();
        // u = 2*(95-80)/20 - 1 = 0.5: T1 = .5, T2 = -.5, T3 = 4u^3-3u = -1
        let (p, v) = pack.state(3, 0, 95.0).unwrap();
        assert!((p[0] - (1.0 + 1.0 - 0.5)).abs() < 1e-15, "{p:?}");
        // d/du = 2 + 0.5*(12u^2-3) = 2; × 2/20
        assert!((v[0] - 0.2).abs() < 1e-15, "{v:?}");
        let p = pack.position(301, 3, 95.0).unwrap();
        assert!((p[0] - (3.0 - 0.125)).abs() < 1e-15);
        let (pe, ve) = pack.state(399, 3, 95.0).unwrap();
        assert!((pe[0] + 0.5 * p[0]).abs() < 1e-15 && ve[1] == 0.0);
        // chained across the file seam at 100: main file from 100 on
        assert_eq!(pack.coverage(3, 0), Some((80.0, 140.0)));
        assert_eq!(pack.coverage(399, 3), Some((80.0, 120.0)));
        assert!((pack.position(3, 0, 107.5).unwrap()[0] - 2.0).abs() < 1e-15);
        assert!(pack.position(301, 399, 90.0).is_ok());
        assert!(matches!(pack.position(399, 3, 79.0), Err(Error::OutOfRange { body: 399, .. })));
        // alias of a missing pair is rejected
        let bad = write_pack_v2(
            &[RowV2 {
                target: 399,
                center: 3,
                n_coeff: 0,
                n_f64: 0,
                kind: 1,
                alias_target: 301,
                start: 80.0,
                interval: 20.0,
                n_records: 1,
                payload: 1f64.to_le_bytes().to_vec(),
            }],
            "{}",
        );
        assert!(Pack::from_bytes(bad).is_err());
    }

    #[test]
    fn many_small_bodies_chain_dynamically() {
        use super::format::tests::{write_pack, TestSeries};
        let series: Vec<TestSeries> = (0..300)
            .map(|k| TestSeries {
                target: 2_000_100 + k,
                center: if k % 2 == 0 { 0 } else { 3 },
                n_coeff: 2,
                start: 100.0,
                interval: 20.0,
                records: vec![vec![k as f64, 0., 0., 0., 0., 0.]],
            })
            .collect();
        let ast = write_pack(&series, r#"{"name":"ast"}"#);
        let pack = Pack::from_bytes_many(vec![synthetic_pack(), ast]).unwrap();
        assert_eq!(pack.extra_small_bodies().len(), 300);
        for k in [0, 1, 150, 299] {
            let p = pack.position(2_000_100 + k, 0, 107.5).unwrap();
            // odd k chains through the EMB (3 wrt 0: x = 2 at 107.5)
            let want = k as f64 + if k % 2 == 0 { 0.0 } else { 2.0 };
            assert!((p[0] - want).abs() < 1e-12, "{k}: {p:?}");
        }
        assert!(matches!(pack.position(2_000_099, 0, 107.5), Err(Error::UnknownBody(2_000_099))));
        // the fixed minor-body ids still work alongside
        assert!(matches!(pack.position(2_002_060, 0, 107.5), Err(Error::UnknownBody(2_002_060))));
    }

    #[test]
    fn rejects_garbage() {
        assert!(Pack::from_bytes(vec![0u8; 10]).is_err());
        let mut b = synthetic_pack();
        b[0] = b'X';
        assert!(Pack::from_bytes(b).is_err());
        let mut b = synthetic_pack();
        b.truncate(b.len() - 8);
        assert!(Pack::from_bytes(b).is_err());
    }
}
