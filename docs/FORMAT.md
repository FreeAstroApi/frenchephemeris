# FEPK — FrenchEphemeris coefficient pack, format v1 / v2

An immutable, memory-mappable file of Chebyshev coefficient series generated
from public JPL Development Ephemerides (DE441 for the early range, DE442 for
1550+). No Swiss Ephemeris data is involved anywhere in its construction.

All integers and floats are little-endian. All offsets are absolute byte
offsets from the start of the file. Coefficient blocks are 8-byte aligned.

## Header (64 bytes, offset 0)

| off | type   | field                                    |
|-----|--------|------------------------------------------|
| 0   | [u8;8] | magic `b"FEPK\0\0\0\x01"`                |
| 8   | u32    | format_version = 1                       |
| 12  | u32    | n_series                                 |
| 16  | u64    | series_table_offset                      |
| 24  | u64    | manifest_offset                          |
| 32  | u64    | manifest_len (bytes, UTF-8 JSON)         |
| 40  | u64    | data_offset (first coefficient byte)     |
| 48  | [u8;16]| reserved, zero                           |

## Series table (`n_series` × 64 bytes)

| off | type | field                                                     |
|-----|------|-----------------------------------------------------------|
| 0   | i32  | target (NAIF id)                                          |
| 4   | i32  | center (NAIF id)                                          |
| 8   | u32  | n_coeff — Chebyshev coefficients per component            |
| 12  | u32  | n_records                                                 |
| 16  | f64  | start_jd_tdb — start of record 0                          |
| 24  | f64  | interval_days — length of each record                     |
| 32  | f64  | end_jd_tdb — = start + n_records·interval (coverage end)  |
| 40  | u64  | coeff_offset                                              |
| 48  | u32  | source_index — index into `manifest.sources`              |
| 52  | [u8;12] | reserved, zero                                         |

Series are sorted by `(target, center, start_jd_tdb)`. Series of the same
`(target, center)` are contiguous in time and do not overlap; where two
sources meet (DE441 → DE442) the boundary is a record boundary.

## Coefficients

For each series: `n_records × 3 × n_coeff` f64 values laid out as
`[record][component x,y,z][k]`. Units: km, ICRF axes, TDB.

Evaluation for time `t` (JD TDB) inside record `i`:

```
t0 = start_jd_tdb + i·interval_days
s  = 2·(t − t0)/interval_days − 1            ∈ [−1, 1]
pos_c = Σ_k c[i][c][k] · T_k(s)               km
vel_c = Σ_k c[i][c][k] · T'_k(s) · 2/interval_days   km/day
```

`t == end_jd_tdb` evaluates in the last record. Times outside every series of
the requested pair are `OutOfRange`. Where two series of a pair meet
(`end` of one == `start` of the next), the epoch itself is evaluated in the
*later* series (the one with the greatest `start_jd_tdb <= t`), i.e. the seam
epoch belongs to DE442.

Readers accept `2 <= n_coeff <= 32` (DE44x use 6..14).

## Format v2 (extended-range refit packs)

v2 keeps the magic, header and 64-byte series rows of v1 and gives meaning to
the row's reserved bytes 52–63 (zero in every v1 file, which a v2 reader checks):

| off | type | field |
|-----|------|-------|
| 52  | u32  | `n_coeff_f64` — leading coefficients per component stored as f64; `0` or `>= n_coeff` = all f64 (v1 layout) |
| 56  | u32  | `kind` — 0 Chebyshev series, 1 alias |
| 60  | i32  | `alias_target` (kind 1) |

Header `format_version` is 2 for files that use these fields; readers accept 1 and 2.

**Mixed precision** (kind 0, `0 < n_coeff_f64 = m < n_coeff`). A record is
`[x,y,z][k < m]` as f64 followed by `[x,y,z][m <= k < n]` as f32
(little-endian IEEE single), the f32 block zero-padded to a multiple of 8
bytes: record size = `24·m + 8·ceil(3(n−m)/2)` bytes. Evaluation is unchanged
(f32 values are widened to f64). The writer picks `m` per series so that
`Σ_{k>=m} max|c_k| · 2^-24` (the worst-case rounding error, |T_k| <= 1) stays
below a quarter of the series' fit tolerance.

**Alias** (kind 1). The series has `n_records = 1`, `n_coeff = 0`,
`interval_days = end − start`, and `coeff_offset` points to one f64 `scale`.
Its value at t is `scale × state(alias_target wrt center, t)`, evaluated from
the stored series of that pair, which must cover the alias series and must
not itself be an alias. Used for Earth wrt EMB = −1/EMRAT × Moon wrt EMB.

Files are written data-first: coefficients from offset 64, then the series
table, then the manifest (all offsets are absolute, so this is legal in v1 too).

### Extended-range packs (`tools/build_ext_pack.py`)

| file | span (JD TDB) | content |
|---|---|---|
| `fepk-de441-ext-past-v1.bin` | DE441 start (−3100015.5, −13200) … main pack start (per pair, 1720624.5 / 1720628.5) | refit of DE441 part 1 |
| `fepk-de441-ext-future-v1.bin` | main pack end (per pair, 2598072.5 … 2598096.5) … ≤ 8000016.5 (17191) | refit of DE441 part 2 |

Same pairs as the main pack; 399 wrt 3 is an alias of 301 wrt 3. Each pair is
one series per chunk of 65,536 days with its own `n_coeff` (adaptive);
record lengths are fixed per body. Series end/start exactly at the main
pack's per-pair coverage limits, so `Pack::open_many([main, past, future])`
chains them as contiguous series of the same pair. Fit method, tolerances
and achieved errors: docs/EXTENDED_RANGE_REPORT.md and the manifest `fit`.

## Bodies (v1)

| target | center | meaning                         |
|--------|--------|---------------------------------|
| 1..9   | 0      | planet system barycenters / SSB |
| 10     | 0      | Sun / SSB                       |
| 301    | 3      | Moon / Earth-Moon barycenter    |
| 399    | 3      | Earth / Earth-Moon barycenter   |

(Mercury 199 and Venus 299 coincide with their barycenters in DE44x; the
reader accepts 199/299 as aliases of 1/2.)

### Asteroid packs (`fepk-asteroids-v1.bin`, on-demand `ast_<n>.bin`)

Same format. Targets are NAIF small-body ids 2000000 + MPC number (20000000 +
n for n >= 1e6), centre 0 (SSB; barycentric = JPL heliocentric + the main
pack's Sun), JD 2305426.5 (1599-12-11) .. 2598072.5. A body's series may be
many short contiguous series where the fitter halved records near perihelion.
Any small-body id found in a loaded pack becomes chainable (the reader assigns
nodes when the packs are opened); a body may come from one file only.
Manifest `bodies` is keyed by MPC number (name, reasons, JPL solution id,
coverage, fit statistics, JPL reference states). docs/ASTEROIDS_REPORT.md.

### Minor-body pack (`fepk-minor-v2.bin`, tools/fit_minor.py)

Same format, separate file, loaded together with the main pack
(`Pack::open_many([main, minor])`; Python `Engine(path, extra_packs=[minor])`,
or `open_default()` which picks up `fepk-minor-v2.bin` next to the main pack).

| target | center | meaning |
|--------|--------|---------|
| 2000001, 2000002, 2000003, 2000004 | 0 | Ceres, Pallas, Juno, Vesta / SSB |
| 2002060, 2005145, 2136199 | 0 | Chiron, Pholus, Eris / SSB |
| 2000433 | 10 | Eros / Sun |

Targets are NAIF small-body ids (2000000 + MPC number). The center is chosen
per body by the fitter (Sun 10 or SSB 0, whichever needs fewer coefficients;
barycentric = JPL heliocentric + the main pack's Sun, which the reader's chain
subtracts again). A pair may
consist of several contiguous series: a grid series plus a one-record tail,
and, for barycentric series spanning 1550, a restart at the main pack's
DE441/DE442 seam. Coefficients are our own Chebyshev-Lobatto fits, not
copies (docs/MINOR_BODIES_REPORT.md). In v2, Chiron, Pholus, Eris and Eros
have further series before JD 2305426.5 (1599-12-11, JPL's first epoch),
fitted to our own integration of the JPL orbit; their `source_index` points
to a manifest source of kind "FrenchEphemeris integration from JPL#…". Series of different files may not
overlap; outside a body's series the reader returns `OutOfRange`.

Earth wrt EMB is stored even though in DE44x it is (to rounding) Moon wrt EMB
scaled by −1/EMRAT (≈ 28 % of a full pack): keeping it preserves bit-identity
with JPL; dropping it is a possible v2 size optimisation.

## Manifest (JSON)

`sources` is listed in time order and `source_index` in the series table
indexes it. The builder (tools/build_pack.py, v1.0.0) additionally writes:
`requested_jd_tdb`, `margin_days`, `pair_coverage_jd_tdb` (`"target/center"` →
`[start, end]`), `seams` (per pair: seam JD, DE441 end record, DE442 start
record, whether the grids were aligned), `series` (one object per table row,
including `first_source_record`, the index of the source SPK record copied
into record 0), per-source `constants`, `size_bytes`, `coverage_jd_tdb` and
`md5_verified` (false for the remotely range-read DE441, whose md5 is the
published one). `coverage_jd_tdb` is the intersection of all pair coverages.
A sidecar `<pack>.manifest.json` repeats the manifest plus the file size and
sha256 (these cannot live inside the file they describe).

Abridged:

```json
{
  "name": "fepk-de441-de442",
  "format_version": 1,
  "created_utc": "...",
  "builder": {"tool": "tools/build_pack.py", "version": "..."},
  "coverage_jd_tdb": [start, end],
  "sources": [
    {"name": "de441_part-1.bsp", "url": "...", "md5": "...", "used_jd_tdb": [a, b]},
    {"name": "de442.bsp", "url": "...", "md5": "...", "used_jd_tdb": [b, c]}
  ],
  "constants": {"AU_km": ..., "EMRAT": ..., "GM_km3_s2": {...}},
  "license": "JPL ephemeris data: public, see https://naif.jpl.nasa.gov/naif/rules.html"
}
```

The coefficients are re-encoded exactly from the source SPK type-2 records
(bit-identical doubles), so dataset reconstruction error versus direct JPL
evaluation is zero up to floating-point evaluation order.
