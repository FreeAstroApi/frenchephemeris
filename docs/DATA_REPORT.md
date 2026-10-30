# FEPK data report — `fepk-de441-de442-v1.bin`

Built 2026-10-28 by `tools/build_pack.py` v1.0.0. Every number below was
produced by the commands listed with it.

## Packs

| pack | size (bytes) | series | coverage (JD TDB, all pairs) | sources |
|---|---:|---:|---|---|
| `fepk-de441-de442-v1.bin` | 247,455,672 (247.5 MB) | 24 | 1720628.5 – 2598072.5 (−2 Oct 29 Julian – 2401 Mar 7) | DE441 part-1 + DE442 |
| `fepk-de442-1900-2100-v1.bin` (dev, `--range 1900 2100`) | 20,750,344 (20.8 MB) | 12 | 2414956.5 – 2488500.5 | DE442 only |

Requested range: JD 1720692.5 (1 Jan −1, Julian) to 2598007.5 (1 Jan 2401) with
a 64-day margin on both sides, widened per series to whole source records (so
per-pair coverage differs slightly: 4-day Moon/Earth records start at
1720628.5, the others at 1720624.5; ends range 2598072.5 – 2598096.5).
Packs live in `data/packs/` (gitignored)
with a `<pack>.manifest.json` sidecar holding size and sha256.

### Sources

| kernel | md5 | how read | used (JD TDB) |
|---|---|---|---|
| `de442.bsp` (119,771,136 B) | `446656322267e7b819a26cb08a0d8718` — matches `aa_checksums.txt` | local file | 2287184.5 – 2598096.5 |
| `de441_part-1.bsp` (1,651,119,104 B) | `7e5fcf9ecb5d08e1ab70c049baa60cd3` published, **not recomputed** | HTTP Range reads, 1 MiB chunk cache in `data/sources/de441_cache/` (183 MB: file record, comments, summaries, segment directories, and the 169 MB of records actually used) | 1720624.5 – 2287184.5 |

Instead of a full-file md5, DE441 data integrity rests on: Content-Length and
Last-Modified recorded in the cache, the DAF/type-2 structure checks in
`tools/spk.py` (`start_i + rsize·n + 4 − 1 == end_i`), every copied record's
MID/RADIUS matching the uniform grid exactly, and the independent jplephem
comparison below.

Constants in the manifest (from kernel comments): GM table (km³/s² and
au³/day²) for both kernels; DE441 prints AU = 149597870.7 km, EMRAT =
81.300568221497215, CLIGHT = 299792.458 km/s; DE442's comments print only the GM
table, so its EMRAT = GMB/GMM − 1 = 81.30056846634166 is derived (the same
derivation reproduces DE441's printed EMRAT to 1e-15), and AU is the IAU 2012
value (the two GMS columns of DE442 imply 149597870.7000 km).

### Series and record counts (full pack)

| pair | record (days) | coeffs/comp | records DE441 + DE442 = total | MB | coverage JD |
|---|---:|---:|---|---:|---|
| 1 wrt 0 Mercury bary | 8 | 14 | 70820 + 38861 = 109681 | 36.85 | 1720624.5 – 2598072.5 |
| 2 wrt 0 Venus bary | 16 | 10 | 35410 + 19431 = 54841 | 13.16 | 1720624.5 – 2598080.5 |
| 3 wrt 0 EMB | 16 | 13 | 35410 + 19431 = 54841 | 17.11 | 1720624.5 – 2598080.5 |
| 4 wrt 0 Mars bary | 32 | 11 | 17705 + 9716 = 27421 | 7.24 | 1720624.5 – 2598096.5 |
| 5 wrt 0 Jupiter bary | 32 | 8 | 17705 + 9716 = 27421 | 5.26 | 1720624.5 – 2598096.5 |
| 6 wrt 0 Saturn bary | 32 | 7 | 17705 + 9716 = 27421 | 4.61 | 1720624.5 – 2598096.5 |
| 7 wrt 0 Uranus bary | 32 | 6 | 17705 + 9716 = 27421 | 3.95 | 1720624.5 – 2598096.5 |
| 8 wrt 0 Neptune bary | 32 | 6 | 17705 + 9716 = 27421 | 3.95 | 1720624.5 – 2598096.5 |
| 9 wrt 0 Pluto bary | 32 | 6 | 17705 + 9716 = 27421 | 3.95 | 1720624.5 – 2598096.5 |
| 10 wrt 0 Sun | 16 | 11 | 35410 + 19431 = 54841 | 14.48 | 1720624.5 – 2598080.5 |
| 301 wrt 3 Moon | 4 | 13 | 141639 + 77722 = 219361 | 68.44 | 1720628.5 – 2598072.5 |
| 399 wrt 3 Earth | 4 | 13 | 141639 + 77722 = 219361 | 68.44 | 1720628.5 – 2598072.5 |

Total 877,452 records.

## The seam

DE441 part-1 starts every segment at JD −3100015.5 and DE442 at 2287184.5; the
difference, 5,387,200 days, is a multiple of 32, and both kernels use the same
record length per body (4/8/16/32 days). So the grids **align exactly** and
the seam is JD 2287184.5 (1549-12-31 0h TDB, as labelled in the DE442 header) for all 12 pairs: DE441 records up to
the one ending there, then DE442 from its record 0. (`seam_epoch()` in the
builder searches for a common boundary generally and fails loudly if none
exists; it was not needed.) The seam epoch itself evaluates in DE442.

`python harness/seam_report.py` — geometric geocentric directions (body − Earth,
both from SSB, no light time) at JD 2287184.5:

| body | angle DE441 vs DE442 (arcsec) | geocentric distance change (km) |
|---|---:|---:|
| Sun | 0.000062 | −0.001 |
| Moon | 0.381108 | −0.011 |
| Mercury (bary) | 0.000217 | +0.058 |
| Venus (bary) | 0.000065 | −0.138 |
| Mars (bary) | 0.007546 | +1.752 |
| Jupiter (bary) | 0.004885 | +1.497 |
| Saturn (bary) | 0.000124 | +1.089 |
| Uranus (bary) | 1.846917 | −2236.886 |
| Neptune (bary) | 0.291460 | +254.751 |
| Pluto (bary) | 0.008252 | −20.667 |

Jump of each stored series at the seam (DE441 last record at u = +1 vs DE442
first record at u = −1, read from the pack):

| series | position jump (km) | velocity jump (km/day) |
|---|---:|---:|
| 1 wrt 0 Mercury bary | 1.494 | 8.6e-03 |
| 2 wrt 0 Venus bary | 1.306 | 3.8e-03 |
| 3 wrt 0 EMB | 1.453 | 1.0e-03 |
| 4 wrt 0 Mars bary | 13.851 | 1.2e-01 |
| 5 wrt 0 Jupiter bary | 16.460 | 2.6e-02 |
| 6 wrt 0 Saturn bary | 2.499 | 1.2e-03 |
| 7 wrt 0 Uranus bary | 24501.348 | 5.37 |
| 8 wrt 0 Neptune bary | 6218.035 | 6.4e-01 |
| 9 wrt 0 Pluto bary | 245.649 | 7.3e-03 |
| 10 wrt 0 Sun | 1.451 | 2.7e-04 |
| 301 wrt 3 Moon | 0.740 | 1.6e-01 |
| 399 wrt 3 Earth | 0.009 | 2.0e-03 |

Reading: below 0.01″ everywhere except the Moon (0.38″, the 0.74 km lunar
orbit difference seen from 380,000 km), Neptune (0.29″) and Uranus (1.85″).
The Uranus difference is expected: DE442 exists to refit the Uranus barycenter
(Uranus occultation data), and in 1550 both kernels are far outside the
observation span. A step this size (1.8″ ≈ 0.0005°) is well under the arcminute
precision astrology output uses, but a search that crosses 1549-12-31 (for
example a Uranus station or an aspect) may see a 1.8″ step there. Nothing in
the pack smooths it: coefficients stay JPL-exact on both sides.

## Reconstruction accuracy (pack vs direct JPL)

The independent reference is jplephem 2.24 (its own DAF reader and Chebyshev
code); DE441 is read through a *sparse* local copy of the range cache
(`data/sources/de441_part-1.sparse.bsp`: 1.65 GB apparent, 183 MB on disk),
so jplephem opens it as an ordinary kernel. jplephem is called with the epoch
split into whole and fractional days: with a single JD float its
seconds-offset arithmetic loses ~2e-6 s, which shows up as 1e-4 km for
Mercury and is a jplephem effect, not a pack error.

`python harness/check_pack_vs_jpl.py --n 50000 --rust` (full pack):

* **Bit identity:** all 877,452 records in 24 series are bit-for-bit equal to
  the source SPK coefficient words read by jplephem's DAF reader.
* **Evaluation:** 50,000 random epochs per pair over the full coverage (about
  32k DE441 and 18k DE442 each, plus the seam and its neighbouring float), for
  the 12 stored pairs and 8 chained pairs (10/399, 301/399, 4/399, 5/399,
  1/399, 399/0, 301/0, 9/301), through the Rust reader:
  * max |Δr| = 9.8e-7 km (Pluto barycenter, |r| ≈ 7e9 km: 1 ulp); stored-pair
    relative error ≤ 8e-16 everywhere;
  * Moon wrt EMB 8.3e-11 km, Moon geocentric 1.3e-10 km, Sun geocentric
    3.1e-8 km, Mars geocentric 6.8e-8 km;
  * max |Δv| = 2.1e-9 km/day.

  The numpy evaluator in the script agrees to the same level (max relative
  error 8e-16).
* `cargo test -p fe-core --test pack`: 3,952 jplephem fixtures (12 stored + 7
  chained pairs, all eras, record boundaries and seam ±1e-6 d), max |Δr|
  6.7e-7 km, max relative 2.2e-16, max |Δv| 1.1e-9 km/day; same on the dev
  pack (1,577 fixtures in range). Tolerance: max(1e-6 km, 1e-15·|r|) and
  max(1e-6 km/day, 1e-13·|v|).

## Evaluation speed

`cargo test --release -p fe-core --test pack -- --ignored --nocapture`
(Apple M4, full pack mmapped and warm, 1e6 calls each):

| call | random epochs over 2400 yr | sequential (0.01 d steps) |
|---|---:|---:|
| `position(301, 3)` Moon wrt EMB (1 series) | 102 ns | 35 ns |
| `position(301, 399)` geocentric Moon (2 series) | 206 ns | 64 ns |
| `position(4, 0)` Mars bary (1 series) | 32 ns | 29 ns |
| `position(4, 399)` geocentric Mars (3 series) | 231 ns | 93 ns |
| `state(301, 399)` | 256 ns | 95 ns |
| `state(4, 399)` | 314 ns | 133 ns |

Random-epoch times are dominated by cache misses on the 247 MB file (each
series lookup touches a different page); with time locality, typical of
chart, transit or search workloads, a geocentric Moon costs about 60 ns.
`state`/`position` do not allocate.

## Reproduce

```sh
PY=.venv/bin/python
$PY tools/spk.py data/sources/de442.bsp                   # list segments
$PY tools/build_pack.py --range 1900 2100                   # dev pack, ~1 s
$PY tools/build_pack.py                                     # full pack (fetches ~183 MB of DE441 once)
$PY crates/fe-core/tests/fixtures/pack/gen_pack_fixtures.py # fixtures
cargo test -p fe-core --test pack
$PY harness/check_pack_vs_jpl.py --n 50000 --rust
$PY harness/seam_report.py
```
