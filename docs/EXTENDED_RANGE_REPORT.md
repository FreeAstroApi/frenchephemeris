# Extended range (gap G1): −13198 … +17190

Date 2026-10-29. Branch `ws/extended-range` (merged with `ws/integration`).
Machine: Apple M4. Everything below was run; reproduce commands are in each section.

## 1. Packs

| file | span (JD TDB) | years | size |
|---|---|---|---|
| `fepk-de441-de442-v1.bin` (main, unchanged) | 1720624.5 … 2598072.5 | −1 … 2400 | 247.5 MB |
| `fepk-de441-ext-past-v1.bin` | −3099343.5 … 1720624.5 | −13198 … −1 | 144.0 MB |
| `fepk-de441-ext-future-v1.bin` | 2598096.5 … 7999696.5 | 2401 … 17190 | 160.9 MB |

The two extended packs add **305 MB**, within the 400 MB target. They are separate files.
`Pack::open_many([main, past, future])` chains each pair across the three files. The series are
contiguous at the main pack's per-pair limits, and outside the loaded span every call returns
`OutOfRange`. Python `open_default()` picks the extended packs up when they sit next to the main
pack (`FRENCHEPHEMERIS_EXT_PACKS` overrides this; `""` disables). Without them the behaviour is
the same as before.

The common span starts 672 days after DE441's first epoch (−13200 Aug 15) and ends 320 days
before its last. Outer planets use 1024-day records on a grid anchored at the main-pack seams.
Mercury and the Moon reach the DE441 limits within one record.

Build: `python tools/build_ext_pack.py --era past|future`. DE441 is streamed by HTTP Range,
one body × one 65,536-day chunk at a time, in memory only, prefetched one chunk ahead. The source
is ssd.jpl.nasa.gov `de441.bsp`; it serves the same DE441 as NAIF's two parts and was 5× faster
here. Streamed: 1.02 GB for the past pack and 1.17 GB for the future pack, about 50 and 63 min.
Nothing of the source touched the disk. The build is resumable: per-body part files are written
next to the output and deleted after assembly.

## 2. Format v2 (docs/FORMAT.md)

v2 uses the reserved bytes 52–63 of a series row; v1 files have zeros there and still load.

- **Mixed precision.** Coefficients of degree ≥ `n_coeff_f64` are stored as float32. The writer
  chooses `n_coeff_f64` so that Σ max|c_k|·2⁻²⁴ ≤ tol/4.
- **Alias series.** Earth wrt EMB is −1/EMRAT × (Moon wrt EMB). In DE441 the coefficients agree
  to 5e-13 km (checked at seven epochs across both parts). Earth is therefore not stored, which
  saves 28 %.

The reader stays allocation-free. Unit tests cover mixed records, aliases and multi-file chaining.

## 3. Refit and error budget

Each record is the Chebyshev interpolant at N Gauss–Lobatto nodes. The end points are nodes, so
neighbouring records and the main pack meet at the source value. Record length is fixed per body;
it was chosen from a size/accuracy study on DE442. N is the smallest value that passes, chosen
**per chunk**; each chunk is one series.

Held-out check, on every record and after float32 rounding: the N−1 midpoints between nodes plus
3N random epochs, about 47 million points in total. A record is accepted only if
max |r_fit − r_DE441| ≤ 0.8·tol.

Tolerances (3-D, km):
- **Planets:** 0.001″ seen from Earth at the body's closest approach, split between the body and
  the EMB (EMB 0.04 km). Venus gets 0.15, Mars 0.2, Mercury 0.3 and the Sun 0.6. The outer
  planets get 0.9, so they stay under 1 km.
- **Moon:** 0.1 km, which is 0.058″ at perigee. A 0.001″ Moon would need 1.7 m; the study
  estimated about 300 MB for the Moon alone, so the spec's 0.1 km governs.

| pair | record (d) | N | n_f64 | tol | max past / future (km) | rms (km) | max Δv (km/d) | MB past + future |
|---|---|---|---|---|---|---|---|---|
| Moon 301/3 | 32 | 28–29 | 4 | 0.10 | 0.052 / 0.080 | 0.008 | 0.14 | 60.3 + 67.2 |
| EMB 3/0 | 128 | 28 | 5–6 | 0.04 | 0.028 / 0.029 | 0.010 | 0.013 | 15.4 + 16.9 |
| Mercury 1/0 | 32 | 19 | 5 | 0.30 | 0.173 / 0.196 | 0.035 | 0.35 | 43.4 + 48.6 |
| Venus 2/0 | 128 | 16–17 | 6 | 0.15 | 0.119 / 0.119 | 0.014 | 0.029 | 10.2 + 11.3 |
| Mars 4/0 | 256 | 18–20 | 5 | 0.20 | 0.160 / 0.160 | 0.020 | 0.023 | 5.4 + 6.3 |
| Sun 10/0 | 256 | 24 | 1 | 0.60 | 0.453 / 0.464 | 0.165 | 0.081 | 5.7 + 6.4 |
| Jupiter–Pluto | 1024 | 6–17 | 2–4 | 0.90 | ≤ 0.719 | ≤ 0.20 | ≤ 0.023 | 3.5 + 4.0 |
| Earth 399/3 | alias −Moon/EMRAT | — | — | — | Moon error /81 | | | 0 |

Full per-chunk statistics are in each pack's manifest (`fit.bodies`, `series`).

Geocentric worst case, adding the body, EMB and Earth terms at closest approach: Venus 0.19 km is
0.001″ at 0.26 au, and Mars 0.36 km is 0.0013″ at 0.37 au. Measured errors (next section) are
well below these bounds.

## 4. Validation of positions

**JPL Horizons (DE441), geometric geocentric ICRF vectors.** Horizons serves only 9999 BC Mar 21
… AD 9999, so −13198 … −9999 and 10000 … 17190 rest on the in-fitter DE441 check alone. Run with
`harness/horizons_ext.py --per-era 30`, 30 epochs per body per era.

| body | −9998 … −1: max km / max ″ | 2401 … 9999: max km / max ″ |
|---|---|---|
| Sun | 0.356 / 0.00050 | 0.411 / 0.00057 |
| Moon | 0.015 / 0.0079 | 0.022 / 0.0126 |
| Mercury | 0.050 / 0.00010 | 0.120 / 0.00018 |
| Venus | 0.077 / 0.00036 | 0.040 / 0.00019 |
| Mars | 0.046 / 0.00010 | 0.071 / 0.00006 |
| Jupiter | 0.396 / 0.00013 | 0.448 / 0.00013 |
| Saturn | 0.187 / 0.00003 | 0.132 / 0.00002 |
| Uranus / Neptune / Pluto | ≤ 0.394 / ≤ 0.00002 | ≤ 0.328 / ≤ 0.00002 |
| EMB − Earth (tests the alias) | 0.0002 km | 0.0003 km |

Rust test `tests/pack_ext.rs` checks 264 Horizons vectors against the budget (worst error/budget
0.64), Earth = −Moon/EMRAT to 1e-9 km, coverage and `OutOfRange`, and seam continuity:
- **Past seam** (DE441 main → DE441 refit): below the body tolerance.
- **Future seam at 2401** (DE442 → DE441): the jump is the DE442−DE441 difference — Uranus
  23,600 km, Neptune 6,075 km, Pluto 349 km, Mars 31 km, Moon 0.05 km, inner bodies about 1.3 km.
  This is the same kind of seam as the main pack's 1550 seam.

**Swiss Ephemeris (black box).** The oracle directory has `sepl_`/`semo_` files for −600 … 2400
only. Outside that span Swiss falls back to Moshier; those calls are detected by their flags and
skipped. So **Swiss can only be compared for −600 … 0** in the extended range (0 … 2400 was
covered before). Run with `harness/extended_range_checks.py --sections positions`, 60 epochs per
body, TT:

| −600 … 0 | Sun | Moon | Mercury | Venus | Mars | Jupiter | Saturn | Uranus | Neptune | Pluto |
|---|---|---|---|---|---|---|---|---|---|---|
| max Δlon ″ | 0.001 | 0.039 | 0.028 | 0.044 | 0.023 | 0.014 | 0.029 | 0.007 | 0.019 | 0.149 |
| max Δlat ″ | 0.58 | 0.54 | 0.49 | 0.56 | 0.59 | 0.56 | 0.55 | 0.54 | 0.56 | 0.54 |

The common 0.5″ in latitude is the ecliptic-of-date definition, not position data (see §5,
obliquity).

## 5. Engine at extreme dates

All results come from `harness/extended_range_checks.py`, compared with pyswisseph 2.10.3.2.

| area | result |
|---|---|
| julday / revjul | 40,000 random dates, −13200 … 17200, both calendars: round trip exact and identical to Swiss. |
| Delta T | Same long-term parabola as Swiss (SMH 2016 −320 + 32.5u²). Ours (Swiss-compatible default) − Swiss `deltat`: +2769 s at −13000, +1228 s at −8000, +194 s at −2000, +47 s at 0, −32 s at 3000, −441 s at 10000, −849 s at 17000. The past part is the tidal-acceleration convention: we apply DE441's −25.936″/cy² to the parabola, while Swiss's bare `deltat` uses −25.80 (0.91·Δṅ·u² = 2790 s at −13000). The future part is the two extrapolations. Both are model choices; nobody knows ΔT to better than hours at ±13 kyr. |
| Precession (Vondrák 2011), nutation (IAU 2000A), ERA | vs ERFA at 32 epochs over −13200 … 17191: ltp/ltpb ≤ 3e-16, nut00a ≤ 2e-15 rad, era ≤ 1.5e-11 rad (`tests/frames_ext.rs`). TDB−TT (Fairhead) stays ≤ 1.8 ms. |
| Sidereal time | **Fixed.** Python `sidtime()` used the IAU 2006 GAST polynomial at all dates; at ±13 kyr it was wrong by 9 h (−32,196 s at −13000). It now uses the engine's documented model: IAU 2006 in 1850–2050, Swiss long-term (Simon 1994 + Vondrák) outside. Remaining difference vs Swiss: ≤ 0.2 s before 1850; +0.1 … −3.2 s after 2050, growing to 17000, which is the future ΔT difference (0.041″ of mean Sun per second of ΔT). |
| Obliquity / houses | **Fixed.** The engine (`API_FRAME_MODEL`) now uses the published Vondrák, Capitaine & Wallace (2011) ε_A series (Eq. 10, Table 3), and so does `SWISS_COMPAT`. That series is Swiss's mean obliquity: mean − Swiss ≤ 5e-11″ and true − Swiss ≤ 3e-10″ over −13000 … +17000. Placidus cusps/Asc at 47°N, with Swiss's ARMC, match to ≤ 6e-5″. Before, they were off by up to 110″, because the old value (the angle between Vondrák's poles) differed by 1.0″ at −1000, 46″ at −13000 and 66″ at 14000. The precession matrix is unchanged (Vondrák poles, as in Swiss). The long-term sidereal time now uses the same ε_A; its residual vs Swiss is a constant on each side (−0.10″ / +1.95″) instead of drifting to 23″. See TIME_FRAMES_REPORT §7.2, `tests/obliquity.rs`. The 0.5″ ecliptic-latitude difference at −600 in §4 had the same cause and should now vanish; it has not been re-run through the Python harness because no maturin build was allowed. |
| Ayanamsas | Lahiri/Fagan/Raman/KP vs Swiss ≤ 0.009″ over −13000 … 17000. `true_citra` works over the whole range; Swiss returns garbage outside its files (178° at −13000). |
| Rise/set | Sun and Moon rise/set succeed at all 30 test years. Sunrise − Swiss ≤ 0.02 s where Swiss has files, 0.15 s at 2100 (ΔT). |
| Eclipses | NASA Five Millennium Canons, −1999 … 3000, chained global searches (`harness/eclipses_canon_ext.py`). **Solar:** 11,898 in the catalogue, 1 missed (grazing partial at −1838) and 1 spurious (grazing partial at −741), 4 type flips, all before −1000. **Lunar:** 12,064, 0 missed, 5 spurious penumbral, 4 type flips (N↔P), all before year 0. Greatest-eclipse residual in TT, median / max per 500 years: 1–3 / ≤ 5 s in 500–2500, 13 s at 3000, 17 / 34 s at 0, 95 / 145 s at −1000, 272 / 358 s at −2000. Swiss (DE431) agrees with us to ≤ 8 s at −550, so the growing residual is NASA's ELP-2000/82 lunar theory vs DE441, not our search. The misses and flips are grazing eclipses that move across a limit with that 2–6 min lunar difference. |

**Fixed: Delta T bug in eclipses and occultations.** The eclipse and occultation code, and the
Python eclipse/occultation/`solcross`/`mooncross` TT wrappers, converted UT↔TT with the
library-default ΔT model instead of the engine's. With the Python default (Swiss-compatible ΔT)
that put UT eclipse times, and local circumstances, off by the model difference: 183 s at year
1000, 178 s at −2000. They now use `Engine::ut_to_tt` / `tt_to_ut`.

## 6. Speed

Measured with `cargo test -p fe-core --test pack_ext -- --ignored --nocapture`, dev profile
(opt-level 1, so only relative numbers matter). Geocentric `position`, 1e6 calls:
- **Extended span:** Moon 415 ns, Mars 528 ns, Sun 567 ns at random epochs.
- **Main span, same build:** 241 / 272 / 293 ns.

The extended path is about 1.8× slower for two reasons. Each pair holds about 80 series, which
costs a binary search. The Earth alias evaluates the Moon series a second time.

## 7. Not done / open

- Obliquity / ecliptic of date vs Swiss outside 1000–3000 (§5).
- Horizons cannot check −13198 … −9999 or 10000 … 17190 (no service there). The refit is checked
  against DE441 itself on every record.
- No release-profile timing (dev builds only, per instructions).
- The extended packs are build products (gitignored), in `data/packs/` with `.manifest.json`
  sidecars (size, sha256).
