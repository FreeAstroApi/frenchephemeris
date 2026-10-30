# FrenchEphemeris — status, accuracy and speed report

Date: 2026-10-28. Branch `ws/integration`. Machine: Apple M4 (macOS arm64).
Workstream details: DATA_REPORT.md, TIME_FRAMES_REPORT.md, GEOMETRY_REPORT.md,
STARS_REPORT.md, EVENTS_REPORT.md, ECLIPSES_REPORT.md, MINOR_BODIES_REPORT.md.

Reproduce (Swiss comparison reports are written to the private oracle corpus,
`private/oracle/harness/…`; see tools/privatize_oracle_fixtures.py):

```bash
.venv/bin/python harness/compare_swiss.py --n 600 --json harness/swiss_report.json
.venv/bin/python harness/compare_swiss.py --n 600 --swiss-delta-t --json harness/swiss_report_compat.json
.venv/bin/python harness/horizons.py --per-body 25 --json harness/horizons_report.json
.venv/bin/python bench/py_vs_swiss.py
cargo run --release -p fe-core --example bench_calc
cargo run --release -p fe-core --example search_demo
```

## 1. What exists

| Capability | Status |
|---|---|
| JPL data pack (DE441 0–1550 + DE442 1550–2400, 247 MB, bit-identical coefficients) | done |
| Extended range −13198 … +17190: optional DE441 refit packs (305 MB, FEPK v2, ≤ 0.001″ planets, Moon ≤ 0.1 km), chained by `open_many` / `open_default` (EXTENDED_RANGE_REPORT.md) | done (2026-10-29) |
| Time: julday/revjul, Delta T (SMH 2020 + USNO observed; Swiss-compatible option), TDB | done |
| Frames: Vondrák 2011 precession, IAU 2000A / 77-term nutation, GMST/GAST (IAU 2006) | done |
| Apparent / true / astrometric, geo / topo / helio / bary, ecliptic / equatorial, of date / J2000 | done |
| Speeds (numerical, consistent with positions) | done |
| 8 house systems, house_pos, polar rule + Swiss fallback helper | done |
| 7 ayanamsas (Pushya ≡ True Pushya, from δ Cnc with JPL states) | done |
| Mean / true node, mean / osculating Lilith | done |
| Searches: ingresses, returns, exact aspects, lunations, stations | done (core) |
| Python bindings (batch numpy, GIL released), `swe_compat` shim with per-task settings | done |
| Fixed stars: 410 stars / 429 names (Hipparcos 2007 + IAU WGSN names + SIMBAD RVs), rigorous space motion | done |
| Rise/set/transit, twilight, Hindu rising, pheno, solcross/mooncross, equation of time, cotrans, azalt, refrac | done |
| Solar & lunar eclipses: global, where, local; NASA canon 0–2400: 0 missed / 0 spurious / 0 type mismatches | done |
| Ceres, Pallas, Juno, Vesta (JPL sb441-n16, 0–2400); Chiron, Pholus, Eris, Eros (JPL Horizons 1600–2400, extended by our own integration from the JPL orbits: Chiron from 675, Pholus from 0, Eris/Eros from 1500 — same coverage as Swiss) | done (own Chebyshev fit ≤ 0.07 km) |
| Output flags XYZ / RADIANS / ICRS / ASTROMETRIC / SPEED3, J2000 frame bias, `ECL_NUT`, planet-centric `calc_pctr`, `helio_cross`, `mooncross_node`, `lun_eclipse_when_loc`, UTC with IERS leap seconds, per-engine Delta T controls (user Delta T, tidal acceleration), calendar / `split_deg` / centisecond helpers, names (TIME_FRAMES §10, GEOMETRY, EVENTS, ECLIPSES "Quick wins") | done (2026-10-29) |

## 2. Physical accuracy (independent truth: JPL Horizons, DE441)

Astrometric ICRF RA/Dec, geocentric, TT (no Delta T, no frame-model dependence),
max error over 25 epochs per body per era, arcsec:

| Era | FrenchEphemeris | Swiss 2.10 |
|---|---|---|
| 0–1550 (same DE441 data) | ≤ 0.0005 | ≤ 0.006 |
| 1550–2400, inner planets / Sun / Moon | ≤ 0.13 (DE442 vs DE441) | ≤ 0.10 |
| 1550–2400, Uranus–Pluto | ≤ 1.76 (DE442 vs DE441) | ≤ 1.24 |

After 1550 FrenchEphemeris uses DE442 (newer Uranus/outer-planet fit) while
Horizons serves DE441, so the post-1550 column measures the DE442−DE441
difference, not an error of either. Horizons' "apparent ecliptic of date" uses
older precession/nutation conventions and differs from both engines by up to
3″ at year 0; frame matrices are instead validated against ERFA to ≤ 4e-16 rad.

## 3. Swiss compatibility (as a chart sees it: UT input, planets then houses)

Decision rule (product owner): where a difference would be visible to an astrologer
comparing with Swiss (displayed 0.001° = 3.6″, times to the second, classifications),
get closer to Swiss using published methods only; where we are only more precise and
the difference is invisible, keep ours. No Swiss-derived numbers ship.

Max |FrenchEphemeris − Swiss|, arcsec, 200 random charts per era, Python defaults:

| Era | Sun | Moon | Mercury | Mars | Jupiter | Uranus | Pluto | Mean node | True node | Mean Lilith | Lahiri Moon | Placidus cusps |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0–1000 | 0.003 | 0.034 | 0.010 | 0.010 | 0.006 | 0.002 | 0.049 | 1.18 | 2.12 | 2.90 | 0.041 | 1.01 |
| 1000–1550 | 0.002 | 0.005 | 0.002 | 0.004 | 0.001 | 0.001 | 0.002 | 0.14 | 0.08 | 2.68 | 0.013 | 0.23 |
| 1550–1800 | 0.001 | 0.47 | 0.002 | 0.033 | 0.004 | 1.95 | 0.012 | 0.10 | 0.06 | 1.33 | 0.46 | 0.30 |
| 1800–2100 | 0.066 | 0.83 | 0.10 | 0.088 | 0.025 | 0.62 | 0.64 | 0.033 | 0.039 | 0.54 | 0.84 | 4.9 (after 2050) |
| 2100–2400 | 2.1 | 32 | 4.3 | 1.6 | 0.45 | 1.4 | 1.6 | 0.047 | 0.63 | 0.58 | 32 | 9.5 |

No sign or retrograde flips in 8,400 samples. What made historical charts match:

1. **Delta T** (default `swiss_delta_t=True`): SMH 2016 spline corrected to the DE441 lunar
   tidal acceleration −25.936″/cy² (JPL Horizons / Giorgini, quoted in Swiss's programmer
   docs) — physically consistent with the DE441 Moon we use, and what Swiss applies inside
   charts before 1800. (Swiss's bare `swe_deltat` call reports −25.80-based values until a
   pre-1800 planet is computed — another call-history dependence.) The newer SMH 2020
   spline remains available (`swiss_delta_t=False`).
2. **Sidereal time outside 1850–2050** (default `swiss_sidereal_time=True`): the method
   Swiss documents (swisseph.htm §6.1: Simon et al. 1994 mean Sun on Vondrák 2011
   precession), implemented from the publications, with Delta T at the DE431 tidal value
   (Swiss's houses use it at all dates). Before 2050: ARMC within 0.1″, cusps ≤ 1″.
3. **Ayanamsas**: published zero points with the precession-reconversion constants and the
   Lahiri refined value documented by Swiss (GEOMETRY_REPORT.md): ≤ 0.008″.
4. **Star names**: six names resolve to the stars Swiss uses; α Cen = A+B barycentre.

Visible residuals that cannot be closed without copying Swiss numbers:

- **After 2100**: Moon up to 32″ (≈ 0.009°) — Swiss's Delta T prediction curve matches no
  published formula; ours uses the IERS Bulletin A predictions for the next 13 months, then
  the last 10 years' rate with the SMH long-term acceleration (FE-DT/2, TIME_FRAMES_REPORT §2).
  Nobody knows the future Earth rotation.
- **Houses after 2050**: up to 10″ (≈ 0.003°) — Swiss adds an unpublished continuity
  constant (1.95″ in ARMC) at 2050; ours is continuous.
- **Mean Lilith**: 0.25″ in 1900–2050 (last digit differs in ~3% of values) — Swiss uses
  Moshier's 1992 fit plus unpublished corrections. A refit of ELP/MPP02 to the DE441 Moon was
  tried and moves *away* from Swiss (0.40″): Swiss's mean elements follow the nominal published
  theory, not JPL (GEOMETRY_REPORT.md). Kept: published ELP/MPP02.
- **Chiron before ~1000 / Pholus before 1600**: Chiron agrees with Swiss to 0.3″ in
  1240–1600 but differs by ~47″ in 900 and >1° before 720, when it passed 0.18 AU from
  Saturn. Before that encounter its position is physically undetermined in every engine
  (our 3σ from JPL's covariance: 20″ at 1000, 2′ at 750, ~1° at 675; Pholus up to ~10′).
  `frenchephemeris.minor_uncertainty` now returns the 1-sigma error and an `is_uncertain`
  flag for Chiron, Pholus, Eris and Eros, so callers can show these points as uncertain.
- **Rise/set**: closed. Swiss's refraction (written with V. Reijs) is Sinclair's 1982
  low-altitude formula as printed in the Explanatory Supplement; implemented from that
  publication it matches Swiss to 0.00″ at 0–16° altitude, and sunrise/sunset to ≤ 0.18 s.

Other measured differences:

- Speeds: ours equal the numerical derivative of our positions (≤ 0.0005″/day);
  Swiss's Moon speed differs from the derivative of its own positions by up to 0.18″/day.
- Topocentric Moon: ≤ 0.4″ — diurnal aberration (included by IERS/SOFA, apparently
  not by Swiss) and ~0.3 km observer-placement convention.
- Heliocentric Mercury/Venus: ≤ 1″ (Swiss heliocentric light-time convention not
  identified; heliocentric output is not used by FreeAstroAPI).
- Swiss is history-dependent: the first `calc_ut` after an ephemeris-file switch
  can differ by up to 26″ (Moon) from the same call repeated. FrenchEphemeris is
  a pure function of its inputs.

## 4. Speed (Python level, as FreeAstroAPI calls it)

| Case | Swiss | FrenchEphemeris |
|---|---|---|
| 10 bodies × 1,000 times, with speeds | 13.9 ms | 12.5 ms per-body calls; **1.87 ms batch (7.4×)** |
| positions only | 9.3 ms | **1.0 ms batch (9×)** |
| random-date 10-body chart, per-body calls | 55 µs | **15 µs (3.6×)** |
| Rust kernel, 10-body chart with speeds (1 core) | — | 10 µs |
| 30 years of equinoxes / new moons | — | 4 ms / 24 ms |

No process-wide lock: the engine is immutable and `Send + Sync`; batches are
split across cores and release the GIL.

**Throughput, FrenchEphemeris 0.1.5** (`python tools/bench_throughput.py --swiss`, Apple M4,
pyswisseph 2.10.3.2; 200,000 calls over 10 bodies at random dates 1900–2100, with speeds):

| positions/s | Swiss, 1 thread | Swiss, 8 threads | FrenchEphemeris, 1 thread | FrenchEphemeris, 8 threads |
|---|---|---|---|---|
| tropical, one body per call | 313,000 | 311,000 | 502,000 | 580,000 |
| Lahiri sidereal, one body per call | 262,000 | 263,000 | 386,000 | 398,000 |
| batch API (10 bodies × 10,000 times) | — | — | 4,460,000 | — |

Single positions keep the GIL (releasing it for a few microseconds of work made threads
contend), so per-call throughput stays flat with more threads; FrenchEphemeris keeps each
thread's sidereal mode and topocentre separate, where Swiss settings are process-wide.

## 4b. Other capabilities (details in the workstream reports)

| | vs Swiss | vs independent truth | speed vs Swiss |
|---|---|---|---|
| Fixed stars | median 0.003″ (p95 0.26″, all eras) | ERFA chain to 2e-4 µas | 25 stars 20 vs 34 µs |
| Rise/set (Sinclair 1982 refraction, Explanatory Supplement §3.28) | ≤ 0.18 s below 60° lat, ≤ 0.58 s overall; same no-event classification in 3276 cases; displayed second differs in 0.8% of sunrises | Horizons ≤ 0.5 s | 1.2–28× faster |
| Eclipses | maxima 1–2 s p50; 33 Swiss solar type errors (NASA agrees with us) | NASA canon: contacts ≤ 2 s from NASA elements | 4–8× faster |
| Minor bodies (1800–2100) | median 0.2–1″, Pholus 3″ (older Swiss orbits) | Horizons ≤ 0.004″ (Chiron/Pholus/Eris/Eros), ≤ 0.25″ (Ceres–Vesta) | — |

## 5. Open decisions

Decided (rule in §3): Swiss-compatible Delta T and long-term sidereal time are the Python
defaults; DE442 kept for 1550+ (Uranus seam 1.85″ is below display precision); polar
Placidus/Koch keep today's API error behaviour.
4. Polar Placidus/Koch: keep today's error (pyswisseph raises) or switch to the
   Porphyry fallback that the C library computes.
5. Ayanamsa zero points: strictly published values (branch ws/provenance; constant
   0.4–0.8″ offset vs Swiss for Lahiri/Raman/KP/Fagan) or an authoritative source matching
   today's results.
6. Chiron/Pholus/Eris/Eros before 1600: JPL publishes no data; Swiss serves Chiron from 675.
   The API requests Chiron by default, so pre-1600 charts need the adapter to drop/flag it.
7. Six star names resolve to different stars than Swiss (IAU vs Swiss naming); add
   continuity aliases or follow the IAU.
8. Licence confirmations: ESA Hipparcos / CDS SIMBAD commercial terms; IAU WGSN (CC BY)
   attribution; JPL attribution wording.
9. FreeAstroAPI bug found: Pholus is passed as 5145 (not a valid id), so it errors at every date.

## 6. Gap analysis steps 1–2 (2026-10-29)

| Area | Status vs Swiss |
|---|---|
| House systems | all 24 Swiss systems ('J' Savard-A from Savard's own published page, ≤ 4e-10°; 'i' Sunshine/Makransky differs in 76 high-latitude rows — method only in his copyrighted 1988 book); cusps ≤ 1e-7° (G/P near polar circle ≤ 3e-6°, Swiss's own iteration); `houses_ex2` speeds (true derivatives; several Swiss speed defects documented); Gauquelin sectors, 6 methods; `house_name`; sidereal projection-mode houses ≤ 2.4e-6° |
| Ayanamsas | 44 of 47 modes + user mode + all SIDBIT options. ≤ 0.001″ for the common modes; visible residuals: Djwhal Khul 4.9″, Eta Piscium 0.93″, Galactic-Centre family 0.10″ modern (now from the catalogue's Galactic Center record; the rest is Swiss's own Sgr A* star record), IAU1958 galactic equator 0.18″, Lahiri ICRC 0.07″. Not implemented (no published definition): J.N. Bhasin, Sheoran, Vettius Valens |
| Fixed stars | 1,176 objects; 1,706 of 1,708 Swiss-accepted names return Swiss's exact string (missing: `Apex`, `Gal.Equ`); positions median 0.0000″ at 2000; magnitudes median 0.000 mag (p95 0.02); user list in STAR_NAMES.md |
| Natural Lilith / Priapus | Swiss's documented construction on ELP 2000-82B (918 terms, authors' IMCCE files): max 0.013–0.023° in every era, rms 0.005–0.007° — visible; the rest is a 27.5-day term in Swiss's undocumented implementation |
| Hypothetical bodies | Uranians (Neely 1980), Harrington: ≤ 0.00002°; Kronos ≤ 0.08° (published a), Le Verrier/Adams/Lowell 0.1–0.6° (historical variants differ); not implemented (no published elements found): Isis-Transpluto, Nibiru, Pickering, Vulcan, Proserpina; White Moon and Waldemath from secondary/reconstructed sources (no oracle) |
| Quick wins | XYZ/RADIANS/ICRS/ASTROMETRIC/SPEED3, ECL_NUT, calc_pctr, UTC & leap seconds, utilities (identical), local lunar eclipses (flags identical), node/helio crossings, Delta T controls (scoped, not global) |

Swiss global-state defects found along the way: `PREC_ORIG` contaminates later ayanamsas;
Sunshine `house_pos` reads the Sun's declination left by the last houses call;
Hindu rising uses the previous caller's `set_topo`; `calc_ut`/`deltat` depend on which data
file the last call loaded.
