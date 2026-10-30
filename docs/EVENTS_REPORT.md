# Events and utilities: rise/set/transit, phenomena, crossings, equation of time

Module: `crates/fe-core/src/events.rs`.
Python: `crates/fe-py/src/events.rs` (module-level functions taking `Engine`),
`python/frenchephemeris/swe_compat.py`.
Tests: `crates/fe-core/tests/events.rs`, fixtures `crates/fe-core/tests/fixtures/events/`
(`gen_events.py` → `events.json`, `gen_horizons.py` → `horizons.json`).
Project rule: no table or constant derived from Swiss output ships; Swiss is used only to
observe semantics and to report deltas (see §5 Provenance).
Bench: `cargo run --release -p fe-core --example bench_events`.

Clean room: Swiss Ephemeris was used only through the public pyswisseph 2.10.3.2 API
to *measure* behaviour; algorithms come from Meeus (Astronomical Algorithms, ch. 13, 15,
16, 28), the Explanatory Supplement, Bennett (1982), Sæmundsson (1986), Mallama & Hilton
(2018) and IAU/SOFA conventions.

## 1. What FreeAstroAPI calls (grep of `core/`, `routers/`, worktrees)

| Call site | Call | Used output |
|---|---|---|
| `core/vedic_astrology.py` sunrise/sunset (panchang) | `rise_trans(jd, SUN, CALC_RISE / CALC_SET, (lng, lat, 0), 0, 0, FLG_SWIEPH)` | `tret[0]` |
| `core/astrology.py` moonrise/moonset | `rise_trans(jd, MOON, CALC_RISE / CALC_SET, (lon, lat, 0), 0, 0, FLG_SWIEPH)`, `res != 0` → None | `res`, `tret[0]` |
| `core/vedic_strength_v2.py` (Shadbala day/night) | `rise_trans(jd, SUN, geopos=(lng,lat,0), rsmi=BIT_HINDU_RISING\|BIT_NO_REFRACTION\|BIT_DISC_CENTER\|CALC_RISE/SET, flags=SWIEPH\|SIDEREAL\|TRUEPOS\|NONUT\|SPEED)` | `tret[0]` |
| `core/astrology.py` moon phase | `pheno_ut(jd, MOON, FLG_SWIEPH)` | `[1]` illuminated fraction |
| `core/vedic_jaimini_v2.py`, `natal_pkg/data.py` | `solcross_ut(natal_sun_lon, jd, FLG_SWIEPH)` | JD UT |
| `core/chinese_astrology.py` true solar time | `time_equ(jd_ut) * 1440` minutes; LMT = UTC + lon/15 h; true solar = LMT + EoT | E (days) |
| primary directions, astrocartography (worktrees) | `cotrans((lon, lat, 1), -eps)`; `sidtime(jd_ut)` | RA/Dec, GST |

So the API needs: default-flag rise/set of Sun and Moon at height 0 with pressure 0 and
temperature 0, Hindu rising of the Sun, the −2 "no event" code, Moon illuminated fraction,
solar crossings, EoT, cotrans. Everything else (planets, P/T, horizon height, twilight,
azalt, refrac, mooncross, TT variants) is implemented for the compat layer.

## 2. Semantics (measured black-box, then implemented)

### rise_trans
* **Place**: topocentric apparent place of date (true equinox), whatever the `flags`
  (SIDEREAL, TRUEPOS, NONUT, TOPOCTR change the Swiss result by exactly 0.0 s).
  Altitude from the hour angle with GAST + east longitude and *geodetic* latitude.
* **Event condition**: the chosen disc point reaches the true altitude
  `h0 = horhgt − R(horhgt)`. Upper limb by default (`BIT_DISC_CENTER` → centre,
  `BIT_DISC_BOTTOM` → lower limb); semidiameter = asin(R/Δ) with the IAU mean radii
  (`body_radius_km`, IAU WGCCRE 2015; Swiss uses 1737.5 km for the Moon, 0.1″ in diameter); `BIT_FIXED_DISC_SIZE`: Sun at
  1 AU, Moon at 384 400 km, planets unchanged.
* **Refraction** (default `RefractionModel::Sinclair1982`, used by rise/set and azalt):
  Sinclair (1982), as printed in the *Explanatory Supplement to the Astronomical Almanac*
  (1992) §3.28 (HM Nautical Almanac Office low-altitude formula), h = **apparent** altitude:

      R0 = (34.46 + 4.23 h + 0.004 h²) / (1 + 0.505 h + 0.0845 h²)  arcmin   (h < hₓ)
      R0 = 0.97′ · cot h                                               (h ≥ hₓ)
      R  = R0 · (P − 80)/930 / (1 + 8·10⁻⁵ (R0 + 39)(T − 10))

  hₓ = 17.904° is where the two expressions meet (computed at run time by bisection, so the
  model is continuous); 0.97′ = 58.2″ is the classical standard-conditions (1010 hPa / 10 °C)
  refraction constant of the R = k·cot h asymptote. The rise/set horizon is R(0) = 34.46′ ×
  scaling = 36.74′ at 1013.25 hPa / 0 °C. True→apparent is the numerical inverse
  (fixed point h_a = h_t + R(h_a)); below the apparent horizon no refraction is applied.
  *How identified*: the classic API documents that its extended refraction was written with
  archaeoastronomer V. Reijs, whose ARCHAEOCOSMO pages describe "the inverse of Sinclair's
  formula"; the coordinator matched the ES §3.28 form black-box, and here it was verified
  against `refrac_extended(h, 0, P, T, 0.0065, APP_TO_TRUE)` on 1547 points (h −20…90°, seven
  P/T pairs from 500 hPa to 1100 hPa, −40…30 °C): **apparent→true ≤ 0.0000″, true→apparent
  ≤ 0.09″** (the oracle's own inversion tolerance). The high-altitude branch is 0.97′·cot h
  with the *same* (P − 80)/930 and temperature factors (the ES 0.00452° P tan z/(273+T) form
  is 0.06″–14″ off), and the switch is the crossing point, not 15°.
  Sæmundsson (1986), Bennett (1982), the simplified ES Sinclair form and the AA fixed 34′
  remain selectable. `pressure = 0` → standard-atmosphere pressure at the observer height;
  temperature as given (0 °C by default). The dip of the horizon is **not** applied
  (observed semantics). `BIT_NO_REFRACTION` → R = 0. A horizon below apparent altitude 0
  gets `min(horizon, −R(0))` (observed semantics at sea level).
* **Hindu rising** (`BIT_GEOCTR_NO_ECL_LAT`, part of `BIT_HINDU_RISING = 896`): geocentric
  apparent place with ecliptic latitude set to 0, rotated with the true obliquity.
* **Transits** (`CALC_MTRANSIT`/`CALC_ITRANSIT`): hour angle 0°/180° of the topocentric
  place; all `BIT_*` options are ignored (measured). Swiss quirk found: with the Hindu bits
  the transit branch uses the process-global `swe_set_topo` location instead of `geopos`
  (it errors "geographic position has not been set" in a fresh process, and otherwise
  silently uses the previous call's location). FrenchEphemeris always uses the given
  observer; the fixture generator calls `set_topo(geopos)` to make the oracle deterministic.
* **Twilight bits**: disc centre at true altitude −6/−12/−18°, no refraction (measured).
* **Search window / "no event"**: the next event after `tjd` within **26 h**. Measured on
  4000 random high-latitude cases: Swiss reports events up to +1.0795 d and returns −2 for
  everything from +1.0840 d on, without exception. Outcome `RiseSet::Event(jd)`,
  `AlwaysAbove` (circumpolar / polar day) or `AlwaysBelow` (polar night / never rises);
  both "no event" categories map to Swiss `res = −2, tret[0] = 0`. Transits always exist.
* **Algorithm**: an `EpochInterp` (3 exact epochs, quadratic interpolation of time scales and
  frame matrices; ≤ 0.25 mas over the window, unit-tested) feeds the altitude function;
  upper/lower transits in the window by secant iteration on the wrapped hour angle; altitude
  at window ends and transits; the first bracketing interval of the right sign is refined
  with `search::brent` (1e-8 d). Grazing events between two consecutive culminations are
  therefore always bracketed.

### pheno_ut
`[phase angle, phase, elongation, diameter, magnitude, horizontal parallax, 0 …]` (20 floats).
* Elongation: angle between the apparent geocentric Sun and body (as Swiss, < 0.3″).
* Phase angle: Moon from the apparent geocentric and apparent heliocentric vectors (Swiss
  definition, < 0.15″); planets geometric (body at t−τ, Sun at t−τ, observer at t) —
  Swiss's planet definition could not be reproduced by any apparent/true combination
  (10–35″ apart; ≤ 50 ppm in the illuminated fraction).
* Phase = (1 + cos i)/2. Diameter = 2 asin(R/Δ). Horizontal parallax (Moon) = asin(6378.137/Δ).
* Sun: Swiss returns phase angle 0, **phase 0** and elongation 0; reproduced for compatibility.
* Magnitude model (documented choice, not Swiss's): Mallama & Hilton (2018) for
  Mercury–Neptune (Saturn with rings from the IAU pole when α ≤ 6.5°, β ≤ 27°, globe
  formula otherwise; Uranus sub-latitude term; Neptune secular term; Mars without the
  orientation terms); Pluto H = −1.01; Moon Allen (1976) `0.21 + 5 log(rΔ) + 0.026α + 4e-9 α⁴`;
  Sun −26.74 + 5 log Δ.

### solcross_ut / mooncross_ut (and TT variants)
Next instant after `jd_ut` when the apparent tropical longitude (or per flags: sidereal,
true, J2000 …) equals x. Newton on `calc_ut` longitude/speed from a mean-motion guess.

### time_equ, lmt_to_lat, lat_to_lmt, true solar time
Measured: Swiss `time_equ` = **GAST − α☉(apparent, true equinox of date) − (UT − 12h)**
(identical to this formula evaluated with Swiss's own sidtime and RA to < 1 µs; *not* the
Meeus 28.3 mean-longitude form). `lmt_to_lat(jd, lon) = jd + E(jd − lon/360)`.
`Engine::true_solar_time(jd_ut, lon) = jd_ut + lon/360 + E(jd_ut)` is the helper for the
Chinese true-solar-time code (its `true_solar_absolute = LMT + EoT`).

### cotrans / cotrans_sp
Rotation about x by `eps` (eps < 0: ecliptic→equatorial). Exact (≤ 1e-12° vs Swiss).

### azalt / refrac (compat helpers)
azalt: azimuth from south through west, true altitude, apparent altitude with
Sinclair 1982 (the oracle's azalt uses the same model: ≤ 1.7″ above the horizon, the residual
being the sidereal-time offset). `refrac` in Rust and the native binding default to
Sinclair 1982 (`model=` selects another). Note: the oracle's *simple* `refrac` is a different
function — Sæmundsson below ~15° (< 0.0001°) and an unidentified form above (Sæmundsson
1.7″ p50 / 8.6″ max there; Sinclair differs from it by up to 121″ near the horizon) — so
`swe_compat.refrac` passes `model="saemundsson"` to keep drop-in behaviour.

## 3. Validation

Fixtures: 3276 rise/set/transit cases (Sun, Moon, Mercury–Pluto; latitudes −88.9…+88.7 with
1000 cases in 60–89°; 399 within 2° of the date line; years 0–2399; all `BIT_*` variants;
P/T/height/horizon variants; the three exact API call patterns), 600 pheno, 300 crossings,
400 EoT, 200 cotrans, 200 azalt. Swiss SE files cover −600…2399.

### rise/set/transit vs Swiss (|FE − Swiss|, seconds)

3276 cases (the sets above plus a dense daily scan of 876 Sun/Moon rise/sets around the
polar-day/night and Moon-circumpolar boundaries at 62–70.5° N/S).

Geometry is unchanged by the refraction model — no refraction, Hindu rising and transits:

| |lat| band | n | p50 | p95 | max |
|---|---:|---:|---:|---:|
| 0–30 | 275 | 0.001 | 0.15 | 0.53 |
| 30–50 | 202 | 0.001 | 0.14 | 0.52 |
| 50–60 | 89 | 0.002 | 0.13 | 0.40 |
| 60–70 | 168 | 0.010 | 0.22 | 0.40 |
| 70–89 | 62 | 0.008 | 0.13 | 0.15 |

Refracted rise/set (Sinclair 1982, 1900–2100):

| |lat| band | n | p50 | p95 | max |
|---|---:|---:|---:|---:|
| 0–30 | 408 | 0.003 | 0.14 | 0.17 |
| 30–50 | 301 | 0.005 | 0.14 | 0.18 |
| 50–60 | 125 | 0.020 | 0.16 | 0.18 |
| 60–70 | 817 | 0.005 | 0.10 | 0.58 |
| 70–89 | 77 | 0.004 | 0.08 | 0.14 |

(Previous published-Sæmundsson default: p95 3.2 / 4.7 / 6.5 / 30 / 44 s, max 225 s.)
By API call pattern: vedic sunrise/sunset p50 0.002 s, p95 0.13 s, max 0.16 s;
moonrise/moonset p50 0.008 s, p95 0.17 s, max 0.18 s; Hindu rising max 0.13 s.
Non-default atmospheres (horizon 0) ≤ 0.20 s; non-zero horizons ≤ 0.17 s; the boundary scan
Sun ≤ 0.008 s, Moon ≤ 0.58 s. Negative horizon for elevated observers ~200–300 s
(dip-aware refraction, not modelled).

**Classification**: 0 mismatches in 3276 cases (256 Swiss −2 results, all reproduced as
`AlwaysAbove`/`AlwaysBelow`; no spurious −2), including the daily boundary scan; with the
horizon now identical to the oracle's, flips can only come from the sub-second geometry
residual.

Long range (years 0–1900) with `EngineConfig::SWISS_TIME` and no correction: refracted
Sun ≤ 0.13 s, planets ≤ 0.27 s, Moon ≤ 3.5 s; geometric Sun ≤ 0.14 s, planets ≤ 0.28 s,
Moon ≤ 2.1 s (lunar ephemeris DE441 vs DE431).

Notes:
* The ~0.13 s ceiling of the geometric groups is Swiss's sidereal time drifting from the
  IAU ERA-based GAST after ~2050 (−0.127 s in 2060) and its root tolerance (one case,
  Delhi 2025-06-21 Hindu sunrise, 0.35 s: Swiss leaves the Sun at +4.2″).
* Swiss `sidtime` differs from IAU 2006 GAST by +2.1 s (year 1050), +14.8 s (340), +20.8 s
  (38) (measured vs ERFA `gst06a`); FrenchEphemeris follows IAU. *Affects houses/ARMC too.*

### Independent reference: JPL Horizons (15 events, 1962–2024, ±70° lat, date line)
Horizons airless apparent topocentric elevation + angular diameter at 1-min steps (DE441,
its own UT1), solved for `elev + sd = −34′` (AA convention, a fixed published horizon;
FE run with `HorizonRefraction::Almanac34`), so this checks geometry and time scales
independently of refraction: **FE − Horizons: Sun p50 0.16 s, max 0.46 s; Moon p50 0.06 s,
max 0.49 s.** The largest are Sydney/McMurdo cases (~0.45 s), probably Horizons' Earth
figure. No case > 1 s.

### pheno vs Swiss (1900–2030; arcsec / ppm / mag)

| quantity | Moon | planets (worst) |
|---|---|---|
| phase angle | 0.13″ | 35″ (Mars), 27″ (Jupiter, Saturn) — definition differs |
| phase (illuminated fraction) | 0.25 ppm | 48 ppm (Mercury) |
| elongation | 0.19″ | 0.27″ |
| diameter | 0.10″ (Moon radius 1737.4 vs Swiss 1737.5 km) | 0.0000″ |
| horizontal parallax | 0.35″ | — |
| magnitude | 0.021 (α < 150°); 0.5–3.5 at α > 150° (thin crescent: Allen vs Swiss's own) | 0.049 (Uranus), 0.042 (Saturn), 0.017 (Neptune), 0.010 (Pluto), ≤ 0.001 (Mercury–Jupiter); Sun 0.12 |

Mercury–Jupiter magnitudes agree to < 0.001 mag, i.e. Swiss also uses Mallama & Hilton
(2018) for those; Saturn/Uranus differ by ring/latitude handling, Sun by the adopted
V☉ (Swiss −26.86). Outside 1900–2030 the Moon differs by up to ~100″ purely from Delta T.

### Other functions vs Swiss

| function | 1900–2100 | 0–1900 |
|---|---|---|
| solcross_ut (UT) | median 0.03 s; max 48 s = ΔT extrapolation after 2026 | ≤ 178 s (ΔT) |
| solcross in TT (ΔT removed) | max 0.32 s (0.013″ of longitude) | 0.05 s |
| mooncross in TT | max 0.40 s | 0.74 s |
| time_equ | max 0.68 s raw; 0.013 s after ST/ΔT offsets | 21 s raw (Swiss ST); 0.15 s corrected |
| lmt_to_lat / lat_to_lmt | as time_equ; round trip < 1e-4 s | |
| cotrans / cotrans_sp | < 1e-12° / < 1e-12°/d | |
| azalt true altitude | ≤ 1.8″ (Swiss ST after 2050) | |

### Timings (release, `bench_events`, Apple Silicon, single thread)

| call | FrenchEphemeris | pyswisseph 2.10 (incl. ~0.2 µs Python call) |
|---|---:|---:|
| Sun rise (default) | 15 µs | 18 µs |
| Sun rise (Hindu) | 16 µs | 15 µs |
| Moon rise | 19 µs | 33 µs |
| Jupiter rise | 15 µs | 18 µs |
| upper transit | 8 µs | 23 µs |
| Sun rise at 69.65°N (polar day/night) | 15 µs | 426 µs |
| pheno_ut Moon / Saturn | 3.7 / 3.3 µs | 6.6 / 6.6 µs |
| solcross_ut / mooncross_ut | 22 / 25 µs | 23 / 18 µs |
| time_equ | 3.6 µs | 4.4 µs |

(Timings measured with the earlier refraction model; Sinclair 1982 is a closed form at the
horizon, no slower.)
Half of a rise/set is the three exact epochs (2.5 µs each: Delta T, Vondrák precession,
nutation, CIO locator); each altitude evaluation is 0.43 µs. All functions are `&self`
on the immutable `Engine` and run without the GIL from Python.

## 4. Open issues

1. **Swiss sidereal time in the past/future**: now available as the documented long-term
   model `SiderealTime::SwissLongTerm` (`EngineConfig::SWISS_TIME`, Python default; see
   TIME_FRAMES_REPORT §7.1). Long-range (0–1900) unrefracted rise/set/transits with it and no
   correction: Sun ≤ 0.14 s, planets ≤ 0.28 s, Moon ≤ 2.5 s (`long_range_rise_set_swiss_time`);
   after 2050 a constant 1.95″ (0.13 s) remains (Swiss's unpublished continuity constant).
   Events now use the engine's ΔT model. The default Rust engine keeps IAU GAST.
2. **Refraction**: resolved — Sinclair (1982) / ES (1992) §3.28 reproduces the oracle's
   extended refraction exactly; refracted rise/set ≤ 0.18 s below 60°, ≤ 0.58 s overall.
   Still open: negative horizons for elevated observers (dip-aware refraction) are not
   modelled; the source text of the 0.97′·cot h high-altitude branch (classical 58.2″
   standard constant) was not located in the ES — it was identified as the oracle's form
   black-box and only matters above 17.9° (azalt/refrac, not rise/set).
3. Moon magnitude near new Moon (α > 150°) differs from Swiss by up to 3.5 mag (Allen's
   formula is invalid there; neither is authoritative). Sun magnitude −26.74 vs Swiss −26.86.
4. Planet phase-angle definition differs from Swiss by ≤ 35″ (≤ 50 ppm in the phase).
5. Fixed stars: `Engine::rise_set_target(&TargetFn, …)` accepts any `Fn(&Epoch) → (RA, Dec,
   semidiameter)`; `swe_compat.rise_trans` raises `NotImplementedError` for star names until
   the stars workstream provides the place function.
6. `swe_compat.rise_trans` ignores `flags` (as Swiss does) and returns `tret` padded to 10.
   The native `rise_trans` additionally returns the category
   (`"event"`, `"always_above"`, `"always_below"`).
7. Python bindings were compile-checked (`cargo build -p fe-py`) and smoke-tested by loading
   the debug `lib_native.dylib` from a scratch copy of the package (no `maturin develop`).

## 5. Provenance of every numerical constant in `events.rs`

No value below was fitted to, or tabulated from, Swiss output.

| constant | value | source |
|---|---|---|
| Sæmundsson refraction | 1.02′/tan(h + 10.3/(h + 5.11)) | Sæmundsson, Sky & Telescope 72, 70 (1986); Meeus, *Astronomical Algorithms* 2nd ed., eq. 16.4 |
| Bennett refraction | 1′/tan(h + 7.31/(h + 4.4)) | Bennett, J. Navigation 35, 255 (1982); Meeus eq. 16.3 |
| P/T scaling | (P/1010)·(283/(273 + T)) | Meeus ch. 16 |
| Sinclair 1982 low-altitude (default) | R0 = (34.46 + 4.23h + 0.004h²)/(1 + 0.505h + 0.0845h²)′ | Sinclair (1982), NAO Technical Note 59; *Explanatory Supplement* (1992) §3.28 |
| Sinclair P/T factors | (P − 80)/930; 1/(1 + 8·10⁻⁵ (R0 + 39)(T − 10)) | same; identified via Reijs (ARCHAEOCOSMO) and verified black-box (§2) |
| high-altitude branch | 0.97′ cot h (= 58.2″ cot h at 1010 hPa / 10 °C), switch at the crossing 17.904° | classical standard refraction constant; branch form identified black-box (source text not located) |
| simplified Sinclair (selectable) | P(0.1594 + 0.0196a + 0.00002a²)/((273+T)(1 + 0.505a + 0.0845a²)) | Explanatory Supplement (1992) eq. 3.283-2 |
| AA horizon | 34′ refraction (+16′ semidiameter = −50′ centre) | Astronomical Almanac, rise/set convention |
| standard pressure | 1013.25·(1 − 0.0065h/288)^5.255 hPa | ICAO / US Standard Atmosphere 1976 troposphere |
| radii (km) | Moon 1737.4, Mercury 2439.4, Venus 6051.8, Mars 3389.5, Jupiter 69911, Saturn 58232, Uranus 25362, Neptune 24622, Pluto 1188.3 | IAU WGCCRE 2015 (Archinal et al. 2018, Celest. Mech. 130:22) mean radii |
| Sun radius | 696 000 km | Allen, *Astrophysical Quantities* (R☉ = 6.960·10⁵ km); IAU 2015 B3 nominal is 695 700 km |
| Earth radius | 6378.137 km | WGS84 / IERS |
| Moon mean distance (fixed disc) | 384 400 km | conventional mean distance (Explanatory Supplement) |
| twilight altitudes | −6°, −12°, −18° (centre, no refraction) | standard civil/nautical/astronomical definitions (AA) |
| planetary magnitudes | coefficients per planet | Mallama & Hilton, Astron. Comput. 25, 10 (2018) |
| Saturn / Uranus poles | (40.589°, 83.537°), (257.311°, −15.175°) | IAU WGCCRE 2015 |
| Neptune secular term | −6.89 → −7.00 (1980–2000) | Mallama & Hilton (2018) |
| Pluto H | −1.01 | Astronomical Almanac |
| Moon magnitude | 0.21 + 5 log(rΔ) + 0.026α + 4·10⁻⁹α⁴ | Allen, *Astrophysical Quantities* (1976) |
| Sun magnitude | −26.74 | Allen (1976) / AA |
| sidereal rate (transit iteration seed) | 1.00273790935 | IERS / Aoki et al. (1982) |
| `RISE_SET_WINDOW_DAYS` | 26 h | *semantic* choice (what counts as "no event"), set to reproduce the classic API's classification; a round value, not fitted |
| `refrac` domain | models are functions of the apparent altitude ≥ 0; none below the apparent horizon | observed semantics for a sea-level observer (matches the oracle's extended refraction) |
| interpolation / tolerances | 3-node epochs, 1e-9 d, 1e-8 d | implementation choices |

Definitions (not numbers) adopted from observed semantics: topocentric apparent place
whatever the flags; upper limb default; horizon < 0 rule `min(horizon, −R(0))`; Hindu
rising = geocentric, ecliptic latitude 0; transits ignore disc/refraction options; EoT =
GAST − α☉ − (UT − 12h); Moon phase angle from apparent geocentric and apparent heliocentric
vectors; Sun pheno returns phase 0.

## Quick wins (2026-10-29): `helio_cross`, `mooncross_node`

Tests: `tests/quickwins.rs::crossings_vs_swiss` (60 / 80 random cases, 1900–2050).

- `mooncross_node(_ut)`: next zero of the Moon's apparent geocentric latitude (1-day scan,
  Brent). vs Swiss: time ≤ 0.034 s, longitude ≤ 0.017″, latitude ≤ 0.0008″. The oracle's
  `mooncross_node_ut` (pyswisseph 2.10.03) returns exactly the TT-variant instant for the same
  number; we keep the documented UT-in/UT-out semantics.
- `helio_cross(_ut)`, forward and backward: monotonic heliocentric longitude scan (≤ 5° steps at
  the starting speed) + Brent. The crossing times differ from Swiss by the pre-existing
  heliocentric longitude difference (≤ 1.01″ inner, ≤ 0.51″ outer planets; Swiss's
  heliocentric light-time convention is still unidentified, see ENGINE_REPORT §3) divided by
  the heliocentric speed: ≤ 11 s for Mercury–Mars, up to 72 min for Pluto (0.004°/day).
