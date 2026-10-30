# Eclipses — method and validation

Code: `crates/fe-core/src/eclipses.rs` (kernel), `crates/fe-py/src/eclipses.rs`
(bindings), `python/frenchephemeris/swe_compat.py` (pyswisseph-shaped wrappers).
Tests: `crates/fe-core/tests/eclipses.rs`; fixtures and their generator in
`crates/fe-core/tests/fixtures/eclipses/` (`gen_eclipses.py`).

## What the API uses

From `astro/core/moon_eclipse.py` and `core/moon_forecast.py` (`routers/v1/moon.py`
calls `get_eclipse_info`):

| call | fields read |
|---|---|
| `lun_eclipse_when(jd, FLG_SWIEPH, 0)` | flags `ECL_TOTAL/PARTIAL/PENUMBRAL`; `tret[0]` max; `tret[2,4,6]` / `[3,5,7]` (earliest start / latest end = P1 / P4); all contacts for visibility |
| `lun_eclipse_how(jd, (lon, lat, 0))` | `attr[0]` umbral magnitude at max; `attr[6]` apparent Moon altitude at each contact (visible if > −0.833°) |
| `sol_eclipse_when_glob(jd, FLG_SWIEPH, 0, False)` | flags `ECL_ANNULAR_TOTAL/TOTAL/ANNULAR/PARTIAL`; `tret[0]` max, `[2]`/`[3]` begin/end, `[4]`/`[5]` totality begin/end |
| `sol_eclipse_where(jd_max)` | `geopos[0..1]` lon/lat; `attr[8]` NASA magnitude; `attr[2]` obscuration |
| `sol_eclipse_when_loc(jd, (lon, lat, 0), …)` | flags type + `ECL_VISIBLE`; `tret[0]` local max, `[1]` C1, `[4]` C4; `attr[8]`, `attr[2]` |

All are now provided by `swe_compat` with the same argument order and return shapes
(ephemeris flags ignored; times JD UT).

## Method

Clean-room, from published methods; Swiss Ephemeris was used only as a black box
through its Python API.

**Ephemeris input.** JPL DE441/DE442 states from the pack: Earth and Sun barycentric,
Moon geocentric, three evaluations per instant. Light time is applied linearly
(Sun: 500 s at 0.01 km/s; Moon: 1.3 s), leaving errors below 1e-5 km.

**Candidates.** Each mean lunation (Meeus ch. 49) is checked against the node.
Lunations with |sin F| > 0.5 are rejected without any ephemeris work. F more than
30° from a node is far outside the ≈ 21° limit, even counting the up-to-0.6-day
offset between mean and true phase.

**Lunar eclipses.**
- Greatest eclipse is the minimum angular distance σ between the Moon and the
  shadow axis. It is found by Gauss–Newton on the chord vector with analytic rates
  (3–5 iterations).
- Shadow geometry: the Moon is seen at t = T − τ☾. The shadow it is in was cast by
  sunlight passing Earth at t − τ☾. So **m = Moon(t) − Earth(t − τ☾)** and the
  axis is **Earth(t − τ☾) − Sun(t − τ☾ − τ☉)**.
- To better than 1″ this equals the classical "apparent Moon against apparent Sun
  + 180°" used by NASA and Swiss. A first version that used astrometric vectors
  (light time only, which is right for the solar case) put the Moon 42″ off and
  made every maximum about 80 s late.
- Shadow radii use **Danjon's rule** (Earth radius enlarged by 1 %):
  ρ_u = 1.01·π☾ + π☉ − s☉ and ρ_p = 1.01·π☾ + π☉ + s☉, with Moon k = 0.2725076.
- Danjon is also what the oracle uses. A black-box fit of the shadow radii implied
  by Swiss magnitudes gives 1.0100 (umbra) and 1.0097 (penumbra). Chauvenet's
  2 % (×1.02 on all terms) fits 1.0138 and 1.0076, which is inconsistent.
- Magnitudes are (ρ − σ + s☾)/(2 s☾). Contacts solve σ = ρ ± s☾ (Brent, 1e-8 d).

**Solar eclipses — Besselian elements evaluated directly** at each instant
(Explanatory Supplement 3rd ed. ch. 11; Meeus, *Elements of Solar Eclipses*), not
from fitted polynomials.
- Inputs are geocentric astrometric Sun and Moon. For an observer, aberration maps
  both directions the same way, so it cancels in every alignment.
- Shadow axis ẑ = (☉ − ☾)/|☉ − ☾|. x̂ and ŷ are taken relative to the true pole of
  date. From these come x, y, z, d, f1, f2, l1 = z·tan f1 + k1/cos f1 and
  l2 = z·tan f2 − k2/cos f2.
- k1 = 0.272488 and k2 = 0.272281, the constants printed on NASA's
  Besselian-element pages.
- Greatest eclipse is the minimum of √(x² + y²), i.e. the axis closest to the
  geocentre (NASA's definition). Gamma is that distance signed by y.
- Earth's silhouette is the ellipse with semi-axes 1 and ρ1 = √(1 − e² cos² d)
  (WGS84). Contacts use the **exact point–ellipse distance**:
  - P1/P4: distance = l1.
  - U1/U4: distance = |l2|.
  - C1/C2 (central line): distance = 0.
- Type:
  - central if the axis meets the ellipsoid;
  - total or annular from the sign of L2 = l2 − ζ·tan f2 at the central point;
  - **hybrid** if L2 changes sign along the central line (transition times h1/h2);
  - non-central total/annular if only the umbra or antumbra touches Earth;
  - otherwise partial.
- Central-line point: the axis is intersected exactly with the WGS84 ellipsoid.
  When the axis misses Earth, the nearest limb point is used instead.

**Local circumstances.**
- Topocentric WGS84 site, with Earth rotation from ERA and fixed NPB (< 0.01″ per
  hour). Contacts are roots of sep − (s☉ + s☾(k1)) and sep − |s☾(k2) − s☉|.
- Lunar `lun_eclipse_how` gives the Moon's topocentric altitude, plus Sæmundsson
  refraction (1013.25 hPa, 10 °C). Like the oracle, refraction is dropped when the
  refracted body would still be below the horizon.

**Constants.**

| quantity | value |
|---|---|
| Earth equatorial radius | a = 6378.137 km |
| Earth flattening | f = 1/298.257223563 (WGS84 = IAU 2015 nominal) |
| Sun radius | R☉ = 696 000 km (959.63″ at 1 au, the NASA canon value; IAU 2015 nominal is 695 700 km) |
| Moon k, lunar eclipses | 0.2725076 |
| Moon k1, solar penumbra | 0.272488 |
| Moon k2, solar umbra | 0.272281 |
| Danjon factor | 1.01 |

### Swiss behaviour reproduced deliberately (measured black-box)

**"Next" eclipse.** The search returns the eclipse whose greatest eclipse is more
than 1e-4 d after the start, or before it when searching backward.
`SEARCH_MARGIN_DAYS` holds this; Swiss skips an eclipse until 1e-4 d before its
maximum.

**Local solar visibility.** The oracle uses two different horizons:
- **Sunrise/sunset** (`tret[5]`, `tret[6]`) is placed where the Sun's true centre
  altitude is s☉ + 0.053° − 0.665°·P/P₀, with the standard-atmosphere pressure at
  the site height. That is −0.349° at sea level in June: roughly the apparent
  *lower* limb on the horizon. The rule was fitted over 0–4000 m and gives
  sunrise/sunset within 1.2 s p50 and 3.4 s p95.
- **Contact and maximum visibility flags** use "apparent centre altitude > 0", i.e.
  true altitude ≳ −0.58°.
- **The eclipse counts as visible**, and the search does not skip it, when the
  maximum is visible or the Sun is above the sunrise horizon at some time in
  [C1, C4].
- **Maximum below the horizon.** The reported maximum moves to the sunrise or
  sunset, and `attr` and the local type describe that visible instant.

**Ordering.** `tret` and `attr` follow the pyswisseph layouts:
- Azimuth is measured from south through west.
- `attr[3]` is the core-shadow diameter, negative for the umbra.
- Magnitudes in `lun_eclipse_how` are clamped at 0.
- Not computed: solar `tret[1]` (local-noon time), the limit coordinates in
  `geopos[2..9]`, and saros numbers. These return 0 and −99999999 respectively.

## Validation

### NASA Five Millennium Canons, years 0–2400 (compared in TT)

Sources (Espenak & Meeus, NASA/TP-2008-214170 and NASA/TP-2009-214173):

- https://eclipse.gsfc.nasa.gov/5MCSE/5MKSEcatalog.txt (solar)
- https://eclipse.gsfc.nasa.gov/5MCLE/5MKLEcatalog.txt (lunar)

Comparing "TD of greatest eclipse" with our TT maximum removes NASA's Delta T. Our
own engine searched the whole period with chained "next eclipse" calls. Every
catalogue entry was matched to within 0.5 d.

| | solar | lunar |
|---|---|---|
| catalogue eclipses 0–2400 | 5694 | 5801 |
| found / missed / spurious | 5694 / **0** / **0** | 5801 / **0** / **0** |
| type mismatches (P/A/T/H, N/P/T) | **0** | **0** |
| greatest eclipse − NASA, p50 / p95 / max | 1.1 / 23.7 / 34.4 s | 1.4 / 23.3 / 34.3 s |
| gamma − NASA, max | 0.0002 | 0.0002 |
| magnitude − NASA, p95 / max | 0.0003 / 0.0005 | umbral 0.0002 / 0.0004, penumbral 0.0003 / 0.0004 |

The time residual grows with distance from the present. It is the ELP2000-82/85 vs
DE441 lunar ephemeris difference (secular acceleration), not Delta T; it is 1 s
around 2000.

NASA's magnitude for *non-central* total/annular eclipses is the covered-diameter
fraction at the limb point, and `SolarEclipse.magnitude` follows that.
`sol_eclipse_where` `attr[8]` keeps the Swiss diameter-ratio convention.

### NASA polynomial Besselian elements (15 eclipses 2016–2035)

Pages such as https://eclipse.gsfc.nasa.gov/SEbeselm/SEbeselm2001/SE2024Apr08Tbeselm.html
(VSOP87/ELP2000-85), sampled at t0 + {−2, −1, 0, 1, 2} h:

| quantity | p50 | p95 | max |
|---|---|---|---|
| x, y | 1e-4 Earth radii (0.6 km) | | |
| l1, l2, tan f1, tan f2 | < 1e-5 | | |
| d | 6″ | 8″ | |
| μ | 20″ | 22″ | |
| greatest-eclipse point vs NASA (lat/long to 0.1′, 13 central eclipses) | 0.33 km | 0.60 km | 2.0 km |

μ is NASA's ephemeris hour angle. The ~20″ in d and μ is NASA building the frame
from aberrated directions. It is a rotation of the fundamental plane that shifts
ground points by < 1 km.

Contacts derived from NASA's own elements with the same exact-ellipse rule
(ephemeris-independent definitions):

| contact | ours − NASA elements, p50 / max | Swiss − NASA elements, p50 / max |
|---|---|---|
| P1/P4 | 0.37 / 0.84 s | 6.1 / 24.4 s |
| U1/U4 | 0.43 / 2.1 s | 3.9 / 12.9 s |
| C1/C2 | 0.40 / 2.1 s | 4.5 / 16.7 s |

Our global contacts are therefore more accurate than the oracle's. A spherical-Earth
variant matched Swiss worse (P1 p50 18 s), so Swiss uses some other approximation.

### Swiss black-box (pyswisseph 2.10.03), same 11 495 eclipses

| | solar | lunar |
|---|---|---|
| max − Swiss (TT), p50 / p95 / max | 1.0 / 7.5 / 10.7 s | 2.2 / 6.8 / 10.9 s |
| magnitudes − Swiss, p95 | `where` mag 0.0009, obscuration 0.0017 | umbral 0.0007, penumbral 0.0005 |
| contacts − Swiss, p50 / p95 | P1/P4 11 / 46 s; U1/U4 4 / 16 s; C1/C2 4 / 13 s (max 330 s in grazing geometry) | 1.3–1.8 / 4–6 s (max 162 s, grazing) |
| `where` point − Swiss (Swiss ΔT used for both) | 5.8 km p50, 9.8 km p95 | — |

Type differences with Swiss:
- **Solar: 33**, and NASA agrees with us in every one.
  - 16 × Swiss "total" where NASA and we say hybrid (H, H2, H3).
  - 14 × Swiss "hybrid" where NASA and we say annular.
  - 3 × Swiss "partial" where NASA and we say non-central A or T.
- **Lunar: 1**, at the umbral-magnitude 1.0000 boundary (ours 1.0002, Swiss 0.9998,
  NASA total).

The `where` offset of 5.8 km is on the Swiss side: we match NASA's greatest-eclipse
points to 0.3 km.

**Local solar (`sol_eclipse_when_loc`).** Ten cities (Paris, New York, Sydney,
Tokyo, Cairo, Santiago, Reykjavik, Dallas, Madrid, Longyearbyen), every eclipse
2000–2100, 388 cases:

| quantity | p50 | p95 |
|---|---|---|
| maximum | 0.37 s | 1.4 s |
| C1 | 0.25 s | 1.0 s |
| C4 | 0.57 s | 1.6 s |
| C2/C3 | 0.33 s | 0.76 s |
| sunrise/sunset | 1.2 s | 3.4 s |
| magnitude | | 0.0007 |
| obscuration | | 0.0003 |

- Same eclipse found in all 388 cases.
- One visibility-flag difference: a maximum at apparent altitude +0.005°.
- Remaining outliers are oracle artefacts on very small eclipses (magnitude
  0.0015–0.01): a C4 before the maximum, or 40–425 s off.

**Lunar local visibility (`lun_eclipse_how` at the contacts).** 60 eclipses × 10
cities with the API rule `attr[6] > −0.833`:
- **0 visibility differences.**
- Moon true altitude within 0.0009°, apparent altitude within 0.003°.

**End-to-end API check.** The API's own `core/moon_eclipse.py` and
`core/moon_forecast.py` were run with `swisseph` swapped for `swe_compat` (scratch
build of the extension, not installed). 3810 calls: 40 years of dates every 23
days × 5 cities + no location, plus `find_next_*`.
- **0 differences** in kind, type, blood-moon flag, visibility, local type,
  `days_*`, or the presence of any field.
- Timestamp differences:

  | field | p50 | p95 |
  |---|---|---|
  | `date` | 1 s | 6 s |
  | `starts_at` / `ends_at` (solar P1/P4: Swiss's less accurate contacts) | 2 s | 26 s |
  | `local_*` | ≤ 1 s | 2 s |

- Magnitude difference ≤ 0.001.
- `maximum_location` differs by 0.03° p50 (each engine at its own instant and
  Delta T).

## Timings

The Swiss figures are from the same Mac. The FrenchEphemeris figures come from the
Rust kernel (`cargo test --release … next_eclipse_timing`); through the Python
binding (dev build), the same calls take 22 / 32 / 197 / 4.3 / 4.3 µs.

| call (random start dates) | Swiss (pyswisseph) | FrenchEphemeris |
|---|---|---|
| next lunar eclipse (`lun_eclipse_when`) | 147 µs | **20 µs** |
| next solar eclipse (`sol_eclipse_when_glob`) | 255 µs | **32 µs** |
| next local solar eclipse, Paris (`sol_eclipse_when_loc`) | 870 µs | **197–215 µs** |
| `lun_eclipse_how` | 9.7 µs | 4.3 µs |
| `sol_eclipse_where` | 11.9 µs | 4.3 µs |

Searching the full 0–2400 range takes about 0.3 s for solar and 0.4 s for lunar.

## Open issues

- The pack ends shortly after 2400, so a "next eclipse" search there returns
  `OutOfRange`.
- Solar `tret[1]` (eclipse at local apparent noon), the umbra/penumbra limit
  coordinates `geopos[2..9]` and saros numbers were added later (see "Eclipse extras"
  below).
- Horizon altitudes ignore aberration and the height dip, at the 0.006° level. The
  oracle's two horizon rules are reproduced empirically, not derived.
- Swiss type flags differ in 33 solar cases (always against NASA). The API's solar
  `type` would change for those dates, e.g. hybrid where Swiss said total.
- `swe_compat`'s `ecltype` filter loops over our types; Swiss's exact filtering
  semantics for combined masks were not measured (the API always passes 0).

## Quick wins (2026-10-29): `lun_eclipse_when_loc`

`lunar_local.rs`. 12 sites (incl. Anchorage, Helsinki, Singapore, Cape Town), eclipses
1950–1960 and 2000–2025 forward, one backward search each: 132 cases.

Rules (documentation + black-box): the Moon is visible while its whole disc is above the
apparent horizon (lower-limb rise/set with refraction, `rise_set` with `DiscPoint::LowerLimb`
— the oracle's moonrise/moonset match it to 0.03 s); contacts while the Moon is not visible
are 0; `tret[0]` is the maximum if visible, else the moonrise/moonset bounding the visible
part nearest to it, and `attr` is evaluated there; the reported type is the deepest phase
visible (a total eclipse seen only in its partial phase is ECL_PARTIAL).

| vs Swiss | max |
|---|---|
| retflag (type and all 8 visibility bits) | identical, 132/132 |
| tret[0] | 6.9 s |
| U1/U4 (tret[2,3]) | 25 s (global contact differences, see above) |
| U2/U3, P1/P4 | 3.3 s, 6.6 s |
| moonrise/moonset (tret[8,9]) | 0.03 s |
| umbral/penumbral magnitude at tret[0] | 0.0008 |
| altitude at tret[0] | 0.023° (follows the tret[0] difference while the Moon rises) |

Saros numbers are not computed (−99999999, as in `lun_eclipse_how`).

## Lunar occultations (E2, 2026-10-29): `lun_occult_when_glob`, `lun_occult_when_loc`, `lun_occult_where`

Code: `crates/fe-core/src/occultations.rs`; bindings `crates/fe-py/src/eclipses.rs`; wrappers in
`swe_compat.py`; tests `crates/fe-core/tests/occultations.rs`, fixtures
`tests/fixtures/occultations/` (`gen_occultations.py`).

**Method.** This is the solar-eclipse machinery above, unchanged, with the occulted body in
place of the Sun (`eclipses::Target`/`Shape`): Besselian elements of the Moon's "shadow", the
exact WGS84 outline, contacts and local circumstances.

- **Stars** are point sources (radius 0). The direction comes from the catalogue (`stars.rs`,
  1176 objects) with space motion and parallax.
- **Planets** are discs of their IAU 2015 mean (volumetric) radius (Archinal et al. 2018), light
  time applied.
- The **lunar limb** is the mean radius k = 0.2725076 for every contact.
- **Candidates** are the Moon–body conjunctions in ecliptic longitude. Those with
  |Δβ| > 1.75° are rejected without further work. Stars with |β| > 7.5° are reported as never
  occulted (Swiss would search until the end of its ephemeris).
- **Local daylight flags** (`ECL_OCC_BEG/END_DAYLIGHT`) use the standard sunrise altitude of the
  Sun's centre, −0.8333°.
- **Body visibility** uses apparent altitude > 0, the same rule as eclipse contacts.
- `tret[5]` / `tret[6]` hold the body's rise/set between C1 and C4.
- For stars, C2 = C1 and C3 = C4, as Swiss documents.

**Return layouts.**

- `ECL_ONE_TRY` is read from `backwards`, as pyswisseph and the C API do
  (`backward |= SE_ECL_ONE_TRY`); it is also accepted in `flags`. Only the next conjunction is
  tested. When that conjunction gives no occultation, the flags are 0 and `tret[0]` is the time
  of that conjunction.
- `attr` in `when_loc` is the covered fraction, at most 1. `where` keeps the unclamped classic
  values; for example, obscuration is (Moon/planet diameter ratio)².

**Validation against independent predictions (IOTA).** IOTA's smooth-Moon, sea-level graze
limit lines come from the Occult program (D. Dunham, E. Riedel). The test places an observer on
each tabulated limit point and computes the closest approach of the limb.

| Event (IOTA file) | Points | Limb gap p50 / max | Graze time − IOTA p50 / max |
|---|---|---|---|
| Spica 2024-11-27 N limit (`20241127SpicaNlimit.txt`) | 159 | 0.056 / 0.074″ | 1.0 / 1.6 s |
| Spica 2024-11-27 S limit (`20241127SpicaSlimit.txt`) | 73 | 0.059 / 0.060″ | 3.0 / 3.9 s |
| Mars 2025-01-14 N, inner / outer (`20250114MarsNLimit.txt`) | 81 + 81 | 0.08 / 0.11″ and 0.06 / 0.17″ | 0.9 / 3.1 s and 0.8 / 2.9 s |
| Mars 2025-01-14 S, inner / outer (`20250114MarsSLimit.txt`) | 70 + 70 | 0.09 / 0.31″ and 0.13 / 0.18″ | 1.3 / 3.9 s and 1.5 / 3.9 s |

Sources: https://occultations.org/publications/rasc/2024/20241127Spica{N,S}limit.txt and
https://occultations.org/publications/rasc/2025/20250114Mars{N,S}Limit.txt.

1″ of limb gap is about 1.9 km across the path. IOTA times are rounded to 1 s and use IOTA's
Delta T.

**Swiss black box (pyswisseph 2.10.03).**

*Global search:*
- 1681 events: Regulus, Spica, Aldebaran and Antares 1950–2050, and Mercury–Pluto 1990–2030.
- Missed 0. "Spurious" 2: Jupiter 2001-05 and Neptune 2007-08, real grazing partial
  occultations (γ ≈ −1.27) that the oracle skips.
- Type/central flags differ in 1 case, a borderline Mercury partial/total.

| Global quantity | p50 | p95 | max |
|---|---|---|---|
| Maximum | 0.48 s | 4.7 s | 8.8 s |
| P1/P4 | 6.2 s | 23.7 s | |
| C1/C2 | 3.8 s | 11.9 s | |
| `where` point | 3.0 km | 9.6 km | |

The oracle's global contacts are less precise than its local ones; see the solar section above.

*Local search:*
- 549 events: 7 cities × 4 stars + Venus, Mars, Jupiter, Saturn, 2000–2040.
- Same events, and identical flags in all 549 cases, including the visibility and daylight bits.

| Local quantity | p50 | p95 | max |
|---|---|---|---|
| Maximum | 0.9 s | 2.9 s | |
| C1/C4 | 1.3 s | 4.1 s | 24 s (a graze) |
| C2/C3 | 1.9 s | 5.1 s | |
| Body rise/set inside the event | 11 s | 17 s | |

The oracle's rise/set rule is not published.

*Other modes:*
- `ECL_ONE_TRY`: the same found/not-found outcome in 239/240 cases; the exception is the Jupiter
  graze above. Found maxima agree within 1.5 s p50.
- When nothing is found, the oracle's `tret[0]` is about 1.0 day after the true conjunction from
  the 4th call on (an oracle artefact); ours is the conjunction.
- Backward search: 120/120 same events.

**Timings.**

| Search | Ours (dev profile, opt-level 1) | Swiss |
|---|---|---|
| Global | 86 µs | 400–950 µs |
| Local | 580 µs | 4–22 ms |

## Eclipse extras (E5, 2026-10-29): saros, local apparent noon, shadow limits

- **Saros series and member** (`eclipses::saros`, `attr[9..10]` of `sol_eclipse_where/how`,
  `sol_eclipse_when_loc`, `lun_eclipse_how`). Built from NASA's published series catalogues
  (series number, member count, date of first eclipse), with members 223 lunations apart:
  https://eclipse.gsfc.nasa.gov/SEsaros/SEsaroscat.html and
  https://eclipse.gsfc.nasa.gov/LEsaros/LEsaroscat.html. The table is generated by
  `tools/gen_saros.py` into `src/eclipses/saros.rs`.
  - Series number equal to NASA's canon column for all 5694 solar and 5801 lunar eclipses of
    0–2400.
  - Series and member identical to Swiss for 911/911 eclipses in 1900–2100.
- **Eclipse at local apparent noon** (`tret[1]` of `sol_eclipse_when_glob`,
  `eclipses::solar_local_noon`): the time between P1 and P4 when the body is on the meridian at
  the point of greatest eclipse of that instant.
  - The upper transit is used if there is one. Otherwise the lower transit (midnight Sun) is
    used, which the oracle also reports.
  - Agreement with Swiss: 0.11 s p50, 0.38 s p95. Presence matches in 446/453 eclipses; the
    other 7 are polar passes where the point sweeps past the pole between scan steps.
- **Umbra/penumbra limits** (`geopos[2..9]` of `sol_eclipse_where` and `lun_occult_where`): the
  points of the shadow outline (radius l − ζ·tan f at the site) where the shadow's motion
  *relative to the rotating Earth* is tangent to the outline (Explanatory Supplement 3rd ed.
  §11.3). The oracle does not implement these and returns zeros.
  - Against NASA's path tables (2017 Aug 21 T, 2023 Oct 14 A, 2024 Apr 08 T; limits every
    2 min), the umbral limits agree within 0.36 / 0.46 km p50 (N / S) and 3.5 km p95.
  - The maximum, 14 km, is at the sunrise/sunset ends of the path.
  - Using the shadow's motion without the Earth's rotation gave 24 km; that is why the rotation
    is included.

## Visibility at a place (2026-10-29, drop-in comparison fixes)

- **`lun_eclipse_how` / `sol_eclipse_how` flags.** As the classic functions do, nothing is
  reported while the body's centre is not above the apparent horizon: refracted altitude ≤ 0
  (Sæmundsson, at the standard pressure for the site height). The switch was measured
  black-box at sea level and at 2000 m. The magnitudes are still given for the Moon below the
  horizon.
  - With no eclipse at all, the flags are 0 and the eclipse attributes are zero.
  - A visible solar eclipse adds `ECL_VISIBLE` and `ECL_CENTRAL` or `ECL_NONCENTRAL`, depending
    on whether the shadow axis meets the Earth at that instant. Core helpers:
    `eclipses::lunar_seen_kind`, `solar_seen_kind`, and `SolarAttrs::axis_on_earth`.
- **`lun_eclipse_when_loc` at polar sites.** The Moon's visibility inside [P1, P4] is now
  sampled every 5 minutes and refined. Before, it was inferred from the next rising and
  setting, which wrongly marked the Moon as up when it rose in the window and did not set
  again (or the reverse).
- **`sol_eclipse_when_loc`** sets `ECL_NONCENTRAL` when the shadow axis misses the Earth at the
  reported maximum. It never sets `ECL_CENTRAL`, as measured.
- **Places harness** (`harness/compare_places.py`, 40 places, 1040 cases, before → after):

  | Check | Before | After |
  |---|---|---|
  | Eclipse visible at place | 605 | 0 |
  | Local eclipse flags | 375 | 15 |
  | Next local eclipse is another event | 43 | 1 |

  - The 15 remaining flag differences are single contact-visibility bits for contacts within
    seconds of moonrise or moonset.
  - The one remaining event (Utqiagvik 1919-10-07) is an eclipse whose end is seen with the
    Moon up to 4° high. The oracle skips it.
- **Regression tests:** `tests/eclipse_visibility.rs`.
  - Public: physical consistency.
  - Private: oracle `eclipses/visibility.json`, 700 `*_how` cases with 0 differences, and 336
    local searches with 1 different event (Longyearbyen 2005-10-17, Moon 2.5° up at P4,
    skipped by the oracle) and 0 NONCENTRAL differences.
