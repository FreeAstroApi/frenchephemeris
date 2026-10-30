# Astrology geometry — definitions, Swiss behaviour, accuracy

Workstream: `crates/fe-core/src/{houses,ayanamsa,points}.rs` (+ `houses/*.rs`, `ayanamsa/{vondrak,star,galactic,orig_precession}.rs`).
Oracle: pyswisseph 2.10.3.2 (`swe.version` = 2.10.03), Swiss files from
the Swiss Ephemeris data directory (`SE_EPHE_PATH`), used strictly as a black box (public
Python API, plus the public C entry point `swe_houses_armc` via ctypes to read
the cusps Swiss fills in when it returns ERR). No Swiss source was consulted.

Fixtures are regenerated with the scripts next to them:

```
.venv/bin/python crates/fe-core/tests/fixtures/houses/gen_houses.py
.venv/bin/python crates/fe-core/tests/fixtures/houses/gen_houses_ext.py
.venv/bin/python crates/fe-core/tests/fixtures/ayanamsa/gen_ayanamsa.py
.venv/bin/python crates/fe-core/tests/fixtures/points/gen_points.py
cargo test -p fe-core -- --nocapture   # prints the max-error tables below
```

The ayanamsa and points generators also *assert* the black-box relations
stated below (they fail if Swiss behaves differently).

---

## 1. Houses (`houses.rs`)

Inputs: ARMC (= local apparent sidereal time × 15), geographic latitude φ,
obliquity ε; degrees. Output cusps 1..12 and the 8 Swiss `ascmc` angles
(`Houses::ascmc()` returns them in Swiss order).

### Definitions used

Notation: `asc(θ, p)` = ecliptic longitude rising on an oblique horizon of pole
height `p` when the meridian has RA θ:
`asc = atan2(cos θ, −(sin θ cos ε + tan p sin ε))`.

| Item | Definition |
|---|---|
| MC | `atan2(sin θ, cos θ cos ε)` |
| Asc | `asc(θ, φ)`, then the hemisphere rule below |
| Vertex | `asc(θ+180, 90−φ)`, forced into the western half `(MC−180, MC)` |
| Equatorial Asc (East point) | `asc(θ, 0)` |
| Co-Asc (W. Koch) | `asc(θ, −φ)` |
| Co-Asc (M. Munkasey) | `asc(θ, 90−φ)` |
| Polar Asc (M. Munkasey) | `asc(θ, −φ) + 180` |
| Placidus | cusp 11/12: `RA = θ + f·DSA(δ)`, f = 1/3, 2/3; cusp 2/3: `RA = θ + 180 − f·NSA(δ)`, f = 2/3, 1/3; δ the declination of the ecliptic point with that RA. Solved by damped Newton on RA (analytic derivative), residual < 1e-13 rad |
| Koch | D = diurnal semi-arc of the MC degree; cusps 11, 12, 2, 3 = `asc(θ − 2D/3)`, `asc(θ − D/3)`, `asc(θ + D/3)`, `asc(θ + 2D/3)` at latitude φ |
| Regiomontanus | equator trisected: cusp at `asc(θ + H − 90, p)`, H = 30, 60, 120, 150, `tan p = tan φ · sin H` |
| Campanus | prime vertical trisected: azimuth-angle k = 30, 60, 120, 150; `sin p = sin φ sin k`, equator crossing `tan H = cos φ tan k`, cusp `asc(θ + H − 90, p)` |
| Topocentric (Polich–Page) | `asc(θ + H − 90, p)`, H = 30, 60, 120, 150, `tan p = tan φ · n/3` (n = 1, 2, 2, 1) |
| Porphyry | trisection in longitude of the MC→Asc and Asc→IC arcs |
| Equal | Asc + 30k |
| Whole sign | `floor(Asc/30)·30 + 30k` |

Sources: spherical-astronomy definitions of the systems as published (e.g.
M. Munkasey, *A Concise Guide to Astrological House Systems*; Polich & Page,
*The Topocentric System of Houses*, 1964; Meeus, *Astronomical Algorithms*,
ch. 13 for coordinate transforms). The Koch below-horizon rule (cusps 2/3 use
the same D as 11/12) and all edge rules were confirmed against the oracle.

### Swiss edge behaviour (measured) and how it is reproduced

1. **Polar rule (Placidus, Koch).** Swiss refuses P and K iff
   `|φ| >= 90 − ε` (bisection on φ to 1e-14° for several ε and ARMC: the
   threshold is exactly 90 − ε, independent of ARMC; `|φ| = 90 − ε` itself is
   refused). We return `Error::PolarHouses`.
   - The C library returns ERR **and** fills cusps and all 8 ascmc with
     **Porphyry** (verified identical through `swe_houses_armc` via ctypes;
     the message from `houses_ex2` is "within polar circle, switched to
     Porphyry").
   - **pyswisseph raises `swe.Error`** for `houses`, `houses_ex`,
     `houses_armc`, `houses_ex2`, so FreeAstroAPI today *fails* the request for
     P/K charts above the polar circle (no try/except around
     `swe.houses_ex` in `core/astrology.py`, `core/vedic_astrology.py`,
     `core/vedic_kp_v2.py`).
   - `houses_swiss_compat()` reproduces the C fallback (returns
     `(Porphyry houses, true)`). The adapter must pick: raise (today's API
     behaviour) or fall back.
   - `house_pos` does **not** refuse P/K inside the polar circle (see below).
2. **Latitude clamp.** Swiss uses `|φ| = 90 − 1e-10` for |φ| beyond that
   (inferred from the Asc at φ = ±90, which matches φ = ±(90 − 1e-10) to
   1e-12). Reproduced (`MAX_LAT_DEG`).
3. **Asc hemisphere rule.** The Asc is always in `(MC, MC+180)`; inside the
   polar circle the raw rising point can lie west of the meridian and is then
   replaced by its opposite (holds for all 8 systems and all grid points).
4. **R/T/C flip.** When rule 3 flips the Asc, Regiomontanus, Topocentric and
   Campanus additionally turn every cusp except 1 and 7 by 180°, **and
   `ascmc[1]` (MC) becomes the IC** for those systems (Porphyry/Equal/Whole
   keep the MC and only use the flipped Asc). Reproduced.
5. **Vertex** always in the western half: if `asc(θ+180, 90−φ)` falls in
   `(MC, MC+180)` exclusive, it is turned by 180. At φ = 0 exactly the raw
   vertex is MC ± 180 and is kept (Swiss gives 180 for ARMC 0).
6. **Exact angles.** Swiss returns exact 90/180/0 for ARMC multiples of 90°;
   we use degree-exact sin/cos and axis-exact atan2 so Whole-sign cusps do
   not jump a sign by one ulp.
7. **Degenerate ecliptic = horizon** (`|φ| = 90 − ε` exactly and ARMC = 90/270):
   Asc := MC + 90 (matches Swiss). In a ~1e-4° neighbourhood of that
   singularity Swiss snaps the Asc to exactly 90/270, we give the
   mathematically continuous value (differences ≤ 5e-5°); excluded from the
   1e-7 test (measure-zero, documented).
8. **Whole sign at a sign boundary**: when the Asc is within ~1e-9° of a
   multiple of 30°, Swiss' Asc rounding (≈4e-10° below ours at φ = ±1e-9)
   can select the other sign. Test skips rows with Asc within 1e-8° of a
   boundary.
9. **Placidus precision.** Swiss' own Placidus is not fully converged:
   ≤ 5.4e-8° for |φ| ≤ 60°, up to 1.04e-6° just inside the polar circle. Our
   cusps satisfy the semi-arc equation to < 1e-13 rad; the residual vs Swiss
   is Swiss' iteration error (its trajectory/stop rule could not be
   reproduced black-box; not a definition difference).

### `house_position` (= `swe.house_pos`)

| System | Definition reproduced |
|---|---|
| P | fraction of the point's own diurnal / nocturnal semi-arc (from RA, δ incl. ecliptic latitude); circumpolar points use clamped semi-arcs (AD = ±90). Works inside the polar circle. |
| K | point's rising (east of meridian) / setting (west) sidereal time relative to the trisected semi-arc D of the MC degree; **returns 0.0** (like Swiss) when that time is outside ±D (points near the meridian with latitude) |
| R | angle ψ of the point around the north–south horizon axis, mapped to the equator `tan H = cos φ tan ψ` |
| C | ψ directly (uniform on the prime vertical) |
| T | solve `RA = θ + H + asin(tan p(H) tan δ)`, `tan p = tan φ · min(H,180−H)/90` (west half by mirror) |
| E/W/O | linear in longitude between the cusps of `houses()` |

Swiss quirks **not** reproduced (provenance rule: no oracle-derived
numbers): Swiss adds **+1 mas** to the "house degree" for P, K, R, C (to the
longitude for E, W, O), returns exactly the cusp number for points within
1 mas of a cusp, and skips the offset for circumpolar P points and for T.
Our `house_position` is pure geometry (a point on a cusp gets that cusp's
number). Difference: ≤ 1 mas (2.8e-7 house-degrees = 9.3e-9 house units);
for Porphyry the 1 mas is in longitude, i.e. up to 1.9e-6 house-degrees in
the narrow houses near the poles. An integer house number can differ only
for a point within 1 mas (2.8e-7°) of a cusp.

FreeAstroAPI calls `swe.house_pos(ascmc[2], lat, 23.44, (lon, 0.0), hsys)`
(fixed ε = 23.44, latitude 0, sidereal longitude when sidereal) — the adapter
must pass the same arguments to reproduce today's house numbers.

### Accuracy vs Swiss (12 576 house rows, 4 320 house_pos rows)

| System | |φ| ≤ 60 | 60 < |φ| < 90−ε | polar region / fallback |
|---|---|---|---|
| O, E, W, R, C, T, K | ≤ 4.1e-10° | ≤ 1.0e-10° | ≤ 8.3e-10° (K fallback = Porphyry, 1.0e-10°) |
| P | ≤ 5.4e-8° | ≤ 1.04e-6° (Swiss iteration error) | fallback ≤ 1.0e-10° |
| error/fallback category | exact on all rows |||

house_pos (in house-degrees, 30 per house): P, K, R, C, E, W ≤ 2.8e-7 and
O ≤ 1.9e-6 (exactly Swiss' 1-mas shift, which we do not apply; with the
shift added the geometry agrees to 3.3e-12);
Koch zero/undefined category exact; T ≤ 5.5e-6 for non-circumpolar points;
T circumpolar points: Swiss picks a different root (up to 87 house-degrees) —
not reproduced (open issue, Topocentric house_pos only).

Speed (release, Apple M-series): R/T/C/K ≈ 0.27–0.37 µs, O/E/W ≈ 0.16 µs,
P ≈ 1.0–1.2 µs per call (all 12 cusps + 8 angles); house_position(P) ≈ 65 ns.

---

### Additional systems, speeds, Gauquelin sectors, names (`houses/*.rs`)

Every remaining code of the Swiss documentation (swisseph.htm §6.2, swephprg.htm §15):
`A V D N B X F M H U Y L Q S I i G`, plus `J` (Savard-A, listed by the
oracle's `house_name`, not in its documents). Codes are case-insensitive except
`I` (Sunshine, variant I) vs `i` (Sunshine, Makransky's solution).

API: `houses_ex(sys, armc, φ, ε, sun_dec)` (12-cusp systems; `houses` = no
Sun declination), `house_cusps` / `house_cusps_swiss_compat` (any system,
`Vec` of 12 or 36 cusps), `gauquelin_cusps`, `house_cusps_speeds(_swiss_compat)`
(`houses_armc_ex2` layout), `house_position_ex` (all systems, Sun declination for
I/i), `house_name`, `Engine::house_cusps_ut` (`houses_ex2`: date, Sunshine Sun
declination = apparent geocentric declination of date, sidereal, speeds) and
`Engine::gauquelin_sector`. Python: `houses_armc(…, ascmc9, swiss_fallback)`,
`houses_armc_ex2`, `houses_ex2(engine, …)`, `house_pos(…, sun_dec)`, `house_name`,
`gauquelin_sector(engine, …)`; `swe_compat` mirrors pyswisseph.

Notation as above; `pc(k)` = the ecliptic point on the circle through the north
and south points of the horizon that meets the equator `k` degrees east of the
meridian: `asc(θ + k − 90, p)`, `tan p = tan φ · sin k` (Regiomontanus' circle
family); the circle through a point at hour-angle offset `d = RA − θ` and
declination δ has `tan k = sin d / (cos d + tan φ tan δ)`.

| Code | Definition implemented | Source |
|---|---|---|
| A | = E (Asc + 30k) | Swiss doc §6.2.5 |
| V | Vehlow: Asc − 15 + 30k (Asc in mid-house 1) | Swiss doc §6.2.5 |
| D | Equal from MC: MC + 90 + 30k (cusp 10 = MC) | Swiss doc §6.2.5 |
| N | 0° Aries + 30k | Swiss doc §6.2.5 |
| B | Alcabitius: Asc's diurnal/nocturnal semi-arcs trisected on the equator from the meridian; cusp = ecliptic point with that RA (hour circles); 5/6/8/9 opposite | Rhetorius / al-Qabisi as described in Swiss doc §6.2.11 |
| X | Meridian (axial rotation): RA = θ + 90 + 30k → ecliptic point with that RA | Swiss doc §6.2.7 |
| F | Carter poli-equatorial: RA = RA(Asc) + 30k → ecliptic | C. E. O. Carter, *Essays on the Foundations of Astrology* (1947) pp. 158–159 |
| M | Morinus: equator point RA = θ + 90 + 30k projected along circles through the ecliptic poles (`tan λ = tan α cos ε`) | Swiss doc §6.2.8 |
| H | Horizon: verticals at azimuth 30k from the south point (cusp 10) through east (cusp 1); vector intersection `(Z × d) × n_ecl` | Swiss doc §6.2.9 |
| U | Krusinski-Pisa-Goelzer: great circle Asc–zenith divided in 12 from the Asc (zenith = cusp 10); points projected to the ecliptic along hour circles | B. Krusinski (1995) as quoted in Swiss doc §6.2.13; G. Goelzer, *Der Ich-Kosmos* (1995) |
| Y | APC (Knegt): parallel of the Asc divided in six above and six below the horizon; cusps `pc(k)` through the division points | Swiss doc §6.2.14 |
| L | Pullen SD: middle house of a quadrant q = 30 + (q − 90)/2, outer ones 30 + (q − 90)/4; negative middle house → 0 (outer = q/2) | W. Pullen (Astrolog, 1994), Jan-2016 rule, both quoted in Swiss doc §6.2.6 |
| Q | Pullen SR: small quadrant s = rx + x + rx, large 180 − s = r³x + r⁴x + r³x, solved for r (bisection) | W. Pullen, 24 Jan 2016, quoted in Swiss doc §6.2.6 |
| S | Sripati: midpoints of consecutive Porphyry cusps | Swiss doc §6.2.6 |
| I | Sunshine: the Sun's diurnal/nocturnal semi-arcs trisected, cusps `pc(k)` through the trisection points on the Sun's parallel; circumpolar Sun → arcs 180/0 and the empty half's cusps on the meridian | B. Makransky, *Primary Directions* (1988), as described in Swiss doc §6.2.15 |
| i | same construction, Makransky's solution (see deltas) | idem |
| G | Gauquelin: 36 sectors, clockwise, Asc = 1, MC = 10, Desc = 19, IC = 28; each semi-arc divided in 9 (Placidus construction with fractions k/9) | M. Gauquelin; Swiss doc §6.2.12 |
| J | Savard-A: parallels of declination divide the prime vertical between the east/west points and the zenith (nadir) in thirds of the zenith's declination φ — the prime-vertical point at altitude h has sin δ = sin φ sin h, so cusp 11/12 lie at sin h = sin(2φ/3)/sin φ, sin(φ/3)/sin φ (limit 2/3, 1/3 at φ = 0); house circles through the north and south points of the horizon (Campanus construction with those points); 5/6/8/9 opposite. Asc flip: Campanus rule (all cusps but 1/7 turned, MC → IC) | J. J. G. Savard, "Astrological House Systems", quadibloc.com/other/as01.htm ("Albategnus"; © 2012); name "Savard-A" (the system is not in al-Battani's works): Astrolog 7.30 change notes, swisseph@groups.io thread "What is the origin of this…" (2021) |

Measured edge rules (pyswisseph 2.10.03, `swe_houses_armc_ex2` via ctypes):

1. **Errors/fallback.** Only P, K, G (|φ| ≥ 90 − ε) and `i` (|tan φ tan δ☉| ≥ 1)
   return ERR; all other new systems are defined everywhere (B clamps the
   Asc's semi-arc). Fallback = Porphyry cusps + Porphyry ascmc ("within polar
   circle, switched to Porphyry"); for G the C library fills sectors 1–12 with
   the 12 Porphyry cusps and 13–36 with 0. Reproduced by `*_swiss_compat`.
   `I` additionally returns ERR in a ≈1e-8-wide, ARMC-dependent band around
   |tan φ tan δ☉| = 1 (numerical); not reproduced (we return the cusps).
2. **Asc flip.** All systems use the hemisphere-ruled Asc. Y and I behave like
   R/T/C when it flips: MC ↔ IC in `ascmc[1]`; I turns every cusp but 1/7 by
   180°, Y relabels the raw construction's houses by 6 (its cusps 2/3 and 8/9
   are not opposite). `i` keeps the raw MC and cusps. All others keep the MC.
3. **H orientation.** Northern hemisphere: azimuths from the south point,
   intersection along `(Z×d)×n_ecl`. For φ ≤ 0 (the equator included) the
   oracle counts from the north point and turns the result by 180°.
4. **Sunshine Sun declination.** `houses_armc` takes it in `ascmc[9]`
   (pyswisseph: `ascmc9`, default 0); dated `houses*` use the apparent
   geocentric declination of date (verified to 4e-14°; the TRUEPOS one differs).
   The oracle's `house_pos('I'/'i')` reads it from the process-global state of
   its last houses call; ours takes it as an argument (`swe_compat` keeps a
   per-context value set by `houses_armc(…, ascmc9)`).

Accuracy vs Swiss (14 076 cusp rows incl. random grid φ ∈ (−89.9, 89.9),
ε 22.5–24.5, all ARMC; 2 000 speed rows; 7 905 house_pos rows; 240 Gauquelin
sectors; `cargo test -p fe-core --test houses_ext -- --nocapture`):

| System | cusps + ascmc max Δ | house_pos max Δ (house-deg) |
|---|---|---|
| A V D N L Q S M X F | ≤ 5.0e-10° | A V D 2.8e-7, N L Q M X F ≤ 1.5e-9, S ≤ 2.2e-6 (1 mas of longitude in narrow polar houses) |
| B | ≤ 1.5e-8° (φ = −89.9), else ≤ 4.2e-10° | 1.5e-9 |
| H | ≤ 9.0e-9° | 2.8e-7 |
| U | ≤ 4.2e-10° | 2.8e-7 |
| Y | ≤ 5.0e-10° | 2.8e-7 (1.1e-6 on the equator, φ = 0 exactly) |
| I | ≤ 6.5e-9° | 2.8e-7 |
| i | ≤ 4.2e-10° outside the cases below | 2.8e-7 |
| G | ≤ 1.2e-7° (|φ| ≤ 60), ≤ 3.2e-6° (60–66.5: the oracle's Placidus iteration error); fallback exact | 8.4e-7 (= the oracle's 1 mas Placidus offset, in sector units ×3) |
| J | ≤ 4.2e-10° (also at the exact poles and φ = 0; no error or fallback anywhere) | 1.5e-9 (linear in the Campanus position between the cusps' Campanus positions — the oracle's form, measured; no 1-mas offset); speeds vs the oracle ≤ 5.5e-8 relative, no speed defect |
| error / fallback category | exact on every row | |

The 2.8e-7 house-degree class is the oracle's +1 mas offset/snap (see above),
which we do not apply; N, L, Q, B, M, X, F have no such offset in the oracle.
`house_pos` forms (chosen to agree with the oracle, which documents Y L Q I S
as "simplified" but actually computes): linear in longitude between cusps for
A V D N L Q; Sripati = Porphyry position + ½; geometric for the others —
X: RA − θ; F: RA − RA(Asc); M: equator point of the ecliptic longitude (ecliptic
latitude ignored); B: linear in RA between the cusp hour circles (west of the
meridian with the Desc's semi-arcs); U: angle along the Asc–zenith circle where
the point's hour circle meets it; H: azimuth from the south point (both
hemispheres); G: 1 + (90 − Placidus house degree)/10 (clockwise); Y/I/i: the
point's `pc` circle cuts the reference parallel (Asc, resp. Sun) at hour-angle
offset d, house degree = d as a fraction of that parallel's semi-arcs.

Oracle behaviours **not** reproduced (documented, measure-zero or defects):

- `i` (Makransky): when 2/3 of a semi-arc of the Sun exceeds 90° (|φ| ≳ 60,
  |tan φ tan δ☉| > ≈ 0.7) the oracle's cusps 8/12 (diurnal) or 2/6 (nocturnal)
  differ from the geometric construction (by up to tens of degrees, 76 of 1 088
  `i` rows). Searched for Makransky's own formulas (2026-10-29): his public
  description (B. Makransky, "The Sunshine House System", astro.nu, 2017;
  essay "Intuition in Astrology: The Sunshine House System" in *Topics in
  Astrology*) defines only the geometric construction implemented here (Sun's
  diurnal/nocturnal arcs in six, house circles through the N/S horizon points);
  the calculation recipe is in *Primary Directions: A Primer of Calculation*
  (1988), whose author-sanctioned free download (dearbrutus.com) is gone and
  which is otherwise only a copyrighted book (2022 paperback) or unauthorised
  uploads, so it was not used. The 76 rows stay open until a licensed copy is
  consulted; the oracle's deviation looks like a hand-calculation formula valid
  only for trisection points within 90° of the meridian. Also a
  singular branch in the oracle's `i` for δ☉ = 0 exactly and for φ = 0 exactly
  (cusps 2/8, 5/11 change half-plane; any |δ☉| ≥ 1e-9 agrees).
- `I` cusps are NaN in the oracle for a few circumpolar-Sun rows (7 rows).
- Exact poles (|φ| = 90, clamped to 90 − 1e-10): tan φ ≈ 6e11 amplifies
  1e-12 differences of the Asc declination: B (17°), U and I (≈ 2e-3°), Y (180°),
  Q (2.8e-7°). Only these rows; |φ| = 89.9 agrees as tabled.
- co-Asc (Munkasey) at φ = 0 exactly (pole height 90°): the oracle returns 0 or
  180 depending on the house system; H at φ = 0 with ARMC = 0/180 (equinox at
  the zenith, every vertical meets the ecliptic there).
- `house_pos`: S — the oracle returns 1.0 for every point whose Porphyry
  position + ½ is ≥ 12 (Sripati house 12 and the first half of house 1); we
  return the value. H — the oracle numbers plain azimuths from the south point
  in both hemispheres, which is the mirror of its own cusps for φ ≤ 0 and when
  the zenith is south of the ecliptic; its cusp snap then gives a point exactly
  on a cusp that cusp's number (15 rows), we give the azimuth position.
  U — the oracle's value is ill-conditioned on the equator (moves 3e-5 house
  between φ = 0 and 1e-9). I/i with a circumpolar Sun — the oracle returns
  degenerate values (4.0, 10.0 or arbitrary); we return 4 (no nocturnal arc) /
  10 (no diurnal arc).

**Speeds (C2).** The oracle's speeds are ∂/∂ARMC × the sidereal rate with ε
(and δ☉) fixed; its ARMC speed is 360.98564733°/day. We use 360 × 1.002737909350795
(Aoki et al. 1982) = 360.98564737°/day (Δ = 3e-8°/day) and central differences
of ±1e-3° ARMC (one-sided across a cusp jump). Against the oracle's own cusps
differentiated the same way, every system agrees to ≤ 2e-8 relative (P/G ≤ 1e-6:
the oracle's Placidus iteration noise). Against the oracle's returned speeds,
where those are consistent with its own cusps: ≤ 9e-6 relative. Oracle speed
defects (its returned speed ≠ the derivative of its own cusps; we return the
derivative): Porphyry (and every Porphyry fallback) cusps 5/6/11/12 (always);
W and N cusps 1/4/7/10 (the oracle returns the Asc/MC speeds for these fixed
cusps); U cusps other than 1/4/7/10 (the oracle returns 0); P, K, G inner
cusps in about half of the rows (wrong by up to 30 % near the polar circle).

**Gauquelin sectors (C3)** (`Engine::gauquelin_sector`, methods 0–5):
0/1 = `house_position('G')` of the body's apparent ecliptic place of date
(topocentric when the flags say so) with ARMC = LAST and true ε, latitude
dropped for 1; 2–5 = previous and next rise/set of the disc centre (2, 3) or
upper limb (4, 5) without (2, 4) or with (3, 5) refraction for the given
atmosphere, through `events::rise_set`: above the horizon
s = 1 + 18 (t − rise)/(set − rise), below s = 19 + 18 (t − set)/(rise − set); a
body that does not rise and set in the search window is an error, as in the
oracle. Agreement 1950–2050 (240 rows, all bodies, |φ| ≤ 80): methods 0/1
≤ 6e-6 sector, 2–5 ≤ 1.3e-5 sector (rise/set instants ~0.1 s); the
circumpolar error category is exact. Fixed stars are not supported (stars
workstream). Dated `houses_ex2` (I, i, G, P, Y, H) agree to 2.2e-6° (ARMC/ε
of date; `i` cusps of the unreproduced cases excluded).

**`house_name` (C4)**: the oracle's labels for all codes ('' for unknown codes).

---

## 2. Ayanamsa (`ayanamsa.rs`)

### Swiss behaviour (measured; asserted in `gen_ayanamsa.py`)

- `get_ayanamsa(jd_tt)` = `get_ayanamsa_ex(jd_tt, FLG_NONUT)` = **mean**
  ayanamsa (no nutation). `get_ayanamsa_ut(jd_ut)` = the same at
  `jd_ut + ΔT` — this is what FreeAstroAPI calls.
- `get_ayanamsa_ex(jd_tt, 0)` = **mean + Δψ** (true ayanamsa), to 1e-12°.
  Other flags (TRUEPOS, NOABERR, SWIEPH/MOSEPH) do not change it.
- **Sidereal positions**: `calc(FLG_SIDEREAL | f)` longitude =
  `calc(f)` longitude − `get_ayanamsa_ex(jd_tt, f)` **exactly** (≤ 1e-13°;
  sampled over bodies incl. nodes/Lilith, modes and f ∈ {0, NONUT, TRUEPOS,
  TRUEPOS|NONUT}). Latitude unchanged; speed = tropical speed − ayanamsa rate.
  No re-projection onto the ecliptic of t0 is involved. Hence:
  - default flags: `sid = apparent tropical (true equinox) − (mean + Δψ)` =
    `tropical referred to the mean equinox − mean ayanamsa`;
  - with `FLG_NONUT`: `sid = tropical(mean equinox) − mean ayanamsa` — the same
    number (Δψ cancels), to the extent the non-nutated tropical differs only
    by Δψ in longitude.
- **Sidereal houses** (`houses_ex(..., FLG_SIDEREAL)`): tropical houses (ARMC
  from GAST, true obliquity, both at `jd_ut + ΔT` for ε) minus the **true**
  ayanamsa for every cusp and every ascmc entry **except ARMC**. With
  `FLG_SIDEREAL | FLG_NONUT` (used by `core/vedic_strength_v2.py`): ARMC from
  GMST (= ARMC − equation of the equinoxes), **mean** obliquity, minus the
  **mean** ayanamsa. ARMC = `swe.sidtime(jd_ut)·15 + lon`.
- FreeAstroAPI maps both `"pushya"` and `"true_pushya"` to
  `SIDM_TRUE_PUSHYA` (`core/vedic_astrology.py`); `Ayanamsa::Pushya` is
  therefore an alias of `TruePushya`. Other modes used: LAHIRI, RAMAN,
  KRISHNAMURTI ("kp"), FAGAN_BRADLEY, YUKTESHWAR.

### Definition implemented (fixed-epoch modes)

The sidereal zero point is fixed on the mean ecliptic of the reference epoch
t0, A0 behind the mean equinox of t0; the ayanamsa at t is measured **on the
ecliptic of t0**:

`ayanamsa(t) = A0 − λ_{ecl(t0)}(mean equinox of t)`

with the **Vondrák, Capitaine & Wallace (2011)** long-term ecliptic and equator
poles (A&A 534, A22; coefficients verified against ERFA `ltpecl/ltpequ` to
1e-16). Tested alternatives: same method with IAU 2006 precession (0.34" off
at year 0), and "precess the zero point to date" (1.7" off) — rejected.

| Mode | t0 (JD TT) | published value (digits as published) | interpretation (Swiss Ephemeris documentation, "Ayanamshas") | A0 used (mean) |
|---|---|---|---|---|
| Lahiri | 2435553.5 (1956-03-21 0h TDT) | 23°15′00.658″ (Indian Astronomical Ephemeris 1985 refinement of the ICRC 23°15′00″, CSIR 1955) | value is TRUE (with nutation); IAU 1976 precession offset −0.13036″ | value + 0.13036″ − Δψ_IAU2000B(t0) |
| Raman | 2415020.0 (J1900.0) | 21.01444° (B. V. Raman) | Newcomb precession offset +0.82800″ | value − 0.828″ |
| Krishnamurti | 2415020.0 (J1900.0) | 22.363889° (K. S. Krishnamurti, table in *KP Reader I*) | Newcomb offset +0.82800″ | value − 0.828″ |
| Fagan–Bradley | 2433282.42345905 (B1950.0) | 24°02′31.36″ (SVP 335°57′28.64″, Fagan & Bradley) | Newcomb offset +0.41256″ | value − 0.41256″ |
| Yukteshwar | 2415020.0 (J1900.0) | 21.082222° | no offset | value |

The interpretation (true vs mean, the precession-model offsets) is taken
from the Swiss Ephemeris public documentation (swisseph.htm, ayanamsha
section; documentation, not code). Black-box comparison was used only to
choose among published alternatives: J1900.0 vs 1 Jan 1900 0h (Raman/KP),
B1950.0 vs 1 Jan 1950 0h (Fagan–Bradley), and 23°15′00″ vs 23°15′00.658″
(Lahiri). The hypothesis "offset = Δψ(t0)" was tested and rejected (Δψ is
17.3″ at 1900 and −3.3″ at 1950; the offsets are 0.828″ and 0.413″).

**True Pushya** (P. V. R. Narasimha Rao): δ Cancri at sidereal 106° (16° Cancer).
Swiss: `get_ayanamsa_ex(jd, f) = fixstar2(",deCnc", jd, f).lon − 106` exactly
(6e-14°), where the star longitude is **apparent** (annual parallax, light
deflection, annual aberration) referred to the mean equinox (f = NONUT → mean
ayanamsa) or the true equinox (f = 0 → true ayanamsa). `true_pushya_deg(lon)`
implements exactly that. The δ Cnc catalogue is **Hipparcos, the New
Reduction** (van Leeuwen 2007, A&A 474, 653; HIP 42911, epoch J1991.25:
α 2.2893709279 rad, δ 0.3168621473 rad, μα* −17.67, μδ −229.26 mas/yr,
ϖ 24.98 mas, CDS catalogue I/311) with RV +16.39 km/s (Massarotti et al.
2008), linear space motion; the oracle's catalogue is the Hipparcos solution
as well.

`ayanamsa_mean_deg(TruePushya | Pushya)` uses a built-in apparent place
(`ayanamsa/star.rs`): catalogue + IAU 2006 frame bias + Meeus ch. 25 Earth +
deflection + first-order aberration + Vondrák ecliptic of date.
`star::delta_cnc_apparent_mean_from_states()` takes exact barycentric Earth
(pos, vel) and Sun states from the ephemeris instead (preferred in the
adapter).

### Published model vs Swiss (1 101 dates; dense 1800–2100 and a conjunction window)

ours − Swiss, mean ayanamsa (identical for the true ayanamsa, since both add Δψ):

| Mode | 0–1000 | 1000–1800 | 1800–2100 | 2100–2400 |
|---|---|---|---|---|
| Lahiri | −0.008″ (2.2e-6°) | same | same | same |
| Raman | −0.00005″ (1.3e-8°) | same | same | same |
| Krishnamurti | −0.00005″ (1.3e-8°) | same | same | same |
| Fagan–Bradley | +0.00003″ (7.4e-9°) | same | same | same |
| Yukteshwar | < 1e-9° | same | same | same |
| True Pushya, built-in Earth | ±0.019″ (5.2e-6°) | ±0.018″ (5.1e-6°) | ±0.020″ (5.4e-6°) | ±0.018″ (5.1e-6°) |
| True Pushya, JPL/oracle states | ≤ 0.004″ (1.2e-6°) | ≤ 0.002″ (4.3e-7°) | ≤ 0.001″ (1.7e-7°) | ≤ 0.001″ (1.5e-7°) |
| True Pushya, definition (star lon − 106°) | 6e-14° | | | |

The remaining Lahiri −0.008″ is unexplained (nutation model choice changes it
by only 0.0004″); it is constant and invisible.

(True Pushya rows exclude the few hours a year when δ Cnc is behind the
solar disc, where the oracle additionally attenuates the deflection:
+1.0–1.3e-4°.) Sidereal longitudes and sidereal house cusps inherit these
deltas 1:1.

### All Swiss sidereal modes 0–46, SE_SIDM_USER and the SIDBIT options

Fixtures: `tests/fixtures/ayanamsa_all/gen_ayanamsa_all.py` (601 dates over
0–2400, yearly in 1800–2100); test `tests/ayanamsa_all.rs` prints every table
below (`cargo test -p fe-core --test ayanamsa_all -- --nocapture`).
Definitions are taken from the Swiss Ephemeris public documentation
(swisseph.htm §2.8, swephprg.htm §12, read 2026-10-29) and the primary
literature cited there; the oracle was used only to choose among published
readings (UT vs TT for a calendar t0; which Galactic pole; projection or not).
No constant below was fitted to oracle output.

**Black-box facts (asserted by the generator or the test):**
- true = mean + Δψ for every mode (also star/galactic ones).
- Star modes: `ayanamsa_ex(f) = fixstar2(star, f).lon − X` exactly (Spica 180°,
  ζ Psc 359°50′, λ Sco 240°, δ Cnc 106°); the apparent place enters (TRUEPOS
  changes it by 7–19″).
- Modes 18 J2000, 19 J1900, 20 B1950 and 34 Skydram/Mardyks are **projection
  modes** by default (swephprg.htm §3.3 note): sidereal = the position
  projected onto the mean ecliptic of t0 (Vondrák 2011), lon − A0, latitude
  on that plane; nutation never enters. Their `get_ayanamsa` is still the
  traditional value. Same rule for any fixed-epoch mode with `SIDBIT_ECL_T0`.
- With `FLG_SIDEREAL|FLG_EQUATORIAL`: traditional modes return RA/Dec on the
  mean equator of date (nutation dropped); projection modes on the mean
  equator of t0.
- Star-anchored and galactic modes ignore every SIDBIT.
- `SIDBIT_PREC_ORIG` leaves state behind inside the oracle library (all later
  ayanamsas in the process are changed): the generator runs each bit
  combination in a fresh process.

**Provenance and delta per mode** (ours − Swiss, mean ayanamsa, max |Δ| in
arcsec; built-in path; star/galactic modes with JPL states via the engine in
brackets). "UT" = calendar instant read as UT and converted with
`DeltaTModel::SWISS_COMPAT` (chosen black-box: TT reading gives 0.008–0.03″).

| # | Mode | Definition used (source) | 0–1000 | 1000–1800 | 1800–2100 | 2100–2400 |
|---|---|---|---|---|---|---|
| 0 | Fagan/Bradley | 24°02′31.36″ at B1950, Newcomb offset 0.41256″ | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | Lahiri | 23°15′00.658″ true, 1956-03-21 0h TT, offset −0.13036″ | 0.0078 | 0.0078 | 0.0078 | 0.0078 |
| 2 | De Luce | 0 at 1 Jan 1 BCE (year 0) 0h UT (birth of Jesus, DeLuce 1963) | 0.0001 | 0.0001 | 0.0001 | 0.0001 |
| 3 | Raman | 21.01444° at J1900, offset 0.828″ | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | Usha/Shashi | 18°39′39.46″ at J1900 | 0.0040 | 0.0040 | 0.0040 | 0.0040 |
| 5 | Krishnamurti | 22.363889° at J1900, offset 0.828″ | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 6 | Djwhal Khul | 30° on 1 July 2117 0h UT (Age of Aquarius 2117, G. Dawson) | 4.81 | 4.87 | 4.88 | 4.88 |
| 7 | Yukteshwar | 21.082222° at J1900 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 8 | J.N. Bhasin | Lahiri − 1°07′ (Bhasin, introductory note to his translation of *Sarvarth Chintamani*); differs from Swiss's mode of this name | 78.15 | 78.17 | 78.17 | 78.17 |
| 9–11 | Kugler 1/2/3 | −5°40′ / −4°16′ / −3°25′ at 1 Jan −100 0h UT | 0.0120 | 0.0120 | 0.0120 | 0.0120 |
| 12 | Huber | −4°28′ at 1 Jan −100 0h UT (Huber 1958, value as corrected in SE 2.05) | 0.0120 | 0.0120 | 0.0120 | 0.0120 |
| 13 | Eta Piscium | −5°04′46″ at 1 Jan −129 0h UT (Mercier) | 0.9307 | 0.9307 | 0.9307 | 0.9307 |
| 14 | Aldebaran 15 Tau | A0 = geometric mean-of-date λ(Aldebaran, 1 Jan −100 UT) − 45° (our Hipparcos catalogue) | 0.0222 | 0.0221 | 0.0221 | 0.0222 |
| 15 | Hipparchos | −9°20′ at JD 1674484.0 TT (Mercier) | 0.0069 | 0.0069 | 0.0069 | 0.0069 |
| 16 | Sassanian | 0 at JD 1927135.8747793 read as UT | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 17 | Gal. Center 0 Sag | apparent λ(Sgr A*, catalogue "Galactic Center") − 240° | 0.801 [0.805] | 0.415 [0.417] | 0.109 [0.099] | 0.149 [0.137] |
| 18–20 | J2000, J1900, B1950 | 0 at the epoch, projection mode | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 21–26 | Suryasiddhanta, …mean Sun, Aryabhata, …mean Sun, SS Revati, SS Citra | 0 / −0.21463395 / 0 / −0.23763238 / −0.79167046 / 2.11070444 at 499-03-21 07:30:31.57 UT (Aryabhata: 06:56:55.57 UT), Julian calendar | 0.0000 | 0.0000 | 0.0001 | 0.0000 |
| 27 | True Citra | apparent λ(Spica) − 180° | 0.026 [0.0132] | 0.018 [0.0016] | 0.020 [0.0000] | 0.019 [0.0001] |
| 28 | True Revati | apparent λ(ζ Psc) − 359°50′ | 0.020 [0.0022] | 0.018 [0.0006] | 0.019 [0.0001] | 0.018 [0.0003] |
| 29 | True Pushya | apparent λ(δ Cnc) − 106° (unchanged) | 0.020 [0.0043] | 0.019 [0.0015] | 0.019 [0.0006] | 0.018 [0.0005] |
| 30 | Gal. Center (Gil Brand) | apparent λ(Sgr A*) − (300° − 90°/φ) (golden section 0° Sco–0° Aqu) | as 17 | | | |
| 31 | Gal. Equator (IAU 1958) | node of the galactic equator on the mean ecliptic of date − 240°; pole IAU 1958 in ICRS (Liu, Zhu & Zhang 2011 eq. 26) | 0.281 | 0.176 | 0.176 | 0.176 |
| 32 | Gal. Equator (true) | same with the Sgr A*-based pole (Liu, Zhu & Zhang 2011 eq. 31) | 0.105 | 0.0037 | 0.0011 | 0.0028 |
| 33 | Gal. Equator mid-Mula | as 32, node at 246°40′ | as 32 | | | |
| 34 | Skydram (Mardyks) | 30° at the apparent autumnal equinox 1998 (computed with this engine, JD TT 2451079.734892), projection mode | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 35 | True Mula | apparent λ(λ Sco) − 240° | 0.037 [0.0270] | 0.019 [0.0005] | 0.020 [0.0001] | 0.019 [0.0001] |
| 36 | Dhruva/GC/Mula (Wilhelm) | polar projection of Sgr A* (great circle through the mean pole) − 246°40′ | 0.818 [0.821] | 0.435 [0.438] | 0.113 [0.103] | 0.155 [0.142] |
| 37 | Aryabhata 522 | 0 at 522-03-21 05:46:44 UT | 0.0000 | 0.0000 | 0.0001 | 0.0000 |
| 38 | Britton | −3.2° at 1 Jan year 0, 0h UT | 0.0001 | 0.0001 | 0.0001 | 0.0001 |
| 39 | "Vedic"/Sheoran | δ Cnc (Asellus Australis) at 103°29′32.9375″, true position (Sheoran, *The Science of Time and Timeline of World History*, 2017, §5.2, Table 1.16) | 0.59 | 0.59 | 0.59 | 0.59 |
| 40 | Cochrane | apparent λ(Sgr A*) − 270° | as 17 | | | |
| 41 | Fiorenza | 25° on 1 Jan 2000 0h UT | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 42 | Vettius Valens | — no published definition available (Holden 1995 not accessible) | n/a | | | |
| 43 | Lahiri 1940 | 22°26′45.50″ at J1900, offset 0.828″ | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 44 | Lahiri VP285 | 0 at 285-03-22 17:54:02 TT | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 45 | KP-Senthilathiban (VP291) | 0 at 291-03-21 06:10:29 TT | 0.0007 | 0.0007 | 0.0007 | 0.0007 |
| 46 | Lahiri ICRC | 23°15′00″ true at 1956-03-21 0h TT, Δψ IAU 2000B, Newcomb reconversion offset 0.36189″ computed from the published models (below) | 0.066 | 0.066 | 0.066 | 0.066 |

Causes of the non-zero rows:
- **Djwhal Khul 4.9″**: the oracle's value is anchored at J1900 (its 1900
  value is not published); the published definition (30° on 1 July 2117) is
  what we implement.
- **Lahiri ICRC 0.066″**: the documentation says the IENA values refer to
  Newcomb precession and that *Woolard* (1953) nutation was subtracted from
  the true value; the Newcomb reconversion offset is now computed (see
  "Reconversion offsets" below), the remaining 0.066″ equals the oracle's
  NO_PREC_OFFSET difference and is the Woolard−IAU 2000B Δψ difference at
  t0 (Woolard's series is not implemented: its coefficient table was not
  accessible).
- **Eta Piscium 0.93″, Kugler/Huber 0.012″, Hipparchos 0.007″,
  Usha/Shashi 0.004″**: constant offsets, i.e. the oracle's internal A0 or t0
  differs slightly from the documented value (rounding of the published
  numbers, or an undocumented conversion). Not reproducible from the
  documentation; not fitted.
- **Galactic Centre modes (17, 30, 36, 40)**: Sgr A* now comes from the star
  catalogue's "Galactic Center" record (SIMBAD place, Petrov et al. 2011;
  proper motion Reid & Brunthaler 2020; distance GRAVITY 2019), one source of
  truth with `fixstar`. Engine vs Swiss, max |Δ| per era: 0.805″ / 0.417″ /
  0.099″ / 0.137″ (Wilhelm 0.821 / 0.438 / 0.103 / 0.142; built-in path
  0.801 / 0.415 / 0.109 / 0.149). Black-box: the oracle's mode-17 ayanamsa is
  exactly its own Galactic Center star − 240° − 0.0197″ at all epochs, so the
  remaining drift (≈ 0.04″/century into the past) is the difference between
  our catalogue record and the oracle's star-file record of Sgr A*, not the
  ayanamsa definition; the constant 0.0197″ is unexplained.
- **Gal. Equator IAU 1958 0.18″**: the oracle converts the 1958 pole to the
  ICRS differently; the published rigorous conversion (Liu et al. 2011 eq. 26)
  is used. Other published conversions tried: Hipparcos (192.85948°,
  27.12825°) 0.19″, FK4→FK5 `fk45z` 0.10″, Reid & Brunthaler 0.12″.
- **True Citra/Revati/Mula built-in path ~0.02″**: low-precision built-in
  Earth; with JPL states (engine path) ≤ 0.0001–0.0006″ in 1800–2400. Before
  1000 the catalogue (Hipparcos vs the oracle's star file) dominates.

**Sidereal positions** (engine, 1800–2100, sidereal − tropical compared so
that ephemeris/ΔT differences cancel): traditional modes 0.0000″; J2000,
J1900, B1950, Mardyks (projection) lon 0.0000″, lat ≤ 0.0015″, RA/Dec ≤ 0.0025″;
Lahiri+ECL_T0 0.0078″ (the Lahiri residual), lat 0.0011″; speeds ≤ 0.005″/day.

**Option bits** (per call: `CalcOptions::sid_bits`, `Engine::ayanamsa_deg_bits`,
`houses_ut_bits`; Python `sidereal="lahiri|256"` or `257`):

| Bit | Semantics implemented | vs Swiss (max over 0–2400) |
|---|---|---|
| `SE_SIDM_USER` (+`SIDBIT_USER_UT`) | (t0, ayan_t0), traditional; t0 UT converted with SWISS_COMPAT ΔT | 0.00000″ (TT and UT) |
| `ECL_T0` 256 | projection onto the mean ecliptic of t0 (above) | = base-mode residual |
| `SSY_PLANE` 512 | projection onto the invariable plane (Souami & Souchay 2012: i = 1.578694°, Ω = 107.582222° on the J2000 ecliptic); lon = λ_ssy − λ_ssy(equinox t0) − A0 | lon = base residual, lat 0.042″ |
| `ECL_DATE` 2048 | zero point (λ = A0 on the ecliptic of t0) re-measured on the mean ecliptic of date | ≤ 0.0004″ (+ base residual; ICRC 0.30″) |
| `NO_PREC_OFFSET` 4096 | drop the documented precession-model offsets (FB, Lahiri, Raman, KP, Lahiri 1940) | 0.0000″ (Lahiri 0.0078″); ICRC 0.065″ |
| `PREC_ORIG` 8192 | ayanamsa with the definition's original precession (IAU 1976 Lieske 1977 for Lahiri; Bessel–Newcomb in Kinoshita's 1975 formulation, coefficients as published in SLALIB/PAL `palPrebn`, for FB, Raman, KP, Lahiri 1940, ICRC), no offset | 1800–2100 ≤ 0.006″, 2100–2400 ≤ 0.013″ (ICRC 0.066″); 0–1000: 0.11″ (Lahiri), 0.60–0.73″ (Newcomb; was 3.5–4.1″ with the classical ESAA 1961 expressions) |

**Reconversion offsets.** The documented offsets are exactly the difference
of the precession from t0 to J2000.0 between the author's model and Vondrák
2011 (`orig_precession::reconversion_offset_arcsec`): IAU 1976 at 1956-03-21
gives 0.13036″ (documented 0.13036″); Kinoshita–Newcomb gives 0.41350″ at
B1950 (documented 0.41256″) and 0.83031″ at J1900 (documented 0.828″). The
documented values are kept for those modes; the same computation supplies the
undocumented Lahiri ICRC offset (0.36189″; the oracle's NO_PREC_OFFSET
behaviour implies 0.361″).

**Fixed stars** use the same sidereal rules as planets
(`Engine::sidereal_direction`, called from `stars.rs`): projection modes,
mean-equator sidereal RA/Dec. Ecliptic sidereal star positions vs the oracle
(sidereal − tropical): traditional 0.0000″, J2000/J1900 ≤ 0.0012″,
Lahiri+ECL_T0 0.0078″, SSY 0.043″. The oracle's SIDEREAL|EQUATORIAL output
for *stars* is not the mean-equator RA/Dec it returns for planets (its RA is
shifted by ≈ the ayanamsa); we keep the planet rule and do not imitate it.

Searched without result (no accessible publication with a usable number;
second round of searches 2026-10-29 included Holden's *A History of
Horoscopic Astrology*, Bhasin's ephemerides, Sheoran's 2017 book/website,
Mercier 1976/77 (Archives Int. d'Histoire des Sciences 26–27; reprinted in
*Studies on the Transmission of Medieval Mathematical Astronomy*, 2004),
Dawson/Solar Fire for Djwhal Khul, and Woolard 1953 / Explanatory Supplement
1961 nutation tables — only bibliographic records or scanned images without
text were reachable):
J.N. Bhasin, Sheoran ("Vedic"), Vettius Valens (Holden 1995, AFA Journal of
Research 7(2)), a J1900 value for Djwhal Khul, Mercier's exact Eta Piscium
epoch/value, an IAU 1958 pole conversion reproducing the oracle (four
published conversions tried, 0.10–0.19″).

Third round (2026-10-30), still without a usable primary definition:
- **J.N. Bhasin**: his *Medical Astrology: A Rational Approach* describes the zodiac only as
  measured "from a fixed point viz the star ashvini" and about a degree below Lahiri; *Astro
  Sutras* (archive.org full text) does not discuss ayanamsa. No epoch or value.
- **"Vedic"/Sheoran**: the originator is Sunil Sheoran, *The Science of Time and Timeline of
  World History* (2017, ISBN 9781548381202). Both free-download links on the author's blog
  return 404; the definition needs a copy of the book.
- **Vettius Valens**: Holden, "The Classical Zodiac", *AFA Journal of Research* 7(2), 1995 is
  not online; Kollerstrom (*Culture and Cosmos* 1(2)) only confirms Valens's longitudes are
  sidereal (after Neugebauer, *Greek Horoscopes* p. 172).
Resolved 2026-10-30 from primary sources supplied by the product owner: Sheoran's 2017 eBook
(Table 1.16) and Bhasin's own note (Lahiri − 1°07′). Both differ from Swiss's modes of the same
name: Sheoran's δ Cnc is 0.57″ below Swiss's value (invisible at the 0.001° display), Bhasin's rule
is 1′18″ below Swiss's J.N. Bhasin mode (visible; we keep the author's rule). Vettius Valens
still needs Holden 1995 (AFA Journal of Research 7(2), p. 12).

Sidereal houses (`Engine::houses_ut_bits` / `house_cusps_ut_bits`, all
systems; measured on `houses_ex2(FLG_SIDEREAL)`, fixture `houses_sidereal.json`,
1900–2050): traditional modes subtract the true ayanamsa from cusps and angles
(not ARMC). Projection modes intersect the house circles with the projection
plane instead of the ecliptic of date: the plane's pole (J2000) is carried to
the true equator of date (Vondrák ecliptic-of-date matrix, then
R1(−ε_true)·R3(−Δψ)); houses are computed for its inclination ε′ with
ARMC′ = ARMC − Ω (Ω: RA of its ascending node on that equator), and the node's
sidereal longitude is added. The oracle returns ARMC′ in `ascmc[2]`, and uses
the zero point **without** the precession-model offset (FB +0.41256″,
Lahiri −0.13036″, KP/Raman +0.828″ — measured as constant cusp offsets of
exactly those amounts), i.e. A0 as with `SIDBIT_NO_PREC_OFFSET`. Zodiac-anchored
systems: 'N' cusps are 0°, 30°… of the sidereal zodiac in every mode; 'W' uses
the sidereal sign of the Asc in traditional modes (in projection modes the
signs of the plane's own construction, as the oracle). Deltas: J2000, J1900,
B1950 ≤ 1.7e-6°; ECL_T0 (FB, Lahiri, KP) ≤ 2.4e-6°; traditional Lahiri ≤ 2.4e-6°;
SSY_PLANE ≤ 1.0e-4° in cusps (ARMC′ 1.2e-5°: the invariable-plane orientation
residual quoted above, amplified at high latitude). Before 1900 / after 2050 the
oracle's sidereal time and ΔT dominate (~1e-3°). Speeds stay tropical (oracle
semantics). Fixed stars follow the same projection rules through
`Engine::sidereal_direction`.

Third-decimal change rates (fraction of values displayed to 0.001° that
differ from Swiss, ≈ δ/0.001°, 1800–2100): Djwhal Khul 100 %; Eta Piscium
26 %; Lahiri ICRC 1.8 %; Gal. Equator IAU 1958 5 %; Galactic-Centre modes ≈ 3 % (≤ 0.10″)
(engine, 1800–2100 ≤ 0.10″ → ≈ 3 %); Aldebaran 0.6 %; Kugler/Huber 0.3 %; Lahiri, Hipparchos 0.2 %;
Usha/Shashi 0.1 %; Gal. Equator true/mid-Mula 0.03 %; all other implemented
modes ≤ 0.02 %.

---

## 3. Lunar points (`points.rs`)

### Mean node and mean apogee (mean Lilith)

- Model (published only): ELP/MPP02 secular polynomials W3 (node) and W2
  (perigee) of Chapront & Francou (2003, A&A 404, 735; Chapront,
  Chapront-Touzé & Francou 2002, A&A 387, 700), which are referred to the
  inertial J2000 ecliptic, plus the IAU 2006 general precession in
  longitude p_A (Capitaine, Wallace & Chapront 2003):
  - W3 = 125°02′40.39816″ − 6967919.5383″T + 6.3590″T² + 0.007625″T³ − 0.00003586″T⁴
  - W2 = 83°21′11.67475″ + 14643420.3171″T − 38.2631″T² − 0.045047″T³ + 0.00021301″T⁴
  - p_A = 5028.796195″T + 1.1054348″T² + 0.00007964″T³ − 0.000023857″T⁴ − 0.0000000383″T⁵

  The ELP/MPP02 LLR/DE405 fit corrections ΔW are not applied.
  Published alternatives measured (max mean-apogee longitude delta,
  1800–2100 / 0–1000): ELP/MPP02 + IAU 2006 p_A **0.54″ / 2.90″ (chosen)**;
  same with the ELP/MPP02 DE405-mode secular corrections 0.52″ / 1.95″ (but
  node 1.39″ in 0–1000); with the LLR-mode ΔW0 and ΔW1 constants (the derived Cw
  rate terms, ≈ 0.01″/cy, omitted) 0.53″ / 3.87″; with the DE405-mode
  ΔW0/ΔW1 0.49″ / 2.57″
  (node 2100–2400 worse); IAU 1976 precession instead of p_A 1.12″ / 5.03″;
  Meeus (1998) ch. 47 (Ω, ϖ = L′ − M′) 1.83″ / 6.98″. Per the Swiss
  Ephemeris documentation the oracle uses Moshier's adjustment of
  ELP2000-85 to DE404 (Moshier 1992, A&A 262, 613) plus its own corrections
  derived from DE431; those coefficients are not available to us in a
  published form we can use, so the oracle's perigee rate (≈0.23″/cy below
  ours) cannot be matched with provenance. This is the one remaining
  visible residual (mean Lilith: 3 % of 0.001°-rounded longitudes in
  1900–2050).
- **Mean apogee on the mean orbit plane**, inclination 5°08′43.427″ =
  5.1453964° (Astronomical Almanac, Moon section p. D2, editions before 2005;
  the USNO 2005 errata replaced the printed value by 5.1566898°, the ELP mean
  inclination 2·asin γ). We keep the classic AA value, which the oracle also
  uses; with 5.1566898° mean-Lilith latitudes would change by up to 0.0113°.
  `λ = Ω + atan2(sin u cos i, cos u)`, `β = asin(sin u sin i)`, `u = ϖ + 180 − Ω`.
- Distances are conventions (not positions): node 384 400 km (AA mean
  distance), apogee 384 400 km × (1 + 0.054900489) (AA mean eccentricity);
  these equal the oracle's values exactly (0.0025695552898 AU,
  0.00271062513172 AU).
- Swiss default flags = NONUT value + Δψ in longitude, latitude unchanged
  (asserted for MEAN_NODE, MEAN_APOG, TRUE_NODE, OSCU_APOG).
- The oracle's mean node/perigee are quartics of its own (fitted by us only
  for this comparison, not shipped); vs our published model, ours − Swiss:

| Era | mean node lon | mean apogee lon | mean apogee lat |
|---|---|---|---|
| 0–1000 | ≤ 1.20″ (3.3e-4°) | ≤ 2.90″ (8.1e-4°) | ≤ 0.25″ |
| 1000–1800 | ≤ 0.14″ | ≤ 2.68″ (7.4e-4°) | ≤ 0.24″ |
| 1800–2100 | ≤ 0.034″ (9.4e-6°) | ≤ 0.54″ (1.5e-4°) | ≤ 0.047″ |
| 2100–2400 | ≤ 0.15″ | ≤ 0.36″ (1.0e-4°) | ≤ 0.027″ |

Speeds: analytic derivative; |Δ| ≤ 1.4e-7°/day (oracle speeds are numeric).

### Osculating node and apogee (true node, true Lilith) — the Moon state

Measured (all alternatives tested and rejected: apparent / light-time
barycentric / date-frame velocity / finite-difference velocity):

1. Moon = **geometric geocentric** state (`FLG_TRUEPOS|FLG_NONUT`), **at the
   light-retarded time `t − |r(t)|/c`** (position *and* velocity at that
   time; a single light-time step). Without this step the node is off by
   ≈ 1.3 s of motion (1e-6° node, 1e-4° apogee).
2. In an **inertial frame** (J2000): h = r × v, e = v × h / GM − r/|r|;
   node direction ẑ × h, apogee = −ê · a(1+e).
3. The oracle's GM is 3.98600448e14 m³ s⁻² × (1 + 1/81.30056)
   (= 403503.24877 km³/s², identified black-box). **Not shipped**: the
   engine passes its ephemeris' GM_E + GM_M (DE440: 403503.235648 km³/s²);
   the GM difference alone moves the osculating apogee by ≤ 5.4e-5° in
   longitude and ≤ 4.6e-6° in latitude (measured on the 701 fixture states)
   and the node not at all. The test fixture comparison uses the
   oracle's GM locally (test-only).
4. Vectors rotated to the **mean ecliptic of date** with Swiss' frame
   rotation, which is `R1(ε_mean,Swiss) · P_Vondrák2011(equator) · R1(−ε0)`
   (verified to 5e-12; note: Swiss' mean obliquity is *not* the V11
   pole-to-pole angle, and equals IAU 2006 `obl06` to 3e-10° at 2000, 4.5e-8°
   at 1800, −1.4e-7° at 2400, 2.8e-6° at year 0, 5.6e-6° at year 300).
   Using ERFA `ltecm` (V11 ecliptic) instead gives node errors up to 3.3e-4°
   (year 271), 1e-5° (≈1100), ≤ 5.5e-6° (1800–2300).
5. Apparent output = mean-equinox longitude + Δψ.

Swiss' date-frame `FLG_XYZ|FLG_SPEED` velocities include the frame-rotation
term; do not feed those to `osculating_*`.

### Osculating formulas vs Swiss (701 dates, years 0–2400; oracle Moon state and GM as inputs)

| Point | lon | lat | other |
|---|---|---|---|
| true node (given Swiss' Moon state + rotation) | 1.2e-11° | 0 | node distance rel. 4.7e-11 |
| osculating apogee (same) | 3.6e-8° | 3.1e-9° | distance rel. 4.8e-11 |

---

### Mean elements refitted to DE441 (tested, not adopted)

A study refitted the ELP/MPP02 mean-argument polynomials W1, W2, W3 (corrections
c·T^k, k = 0..4, plus a constant frame rotation) to the geocentric DE441/DE442 Moon
over years 0–2400 (30 000 epochs, Gauss–Newton with the analytical partials of the
full series; no oracle data entered the fit). It reduced the Moon-position rms of the
theory before 1800 by one to two orders of magnitude (longitude 11.8″ → 0.08″ in
0–1000), but moved the mean node and mean apogee *away* from the oracle: in 1900–2050
the node went from 0.012″ to 0.15″ and the apogee from 0.25″ to 0.40″, and before 1000
the node from 1.2″ to 10″.

**Decision: not adopted, and the table has been removed** (it also depended on a
GPL-licensed reformatting of the series files). The oracle's node agrees with the
nominal published ELP/MPP02 W3 + IAU 2006 p_A to 0.012″ today, while DE441 prefers a
node 0.15″ different; the oracle's perigee keeps the nominal ELP constant
83°21′11.67475″, while DE441 (and ELP's own LLR fit, ΔW2 = +0.168″) put it ~0.1–0.26″
higher. So the oracle's mean elements are not tied to DE431 at the 0.1″ level, and no
fit to JPL data can reproduce them; the residual mean-Lilith difference (≤ 0.25″ in
1900–2050, 2.9 % of 0.001° values) remains. Long-span (pre-1000) disagreements of such a
fit also reflect that "mean" elements over two millennia depend on the
ecliptic-precession convention (ELP's Laskar P,Q) more than on the Moon itself.

## Provenance (every shipped constant → its source)

Rule: no constant or table fitted to Swiss output ships; Swiss is used only
for semantics and for the compatibility deltas above.

| Constant (file) | Value | Source |
|---|---|---|
| Vondrák long-term precession series (`ayanamsa/vondrak.rs`) | tables PQPOL/PQPER/XYPOL/XYPER, ε0 = 84381.406″ | Vondrák, Capitaine & Wallace 2011, A&A 534, A22 (+2012 erratum) |
| Lahiri t0, value (`ayanamsa.rs`) | 1956-03-21 0h TDT, true 23°15′00.658″ | Indian Astronomical Ephemeris (1985) refinement of the ICRC value (CSIR 1955); interpretation per Swiss Ephemeris documentation |
| `DPSI_LAHIRI_T0_DEG` | 0.004660197838803416° | computed: IAU 2000B nutation (McCarthy & Luzum 2003, `frames::nut00b`) at t0; unit-tested |
| Precession-model offsets | Lahiri −0.13036″ (IAU 1976), FB +0.41256″, Raman/KP +0.82800″ (Newcomb) | Swiss Ephemeris documentation, ayanamsha section (published text) |
| Raman | 21.01444° at J1900.0 | B. V. Raman (published 1900 value) |
| Krishnamurti | 22.363889° at J1900.0 | K. S. Krishnamurti (*KP Reader I* table) |
| Yukteshwar | 21.082222° at J1900.0 | Sri Yukteshwar ayanamsa as published |
| Fagan–Bradley | 24°02′31.36″ (SVP 335°57′28.64″) at B1950.0 | C. Fagan & D. Bradley |
| True Pushya star longitude | 106° (16° Cancer) | P. V. R. Narasimha Rao (definition of True Pushya paksha) |
| δ Cnc astrometry (`ayanamsa/star.rs`) | α, δ, μ, ϖ at J1991.25 | Hipparcos New Reduction, van Leeuwen 2007, A&A 474, 653 (HIP 42911) |
| δ Cnc radial velocity | +16.39 km/s | Massarotti et al. 2008, AJ 135, 209 |
| Frame bias | ξ0 −16.617 mas, η0 −6.8192 mas, dα0 −14.6 mas | IAU 2006 (Capitaine et al. 2003; IERS Conventions 2010) |
| All other sidereal modes (`ayanamsa.rs` `kind()`) | t0 / A0 / anchor per mode | §2 "All Swiss sidereal modes" table (Swiss documentation swisseph.htm §2.8 and the authors cited there) |
| Sgr A* (`ayanamsa/galactic.rs`) | α 17h45m40.0409s, δ −29°00′28.118″, μ (−3.151, −5.547) mas/yr | Reid & Brunthaler 2004, ApJ 616, 872 (parallax 0.125 mas = 8 kpc, RV 0: our choice, immaterial) |
| Galactic poles | IAU 1958 in ICRS 12h51m26.27469s +27°07′41.7087″; Sgr A*-based 12h51m36.7151981s +27°06′11.193172″ | Liu, Zhu & Zhang 2011, A&A 526, A16, eqs. 26 and 31 |
| Invariable plane | i 1.578694°, Ω 107.582222° (J2000 ecliptic) | Souami & Souchay 2012, A&A 543, A133 |
| IAU 1976 / Newcomb precession (`ayanamsa/orig_precession.rs`, PREC_ORIG only) | ζ, z, θ, ε polynomials | Lieske et al. 1977, A&A 58, 1; Explanatory Supplement 1961 |
| `MARDYKS_T0_TT` | 2451079.734892 | computed with this engine (apparent Sun at 180°); unit-tested |
| Low-precision Sun | Meeus series | Meeus, *Astronomical Algorithms* 2nd ed., ch. 25 |
| GM_Sun, c | 1.32712440041e20 m³ s⁻², 299 792 458 m/s | IAU 2009/2015 / SI |
| Mean node W3, perigee W2 (`points.rs`) | see §3 | Chapront & Francou 2003 (ELP/MPP02), Chapront, Chapront-Touzé & Francou 2002 |
| p_A | 5028.796195″T … | Capitaine, Wallace & Chapront 2003 (IAU 2006 P03) |
| `MEAN_INCLINATION_DEG` | 5.1453964° | Astronomical Almanac, p. D2 (pre-2005 editions) |
| `MEAN_MOON_DIST_KM`, `MEAN_ECCENTRICITY` | 384 400 km, 0.054900489 | Astronomical Almanac lunar constants (conventions, used only as reported distances) |
| `points/elp_main.rs` (natural apogee/perigee) | ELP 2000-82B series, 918 terms | generated by `tools/gen_intp_series.py` from the authors' files ELP1–ELP36, IMCCE ftp://ftp.imcce.fr/pub/ephem/moon/elp82b/ (checksums in docs/FICTITIOUS_BODIES.md) |
| `MAX_LAT_DEG` (`houses.rs`) | 90 − 1e-10° | numerical guard (our choice; behaviour-level; coincides with the oracle's) |
| Placidus/topocentric solver tolerances, degeneracy thresholds | 1e-13 rad, 1e-12 | numerical choices |

Removed in this pass (were oracle-derived): fitted mean node/perigee
quartic rates (`NODE_C`, `PERI_C`), oracle A0 values for all five fixed
ayanamsas, the oracle's δ Cnc catalogue, `GM_EARTH_MOON_SWISS_KM3_DAY2`
(now test-only), and the 1-mas `house_pos` offset/snap
(`HOUSE_POS_OFFSET_DEG`).

### Can the deltas flip an API classification (values shown to 0.001°)?

A constant longitude delta δ flips a sign (30°) or nakshatra (13°20′)
assignment only for points within δ of a boundary: probability per point
≈ δ/segment; and changes the third displayed decimal for ≈ δ/0.001 of values.

| Item (era 1800–2100 unless noted) | δ | 3rd decimal changes | nakshatra flip / point | sign flip / point |
|---|---|---|---|---|
| Lahiri sidereal longitudes & cusps | 2.2e-6° | ~0.2 % | 1.6e-7 | 7e-8 |
| Raman, KP, Fagan–Bradley, Yukteshwar | ≤ 1.3e-8° | ~0.001 % | 1e-9 | 4e-10 |
| True Pushya (engine path: JPL states) | ≤ 1.7e-7° (≤ 1.2e-6° before 1000) | ~0.02 % (0.1 %) | 1e-8 (9e-8) | 6e-9 (4e-8) |
| True Pushya (built-in Earth fallback) | ≤ 5.4e-6° | ~0.5 % | 4e-7 | 2e-7 |
| Mean node (Rahu/Ketu mean) | 9.4e-6° (3.3e-4° before 1000) | 0.2 % (12 % in 0–1000), measured | 7e-7 (2.5e-5) | 3e-7 (1.1e-5) |
| Mean Lilith lon | 7.0e-5° in 1900–2050, 1.5e-4° in 1800–2100, 8.1e-4° in 0–1000 | 3.0 % / 5.9 % / 74 %, measured on 4000 random dates | 5e-6 / 1.1e-5 / 6e-5 | 2e-6 / 5e-6 / 2.7e-5 |
| house_pos | 2.8e-7° (1 mas) | — | integer house differs only within 1 mas of a cusp (~1e-8 per point) | |

All ayanamsas are now at or below 0.2 % third-decimal changes. The one
visible residual is **mean Lilith** (3 % of displayed values in 1900–2050,
6 % in 1800–2100, most values in antiquity), see §3 — no published
mean-perigee theory reproduces the oracle's (Moshier's DE404 fit plus
DE431 corrections, per the Swiss documentation, whose coefficients are not
published in an accessible form).

## 4. What the integration layer must do

- **Ayanamsa**: `true = ayanamsa_mean_deg(mode, jd_tt) + Δψ(jd_tt)`.
  `get_ayanamsa_ut(jd_ut)` corresponds to `ayanamsa_mean_deg(mode, jd_ut + ΔT)`
  (residual ≤ 2.2e-6°, Lahiri).
  Sidereal longitude = tropical (same flags) − (mean + Δψ if nutation is
  applied, else mean). Sidereal speed = tropical speed − ayanamsa rate
  (`ayanamsa_mean_rate_deg_per_day` + dΔψ/dt when apparent).
- **True Pushya / Pushya**: prefer
  `true_pushya_deg(star::delta_cnc_apparent_mean_from_states(jd_tt, E_bary, V_bary, Sun_bary).0)`
  with JPL states (≤ 1.7e-7° vs the oracle in 1800–2100, ≤ 1.2e-6° before
  1000); the built-in Earth fallback is ≤ 5.4e-6°.
- **Houses**: ARMC = GAST·15 + λ_geo (apparent) and ε = true obliquity at
  `jd_ut + ΔT`; for `FLG_NONUT`, GMST and mean ε. Sidereal: subtract the true
  (or, with NONUT, mean) ayanamsa from cusps and ascmc[0,1,3..7], never from
  ARMC. Placidus/Koch above the polar circle: today's API raises; use
  `houses_swiss_compat` only if the fallback is wanted.
  Placidus cusps amplify an obliquity error by (measured, max over ARMC)
  0.2× at φ=0, 1.6× at 45°, 4.8× at 60°, ~50× at 66°: frames must match
  Swiss' ε (see §3 item 4) to reach 1e-7° charts.
- **house_pos**: call with `(armc, lat, 23.44, lon, 0.0)` to match FreeAstroAPI.
- **True node / true Lilith**: geometric geocentric Moon at `t − |r|/c`,
  rotated with the frames module's J2000→mean-ecliptic-of-date matrix
  (which must equal Swiss' `R1(ε_mean) · P_V11` for the < 1e-6° targets),
  the ephemeris' GM_E + GM_M, then + Δψ. Our Moon (JPL DE441/442 direct)
  will differ from Swiss' compressed `semo` files; the osculating apogee
  amplifies a relative Moon-velocity error δ by ≈ δ/e ≈ 18 δ rad — to be
  measured in the harness once `ephem` exists.
- **Mean points**: `mean_node_deg(jd_tt)`, `mean_apogee_deg(jd_tt)` (+ Δψ for
  apparent); distances `MEAN_NODE_DIST_AU`, `MEAN_APOGEE_DIST_AU`.

## 5. Open issues

1. Swiss' Placidus cusps carry up to 1.0e-6° iteration error near the polar
   circle (5e-8° at |φ| ≤ 60); we are exact. Not reproducible black-box.
2. Topocentric `house_pos`: ≤ 5.5e-6 house-degrees difference; circumpolar
   points take a different root in Swiss.
3. Near-degenerate Asc (ecliptic ≈ horizon, ~1e-4° neighbourhood) and Asc on
   an exact sign boundary (Whole sign) differ at the 1e-10..5e-5° level.
4. True Pushya while δ Cnc is behind the solar disc (≈ ±6 h/yr): 1.3e-4°.
5. Published-definition deltas (Provenance section) are intentional; if
   exact oracle parity of sidereal output is ever required it must come
   from a compatibility layer outside the shipped constants.
6. Swiss mean obliquity model (needed by frames for houses/true node to
   < 1e-7°) is not IAU 2006 outside ≈1800–2200 — characterised above, not
   implemented here.

## Quick wins (2026-10-29): output flags and planet-centric positions

Tests: `tests/quickwins.rs` (`output_flags_vs_swiss`, `calc_pctr_vs_swiss`), 60 dates
1900–2100 × 8 bodies × 16 flag sets, TT input.

- `CalcOptions` gained `xyz`, `radians`, `icrs` (Python flags `FLG_XYZ`, `FLG_RADIANS`,
  `FLG_ICRS`; `FLG_ASTROMETRIC` = NOABERR|NOGDEFL; `FLG_SPEED3` accepted as SPEED).
  XYZ is the cartesian form of the final polar output (sidereal included) with chain-rule
  velocities; RADIANS converts angles and angular rates. Neither adds error: max vs Swiss
  equals the plain apparent place (planets 0.23″ — Pluto after 2050 —, Moon 0.005″).
- **J2000 equatorial now includes the IAU 2006 frame bias**, as Swiss does without
  `SEFLG_ICRS` (Moon 0.024″ → 0.005″). `icrs` removes the bias everywhere (the ICRS vector is
  precessed as if it were mean J2000), matching Swiss's ICRS outputs to the same 0.005″.
  Lunar points (nodes, apogees) still ignore J2000/equatorial frames (pre-existing).
- SPEED3: Swiss's three-point Moon speed differs from its own SPEED by up to 16″/day; ours is
  the same derivative for both flags.
- `calc_pctr` (`Origin::Body(naif)`): light time, solar deflection and aberration from the
  centre's barycentric velocity (Swiss confirms all three by flag tests). Positions ≤ 0.49″
  vs Swiss (all flag variants, Moon ≤ 0.22″), longitude speeds ≤ 0.07″/day. The oracle's
  distance speeds for fast centres (e.g. Neptune from Mercury) differ from the derivative of
  its positions by up to 7e-5 AU/day; ours are the derivative. Like the oracle, centre flags
  HELCTR/BARYCTR force an astrometric place and TOPOCTR is ignored; lunar points as target or
  centre are refused (Swiss accepts them).
