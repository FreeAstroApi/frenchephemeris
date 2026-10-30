# Minor bodies: Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Eros

Status 2026-10-28. Pack `fepk-minor-v2.bin` (6,684,336 B), built by
`tools/fit_minor.py` from JPL data and, before JPL's 1599-12-11 limit, from
our own integration of the JPL orbits (section "Extension before 1600").
v1 (JPL spans only) is reproduced by `fit_minor.py --no-integration`. Validation numbers below come
from `harness/check_minor.py` (full output: `harness/minor_report.json`) and
`cargo test -p fe-core --test minor`.

## Ids

| body | API / Swiss id | also accepted | NAIF id in pack | MPC |
|---|---:|---|---:|---:|
| Chiron | 15 | 12060 | 2002060 | 2060 |
| Pholus | 16 | 15145 | 2005145 | 5145 |
| Ceres | 17 | 10001 | 2000001 | 1 |
| Pallas | 18 | 10002 | 2000002 | 2 |
| Juno | 19 | 10003 | 2000003 | 3 |
| Vesta | 20 | 10004 | 2000004 | 4 |
| Eris | 146199 (`AST_OFFSET` + 136199) | | 2136199 | 136199 |
| Eros | 10433 (`AST_OFFSET` + 433) | | 2000433 | 433 |

How FreeAstroAPI asks for them today (read from the FreeAstroAPI code):
`core/astrology.py` `EXTRA_BODIES` uses `swe.CHIRON`, `swe.PHOLUS`, `swe.CERES`,
`swe.PALLAS`, `swe.JUNO`, `swe.VESTA` (15–20) with `FLG_SWIEPH|FLG_SPEED`;
the natal request (`backend/routers/v1/natal.py`) has `chiron` (default **true**),
`pholus`, `ceres`…`vesta`, `eris`, `eros` toggles; the chart package
`dependencies/natal_pkg/natal/const.py` passes Eris = 146199 and Eros = 10433
but **Pholus = 5145**, which is not a valid Swiss id (Swiss raises an error at
every date; Pholus is 16 or 15145). That is a bug in the API, not mapped here:
`Body::from_id(5145)` is `UnknownBody`, as in Swiss.

Python: `frenchephemeris.CHIRON … VESTA`, `ERIS`, `EROS`, `AST_OFFSET`
(also in `swe_compat`). `Engine(path, swiss_delta_t, extra_packs=[minor])`;
`open_default()` loads `fepk-minor-v2.bin` when it sits next to the main pack
(override with `FRENCHEPHEMERIS_MINOR_PACK`, `""` disables).

## Data sources (provenance ledger: `docs/MINOR_PROVENANCE.json`)

| body | JPL source | solution | integration / notes |
|---|---|---|---|
| Ceres, Pallas, Juno, Vesta | `sb441-n16.bsp` (645,727,232 B, Last-Modified 2021-03-31), https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/ — read with HTTP Range requests (≈90 MB cache incl. directory reads) | SB441-N16: the 16 asteroid perturbers integrated with DE441 (JPL IOM 392R-21-005) | type 2, centre Sun, 32-day records, JD −1200525.5 … 5008242.5 |
| Chiron | Horizons SPK `20002060`, 14.2 MB | JPL#171 (soln 2026-Jun-05, 6026 obs 1895–2026) | type 21, centre Sun; planets+Moon, SB441-N16 asteroids, PPN, Earth J2 |
| Pholus | Horizons SPK `20005145`, 14.3 MB | JPL#24 (soln 2024-Oct-11, 397 obs 1977–2018) | idem |
| Eris | Horizons SPK `20136199`, 14.4 MB | JPL#103 (soln 2026-Jun-06, 4064 obs 1954–2026) | idem |
| Eros | Horizons SPK `20000433`, 14.4 MB | JPL#659 (soln 2021-May-24, 9130 obs 1893–2021) | idem |
| Chiron, Pholus, Eris, Eros before 1599-12-11 | "FrenchEphemeris integration from JPL#171 / #24 / #103 / #659" | the JPL SPK state at 1599-12-11 | our integration (`tools/integrate_minor.py`); covariance from the SBDB API; ledger key `integrations` |

Horizons SPKs requested 2026-10-28 via the Horizons API (`EPHEM_TYPE='SPK'`,
`START_TIME='1599-12-11'`, `STOP_TIME='2401-04-01'`; one request per body,
the API accepted the 801-year span). Files and sha256 are in the ledger and
the pack manifest. Downloads total ≈150 MB under
`data/sources/minor/` (gitignored). `tools/fetch_minor.py` re-fetches missing
files and rebuilds the ledger.

Licence / attribution: JPL Solar System Dynamics data (NASA/JPL-Caltech); no
licence terms are attached to the files. Credit "JPL Horizons On-Line
Ephemeris System, Solar System Dynamics Group, JPL/Caltech" and, for the four
big asteroids, Park et al. 2021 (DE440/DE441, AJ 161:105; IOM 392R-21-005).
The Horizons SPK comments call the files informal releases not for
spacecraft-safety use. Confirm attribution wording with JPL SSD before
public redistribution of the derived pack (open item).

Clean room: Swiss Ephemeris is only called through `pyswisseph` in the
harness, as a black-box oracle; none of its files or outputs feed the pack.

## Supported ranges (explicit; outside → `Error::OutOfRange`)

| body | JD TDB | dates | why it stops there |
|---|---|---|---|
| Ceres, Pallas, Juno, Vesta | 1720624.5 – 2598072.5 (geocentric from 1720628.5, the main pack's Earth) | −2 Oct 29 (Jul.) – 2401 Mar 7 | same as the planets; sb441-n16 covers far beyond |
| Chiron | 1967601.5 – 2598072.5 | 675 Jan 1 (Jul.) – 2401 Mar 7 | Swiss's first date; JPL from 1599-12-11, our integration before |
| Pholus | 1720624.5 – 2598072.5 (geocentric from 1720628.5) | −2 Oct 25 – 2401 Mar 7 | main pack start; integration before 1599-12-11 |
| Eris, Eros | 2268932.5 – 2598072.5 | 1500 Jan 1 (Jul.) – 2401 Mar 7 | Swiss's first date; integration before 1599-12-11 |

**Horizons has no small-body ephemeris before A.D. 1599-DEC-11 TDB** (it
answers "No ephemeris for target … prior to") and no NAIF kernel covers
these four earlier, so everything before that epoch is our model
extrapolation of the JPL orbit (next section), with growing uncertainty.

Tests (`supported_ranges_are_explicit`, `integrated_extension_joins_jpl_continuously`)
check that the pack's chain coverage equals these ranges, that the joint at
1599-12-11 is continuous (< 0.05 km, < 0.5 km/day), that 1e-6 d outside either end is `OutOfRange` (with the
NAIF id), that repeated calls give the same error, and that an engine without
the minor pack reports `UnknownBody` rather than a wrong value.

Reliability inside the range (Horizons 3-sigma RA/Dec uncertainty, QUANTITIES
36, max over the sampled epochs): Ceres ≤ 0.003″, Juno ≤ 0.7″, Pallas ≤ 0.5″,
Vesta ≤ 0.9″, Eros ≤ 0.19″, Chiron ≤ 3.0″ (1600–1800), **Pholus ≤ 114″**,
**Eris ≤ 151″** (1600–1800; 27″ in 1800–2100). Pholus (short arc, 1977–2018)
and Eris (very slow orbit) are intrinsically uncertain at the arcminute level
away from the present; Chiron's close Saturn encounters lie outside
1600–2400. These are JPL's own formal uncertainties, not our error.

## Extension before 1600: our integration of the JPL orbits

`tools/integrate_minor.py` drives `crates/fe-core/examples/integrate_minor.rs`:
a massless particle, barycentric ICRF, integrated with a 6-stage
Gauss-Legendre collocation method (implicit RK, order 12, symplectic; stage
equations iterated to round-off; Kahan-compensated update), fixed step
(Chiron, Pholus 4 d, Eris 8 d, Eros 0.5 d). Forces: Sun, Mercury and Venus
barycentres, EMB (Earth + Moon separately for Eros), Mars–Pluto system
barycentres, all from the main pack (DE441 before 1550, DE442 after), GMs
from the main pack manifest (DE442 table); Ceres, Pallas, Juno, Vesta from
the minor pack (sb441-n16) with DE440 GMs rounded to 4–5 digits; 1PN
(Schwarzschild) term of the Sun. No Swiss data anywhere.

Initial state: the JPL Horizons SPK state (position and velocity) at its
first epoch, JD 2305426.5 = 1599-12-11 TDB, + the main pack's Sun;
integrated backwards. The joint is continuous by construction (pack: jump
< 0.05 km, i.e. one double of time × velocity).

Numerical checks (`docs/MINOR_INTEGRATION.json`):

| body | span | steps | h vs h/2, max km | round trip to 1599, km |
|---|---|---:|---:|---:|
| Chiron | 675–1600 | 84,456 | 9 (after the 720 Saturn encounter; 0.002 before 1200) | 4.3 |
| Pholus | −2–1600 | 146,200 | 0.0005 | 0.0002 |
| Eris | 1500–1600 | 4,562 | 6e-6 | 2e-6 |
| Eros | 1500–1600 | 72,988 | 2e-5 | 7e-5 |

Model check against JPL over 1600–2401 (same model, JPL state at 2401,
integrated backwards 800 years — the production situation mirrored in time):

| body | 100 y | 200 y | 400 y | 600 y | 800 y |
|---|---:|---:|---:|---:|---:|
| Chiron | 0.016″ | 0.10″ | 0.54″ | 1.6″ | 3.7″ |
| Pholus | 0.004″ | 0.08″ | 0.11″ | 0.19″ | 0.19″ |
| Eris | 0.003″ | 0.016″ | 0.020″ | 0.07″ | 0.22″ |
| Eros | 0.10″ | 0.29″ | 0.59″ | 0.92″ | 0.94″ |

(max geocentric angle vs JPL). This bounds our model error; part of it is
not ours: Horizons integrates with DE441 planets, while 1550–2401 of our main
pack is DE442, and the two differ by ~2,200 km for Uranus in 1550
(DATA_REPORT seam table). That is why Chiron, which passes 2.2 AU from
Uranus in 1679, shows ~6,000 km over 1600–1700 in forward tests from 1600
(28″ by 2400) while the other windows agree to ~100 km per century. Before
1550 our planets are DE441, as in JPL's integrations. A fitted solar-GM scale
(standing in for the 12 SB441 asteroids we omit) came out with inconsistent
signs (Chiron +3.6e-9, Pholus −1.2e-8, Eris −8.8e-9, Eros +6.4e-9) and made
Eros worse, so production uses the plain model.

Orbit uncertainty (the dominant term): JPL's SBDB covariance (elements at
the solution epoch) is mapped to Cartesian, and its Cholesky columns are
integrated as 6 clones from the solution epoch to the start; 3σ of the
transverse position, seen from Earth:

| body | 1600 | 1500 | 1200 | 1000 | 900 | 800 | 750 | 700 | ~675 | 500 | 200 | 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Chiron × 4 (see note) | 0.8″ | 1″ | 0.6″ | 20″ | 49″ | 49″ | 120″ | 900″ | 1.2° | – | – | – |
| Pholus | 12″ | 7″ | 43″ | 68″ | 41″ | 38″ | 340″ | 34″ | 261″ | 51″ | 2″ | 18″ |
| Eris | 34″ | 11″ | | | | | | | | | | |
| Eros | 0.1″ | 0.23″ | | | | | | | | | | |

(values at the nearest 10-year sample; they oscillate with the viewing
geometry, max before 1600: Pholus 590″, Eris 34″, Eros 0.23″). The mapping
reproduces Horizons' own 3σ ellipse (QUANTITIES 37) within 5–20 % for Pholus
(1630: 99″ vs 119″; 1910: 15.6″ vs 18.7″), Eris (1719: 138″ vs 144″) and
Eros (1755: 0.12″ vs 0.13″), but is 3.5–4.7× smaller than Horizons for
Chiron (1797: 0.67″ vs 2.3″; 1600: 0.2″ vs 0.94″), so the Chiron row above
is multiplied by 4. **Chiron passes 0.18 AU from Saturn in 720 AD**: before
that its position is essentially undetermined (3σ ~1° at 679), and the
step-size check shows the encounter amplifying even round-off (h vs h/2:
0.002 km before 1200, 9 km after 720). Pholus's uncertainty comes from its
short observed arc (1977–2018), not from encounters (closest: Saturn 1.7 AU
in 1447).

Fit into the pack: the integrated spans are separate series of the same
pairs (same centre as the JPL span), fitted like the JPL spans (local record
halving available for encounters, not needed): Chiron 334 records × 1024 d ×
19 coefficients (max 0.036 km vs the integration), Pholus 572 × 1024 d × 18
(0.044 km), Eris 36 × 1024 d × 6 (0.013 km), Eros 286 × 128 d × 19 (0.040 km;
heliocentric, restarting at the 1550 seam); +546 KB in total. Source entries
in the manifest are "FrenchEphemeris integration from JPL#…".

## Chiron before ~1000 AD vs Swiss: where the difference comes from

Swiss's documentation (swisseph.htm §2.3.2, §2.3.4, read as documentation
only): its asteroid ephemerides come from its authors' own integration
(Moshier's integrator, built to reproduce DE200) started from Bowell/Lowell
`astorb.dat` elements, with all planets, the Moon, Ceres, Pallas and Vesta;
"Positions of Chiron can be well computed for the time between 700 CE and
4650 CE … close encounters with Saturn in Sept. 720 CE … earlier Chiron
ephemerides … are meaningless." The documentation names neither the element
epoch nor the astorb.dat version.

Controlled experiments (`tools/chiron_study.py` → `docs/CHIRON_STUDY.json`;
our integrator, geocentric angle in arcsec relative to production; FE − Swiss
is the black-box column from `harness/check_minor.py`):

| year | FE − Swiss | JPL 3σ (×4, ≈ Horizons) | MPC orbit − JPL#171 | no asteroids | no GR | DE200 masses | Saturn +6,500 km | Uranus +2,200 km | our model from 2026 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1600 | 0.66 | 0.80 | 0.10 | 0 | 0 | 0 | 0 | 0 | 0.89 |
| 1500 | 0.61 | 0.97 | 0.14 | 0.004 | 0.03 | 1.2 | 0.01 | 0 | 0.81 |
| 1200 | 0.58 | 0.59 | 0.06 | 0.007 | 0.06 | 15 | 1.1 | 0.02 | 1.4 |
| 1100 | 2.5 | 2.4 | 0.16 | 0.03 | 0.17 | 51 | 4.3 | 0.08 | 5.7 |
| 1000 | 19 | 20 | 1.2 | 0.23 | 1.2 | 367 | 31 | 0.6 | 45 |
| 900 | 47 | 49 | 3.1 | 0.59 | 3.0 | 924 | 78 | 1.5 | 114 |
| 800 | 49 | 49 | 3.0 | 0.58 | 2.9 | 803 | 76 | 1.5 | 112 |
| 760 | 213 | 216 | 15 | 2.9 | 15 | 4,220 | 386 | 7.5 | 563 |
| 721 | 511 | 580 | 24 | 4.7 | 23 | 6,803 | 620 | 12 | 905 |
| 700 | 5,555 | 909 | 53 | 53 | 268 | 16,608 | 9,880 | 136 | 5,658 |
| 680 | 15,090 | 4,499 | 280 | 151 | 764 | 11,286 | 26,464 | 390 | 18,889 |

Columns: "MPC orbit − JPL#171" = the MPC's independent current orbit
(orbfit, DE431, 17 asteroids, 1895–2026) vs JPL#171, both integrated with our
model from the 2026 epoch; "our model from 2026" = JPL#171's 2026 state
integrated with our model instead of starting from JPL's 1599 state (it
inherits our 1600–2026 model error of 0.9″, mostly the DE441/DE442 Uranus
difference). Step halving changes nothing (< 0.005″). Chiron passes 0.182 AU
from Saturn on JD 1984177.5 (720.4 AD) in our integration, as Swiss's docs say.

Findings:

1. **FE − Swiss equals JPL's 3σ to within 5 % from 1600 down to 760 AD**
   (ratio 0.95–0.99 at every 10-year sample, through the geometry
   oscillations: 26.5″ vs 27.9″ at 850, 121″ vs 123″ at 740). Swiss's Chiron
   is therefore today's orbit displaced along the least-determined direction
   (essentially the mean motion), by about JPL's (Horizons-scaled) 3σ: what one
   expects of an older element set fitted to a shorter arc. After the 720
   encounter the difference grows non-linearly (6× the 3σ at 700).
2. **(a) initial orbit is the only candidate of the right size and shape.** Two
   current published solutions (JPL#171, MPC) agree to 1.2″ at 1000 and 15″ at
   760, i.e. at ~6 % of FE − Swiss: today's orbit is settled, and Swiss's
   deviation is what an older Bowell element set would give. We could not
   test that directly: the historical `astorb.dat` snapshots are not online at
   Lowell (only the current file), the Internet Archive was offline, and the
   PDS SBN archive is behind a bot challenge (not bypassed). MPC/MPEC element
   histories were not searched further.
3. **(b) model differences are either far too small or far too large.** The
   asteroids (0.6″ at 900), relativity (3″), Uranus at the DE441/DE442 level
   (1.5″) are 15–80× smaller than FE − Swiss. DE200 planetary masses applied to
   today's orbit give 367″ at 1000 (20× too large) and a 1″ Saturn ephemeris
   error gives 31″ at 1000 — so Swiss cannot be "today's orbit with a DE200
   model", and its result depends on how its starting elements were fitted,
   which only the element set itself could settle. All along-track effects
   share the 3σ shape, so shape alone cannot separate (a) from (b).
4. **(c) chaos**: before ~760 every one of these choices moves Chiron by
   arc-minutes to degrees (Saturn +6,500 km: 2.7° at 700; MPC vs JPL: 53″ at 700,
   280″ at 680). Values there are meaningless in every engine, as Swiss says.

**Conclusion: not legitimately reducible now.** Reproducing Swiss would need
its starting elements; the only other route, adjusting our orbit to match
Swiss, is fitting to Swiss output and is excluded. Our Chiron uses the current
JPL solution, which by every comparison above is the better-determined orbit;
the Swiss differences down to 760 AD are within JPL's 3σ. If the old
`astorb.dat` (late 1990s) becomes retrievable, `tools/chiron_study.py` can
integrate it with our model and the DE200 masses as a `chiron_compat`
candidate in an hour.

**API semantics recommended** (consistent with Swiss's own documentation):
keep returning Chiron from 675 AD like Swiss, and attach an uncertainty to
pre-1600 positions: 3σ (Horizons-scaled, `docs/MINOR_INTEGRATION.json` ×4):
≤ 1″ in 1200–1600, ~20″ at 1000, ~50″ at 800–900, ~4′ at 760, ~15′ at 700,
~1.2° at 680. Flag Chiron as `uncertain` before 760 AD (3σ > 2′; differences
between engines reach degrees) and say "not meaningful before 720 AD (Saturn
encounter)" in the response metadata rather than refusing the chart. Before
675 AD return `OutOfRange` as now (Swiss errors there too).

## Fit (our own dataset)

Per body, heliocentric (centre 10) or barycentric (centre 0) ICRF positions
are interpolated per record at Chebyshev-Gauss-Lobatto nodes (record ends are
nodes, so positions are continuous across records). Record length (power of 2
days) and degree are chosen to minimise size; then every record is checked at
out-of-sample points (midpoints between nodes + 8 random). Barycentric =
JPL heliocentric + the main pack's Sun (the chain subtracts it again), which
removes the Sun's reflex wobble (Mercury ~10 km/88 d, Venus ~270 km, EMB
~450 km) that makes heliocentric series 3–20× longer for distant bodies;
barycentric series restart at the main pack's DE441→DE442 seam (1549-12-31),
where the Sun jumps by ~0.5 km. Tolerance: the spec allows max(1 km, 1e-5″ ×
geocentric distance); 1 km is the looser bound for all eight bodies (1e-5″ =
1 km only at 138 AU). We fitted to 0.1 km.

| body | centre | record (d) × N | records | KB | max km (fit check) | dense check (a): max / p99 km | ≤ ″ seen from Earth |
|---|---:|---|---:|---:|---:|---|---:|
| Ceres | 0 | 512 × 26 | 1715 | 1070 | 0.089 | 0.054 / 0.021 | 8e-5 |
| Pallas | 0 | 512 × 31 | 1715 | 1276 | 0.054 | 0.060 / 0.015 | 9e-5 |
| Juno | 0 | 512 × 30 | 1715 | 1235 | 0.050 | 0.048 / 0.021 | 7e-5 |
| Vesta | 0 | 512 × 29 | 1715 | 1194 | 0.055 | 0.055 / 0.026 | 7e-5 |
| Chiron | 0 | 1024 × 18 | 286 | 124 | 0.051 | 0.051 / 0.026 | 1e-5 |
| Pholus | 0 | 1024 × 19 | 286 | 130 | 0.049 | 0.049 / 0.014 | 9e-6 |
| Eris | 0 | 1024 × 7 | 286 | 48 | 0.019 | 0.019 / 0.010 | 7e-7 |
| Eros | 10 | 128 × 19 | 2287 | 1043 | 0.071 | 0.070 / 0.034 | 8e-4 |

(a) = `check_minor.py` section a: ~200,000 random epochs per body plus every
record boundary (and the epoch just before it), our pack through the
engine's chaining rule vs the source SPK. Velocity jumps at record
boundaries ≤ 0.024 km/day. The record counts include tail/seam series
(Ceres–Vesta: 4 series each, others 2). Total pack 6.1 MB.

## SPK type 21 reader (`tools/spk.py`)

Own implementation of the modified-difference-array evaluation from the NAIF
SPK Required Reading description. jplephem 2.24 does not support type 21,
so it was validated against Horizons: (1) adjacent difference lines agree at
their shared epoch to ≤ 3e-5 km (Chiron, all 39 segments); (2) Horizons
VECTORS (centre 500@10, ICRF, TDB) at single-epoch requests agree with our
evaluation to ≤ 0.1 m. Note: for multi-epoch TLIST requests Horizons
re-integrates from the first epoch, and its later epochs drift from its own
SPK by up to 26 km (Chiron, 2380); all Horizons comparisons here therefore
use one epoch per request. Layout facts found: the directory holds
floor(N/100) epochs; TL is each line's reference epoch.

## Accuracy vs JPL Horizons (b; geocentric, TT, 5 epochs per era, one per request)

Astrometric ICRF RA/Dec (QUANTITIES 1) separation, max, arcsec:

| body | 1600–1800 | 1800–2100 | 2100–2400 | Swiss (worst era) |
|---|---:|---:|---:|---:|
| Chiron | 0.0000 | 0.0000 | 0.0000 | 2.6 |
| Pholus | 0.0002 | 0.0004 | 0.0000 | 46 |
| Eris | 0.0000 | 0.0001 | 0.0000 | 185 |
| Eros | 0.0035 | 0.0006 | 0.0019 | 9.9 |
| Ceres | 2.85 | 0.09 | 1.53 | 2.5 |
| Pallas | 0.46 | 0.20 | 0.03 | 5.1 |
| Juno | 0.30 | 0.06 | 2.26 | 14.3 |
| Vesta | 1.12 | 0.25 | 0.28 | 7.2 |

For the Horizons-sourced bodies we reproduce Horizons to milliarcseconds (the
residual is the fit plus light-time/frame details). For Ceres–Vesta the
differences are between two JPL orbit solutions: sb441-n16 (2020, the DE441
perturbers) vs Horizons' current SBDB solutions (Ceres JPL#48, Pallas #74,
Juno #144, Vesta #36). They exceed Horizons' formal 3σ (Ceres 0.003″), i.e.
the formal covariance is optimistic about model differences; both are JPL.
Apparent ecliptic longitude of date (QUANTITIES 31) differs from ours by up
to ~1″ for every body including the Horizons-sourced ones: the same ~1″ at
1600 appears for the planets (harness/horizons_report.json), so it is
Horizons' older precession-nutation model for that quantity, not the
ephemeris.

Rust: `astrometric_places_match_horizons` checks 120 Horizons epochs with
per-body limits (1.5 × observed + 0.002″: Chiron/Pholus/Eris 0.003″, Eros
0.008″, Pallas 0.69″, Vesta 1.68″, Juno 3.4″, Ceres 4.3″).

## Coverage and differences vs Swiss (c; black box, `FLG_SWIEPH|FLG_SPEED`, TT, 200 epochs per body and era)

Swiss ranges found (SWIEPH result, not error or Moshier fallback; bisected to
1 day): Ceres, Pallas, Juno, Vesta, Pholus: before −11 … 2400-01-22 (then
Moshier, which has no asteroids); Chiron 675-01-01 … 2400-01-22 (error
before 675); Eris 1500-01-01 … 2102-05-13; Eros 1500-01-01 … 2100-03-04
(the `s136199s.se1` / `se00433s.se1` short files).

| body | 0–600 | 600–1200 | 1200–1600 | 1600–1800 | 1800–2100 | 2100–2400 |
|---|---|---|---|---|---|---|
| Ceres | both; 24″ | both; 20″ | both; 3.3″ | both; 0.6″ | both; 1.2″ | both; 5.5″ |
| Pallas | both; 110″ | both; 16″ | both; 15″ | both; 12″ | both; 8.2″ | both; 1.6″ |
| Juno | both; 248″ | both; 65″ | both; 28″ | both; 4.2″ | both; 4.7″ | both; 19″ |
| Vesta | both; 29″ | both; 15″ | both; 11″ | both; 11″ | both; 1.6″ | both; 3.5″ |
| Chiron | neither | both from 675; median 44″, max 5.6° (before 720) | both; 1.9″ | both; 4.3″ | both; 2.4″ | both; 4.2″ |
| Pholus | both; 127″ | both; 243″ | both; 89″ | both; 71″ | both; 55″ | both; 56″ |
| Eris | neither | neither | both from 1500; 65″ | both; 251″ | both; 41″ | FE only after 2102 |
| Eros | neither | neither | both from 1500; 59″ | both; 71″ | both; 24″ | FE only after 2100 |

Before 1599-12-11 (our integration; `e_integrated_span_vs_swiss`, every 10
years): FE − Swiss separation vs our 3σ. Chiron: 0.2–0.7″ in 1240–1600,
2–11″ in 1060–1180, 19–57″ in 820–1000, 213″ at 759, 1.5° at 699, up to 5.6°
just after 675 — the ratio to our (un-scaled) 3σ is a steady ~4, the same
factor by which Horizons' Chiron 3σ exceeds ours, so Swiss's Chiron is
about one Horizons-sigma away from JPL's orbit until the 720 Saturn encounter
amplifies it. Pholus: median 13″, max 318″, 99 % of epochs inside our 3σ.
Eris: 43–77″ (Swiss's old orbit is outside JPL's 3σ of 12–34″, as it already
is in 1600–1800: 185″ vs Horizons). Eros: 1–10″ (Swiss's orbit; JPL 3σ 0.2″,
and Swiss is 9.9″ off Horizons in 1600–1800). So the differences before 1600
have the size of the Swiss–JPL differences after 1600: they are Swiss's orbit
solutions, not our extrapolation.

Cells: who returns a position; max |Δλ| FE − Swiss (apparent, of date).
Median |Δλ| 1800–2100: Ceres 0.4″, Pallas 0.9″, Juno 0.8″, Vesta 0.3″, Chiron
0.2″, Pholus 3.2″, Eris 1.0″, Eros 0.6″. Swiss's asteroid files are older
orbit solutions: against Horizons Swiss is off by up to 46″ (Pholus), 185″
(Eris), 14″ (Juno), while we match Horizons as above. All Swiss–FE
differences are below Horizons' own 3σ for Pholus/Eris and are orbit-solution
differences, not computation errors. Speed differences are reported in the
JSON (`dlon_speed_arcsec_per_day`).

## Root-finding bounds (d)

Max |dλ/dt| (apparent geocentric) over each range: Pallas 0.69, Juno 0.60,
Vesta 0.54, Ceres 0.46, Eros 0.95, Chiron 0.15, Pholus 0.14, Eris 0.04 °/day;
`search::max_speed` uses 1.0 / 3.5 / 0.2 / 0.1.

## Reproduce

```
.venv/bin/python tools/fetch_minor.py                # sources + docs/MINOR_PROVENANCE.json
cargo build -p fe-core --example integrate_minor
.venv/bin/python tools/integrate_minor.py validate   # model checks vs JPL (~5 min)
.venv/bin/python tools/integrate_minor.py build      # integrations + clones -> data/sources/minor/int_*.npy
.venv/bin/python tools/fit_minor.py                  # data/packs/fepk-minor-v2.bin (~6.5 min)
.venv/bin/python crates/fe-core/tests/fixtures/minor/gen_minor_fixtures.py
cargo test -p fe-core --test minor
PYTHONPATH=<frenchephemeris build with this branch> .venv/bin/python harness/check_minor.py --fixture
```

The harness was run against a debug build of `fe-py` (`cargo build -p fe-py`,
`lib_native.dylib` copied as `frenchephemeris/_native.abi3.so` into a scratch
package dir), since the shared venv's editable install must not be rebuilt.

## Open issues

1. **Chiron vs Swiss before 1000 AD** (section above): equal to JPL's 3σ down
   to 760, degrees before 720; not reducible without Swiss's starting
   elements (historical astorb.dat not retrievable here). API: flag Chiron
   before 760 AD as uncertain.
1b. **Before 1600 = model extrapolation.** Coverage now matches Swiss (Chiron
   from 675, Pholus from year 0, Eris/Eros from 1500), but these positions are
   our integration of today's JPL orbits, with the uncertainties above. Chiron
   before ~760 (3σ > 2′) and especially before the 720 Saturn encounter
   (~1°) is not physically meaningful in either engine; the API may want to
   flag Chiron before 760 and Pholus generally (3σ up to ~10′). Chiron's
   mapped covariance is 4× below Horizons' (unexplained; scaled in the table).
2. **Ceres–Vesta solution choice.** sb441-n16 (one seamless JPL integration,
   0–2400) differs from Horizons' current solutions by ≤ 0.1″ around 2000 but
   up to 2.9″ (Ceres, 1600) / 2.3″ (Juno, 2400). Alternative: Horizons SPKs
   for 1600–2400 plus sb441 before, blended over a transition window to avoid
   a ~3″ step at 1600. Below the API's 0.001° rounding either way.
3. API bug: natal_pkg passes 5145 for Pholus (always a Swiss error).
4. Horizons solutions are updated over time; the ledger pins the solution IDs
   used. Rebuilding later fetches newer solutions (different numbers).
5. Licence/attribution wording to confirm with JPL SSD (see above).
6. The pack manifest's `builder.argv` records the scratch output path used
   for this build (the file was moved into `data/packs/` afterwards).
7. Integration model: 12 of the 16 SB441 asteroids are omitted (no data for
   them in our packs); DE442 instead of DE441 planets over 1550–1600 (Uranus
   ~2,000 km apart) — both far below the orbit uncertainty.
