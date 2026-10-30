# FrenchEphemeris user guide

This guide covers everything needed to use FrenchEphemeris from Python and Rust: installing
the library and its data packs, the pyswisseph-style `swe_compat` module, the native engine,
bodies, houses, sidereal zodiacs, fixed stars, events, time scales, threads and errors.
Measurements behind every number are in the [reports](README.md).

- [1. Installation](#1-installation)
- [2. Data packs](#2-data-packs)
- [3. The `swe_compat` module](#3-the-swe_compat-module)
- [4. Migrating from pyswisseph](#migrating-from-pyswisseph)
- [5. The native engine](#5-the-native-engine)
- [6. Bodies](#6-bodies)
- [7. Houses](#7-houses)
- [8. Sidereal zodiacs (ayanamsas)](#8-sidereal-zodiacs-ayanamsas)
- [9. Fixed stars](#9-fixed-stars)
- [10. Rise, set, crossings, eclipses, heliacal events](#10-rise-set-crossings-eclipses-heliacal-events)
- [11. Time scales and Delta T](#11-time-scales-and-delta-t)
- [12. Minor-body uncertainty](#12-minor-body-uncertainty)
- [13. Threads, async code and servers](#13-threads-async-code-and-servers)
- [14. Errors](#14-errors)
- [15. Rust](#15-rust)
- [16. Building packs yourself](#16-building-packs-yourself)

## 1. Installation

```bash
pip install frenchephemeris        # CPython 3.9+; wheels for Linux, macOS, Windows
python -m frenchephemeris.packs download
```

From source: Rust 1.88+, then `pip install maturin && maturin build --release -o dist &&
pip install dist/*.whl`.

## 2. Data packs

The engine reads FEPK packs: JPL Chebyshev coefficients, memory-mapped from disk
([format](FORMAT.md)). Download a set with `python -m frenchephemeris.packs download
--set <name>` (each file is checked against a SHA-256 sum; files already present are verified
and skipped):

| Set | Files | Size | Years | Fidelity |
|---|---|---|---|---|
| `standard` (default) | main DE441/DE442 pack, minor bodies, 2401–2581 refit | ~256 MB | −2 … 2581 | bit-identical to JPL |
| `modern` | DE442 1800–2401, minor bodies | ~69 MB | 1800 … 2401 | bit-identical |
| `dev` | DE442 1900–2100, minor bodies | ~27 MB | 1900 … 2100 | bit-identical |
| `extended` | DE441 refits before −2 and after 2401 | ~305 MB | −13198 … 17190 | planets ≤ 0.001″, Moon ≤ 0.1 km |
| `asteroids` | 326-asteroid bundle | ~108 MB | per body | JPL small-body fits |
| `compact` | lossy refit, minor bodies | ~13 MB | 1800 … 2400 | planets ≤ 0.001″, Moon ≤ 0.04″ |
| `mobile` | quantized refits of planets and minor bodies | ~2 MB | 1800 … 2400 | planets ≤ 0.007″, Moon ≤ 0.4″ |

`python -m frenchephemeris.packs list` prints every pack; `path` prints the data directory.

**Where packs are found**, in order: `$FRENCHEPHEMERIS_PACK` (the main pack file), then
`$FRENCHEPHEMERIS_DATA` or the per-user data directory (`~/.local/share/frenchephemeris/packs`
on Linux, `~/Library/Application Support/frenchephemeris/packs` on macOS,
`%LOCALAPPDATA%\frenchephemeris\packs` on Windows), then `data/packs/` in a source
checkout. Bit-identical packs are preferred over lossy ones when several are installed.
`$FRENCHEPHEMERIS_MINOR_PACK` overrides the minor-body pack ("" disables it).

**Docker and servers:** run `python -m frenchephemeris.packs download --dest /app/packs`
during the image build and set `FRENCHEPHEMERIS_DATA=/app/packs`. Packs are memory-mapped,
so all worker processes share the same pages.

## 3. The `swe_compat` module

`import frenchephemeris.swe_compat as swe` exposes pyswisseph's names, flags and return
shapes. Angles are degrees, times Julian days.

| Group | Functions |
|---|---|
| Positions | `calc_ut`, `calc`, `calc_pctr`, `get_planet_name`, `pheno_ut`, `pheno`, `azalt`, `cotrans`, `cotrans_sp`, `refrac` |
| Houses | `houses`, `houses_ex`, `houses_ex2`, `houses_armc`, `houses_armc_ex2`, `house_pos`, `house_name`, `gauquelin_sector` |
| Sidereal | `set_sid_mode`, `get_ayanamsa_ut`, `get_ayanamsa`, `get_ayanamsa_ex_ut`, `get_ayanamsa_ex`, `get_ayanamsa_name` |
| Fixed stars | `fixstar2_ut`, `fixstar2`, `fixstar2_mag`, `fixstar_mag` |
| Rise and set | `rise_trans`, `rise_trans_true_hor` |
| Crossings | `solcross_ut`, `solcross`, `mooncross_ut`, `mooncross`, `mooncross_node_ut`, `mooncross_node`, `helio_cross_ut`, `helio_cross` |
| Eclipses, occultations | `sol_eclipse_when_glob`, `sol_eclipse_when_loc`, `sol_eclipse_where`, `sol_eclipse_how`, `lun_eclipse_when`, `lun_eclipse_when_loc`, `lun_eclipse_how`, `lun_occult_when_glob`, `lun_occult_when_loc`, `lun_occult_where` |
| Heliacal | `heliacal_ut`, `heliacal_pheno_ut`, `vis_limit_mag` |
| Time | `julday`, `revjul`, `date_conversion`, `day_of_week`, `utc_to_jd`, `jdet_to_utc`, `jdut1_to_utc`, `utc_time_zone`, `deltat`, `deltat_ex`, `set_delta_t_userdef`, `set_tid_acc`, `get_tid_acc`, `sidtime`, `sidtime0`, `time_equ`, `lmt_to_lat`, `lat_to_lmt` |
| Settings | `set_topo`, `set_ephe_path` (accepted, ignored), `get_library_path`, `get_current_file_data` |
| Formatting | `split_deg`, `cs2degstr`, `cs2lonlatstr`, `cs2timestr` |

Flags: `FLG_SPEED`, `FLG_SIDEREAL`, `FLG_TOPOCTR`, `FLG_HELCTR`, `FLG_BARYCTR`,
`FLG_EQUATORIAL`, `FLG_XYZ`, `FLG_RADIANS`, `FLG_J2000`, `FLG_ICRS`, `FLG_NONUT`,
`FLG_TRUEPOS`, `FLG_NOABERR`, `FLG_NOGDEFL`, `FLG_ASTROMETRIC`, `FLG_SPEED3`; the ephemeris
selectors `FLG_SWIEPH`, `FLG_JPLEPH`, `FLG_MOSEPH` are accepted and ignored (positions always
come from the packs).

```python
import frenchephemeris.swe_compat as swe

jd = swe.julday(1990, 5, 17, 14.5)
xx, flags = swe.calc_ut(jd, swe.MARS, swe.FLG_SPEED | swe.FLG_EQUATORIAL)   # RA, Dec, dist, speeds
swe.set_topo(2.3522, 48.8566, 35)                                         # lon, lat, metres
xx, _ = swe.calc_ut(jd, swe.MOON, swe.FLG_TOPOCTR | swe.FLG_SPEED)
cusps, ascmc = swe.houses_ex(jd, 48.8566, 2.3522, b"K")                    # Koch
tret = swe.rise_trans(jd, swe.SUN, swe.CALC_RISE, (2.3522, 48.8566, 35))   # next sunrise
ecl = swe.sol_eclipse_when_glob(jd, swe.ECL_ALLTYPES_SOLAR)               # next solar eclipse
```

## Migrating from pyswisseph

| Topic | What changes |
|---|---|
| Import | `import frenchephemeris.swe_compat as swe` instead of `import swisseph as swe` |
| Ephemeris path | `set_ephe_path()` is ignored; packs are found as in [§2](#2-data-packs) |
| Settings scope | `set_sid_mode`, `set_topo` and the Delta T settings apply to the current thread or asyncio task, not the process: set them where you compute |
| Date range | outside the loaded packs, calls raise `swe_compat.Error` instead of falling back to a lower-precision model |
| Results | positions come from JPL DE441/DE442 and differ from Swiss Ephemeris by fractions of an arcsecond for modern dates; after 2025 both engines forecast Delta T and differ by up to a few arcseconds for the Moon by 2100 |
| Ayanamsas | 46 of the 47 modes; J.N. Bhasin and Sheoran follow their authors' publications; Vettius Valens raises an error |
| Speeds | speeds are the true rate of change of the positions |

## 5. The native engine

```python
import numpy as np
import frenchephemeris as fe

eng = fe.open_default(swiss_delta_t=True)    # main pack + every optional pack found
eng = fe.Engine("fepk-de441-de442-v1.bin", True, ["fepk-minor-v2.bin"])   # explicit files
```

| Method | Returns |
|---|---|
| `calc_ut(jd_ut, body, flags=FLG_SPEED, sidereal=None, topo=None)` | `[lon, lat, dist, lon_speed, lat_speed, dist_speed]` |
| `calc_tt(jd_tt, body, flags, sidereal, topo)` | the same, from TT |
| `calc_batch_ut(jd_ut_array, bodies, flags=FLG_SPEED, sidereal=None, topo=None)` | NumPy array `(len(jd), len(bodies), 6)`, computed without the GIL |
| `houses_ut(jd_ut, lat, lon, hsys="P", sidereal=None)` | `(cusps[12], ascmc[8])` |
| `ayanamsa_ut(jd_ut, mode, true_equinox=True)`, `ayanamsa_tt(...)` | degrees |
| `deltat(jd_ut)` | TT − UT1 in days |
| `manifest()` | the packs' JSON manifests |

`sidereal` takes an ayanamsa identifier (`"lahiri"`, `"true_citra"`, …, see
`fe.ayanamsa_idents()`), `topo` a `(lon, lat, alt_m)` tuple. The module also offers the
functional API used by `swe_compat` (`fe.fixstar_ut(engine, name, jd)`,
`fe.sol_eclipse_when_glob(engine, jd)`, `fe.rise_trans(engine, ...)`, …), and
`fe.engine_with_delta_t(engine, user_days=None, tid_acc=None)` for a custom Delta T.

## 6. Bodies

| Id | Body |
|---|---|
| 0–9 | Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto |
| 10, 11 | mean node, true node |
| 12, 13 | mean Lilith (mean apogee), osculating apogee |
| 21, 22 | interpolated apogee and perigee |
| 15, 16 | Chiron, Pholus |
| 17–20 | Ceres, Pallas, Juno, Vesta |
| 40 … | Uranian and other hypothetical bodies ([FICTITIOUS_BODIES.md](FICTITIOUS_BODIES.md)) |
| `AST_OFFSET + n` (10000 + n) | numbered asteroid n: Eris (136199) and Eros (433) in the minor pack, 326 more in the asteroid bundle, any other via `frenchephemeris.asteroids.fetch(n)` from JPL |

## 7. Houses

24 systems, selected by the classic one-letter codes:

| Code | System | Code | System |
|---|---|---|---|
| P | Placidus | K | Koch |
| W | Whole Sign | A / E | Equal (Ascendant) |
| D | Equal (MC) | V | Vehlow equal |
| N | Equal, 1st = 0° Aries | C | Campanus |
| R | Regiomontanus | O | Porphyry |
| B | Alcabitius | T | Topocentric (Polich/Page) |
| M | Morinus | X | Meridian (axial rotation) |
| H | Horizon / azimuth | F | Carter poli-equatorial |
| G | Gauquelin sectors (36) | U | Krusinski-Pisa-Goelzer |
| Y | APC | S | Sripati |
| L | Pullen sinusoidal delta | Q | Pullen sinusoidal ratio |
| I / i | Sunshine (Treindl / Makransky) | J | Savard-A |

`house_pos(armc, lat, eps, (lon, lat), hsys)` gives the house position of any point. At
polar latitudes, Placidus, Koch and related systems raise `PolarHousesError`; catch it and
fall back to another system.

**House variants** (`frenchephemeris.house_variants.houses(jd_ut, lat, lon, variant,
flags=0)`): `sunrise`, `moonrise`, `north_node`, `south_node` (the body on the 1st cusp,
equal houses), `sunset`, `moonset` (on the 7th), `equal_vehlow_mc` (MC in the middle of the
10th), `whole_sign_aries` … `whole_sign_pisces`, `asc_mc` (angles only).

## 8. Sidereal zodiacs (ayanamsas)

46 built-in modes plus user-defined values (`SIDM_USER` with `set_sid_mode(mode, t0,
ayan_t0)`), with every `SIDBIT_*` option (projection on the ecliptic of t0, on the solar
system plane, original precession models). Identifiers for the native engine:

`fagan_bradley`, `lahiri`, `deluce`, `raman`, `ushashashi`, `krishnamurti`, `djwhal_khul`,
`yukteshwar`, `jn_bhasin`, `babyl_kugler1`, `babyl_kugler2`, `babyl_kugler3`,
`babyl_huber`, `babyl_etpsc`, `aldebaran_15tau`, `hipparchos`, `sassanian`, `galcent_0sag`,
`j2000`, `j1900`, `b1950`, `suryasiddhanta`, `suryasiddhanta_msun`, `aryabhata`,
`aryabhata_msun`, `ss_revati`, `ss_citra`, `true_citra`, `true_revati`, `true_pushya`,
`galcent_rgilbrand`, `galequ_iau1958`, `galequ_true`, `galequ_mula`, `galalign_mardyks`,
`true_mula`, `galcent_mula_wilhelm`, `aryabhata_522`, `babyl_britton`, `true_sheoran`,
`galcent_cochrane`, `galequ_fiorenza`, `lahiri_1940`, `lahiri_vp285`, `krishnamurti_vp291`,
`lahiri_icrc` (and `valens_moon`, not yet available). Each mode's definition and source is in
[GEOMETRY_REPORT.md §2](GEOMETRY_REPORT.md).

## 9. Fixed stars

1,175 stars from the Hipparcos New Reduction (van Leeuwen 2007) with IAU WGSN names,
SIMBAD radial velocities and rigorous space motion. Look them up by traditional name,
Bayer/Flamsteed designation or catalogue number: `swe.fixstar2_ut("Spica", jd, 0)` returns
`(position, "Spica,alVir", flags)`. See [STARS_REPORT.md](STARS_REPORT.md) and
[STAR_NAMES.md](STAR_NAMES.md).

## 10. Rise, set, crossings, eclipses, heliacal events

- `rise_trans` (rise, set, upper and lower transit, civil/nautical/astronomical twilight,
  disc centre or limb, with or without refraction, Hindu rising), `rise_trans_true_hor`.
- `solcross_ut` / `mooncross_ut`: the next time the Sun or Moon reaches a longitude
  (ingresses, solar and lunar returns); `mooncross_node_ut`: the next node crossing;
  `helio_cross_ut` for any planet.
- Eclipses: global and local solar eclipses (type, times of contacts, magnitude,
  obscuration, path), lunar eclipses (umbral and penumbral magnitudes, contacts), lunar
  occultations of planets and stars. Validated against the NASA Five Millennium Canons
  ([ECLIPSES_REPORT.md](ECLIPSES_REPORT.md)).
- Heliacal rising and setting of planets, stars and the Moon (Schaefer's visibility model),
  limiting magnitude (`vis_limit_mag`), see [HELIACAL_REPORT.md](HELIACAL_REPORT.md).
- In Rust, `fe_core::search` also finds exact aspects, lunations and stations.

## 11. Time scales and Delta T

Delta T (TT − UT1) follows model FE-DT/2: Stephenson, Morrison, Hohenkerk and Zawilski's
2020 spline before 1955, USNO observed values from 1955, IERS Bulletin A predictions for
the following year, then the last 10 years' rate with the long-term acceleration
([TIME_FRAMES_REPORT.md §2](TIME_FRAMES_REPORT.md)). `open_default(swiss_delta_t=True)`
uses the tidal acceleration that matches the DE441 Moon; `set_delta_t_userdef` and
`set_tid_acc` override it per thread. UTC conversions use the IERS leap-second table.

## 12. Minor-body uncertainty

Before JPL's tabulation, Chiron, Pholus, Eris and Eros come from our own integration of
JPL's orbits, and their positions are uncertain far from the present (Chiron passed 0.18 AU
from Saturn in 720).

```python
from frenchephemeris import minor_uncertainty as mu
mu.sigma_arcsec(swe.CHIRON, jd_ut)        # 1-sigma error of the geocentric direction
mu.is_uncertain(swe.CHIRON, jd_ut)        # True when 3 sigma exceeds 0.001 degree
xx, sigma, uncertain = mu.calc_ut(jd_ut, swe.CHIRON)
```

Chiron's 1σ is about 0.06″ in 2000, 6″ in 1000 and 2.6′ in 720.

## 13. Threads, async code and servers

The engine is immutable and has no global mutable state. In `swe_compat`, the sidereal mode,
observer position and Delta T settings live in context variables, so each thread and each
asyncio task has its own. Throughput is the same from 1 to 8 threads (~500,000–580,000
positions/s); batch calls release the GIL and run in parallel.

```python
from concurrent.futures import ThreadPoolExecutor

def sun_sidereal(jd, mode):
    swe.set_sid_mode(mode)                    # this thread only
    return swe.calc_ut(jd, swe.SUN, swe.FLG_SIDEREAL)[0][0]

with ThreadPoolExecutor(8) as pool:
    a = pool.submit(sun_sidereal, 2451545.0, swe.SIDM_LAHIRI)
    b = pool.submit(sun_sidereal, 2451545.0, swe.SIDM_RAMAN)
```

## 14. Errors

| Exception | When |
|---|---|
| `swe_compat.Error` | any engine error through `swe_compat`, including dates outside the loaded packs and unavailable ayanamsas |
| `frenchephemeris.OutOfRangeError` | a date outside the packs (native API) |
| `frenchephemeris.PolarHousesError` | Placidus, Koch and similar systems at polar latitudes |
| `FileNotFoundError` | no pack installed: run `python -m frenchephemeris.packs download` |

## 15. Rust

```toml
[dependencies]
fe-core = { git = "https://github.com/FreeAstroApi/frenchephemeris", tag = "v1.2.0" }
```

```rust
use fe_core::context::{Body, CalcOptions, Engine};

let eng = Engine::open("packs/fepk-de441-de442-v1.bin".as_ref())?;
let opt = CalcOptions::default();
let mars = eng.calc_ut(Body::Mars, 2451545.0, &opt)?;          // lon, lat, dist and speeds
let eq = fe_core::search::longitude_crossings(&eng, Body::Sun, 0.0, 2451544.5, 2462502.5, &opt)?;
```

The engine is `Send + Sync`: share one `Engine` across threads. `Pack::open_many` chains
several packs. See [ARCHITECTURE.md](ARCHITECTURE.md).

## 16. Building packs yourself

Every published pack can be rebuilt from the public JPL files:

```bash
python tools/build_pack.py                     # DE441+DE442, years -2..2401, bit-identical
python tools/build_pack.py --range 1800 2400   # bit-identical subset (modern set)
python tools/build_ext_pack.py --era past      # refits before -2; --era future after 2401
python tools/fetch_minor.py && python tools/fit_minor.py   # minor bodies
```

The compact and mobile commands are listed at the top of `tools/build_ext_pack.py`;
`docs/FORMAT.md` documents the file format, including the quantized series of the mobile
packs.
