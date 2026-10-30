# French Ephemeris

[![PyPI](https://img.shields.io/pypi/v/frenchephemeris.svg)](https://pypi.org/project/frenchephemeris/)
[![Python](https://img.shields.io/pypi/pyversions/frenchephemeris.svg)](https://pypi.org/project/frenchephemeris/)
[![CI](https://github.com/FreeAstroApi/frenchephemeris/actions/workflows/ci.yml/badge.svg)](https://github.com/FreeAstroApi/frenchephemeris/actions/workflows/ci.yml)
[![License: AGPL-3.0 or commercial](https://img.shields.io/badge/license-AGPL--3.0%20or%20commercial-blue.svg)](https://github.com/FreeAstroApi/frenchephemeris/blob/main/COMMERCIAL-LICENSE.md)

**FrenchEphemeris** is a fast, thread safe **ephemeris and astrology calculation library**
written in Rust, with Python bindings. It computes planetary, lunar and asteroid positions
from **NASA/JPL DE441/DE442**, the same data JPL Horizons serves, together with house
cusps, ayanamsas, lunar nodes and Lilith, fixed stars, rise and set times, solar and lunar
eclipses, occultations and heliacal events. A **pyswisseph style API** (`swe_compat`) lets
existing Swiss Ephemeris code switch with one import change.

It is the engine behind every calculation of [FreeAstroAPI](https://www.freeastroapi.com),
and an independent, clean room alternative to the Swiss Ephemeris.

```bash
pip install frenchephemeris          # works immediately: years 0 to 2400 are built in
python -m frenchephemeris.packs download   # optional: JPL's exact numbers for 1800 to 2401
```

- **Accurate:** within 0.0005″ of JPL Horizons for years 0–1550; every eclipse of the NASA
  Five Millennium Canons found, years 0–2400.
- **Fast:** ~500,000 positions/s per call from Python, ~4.5 million/s in batches.
- **Thread safe:** no global state; settings are per thread and per asyncio task.
- **Complete:** 24 house systems, 46 ayanamsas, 1,175 fixed stars, Chiron, Pholus, main
  asteroids, any numbered asteroid, years −13198 … +17190.
- **No setup:** the package ships a compact data set for years 0 to 2400; download JPL's
  exact numbers when you want them, or a 2 MB set for mobile apps.

## Contents

- [Installation](#installation)
- [Quickstart](#quickstart)
- [Data packs](#data-packs)
- [Features](#features)
- [Accuracy and validation](#accuracy-and-validation)
- [Performance](#performance)
- [FrenchEphemeris and the Swiss Ephemeris](#frenchephemeris-and-the-swiss-ephemeris)
- [Documentation](#documentation)
- [FAQ](#faq)
- [Licensing](#licensing) · [Data sources](#data-sources-and-attribution) · [Trademarks](#trademarks) · [Disclaimer](#disclaimer) · [Security](#security)

## Installation

Prebuilt wheels (CPython 3.9 and later, one abi3 wheel per platform):

| Platform | Architectures |
|---|---|
| Linux (manylinux 2.28) | x86_64, aarch64 |
| macOS | Apple Silicon, Intel |
| Windows | x64 |

```bash
pip install frenchephemeris
python -m frenchephemeris.packs download          # exact JPL data for 1800 to 2401, ~69 MB, SHA-256 checked
```

The engine reads **data packs** (JPL coefficients in the FEPK format). The wheel includes a
compact lossy set for years 0 to 2400 (7.7 MB), so the library works right after
`pip install`. `packs download` adds JPL's exact numbers from the
[packs v1 release](https://github.com/FreeAstroApi/frenchephemeris/releases/tag/packs-v1);
they are then used wherever they apply, and the built in set covers the other dates. See
[Data packs](#data-packs) for every set. To build from source you need Rust 1.88+ and
[maturin](https://www.maturin.rs): `pip install maturin && maturin build --release`.

## Quickstart

### Python, pyswisseph style

```python
import frenchephemeris.swe_compat as swe

jd = swe.julday(1990, 5, 17, 14.5)                       # UT

(lon, lat, dist, speed, *_), _ = swe.calc_ut(jd, swe.MOON, swe.FLG_SPEED)
cusps, ascmc = swe.houses(jd, 48.8566, 2.3522, b"P")     # Placidus, Paris; ascmc[0] = Ascendant

swe.set_sid_mode(swe.SIDM_LAHIRI)                         # this thread / async task only
(sid_lon, *_), _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SIDEREAL)

xx, name, _ = swe.fixstar2_ut("Spica", jd, 0)             # fixed stars by name
```

Most code written for `pyswisseph` runs after replacing `import swisseph as swe` with
`import frenchephemeris.swe_compat as swe`. See the
[migration notes](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/GUIDE.md#migrating-from-pyswisseph).

### Python, native engine and NumPy batches

```python
import numpy as np
import frenchephemeris as fe

eng = fe.open_default()                                   # every installed pack
times = 2448029.1 + np.arange(10_000.0)                   # JD UT
pos = eng.calc_batch_ut(times, list(range(10)))           # Sun .. Pluto, shape (10000, 10, 6)
sid = eng.calc_batch_ut(times, [1], sidereal="lahiri", topo=(2.3522, 48.8566, 35.0))
```

Batches run without the GIL and return NumPy arrays: longitude, latitude, distance and
their speeds.

### Rust

```toml
[dependencies]
fe-core = { git = "https://github.com/FreeAstroApi/frenchephemeris", tag = "v1.3.0" }
```

```rust
use fe_core::context::{Body, CalcOptions, Engine};

fn main() -> fe_core::Result<()> {
    let eng = Engine::open("packs/fepk-de441-de442-v1.bin".as_ref())?;
    let moon = eng.calc_ut(Body::Moon, 2448029.104166667, &CalcOptions::default())?;
    println!("Moon {:.6}°, {:.4}°/day", moon.lon, moon.lon_speed);
    Ok(())
}
```

## Data packs

`python -m frenchephemeris.packs download --set <name>` fetches a set and checks every file
against published SHA-256 sums. `list` shows every pack, `path` prints the data directory.

| Set | Size | Years | Fidelity | Use |
|---|---|---|---|---|
| built in (`lite`) | 7.7 MB, in the wheel | 0 … 2400 | lossy, planets ≤ 0.007″, Moon ≤ 0.4″ | works with no download |
| `modern` (the default download) | ~69 MB | 1800 … 2401 | JPL coefficients, bit identical | exact results for modern charts |
| `full` | ~256 MB | −2 … 2581 | bit identical | exact results for every era |
| `dev` | ~27 MB | 1900 … 2100 | bit identical | development, CI |
| `extended` | ~305 MB | −13198 … 17190 | refit, ≤ 0.001″ planets | ancient and far future dates |
| `asteroids` | ~108 MB | varies | refit | 326 asteroids used in astrology |
| `compact` | ~13 MB | 1800 … 2400 | lossy, planets ≤ 0.001″, Moon ≤ 0.04″ | small installs |
| `mobile` | ~2 MB | 1800 … 2400 | lossy, planets ≤ 0.007″, Moon ≤ 0.4″ | mobile and embedded apps |

Packs go to `$FRENCHEPHEMERIS_DATA` if set, otherwise a per user data directory
(`~/.local/share/frenchephemeris/packs`, `~/Library/Application Support/frenchephemeris/packs`
or `%LOCALAPPDATA%\frenchephemeris\packs`). `$FRENCHEPHEMERIS_PACK` names the main pack
directly. Packs are memory mapped, so worker processes share them. Exact packs are always
preferred; the built in set answers only for dates or bodies they do not cover
(`FRENCHEPHEMERIS_NO_FALLBACK=1` turns that off). Outside every loaded pack, calls raise an
error instead of returning an approximation.

## Features

| Area | What is included |
|---|---|
| Bodies | Sun, Moon, Mercury … Pluto; mean and true lunar node; mean, osculating and interpolated Lilith (apogee and perigee); Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Eros; a 326-asteroid bundle; any numbered asteroid on demand from JPL; Uranian and other hypothetical bodies |
| Coordinates | apparent, true or astrometric; geocentric, topocentric, heliocentric, barycentric, planet centric; ecliptic or equatorial; of date, J2000 or ICRS; Cartesian; speeds |
| Houses | 24 systems (Placidus, Koch, Whole Sign, Equal, Campanus, Regiomontanus, Porphyry, Alcabitius, Topocentric, Morinus, Meridian, Gauquelin sectors, Krusinski, APC, Sunshine, …), house position of a point, plus variants anchored on the Sun, Moon or nodes (`house_variants`) |
| Sidereal | 46 ayanamsas (Lahiri, Raman, Krishnamurti, Fagan/Bradley, True Citra, True Pushya, J.N. Bhasin, Sheoran, galactic centre and Babylonian modes, …), user defined, all `SIDBIT` options |
| Stars | 1,175 fixed stars from the Hipparcos New Reduction with IAU names, rigorous space motion |
| Events | rise, set and meridian transit, twilight; Sun, Moon and heliocentric crossings (ingresses, returns); solar and lunar eclipses (global and local); lunar occultations; heliacal rising and setting; in the Rust core also exact aspects, lunations and stations (`fe_core::search`) |
| Time | Delta T from observations, IERS predictions and a published long term model; UTC with leap seconds; TT, TDB, UT1; calendar conversions |
| Uncertainty | 1σ position error and an "uncertain" flag for Chiron, Pholus, Eris and Eros from JPL's orbit covariance (`minor_uncertainty`) |

The complete function list is in the [user guide](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/GUIDE.md).

## Accuracy and validation

| Check | Result |
|---|---|
| Positions vs JPL Horizons (astrometric, years 0–1550, same DE441) | ≤ 0.0005″ |
| After 1550 (DE442 vs Horizons' DE441) | ≤ 0.13″ inner planets, Sun, Moon; ≤ 1.76″ Uranus–Pluto |
| Agreement with Swiss Ephemeris 2.10, 1800–2100 charts | Sun ≤ 0.07″, Moon ≤ 0.83″, nodes ≤ 0.04″; no sign or retrograde changes in 8,400 samples |
| 10,000 random births through a full astrology API (Western, Vedic, dasha, KP, panchang) | no sign, nakshatra, pada or dasha lord changed for 1900–2025 |
| Eclipses vs the NASA Five Millennium Canons, years 0–2400 | 0 missed, 0 spurious, 0 type mismatches |
| Eclipse contacts vs NASA's Besselian elements | 0.4 s median (Swiss: 4–6 s) |

Every figure has its method, sample and reproduction command in the
[engine report](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/ENGINE_REPORT.md)
and the topic reports listed under [Documentation](#documentation).

## Performance

Measured on an Apple M4 with `tools/bench_throughput.py --swiss` (Python, 10 bodies, random
dates 1900–2100, with speeds):

| Positions per second | FrenchEphemeris | Swiss Ephemeris 2.10 |
|---|---|---|
| One body per call, 1 thread | ~500,000 | ~313,000 |
| One body per call, 8 threads | ~580,000 | ~311,000 |
| Sidereal (Lahiri), 1 thread | ~386,000 | ~262,000 |
| Batch API, 10 bodies × 10,000 times | ~4,460,000 |, |

The engine is immutable and `Send + Sync`; single positions keep the GIL (a few
microseconds of work), batches and searches release it.

## FrenchEphemeris and the Swiss Ephemeris

FrenchEphemeris is an independent implementation, written from published sources only.
Where an independent reference exists (JPL
Horizons, NASA eclipse canons, IERS Delta T), FrenchEphemeris is the closer of the two;
where none exists (future Earth rotation, Mean Lilith conventions), the reports say so. Data
files are a choice: Swiss ships small lossy files, FrenchEphemeris ships JPL's numbers
bit identical, plus a lossy `mobile` set that is smaller than Swiss's files for the same
years. The full point by point comparison is in
[docs/ENGINE_REPORT.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/ENGINE_REPORT.md)
and [docs/SWISS_GAP_ANALYSIS.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/SWISS_GAP_ANALYSIS.md).

## Documentation

| Document | Contents |
|---|---|
| [User guide](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/GUIDE.md) | installation, packs, `swe_compat`, native engine, bodies, houses, sidereal, stars, events, threads, errors, Rust |
| [Documentation index](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/README.md) | every report, with what it validates |
| [Online docs](https://www.freeastroapi.com/docs/utilities/frenchephemeris) | the same quickstart on freeastroapi.com |
| [Engine report](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/ENGINE_REPORT.md) | accuracy, Swiss compatibility, speed |
| [Pack format](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/FORMAT.md) | the FEPK file format |
| [Changelog](https://github.com/FreeAstroApi/frenchephemeris/blob/main/CHANGELOG.md) | releases |

## FAQ

**What is FrenchEphemeris?**
An open source ephemeris library for Python and Rust that computes the positions of the Sun,
Moon, planets, asteroids and stars from NASA/JPL data, with the geometry astrology needs:
houses, ayanamsas, nodes, crossings, eclipses and rise/set times.

**Is it a replacement for pyswisseph and the Swiss Ephemeris?**
Yes for most code. `frenchephemeris.swe_compat` accepts pyswisseph's function names, flags
and return shapes (`calc_ut`, `houses`, `set_sid_mode`, `fixstar2_ut`, `rise_trans`,
eclipse functions and more). Results agree with Swiss Ephemeris 2.10 to under 1″ for modern
charts. It is an independent implementation, not affiliated with the Swiss Ephemeris.

**How accurate is it?**
It reads JPL's DE441/DE442 coefficients directly: within 0.0005″ of JPL Horizons for years
0–1550, and within 0.13″ for the Sun, Moon and inner planets after 1550 (where it uses the
newer DE442). Charts display 0.001° = 3.6″, far coarser than either figure.

**Which data files do I need, and how large are they?**
None to start: the package includes a compact set for years 0 to 2400 (planets within
0.007″ of JPL, Moon within 0.4″). `python -m frenchephemeris.packs download` adds JPL's exact
numbers for 1800 to 2401 (69 MB), `--set full` for −2 to 2581 (256 MB), and `--set mobile`
(2 MB) is a lossy set for apps.

**Can I use it in a mobile or embedded app?**
Yes. The `mobile` set is 2.03 MB for 1800–2400, planets, Moon and minor bodies included,
with planets within 0.007″ and the Moon within 0.4″ of JPL. The Rust crate depends only on
`memmap2` and `thiserror`.

**Is it thread safe? Can I use it in FastAPI, Django or Celery?**
Yes. There is no global mutable state; the sidereal mode, observer position and Delta T
settings in `swe_compat` are stored per thread and per asyncio task, so concurrent requests
never see each other's settings.

**Does it support Vedic and KP astrology?**
Yes: 46 ayanamsas including Lahiri, Raman, Krishnamurti and True Chitra, sidereal houses,
nakshatra ready longitudes, the true node, and sunrise times for panchang. FreeAstroAPI runs
its Vedic, KP and dasha endpoints on it.

**Which date range is covered?**
Years 0 to 2400 out of the box, −2 to 2581 with `--set full`, and −13198 to +17190 with the
extended packs. Outside the loaded range, calls raise an error.

**Is it free for commercial use?**
Yes under the AGPL-3.0, provided you meet its conditions (including publishing the source of
modified versions you offer over a network). For closed source products, the commercial
licence costs USD 200 once per company, with every future version and unlimited end users.

**Who maintains it, and is it used in production?**
Gabriel Naulleau (FreeAstroAPI). It has computed every FreeAstroAPI response since
October 2026.

**How do I report a bug or a security issue?**
Open an issue for wrong results or crashes; report vulnerabilities privately as described in
[SECURITY.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/SECURITY.md).

## Licensing

FrenchEphemeris is **dual licensed**:

- **AGPL-3.0-only** ([LICENSE](https://github.com/FreeAstroApi/frenchephemeris/blob/main/LICENSE)), free, with the AGPL's obligations, including
  offering the source of modified versions to users who interact with them over a network
  (§ 13); or
- a **commercial licence**, USD 200 one time, per company, for the full term of
  copyright with all updates,
  including use behind your own APIs and apps for unlimited end users, for closed source
  or otherwise non AGPL use; see
  [the commercial licence terms](https://github.com/FreeAstroApi/frenchephemeris/blob/main/COMMERCIAL-LICENSE.md), contact
  [contact@freeastroapi.com](mailto:contact@freeastroapi.com).

SPDX: `AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial`.
Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI.
Contributions require the contributor licence agreement in [CONTRIBUTING.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/CONTRIBUTING.md).

## Data sources and attribution

Third party data and their terms are listed in [NOTICE](https://github.com/FreeAstroApi/frenchephemeris/blob/main/NOTICE): NASA/JPL DE441/DE442 and
small body ephemerides (Park et al. 2021), ERFA (BSD-3) coefficient tables, the IAU WGSN
star names (CC BY 4.0), Hipparcos (ESA; van Leeuwen 2007), SIMBAD (CDS), IERS EOP and leap
seconds, USNO and Stephenson–Morrison–Hohenkerk Delta T, NASA eclipse catalogues
(F. Espenak), the ELP 2000-82B and ELP/MPP02 lunar theories (Chapront-Touzé, Chapront,
Francou), and the published astronomy and astrology
definitions behind each feature. Per constant provenance is in the `docs/*REPORT*.md` files.

The Swiss Ephemeris was used only as a black box oracle through its public API to measure
compatibility. No Swiss code, data files or outputs are part of this repository; the
oracle comparison fixtures are kept in a private corpus (`FE_ORACLE_DIR`), and the tests
that need them are skipped without it.

## Trademarks

"Swiss Ephemeris" is a trademark of its owner. FrenchEphemeris is an independent
implementation; it is not affiliated with, sponsored or endorsed by the owner or the
authors of the Swiss Ephemeris. The
name is used here only to describe compatibility and to report measured comparisons. The
`swe_compat` Python module is an interface compatibility layer that accepts the same
function names and argument shapes; it contains no Swiss Ephemeris code.

"FrenchEphemeris" and "FreeAstroAPI" are trade names of Gabriel Naulleau EI; see
[TRADEMARKS.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/TRADEMARKS.md) for how forks and derived works may refer to them.

## Disclaimer

FrenchEphemeris is provided "as is", without warranty of any kind, as stated in the
AGPL-3.0 (sections 15–16) and in the commercial licence. Results are intended for
astrology, education and general astronomy. It is **not** designed or certified for
navigation, surveying, spacecraft operations, legal or civil timekeeping, or any
safety critical use. Positions far from the present (before about 1600 and after about
2100) inherit the uncertainties of the underlying models (Delta T, minor body orbits),
which are documented in the reports.

## Security

See [SECURITY.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/SECURITY.md) to report a vulnerability. Community rules are in
[CODE_OF_CONDUCT.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/CODE_OF_CONDUCT.md).
