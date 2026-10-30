# FrenchEphemeris

A fast, thread-safe ephemeris and astrology-geometry engine in Rust, with Python bindings
(PyO3). Its planetary, lunar and minor-body positions come directly from NASA/JPL
ephemerides (DE441/DE442, sb441-n16, Horizons). It implements the classic astrology API:
houses, ayanamsas, nodes and Lilith, fixed stars, rise/set, eclipses, occultations,
heliacal events and hypothetical bodies. It was written as an independent replacement for
the Swiss Ephemeris in [FreeAstroAPI](https://freeastroapi.com).

Repository: <https://github.com/FreeAstroApi/frenchephemeris>

## Headline numbers

Measured in October 2026 on an Apple M4 (macOS, single machine) against Swiss Ephemeris
2.10.03 via pyswisseph 2.10.3.2 with its standard `.se1` files: accuracy with
FrenchEphemeris 0.1.0 (planetary positions have not changed since), throughput with 0.1.4.
Method, sample sizes and the commands to reproduce every figure are in
[docs/ENGINE_REPORT.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/ENGINE_REPORT.md); results on other hardware will differ.

| | |
|---|---|
| Physical accuracy vs JPL Horizons (astrometric ICRF, years 0–1550, same DE441) | ≤ 0.0005″ |
| After 1550 (DE442 vs Horizons' DE441) | ≤ 0.13″ inner planets/Sun/Moon, ≤ 1.76″ Uranus–Pluto |
| Range | −13198 … +17190 with the optional extended packs |
| Agreement with Swiss Ephemeris 2.10, 1800–2100 charts | Sun ≤ 0.07″, Moon ≤ 0.83″, nodes ≤ 0.04″, no sign or retrograde flips in 8,400 samples |
| Per-call throughput through `swe_compat` (Python, with speeds, 1–8 threads) | ~500,000–580,000 tropical / ~390,000 sidereal positions/s vs Swiss ~310,000 / ~260,000 |
| Batch API, 10 bodies × 10,000 times (Python, numpy out) | ~4.5 million positions/s |
| Rust kernel, 10-body chart with speeds, one core | 10 µs |
| Eclipses vs the NASA Five Millennium Canons, years 0–2400 | 0 missed, 0 spurious, 0 type mismatches |

No global mutable state: the engine is immutable and `Send + Sync`, and batch calls release
the GIL and run across cores.

## Layout

```
crates/fe-core   Rust kernel
crates/fe-py     PyO3 bindings -> Python package `frenchephemeris` (python/)
tools/           pack builders and table generators (Python)
harness/         differential checks vs JPL DE files and Horizons
docs/            design, data format and validation reports
```

See [docs/ARCHITECTURE.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/ARCHITECTURE.md) and [docs/FORMAT.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/docs/FORMAT.md).

Comparisons with the Swiss Ephemeris were run with pyswisseph as a black-box oracle. Those
scripts and their outputs (fixture generators `gen_*.py`, `harness/compare_swiss.py` and
similar) are kept in a private validation corpus and are not part of this repository; the
reports in `docs/` describe how each comparison was made.

## Quickstart (Python)

```bash
pip install maturin && maturin build --release -o dist && pip install dist/*.whl
python -m frenchephemeris.packs download          # standard packs, ~256 MB, checksummed
python -m frenchephemeris.packs download --set dev   # or: 1900–2100 only, ~27 MB
```

```python
import numpy as np
import frenchephemeris as fe
import frenchephemeris.swe_compat as swe   # pyswisseph-style calls

jd = swe.julday(1990, 5, 17, 14.5)
(lon, lat, dist, speed, *_), _ = swe.calc_ut(jd, swe.MOON, swe.FLG_SPEED)
cusps, ascmc = swe.houses(jd, 48.8566, 2.3522, b"P")
swe.set_sid_mode(swe.SIDM_LAHIRI)            # per thread / async task, not process-wide
(sid_lon, *_), _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SIDEREAL)

eng = fe.open_default()                      # native engine
pos = eng.calc_batch_ut(jd + np.arange(1000.0), list(range(10)))   # shape (1000, 10, 6)
```

`python -m frenchephemeris.packs list` shows every pack and set (`standard`, `dev`,
`extended` for −13200 … +17191, `asteroids`, `all`). Packs go to `$FRENCHEPHEMERIS_DATA`
or a per-user data directory (`python -m frenchephemeris.packs path`).

## Building

```bash
cargo build --release --workspace
FE_PACK_DIR=~/path/to/packs cargo test --release --workspace --exclude fe-py   # pack tests skip without packs
pytest tests                                # after installing the wheel
```

### Building the packs yourself

The published packs can be rebuilt from the public JPL files into `data/packs/`
(Python 3.9+, `numpy`):

```bash
python tools/build_pack.py                    # DE441+DE442, years -1..2400 (~247 MB)
python tools/build_pack.py --range 1900 2100  # small DE442-only development pack
python tools/build_ext_pack.py --era past     # optional: -13200 .. -1
python tools/build_ext_pack.py --era future   # optional: 2401 .. 17191
python tools/fetch_minor.py && python tools/fit_minor.py   # Chiron, Pholus, Ceres–Vesta, Eris, Eros
```

`build_pack.py` expects `de442.bsp` in `data/sources/` and reads DE441 remotely with HTTP
range requests. `frenchephemeris.open_default()` looks for packs in `$FRENCHEPHEMERIS_DATA`
(or the per-user data directory), then `data/packs/` in a source checkout;
`$FRENCHEPHEMERIS_PACK` names the main pack directly. Rust tests read `$FE_PACK_DIR`.

## Licensing

FrenchEphemeris is **dual-licensed**:

- **AGPL-3.0-only** ([LICENSE](https://github.com/FreeAstroApi/frenchephemeris/blob/main/LICENSE)) — free, with the AGPL's obligations, including
  offering the source of modified versions to users who interact with them over a network
  (§ 13); or
- a **commercial licence** — USD 200 one-time, per company, for the full term of
  copyright with all updates,
  including use behind your own APIs and apps for unlimited end users — for closed-source
  or otherwise non-AGPL use; see
  [COMMERCIAL-LICENSE.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/COMMERCIAL-LICENSE.md), contact
  [contact@freeastroapi.com](mailto:contact@freeastroapi.com).

SPDX: `AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial`.
Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI.
Contributions require the contributor licence agreement in [CONTRIBUTING.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/CONTRIBUTING.md).

## Data sources and attribution

Third-party data and their terms are listed in [NOTICE](https://github.com/FreeAstroApi/frenchephemeris/blob/main/NOTICE): NASA/JPL DE441/DE442 and
small-body ephemerides (Park et al. 2021), ERFA (BSD-3) coefficient tables, the IAU WGSN
star names (CC BY 4.0), Hipparcos (ESA; van Leeuwen 2007), SIMBAD (CDS), IERS EOP and leap
seconds, USNO and Stephenson–Morrison–Hohenkerk Delta T, NASA eclipse catalogues
(F. Espenak), the ELP 2000-82B and ELP/MPP02 lunar theories (Chapront-Touzé, Chapront,
Francou), and the published astronomy and astrology
definitions behind each feature. Per-constant provenance is in the `docs/*REPORT*.md` files.

The Swiss Ephemeris was used only as a black-box oracle through its public API to measure
compatibility. No Swiss code, data files or outputs are part of this repository; the
oracle comparison fixtures are kept in a private corpus (`FE_ORACLE_DIR`), and the tests
that need them are skipped without it.

## Trademarks

"Swiss Ephemeris" is a trademark of its owner. FrenchEphemeris is an independent
implementation; it is not affiliated with, sponsored or endorsed by the owner or the
authors of the Swiss Ephemeris. The
name is used here only to describe compatibility and to report measured comparisons. The
`swe_compat` Python module is an interface-compatibility layer that accepts the same
function names and argument shapes; it contains no Swiss Ephemeris code.

"FrenchEphemeris" and "FreeAstroAPI" are trade names of Gabriel Naulleau EI; see
[TRADEMARKS.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/TRADEMARKS.md) for how forks and derived works may refer to them.

## Disclaimer

FrenchEphemeris is provided "as is", without warranty of any kind, as stated in the
AGPL-3.0 (sections 15–16) and in the commercial licence. Results are intended for
astrology, education and general astronomy. It is **not** designed or certified for
navigation, surveying, spacecraft operations, legal or civil timekeeping, or any
safety-critical use. Positions far from the present (before about 1600 and after about
2100) inherit the uncertainties of the underlying models (Delta T, minor-body orbits),
which are documented in the reports.

## Security

See [SECURITY.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/SECURITY.md) to report a vulnerability. Community rules are in
[CODE_OF_CONDUCT.md](https://github.com/FreeAstroApi/frenchephemeris/blob/main/CODE_OF_CONDUCT.md).
