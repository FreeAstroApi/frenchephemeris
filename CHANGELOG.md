# Changelog

All notable changes to FrenchEphemeris are listed here. Versions follow
[Semantic Versioning](https://semver.org/).

## 1.1.0 — 2026-10-30

### Changed
- **Delta T FE-DT/2.** IERS Bulletin A predictions for the 13 months after the last observed
  value (`tools/update_delta_t.py`), then extrapolation with the mean rate of the last 10 years
  instead of 12 months. Future positions move: Moon about 0.07″ in 2027, 0.85″ in 2050, 2.8″ in
  2100. Two Swiss-comparison test tolerances at future dates were widened accordingly.

### Added
- `--set modern`: DE442 for 1800–2401, bit-identical to JPL (~69 MB with the minor bodies).
- `--set compact`: lossy 6.5 MB refit for 1800–2400 (`tools/build_ext_pack.py --era compact`),
  used only when no bit-identical pack is installed.
- `frenchephemeris.house_variants`: house variants anchored on a point other than the
  Ascendant — Sun, Moon or a node on the 1st cusp (sunrise, moonrise, nodal charts), Sun or
  Moon on the 7th (sunset, moonset), equal houses with the MC in the middle of the 10th
  (Vehlow-MC), whole-sign houses starting at any sign, and angles only.
- `frenchephemeris.minor_uncertainty`: the 1-sigma error of Chiron, Pholus, Eris and Eros
  positions (arcsec) and an `is_uncertain` flag (3 sigma above the 0.001° display step),
  from JPL's orbit covariance propagated forward and backward (`tools/minor_uncertainty.py`),
  never below JPL Horizons' own 3-sigma at the dates checked. Chiron: about 0.06″ in 2000,
  6″ in 1000, 2.6′ in 720.

## 1.0.0 — 2026-10-30

First public release. FrenchEphemeris has served all of FreeAstroAPI's production traffic
since 0.1.4. From this version on, the public API (`frenchephemeris`, `swe_compat`, the
Rust `fe-core` crate and the FEPK pack format) follows Semantic Versioning: breaking
changes need a new major version.

- Published on PyPI (`pip install frenchephemeris`) with wheels for Linux, macOS and
  Windows.
- The Swiss Ephemeris comparison scripts are kept in a private validation corpus; the
  public test suite runs without them.

No changes to computed positions since 0.1.4.

## 0.1.6 — 2026-10-30

### Licensing
- **Commercial Licence Agreement v1.0** (COMMERCIAL-LICENSE.md) replaces the draft terms:
  B2B licence for the full term of copyright, all updates, OEM use excluded, warranty of
  title with an IP remedy, liability capped at the price paid, French law, Paris courts.
- **Contributor Licence Agreement v1.0** (CONTRIBUTING.md): non-exclusive licence of each
  right with commitment to keep contributions available under an open-source licence.
- SECURITY.md states a support period of at least five years for security updates.

No code changes.

## 0.1.5 — 2026-10-30

### Added
- **`python -m frenchephemeris.packs download`** fetches the data packs from the `packs-v1`
  GitHub release and checks each SHA-256 (`--set standard|dev|extended|asteroids|all`).
  Packs go to `$FRENCHEPHEMERIS_DATA` or a per-user data directory, which
  `open_default()` now searches before a source checkout's `data/packs/`.
- CI (rustfmt, clippy, Rust and Python tests on the development packs) and wheel builds
  for Linux x86_64/aarch64, macOS arm64/x86_64 and Windows x64.
- README quickstart.

### Changed
- Tests, tools and harness scripts no longer contain machine-specific paths: Rust tests
  read `FE_PACK_DIR`, Swiss-oracle checks run only when `SE_EPHE_PATH` is set, and Python
  tests skip when no pack is installed.
- Code formatted with rustfmt and clean under `clippy -D warnings`; minimum Rust 1.88.

## 0.1.4 — 2026-10-30

### Performance
- **Multi-threaded callers.** `calc_ut`, `calc_tt`, `houses_ut` and `houses_ex2` keep the GIL:
  a single position takes a few microseconds, and releasing the GIL for it made concurrent
  threads contend (8 threads ran at a third of the single-thread rate). Batch, search,
  eclipse and other long calls still run without the GIL.
- **`swe_compat`** passes the sidereal mode only for `FLG_SIDEREAL` calls, as an int where
  the binding accepts one (a string made the binding fail an int conversion on every call).
  Call-by-call throughput: ~600,000 tropical / ~400,000 sidereal positions per second, the
  same with 1 or 8 threads.

## 0.1.3 — 2026-10-30

### Fixed
- **Star-anchored ayanamsas with true positions.** True Citra, True Revati, True Pushya,
  True Mula and the galactic-centre modes now place their reference star (or Sgr A*) with
  the body's own aberration and light-deflection settings. With `FLG_TRUEPOS`,
  `FLG_NOABERR` or `FLG_NOGDEFL` the sidereal position was measured against the apparent
  star instead, up to ~19″ off (True Pushya: dasha dates moved by up to two days).

### Added
- `tools/build_ext_pack.py --assemble-partial`: build a pack of just the chunks next to the
  main pack, e.g. `fepk-de441-ext-near-future-v1.bin` (2401–2580, 2 MB).

## 0.1.2 — 2026-10-29

### Fixed
- **Lunar points in equatorial and J2000 coordinates.** Mean and true node, mean and
  osculating apogee and the interpolated apsides ignored `FLG_EQUATORIAL` / `FLG_J2000` and
  returned ecliptic-of-date values (e.g. a declination of 0 for the node). They are now
  rotated to the requested frame like every other body.
- The Python package imports on Python 3.9 again (as declared in `requires-python`).
- **House position on a cusp.** `house_pos` placed a point lying exactly on a cusp (an
  angle or cusp evaluated at its own longitude) just below the cusp because of rounding,
  so `int(house_pos)` gave the previous house (e.g. the Ascendant in house 12). A point
  less than 1 mas before a cusp is now placed on it.
- **Topocentric houses: point on the IC.** `house_pos` with `T` returned 10 (the MC) for a
  point exactly on the IC; it now returns 4.

## 0.1.1 — 2026-10-29

### Fixed
- **Eclipse visibility at a place.** `lun_eclipse_how` and `sol_eclipse_how` now report
  no eclipse when the Moon (lunar) or the Sun (solar) is below the apparent horizon at the
  observer, and zero attributes when there is no eclipse at all. Previously an eclipse
  could be reported as visible from places where the body was below the horizon.
- **Local lunar eclipses at high latitude.** `lun_eclipse_when_loc` now follows the Moon's
  altitude through the eclipse instead of inferring visibility from the next moonrise and
  moonset, which failed when the Moon stayed up or down for long periods (|latitude| ≥ 66°).
- **Local solar eclipses.** `ECL_NONCENTRAL` is set for local partial phases when the
  shadow axis misses the Earth.
- **Planets behind the solar disc.** The gravitational light deflection of a body seen
  behind the Sun is now computed from the solar mass enclosed within the ray's impact
  parameter (solar Model S, Christensen-Dalsgaard et al. 1996): continuous at the limb and
  zero at the centre. Previously it grew without bound near exact superior conjunction,
  shifting positions by up to ~100″ and speeds by up to ~3°/day for a few minutes.

### Changed
- The fixed-star catalogue no longer contains values derived from Gaia data releases
  (their licence carries share-alike terms). Radial velocities and a few cluster positions,
  motions and distances now come from non-Gaia published catalogues; the Hipparcos New
  Reduction is read from CDS catalogue I/311 (values identical).
- IOTA occultation predictions are used for validation only and are no longer distributed.

## 0.1.0 — 2026-10-29

First release.
