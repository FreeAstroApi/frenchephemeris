# FrenchEphemeris documentation

Start with the **[user guide](GUIDE.md)**: installation, data packs, the pyswisseph style
`swe_compat` module, the native engine, bodies, houses, ayanamsas, stars, events, threads,
errors and Rust. The reports below hold the method, sources and measurements behind every
feature; each states what was compared, on how many samples, and how to reproduce it.

## Using FrenchEphemeris

| Document | Contents |
|---|---|
| [GUIDE.md](GUIDE.md) | user guide and API reference |
| [FORMAT.md](FORMAT.md) | FEPK pack format: series table, mixed precision, aliases, quantized (mobile) series |
| [ARCHITECTURE.md](ARCHITECTURE.md) | crate layout, data flow from packs to positions, design choices |
| [../CHANGELOG.md](../CHANGELOG.md) | what changed in each release |

## Validation reports

| Report | What it validates |
|---|---|
| [ENGINE_REPORT.md](ENGINE_REPORT.md) | overall accuracy against JPL Horizons, agreement with Swiss Ephemeris by era, speed |
| [TIME_FRAMES_REPORT.md](TIME_FRAMES_REPORT.md) | Delta T (FE-DT/2), precession, nutation, sidereal time, frames, against ERFA and IERS |
| [GEOMETRY_REPORT.md](GEOMETRY_REPORT.md) | houses, ayanamsas (definition and source of each of the 46 modes), nodes, Lilith, apsides |
| [EVENTS_REPORT.md](EVENTS_REPORT.md) | rise, set, transits, twilight, crossings, searches |
| [ECLIPSES_REPORT.md](ECLIPSES_REPORT.md) | solar and lunar eclipses against the NASA Five Millennium Canons, occultations |
| [HELIACAL_REPORT.md](HELIACAL_REPORT.md) | heliacal events and visibility (Schaefer's model) |
| [STARS_REPORT.md](STARS_REPORT.md), [STAR_NAMES.md](STAR_NAMES.md) | fixed star catalogue, names, astrometry |
| [MINOR_BODIES_REPORT.md](MINOR_BODIES_REPORT.md) | Chiron, Pholus, Ceres–Vesta, Eris, Eros: sources, integration, uncertainty |
| [ASTEROIDS_REPORT.md](ASTEROIDS_REPORT.md) | the 326-asteroid bundle and on demand asteroids |
| [FICTITIOUS_BODIES.md](FICTITIOUS_BODIES.md) | Uranian and other hypothetical bodies, the ELP based lunar apsides |
| [EXTENDED_RANGE_REPORT.md](EXTENDED_RANGE_REPORT.md) | the −13198 … +17190 refit packs |
| [DATA_REPORT.md](DATA_REPORT.md) | the main pack: JPL sources, checksums, DE441/DE442 seams, bit identity checks |
| [SWISS_GAP_ANALYSIS.md](SWISS_GAP_ANALYSIS.md) | feature by feature coverage of the classic Swiss Ephemeris API |

## Provenance data

[ASTEROID_LIST.json](ASTEROID_LIST.json), [ASTEROIDS_PROVENANCE.json](ASTEROIDS_PROVENANCE.json),
[MINOR_PROVENANCE.json](MINOR_PROVENANCE.json), [MINOR_INTEGRATION.json](MINOR_INTEGRATION.json)
and [CHIRON_STUDY.json](CHIRON_STUDY.json) record the sources and fits of every minor body.
