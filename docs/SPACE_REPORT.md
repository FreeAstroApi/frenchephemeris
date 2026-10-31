# Charts away from Earth: the `space` extension

The `space` extension casts a chart for an observer on the surface of another world (Mercury,
Venus, the Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, and Earth itself), at the
centre of a world, or at any point in space. It is optional: the cargo feature `space` in Rust,
`frenchephemeris.ext.space` in Python. The core engine's results do not change.

## What it computes

| Observer | Positions | Horizon (azimuth, elevation) | Houses and angles |
|---|---|---|---|
| Site on a world's surface | yes | yes | yes |
| Centre of a world | yes | no | no |
| Any point in space | yes | no | no |

Positions go through the engine's usual pipeline (light time, solar light deflection,
aberration) with the observer's own barycentric velocity, rotation of the world included.
Earth appears in the chart as a planet; the observer's own world is left out.

## Sources

| Concern | Source |
|---|---|
| World positions | the engine's JPL packs (DE441/DE442); outer planets as system barycentres |
| Orientation (pole, prime meridian) | IAU WGCCRE 2015 rotation models, Archinal et al. 2018, Celest. Mech. Dyn. Astr. 130:22; Earth uses the engine's own Earth rotation |
| Shape | IAU 2015 reference ellipsoids; WGS84 for Earth |
| Mars time | Mars Sol Date and Coordinated Mars Time, Allison and McEwen 2000 (the Mars24 algorithm) |

## Zodiac

- **Earth** (default): Earth's true ecliptic and equinox of date. Signs mean what they mean in
  an ordinary chart.
- **Native**: the world's orbital plane (instantaneous heliocentric orbit normal; for the Moon
  the Earth Moon barycentre's orbit), with 0° at the world's vernal equinox, the ascending node
  of its equator (IAU north pole) on that plane. The angle between the pole and the orbit
  normal is the world's axial tilt; the tests check it against the published values
  (Mars 25.19°, Venus 177.36°, Uranus 97.77°, Pluto 119.61°, ...).

## Houses

The house circle is the world's equator, oriented by its spin axis, so the Ascendant is the
rising point of the zodiac on every world. The house routines receive, in that frame, the
sidereal angle of the local meridian (the ARMC analogue), the latitude from the spin axis and
the obliquity against the chosen zodiac. On worlds spinning against the zodiac's sense
(Venus, Uranus, Pluto: obliquity above 90°), the cusps are computed with the opposite zodiac
pole and mapped back as directions, so their numbering follows the sky's daily motion.

Placidus and Koch are undefined inside a world's polar circle, as on Earth, and raise
`PolarHousesError`. With Uranus tilted 98°, its polar circle reaches 8° from the equator.
Sunshine and Gauquelin sectors are not available away from Earth.

## Validation against JPL Horizons

`crates/fe-core/tests/fixtures/space/gen_space.py` asks JPL Horizons for sites on every world,
four dates from 1950 to 2073 and six or seven targets each: the astrometric ICRF right
ascension and declination (independent check of the observer's position) and the apparent
azimuth and elevation (independent check of the rotation model and of the site on the
ellipsoid). Largest gaps, arcseconds:

| Observer on | Astrometric position | Azimuth and elevation |
|---|---|---|
| Mercury | 0.158 | 0.004 |
| Venus | 0.030 | 0.003 |
| Moon | 0.058 | 11.1 |
| Mars | 0.019 | 0.003 |
| Jupiter | 0.041 | 0.039 |
| Saturn | 0.051 | 0.050 |
| Uranus | 0.355 | 0.355 |
| Neptune | 0.073 | 0.072 |
| Pluto | 0.107 | 0.111 |
| Earth (control) | 0.544 | 73.3 |

- **Uranus, Neptune, Pluto**: the packs use DE442 after 1550, Horizons serves DE441; the two
  JPL editions place Uranus 0.35″ apart, and an observer standing there inherits the gap.
- **Moon**: Horizons orients the Moon with JPL's libration angles; the IAU series approximates
  them to about 100 m on the surface, 11″ on the horizon.
- **Earth**: Horizons rotates the Earth with measured UT1, the engine with its Delta T model;
  a second of time is 15″ on the horizon. Positions and houses of an Earth site otherwise
  match the ordinary chart (houses within 0.02°, positions within 0.01″).

Geometric tests check, on every world, several latitudes and both zodiacs, that the Ascendant
lies on the horizon and in the east, and the MC on the local meridian (outside polar circles,
where the house convention may take the setting point).

## Not included yet

- Moons other than Earth's (Titan, Europa, ...) need satellite ephemerides; true planet centres
  (instead of system barycentres) for Jupiter to Pluto likewise.
- Sidereal zodiacs (ayanamsas) are defined against Earth's equinox and are not applied away
  from Earth.
