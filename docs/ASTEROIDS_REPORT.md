# Asteroids beyond the minor pack

Two ways to get any numbered asteroid, both JPL-only (clean room: Swiss is a
black-box oracle in the harness, never a data source):

1. **Bundle** `fepk-asteroids-v1.bin`: the asteroids astrologers actually use
   (list and reasons: `docs/ASTEROID_LIST.json`, built by `tools/asteroid_list.py`).
2. **On demand**: `frenchephemeris.asteroids.fetch(mpc_number)` fetches any
   numbered asteroid from JPL Horizons, fits it and caches a small pack that
   `open_default()` loads.

The 8 bodies of the minor pack (Ceres, Pallas, Juno, Vesta, Chiron, Pholus,
Eris, Eros) are not in the bundle; they keep their ids and their longer,
integration-extended coverage (docs/MINOR_BODIES_REPORT.md).

Status: `fepk-asteroids-v1.bin` = 107,936,784 bytes, **325 of 326 listed
bodies**, 64,797 series, all covering JD 2305426.5–2598072.5
(1599-12-11 – 2401-03-07). Not included: 2101 Adonis (see Open issues).
Per body: median 270 kB (main belt 190–330 kB, TNOs 55 kB, eccentric
NEAs up to 2.5 MB: Icarus, Talos, Phaethon, Hathor, Ra-Shalom, Khufu).

## Ids

Swiss convention: engine id = `AST_OFFSET` (10000) + MPC number, for TNOs and
dwarf planets too (Sedna 100377, Haumea 146108, Makemake 146472, Gonggong
235088). Rust: `Body::from_id(10000 + n)` → `Body::Asteroid(n)` (the 8 minor
bodies map to their own variants; 10001–10004, 12060, 15145, 10433, 146199
are aliases of them). NAIF id in packs: 2000000 + n (20000000 + n for
n ≥ 1,000,000). Python: `fe.asteroid_id(n)`, `fe.AST_OFFSET`.

An asteroid whose pack is not loaded is `UnknownBody` (NAIF id), never a wrong
value; outside its range `OutOfRange`.

## Selection (`docs/ASTEROID_LIST.json`, 326 bodies)

Every entry lists why it is included (several reasons possible):

| reason | count | what |
|---|---:|---|
| `swiss_docs: <section>` | 99 | listed in the Swiss Ephemeris documentation's asteroid introduction (main-belt giants, Aten/Apollo/Amor families, Hirayama families, Trojans, Jupiter-crossers, centaurs, "transplutonians" Isis/Cupido/Poseidon/Vulcano/Zeus/Apollo/Admete/Lilith, love/wisdom asteroids, Beatles, composers, frequent names) |
| `requested` | 40 | the astrological asteroids requested for FreeAstroAPI (Hygiea, Astraea, Psyche, Lilith, Sappho, Amor, Hidalgo, Nessus, Chariklo, Asbolus, Sedna, Quaoar, Orcus, Haumea, Makemake, Varuna, Ixion, Gonggong, Salacia, Varda, Hebe, Iris, Flora, Metis, Kassandra, Urania, Lucifer, Bacchus, Apollo, Icarus, Toro, Pandora, Ariadne, Eurydike, Siva, Medusa, Persephone, Hekate, Nemesis, Dembowska) |
| `classical` | 196 | MPC numbers 5–200: the classical, mostly mythologically named asteroids that asteroid-astrology software ships by default |
| `tno_centaur` | 28 | further dwarf planets, large TNOs and centaurs (Huya, Ceto, Typhon, Chaos, Okyrhoe, Bienor, Thereus, Elatus, Hylonome, Echeclus, Crantor, Dziewanna, Sila-Nunam, Albion, Altjira, …) |
| `neo_myth` | 27 | mythologically named near-Earth asteroids (Adonis, Hephaistos, Phaethon, Daedalus, Quetzalcoatl, Toutatis, Nereus, Apophis, Geographos, …) |

## Data and fit

Per body (`tools/build_asteroids.py fetch`, one at a time; the SPK is deleted
right after its fit):

* JPL Horizons small-body SPK (type 21, centre Sun), requested for
  1599-12-11 .. 2401-04-01 — the longest span Horizons serves (it has no
  small-body ephemeris before 1599-12-11). Narrower fallbacks (1700–2300,
  1800–2200, 1900–2100) are tried only if Horizons refuses.
* Fit (`frenchephemeris/asteroids.py`, the same code serves on-demand
  requests): barycentric ICRF positions (JPL heliocentric + the main pack's
  Sun), Chebyshev-Gauss-Lobatto interpolation per record; record length and
  degree chosen per body for the smallest pack over a small grid (base length
  from the orbital period, × ½…16, degree 8–28), and records failing the
  tolerance halved locally (near perihelion of eccentric NEAs). Tolerance
  0.5 km on out-of-sample points (node midpoints + 4 random per record),
  within the 1 km requirement. The SPK type-21 reader is vectorised over
  epochs (15× faster than tools/spk.py; agrees with it to 3e-8 km).
* Coverage per body = JD 2305426.5 .. 2598072.5 (1599-12-11 .. 2401-03-07),
  explicit; outside → `OutOfRange`. No integration extension before 1600
  (unlike Chiron/Pholus/Eris/Eros): JPL's own data stop there and the
  bundle is meant to be JPL-faithful; `tools/integrate_minor.py` could extend
  individual bodies where an uncertainty study justifies it.
* Provenance: `docs/ASTEROIDS_PROVENANCE.json` (per body: JPL solution id,
  solution date, observation arc, SPK sha256, request, retrieval time) and the
  pack manifest (adds fit statistics and 8 JPL reference states per body,
  which become `crates/fe-core/tests/fixtures/asteroids/asteroid_states.json`).

## Engine changes

* `Pack`: any small-body id (≥ 2000000) present in a loaded pack becomes a
  chain node at `Pack::open_many` time (sorted id table, binary search; no
  allocation at evaluation). Test `many_small_bodies_chain_dynamically`
  (300 bodies).
* `Body::Asteroid(u32)`; `Body::id()` now matches explicitly (`#[repr(i32)]`).
  `search::max_speed(Asteroid)` = 20°/day (orbit class unknown; NEAs near
  Earth), so root scans are conservative for asteroids.
* Python: `frenchephemeris.asteroids` (numpy only, no native import),
  `fe.asteroid_id`, `fe.default_extra_packs()`; `open_default()` loads the
  minor pack, the bundle (`FRENCHEPHEMERIS_ASTEROID_PACK`, "" disables) and
  cached on-demand packs (`FRENCHEPHEMERIS_ASTEROID_CACHE`, default
  `~/.cache/frenchephemeris/asteroids`) for bodies not already provided.

## On demand

```python
from frenchephemeris import asteroids
p = asteroids.fetch(3753)                  # Cruithne -> <cache>/ast_3753.bin (+ .json provenance)
eng = fe.open_default()                    # loads it
eng.calc_ut(jd, fe.asteroid_id(3753))
asteroids.fetch(99942, "1900-01-01", "2100-01-01")   # custom span
```

Network access (curl when available, else urllib) happens only in this
Python function; the Rust core never touches the network. The SPK is
downloaded into a temporary directory inside the cache and deleted after the
fit. A cached pack is reused unless `force=True`. Throughput is set by
Horizons: an 800-year SPK is ~19 MB of base64; at JPL's bandwidth during
this work (20 KB/s – 1 MB/s) a fetch took 45 s – 9 min.

## Validation

(`harness/check_asteroids.py` → `harness/asteroids_report.json`;
`cargo test -p fe-core --test asteroids`; `tests/test_asteroids.py`.)

* **Fit vs JPL**: max error per body (out-of-sample, at fit time) median
  0.494 km, max 0.500 km (tolerance 0.5 km). Independently, 2,600 JPL SPK
  states recorded at fit time (8 per body, including both range ends) are
  reproduced by the pack through the engine's chaining to ≤ 0.433 km (Rust
  test).
* **vs JPL Horizons**, geocentric astrometric RA/Dec, one epoch per request,
  60 bodies spread over the list (every 5th by number: main belt, NEAs,
  Trojans, centaurs, TNOs), one epoch per era:

  | era | median | p99 | max |
  |---|---:|---:|---:|
  | 1600–1800 | 0.00002″ | 0.0038″ | 0.0092″ (1620 Geographos, an NEA near Earth) |
  | 1800–2100 | 0.00001″ | 0.0001″ | 0.00014″ |
  | 2100–2400 | 0.00002″ | 0.0001″ | 0.00013″ |

  (The engine's J2000 equatorial output now includes the IAU 2006 frame
  bias unless `FLG_ICRS` is passed; the harnesses compare with `FLG_ICRS`.)
* **vs Swiss** (black box): Swiss needs one `se<number>s.se1` file per
  asteroid; the API's ephe directory has only Eris and Eros, so Swiss errors
  ("file not found") for all 325 bundle bodies. No comparison is possible
  with the files the API ships; the API cannot compute any of these bodies
  today.

## Open issues

1. **2101 Adonis is excluded.** Its Horizons SPK (1599–2401) is
   discontinuous at type-21 line boundaries by up to 232 km (2036 Feb, its
   close Earth approach), 55 km (1894) and 6 km (2134), and inside the lines
   around the approach, so no smooth fit reaches 1 km even with 3-minute
   records. The fitter cuts series at line-boundary jumps (> 0.1 km) so such
   steps are reproduced, but the remaining intra-line error was still 83 km;
   `fetch(2101)` raises instead of returning a bad pack.
2. Rebuilding the bundle is bandwidth-bound (15 s – 9 min per body depending
   on JPL's throughput, ~6 h for 326 bodies here); `tools/build_asteroids.py
   fetch` is resumable, then `bundle` rewrites the pack, ledger and fixture.
   The per-body work packs were deleted after bundling (disk).
3. Swiss cannot be used as an oracle for these bodies with the API's ephe
   directory (one `se?????s.se1` file per asteroid; only Eris and Eros are
   present): today's API cannot compute any of them, so the bundle is a pure
   addition.
4. Root scanning uses a conservative 20°/day bound for unknown asteroids; a
   per-body bound from the manifest would speed up searches.
5. Horizons solutions change over time: the ledger pins the solution ids used.
