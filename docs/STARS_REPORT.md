# Fixed stars — contract, catalogue, pipeline, validation

Status: `fe-core::stars` + `Engine::fixstar_ut` / `fixstars_ut` / `fixstar_prepared`,
`Engine::fixstar_mag` / `fixstars_mag` / `fixstar_lookup`, `stars::lookup` / `lookup_compat`,
`Engine::with_star_compat`; Python `fixstar_ut`, `fixstars_ut`, `fixstar_mag`, `fixstars_mag`,
`fixstar_names`, `fixstar_info` (module-level functions taking the `Engine` handle),
`Engine(..., star_compat=path)`; `swe_compat.fixstar2_ut` / `fixstar2_mag` / `fixstar_mag`.

**Default catalogue (public):** 1175 objects, 504 name keys and 6691 designation keys.
- Names come only from public sources: IAU WGSN, SIMBAD `NAME` identifiers, and traditional names from Allen 1899.
- Designations are in standard form (`alf Vir`, `α Vir`, `61 Cyg`, HR, HD, HIP, M/NGC/IC).
- The returned string is `"Name,designation"`, e.g. `"Spica,alf Vir"`.
- The user-facing list is `docs/STAR_NAMES.md`.

**Optional compatibility alias table (private, not in this repository):** loaded at runtime. It holds the classic software's names and `"Name,nomenclature"` strings, and the names that classic software assigns to a different star. With it loaded, 1706 of the 1708 queries the classic API accepts return the identical string and the same object (§2b).

## 1. FreeAstroAPI contract (read from `astro/backend`)

Fixed stars are used in `core/astrology.py` (`AstrologyEngine.calculate_fixed_stars`) and
`core/ephemeris_service.py` (`EphemerisService.calculate_fixed_stars`). They are exposed by
`routers/v1/natal.py` and `routers/v1/ephemeris.py` (`include_fixed_stars`, `fixed_stars`).

- **Request:** `fixed_stars` is a list, or a comma-separated string, of names. They are passed
  **verbatim** to `swe.fixstar2_ut`, so the accepted names are exactly the names Swiss accepts.
  The API has no alias table of its own. `include_fixed_stars=true` with no list uses
  `DEFAULT_ASTROLOGY_FIXED_STARS` (25): Algol, Alcyone, Aldebaran, Rigel, Betelgeuse, Sirius,
  Canopus, Procyon, Castor, Pollux, Regulus, Denebola, Alkaid, Spica, Arcturus, Agena, Toliman,
  Unukalhai, Antares, Vega, Altair, Deneb, Fomalhaut, Scheat, Markab.
- **Normalisation:** names are stripped and de-duplicated case-insensitively. More than
  `MAX_FIXED_STARS = 100` names is an error. An unknown name goes to `fixed_star_errors`
  (`{"name", "error"}`) and does not fail the request.
- **Calls:** `fixstar2_ut(name, jd_ut, SWIEPH | sidereal_flags [| SPEED])`. Sidereal uses the
  same global `set_sid_mode` as the planets. A second call with `SWIEPH | EQUATORIAL` gives
  RA/Dec, always tropical and always with speed off.
- **Output per star**, keyed by `name`:
  - `name` = the part of Swiss's returned `"Name,nomenclature"` before the comma
  - `id` = `fixed_star_<slug(name)>`
  - `requested_name`, `swiss_name` (the full `"Name,nomenclature"` string)
  - `sign`, `sign_id`, `pos`, `abs_pos`, `longitude_deg`, `latitude_deg`
  - `distance_au` (Swiss returns parsec-scale distances in AU)
  - `degree_in_sign`, `right_ascension_deg`, `declination_deg`
  - `speed` (deg/day) and `is_stationary` or `motion_state`
  - `house` (`swe.house_pos`, natal only)
- **Swiss name matching (black-box):**
  - case-insensitive, and spaces are ignored (`"asellusaustralis"`)
  - a leading comma looks the star up by its compact nomenclature; a bare nomenclature is rejected
  - aliases are separate entries returning their own name with the shared nomenclature
  - without `sefstars.txt`, Swiss still knows a small built-in set (Spica, δ Cnc, …). Production has the file in `backend/ephe`.

`fe-core` reproduces this contract:

| Contract item | `fe-core` behaviour |
|---|---|
| Lookup | Case-insensitive, whitespace removed |
| `"Name,xx"` | The name part wins |
| `",designation"` | Accepted. Default: standard designations (`",alf Vir"`). The compatibility table adds the classic compact forms. |
| Also accepted | `"alf Vir"`, `"α Vir"`, `"alf1 Cap"`/`"alf01 Cap"`, `"61 Cyg"`, `"HR 5056"`, `"HD 116658"`, `"HIP 65474"`, `"M 31"` |
| Name semantics | **Default:** IAU WGSN wins, then traditional names, then SIMBAD names. **With the compatibility table:** the classic meaning wins (§2a). |
| Returned name | **Default:** `"Name,designation"` (`"Spica,alf Vir"`; `",109 Vir"` without a public name). The API's `name`/`id` derivation (the part before the comma) keeps working. **With the table:** the classic strings, identical to today's API output. |
| Errors | `InvalidArgument("unknown fixed star …")`, which is `ValueError` in Python. In the batch it is a per-star error. |
| 6-vector | Same as `calc_ut`: lon/lat (or RA/Dec) in degrees, distance in AU (from parallax), daily rates. Flags work as for planets: `EQUATORIAL`, `J2000`, `NONUT`, `TRUEPOS`, `NOABERR`, `NOGDEFL`, `TOPOCTR`, `HELCTR`, `BARYCTR`, `SIDEREAL`. |
| Sidereal | Subtracts `Engine::ayanamsa_deg` exactly as for planets. Tested: Lahiri Spica = tropical − `ayanamsa_deg` (`tests/stars.rs`). |

In default mode, 23 of the 25 API default names resolve to the same object as today. The exceptions are the two α Cen names:
- **Toliman** is α Cen B (IAU);
- **Rigil Kentaurus** is α Cen A (IAU);
- the classic software uses one α Cen point for both, which the compatibility table restores.

With the table, all 25 return today's strings.

## 2. Catalogue: sources, provenance, terms

Generated by `tools/build_stars.py`:
- `fetch` downloads about 260 kB into the gitignored `data/sources/stars/`, using `curl`.
- `build` writes `crates/fe-core/src/stars/catalog.rs` (141 kB, about 30 kB in the binary).

`build` is deterministic for a given cache.

| Content | Source | Provenance / terms |
|---|---|---|
| Names, Bayer/Flamsteed IDs, HIP numbers (every IAU-named star that has a HIP number) | **IAU Catalog of Star Names (IAU-CSN)**, WGSN, machine-readable edition by E. Mamajek, `https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt` (2022-04-04) | The file header states that IAU products are released under **Creative Commons Attribution** (free to use, source must be credited). Cite the IAU WGSN. |
| Positions (ICRS, epoch **J1991.25 TT**), proper motions μα\*, μδ, parallax | **Hipparcos, the New Reduction** (van Leeuwen 2007, A&A 474, 653), **CDS catalogue I/311** (`hip2`), retrieved from a CDS VizieR mirror (`vizier.idia.ac.za`) as FITS binary tables at full precision. The values are identical (all 1150 stars, all five fields) to the copy previously taken from the ESA Gaia archive table `public.hipparcos_newreduction`, which is no longer used. | ESA scientific data, free to use with acknowledgement of ESA/Hipparcos and citation of van Leeuwen 2007. See open issue 1. |
| Radial velocities | **SIMBAD** (CDS) `basic.rvz_radvel`, via the CfA mirror `simbad.harvard.edu` TAP. The bibcode of every value is embedded (`StarRecord::rv_ref`). Where SIMBAD's adopted value is Gaia-derived, the best non-Gaia measurement from SIMBAD's `mesVelocities` table is used instead (§2c). The sources are mostly Gontcharov 2006 (Pulkovo, CDS III/252), Nordström 2004, Famaey 2005, Kharchenko 2007, Wilson 1953 / Evans 1967, and spectroscopic surveys. | CDS asks for the acknowledgement "This research has made use of the SIMBAD database, operated at CDS, Strasbourg, France". The values are factual measurements credited to the cited papers. |
| α Cen AB barycentre record (`SYSTEMS`; `hip = 0`, designation `alf Cen AB`) | **Position:** the mass-weighted mean of the Hipparcos 2007 A and B places at J1991.25, with masses 1.1055 / 0.9373 M☉ (Kervella et al. 2016, A&A 594, A107). **Motion:** SIMBAD `alf Cen` system proper motion (−3608, +686 mas/yr) and RV −22.3 km/s (1979IAUS...30...57E). **Parallax:** Akeson et al. 2021 (AJ 162, 14), 750.81 mas. | All values are published; see §2a. |
| α Cen A and B motion override | Both components get the same system motion | Their Hipparcos proper motions include orbital motion (≈ 450″ after 2000 yr). With the override, each keeps its J1991.25 offset from the barycentre, so the error is bounded by the orbital amplitude (≲ 10″). |
| Additional objects (`tools/star_extra_objects.tsv`: 1088 HIP numbers / SIMBAD identifiers, **no names**) | Objects identified, during the compatibility work, behind names the classic API accepts: designation → SIMBAD → HIP. They are real Hipparcos/SIMBAD objects. Their default names come only from the public sources in this table; objects without a public name are reachable by designation. | SIMBAD (CDS), ESA |
| Names and designations of every object | SIMBAD (CDS) identifiers: `NAME …` names, Bayer/Flamsteed (`* alf Vir`), variable-star names, HR, HD, HIP, M/NGC/IC. Greek letters are also accepted in Unicode, and superscripts without leading zeros. | SIMBAD acknowledgement |
| Optional compatibility alias table (private: `private/stars/star_compat_aliases.tsv`; env `FRENCHEPHEMERIS_STAR_COMPAT`) | Built by the public tool `tools/build_star_compat.py` from private probe outputs: the oracle's accepted queries and returned strings (`tools/probe_swiss_stars.py`, public API only, the oracle's data file never opened), mapped to catalogue records by `tools/identify_swiss_stars.py`. | Not distributed |
| Objects without a Hipparcos 2007 solution (clusters, galaxies, η Car, γ¹ Leo, α¹ Her, 61 Cyg, PSR B1257+12, the Trapezium, …) | SIMBAD ICRS J2000 place, proper motion, parallax and RV, propagated to J1991.25 with ERFA `pmsafe`. The bibcode of every field is stored in the record's `rv_ref`. Open clusters whose SIMBAD values are Gaia-derived take place, motion and distance from Dias et al. 2002 (CDS B/ocl v3.5). For multiple-star entries that SIMBAD lists without motion (61 Cyg), we keep the system place and use component A's motion. ζ¹ Aqr is represented by the Hipparcos solution of the ζ Aqr pair (HIP 110960). | SIMBAD (CDS), CDS VizieR, and the cited papers |
| Galactic Center (lookup `"Galactic Center"`, `"Sgr A*"`, SIMBAD names; record id `SIMBAD NAME Sgr A*`) | **Place:** SIMBAD `NAME Sgr A*` (2011AJ....142...35P). **Proper motion:** Reid & Brunthaler 2020 (ApJ 892, 39), μ_l = −6.411, μ_b = −0.219 mas/yr, converted to ICRS. **Distance:** R₀ = 8.178 kpc (GRAVITY 2019, A&A 625, L10). | Published |
| Great Attractor, Virgo Cluster | These are conventional points, not catalogue objects. **GA:** galactic l = 325.3°, b = −7.2° (cf. Kraan-Korteweg et al. 1996), distance 4350 km/s (Lynden-Bell et al. 1988) / H₀ = 70. **Virgo:** 12h24m, +13°00′ (FK4 B1950), distance 16.5 Mpc (Mei et al. 2007). Both round-number definitions agree with the classic software's place to < 0.001″. | Conventional definitions |
| V magnitudes (`fixstar_mag`) | Johnson V from the SIMBAD `flux` table, with its bibcode, for 1118 objects. For 45 HIP stars without a SIMBAD V, we use the Hipparcos Catalogue (ESA 1997) field H5 V from CDS catalogue I/239 (VizieR). α Cen AB is −0.27, combined from A and B. 13 objects with no published V return 999.99, the classic sentinel. | SIMBAD (CDS); ESA |
| Traditional names (`EXTRA_ALIASES`, 23) | R. H. Allen, *Star Names: Their Lore and Meaning* (1899, public domain), and common astrological usage: Agena = β Cen, Bungula = α Cen (AB barycentre), Pushya = δ Cnc, Benetnasch, Alnath, Wega, Gemma, Sirrah, Graffias, Cor Scorpii, Galactic Center = Sgr A* … IAU names take precedence, so "Rigil Kentaurus" / "Toliman" in this list only matter with the compatibility table. | Public domain / usage |

### 2c. No Gaia-derived values

Gaia data are released under CC BY-SA 3.0 IGO (share-alike), so the catalogue contains no Gaia-derived value.
- **Rule:** `tools/build_stars.py` treats a value as Gaia-derived when its bibcode is:
  - a Gaia release (VizieR I/337, I/345–I/358);
  - in an explicit list of Gaia-based publications (`GAIA_BIBCODES`: cluster and globular-cluster catalogues built on Gaia, astrometric radial velocities from Gaia, the Gaia RV-standard catalogue, …, plus GALAH DR3 and one cluster-distance paper excluded out of caution);
  - or a publication whose title names Gaia or TGAS. Titles are taken from SIMBAD's `ref` table.
- **Replacement:** each such value is replaced by the best non-Gaia published value, or dropped.
- **Enforcement:** the generator asserts that no record cites a Gaia-derived source, and `tests/stars.rs::no_gaia_derived_values_in_catalogue` checks the shipped table.
- **Log:** every replaced value is listed in `tools/star_gaia_replacements.tsv`.

| Field | Replaced | Dropped | Replacement sources |
|---|---:|---:|---|
| Radial velocity | 426 | 1 (Trapezium: no non-Gaia value; RV = 0) | SIMBAD `mesVelocities`, best quality grade first, then the compilations Gontcharov 2006 (231), Nordström 2004 (60), Famaey 2005 (46), Evans 1967 (15), Kharchenko 2007, and others, then the newest. Open clusters: Dias 2002 / Mermilliod 2008. |
| Position | 5 | 0 | Open clusters NGC 869, NGC 1039 (M34), NGC 2632 (Praesepe), NGC 6405 (M6), NGC 6475 (M7): Dias et al. 2002 (B/ocl v3.5) |
| Proper motion | 8 | 0 | the same five clusters (Dias 2002 / Dias 2014, UCAC4-based); M22 (Cudworth 2013); M87 (VLBA, Truebenbach & Darling 2017); 61 Cyg A (Hipparcos 1997) |
| Parallax | 9 | 0 | five clusters (Dias 2002 distances); M22 (Baumgardt & Hilker 2018 distance); NGC 6530 and the Trapezium (published distances); 61 Cyg A (Hipparcos 1997) |
| V magnitude | 0 | 0 | none were Gaia-derived |

ζ¹ Aqr's SIMBAD place and motion were Gaia DR2 values, so its record is now the Hipparcos 2007 solution of the ζ Aqr pair (HIP 110960).

**Hipparcos table re-sourced.** The Hipparcos 2007 astrometry had been read from the ESA Gaia archive. To avoid any doubt about that archive's terms, it is now read from CDS catalogue I/311 through a VizieR mirror, and the Hipparcos 1997 V magnitudes from I/239. Every value is identical:
- I/311: 1150 stars × RA, Dec, parallax and both proper-motion components;
- I/239: 1152 V magnitudes.

The rebuilt catalogue before the Gaia replacements was byte-for-byte the same as the previous one.

**Accuracy impact** (positions vs the oracle, apparent ecliptic of date; before → after):

| Set | Year 0 median / p95 / max | Year 2000 median / p95 / max | Year 2400 median / p95 / max |
|---|---|---|---|
| Public: 119 IAU names | 0.0023 / 1.85 / 17.1″ → 0.0018 / 1.85 / 17.1″ | 0.0000 / 0.0034 / 2.35″ → unchanged | 0.0004 / 0.24 / 4.0″ → 0.0004 / 0.34 / 4.0″ |
| Private: all 1188 oracle objects (compat table) | 0.0013 / 0.37 / 76.5″ → 0.0009 / 0.45 / 76.5″ | 0.0000 / 0.0005 / 1.89″ → 0.0000 / 0.0004 / 0.48″ | 0.0003 / 0.057 / 15.9″ → 0.0002 / 0.065 / 15.3″ |

How to read this:
- Medians improve slightly.
- p95 away from 2000 grows a little: the replacement RVs are older ground-based values, and the oracle's own catalogue partly uses Gaia RVs.
- The API's 25 default stars are unchanged within 0.001″ (max 1.8″ at year 0, 0.09″ after 1700).
- The year-0 maximum (M31, 76″) is the oracle's own drift and is unaffected.
- Cluster centres now come from Dias 2002. The Praesepe offset from the oracle drops from 1308″ to 1168″; it remains an exception.
- ζ¹ Aqr becomes a 3.1″ exception at 2000: the pair's Hipparcos solution vs the component.

### 2a. Name semantics: default (IAU) vs compatibility table

**Default.** IAU WGSN names win. Six names (Menkar, Deneb Kaitos, Kurhah, Sadalbari, Ukdah, Algedi) denote a *different star* in classic astrology software; in default mode they are the IAU stars. Rigil Kentaurus and Toliman are α Cen A and B, and Bungula is the AB barycentre.

**Compatibility table.** The six names and the α Cen names resolve to the star classic software uses, and the returned strings are the classic ones. The IAU stars stay reachable by designation or HIP. This mapping lives only in the private table.

**α Cen.** Swiss's Rigil Kentaurus / Toliman place (black-box, barycentric astrometric ICRS) has:
- an implied proper motion of exactly (−3607.8, +686.0) mas/yr between 1990 and 2010, which is the SIMBAD system value;
- a distance that implies a parallax of 742.0 mas;
- a position 2.4″ from our mass-weighted A/B barycentre at J1991.25, versus 8.2″ from A and about 4″ from the photocentre.

That is consistent with a published system (barycentre) position, so the catalogue has an AB barycentre record (`SYS alf Cen AB`).
- **Default:** it is reached as "Bungula" or `alf Cen AB`.
- **With the compatibility table:** the classic α Cen names resolve to it.
- **Components:** reachable as `alf Cen A` / `HIP 71683` and `alf Cen B` / `HIP 71681`.

Sensitivity of our barycentre vs Swiss (astrometric, arcsec):

| Model | Year 0 | 1000 | 2000 | 2400 |
|---|---:|---:|---:|---:|
| Before: A place + system pm, RV −18.6 (mass-weighted A/B) | 39.7 | 13.5 | 8.2 | 8.9 |
| Barycentre + mass-weighted HIP2 pm | 159 | 75 | 2.7 | 30 |
| Barycentre + SIMBAD system pm, RV −18.6 | 35.7 | 7.9 | 2.4 | 1.3 |
| **Barycentre + SIMBAD system pm/RV −22.3, ϖ 750.81 (adopted)** | **3.5** | **2.2** | **2.4** | **3.0** |

The mass-weighted A/B RV is −18.6 km/s, but the literature systemic RV is about −22.4 km/s; the A and B RVs were measured at different orbital phases. The −22.3 km/s system RV removes the year-0 drift. The remaining 2.4″ is the uncertainty of the barycentre place: the masses and the Hipparcos solutions for a 19″ pair.

### 2b. Compatibility table: black-box discovery (private outputs)

`tools/probe_swiss_stars.py` asked the licensed oracle's public API (`swe.fixstar2_ut`) about candidates from public sources:

- **Names (9073):**
  - IAU-CSN;
  - Wikipedia star-name pages (link texts and runs of capitalised words);
  - the API's 25 default stars and about 280 traditional astrological names;
  - a deep-sky / reference-point list.
- **Designations (712 145):** in the classic compact form:
  - Greek letters with superscripts;
  - Latin letters;
  - Flamsteed numbers;
  - variable-star names;
  - HR, HD, M, NGC and IC numbers.
- Queries starting with a digit are star-number lookups in the oracle. They were excluded and are not supported, and star numbers were not enumerated, so the oracle's data file is never reproduced.

The outputs (accepted queries, returned strings, identifications) are written only to the gitignored `private/stars/`.

| Count | What |
|---:|---|
| 9 073 | name candidates probed |
| 712 145 | designation candidates probed |
| 1 708 | accepted queries |
| 615 | — of which names |
| 1 093 | — of which designations |
| 1 207 | distinct returned strings |
| 583 | — of which designation-only |
| 1 205 | identified in SIMBAD / Hipparcos |
| 2 | entries with no published source (a solar-apex point and a galactic-equator point whose definitions could not be matched). These are not reproduced. |
| 1 706 / 1 708 | queries that resolve identically with the compatibility table (asserted when the private files are present) |

In default mode:
- 440 of the 1175 objects have a public name.
- The others are reachable by designation, HR/HD/HIP or catalogue number, and return `",designation"`.
- In 14 names the classic spelling differs from the IAU spelling; default mode returns the IAU spelling.

Why Hipparcos 2007 and not Gaia DR3:
- Most astrologically used stars are brighter than G ≈ 3 and are saturated or missing in Gaia DR3.
- Gaia DR3 astrometry is unreliable for stars brighter than G ≈ 6 (for the ones it has at all).
- Hipparcos 2007 is the reference for this magnitude range, and one catalogue keeps the build simple and consistent.
- Gaia DR3 could replace the fainter WGSN stars (V > 6, mostly exoplanet hosts) later. That changes nothing for the API's stars.

Clean-room:
- `sefstars.txt` was never opened, copied or used.
- The oracle was used only through its public API: to compare positions, and to build the private compatibility table.
- No oracle-derived string is in the public catalogue, tests or documentation.

## 3. Pipeline (`crates/fe-core/src/stars.rs`)

1. **Space motion.** We implement SOFA `starpv` → `starpm` → `pvstar` → `pmsafe` ourselves, from the published algorithms. This includes the relativistic Doppler factor and the catalogue light-time (Roemer) term. It runs from the catalogue epoch J1991.25 to TDB, and gives the barycentric position in AU.
2. **Parallax.** We take the direction from the observer's barycentric position: `ephem::ObsState` for the geocentre, topocentre, Sun or SSB.
3. **Solar deflection.** SOFA `ldsun`, which is `ld` with q = p and the `dlim` limiter. It applies only to geocentric or topocentric observers, and can be switched off with `NOGDEFL`.
4. **Aberration.** SOFA `ab` (relativistic), with the same `ObsState` velocity and Sun distance as the planets. Switched off by `NOABERR` / `TRUEPOS`.
5. **Frame.** `ephem::to_frame` with the epoch's `FrameSet` (`ecl_true`, `bpn`, `bp`, `ecl_mean`, J2000, ICRF). The flag mapping is the same as `CalcOptions::frame`.
6. **Sidereal.** `lon -= Engine::ayanamsa_deg(mode, ep, !mean_equinox)`, applied only to ecliptic output.
7. **Speeds.** Central differences at ±1e-3 d on `Prepared::shifted` epochs, with observer states recomputed there, exactly as in `calc_prepared`.

Batching:
- `fixstars_ut` prepares the instant, the frames and the three observer states once per call.
- The catalogue-epoch state `starpv` is evaluated once per star and reused for all three epochs.

## 4. Validation (`crates/fe-core/tests/stars.rs`)

Fixtures live in `tests/fixtures/stars/`:
- `gen_erfa.py` → `erfa.json`
- `gen_swiss.py` → `swiss.json`
- `diag_swiss.py`, a diagnostic that prints and writes nothing

Tests that need the pack skip when it is missing.

### (a) Against ERFA (pyerfa 2.0.1.5)

Fixture: 146 stars × 9 dates (0–2400 CE), 1314 rows.

| Quantity | Max difference |
|---|---|
| `pmsafe` position | 9e-5 µas |
| proper motion | 4e-9 µas/yr |
| parallax | 1e-10 µas |
| RV | 2e-13 km/s |
| Astrometric direction (rigorous parallax from the observer: `starpv` position − `epv00` Earth) | 1.8e-4 µas |
| Apparent direction (`ldsun` + `ab`) | 1.8e-4 µas |
| Distance (relative) | 5e-16 |

The test asserts < 1 µas.

### (b) Against Swiss Ephemeris 2.10.03

Setup:
- `swe.fixstar2_ut(name, jd_ut, SWIEPH|SPEED)`, which are the API's flags.
- Production star file directory, plus the planetary files.
- 15 dates from year 0 to 2400, on 1 January.

There are two tests:
- **Public** (`apparent_places_vs_swiss_oracle`): 119 IAU names (the API defaults plus every 3rd IAU name the oracle accepts). It excludes the six names that denote another star and the two α Cen names. The test asserts the same object for each name and prints the per-era statistics; the current run gives an overall median of 0.0003″ and p95 of 0.24″.
- **Private** (`every_swiss_accepted_name_resolves_identically_with_compat_table`): runs when the compatibility table and its oracle fixture are present. It asserts the identical string for every accepted query.

The table below is the earlier detailed set of 328 names × 15 dates.

**Apparent ecliptic-of-date separation (arcsec), 328 names × 15 dates = 4920 rows, API frame model (`Engine::new`):**

| Year | Median | p95 | Max |
|---:|---:|---:|---:|
| 0 | 0.182 | 1.56 | 26.5 |
| 300 | 0.088 | 1.26 | 22.1 |
| 600 | 0.033 | 0.99 | 17.8 |
| 900 | 0.0056 | 0.74 | 13.7 |
| 1200 | 0.0048 | 0.52 | 9.7 |
| 1500 | 0.0048 | 0.32 | 5.9 |
| 1700 | 0.0026 | 0.32 | 4.0 |
| 1800 | 0.0014 | 0.21 | 3.1 |
| 1900 | 0.0005 | 0.10 | 2.3 |
| 1950 | 0.0002 | 0.047 | 2.3 |
| 2000 | 0.0000 | 0.0093 | 2.4 |
| 2025 | 0.0001 | 0.037 | 2.4 |
| 2100 | 0.0001 | 0.12 | 2.5 |
| 2200 | 0.0008 | 0.23 | 2.6 |
| 2400 | 0.0049 | 0.32 | 4.3 |
| **all** | **0.0028** | **0.25** | **26.5** |

The 2000–2200 maximum is α Cen (2.4–2.6″). The early maxima are Barnard's Star, Proxima, Algieba and Acamar.

The 25 API default stars:

| Years | Median | p95 | Max |
|---|---:|---:|---:|
| 0–1500 | 0.036″ | 1.46″ | 3.6″ (Toliman at year 0; was 39.6″) |
| 1700–2400 | 0.0008″ | 0.044″ | 3.0″ (Toliman; was 8.9″) |

Maximum over 0–2400 for the API-continuity names:

| Name | Max |
|---|---:|
| Menkar | 0.07″ |
| Deneb Kaitos | 0.10″ |
| Sadalbari | 0.14″ |
| Kurhah | 0.18″ |
| Ukdah | 0.18″ |
| Algedi | 0.21″ |
| Rigil Kentaurus / Toliman / Bungula | 3.6″ |

With `SWISS_COMPAT` frames, the six names drop to ≤ 0.02″ and α Cen stays at 3.5″.

The API rounds `longitude_deg` to 1e-6° (0.0036″) and `pos` to 0.001° (3.6″). So the default stars agree with Swiss at the displayed `pos` precision everywhere, α Cen included.

**Causes, separated with `FrameModel::SWISS_COMPAT` and `diag_swiss.py`:**

1. **Frame model (year-0 median).** The median at year 0 drops from 0.18″ to 0.009″ with `SWISS_COMPAT` frames. Our ecliptic of date uses the Vondrák pole angle for the obliquity, while Swiss uses its own ε_A, which is 0.2″ away at year 0 (see TIME_FRAMES_REPORT §7). This is the same difference the planets show. It is a modelling choice, not an error.
2. **Catalogue proper motions and RVs (the p95/max growth away from 2000).** We compared barycentric astrometric ICRS places, which involve no frames, aberration or deflection:
   - median 0.0000″ / p95 0.009″ at 2000;
   - median 0.0017″ / p95 1.42″ at year 0.

   This growth is linear in |t − 2000|, which points to the input proper motions and RVs, not the propagation algorithm. Swiss also propagates rigorously with RV: Barnard's Star has a 2239″ RV/perspective term but differs by only 26″. The stars that differ most have:
   - high proper motion or a different RV source: Barnard's Star 26″, Proxima 22″;
   - a close binary, where the Hipparcos solution differs from Swiss's source: Algieba (γ¹ Leo) 20″, Acamar (θ¹ Eri) 17″, Almach 3.4″.
3. **α Cen** (Rigil Kentaurus / Toliman / Bungula): now the AB barycentre (§2a).
   - 2.4″ at 2000 and ≤ 3.6″ over 0–2400; before this change it was 8.2″ and 39.6″.
   - The residual is the barycentre-place uncertainty. If A's place were used with the instantaneous Hipparcos proper motion, it would be 457″ at year 0.
4. **Speeds.**
   - Our rates match Swiss's own positions differentiated over ±0.05 d to a median of 7e-7″/day and a p95 of 3e-6″/day.
   - Swiss's *reported* star speeds differ from its own positions by a median of 0.017″/day and up to 0.046″/day, so they are approximate.
   - As a result, our `lon_speed` differs from Swiss-reported speeds by up to 0.92″/day for Fafnir, which is near the ecliptic pole where 1/cos β amplifies everything. That is 2.5e-4°/day, far below `is_stationary`'s 0.001°/day.
5. **Distance.** The median relative difference is 3e-6. Stars where Swiss's parallax source differs reach up to 38% (Acamar) and 11% (Almach).

### (b2) All Swiss-accepted objects (`swiss_names.json`, 1 205 entries, years 0 / 2000 / 2400)

Apparent ecliptic of date, our place vs `swe.fixstar2_ut` (arcsec), excluding the exceptions below:

| Year | Median | p95 | Max |
|---:|---:|---:|---:|
| 0 | 0.18 | 0.42 | 37.8 (61 Cyg) |
| 2000 | 0.0000 | 0.0004 | 0.48 |
| 2400 | 0.0048 | 0.047 | 8.2 (61 Cyg) |

Away from 2000 the largest differences are catalogue proper motions and RVs, as in §4b(2): 61 Cyg 38″, ω Phe 33″, Barnard's star 26″, υ Phe 25″, Proxima 22″, Acamar 17″ (year 0).

**Every object more than 1″ from the oracle at 2000 was investigated.** The test asserts these as exceptions, by record id:

| Object | Δ at 2000 | Cause |
|---|---:|---|
| η Ori (one oracle entry) | 4.5° | The oracle's RA for that entry is 18 min off η Ori; its magnitude is η Ori's. |
| Praesepe (NGC 2632) | 1308″ | Cluster centre: SIMBAD (Tarricq 2021) vs the oracle's older centre |
| M7 | 178″ | Cluster centre |
| M34 | 145″ | Cluster centre |
| M6 | 60″ | Cluster centre |
| NGC 869 | 27″ | Cluster centre |
| M42 | 14.5″ | Nebula centre |
| Trapezium | 8.9″ | SIMBAD centre vs the oracle's θ¹ Ori point |
| α Cen AB | 2.4″ | barycentre vs the oracle's system place (§2a) |
| ζ Vir | 2.3″ | close visual binary: component / photocentre |
| ζ¹ Aqr | 1.9″ | same |
| ζ Cnc | 1.4″ | same |
| γ¹ Leo | 1.3″ | same |
| ξ Sco | 1.0″ | same |

Defined points agree to < 0.001″: Great Attractor, Virgo Cluster, Galactic Center = SIMBAD Sgr A*. M31 and NGC 6530 agree to < 0.001″. ε Mon agrees to 0.00006″ once it is mapped to component A.

### (b3) Magnitudes (`fixstar_mag` vs `swe.fixstar2_mag`)

1678 comparable entries:

| Statistic | Value |
|---|---|
| Median \|Δ\| | 0.000 mag |
| p95 \|Δ\| | 0.020 mag |
| More than 0.1 mag | 25 |

The largest differences:

| Object | Ours (V) | Oracle | Note |
|---|---:|---:|---|
| NGC 4194 | 13.3 | 0.0 | the oracle has no magnitude |
| ξ Cep | 4.29 | 6.5 | the oracle's value looks like the B component |
| T CrB | 10.25 | 10.8 | variable star |
| ρ Her | 4.17 | 4.56 | |
| 61 Cyg | 5.21 | 4.84 | |
| α Her | 2.78 | 3.06 | |
| η Car | 6.48 | 6.21 | variable star |
| γ¹ Leo | 2.23 | 1.98 | |
| α Cen | −0.27 | −0.10 | |

Objects without a V magnitude return 999.99, as Swiss does: Galactic Center, Great Attractor, Virgo Cluster and 10 clusters/nebulae/pulsars.

### (c) JPL Horizons

Not feasible. Horizons does not offer arbitrary catalogue stars as targets, so there is no independent JPL truth for stellar places. ERFA serves as the independent reference instead.

### Timings (Apple Silicon)

| Build | `fixstar_ut` with speed | without speed | 25-star batch (with speed) |
|---|---:|---:|---:|
| debug (opt-level 1) | 5.8 µs | 3.2 µs | 19 µs (0.76 µs/star) |
| release | 5.5 µs | 2.9 µs | 25 µs (1.0 µs/star) |

The release timings are noisy because they were measured with other tests running. A single call is dominated by preparing the instant: frames plus three pack evaluations. Swiss `fixstar2_ut` (with `sefstars.txt`, Python) takes 8.7 µs for one star and 33.6 µs for the 25 defaults.

## 5. Open issues

1. **Licence confirmation for redistribution in a commercial binary.**
   - IAU-CSN (CC BY) needs attribution in the product docs.
   - Hipparcos values come from CDS catalogues I/311 and I/239, and the other values from SIMBAD and VizieR (CDS). These are free-use scientific data with acknowledgement. We embed 1175 rows of factual values, but explicit confirmation of commercial-redistribution terms from ESA and CDS is recommended. No Gaia-derived value is used (§2c).
   - Add the acknowledgement lines to the product's notices.
2. **Default vs classic naming.**
   - Without the private compatibility table, returned strings use IAU names and standard designations.
   - Six names, plus Toliman and Rigil Kentaurus, denote the IAU star instead of the classic one.
   - Deployments that need today's output (FreeAstroAPI) must ship the private table and set `FRENCHEPHEMERIS_STAR_COMPAT`; `open_default()` also picks up `private/stars/star_compat_aliases.tsv`.
   - Earlier commits on the `ws/stars*` branches contain oracle-derived strings in `catalog.rs` / `tools/*.tsv`. History rewriting is needed if those branches are ever published.
3. **Names Swiss accepts that we do not.** We cannot enumerate them without reading Swiss's data file, which is not allowed. Recommendations:
   - log `fixed_star_errors` / requested names in production;
   - extend `EXTRA_ALIASES` from real traffic;
   - Swiss's `%` wildcard search and star numbers are not supported.
4. **Binary orbital motion is not modelled.** Hipparcos proper motions of astrometric binaries (Sirius, Procyon, Algieba, …) include orbital motion. This is mitigated only for α Cen (barycentre record; A and B carry the system motion). Over millennia it amounts to tens of arcseconds for a few stars.
5. **Alula Australis (ξ UMa, HIP 55203)** has no Hipparcos 2007 solution. It is now built from SIMBAD's place and motion.
6. **Sidereal output** works through `Engine::ayanamsa_deg`. The ayanamsa module has its own δ Cnc model (`ayanamsa::star`). If it should use this catalogue's δ Cnc instead, wire it to `stars::lookup("Pushya")` + `fixstar_prepared`.
7. **Clusters and galaxies.**
   - The catalogue uses SIMBAD centres. These differ from the classic software's by up to 22′ (Praesepe), because centre definitions vary between catalogues.
   - Distances for objects without a parallax are not meaningful: parallax is floored.
   - The oracle's M31 drifts by 38 mas/yr; ours is fixed (no published proper motion). That gives 76″ at year 0.
8. **Two reference points** that the classic API accepts are not reproduced, because no published definition matched.
9. **Swiss star-number queries** (a leading digit) are not supported. FreeAstroAPI passes user names through, so a request such as `"4 Sagittarii"` currently returns an unrelated star from Swiss; with fe-core it returns an error.
10. **RV = 0 stars:** SOFA's `starpv` is ill-conditioned when RV is exactly 0. Our results and ERFA's differ by up to 14 µas over 2000 years (rounding only).
