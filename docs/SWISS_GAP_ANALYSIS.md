# What FrenchEphemeris still lacks to replace Swiss Ephemeris in full

Source: the current public Swiss Ephemeris documentation
(https://www.astro.com/swisseph/swephprg.htm, https://www.astro.com/swisseph/swisseph.htm),
read 2026-10-29. Documentation only, no Swiss source code.

FrenchEphemeris already covers everything FreeAstroAPI calls. This list is what a general
astrology program could call that we do not offer yet. Size: S ≈ days, M ≈ 1–2 weeks,
L ≈ several weeks.

## Already covered

| Area | Swiss API | FrenchEphemeris |
|---|---|---|
| Positions | `calc`, `calc_ut` | yes. Sun–Pluto, mean/true node, mean/osculating Lilith, Earth, Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Eros |
| Flags | TRUEPOS, NONUT, EQUATORIAL, J2000, NOABERR, NOGDEFL, HELCTR, BARYCTR, TOPOCTR, SIDEREAL, SPEED, XYZ, RADIANS, ICRS, ASTROMETRIC, SPEED3 | yes (B1–B3 done 2026-10-29) |
| Houses | `houses`, `houses_ex`, `houses_ex2`, `houses_armc`, `houses_armc_ex2`, `house_pos`, `gauquelin_sector`, `house_name` | yes, all documented systems: P K O R C T E W A V D N B X F M H U Y L Q S I i G J |
| Sidereal | `set_sid_mode`, `get_ayanamsa(_ex)(_ut)`, `get_ayanamsa_name`, `SE_SIDM_USER`, SIDBIT options | 44 of 47 modes (all but J.N. Bhasin, Sheoran, Valens), all six SIDBITs, per call; see GEOMETRY_REPORT §2 |
| Stars | `fixstar(2)(_ut)` | yes, 410+ stars, Swiss-compatible names in progress |
| Events | `rise_trans`, `rise_trans_true_hor`, `pheno(_ut)`, `solcross(_ut)`, `mooncross(_ut)`, `time_equ`, `lmt_to_lat`, `lat_to_lmt` | yes |
| Eclipses | `lun_eclipse_when/how`, `sol_eclipse_when_glob/when_loc/where/how`, `lun_occult_when_glob/when_loc/where` | yes, all fields (saros, local noon, shadow limits added 2026-10-29) |
| Coordinates | `cotrans`, `cotrans_sp`, `azalt`, `refrac`, `sidtime` | yes |
| Time | `julday`, `revjul`, `deltat`, `deltat_ex`, `utc_to_jd`, `jdet_to_utc`, `jdut1_to_utc`, `utc_time_zone`, `date_conversion`, `day_of_week`, `sidtime0`, `set_delta_t_userdef`, `get/set_tid_acc` | yes (F1–F3 done 2026-10-29) |
| Quick wins (2026-10-29, branch ws/quickwins) | `ECL_NUT`, `calc_pctr`, `lun_eclipse_when_loc`, `helio_cross(_ut)`, `mooncross_node(_ut)`, `split_deg`, `cs2*str`, `csnorm`, `csroundsec`, `d2l`, `degnorm`, `radnorm`, `difdeg*`, `difcs*`, `difrad2n`, `get_planet_name`, `version`, `get_library_path`, `get_current_file_data`, `get_ayanamsa_name` | yes; deltas vs Swiss in TIME_FRAMES_REPORT.md §10, EVENTS_REPORT.md and ECLIPSES_REPORT.md ("Quick wins" sections) |

## Missing

### A. Bodies and points

| # | Item | Swiss | Size | Notes |
|---|---|---|---|---|
| A1 | Any numbered asteroid (`AST_OFFSET + n`) | separate files for ~all MPC-numbered objects | L | We have 8. Options: a bundle of the ~300 astrologically used asteroids fitted from JPL (same pipeline as the minor pack), plus on-demand fetch-and-fit from JPL Horizons for others. |
| A2 | Interpolated ("natural") apogee and perigee (`INTP_APOG`, `INTP_PERG`: natural Lilith, Priapus) | yes | M | Definition documented in swisseph.htm §2.2.4; must be implemented from that description. |
| A3 | Hypothetical / fictitious bodies (`FICT_OFFSET`): Uranian planets Cupido–Poseidon, Transpluto/Isis, Harrington, Nibiru, Vulcan, White Moon/Selena, Proserpina, Waldemath, Leverrier/Adams Neptune, Lowell/Pickering Pluto | orbital elements in a data file | M | Elements must come from the original publications (Witte/Sieggrün, Landscheidt, …), not Swiss's file. |
| A4 | Planetary moons and planet centres (`PLMOON_OFFSET`, `SEFLG_CENTER_BODY`) | yes (2.10) | M–L | Needs JPL satellite kernels (large). Planet centre vs barycentre matters below 0.001° except for the moons themselves. |
| A5 | Comets / interstellar objects | via files | M | Same pipeline as A1, from JPL Horizons. |
| A6 | `ECL_NUT` pseudo-body (obliquity and nutation) | yes | S | **Done** (2026-10-29): Δψ, Δε identical; obliquity ≤ 1 mas in 1800–2200, ≤ 0.2″ at year 0 (engine's Vondrák obliquity). |
| A7 | Planetary nodes and apsides, `nod_aps(_ut)`: mean and osculating nodes, perihelion/aphelion, focal points, all planets | yes | M | Osculating elements from our JPL states; mean ones from published planetary theories. |
| A8 | `get_orbital_elements`, `orbit_max_min_true_distance` | yes | S–M | From JPL states. |
| A9 | Planet-centric positions `calc_pctr` | yes | S | **Done** (2026-10-29): `Origin::Body`; positions ≤ 0.5″ vs Swiss. |

### B. Calculation flags

| # | Flag | Size | Notes |
|---|---|---|---|
| B1 | `SEFLG_XYZ` (cartesian), `SEFLG_RADIANS` | S | **Done** (2026-10-29). |
| B2 | `SEFLG_ICRS`, `SEFLG_ASTROMETRIC` | S | **Done** (2026-10-29). J2000 equatorial output now applies the IAU 2006 frame bias unless ICRS is requested (as Swiss; was ICRS before, 0.024″ off). |
| B3 | `SEFLG_SPEED3` | S | **Done** (2026-10-29): accepted as SPEED (ours is always the derivative). |
| B4 | `SEFLG_DPSIDEPS_1980`, `JPLHOR`, `JPLHOR_APPROX` (reproduce Horizons' IAU 1976/1980 frame) | M | **Done** (2026-10-29, ws/occultations): IAU 1976/1980 + pinned IERS EOP 14 C04 pole offsets (or modelled ones); vs Horizons ≤ 12 mas frame residual 1700–2400; before ~1600 Horizons departs from published IAU 1976/1980 (4″ at year 0, not reproduced). TIME_FRAMES_REPORT §11. |
| B5 | `SEFLG_CENTER_BODY` | see A4 | |

### C. Houses

| # | Item | Size | Notes |
|---|---|---|---|
| C1 | **Done** (ws/houses-all): all 18 systems (incl. 'J' Savard-A from J. Savard's published description, ws/gaps-houses) with cusps, ascmc, `house_pos` and the oracle's edge rules; deltas in GEOMETRY_REPORT §1. Open: 76 high-latitude 'i' (Makransky) rows, his calculation recipe is only in a copyrighted book | M | |
| C2 | **Done**: `house_cusps_speeds`, `Engine::house_cusps_ut`; oracle speed defects documented | S | |
| C3 | **Done**: `Engine::gauquelin_sector`, methods 0–5 (planets; fixed stars not yet) | S–M | |
| C4 | **Done** | S | |

### D. Sidereal zodiac

Done (branch ws/ayanamsas-all): D1 44 of the 47 Swiss modes from their published
definitions, D2 `SE_SIDM_USER` and SIDBIT_ECL_T0 / SSY_PLANE / USER_UT /
ECL_DATE / NO_PREC_OFFSET / PREC_ORIG as per-call options, D3
`get_ayanamsa_name`. Per-mode deltas: GEOMETRY_REPORT.md §2. Remaining:

| # | Item | Size | Notes |
|---|---|---|---|
| D1a | J.N. Bhasin (8), "Vedic"/Sheoran (39), Vettius Valens (42) | S | No published numeric definition was accessible (the documentation gives none; the oracle shows Sheoran is δ Cnc-anchored). Need the primary sources. |
| D1b | Residuals above 0.001″: Djwhal Khul 4.9″ (oracle anchored at J1900), Eta Piscium 0.93″, Lahiri ICRC 0.066″ (Woolard nutation), Gal. Equator IAU1958 0.18″, Galactic-Centre modes ≤ 0.10″ in 1800–2100, 0.8″ before 1000 (now the catalogue's Sgr A* record; residual = catalogue vs the oracle's star file), Kugler/Huber 0.012″, Hipparchos 0.007″ | — | Constant offsets between published values and the oracle's internal ones; documented, not fitted. |
| D2a | Sidereal houses for projection modes (intersection with the plane of t0) | S | Swiss differs from "tropical − ayanamsa" by 6–43″; we subtract. |
| D2b | `PREC_ORIG` residual | — | Kinoshita (1975) Newcomb now used (SLALIB/PAL coefficients): 0.6–0.7″ before 1000, ≤ 0.013″ after 1800. |

### E. Events

| # | Item | Size | Notes |
|---|---|---|---|
| E1 | Heliacal events: `heliacal_ut`, `heliacal_pheno_ut`, `vis_limit_mag` (heliacal rising/setting, evening first/morning last, visibility limit) | L | **Done 2026-10-29** (Schaefer 1998 model; dates ±1–3 d vs Swiss, see HELIACAL_REPORT.md). |
| E2 | Lunar occultations of planets and stars: `lun_occult_when_glob/when_loc/where` | M | **Done** (2026-10-29, ws/occultations): IOTA graze limits ≤ 0.31″ / ≤ 4 s; Swiss 1681 global (0 missed) and 549 local events (flags identical). ECLIPSES_REPORT "Lunar occultations". |
| E3 | Local lunar eclipse search `lun_eclipse_when_loc` | S | **Done** (2026-10-29): same flags as Swiss in 132/132 cases. |
| E4 | `helio_cross(_ut)` (heliocentric longitude crossing), `mooncross_node(_ut)` | S | **Done** (2026-10-29). |
| E5 | Eclipse extras: saros series and member numbers, local apparent noon, umbra/penumbra limits in `sol_eclipse_where` | S–M | **Done** (2026-10-29): saros = NASA for 11 495 eclipses, = Swiss 911/911; noon 0.11 s p50; umbral limits 0.4 km p50 vs NASA path tables (Swiss returns none). |
| E6 | Star magnitudes `fixstar(2)_mag`, `fixstars_mag` | S | V magnitudes are in the Hipparcos data. |

### F. Time and utilities

| # | Item | Size | Notes |
|---|---|---|---|
| F1 | UTC with leap seconds: `utc_to_jd`, `jdet_to_utc`, `jdut1_to_utc`, `utc_time_zone` | S | **Done** (2026-10-29): IERS Bulletin C 72 table. |
| F2 | `deltat_ex`, `set_delta_t_userdef`, `get_tid_acc`, `set_tid_acc` | S | **Done** (2026-10-29): `DeltaTModel.user_seconds` / `lunar_ndot` + `Engine::with_delta_t`; ContextVar-scoped in `swe_compat`. |
| F3 | `date_conversion`, `day_of_week`, `sidtime0` | S | **Done** (2026-10-29). |
| F4 | Formatting and normalisation helpers: `split_deg`, `cs2*str`, `csnorm`, `degnorm`, `radnorm`, `difdeg*`, `difcs*`, `d2l`, `csroundsec` | S | **Done** (2026-10-29): identical to Swiss on 6 800 `split_deg` and 200 centisecond cases. |
| F5 | `get_planet_name`, `version`, `get_library_path`, `get_current_file_data` | S | **Done** (2026-10-29): FrenchEphemeris equivalents (pack path and coverage, DE number 441). |

### G. Ephemeris range and sources

| # | Item | Size | Notes |
|---|---|---|---|
| G1 | Time range: Swiss covers about −13000 to +17000 with DE431 | M | **Done** (2026-10-29, ws/extended-range): optional refit packs `fepk-de441-ext-past/future-v1.bin` (305 MB, FEPK v2) extend coverage to −13198 … +17190; ≤ 0.001″ geocentric planets, Moon ≤ 0.1 km; see EXTENDED_RANGE_REPORT.md. Open: obliquity/ecliptic-of-date differs from Swiss by 1″ at −1000 … 46″ at −13000. |
| G2 | Choosing another JPL file (`set_jpl_file`) or the Moshier analytical ephemeris (`SEFLG_MOSEPH`) | M | Our builder can pack any DE kernel. An analytical fallback is optional. |

## Suggested order for a "full replacement" claim

1. **Quick wins (≈ 1–2 weeks):** A6, A9, B1–B3, C2, C4, D3, E3, E4, E6, F1–F5. These close most remaining API names cheaply. A6, A9, B1–B3, D3 (current modes), E3, E4, F1–F5 done on ws/quickwins; C2/C4 (ws/houses-all) and E6 (ws/star-names) elsewhere.
2. **What astrologers use most (≈ 3–5 weeks):** C1 + C3 (all house systems), D1a/D2a, A2 (natural Lilith/Priapus), A3 (Uranian and hypothetical bodies), A7–A8, E5.
3. **Bigger items (≈ 4–8 weeks):** A1 (asteroid bundle and on-demand fetch), E2 (occultations), E1 (heliacal events), B4, G1.
4. **Optional:** A4/A5 (moons, comets), G2 (Moshier fallback).
