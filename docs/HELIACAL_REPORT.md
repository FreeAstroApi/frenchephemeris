# Heliacal events and visual limiting magnitude (gap item E1)

Module `crates/fe-core/src/heliacal.rs`; Python `crates/fe-py/src/heliacal.rs`
(`heliacal_ut`, `heliacal_pheno_ut`, `vis_limit_mag`, `HELFLAG_*` / event constants) and
`swe_compat.heliacal_ut / heliacal_pheno_ut / vis_limit_mag`. Tests `crates/fe-core/tests/heliacal.rs`,
fixtures `crates/fe-core/tests/fixtures/heliacal/` (`gen_heliacal.py` → `heliacal.json`).
Clean room: the classic API was used only as a black box; no table or constant from it ships.

## 1. Model and provenance

| Part | Source |
|---|---|
| Extinction (Rayleigh 0.1066·e^(−h/8200)(λ/0.55)⁻⁴; aerosol 0.1(λ/0.55)^−1.3 e^(−h/1500)(1 − 0.32/ln RH)^1.33 (1 + sgn(φ) sin((M−3)·30°)); ozone; water vapour 0.031·0.94·RH·e^(T/15)·e^(−h/8200)), gas/aerosol/ozone airmasses | B. E. Schaefer, "To the visual limits", *Sky & Telescope* 95(5), 57 (1998), program VISLIMIT (V-band constants; published listing, as transcribed in L. Bogan's JavaScript port) |
| Sky brightness: dark sky with 11-yr solar cycle, moonlight (Moon magnitude −12.73 + 0.026α + 4e-9α⁴, scattering function), twilight, daylight; min(daylight, twilight) | same |
| Threshold (Blackwell–Hecht form, scotopic/photopic constants at 1500 nL), Snellen ratio term 5 log SN | same |
| Meteorological range VR (`datm[3] ≥ 1`) | Koschmieder β = 3.912/VR (km⁻¹) over Schaefer's 1.5 km aerosol scale height (our choice) |
| `datm[3]` in (0,1): total k | classic semantics (k used with the gas airmass) |
| Optical aids (`HELFLAG_OPTICAL_PARAMS`): monocular factor 1.41, aperture/pupil light gain, sky dilution by magnification; pupil 7 mm·exp(−½(age/100)²) | Schaefer, PASP 102, 212 (1990); Schaefer, Vistas Astron. 36, 311 (1993) |
| Crescent width W = SD(1 − cos ARCL), q = (ARCV − (11.8371 − 6.3226W + 0.7319W² − 0.1018W³))/10, criteria A–E | B. D. Yallop, NAO Technical Note 69 (1997) |
| Defaults: P/T/RH = 1013.25 hPa / 15 °C / 40 % (height-adjusted: standard-atmosphere pressure, −6.5 K/km), age 36, Snellen 1 | classic documented semantics |
| Positions, magnitudes, refraction | our JPL pipeline; `pheno` magnitudes (Mallama & Hilton 2018), SIMBAD V for stars; Sinclair (1982) refraction for AppAltO |

Semantics (from the classic documentation and black-box observation): event types 1 morning
first (heliacal rising), 2 evening last (heliacal setting), 3 evening first and 4 morning last
(Mercury, Venus, Moon; the Moon has only 3/4), 5/6 acronychal → "not provided"; returns
(first visible, optimum, last visible) JD UT of the event day. `HELFLAG_VISLIM_DARK` puts both
Sun and Moon at the nadir (observed: the classic API reports both at −90°),
`HELFLAG_VISLIM_NOMOON` the Moon; for the Moon as target its own light is not added to the sky.
The heliacal search ("VR" method): per local morning/evening, the window from astronomical
twilight (or 4 h) to sunrise / from sunset; visible when the margin (limiting − object
magnitude) > 0 with the object above the true horizon on the rising (morning, east) or
setting (evening, west) side of the meridian and the Sun below the horizon. The event day is
the first visible day after invisible ones (first) or the last visible day before
invisibility (last), found by a 4-day (2 for Mercury/Venus, 1 for the Moon) coarse scan then
day by day; times by bisection (first/last) and golden section (optimum).
`heliacal_pheno_ut` fills 0–14, 16–25, 27 of the 50-value layout (others 0, see the code).
Not implemented: `HELFLAG_AV*` (arcus-visionis methods; VR is always used),
`HELFLAG_NO_DETAILS` returns (first, 0, 0) of the VR result, "mixed-opic" (res 2).

## 2. Validation vs the classic API (black box, pyswisseph 2.10.3.2)

954 `heliacal_ut` cases: Sirius at Memphis / Alexandria / Babylon every 150 years 0–2250
(events 1, 2); Venus and Mercury (1–4), Mars/Jupiter/Saturn (1, 2), Aldebaran, Spica, Regulus,
Antares, Arcturus, Canopus (1, 2), Moon (3, 4) at Memphis, the equator, 45° N, 35° S, 52° N in
years 100, 800, 1500, 2000, 2024, 2300; NOMOON and non-default atmospheres at up to 58° N.
Engine `EngineConfig::SWISS_TIME`.

Event day, FE − classic (days), events in the same apparition:

| group | n | median | mean |Δ| | p90 |Δ| | max |Δ| |
|---|---:|---:|---:|---:|---:|
| Sirius rising, Memphis | 18 | +1 | 0.83 | 1 | 2 |
| Sirius rising, Alexandria | 16 | +1 | 0.94 | 2 | 2 |
| Sirius rising, Babylon | 16 | +1 | 0.94 | 2 | 2 |
| stars, rising | 220 | 0 | 1.76 | 3 | 7 |
| stars, setting | 216 | −1 | 2.26 | 5 | 8 |
| Mars/Jupiter/Saturn rising | 93 | 0 | 2.62 | 5 | 15 |
| Mars/Jupiter/Saturn setting | 89 | −1 | 2.72 | 5 | 20 |
| Mercury/Venus morning first | 63 | 0 | 1.11 | 2 | 7 |
| Mercury/Venus evening last | 58 | 0 | 0.95 | 2 | 8 |
| Mercury/Venus evening first | 58 | −1 | 2.81 | 6 | 9 |
| Mercury/Venus morning last | 58 | +2 | 3.26 | 6 | 8 |
| Moon evening first | 30 | 0 | 0.13 | 1 | 1 |
| Moon morning last | 30 | 0 | 0.17 | 1 | 1 |

When the day agrees, the start of visibility differs by ≤ 0.17 h for planets/stars (median
0.02–0.09 h) and ≤ 0.41 h for the Moon.

Classification: 0 event/no-event differences (24 cases have no event in either engine);
15 / 954 cases (1.6 %) land in a different apparition (|Δ| > 20 d: 110–507 d), all marginal:
Mercury at 45–52° and 35° S, Venus at 52° N, Mars near conjunction, Antares at 45–52° N,
Jupiter at 35° S — the two models disagree on whether that one apparition is visible at all.

Arcus visionis: Sirius rising, topocentric object − Sun altitude at first visibility: FE
12.75° median (12.0–13.4), classic 11.89° (11.3–12.3). Published estimates for Sirius are
~9–11° for the classical (horizon) arcus visionis, depending on transparency (de Jong 2006;
Schaefer 2000); the two definitions differ by the star's altitude at first sighting.

`vis_limit_mag` (400 random cases, 0 above/below-horizon differences; positions agree to
0.015° in altitude, star magnitudes exactly, planets ≤ 0.04 mag, Moon ≤ 3 mag near new Moon):
limiting magnitude FE − classic: dark sky median −0.53 (p90 0.83), twilight +0.21 (p90 0.71),
daylight +0.26 (p90 1.16); outliers of 4–7 mag for objects < 1° above the horizon (airmass
treatment), 4 photopic/scotopic flag differences. The classic implementation's extinction
coefficient is ~0.30–0.40 at sea level with the defaults versus ~0.19–0.47 here (Schaefer's
seasonal aerosol term), and its dark-sky limit is ~0.5 mag fainter — it evidently extends the
published 1998 program (observer and extinction refinements documented by Reijs that we could
not access); we implement the publication as printed.

Speed (debug build): ~9 ms per `heliacal_ut` (one year search), ~0.1 ms per `vis_limit_mag`.

## 3. Open issues
1. Reijs's extensions (Swiss's extinction variant, observer experience/age, AV methods,
   mixed-opic flag) are not implemented; with them the dark-sky 0.5 mag offset and the
   1–3 day date scatter would likely shrink. Needs the ARCHAEOCOSMO documentation text.
2. Acronychal/cosmical events (types 5/6) are rejected as by the classic API.
3. No independent published heliacal table in the 0–2400 range was reachable for an
   external check beyond the Sirius arcus-visionis literature.
