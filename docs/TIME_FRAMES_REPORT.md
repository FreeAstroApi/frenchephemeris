# Time scales and reference frames: report

Workstream: fundamental astronomy (`crates/fe-core/src/time.rs`, `crates/fe-core/src/frames.rs`).
Every number below was produced by the code and scripts in this repository. Swiss Ephemeris
(`pyswisseph 2.10.3.2`) was used only as a black-box oracle through its public Python API.
No Swiss source code was read.

## 1. What was built

| Item | Where |
|---|---|
| `julday` / `revjul` (Gregorian and Julian, astronomical years) | `time.rs` |
| Delta T, versioned (`DELTA_T_MODEL = "FE-DT/1 …"`), `DeltaTModel::{DEFAULT, SWISS_COMPAT}`, `delta_t_seconds_with`, `tidal_correction_seconds`, `ut_to_tt[_with]`, `tt_to_ut[_with]` | `time.rs`, tables in `time/tables.rs` |
| TDB − TT: Fairhead & Bretagnon truncated to 87 terms (`tdb_minus_tt_seconds`) and the full 787-term series (`tdb_minus_tt_seconds_full`) | `time.rs` |
| Precession: IAU 2006 (Fukushima–Williams `pfw06`/`fw2m`) and Vondrák 2011 (`ltp_ecliptic_pole`, `ltp_equator_pole`, `ltp`, `ltpb`); frame bias (`frame_bias`) | `frames.rs` |
| Nutation: IAU 2000A (`nut00a`), IAU 2000A + IAU 2006 adjustments (`nut06a`), IAU 2000B (`nut00b`), Swiss-compatible variant (`nut_swiss_compat`), `nutation_matrix` | `frames.rs`, tables in `frames/tables.rs` |
| Obliquity: `obl06`, Vondrák pole angle (`ltp(..).1`), Vondrák 2011 published ε_A series (`vondrak_obliquity`, `Obliquity::Vondrak2011`), `obliquity_mean[_with]` | `frames.rs` |
| Earth rotation: `era`, `gmst` (IAU 2006), `gast` (IAU 2006/2000A, ERA − EO), `s06`, `eors`, `eect00`, `equation_of_equinoxes`; long-term model `SiderealTime::SwissLongTerm` (`gmst_long_term`, `emb_mean_longitude_j2000`, `long_term_offsets`) | `frames.rs` |
| `FrameSet::new` (default) / `FrameSet::with_model(jd, FrameModel)`: `bp`, `bpn`, `ecl_mean`, `ecl_true`, `eps_*`, `dpsi/deps`, `eo`, `eqeq`, `gast(jd_ut)`, `gmst(jd_ut)`, `nutation_matrix()`, `bpn_rate()`, `ecl_true_rate()` | `frames.rs` |
| `gcrs_to_ecliptic_j2000` = R1(84381.406″) · frame bias | `frames.rs` |
| Table generator (ERFA coefficients pinned to commit `1a8044cd`, Delta T sources) | `crates/fe-core/tables/gen_tables.py`, `tables/sources/` |
| Fixture generators and fixtures | `crates/fe-core/tests/fixtures/{time,frames}/gen_*.py` and `*.json` |
| Integration tests | `crates/fe-core/tests/time.rs`, `crates/fe-core/tests/frames.rs` |

Model presets (`FrameModel`):

| preset | precession | nutation | mean obliquity |
|---|---|---|---|
| `DEFAULT` | Vondrák 2011 + IAU 2006 bias | IAU 2000A + IAU 2006 adjustments | angle between the Vondrák poles |
| `VONDRAK_2000B` | Vondrák 2011 + bias | IAU 2000B | Vondrák poles |
| `IAU2006_2000A` | IAU 2006 (equivalent to ERFA `pnm06a`) | IAU 2000A/2006 | `obl06` |
| `IAU2006_2000B` | IAU 2006 | IAU 2000B | `obl06` |
| `SWISS_COMPAT` | Vondrák 2011 + bias | Swiss 2.10 variant (below) | Vondrák 2011 ε_A series |
| `context::API_FRAME_MODEL` (engine/API default) | Vondrák 2011 + bias | Swiss 2.10 variant | Vondrák 2011 ε_A series |

The free functions `nutation` and `obliquity_mean` use the `DEFAULT` model. `gmst` and `gast`
are the IAU 2006 / 2006-2000A definitions (ERFA `gmst06` and `gst06a`). `FrameSet::gast` and
`FrameSet::gmst` follow the frame's own model.

## 2. Delta T

### Model FE-DT/1 (default)

| range (UT) | source |
|---|---|
| before −720 | SMH 2016 long-term parabola −320 + 32.5·u², u = (Y − 1825)/100, offset so it joins the spline continuously at −720 |
| −720 … 1950.0 | **Addendum 2020** spline: Morrison, Stephenson, Hohenkerk & Zawilski 2021, Proc. R. Soc. A 477:20200776, Table S15.2020 |
| 1950.0 … 1955.0 | the spline plus a linear ramp of the 1955.0 (observed − spline) offset (+0.66 s), so the curve is continuous |
| 1955.0 … 1973.0 | USNO `historic_deltat.data`, semi-annual observed values |
| 1973-02-01 … 2026-10-01 | USNO `deltat.data`, monthly observed values, linearly interpolated |
| after 2026-10-01 | ΔT_last + v·τ + 32.5·τ² s. τ is Julian centuries after the last value. v is the mean rate over the last 12 months (removes the annual term). |

The future term keeps today's ΔT and rate. It adds the long-term acceleration of the SMH
parabola, 65 s/cy², which is a LOD increase of +1.78 ms/cy and matches the +1.8 ms/cy measured
by SMH 2016. The table has no ΔT jump at any boundary. Unit tests check continuity at −720,
1950, 1955 and 2026-10.

Y = 2000 + (JD − 2451545)/365.2425. This decimal-year mapping reproduces Swiss's behaviour to
≤ 0.014 s. With a Julian-year mapping the error is 1.0 s, so this is the mapping Swiss uses.

The spline refers to ṅ = −25.85″/cy², which HMNAO says goes with DE430/LE430. The option
`DeltaTModel::lunar_ndot` applies the standard correction −0.91072·(ṅ − (−25.85))·u²,
with u = (Y − 1955)/100, to the pre-1955 part.

Sources (the files are copied into `crates/fe-core/tables/sources/`):
- SMH 2016 Table S15: figshare doi:10.6084/m9.figshare.4290866.v1 (`https://ndownloader.figshare.com/files/6998699`). Paper doi:10.1098/rspa.2016.0404.
- Addendum 2020 Table S15.2020: figshare doi:10.6084/m9.figshare.13885863.v2 (`https://ndownloader.figshare.com/files/26513569`). Paper doi:10.1098/rspa.2020.0776.
- USNO: `https://maia.usno.navy.mil/ser7/deltat.data`, `https://maia.usno.navy.mil/ser7/historic_deltat.data` (downloaded 2026-10-28).
- HMNAO data page (blocked from this machine): `http://astro.ukho.gov.uk/nao/lvm/`. The ṅ = −25.85 statement is also in `https://hpiers.obspm.fr/eoppc/series/longterm/hmnao_a_2021.txt`.

### What Swiss 2.10 does (black box)

- **Before 1950:** `swe_deltat` = SMH **2016** spline + tidal correction from ṅ = −25.85 to Swiss's default −25.80 (`swe.get_tid_acc()`). A least-squares fit of swe − SMH2016 gives −0.04551·u². That is −0.91072 × (+0.04997), so Δṅ = 0.050. The residual is ≤ 0.012 s.
  - **`DeltaTModel::SWISS_COMPAT`** (SMH 2016, ṅ = −25.80) reproduces Swiss to **0.0057 s max before 1950**. This is tested.
- **1955 … 2022:** Swiss uses observed values. They agree with the USNO tables to ≤ 0.04 s (1955–1972) and ≤ 0.065 s (1973–2022).
- **From about 2023:** Swiss uses its own stale predictions. For example 68.80 s for 2027–2028, against an observed 69.2 s.
- **Further out:** Swiss extrapolates with a curve that grows more slowly than ours (details below).

### Swiss after its last tabulated value (black box, 2026-10 investigation)

The Swiss documentation (swisseph.htm §7 "ΔT") says: *"For the time after the last tabulated
value, we use the formula of Stephenson (1997; p. 507), with a modification that avoids a
discontinuity at the end of the tabulated period. A linear term is added that makes a slow
transition from the table to the formula over a period of 100 years."* Measured (`swe.deltat`,
Y in Julian years from J2000, residuals of each fit ≤ 1e-11 s):

| interval | Swiss ΔT |
|---|---|
| to 2028 | its table (incl. its own predictions, 68.8 s flat 2027–2028) |
| 2028 … 2128 | the 2128–2500 cubic + a linear term that vanishes at 2128 (the documented 100-year transition) |
| 2128 … 2500 | a cubic in x = (Y − 2000)/100: 64.0002 + 17.3669x + 8.0002x² + 4.0333x³ s |
| after 2500 | a parabola 42.5 + 32.5x² s (32.5 s/cy² is the SMH 2016 acceleration; the apex at 2000 is not SMH's 1825) |

The structure (table → 100-year linear blend → formula) matches the documentation, but the
formula is not any published one available here. Residuals of the published candidates (each
blended linearly from the last USNO value, 69.21 s at 2026.75, over 100 years, as documented)
against Swiss:

| year | Swiss (s) | FE (last value + rate + 32.5 s/cy²) | SMH 2016 parabola −320 + 32.5u² | Morrison & Stephenson 2004 −20 + 32u² | −20 + 31u² (the earlier 31 s/cy² of Stephenson & Morrison, as cited by SMH 2016) | Espenak & Meeus 2006 |
|---:|---:|---:|---:|---:|---:|---:|
| 2050 | 74.6 | −0.9 | −33 | +38 | +36 | +18 |
| 2100 | 93.2 | +1.9 | −99 | +125 | +118 | +110 |
| 2200 | 163.0 | +23.8 | −26 | +279 | +265 | +279 |
| 2300 | 297.0 | +46.4 | +116 | +420 | +397 | +420 |
| 2400 | 519.6 | +45.5 | +235 | +537 | +503 | +537 |

(u = (Y − 1825)/100 resp. (Y − 1820)/100.) Stephenson (1997, *Historical Eclipses and Earth's
Rotation*, CUP) p. 507 could not be consulted, and no published formula reproduces the
measured cubic/parabola. Our existing extrapolation is already the closest published-method
curve (≤ 46 s, Moon ≤ 25″ to 2400), so it is kept and **`SWISS_COMPAT` is unchanged after
2026**: reproducing Swiss there would mean shipping the six coefficients above, which were
obtained by fitting Swiss output (forbidden by the project rule).

### FE default vs Swiss

Moon Δλ ≈ ΔΔT × 0.000152°/s, shown in arcsec. Every 10 years from 0 to 2400 is produced by
`cargo test --release --test time delta_t_report -- --ignored --nocapture`. Excerpt:

| year | swe 2.10 (s) | FE default (s) | default − swe (s) | Moon Δλ (″) | SWISS_COMPAT − swe (s) |
|---:|---:|---:|---:|---:|---:|
| −1 | 10567.279 | 10452.039 | −115.240 | −63.06 | −0.006 |
| 0 | 10556.896 | 10441.313 | −115.584 | −63.25 | −0.006 |
| 100 | 9540.369 | 9405.044 | −135.325 | −74.05 | −0.003 |
| 200 | 8553.871 | 8424.699 | −129.173 | −70.68 | −0.001 |
| 300 | 7578.009 | 7476.111 | −101.898 | −55.76 | +0.002 |
| 400 | 6593.388 | 6535.116 | −58.272 | −31.89 | +0.005 |
| 500 | 5590.100 | 5586.601 | −3.500 | −1.92 | +0.004 |
| 600 | 4596.150 | 4651.655 | +55.505 | +30.37 | +0.002 |
| 700 | 3649.026 | 3760.420 | +111.393 | +60.95 | +0.001 |
| 800 | 2786.218 | 2943.037 | +156.819 | +85.81 | −0.001 |
| 900 | 2045.215 | 2229.648 | +184.433 | +100.92 | −0.002 |
| 1000 | 1463.504 | 1650.393 | +186.889 | +102.27 | −0.003 |
| 1100 | 1062.857 | 1222.881 | +160.024 | +87.57 | −0.002 |
| 1200 | 802.173 | 914.611 | +112.438 | +61.53 | −0.001 |
| 1300 | 624.631 | 681.149 | +56.518 | +30.93 | −0.001 |
| 1400 | 473.411 | 482.288 | +8.877 | +4.86 | −0.000 |
| 1500 | 291.692 | 292.343 | +0.651 | +0.36 | +0.001 |
| 1600 | 88.805 | 109.127 | +20.322 | +11.12 | +0.001 |
| 1700 | 14.287 | 14.100 | −0.187 | −0.10 | +0.000 |
| 1800 | 18.605 | 18.367 | −0.238 | −0.13 | −0.000 |
| 1900 | −1.990 | −1.977 | +0.013 | +0.01 | −0.001 |
| 1950 | 28.932 | 28.932 | −0.000 | −0.00 | −0.000 |
| 2000 | 63.829 | 63.829 | +0.000 | +0.00 | +0.000 |
| 2020 | 69.361 | 69.361 | +0.000 | +0.00 | +0.000 |
| 2026 | 68.900 | 69.110 | +0.210 | +0.11 | +0.210 |
| 2028 | 68.800 | 69.356 | +0.556 | +0.30 | +0.556 |
| 2050 | 74.580 | 73.650 | −0.931 | −0.51 | −0.931 |
| 2100 | 93.182 | 95.107 | +1.926 | +1.05 | +1.926 |
| 2200 | 162.997 | 186.771 | +23.773 | +13.01 | +23.773 |
| 2300 | 296.992 | 343.431 | +46.439 | +25.41 | +46.439 |
| 2400 | 519.581 | 565.089 | +45.508 | +24.90 | +45.508 |

Summary:
- **Before 1955:** the default differs from Swiss only because it uses the newer **2020 addendum** spline and does not rescale for ṅ. The largest difference is **189 s near 960 AD**, which moves the Moon by **≈ 0.029° (104″)**. Typical differences are 50–135 s in 0–1300, and 20 s around 1600.
  - The tidal part alone is −0.0455·u² s (−17 s at year 0, −4 s at 1000).
  - `SWISS_COMPAT` removes the whole pre-1950 difference (≤ 0.006 s).
- **1955–2022:** both use observed data. They agree to ≤ 0.07 s, and to ≤ 0.02 s at most dates.
- **After 2023:** both are predictions or extrapolations. The difference is +0.2 to +0.6 s in 2026–2028 (Swiss's predictions are stale), +1.9 s at 2100, and +46 s at 2300–2400. That is at most 0.007° of Moon longitude. `SWISS_COMPAT` does not reproduce Swiss after 1955.

## 3. TDB − TT

Fairhead & Bretagnon (1990) series with the JPL-mass adjustments of SOFA/ERFA `dtdb`,
geocentric (no topocentric terms). Compared with `erfa.dtdb` at 600 epochs across years 0–2400:

| function | terms | max error | time (release) |
|---|---:|---:|---:|
| `tdb_minus_tt_seconds_full` | 787 | 1.3e-14 s | – |
| `tdb_minus_tt_seconds` | 87 (amplitude × 2^power ≥ 3e-8 s) | 0.32 µs (0.42 µs on a 20 000-point scan) | ≈ 0.3–0.5 µs |

## 4. Calendars

`julday` and `revjul` use exact integer arithmetic for the day number. On 3408 cases (years
−50…2400, both calendars, Feb-29 overflow, fractional hours) they are **bit-identical to
`swe.julday` / `swe.revjul`**, for the JD and for the returned hour. They also match the
Meeus ch. 7 reference dates.

Note: serde_json's default float parser is not always correctly rounded, so fixtures carry
exact values as u64 bit patterns (`*_bits`).

## 5. Frames: validation against ERFA (pyerfa 2.0.1.5)

42 epochs across years −1 … 2400 (`tests/frames.rs::erfa_reference_values`):

| quantity | max error |
|---|---:|
| `obl06`, `pfw06`, frame bias `rb`, `bp06` rbp | ≤ 2.2e-16 rad |
| `nut00a`, `nut06a`, `nut00b` | ≤ 4.0e-16 rad (0.0001 µas) |
| `pnm06a` (IAU2006_2000A bpn) | 3.9e-16 |
| `fw2m(pfw06 + nut00b)` (IAU2006_2000B bpn) | 3.9e-16 |
| `ltpecl`, `ltpequ`, `ltp`, `ltpb` (Vondrák) | ≤ 3.3e-16 |
| `s06` | 1.2e-17 rad |
| `eect00` | 3e-21 rad |
| `ee06a` | 8.9e-16 rad |
| `era00`, `gmst06`, `gst06a`, `gst06` | 4.5e-13 rad (0.09 µas). This is the f64 rounding of ERA's large intermediate value, not a model difference. |

ERFA `pnm00b` uses IAU 2000 (not 2006) precession. It is stored in the fixture, but there is
no IAU 2000 precession model here. `IAU2006_2000B` is validated against
`fw2m(pfw06 + nut00b)` instead.

`bpn_rate` (central difference, h = 0.05 d, cheap 77-term nutation) agrees with a numerical
derivative of the full 2000A matrix to < 1e-9 rad/day. The precession rate is about 6.7e-7 rad/day.

## 6. Precession: IAU 2006 vs Vondrák 2011

The IAU 2006 polynomial is intended for a few centuries around J2000. Vondrák, Capitaine &
Wallace (2011, A&A 534, A22; corrigendum A&A 541, C1) fit ±200 kyr and agree with IAU 2006
near J2000. Differences (bias included in both):

| year | total rotation (″) | equinox longitude on ecliptic (″) | obliquity obl06 − Vondrák poles (″) | equator pole (″) |
|---:|---:|---:|---:|---:|
| 0 | 0.370 | −0.341 | +0.196 | 0.211 |
| 500 | 0.156 | −0.154 | +0.035 | 0.080 |
| 1000 | 0.057 | −0.055 | −0.009 | 0.025 |
| 1500 | 0.013 | −0.011 | −0.007 | 0.000 |
| 1800 | 0.002 | −0.001 | −0.002 | 0.000 |
| 2000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 2100 | 0.0002 | −0.0002 | −0.0001 | 0.000 |
| 2400 | 0.005 | −0.001 | −0.005 | 0.000 |

At most 0.37″ (1e-4°) over the API range, so either model is acceptable astrologically.
Vondrák is the one designed for −20 centuries and the one Swiss uses.

## 7. Which frames Swiss 2.10 uses (black box)

Method: `swe.calc` with `FLG_XYZ | FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL` for 10 bodies in
several frames:
- ICRS: `FLG_J2000 | FLG_EQUATORIAL | FLG_ICRS`
- J2000: `FLG_J2000 | FLG_EQUATORIAL`
- mean of date: `FLG_EQUATORIAL | FLG_NONUT`
- true of date: `FLG_EQUATORIAL`
- ecliptic of date: default flags
- J2000 ecliptic: `FLG_J2000`

The orthogonal Procrustes fit recovers each frame matrix with residual ≤ 8e-16. The fits were
done at years 0…2400 with Moshier, and checked with SWIEPH at 2000. Findings:

| Swiss component | identification | agreement |
|---|---|---|
| J2000 (no `FLG_ICRS`) vs ICRS | IAU 2006 frame bias (= ERFA `bp06` rb) | 1 µas |
| precession (ICRS → mean of date) | **Vondrák 2011 with bias, exactly ERFA `ltpb`** | 0.8 µas at all epochs. IAU 2006 is 370 mas away at year 0. |
| nutation | **the 77 IAU 2000B luni-solar terms evaluated with the full IERS 2003 (quartic) Delaunay arguments, with no 2000B planetary offsets (−0.135/+0.388 mas) and no IAU 2006 adjustment** | 3e-8 ″ (0.03 µas). ERFA `nut00b` itself is 130 mas away at year 500, because its linear arguments are only valid near J2000. |
| J2000 ecliptic | R1(84381.406″) · bias | 0.0002 µas |
| mean obliquity of date | **the published Vondrák, Capitaine & Wallace (2011) ε_A series** (Eq. 10 + Table 3), not the pole angle (§7.2) | 5e-11″ over −13000 … +17000 |
| `swe.sidtime` | equals `gst06a` to ≤ 1.2 mas inside **1850-01-01 … 2049-12-31**. Outside, the documented long-term method (§7.1), reproduced to ≤ 0.0007″ up to two constants. It jumps by 0.1″ at 1850 and **1.9″ at 2050-01-01**, and depends on ΔT at 0.0411″/s (the mean Sun is evaluated at TT). | Swiss − IAU: +327″ at year 0, +172″ at 500, +37″ at 1000, +7.7″ at 1500, −1.8″ at 2100, +8.6″ at 2400 |

`FrameModel::SWISS_COMPAT` reproduces the Swiss frames (tested in `tests/frames.rs::swiss_compat_frames`):
- `bp`: 0.8 µas
- `bpn`: 1.4 µas
- J2000 ecliptic: exact
- true ecliptic of date: 0.74 µas (was 18 mas with `obl06` before the ε_A series was implemented)

`DEFAULT` vs Swiss:
- `bp`: identical (0.8 µas)
- `bpn`: ≤ 2.2 mas (nutation: Δψ ≤ 2.5 mas, Δε ≤ 1.6 mas over 0–2400)
- true ecliptic of date: ≤ 0.19″ (0.2″ obliquity difference at year 0, almost all of it in latitude)

Sidereal time: FE `gast` / `FrameSet::gast` are the IAU 2006 values by default
(`SiderealTime::Iau2006`). Vondrák-model GAST differs from `gst06a` by ≤ 0.30″ at year 0.
`SiderealTime::SwissLongTerm` (§7.1) reproduces Swiss's long-term sidereal time.

### 7.1 Long-term sidereal time (`SiderealTime::SwissLongTerm`)

**What Swiss documents** (swisseph.htm §6.1 "Sidereal Time"): *"Since version 2.00 it uses
sidereal time based on the IAU2006/2000 precession/nutation model. As this solution is not
good for the whole time range of JPL Ephemeris DE431, we only use it between 1850 and 2050.
Outside this period, we use a solution based on the long term precession model Vondrak 2011,
nutation IAU2000 and the mean motion of the Earth according to Simon & alii 1994. To make the
function continuous we add some constant values to our long-term function before 1850 and
after 2050."* (Changelog 2.01/2.02 in swephprg.htm: bug fixes of this long-term function.)

**Method implemented** (`frames::gmst_long_term`, `FrameSet::gmst/gast`), chosen among
constructions of that description by black-box comparison:

1. Mean Sun: λ☉ = λ_EMB + 180°, with the Earth–Moon-barycentre mean longitude referred to the
   fixed J2000 ecliptic and equinox of Simon et al. (1994, A&A 282, 663, §3.5):
   λ_EMB = 100°27′59.22059″ + 1295977422.83429″t − 2.04411″t² − 0.00523″t³, t in Julian
   millennia from J2000, **evaluated at TT**.
2. The point (λ☉, β = 0) on the J2000 ecliptic is rotated to the mean equator of J2000
   (ε₀ = 84381.406″), precessed with the Vondrák, Capitaine & Wallace (2011) matrix (`ltp`),
   and rotated to the mean ecliptic of date with the Vondrák 2011 ε_A series
   (`vondrak_obliquity`; originally `obl06`, see §7.2); its longitude L is taken.
3. GMST = L + 360°·frac(UT − J2000) + C; GAST = GMST + Δψ cos ε (true obliquity, no
   complementary terms — measured: `swe.sidtime0(jd, ε, Δψ) = swe.sidtime(jd)` exactly and
   `sidtime0(jd, ε, 0)` differs by Δψ cos ε).
4. C: one constant before 1850-01-01 0h UT and one from 2050-01-01 0h UT, **computed here**
   as IAU 2006 GMST − (1+2) at the boundary instant with ΔT from `DeltaTModel::SWISS_COMPAT`
   (−21.9295″ and −23.5094″). Swiss's constants are not published and are not used.

Discriminating checks (black box, ΔT set to 0 with `set_delta_t_userdef`): Swiss − (1+2+3)
is constant to 0.0001″ over −500…1850 and 2050…4000 (0.03″ at −3000, where Swiss's mean
obliquity departs from `obl06` — resolved by §7.2). Alternatives fail: Simon's mean longitude *of date*
(3″ drift), Simon J2000 + Vondrák general precession p_A added in longitude (no ecliptic
rotation; 1.3″ drift plus a semi-annual term growing as π_A², ±4″ at −2000), the same
with IAU 2006 p_A (0.6″ drift), or a UT-evaluated mean Sun (433″ at year 0: the measured
0.0411″/s ΔT dependence is the solar mean motion, i.e. the mean Sun is taken at TT).

**Residuals vs `swe.sidtime`** (`tests/sidtime.rs`, fixture `tests/fixtures/sidtime/gen_sidtime.py`,
709 epochs −600…2400 plus boundary straddles, both at Swiss's ΔT):

| years | n | SwissLongTerm − swe (″) | IAU 2006 − swe (″, max) |
|---|---:|---:|---:|
| −600 … 0 | 138 | −0.0991 … −0.0984 | 562 |
| 0 … 1000 | 237 | ≤ 0.099 | 326 |
| 1000 … 1850 | 206 | ≤ 0.099 | 37 |
| 1850 … 2050 | 49 | ≤ 0.0006 (both IAU) | 0.001 |
| 2050 … 2400 | 79 | **+1.948** (constant) | 8.7 |

The shape of Swiss's function is reproduced to ≤ 0.0007″ on each side; what remains is the
difference between our continuity constants and Swiss's unpublished ones: −0.099″ before 1850
and +1.948″ from 2050. Ours make the function continuous at both limits (tested, < 0.005″);
Swiss's are not (its 0.1″ and 1.9″ jumps), so matching them would need two numbers read off
Swiss output. 1.95″ of ARMC is ≤ 0.0006° in MC/ASC.

**Houses end to end** (`Engine::houses_ut`, Placidus, 120 random charts −600…2400, |lat| ≤ 60°,
each engine's own ΔT):

| config | ARMC max (before 2050) | cusps/ASC/MC max (before 2050) | ARMC / cusps (after 2050) |
|---|---:|---:|---:|
| `SWISS_DELTA_T` (IAU sidereal time) | 543″ | 0.200° | 8.0″ / 0.0023° |
| `SWISS_TIME` (long-term sidereal time) | 0.10″ | 0.00023° | 4.0″ / 0.0014° (1.95″ constant + ΔT extrapolation difference) |

**Rise/set** (long-range events fixture, years 0–1900, `SWISS_TIME`, no correction applied):
unrefracted Sun ≤ 0.14 s, planets ≤ 0.28 s, Moon ≤ 2.5 s (lunar ephemeris at ancient dates);
before, the raw differences were up to 27 s. Refracted events keep the refraction-model
residual (≤ 11 s here; docs/EVENTS_REPORT.md).

Wiring: `FrameModel.sidereal` (`SiderealTime::{Iau2006, SwissLongTerm}`), so every user of
`FrameSet::gast` follows it — ARMC/houses, topocentric observer positions (`ephem`), rise/set,
transits, equation of time, eclipse local circumstances. `FrameModel::SWISS_COMPAT` and
`EngineConfig::SWISS_TIME` (= `SWISS_DELTA_T` + long-term ST) select it; `EngineConfig::DEFAULT`
and `API_FRAME_MODEL` stay IAU 2006. Python: `Engine(path, swiss_delta_t=True,
swiss_sidereal_time=True)` (both default True). Events now also use the engine's ΔT model
(`EpochInterp::new_with`) instead of the default one. Cost: ~0.3 µs per GAST outside 1850–2050.

### 7.2 Mean obliquity: Vondrák et al. (2011) ε_A series

Source: Vondrák, Capitaine & Wallace 2011, A&A 534, A22, §4.3 Eq. (10) and Table 3 (ε_A
columns). The open-access copy is at https://hal.science/hal-03785916 (the A&A site serves the
same article behind a bot check). The 2012 corrigendum (A&A 541, C1) changes Table 1 only.

ε_A = 84028.206305″ + 0.3624445″T − 0.00004039″T² − 110×10⁻⁹″T³ + Σ (C cos 2πT/P + S sin 2πT/P),
with ten (P, C, S) terms (P from 203 to 4043 cy) and T in Julian centuries TT from J2000.
It is implemented as `frames::vondrak_obliquity` / `Obliquity::Vondrak2011`. At T = 0 it gives
84381.406″ to the 1e-6″ rounding of the coefficients. It stays within 0.025″ of `obl06` over
years 0 … 3000 (`tests/obliquity.rs`).

Black-box test against Swiss over −13000 … +17000: 193 dates, `swe.calc(SE_ECL_NUT)`, `swe.sidtime`
and `swe.houses` 'P' at 47°N; fixture `tests/fixtures/obliquity/gen_obliquity.py`.

| quantity | Vondrák pole angle (before) | ε_A series (now) |
|---|---:|---:|
| mean obliquity − Swiss, max | 66.6″ (−1.04″ at −1000, −31.7″ at −8000, +66.2″ at 14000) | 5e-11″ |
| true obliquity − Swiss, max | — | 2.6e-10″ |
| Placidus cusps/Asc, Swiss ARMC and ΔT, max | 110″ | 6e-5″ |

**Swiss 2.10's mean obliquity is exactly this series.** It is now used by `API_FRAME_MODEL`,
which is the engine/API default for the ecliptic of date, houses and ECL_NUT, and by `SWISS_COMPAT`.
The precession matrix stays the Vondrák equator/ecliptic pole construction, as in Swiss, so only
the obliquity used for the ecliptic of date changes. `FrameModel::DEFAULT` keeps the pole angle,
which is exactly consistent with its matrix. Both are the authors' representations of the same
quantity, and their difference is a fit difference.

The long-term sidereal time (§7.1) projects the mean Sun with the same ε_A. With `obl06` its
residual vs `swe.sidtime` drifted from −0.10″ to −9.5″ at −12000 and from +1.95″ to +19.5″ at
+16000. With ε_A it is constant over the whole range: −0.111 … −0.088″ before 1850 and
+1.944 … +1.952″ from 2050, which are Swiss's unpublished continuity constants. The 0–2400
results of `tests/sidtime.rs` are unchanged.

### Why this default

The choice is (a) scientifically sound across 0–2400 and (b) as close to Swiss as justified:
- **Precession:** Vondrák 2011 is the published long-term IAU-consistent model, and Swiss uses it exactly.
- **Nutation:** IAU 2000A (+2006) is the current IAU/IERS standard. Swiss's 77-term variant differs by ≤ 2.5 mas (7e-7°), so there is no reason to give up the 2000A accuracy. `VONDRAK_2000B` and `SWISS_COMPAT` cost about 1.8 µs per `FrameSet`, against 9 µs for 2000A, for hot loops.
- **Obliquity:** `DEFAULT` keeps the Vondrák pole angle, which is exactly consistent with its precession matrix. The engine/API model (`API_FRAME_MODEL`) and `SWISS_COMPAT` use the authors' separately published ε_A series, which is what Swiss reports. The two Vondrák representations of the same quantity differ by 0.2″ at year 0, 1″ at −1000 and 66″ at +14000 (§7.2).

## 8. Performance (release, Apple Silicon, `cargo test --release --test frames -- --ignored --nocapture`)

| model | `FrameSet::with_model` | `bpn_rate` |
|---|---:|---:|
| DEFAULT (Vondrák + 2000A) | ≈ 9.4 µs | ≈ 1.8 µs |
| VONDRAK_2000B | ≈ 1.8 µs | ≈ 1.8 µs |
| IAU2006_2000A | ≈ 8.6 µs | ≈ 1.4 µs |
| IAU2006_2000B | ≈ 1.3 µs | ≈ 1.4 µs |
| SWISS_COMPAT | ≈ 1.7 µs | ≈ 1.6 µs |

Other timings: `tdb_minus_tt_seconds` ≈ 0.3 µs; `delta_t_seconds` ≈ 4 ns (table lookup).

The 2000A series builds each term's sine/cosine from tables of e^{ikα} by complex products, so
it calls no trig function per term. This made it 2× faster than one `sin_cos` per term, with
identical results (4e-16 rad vs ERFA). About 1 µs of each `FrameSet` is `s06` (needed for GAST).

## 9. Open items

- Swiss's long-term sidereal time is emulated by its documented method (§7.1) except for its two unpublished continuity constants (residual −0.099″ before 1850, +1.948″ from 2050). Swiss's mean obliquity is now reproduced exactly (§7.2).
- `SWISS_COMPAT` Delta T does not reproduce Swiss's post-2023 predictions or its far-future curve (§2: the measured cubic/parabola is not a published formula available here).
- No IAU 2000 (Lieske + rate corrections) precession, so ERFA `pnm00b` is not matched (it is not needed by the pipeline).
- The USNO monthly table ends 2026-10-01. Refresh `tables/sources/USNO-deltat.data`, rerun `tables/gen_tables.py`, and bump `DELTA_T_MODEL` when updating.

## 10. Quick wins (2026-10-29): UTC, Delta T controls, calendar, ECL_NUT, sidtime0

Code: `utc.rs`, `extras.rs`, `time.rs` (`DeltaTModel::user_seconds`), `context.rs`
(`Engine::with_delta_t`). Tests: `tests/quickwins.rs`, fixtures from
`tests/fixtures/quickwins/gen_quickwins.py` (Swiss 2.10 oracle).

**Leap seconds.** TAI − UTC from IERS Bulletin C 72, `Leap_Second.dat` (July 2026, expires
2027-06-28, SHA-256 `6cb6f5d4…68d7`), 28 entries 1972-01-01 … 2017-01-01 (37 s), embedded in
`utc::LEAP_SECONDS` with the citation. Semantics as documented by Swiss (§9.2–9.3): UT1 before
1972; TT = UTC + (TAI − UTC) + 32.184 s from 1972; time treated as UT1 when
ΔT − (TAI − UTC) − 32.184 s > 1 s (an unknown future leap second); 23:59:60.x accepted only on
leap-second days (error otherwise, like the oracle); `jdet_to_utc` returns 60.x inside a leap
second; `utc_time_zone` carries a leap second through unchanged.

| vs Swiss (max) | ≥ 1972 | < 1972 |
|---|---|---|
| `utc_to_jd` TT | 0.00004 s (JD rounding) | 0.24 s (Delta T models; input is UT1) |
| `utc_to_jd` UT1 | 0.32 s (Delta T models, both observed-value based) | 0.00004 s |
| `jdet_to_utc` | 0.00004 s | 0.29 s (Delta T) |
| `jdut1_to_utc` | 0.12 s (Delta T) | 0.00004 s |
| `utc_time_zone`, `date_conversion` (incl. invalid dates, 25 h, month 13), `day_of_week` | identical | |

Valid/invalid classification identical in every `utc_to_jd` fixture case (impossible dates, 24:00, 60 min, 61 s, 23:59:60 on non-leap days). The "treat as UT1
after the last known leap second" switch follows our Delta T extrapolation, so its year
differs from Swiss's (both > 2030).

**Delta T controls (F2), never global state.** `DeltaTModel` gained `user_seconds`
(fixed Delta T for every UT computation); `lunar_ndot` is the tidal acceleration. 
`Engine::with_delta_t(model)` returns an engine sharing the same `Arc<Pack>` (cheap), so
`set_delta_t_userdef` / `set_tid_acc` in `swe_compat` are ContextVar-scoped derived engines.
The tidal correction is the published −0.91072 Δṅ u² (Morrison & Stephenson 2004); Swiss's
changes are 0.07 % smaller (max 0.83 s for DE200 at −500, 0.06 s for DE441). `deltat_ex`
returns the engine's Delta T (ephemeris flag ignored: one JPL ephemeris); it differs from the
oracle's bare `deltat_ex` by up to 75 s at −500 because the oracle uses ṅ = −25.80 (DE431)
there while our engine uses the DE441 value −25.936 (TIME_FRAMES §2).
`get_tid_acc` = −25.936 by default.

**ECL_NUT** (`Engine::ecl_nut`, body −1 in `calc_ut`/`calc_tt`): `[ε, ε̄, Δψ, Δε, 0, 0]` from the
engine's frames. Δψ, Δε identical to Swiss (0.000″); ε and ε̄ identical to Swiss (≤ 3e-10″, −13000 … +17000) since the engine uses the Vondrák ε_A series (§7.2; was 0.2″ at year 0 with the pole angle). RADIANS
converts; XYZ and EQUATORIAL give zeros, as the oracle does.

**sidtime0** (GMST of the engine's model + Δψ cos ε): ≤ 0.0002 s of time vs Swiss in
1850–2050, ≤ 0.007 s before 1850, 0.26 s after 2050 (Swiss's unpublished continuity constant, §9).

**Helpers (F3–F5)** — `split_deg` (all rounding, zodiac, nakshatra and keep flags),
`cs2degstr`, `cs2lonlatstr`, `cs2timestr`, `csnorm`, `csroundsec`, `d2l`, `degnorm`,
`radnorm`, `difdegn`, `difdeg2n`, `difrad2n`, `difcsn`, `difcs2n` (`util.rs`): identical
integer fields / strings on 6 802 `split_deg` and 188 centisecond cases (secfr ≤ 2e-10 s,
degnorm ≤ 1e-12°), including 32-bit wrap-around for out-of-range centiseconds. The
nakshatra width is 13.33333333333333° (40/3 to 16 digits): the oracle's results at exact
boundaries (360° → nakshatra 0) need it. Names: `get_planet_name` for ids 0–22 and the
served asteroids, `get_ayanamsa_name` for the implemented modes.

## 11. IAU 1976/1980 and JPL Horizons modes (B4, 2026-10-29)

Code: `crates/fe-core/src/frames/iau1980.rs`; `Engine::prepare_horizons`; flags `FLG_JPLHOR`
(= `FLG_DPSIDEPS_1980`) and `FLG_JPLHOR_APPROX` in `calc_ut`/`calc`/`fixstar_ut` and
`swe_compat`. Tests: `crates/fe-core/tests/iau1980.rs`; fixtures `tests/fixtures/horizons/`
(`gen_horizons.py`).

- **Models.** IAU 1976 precession (Lieske 1979, eqs. 6–7), IAU 1980 mean obliquity, and IAU 1980
  nutation (Wahr; Seidelmann 1982). The 106-term table comes from the Explanatory Supplement
  1992 Table 3.222.1, transcribed from ERFA `nut80.c` (BSD-3) by `tools/gen_iau1980.py`. Earth
  rotation uses GMST 1982 and the 1994 equation of the equinoxes.
- **ERFA agreement.** nut80 1e-18 rad, pmat76 2e-16, obl80 1e-16, gmst82 2e-12, eqeq94 1e-18.
- **Horizons convention** (IERS Conventions 1996, as Horizons states in its output header:
  "EOP-corrected IAU76/80"). Of-date matrices are N80(Δψ80 + δψ, Δε80 + δε)·P76, with no ICRS
  frame bias.
  - `JPLHOR`: δψ, δε come from the pinned **IERS EOP 14 C04 (IAU 1980)** daily series,
    1962-01-01 … 2026-01-05, from
    https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt
    (sha256 40347f18…a79df). It is stored losslessly in µas as `frames/eop_c04_iau1980.bin`
    (187 kB).
  - Outside the series the first or last value is held, as Horizons does (Swiss manual §2).
  - `JPLHOR_APPROX`: no data file. δψ, δε are the offsets that put the IAU 1976/1980 pole on
    the IAU 2006/2000A CIP; observed offsets differ from these only by the free core nutation.
    They are held constant outside the same span.
- **Validation against JPL Horizons** (quantities 1, 2, 31; Sun, Moon and the Mars/Jupiter/Saturn
  barycentres; 40 epochs over 0–2400; maxima):

  | Frames | 1700–2400 | before 1700 |
  |---|---|---|
  | Default (Vondrák + 2000A) | 1.2″ | 3.0″ |
  | `JPLHOR` / `JPLHOR_APPROX` | 82 / 83 mas | 4.2″ |

  - The 1700–2400 figures include the DE442 (ours) vs DE441 (Horizons) ephemeris difference.
    With the astrometric ICRF difference removed, the frame residual is **≤ 12 mas**.
  - Around the present the agreement is ≤ 2 mas. Example: Mars 2016-10-25 is 288.9700707° /
    −2.0525054°, while the Swiss manual's Horizons-matching values are 288.9700710 / −2.0525051.
- **Before ~1600 Horizons is not plain IAU 1976/1980.**
  - Its ecliptic obliquity departs from obl80 + Δε by an amount growing roughly as T³: 1.5″
    at year 0.
  - Its longitudes depart by an amount growing roughly as T⁴: 4″ at year 0.
  - This matches neither IAU 2006 nor Vondrák. Horizons' method there is not published, and it
    was not fitted.
- **Swiss comparison.** It is not possible with the oracle as installed. pyswisseph with its own
  `.se1` files ignores `FLG_JPLHOR*`: its output equals its default frames, because the modes
  need an original JPL file. So "JPLHOR_APPROX vs Swiss" (2″) only measures Swiss's default
  frames.
