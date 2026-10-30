# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Drop-in subset of the `swisseph` module surface used by FreeAstroAPI.

Migration aid: `import frenchephemeris.swe_compat as swe` lets existing route
code run on the new engine. Unlike the original, `set_sid_mode` / `set_topo`
store their value in a ContextVar, so they are scoped to the current thread or
asyncio task and never leak between concurrent requests. New code should call
`Engine` directly with explicit arguments.
"""

from __future__ import annotations

import contextvars
import functools

from . import _native
from ._native import (  # noqa: F401  (re-exported constants)
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT,
    FLG_SIDEREAL, FLG_SPEED, FLG_TOPOCTR, FLG_TRUEPOS,
    ECL_1ST_VISIBLE, ECL_2ND_VISIBLE, ECL_3RD_VISIBLE, ECL_4TH_VISIBLE, ECL_ALLTYPES_LUNAR,
    ECL_ALLTYPES_SOLAR, ECL_ANNULAR, ECL_ANNULAR_TOTAL, ECL_CENTRAL, ECL_HYBRID, ECL_MAX_VISIBLE,
    ECL_NONCENTRAL, ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL, ECL_VISIBLE,
)

FLG_JPLEPH = 1
FLG_SWIEPH = 2
FLG_MOSEPH = 4
_IGNORED_FLAGS = FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH

SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO = range(10)
MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, EARTH = range(10, 15)
CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA = range(15, 21)
AST_OFFSET = 10000  # Eris = AST_OFFSET + 136199, Eros = AST_OFFSET + 433
INTP_APOG, INTP_PERG = 21, 22
FICT_OFFSET = 40
CUPIDO, HADES, ZEUS, KRONOS, APOLLON, ADMETOS, VULKANUS, POSEIDON = range(40, 48)
HARRINGTON, NEPTUNE_LEVERRIER, NEPTUNE_ADAMS, PLUTO_LOWELL = range(50, 54)
WHITE_MOON, WALDEMATH = 56, 58  # 48 Isis, 49 Nibiru, 54 Pickering, 55 Vulcan, 57 Proserpina: not provided

GREG_CAL, JUL_CAL = 1, 0

# Every sidereal mode by its Swiss number (SIDM_LAHIRI = 1, ...; 0-46), the
# user-defined mode and the option bits.
for _i, _ident in enumerate(_native.ayanamsa_idents()):
    globals()["SIDM_" + _ident.upper()] = _i
SIDM_GALCENT_RGBRAND = 30  # C-header spelling (pyswisseph: SIDM_GALCENT_RGILBRAND)
SIDM_USER = 255
SIDBIT_ECL_T0 = _native.SIDBIT_ECL_T0
SIDBIT_SSY_PLANE = _native.SIDBIT_SSY_PLANE
SIDBIT_USER_UT = _native.SIDBIT_USER_UT
SIDBIT_ECL_DATE = _native.SIDBIT_ECL_DATE
SIDBIT_NO_PREC_OFFSET = _native.SIDBIT_NO_PREC_OFFSET
SIDBIT_PREC_ORIG = _native.SIDBIT_PREC_ORIG
_SIDBITS_ALL = (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE | SIDBIT_USER_UT | SIDBIT_ECL_DATE
                | SIDBIT_NO_PREC_OFFSET | SIDBIT_PREC_ORIG)

version = f"FrenchEphemeris {_native.VERSION}"


class Error(Exception):
    pass


# (mode | SIDBIT_* bits, t0, ayan_t0), as passed to set_sid_mode.
_sid_mode: contextvars.ContextVar[tuple[int, float, float]] = contextvars.ContextVar(
    "fe_sid_mode", default=(SIDM_FAGAN_BRADLEY, 0.0, 0.0))
_topo: contextvars.ContextVar[tuple[float, float, float] | None] = contextvars.ContextVar("fe_topo", default=None)


@functools.lru_cache(maxsize=1)
def _base_engine() -> _native.Engine:
    from . import open_default

    return open_default()


def _engine() -> _native.Engine:
    """The shared engine, or a derived one (same pack) carrying this context's Delta T settings."""
    user, tid = _delta_t_user.get(), _tid_acc.get()
    if user is None and tid is None:
        return _base_engine()
    return _derived_engine(user, tid)


def set_ephe_path(path: str | None = None) -> None:
    """Accepted for compatibility; the engine reads its own pack."""


def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None:
    """Sidereal mode (+ SIDBIT_* option bits) for the current thread / task.
    `t0` (JD TT, or UT with SIDBIT_USER_UT) and `ayan_t0` are used only with
    SIDM_USER, as in the classic API."""
    base = mode & 0xFF
    if mode & ~(0xFF | _SIDBITS_ALL):
        raise Error(f"unsupported sidereal bits in {mode}")
    if base != SIDM_USER and not get_ayanamsa_name(base):
        raise Error(f"unsupported sidereal mode {mode}")
    _sid_mode.set((mode, float(t0), float(ayan_t0)))


def set_topo(lon: float, lat: float, alt: float = 0.0) -> None:
    _topo.set((lon, lat, alt))


def _sid_name() -> str:
    """Sidereal spec string understood by the native `sidereal=` argument."""
    mode, t0, ayan_t0 = _sid_mode.get()
    bits = mode & ~0xFF
    if mode & 0xFF == SIDM_USER:
        return f"user:{t0!r}:{ayan_t0!r}|{bits}"
    return str(mode)


def _sid_arg():
    """Sidereal mode for the native functions that accept an int: the int itself (mode +
    SIDBIT_* bits), which the binding reads without a failed conversion; a spec string only
    for SIDM_USER. Much cheaper per call, especially with several threads."""
    mode = _sid_mode.get()[0]
    return _sid_name() if mode & 0xFF == SIDM_USER else mode


def get_ayanamsa_name(sid_mode: int) -> str:
    """Name of a sidereal mode ("" for unknown numbers and SIDM_USER)."""
    return _native.get_ayanamsa_name(sid_mode)


def _wrap(fn, *args, **kw):
    try:
        return fn(*args, **kw)
    except (_native.EphemerisError, ValueError) as e:
        raise Error(str(e)) from e


def calc_ut(jd_ut: float, body: int, flags: int = FLG_SWIEPH | FLG_SPEED):
    f = flags & ~_IGNORED_FLAGS
    sid = _sid_arg() if f & FLG_SIDEREAL else None
    r = _wrap(_engine().calc_ut, jd_ut, body, f, sidereal=sid, topo=_topo.get())
    return tuple(r), flags


def calc(jd_tt: float, body: int, flags: int = FLG_SWIEPH | FLG_SPEED):
    f = flags & ~_IGNORED_FLAGS
    sid = _sid_arg() if f & FLG_SIDEREAL else None
    r = _wrap(_engine().calc_tt, jd_tt, body, f, sidereal=sid, topo=_topo.get())
    return tuple(r), flags


def _hsys(h) -> str:
    return h.decode() if isinstance(h, (bytes, bytearray)) else str(h)


# Sunshine house_pos needs the Sun's declination; the original keeps the one of
# its last houses call in a process global. Here: a ContextVar set by
# houses_armc(..., ascmc9) / houses_armc_ex2 (per thread / asyncio task).
_sun_dec: contextvars.ContextVar[float] = contextvars.ContextVar("fe_sun_dec", default=0.0)


def _houses_ex2(jd_ut, lat, lon, hsys, sid):
    return _wrap(_native.houses_ex2, _engine(), jd_ut, lat, lon, _hsys(hsys), sid)


def houses(jd_ut: float, lat: float, lon: float, hsys=b"P"):
    cusps, ascmc, _, _ = _houses_ex2(jd_ut, lat, lon, hsys, None)
    return tuple(cusps), tuple(ascmc)


def houses_ex(jd_ut: float, lat: float, lon: float, hsys=b"P", flags: int = 0):
    sid = _sid_name() if flags & FLG_SIDEREAL else None
    cusps, ascmc, _, _ = _houses_ex2(jd_ut, lat, lon, hsys, sid)
    return tuple(cusps), tuple(ascmc)


def houses_ex2(jd_ut: float, lat: float, lon: float, hsys=b"P", flags: int = 0):
    """(cusps, ascmc, cusp_speeds, ascmc_speeds); speeds in deg/day."""
    sid = _sid_name() if flags & FLG_SIDEREAL else None
    return tuple(tuple(x) for x in _houses_ex2(jd_ut, lat, lon, hsys, sid))


def houses_armc(armc: float, lat: float, eps: float, hsys=b"P", ascmc9: float = 0.0):
    _sun_dec.set(ascmc9)
    cusps, ascmc = _wrap(_native.houses_armc, armc, lat, eps, _hsys(hsys), ascmc9)
    return tuple(cusps), tuple(ascmc)


def houses_armc_ex2(armc: float, lat: float, eps: float, hsys=b"P", ascmc9: float = 0.0):
    _sun_dec.set(ascmc9)
    return tuple(tuple(x) for x in _wrap(_native.houses_armc_ex2, armc, lat, eps, _hsys(hsys), ascmc9))


def house_pos(armc: float, lat: float, eps: float, objcoord, hsys=b"P") -> float:
    lon, lat_ecl = objcoord[0], objcoord[1]
    return _wrap(_native.house_pos, armc, lat, eps, lon, lat_ecl, _hsys(hsys), _sun_dec.get())


def house_name(hsys) -> str:
    return _native.house_name(_hsys(hsys))


def gauquelin_sector(jd_ut: float, body: int, method: int, geopos, atpress: float = 0.0, attemp: float = 0.0,
                     flags: int = FLG_SWIEPH | FLG_TOPOCTR) -> float:
    """Fixed stars (string body) are not supported here."""
    if isinstance(body, str):
        raise Error("gauquelin_sector for fixed stars is not supported")
    return _wrap(_native.gauquelin_sector, _engine(), jd_ut, body, method, tuple(geopos), atpress, attemp,
                 flags & ~_IGNORED_FLAGS)


def get_ayanamsa_ut(jd_ut: float) -> float:
    """Mean ayanamsa (no nutation), as pyswisseph's get_ayanamsa_ut."""
    return _wrap(_engine().ayanamsa_ut, jd_ut, _sid_arg(), False)


def get_ayanamsa_ex_ut(jd_ut: float, flags: int = FLG_SWIEPH):
    """(retflag, ayanamsa): mean with FLG_NONUT, else mean + nutation in longitude."""
    return flags, _wrap(_engine().ayanamsa_ut, jd_ut, _sid_arg(), not flags & FLG_NONUT)


def get_ayanamsa(jd_tt: float) -> float:
    """Mean ayanamsa (no nutation) at JD TT."""
    return _wrap(_engine().ayanamsa_tt, jd_tt, _sid_arg(), False)


def get_ayanamsa_ex(jd_tt: float, flags: int = FLG_SWIEPH):
    """(retflag, ayanamsa) at JD TT: mean with FLG_NONUT, else true."""
    return flags, _wrap(_engine().ayanamsa_tt, jd_tt, _sid_arg(), not flags & FLG_NONUT)


def julday(year: int, month: int, day: int, hour: float = 12.0, cal: int = GREG_CAL) -> float:
    return _native.julday(year, month, day, hour, cal == GREG_CAL)


def revjul(jd: float, cal: int = GREG_CAL):
    return _native.revjul(jd, cal == GREG_CAL)


def deltat(jd_ut: float) -> float:
    """Delta T (days) as used by the engine's UT functions (honours set_delta_t_userdef /
    set_tid_acc in this context)."""
    return _engine().deltat(jd_ut)


def sidtime(jd_ut: float) -> float:
    return _native.sidtime(jd_ut)


# --- eclipses -------------------------------------------------------------------------------
# Same argument order and return-tuple shapes as pyswisseph; ephemeris flags are ignored
# (the engine always uses its JPL pack). Times are JD UT. Fields not computed are 0.0
# (lunar tret[1], solar tret[1] "local noon", where() limit coordinates) or -99999999
# (saros numbers). See docs/ECLIPSES_REPORT.md for the semantics and validation.


def _geo(geopos):
    lon, lat, *rest = geopos
    return (float(lon), float(lat), float(rest[0]) if rest else 0.0)


def lun_eclipse_when(tjdut: float, flags: int = FLG_SWIEPH, ecltype: int = 0, backwards: bool = False):
    f, tret = _wrap(_native.lun_eclipse_when, _engine(), tjdut, ecltype, bool(backwards))
    return f, tuple(tret)


def lun_eclipse_how(tjdut: float, geopos, flags: int = FLG_SWIEPH):
    f, attr = _wrap(_native.lun_eclipse_how, _engine(), tjdut, _geo(geopos))
    return f, tuple(attr)


def sol_eclipse_when_glob(tjdut: float, flags: int = FLG_SWIEPH, ecltype: int = 0, backwards: bool = False):
    f, tret = _wrap(_native.sol_eclipse_when_glob, _engine(), tjdut, ecltype, bool(backwards))
    return f, tuple(tret)


def sol_eclipse_where(tjdut: float, flags: int = FLG_SWIEPH):
    f, geopos, attr = _wrap(_native.sol_eclipse_where, _engine(), tjdut)
    return f, tuple(geopos), tuple(attr)


def sol_eclipse_how(tjdut: float, geopos, flags: int = FLG_SWIEPH):
    f, attr = _wrap(_native.sol_eclipse_how, _engine(), tjdut, _geo(geopos))
    return f, tuple(attr)


def sol_eclipse_when_loc(tjdut: float, geopos, flags: int = FLG_SWIEPH, backwards: bool = False):
    f, tret, attr = _wrap(_native.sol_eclipse_when_loc, _engine(), tjdut, _geo(geopos), bool(backwards))
    return f, tuple(tret), tuple(attr)


# --- lunar occultations ----------------------------------------------------------------------
# `body` is a planet id (SUN, MERCURY..PLUTO) or a fixed-star name. ECL_ONE_TRY tests only
# the next conjunction (result flags 0 and tret[0] = conjunction if none); like pyswisseph
# (and the C API's `backward |= SE_ECL_ONE_TRY`) it is read from `backwards`, and also
# accepted in `flags`.

from ._native import ECL_OCC_BEG_DAYLIGHT, ECL_OCC_END_DAYLIGHT, ECL_ONE_TRY  # noqa: E402,F401
# Horizons-compatible frames for calc/calc_ut/fixstar*: IAU 1976/1980 + IERS EOP pole offsets.
from ._native import FLG_DPSIDEPS_1980, FLG_JPLHOR, FLG_JPLHOR_APPROX  # noqa: E402,F401


def _body(body):
    return body.decode() if isinstance(body, (bytes, bytearray)) else body


def lun_occult_when_glob(tjdut: float, body, flags: int = FLG_SWIEPH, ecltype: int = 0, backwards: bool = False):
    bk = int(backwards)
    f, tret = _wrap(_native.lun_occult_when_glob, _engine(), tjdut, _body(body), (int(flags) | bk) & ECL_ONE_TRY,
                    int(ecltype), bool(bk & 1))
    return f, tuple(tret)


def lun_occult_when_loc(tjdut: float, body, geopos, flags: int = FLG_SWIEPH, backwards: bool = False):
    bk = int(backwards)
    f, tret, attr = _wrap(_native.lun_occult_when_loc, _engine(), tjdut, _body(body), _geo(geopos),
                          (int(flags) | bk) & ECL_ONE_TRY, bool(bk & 1))
    return f, tuple(tret), tuple(attr)


def lun_occult_where(tjdut: float, body, flags: int = FLG_SWIEPH):
    f, geopos, attr = _wrap(_native.lun_occult_where, _engine(), tjdut, _body(body))
    return f, tuple(geopos), tuple(attr)
# --- events (rise/set/transit, phenomena, crossings, equation of time, coordinates) ---------

from ._native import (  # noqa: E402,F401  (re-exported constants)
    APP_TO_TRUE, BIT_ASTRO_TWILIGHT, BIT_CIVIL_TWILIGHT, BIT_DISC_BOTTOM, BIT_DISC_CENTER,
    BIT_FIXED_DISC_SIZE, BIT_FORCE_SLOW_METHOD, BIT_GEOCTR_NO_ECL_LAT, BIT_HINDU_RISING,
    BIT_NAUTIC_TWILIGHT, BIT_NO_REFRACTION, CALC_ITRANSIT, CALC_MTRANSIT, CALC_RISE, CALC_SET,
    ECL2HOR, EQU2HOR, TRUE_TO_APP,
)


def _geopos(geopos):
    g = tuple(float(x) for x in geopos)
    return g if len(g) == 3 else (g[0], g[1], 0.0)


def rise_trans_true_hor(tjdut: float, body, rsmi: int, geopos, atpress: float = 0.0, attemp: float = 0.0,
                        horhgt: float = 0.0, flags: int = FLG_SWIEPH):
    """(res, tret[10]): res 0 = found (tret[0] = JD UT), -2 = circumpolar / no event in the window.

    `flags` is accepted for compatibility; like the original, the computation always uses the
    apparent topocentric place for `geopos` (the global `set_topo` value is not consulted).
    """
    if isinstance(body, (str, bytes)):
        raise NotImplementedError("rise_trans for fixed stars is not implemented yet in FrenchEphemeris")
    res, t, _cat = _wrap(_native.rise_trans, _engine(), tjdut, body, rsmi, _geopos(geopos), atpress, attemp, horhgt)
    return res, (t,) + (0.0,) * 9


def rise_trans(tjdut: float, body, rsmi: int, geopos, atpress: float = 0.0, attemp: float = 0.0,
               flags: int = FLG_SWIEPH):
    return rise_trans_true_hor(tjdut, body, rsmi, geopos, atpress, attemp, 0.0, flags)


def pheno_ut(tjdut: float, planet: int, flags: int = FLG_SWIEPH):
    f = flags & ~_IGNORED_FLAGS
    return tuple(_wrap(_native.pheno_ut, _engine(), tjdut, planet, f, _topo.get() if f & FLG_TOPOCTR else None))


def pheno(tjdet: float, planet: int, flags: int = FLG_SWIEPH):
    f = flags & ~_IGNORED_FLAGS
    return tuple(_wrap(_native.pheno, _engine(), tjdet, planet, f, _topo.get() if f & FLG_TOPOCTR else None))


def _cross(fn, x2cross, tjd, flags):
    f = flags & ~_IGNORED_FLAGS
    return _wrap(fn, _engine(), x2cross, tjd, f, _sid_name() if f & FLG_SIDEREAL else None)


def solcross_ut(x2cross: float, tjdut: float, flags: int = FLG_SWIEPH) -> float:
    return _cross(_native.solcross_ut, x2cross, tjdut, flags)


def solcross(x2cross: float, tjdet: float, flags: int = FLG_SWIEPH) -> float:
    return _cross(_native.solcross, x2cross, tjdet, flags)


def mooncross_ut(x2cross: float, tjdut: float, flags: int = FLG_SWIEPH) -> float:
    return _cross(_native.mooncross_ut, x2cross, tjdut, flags)


def mooncross(x2cross: float, tjdet: float, flags: int = FLG_SWIEPH) -> float:
    return _cross(_native.mooncross, x2cross, tjdet, flags)


def time_equ(tjdut: float) -> float:
    """Equation of time in days (local apparent minus local mean time)."""
    return _wrap(_native.time_equ, _engine(), tjdut)


def lmt_to_lat(tjdlmt: float, geolon: float) -> float:
    return _wrap(_native.lmt_to_lat, _engine(), tjdlmt, geolon)


def lat_to_lmt(tjdlat: float, geolon: float) -> float:
    return _wrap(_native.lat_to_lmt, _engine(), tjdlat, geolon)


def cotrans(coord, eps: float):
    c = tuple(float(x) for x in coord)
    return _native.cotrans((c[0], c[1], c[2] if len(c) > 2 else 1.0), eps)


def cotrans_sp(coord, eps: float):
    return tuple(_native.cotrans_sp([float(x) for x in coord][:6], eps))


def azalt(tjdut: float, flag: int, geopos, atpress: float, attemp: float, xin):
    return _wrap(_native.azalt, _engine(), tjdut, flag, _geopos(geopos), atpress, attemp, [float(x) for x in xin])


def refrac(alt: float, atpress: float, attemp: float, flag: int, model: str = "saemundsson") -> float:
    """Simple refraction conversion. The original simple `refrac` follows Saemundsson (1986)
    below ~15 deg (measured, < 0.0001 deg), so that is the compat default here; the engine
    default elsewhere (rise/set, azalt) is Sinclair (1982) / Explanatory Supplement 1992,
    which is what the original's extended refraction uses. Pass model="sinclair1982" for it.
    """
    return _wrap(_native.refrac, alt, atpress, attemp, flag, model)


def fixstar2_ut(star: str, tjdut: float, flags: int = FLG_SWIEPH):
    """Returns (xx[6], "Name,nomenclature", retflag) like pyswisseph."""
    f = flags & ~_IGNORED_FLAGS
    sid = _sid_name() if f & FLG_SIDEREAL else None
    xx, name = _wrap(_native.fixstar_ut, _engine(), star, tjdut, f, sidereal=sid, topo=_topo.get())
    return tuple(xx), name, flags


def fixstar2(star: str, tjd_tt: float, flags: int = FLG_SWIEPH):
    return fixstar2_ut(star, _native.tt_to_ut(tjd_tt) if hasattr(_native, "tt_to_ut") else tjd_tt - deltat(tjd_tt), flags)


def fixstar2_mag(star: str):
    """Returns (mag, name) like pyswisseph; 999.99 if no V magnitude is known. Names follow the
    engine's star naming (the compatibility table is loaded by open_default when available)."""
    return _wrap(_native.fixstar_mag, _engine(), star)


fixstar_mag = fixstar2_mag


# --- heliacal events (Schaefer 1998 visibility model) -------------------------------------

from ._native import (  # noqa: E402,F401  (re-exported constants)
    ACRONYCHAL_RISING, ACRONYCHAL_SETTING, COSMICAL_SETTING, EVENING_FIRST, EVENING_LAST,
    HELFLAG_AV, HELFLAG_HIGH_PRECISION, HELFLAG_LONG_SEARCH, HELFLAG_NO_DETAILS, HELFLAG_OPTICAL_PARAMS,
    HELFLAG_SEARCH_1_PERIOD, HELFLAG_VISLIM_DARK, HELFLAG_VISLIM_NOMOON, HELFLAG_VISLIM_PHOTOPIC,
    HELIACAL_RISING, HELIACAL_SETTING, MIXEDOPIC_FLAG, MORNING_FIRST, MORNING_LAST, PHOTOPIC_FLAG,
    SCOTOPIC_FLAG,
)

def _hel_args(geopos, atmo, observer):
    g = tuple(float(x) for x in geopos)
    g = g if len(g) == 3 else (g[0], g[1], 0.0)
    return g, [float(x) for x in atmo], [float(x) for x in observer]


def heliacal_ut(tjdut: float, geopos, atmo, observer, objname: str, eventtype: int, flags: int = FLG_SWIEPH):
    """(start, optimum, end) of visibility on the event day, JD UT. The visibility model is
    Schaefer (1998) as published; dates typically agree with the original within 1-3 days."""
    g, a, o = _hel_args(geopos, atmo, observer)
    return _wrap(_native.heliacal_ut, _engine(), tjdut, g, a, o, objname, eventtype, flags & ~_IGNORED_FLAGS)


def heliacal_pheno_ut(tjdut: float, geopos, atmo, observer, objname: str, eventtype: int, flags: int = FLG_SWIEPH):
    g, a, o = _hel_args(geopos, atmo, observer)
    return tuple(_wrap(_native.heliacal_pheno_ut, _engine(), tjdut, g, a, o, objname, eventtype, flags & ~_IGNORED_FLAGS))


def vis_limit_mag(tjdut: float, geopos, atmo, observer, objname: str, flags: int = FLG_SWIEPH):
    g, a, o = _hel_args(geopos, atmo, observer)
    res, dret = _wrap(_native.vis_limit_mag, _engine(), tjdut, g, a, o, objname, flags & ~_IGNORED_FLAGS)
    return float(res), tuple(dret)


def _unsupported(name):
    def f(*_a, **_k):
        raise NotImplementedError(f"swe_compat.{name} is not implemented yet in FrenchEphemeris")

    f.__name__ = name
    return f




# --- gap-analysis quick wins: flags, ECL_NUT, calc_pctr, crossings, local lunar eclipses,
# UTC, Delta T controls, calendar and formatting helpers, names ------------------------------
# Delta T settings are ContextVar-scoped like set_sid_mode: a derived engine (sharing the
# pack) is used while they are set, so they never leak between threads or asyncio tasks.

from ._native import (  # noqa: E402,F401  (re-exported constants)
    ECL_NUT, ECL_PARTBEG_VISIBLE, ECL_PARTEND_VISIBLE, ECL_PENUMBBEG_VISIBLE, ECL_PENUMBEND_VISIBLE,
    ECL_TOTBEG_VISIBLE, ECL_TOTEND_VISIBLE, FLG_ASTROMETRIC, FLG_ICRS, FLG_RADIANS, FLG_SPEED3, FLG_XYZ,
    SPLIT_DEG_KEEP_DEG, SPLIT_DEG_KEEP_SIGN, SPLIT_DEG_NAKSHATRA, SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_SEC, SPLIT_DEG_ZODIACAL,
)

DELTAT_AUTOMATIC = -1e-10
TIDAL_AUTOMATIC = 999999.0
TIDAL_DE200, TIDAL_DE403, TIDAL_DE404, TIDAL_DE405, TIDAL_DE406 = -23.8946, -25.580, -25.580, -25.826, -25.826
TIDAL_DE421, TIDAL_DE422, TIDAL_DE430, TIDAL_DE431, TIDAL_DE441 = -25.85, -25.85, -25.82, -25.80, -25.936
TIDAL_26, TIDAL_STEPHENSON_2016 = -26.0, -25.85
TIDAL_DEFAULT = TIDAL_SWIEPH = TIDAL_JPLEPH = TIDAL_DE431
TIDAL_MOSEPH = TIDAL_DE404

_delta_t_user: contextvars.ContextVar[float | None] = contextvars.ContextVar("fe_delta_t_user", default=None)
_tid_acc: contextvars.ContextVar[float | None] = contextvars.ContextVar("fe_tid_acc", default=None)


@functools.lru_cache(maxsize=64)
def _derived_engine(user: float | None, tid: float | None) -> _native.Engine:
    return _native.engine_with_delta_t(_base_engine(), user, tid)


def set_delta_t_userdef(dt: float) -> None:
    """Fixed Delta T (days) for every UT computation in this context; DELTAT_AUTOMATIC resets."""
    _delta_t_user.set(None if dt == DELTAT_AUTOMATIC else float(dt))


def set_tid_acc(acc: float) -> None:
    """Lunar tidal acceleration (″/cy²) for Delta T in this context; TIDAL_AUTOMATIC resets to
    the engine's (DE441, −25.936)."""
    _tid_acc.set(None if acc == TIDAL_AUTOMATIC else float(acc))


def get_tid_acc() -> float:
    return _native.tid_acc(_engine())


def deltat_ex(tjdut: float, flags: int = FLG_SWIEPH) -> float:
    """Delta T (days) for the engine's ephemeris; the ephemeris flag is accepted and ignored
    (there is one, JPL-derived, ephemeris)."""
    return _engine().deltat(tjdut)


def calc_pctr(tjdet: float, planet: int, center: int, flags: int = FLG_SWIEPH | FLG_SPEED):
    f = flags & ~_IGNORED_FLAGS
    r = _wrap(_native.calc_pctr, _engine(), tjdet, planet, center, f, _sid_name() if f & FLG_SIDEREAL else None)
    return tuple(r), flags


def helio_cross(planet: int, x2cross: float, tjdet: float, flags: int = FLG_SWIEPH, backwards: bool = False) -> float:
    f = flags & ~_IGNORED_FLAGS
    sid = _sid_name() if f & FLG_SIDEREAL else None
    return _wrap(_native.helio_cross, _engine(), planet, x2cross, tjdet, f, sid, bool(backwards))


def helio_cross_ut(planet: int, x2cross: float, tjdut: float, flags: int = FLG_SWIEPH, backwards: bool = False) -> float:
    f = flags & ~_IGNORED_FLAGS
    sid = _sid_name() if f & FLG_SIDEREAL else None
    return _wrap(_native.helio_cross_ut, _engine(), planet, x2cross, tjdut, f, sid, bool(backwards))


def mooncross_node(tjdet: float, flags: int = FLG_SWIEPH):
    f = flags & ~_IGNORED_FLAGS
    return tuple(_wrap(_native.mooncross_node, _engine(), tjdet, f, _sid_name() if f & FLG_SIDEREAL else None))


def mooncross_node_ut(tjdut: float, flags: int = FLG_SWIEPH):
    """UT in, UT out (as documented; pyswisseph 2.10.03 returns the TT-variant instant)."""
    f = flags & ~_IGNORED_FLAGS
    return tuple(_wrap(_native.mooncross_node_ut, _engine(), tjdut, f, _sid_name() if f & FLG_SIDEREAL else None))


def lun_eclipse_when_loc(tjdut: float, geopos, flags: int = FLG_SWIEPH, backwards: bool = False):
    f, tret, attr = _wrap(_native.lun_eclipse_when_loc, _engine(), tjdut, _geo(geopos), bool(backwards))
    return f, tuple(tret), tuple(attr)


def _cal(cal: int) -> bool:
    if cal not in (GREG_CAL, JUL_CAL):
        raise ValueError("cal must be GREG_CAL or JUL_CAL")
    return cal == GREG_CAL


def utc_to_jd(year: int, month: int, day: int, hour: int, minutes: int, seconds: float, cal: int = GREG_CAL):
    """(jd_et, jd_ut1) from UTC (leap seconds per IERS Bulletin C; UT1 before 1972)."""
    return tuple(_wrap(_native.utc_to_jd, _engine(), year, month, day, hour, minutes, float(seconds), _cal(cal)))


def jdet_to_utc(tjdet: float, cal: int = GREG_CAL):
    return tuple(_native.jdet_to_utc(_engine(), tjdet, _cal(cal)))


def jdut1_to_utc(tjdut: float, cal: int = GREG_CAL):
    return tuple(_native.jdut1_to_utc(_engine(), tjdut, _cal(cal)))


def utc_time_zone(year: int, month: int, day: int, hour: int, minutes: int, seconds: float, offset: float):
    return tuple(_native.utc_time_zone(year, month, day, hour, minutes, float(seconds), float(offset)))


def date_conversion(year: int, month: int, day: int, hour: float = 12.0, cal: bytes = b"g"):
    c = cal.decode() if isinstance(cal, (bytes, bytearray)) else str(cal)
    if len(c) != 1:
        raise TypeError("cal must be exactly 1 byte")
    if c not in "gj":
        raise ValueError("cal must be b'g' or b'j'")
    valid, jd, dt = _native.date_conversion(year, month, day, float(hour), c == "g")
    return valid, jd, tuple(dt)


def day_of_week(jd: float) -> int:
    return _native.day_of_week(jd)


def sidtime0(tjdut: float, eps: float, nutation: float) -> float:
    return _native.sidtime0(_engine(), tjdut, eps, nutation)


def split_deg(degree: float, roundflag: int):
    return tuple(_native.split_deg(degree, roundflag))


def _ch(x) -> str:
    c = x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
    if len(c) != 1:
        raise TypeError("parameter must be exactly 1 byte")
    return c


def cs2degstr(cs: int) -> str:
    return _native.cs2degstr(cs)


def cs2lonlatstr(cs: int, plus, minus) -> str:
    return _native.cs2lonlatstr(cs, _ch(plus), _ch(minus))


def cs2timestr(cs: int, sep, suppresszero: bool = False) -> str:
    return _native.cs2timestr(cs, _ch(sep), bool(suppresszero))


from ._native import (  # noqa: E402,F401  (pure helpers, same signatures)
    csnorm, csroundsec, d2l, degnorm, difcs2n, difcsn, difdeg2n, difdegn, difrad2n, radnorm,
)


def get_planet_name(planet: int) -> str:
    return _native.planet_name(planet)



def get_library_path() -> str:
    return _native.__file__


def get_current_file_data(ifno: int):
    """(path, start, end, denum) of the data pack used: 0 planets / 1 Moon (main JPL pack,
    DE441/DE442 -> 441), 2 main asteroids / 3 other minor bodies (minor pack, if loaded),
    4 stars (built in: empty, as the original reports for star files). Start/end are the
    pack's TDB coverage."""
    from . import default_minor_pack_path, default_pack_path

    naif = {0: 399, 1: 301, 2: 2000001, 3: 2002060}.get(ifno)
    if naif is None:
        return "", 0.0, 0.0, 0
    path = default_pack_path() if ifno < 2 else default_minor_pack_path()
    cov = _native.coverage(_base_engine(), naif)
    if path is None or cov is None:
        return "", 0.0, 0.0, 0
    return str(path), cov[0], cov[1], 441
