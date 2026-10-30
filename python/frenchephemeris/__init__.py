# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""FrenchEphemeris: JPL-derived, thread-safe ephemeris engine."""

from __future__ import annotations

import os
from pathlib import Path

from ._native import *  # noqa: F401,F403
from ._native import Engine

SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO = range(10)
MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, EARTH = range(10, 15)
# Minor bodies (need the minor-body pack; JPL-covered ranges only)
CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA = range(15, 21)
AST_OFFSET = 10000
ERIS = AST_OFFSET + 136199
EROS = AST_OFFSET + 433

# Main pack first; the 1900-2100 development pack is used only when it is absent.
_DEFAULT_PACK_NAMES = ("fepk-de441-de442-v1.bin", "fepk-de442-1900-2100-v1.bin")
_MINOR_PACK_NAME = "fepk-minor-v2.bin"
# Extended-range refits of DE441 (-13200 .. -1 and 2401 .. 17191), optional
_EXT_PACK_NAMES = ("fepk-de441-ext-past-v1.bin", "fepk-de441-ext-future-v1.bin")
# 2401 .. 2581 only; used when the full future pack is absent (they overlap).
_NEAR_FUTURE_PACK_NAME = "fepk-de441-ext-near-future-v1.bin"
_ASTEROID_PACK_NAME = "fepk-asteroids-v1.bin"
# Bodies of the minor pack (MPC numbers): never taken from asteroid packs.
_MINOR_NUMBERS = {1, 2, 3, 4, 433, 2060, 5145, 136199}


def asteroid_id(mpc_number: int) -> int:
    """Engine / Swiss id of a numbered asteroid (AST_OFFSET + number)."""
    return AST_OFFSET + int(mpc_number)


def pack_dirs() -> list[Path]:
    """Directories searched for packs: `FRENCHEPHEMERIS_DATA` (or the per-user data
    directory filled by `python -m frenchephemeris.packs download`), then a source
    checkout's data/packs."""
    from .packs import data_dir

    return [data_dir(), Path(__file__).resolve().parents[2] / "data" / "packs"]


def default_pack_path() -> Path:
    """`FRENCHEPHEMERIS_PACK` env var, else the first main pack found in `pack_dirs()`."""
    env = os.environ.get("FRENCHEPHEMERIS_PACK")
    if env:
        return Path(env)
    for name in _DEFAULT_PACK_NAMES:
        for d in pack_dirs():
            if (d / name).exists():
                return d / name
    raise FileNotFoundError(
        "no FEPK pack found: run `python -m frenchephemeris.packs download` "
        "(or set FRENCHEPHEMERIS_PACK / FRENCHEPHEMERIS_DATA)")


def default_minor_pack_path() -> Path | None:
    """`FRENCHEPHEMERIS_MINOR_PACK` env var ("" disables), else `fepk-minor-v2.bin`
    next to the main pack, if present."""
    env = os.environ.get("FRENCHEPHEMERIS_MINOR_PACK")
    if env is not None:
        return Path(env) if env else None
    p = default_pack_path().parent / _MINOR_PACK_NAME
    return p if p.exists() else None


def default_ext_pack_paths() -> list[Path]:
    """`FRENCHEPHEMERIS_EXT_PACKS` env var (os.pathsep-separated; "" disables), else the
    extended-range packs next to the main pack that exist."""
    env = os.environ.get("FRENCHEPHEMERIS_EXT_PACKS")
    if env is not None:
        return [Path(p) for p in env.split(os.pathsep) if p]
    d = default_pack_path().parent
    found = [d / n for n in _EXT_PACK_NAMES if (d / n).exists()]
    if not (d / _EXT_PACK_NAMES[1]).exists() and (d / _NEAR_FUTURE_PACK_NAME).exists():
        found.append(d / _NEAR_FUTURE_PACK_NAME)
    return found


def default_star_compat_path() -> Path | None:
    """Optional fixed-star compatibility alias table: `FRENCHEPHEMERIS_STAR_COMPAT` env var
    ("" disables), else `private/stars/star_compat_aliases.tsv` in the repository if present.
    It is not distributed with FrenchEphemeris (see docs/STAR_NAMES.md)."""
    env = os.environ.get("FRENCHEPHEMERIS_STAR_COMPAT")
    if env is not None:
        return Path(env) if env else None
    p = Path(__file__).resolve().parents[2] / "private" / "stars" / "star_compat_aliases.tsv"
    return p if p.exists() else None


def default_asteroid_pack_path() -> Path | None:
    """`FRENCHEPHEMERIS_ASTEROID_PACK` env var ("" disables), else
    `fepk-asteroids-v1.bin` next to the main pack, if present."""
    env = os.environ.get("FRENCHEPHEMERIS_ASTEROID_PACK")
    if env is not None:
        return Path(env) if env else None
    p = default_pack_path().parent / _ASTEROID_PACK_NAME
    return p if p.exists() else None


def default_extra_packs() -> list:
    """Minor pack, asteroid bundle, then on-demand packs from the asteroid cache
    (`FRENCHEPHEMERIS_ASTEROID_CACHE`, see `frenchephemeris.asteroids.fetch`)
    for bodies not already provided (a body may come from one pack only)."""
    from . import asteroids as _ast

    extra, have = [], set(_MINOR_NUMBERS)
    minor = default_minor_pack_path()
    if minor:
        extra.append(str(minor))
    bundle = default_asteroid_pack_path()
    if bundle:
        extra.append(str(bundle))
        have |= _ast.pack_numbers(str(bundle))
    extra += _ast.cached_packs(exclude=have)
    return extra


def open_default(swiss_delta_t: bool = True) -> Engine:
    """Main pack + extended-range packs + minor-body pack + asteroid bundle + cached
    on-demand asteroid packs, whichever are present, and the optional fixed-star
    compatibility table (`default_star_compat_path`).
    Outside the loaded coverage every call raises (deterministic out-of-range)."""
    extra = [str(p) for p in default_ext_pack_paths()] + default_extra_packs()
    compat = default_star_compat_path()
    return Engine(str(default_pack_path()), swiss_delta_t, extra or None,
                  star_compat=str(compat) if compat else None)
