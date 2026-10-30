# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Offline tests of frenchephemeris.asteroids (no network, no native module needed).

    .venv/bin/python -m pytest tests/test_asteroids.py
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("fe_asteroids", os.path.join(REPO, "python", "frenchephemeris",
                                                                             "asteroids.py"))
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)

DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
EROS_SPK = os.path.join(DATA, "sources", "minor", "hz_433_eros.bsp")


def kepler(jd, a=2.7 * A.AU_KM, e=0.25, period=1600.0):
    """Planar Keplerian ellipse, km."""
    M = 2 * np.pi * (np.asarray(jd) - 2451545.0) / period
    E = M.copy()
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    x = a * (np.cos(E) - e)
    y = a * np.sqrt(1 - e * e) * np.sin(E)
    return np.c_[x, y, 0.1 * y]


def test_fit_meets_tolerance_and_tiles_span():
    series, st = A.fit_body(kepler, 2400000.5, 2410000.5, 1024.0, n=12, tol_km=0.5)
    assert st["max_km"] <= 0.5
    assert st["min_interval_days"] < 1024.0  # refinement happened
    t = 2400000.5
    for s in series:
        assert s["start"] == t
        t = s["start"] + len(s["coef"]) * s["interval"]
    assert t == pytest.approx(2410000.5, abs=1e-9)
    jd = np.random.default_rng(0).uniform(2400000.5, 2410000.5, 2000)
    # evaluate the series directly
    for s in series:
        m = (jd >= s["start"]) & (jd < s["start"] + len(s["coef"]) * s["interval"])
        i = ((jd[m] - s["start"]) // s["interval"]).astype(int)
        u = 2 * (jd[m] - s["start"] - i * s["interval"]) / s["interval"] - 1
        p = np.einsum("rck,rk->rc", s["coef"][i], A._cheb(u, s["coef"].shape[2]))
        assert np.linalg.norm(p - kepler(jd[m]), axis=1).max() <= 0.5


def test_fepk_roundtrip_and_cache_listing(tmp_path):
    series, _ = A.fit_body(kepler, 2400000.5, 2401000.5, 128.0, n=16, tol_km=0.5)
    path = tmp_path / "ast_1181.bin"
    A.write_fepk(str(path), [dict(target=A.naif_id(1181), center=0, start=s["start"], interval=s["interval"],
                                  coef=s["coef"]) for s in series], {"name": "t"})
    f = A.Fepk(str(path))
    jd = np.linspace(2400000.5, 2401000.5, 500)
    assert np.linalg.norm(f.pos(A.naif_id(1181), 0, jd) - kepler(jd), axis=1).max() <= 0.5
    assert A.pack_numbers(str(path)) == {1181}
    (tmp_path / "ast_5.bin").write_bytes(path.read_bytes())
    assert [os.path.basename(p) for p in A.cached_packs(str(tmp_path))] == ["ast_1181.bin", "ast_5.bin"]
    assert [os.path.basename(p) for p in A.cached_packs(str(tmp_path), exclude={5})] == ["ast_1181.bin"]


def test_ids():
    assert A.naif_id(1181) == 2001181
    assert A.naif_id(1_234_567) == 21_234_567
    assert A.base_interval(2.7) in (128.0, 256.0)
    assert A.base_interval(40.0) == 1024.0


@pytest.mark.skipif(not os.path.exists(EROS_SPK), reason="Horizons SPK not present")
def test_type21_reader_matches_tools_spk():
    sys.path.insert(0, os.path.join(REPO, "tools"))
    import spk as S
    a = A.SPK21(EROS_SPK)
    b = S.SPK.open(EROS_SPK)
    jd = np.random.default_rng(1).uniform(2305427, 2598000, 5000)
    p2 = b.state_et(20000433, 10, (jd - 2451545.0) * 86400.0)[0]
    assert np.abs(a.pos(jd) - p2).max() < 1e-6
    assert A.spk_metadata(a)["solution_id"].startswith("JPL#")
