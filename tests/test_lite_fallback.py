# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
import pytest

import frenchephemeris as fe

lite = fe.lite_pack_paths()
needs_lite = pytest.mark.skipif(not lite, reason="lite packs not present (run tools/fetch_bundled_packs.py)")


@needs_lite
def test_open_default_covers_0_to_2400():
    eng = fe.open_default()
    for jd in (1830000.5, 2086308.0, 2451545.0):  # 298, 1000, 2000
        lon = eng.calc_ut(jd, 1)[0]
        assert 0.0 <= lon < 360.0


@needs_lite
def test_lite_close_to_exact_where_both_exist():
    exact = fe.open_default()
    only_lite = fe.Engine(str(lite[0]), True, [str(p) for p in lite[1:]])
    jd = 2451545.0
    for body in range(10):
        a, b = exact.calc_ut(jd, body)[0], only_lite.calc_ut(jd, body)[0]
        assert abs(((a - b) + 180) % 360 - 180) * 3600 < 0.5
