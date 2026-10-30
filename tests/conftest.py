# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Without a data pack, tests that need the engine are skipped (not failed).
Get the packs with `python -m frenchephemeris.packs download --set dev`."""

import pytest


def _packs_available() -> bool:
    try:
        import frenchephemeris

        frenchephemeris.default_pack_path()
        return True
    except FileNotFoundError:
        return False


def pytest_configure(config):
    config.addinivalue_line("markers", "no_packs: the test does not need a data pack")


def pytest_collection_modifyitems(config, items):
    if _packs_available():
        return
    skip = pytest.mark.skip(reason="no FEPK pack: run `python -m frenchephemeris.packs download --set dev`")
    for item in items:
        if "no_packs" not in item.keywords:
            item.add_marker(skip)
