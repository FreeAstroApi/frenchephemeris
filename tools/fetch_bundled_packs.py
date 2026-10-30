#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Download the lossy 0..2400 `lite` packs into python/frenchephemeris/_data/ before a wheel
build (they are bundled in the wheel but not kept in git). Standard library only, so it runs
before the package is built; every file is checked against the release's SHA256SUMS.

    python tools/fetch_bundled_packs.py
"""

import hashlib
import os
import sys
import urllib.request

BASE = "https://github.com/FreeAstroApi/frenchephemeris/releases/download/packs-v1/"
NAMES = ("fepk-de441-lite-0-2400-v1.bin", "fepk-minor-lite-0-2400-v1.bin")
DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python", "frenchephemeris", "_data")


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "frenchephemeris-build"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def main() -> int:
    sums = dict(reversed(line.split()) for line in get(BASE + "SHA256SUMS").decode().splitlines() if line.strip())
    os.makedirs(DEST, exist_ok=True)
    for name in NAMES:
        path = os.path.join(DEST, name)
        if os.path.exists(path) and hashlib.sha256(open(path, "rb").read()).hexdigest() == sums[name]:
            print(f"{name}: present")
            continue
        data = get(BASE + name)
        if hashlib.sha256(data).hexdigest() != sums[name]:
            sys.exit(f"{name}: checksum mismatch")
        with open(path, "wb") as f:
            f.write(data)
        print(f"{name}: {len(data)} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
