# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Download the FEPK data packs (JPL-derived coefficient files) the engine reads.

    python -m frenchephemeris.packs download                # modern set: 1800-2401, exact JPL, ~69 MB
    python -m frenchephemeris.packs download --set full     # -2..2581, exact JPL, ~256 MB
    python -m frenchephemeris.packs download --set dev      # 1900-2100 only, exact, ~27 MB
    python -m frenchephemeris.packs download --set mobile   # 1800-2400, lossy, ~2 MB (apps)
    python -m frenchephemeris.packs download --set compact  # 1800-2400, lossy refit, ~13 MB
    python -m frenchephemeris.packs list
    python -m frenchephemeris.packs path

Packs are published as assets of the `packs-v1` GitHub release, each with a SHA-256
checked after download. They go to `FRENCHEPHEMERIS_DATA` if set, else a per-user data
directory, where `frenchephemeris.open_default()` finds them. Without any download, the
lossy 0..2400 `lite` packs bundled in the package are used; once exact packs are installed
they take over, and the lite packs only cover dates or bodies outside them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO = "FreeAstroApi/frenchephemeris"
RELEASE = "packs-v1"

# name: (bytes, sha256, description)
PACKS = {
    "fepk-de441-de442-v1.bin": (247455672, "9c0a5ddd68fef80d19a4e51771b7aabdfd57ce37e7ea803430fe823d89ee97da",
                                "DE441/DE442 planets, Sun and Moon, years -2 .. 2401 (bit-identical to JPL)"),
    "fepk-minor-v2.bin": (6684336, "752734429253c038deeb1e73a9f0507cc25984034438f7dfe6edeacad919c736",
                          "Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Eros"),
    "fepk-de441-ext-near-future-v1.bin": (1968949, "066c7dd5e374602972da6f2261add6c6bfd56c437037319e06acb3ac72c2808a",
                                          "DE441 refit, years 2401 .. 2581"),
    "fepk-de442-1800-2400-v1.bin": (61947800, "23098485e262eb527923721270936f14ae5e9760a10574e361be02e33b86a68b",
                                    "DE442, years 1800 .. 2401, bit-identical to JPL"),
    "fepk-de442-compact-1800-2400-v1.bin": (6546643, "866cc7c3ceb79943bff7d8d204f8ad1d61fdea94987813a99d90a9c3896c87fc",
                                            "DE442 compact refit, years 1800 .. 2400, LOSSY (planets <= 0.001\", Moon <= 0.04\" of JPL)"),
    "fepk-de441-lite-0-2400-v1.bin": (6890382, "041a66c17ebe48982f432f33ce6f768662436d385ed323bfe6d701e0aaf87e97",
                                      "lite refit, 0 .. 2400, LOSSY (planets <= 0.007\", Moon <= 0.4\"; bundled in the wheel)"),
    "fepk-minor-lite-0-2400-v1.bin": (844755, "e888453939dd3f86f440e4d38167e5ec1ccf5968f0063f6bc9dfe8a718181101",
                                      "minor bodies lite refit, 0 .. 2400, LOSSY (<= 0.013\", Eros <= 0.03\"; bundled)"),
    "fepk-de442-mobile-1800-2400-v1.bin": (1736702, "72b51aa03a31bd62bf3787d58dc9cb7a9bf94d26d908b27b373924e9e8c50ae5",
                                           "DE442 mobile refit, 1800 .. 2400, LOSSY (planets <= 0.007\", Moon <= 0.4\" of JPL)"),
    "fepk-minor-mobile-1800-2400-v1.bin": (291426, "3ca25e4737301d8cb08929c240041c05336b6d253976233e4a1bdcdb533fbfcc",
                                           "minor bodies mobile refit, 1800 .. 2400, LOSSY (<= 0.007\", Eros <= 0.03\")"),
    "fepk-de442-1900-2100-v1.bin": (20750344, "755ef97fbfaa866fe77537584eb29214136146baba9d72c766e7f71c09bc0f5e", "DE442, years 1900 .. 2100 (development pack)"),
    "fepk-de441-ext-past-v1.bin": (144031527, "3894de83aab7d5eb711144b224195e81e5edd07e99cedd4f3904a2730dba2bbf", "DE441 refit, years -13200 .. -2"),
    "fepk-de441-ext-future-v1.bin": (160945165, "c2dc933a66c82203a12e4ff328599cf8c4e079dde2e458218c9b19529014db4c", "DE441 refit, years 2401 .. 17191"),
    "fepk-asteroids-v1.bin": (107936784, "afe56f112fa70566ece7985d3a84b70c0b5e0eb3259ba90d698af4a131da5079", "numbered asteroids bundle"),
}
SETS = {
    "modern": ["fepk-de442-1800-2400-v1.bin", "fepk-minor-v2.bin"],
    "full": ["fepk-de441-de442-v1.bin", "fepk-minor-v2.bin", "fepk-de441-ext-near-future-v1.bin"],
    "compact": ["fepk-de442-compact-1800-2400-v1.bin", "fepk-minor-v2.bin"],
    "mobile": ["fepk-de442-mobile-1800-2400-v1.bin", "fepk-minor-mobile-1800-2400-v1.bin"],
    "lite": ["fepk-de441-lite-0-2400-v1.bin", "fepk-minor-lite-0-2400-v1.bin"],
    "dev": ["fepk-de442-1900-2100-v1.bin", "fepk-minor-v2.bin"],
    "extended": ["fepk-de441-ext-past-v1.bin", "fepk-de441-ext-future-v1.bin"],
    "asteroids": ["fepk-asteroids-v1.bin"],
}
SETS["all"] = sorted({n for names in SETS.values() for n in names})
SETS["standard"] = SETS["full"]  # name before 1.3
CHECKSUMS = "SHA256SUMS"


def data_dir() -> Path:
    """Where downloaded packs live: `FRENCHEPHEMERIS_DATA`, else the platform's per-user
    data directory (~/.local/share, ~/Library/Application Support or %LOCALAPPDATA%)."""
    env = os.environ.get("FRENCHEPHEMERIS_DATA")
    if env:
        return Path(env)
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
    return base / "frenchephemeris" / "packs"


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def _open(url: str, token: str | None, accept: str):
    """GET, with the token only on GitHub's own host (redirects to storage go without it)."""
    headers = {"Accept": accept, "User-Agent": "frenchephemeris"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    try:
        return urllib.request.build_opener(_NoRedirect).open(urllib.request.Request(url, headers=headers), timeout=600)
    except urllib.error.HTTPError as exc:
        if exc.code in (301, 302, 303, 307, 308):
            return urllib.request.urlopen(urllib.request.Request(exc.headers["Location"],
                                                                 headers={"User-Agent": "frenchephemeris"}), timeout=600)
        raise


def _asset_urls(token: str | None) -> dict[str, str]:
    if not token:
        base = f"https://github.com/{REPO}/releases/download/{RELEASE}/"
        return {name: base + name for name in [*PACKS, CHECKSUMS]}
    with _open(f"https://api.github.com/repos/{REPO}/releases/tags/{RELEASE}", token, "application/vnd.github+json") as r:
        release = json.load(r)
    return {a["name"]: f"https://api.github.com/repos/{REPO}/releases/assets/{a['id']}" for a in release["assets"]}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _published_checksums(urls: dict[str, str], token: str | None) -> dict[str, str]:
    if CHECKSUMS not in urls:
        return {}
    with _open(urls[CHECKSUMS], token, "application/octet-stream") as r:
        lines = r.read().decode().splitlines()
    return {line.split()[1]: line.split()[0] for line in lines if len(line.split()) == 2}


def download(pack_set: str = "modern", dest: str | Path | None = None, token: str | None = None,
             quiet: bool = False) -> list[Path]:
    """Download (or verify) the packs of `pack_set`; returns their paths."""
    names = SETS[pack_set]
    dest = Path(dest) if dest else data_dir()
    dest.mkdir(parents=True, exist_ok=True)
    token = token if token is not None else os.environ.get("GITHUB_TOKEN") or None
    urls = _asset_urls(token)
    published = _published_checksums(urls, token)
    out = []
    for name in names:
        size, sha, _ = PACKS[name]
        sha = sha or published.get(name)
        path = dest / name
        if path.exists() and path.stat().st_size == size and (sha is None or _sha256(path) == sha):
            if not quiet:
                print(f"{name}: already present")
            out.append(path)
            continue
        if name not in urls:
            raise RuntimeError(f"{name} is not published in release {RELEASE}")
        tmp = path.with_suffix(path.suffix + ".part")
        digest, done = hashlib.sha256(), 0
        with _open(urls[name], token, "application/octet-stream") as r, tmp.open("wb") as f:
            while chunk := r.read(1 << 20):
                digest.update(chunk)
                f.write(chunk)
                done += len(chunk)
                if not quiet:
                    print(f"\r{name}: {done / 1e6:7.1f} / {size / 1e6:.1f} MB", end="", flush=True)
        if not quiet:
            print()
        if sha and digest.hexdigest() != sha:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"{name}: checksum mismatch")
        tmp.replace(path)
        out.append(path)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m frenchephemeris.packs", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("download", help="download and verify a set of packs")
    d.add_argument("--set", default="modern", choices=sorted(SETS))
    d.add_argument("--dest", help=f"directory (default: {data_dir()})")
    sub.add_parser("list", help="list the published packs and sets")
    sub.add_parser("path", help="print the data directory")
    args = ap.parse_args(argv)
    if args.cmd == "path":
        print(data_dir())
    elif args.cmd == "list":
        for name, (size, _, desc) in PACKS.items():
            print(f"{name:38} {size / 1e6:7.1f} MB  {desc}")
        for s, names in SETS.items():
            print(f"set {s:9} {sum(PACKS[n][0] for n in names) / 1e6:7.1f} MB: {', '.join(names)}")
    else:
        for p in download(args.set, args.dest):
            print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
