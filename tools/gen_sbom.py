# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Generate a software bill of materials (SPDX 2.3, JSON) for a FrenchEphemeris release.

Lists:
  - FrenchEphemeris itself (Rust crates fe-core / fe-py, Python package);
  - every Rust crate in the resolved dependency graph (`cargo metadata`), with version,
    licence and crates.io source;
  - the Python runtime dependencies declared in pyproject.toml, with installed version and
    licence when available;
  - the ephemeris data packs found in data/packs/ (with SHA-256), as data components.

Usage:
    python tools/gen_sbom.py [--out dist/frenchephemeris.spdx.json] [--version 0.1.0]
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.metadata as im
import json
import re
import subprocess
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUPPLIER = "Organization: Gabriel Naulleau EI, trading as FreeAstroAPI (contact@freeastroapi.com)"
PROJECT_LICENSE = "AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial"


def spdx_id(name: str) -> str:
    return "SPDXRef-" + re.sub(r"[^A-Za-z0-9.-]", "-", name)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rust_packages() -> list[dict]:
    meta = json.loads(subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--locked"],
        cwd=ROOT, check=True, capture_output=True, text=True).stdout)
    out = []
    for p in sorted(meta["packages"], key=lambda p: (p["name"], p["version"])):
        local = p.get("source") is None
        out.append({
            "SPDXID": spdx_id(f"rust-{p['name']}-{p['version']}"),
            "name": p["name"],
            "versionInfo": p["version"],
            "supplier": SUPPLIER if local else "NOASSERTION",
            "downloadLocation": "NOASSERTION" if local else f"https://crates.io/crates/{p['name']}/{p['version']}",
            "licenseConcluded": PROJECT_LICENSE if local else (p.get("license") or "NOASSERTION"),
            "licenseDeclared": p.get("license") or "NOASSERTION",
            "copyrightText": "NOASSERTION",
            "filesAnalyzed": False,
            "externalRefs": [] if local else [{
                "referenceCategory": "PACKAGE-MANAGER",
                "referenceType": "purl",
                "referenceLocator": f"pkg:cargo/{p['name']}@{p['version']}",
            }],
        })
    return out


def python_packages() -> list[dict]:
    text = (ROOT / "pyproject.toml").read_text()
    deps = re.findall(r'"([A-Za-z0-9_.-]+)\s*[<>=!~]', text.split("dependencies", 1)[1].split("]", 1)[0])
    out = []
    for name in deps:
        try:
            md = im.metadata(name)
            version = md["Version"]
            lic = md.get("License-Expression") or md.get("License") or "NOASSERTION"
        except im.PackageNotFoundError:
            version, lic = "NOASSERTION", "NOASSERTION"
        out.append({
            "SPDXID": spdx_id(f"py-{name}-{version}"),
            "name": name,
            "versionInfo": version,
            "supplier": "NOASSERTION",
            "downloadLocation": f"https://pypi.org/project/{name}/",
            "licenseConcluded": "NOASSERTION",
            "licenseDeclared": lic if len(lic) < 200 else "NOASSERTION",
            "copyrightText": "NOASSERTION",
            "filesAnalyzed": False,
            "externalRefs": [{
                "referenceCategory": "PACKAGE-MANAGER",
                "referenceType": "purl",
                "referenceLocator": f"pkg:pypi/{name}@{version}",
            }],
        })
    return out


def data_packages() -> list[dict]:
    out = []
    for pack in sorted((ROOT / "data" / "packs").glob("*.bin")):
        out.append({
            "SPDXID": spdx_id(f"data-{pack.stem}"),
            "name": pack.name,
            "versionInfo": "NOASSERTION",
            "supplier": SUPPLIER,
            "downloadLocation": "NOASSERTION",
            "checksums": [{"algorithm": "SHA256", "checksumValue": sha256(pack)}],
            "licenseConcluded": "NOASSERTION",
            "licenseDeclared": "NOASSERTION",
            "copyrightText": "NOASSERTION",
            "filesAnalyzed": False,
            "comment": "FEPK coefficient pack built from public NASA/JPL and other sources; see NOTICE and "
                       f"{pack.name}.manifest.json for provenance.",
            "primaryPackagePurpose": "DATA",
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / "dist" / "frenchephemeris.spdx.json")
    ap.add_argument("--version", default="0.1.0")
    ap.add_argument("--no-data", action="store_true", help="omit data packs (no checksums computed)")
    args = ap.parse_args()

    root = {
        "SPDXID": "SPDXRef-FrenchEphemeris",
        "name": "FrenchEphemeris",
        "versionInfo": args.version,
        "supplier": SUPPLIER,
        "downloadLocation": "https://github.com/FreeAstroApi/frenchephemeris",
        "licenseConcluded": PROJECT_LICENSE,
        "licenseDeclared": PROJECT_LICENSE,
        "copyrightText": "Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI",
        "filesAnalyzed": False,
        "primaryPackagePurpose": "LIBRARY",
    }
    packages = [root] + rust_packages() + python_packages() + ([] if args.no_data else data_packages())
    rels = [{"spdxElementId": "SPDXRef-DOCUMENT", "relationshipType": "DESCRIBES",
             "relatedSpdxElement": root["SPDXID"]}]
    rels += [{"spdxElementId": root["SPDXID"], "relationshipType": "DEPENDS_ON",
              "relatedSpdxElement": p["SPDXID"]} for p in packages[1:]]
    doc = {
        "spdxVersion": "SPDX-2.3",
        "dataLicense": "CC0-1.0",
        "SPDXID": "SPDXRef-DOCUMENT",
        "name": f"FrenchEphemeris-{args.version}",
        "documentNamespace": f"https://github.com/FreeAstroApi/frenchephemeris/spdx/{args.version}-{uuid.uuid4()}",
        "creationInfo": {
            "created": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "creators": ["Tool: frenchephemeris-gen_sbom", SUPPLIER],
        },
        "hasExtractedLicensingInfos": [{
            "licenseId": "LicenseRef-FrenchEphemeris-Commercial",
            "name": "FrenchEphemeris commercial licence",
            "extractedText": "See COMMERCIAL-LICENSE.md in the FrenchEphemeris distribution.",
        }],
        "packages": packages,
        "relationships": rels,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=1))
    print(f"wrote {args.out} ({len(packages)} packages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
