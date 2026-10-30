#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Build fepk-asteroids-v1.bin from docs/ASTEROID_LIST.json.

    python tools/build_asteroids.py fetch     # per body: Horizons SPK -> fit -> ast_<n>.bin, SPK deleted
    python tools/build_asteroids.py bundle    # merge the per-body packs into data/packs/fepk-asteroids-v1.bin
                                              # + docs/ASTEROIDS_PROVENANCE.json

One SPK at a time (deleted right after its fit); the per-body packs live in
data/sources/asteroids/work/ (resumable). Stops if free disk < 2 GB.
Fitting is frenchephemeris/asteroids.py (the same code as the on-demand API).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("fe_asteroids", os.path.join(REPO, "python", "frenchephemeris",
                                                                             "asteroids.py"))
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)

DATA = os.environ.get("FE_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
MAIN = os.path.join(DATA, "packs", "fepk-de441-de442-v1.bin")
WORK = os.path.join(DATA, "sources", "asteroids", "work")
OUT = os.path.join(DATA, "packs", "fepk-asteroids-v1.bin")
LIST = os.path.join(REPO, "docs", "ASTEROID_LIST.json")
LEDGER = os.path.join(REPO, "docs", "ASTEROIDS_PROVENANCE.json")
TOL_KM = 0.5
MIN_FREE = 2 << 30


def free_bytes(path):
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def fetch_all(args):
    os.makedirs(WORK, exist_ok=True)
    bodies = json.load(open(LIST))["bodies"]
    # priority: explicitly requested, then TNOs/centaurs, Swiss-documented, NEOs, classical 1-200
    rank = {"requested": 0, "tno_centaur": 1, "swiss_docs": 2, "neo_myth": 3, "classical": 4}
    # a diverse validation sample first: NEAs, centaurs, TNOs/dwarf planets, Lilith
    first = [1181, 1566, 1862, 7066, 10199, 90377, 136108, 50000, 944, 16, 10]
    bodies.sort(key=lambda b: (b["number"] not in first, first.index(b["number"]) if b["number"] in first else 0,
                               min(rank[r.split(":")[0]] for r in b["reasons"]), b["number"]))
    if args.only:
        want = {int(x) for x in args.only.split(",")}
        bodies = [b for b in bodies if b["number"] in want]
    failures = {}
    for k, b in enumerate(bodies):
        n = b["number"]
        out = os.path.join(WORK, f"ast_{n}.bin")
        if os.path.exists(out):
            continue
        if free_bytes(WORK) < MIN_FREE:
            print(f"STOP: free disk {free_bytes(WORK) / 2**30:.2f} GB < 2 GB", flush=True)
            break
        t = time.time()
        for attempt in range(3):
            try:
                A.fetch(n, cache_dir=WORK, main_pack=MAIN, tol_km=TOL_KM)
                break
            except Exception as e:  # network hiccup or refusal
                err = f"{type(e).__name__}: {e}"
                time.sleep(5)
        else:
            failures[n] = err
            print(f"[{k + 1}/{len(bodies)}] {n} FAILED {err}", flush=True)
            continue
        m = json.load(open(out + ".json"))["bodies"][str(n)]
        print(f"[{k + 1}/{len(bodies)}] {n} {m['spk']['target_line']} {m['spk']['solution_id']} "
              f"{m['fit']['bytes'] / 1e3:.0f} kB max {m['fit']['max_km']:.3f} km  {time.time() - t:.0f}s  "
              f"free {free_bytes(WORK) / 2**30:.2f} GB", flush=True)
    if failures:
        with open(os.path.join(WORK, "failures.json"), "w") as f:
            json.dump(failures, f, indent=1)


def bundle(args):
    bodies = {b["number"]: b for b in json.load(open(LIST))["bodies"]}
    main = A.Fepk(MAIN)
    rows, man_bodies, sources = [], {}, []
    missing = []
    for n, b in sorted(bodies.items()):
        p = os.path.join(WORK, f"ast_{n}.bin")
        if not os.path.exists(p):
            missing.append(n)
            continue
        side = json.load(open(p + ".json"))["bodies"][str(n)]
        pk = A.Fepk(p)
        si = len(sources)
        sources.append({"number": n, "kind": "JPL Horizons small-body SPK (type 21)",
                        "solution_id": side["spk"]["solution_id"], "solution_date": side["request"]["solution_date"],
                        "observations": side["request"]["observations"], "spk_sha256": side["spk"]["sha256"],
                        "spk_target_id": side["spk"]["spk_target_id"], "request": side["request"]["params"],
                        "retrieved_utc": side["request"]["retrieved_utc"],
                        "execution_date_pasadena": side["spk"]["execution_date_pasadena"],
                        "small_body_perturbers": side["spk"]["small_body_perturbers"]})
        for s in pk.series:
            coef = A.np.asarray(pk.mm[s["offset"] // 8:s["offset"] // 8 + s["n_records"] * 3 * s["n_coeff"]]
                                ).reshape(s["n_records"], 3, s["n_coeff"]).copy()
            rows.append(dict(target=s["target"], center=s["center"], start=s["start"], interval=s["interval"],
                             coef=coef, source_index=si))
        man_bodies[n] = {"name": side["spk"]["target_line"].split("{")[0].strip(), "id": A.AST_OFFSET + n,
                         "naif_id": A.naif_id(n), "reasons": b["reasons"], "source_index": si,
                         "solution_id": side["spk"]["solution_id"],
                         "coverage_jd_tdb": side["fit"]["coverage_jd_tdb"], "fit": side["fit"],
                         "jpl_helio_samples_km": side.get("jpl_helio_samples_km", [])}
    manifest = A.pack_manifest("fepk-asteroids", man_bodies, TOL_KM, main)
    manifest["sources"] = sources
    manifest["list"] = "docs/ASTEROID_LIST.json"
    manifest["provenance_ledger"] = "docs/ASTEROIDS_PROVENANCE.json"
    manifest["missing"] = missing
    total = A.write_fepk(OUT, rows, manifest)
    sha = hashlib.sha256(open(OUT, "rb").read()).hexdigest()
    side = dict(manifest)
    side["file"] = {"path": os.path.basename(OUT), "size_bytes": total, "sha256": sha}
    json.dump(side, open(OUT + ".manifest.json", "w"), indent=1)
    ledger = {"generated_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "tool": "tools/build_asteroids.py", "pack": side["file"],
              "clean_room": "All data from JPL Horizons; Swiss is only a black-box oracle in harness/check_asteroids.py",
              "license": manifest["license"], "missing": missing,
              "bodies": [{**{"number": n, "name": v["name"], "solution_id": v["solution_id"],
                             "coverage_jd_tdb": v["coverage_jd_tdb"], "reasons": v["reasons"]},
                          **{k: sources[v["source_index"]][k] for k in ("solution_date", "observations", "spk_sha256",
                                                                         "retrieved_utc")}}
                         for n, v in sorted(man_bodies.items())]}
    json.dump(ledger, open(LEDGER, "w"), indent=1)
    # Rust test fixture: JPL SPK states recorded at fit time (the SPKs themselves are deleted)
    fx = {"generator": "tools/build_asteroids.py bundle", "units": "km, ICRF, heliocentric (wrt 10), JD TDB",
          "tolerance_km": TOL_KM,
          "rows": [[A.naif_id(n), A.AST_OFFSET + n, *row] for n, v in sorted(man_bodies.items())
                   for row in v["jpl_helio_samples_km"]],
          "coverage": {str(A.AST_OFFSET + n): v["coverage_jd_tdb"] for n, v in sorted(man_bodies.items())}}
    fxp = os.path.join(REPO, "crates", "fe-core", "tests", "fixtures", "asteroids", "asteroid_states.json")
    os.makedirs(os.path.dirname(fxp), exist_ok=True)
    json.dump(fx, open(fxp, "w"), indent=0)
    print(f"wrote {OUT}: {total / 1e6:.1f} MB, {len(man_bodies)} bodies, {len(rows)} series; missing {missing}")
    if args.clean:
        shutil.rmtree(WORK)
        print("removed", WORK)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fetch", "bundle"])
    ap.add_argument("--only", help="comma-separated MPC numbers")
    ap.add_argument("--clean", action="store_true", help="bundle: delete the per-body work packs afterwards")
    args = ap.parse_args()
    fetch_all(args) if args.cmd == "fetch" else bundle(args)


if __name__ == "__main__":
    sys.exit(main())
