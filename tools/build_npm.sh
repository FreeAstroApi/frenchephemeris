#!/usr/bin/env bash
# Build the npm package in js/: the WebAssembly engine (js/pkg/) and the bundled mobile data
# set (js/data/). Needs the wasm32-unknown-unknown target and wasm-bindgen-cli 0.2.129.
# Packs come from data/packs/ (build them with tools/build_ext_pack.py, or download a release).
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build -p fe-wasm --release --target wasm32-unknown-unknown
rm -rf js/pkg && wasm-bindgen --target web --out-dir js/pkg --out-name engine \
    target/wasm32-unknown-unknown/release/fe_wasm.wasm
rm -f js/pkg/.gitignore js/pkg/package.json
mkdir -p js/data
for f in fepk-de442-mobile-1800-2400-v1.bin fepk-minor-mobile-1800-2400-v1.bin; do
    cp "data/packs/$f" js/data/
done
cp LICENSE COMMERCIAL-LICENSE.md NOTICE js/
ls -la js/pkg/engine_bg.wasm js/data/*.bin
