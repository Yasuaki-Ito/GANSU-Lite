#!/usr/bin/env bash
# Build both standard and SIMD WASM binaries.
# Usage: cd wasm-eri && ./build-wasm.sh
#
# Prerequisites:
#   cargo install wasm-pack
#   rustup target add wasm32-unknown-unknown

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

OUT_DIR="../public/wasm"
mkdir -p "$OUT_DIR"

echo "=== Building standard WASM binary ==="
wasm-pack build --target web --release --out-dir pkg
cp pkg/wasm_eri_bg.wasm "$OUT_DIR/wasm_eri_bg.wasm"
cp pkg/wasm_eri.js "$OUT_DIR/wasm_eri.js"
cp pkg/wasm_eri.d.ts "$OUT_DIR/wasm_eri.d.ts"
cp pkg/wasm_eri_bg.wasm.d.ts "$OUT_DIR/wasm_eri_bg.wasm.d.ts"
echo "  -> $(wc -c < "$OUT_DIR/wasm_eri_bg.wasm" | tr -d ' ') bytes"

echo ""
echo "=== Building SIMD WASM binary ==="
RUSTFLAGS="-C target-feature=+simd128" wasm-pack build --target web --release --out-dir pkg-simd
cp pkg-simd/wasm_eri_bg.wasm "$OUT_DIR/wasm_eri_simd_bg.wasm"
echo "  -> $(wc -c < "$OUT_DIR/wasm_eri_simd_bg.wasm" | tr -d ' ') bytes"

echo ""
echo "Done. Files in $OUT_DIR:"
ls -lh "$OUT_DIR"/wasm_eri*
