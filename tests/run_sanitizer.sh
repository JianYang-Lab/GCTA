#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-${ROOT}/build/Sanitizer}"
JOBS="${JOBS:-4}"
CMAKE_ARGS=()
if [[ -n "${MKL_DIR:-}" ]]; then
    CMAKE_ARGS+=("-DMKL_DIR=${MKL_DIR}")
fi

cmake -S "${ROOT}" -B "${BUILD}" -G Ninja \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DGCTA_SANITIZE=address,undefined \
    "${CMAKE_ARGS[@]}"
cmake --build "${BUILD}" -j "${JOBS}"

ASAN_LIBRARY="$(gcc -print-file-name=libasan.so.8)"
if [[ ! -e "${ASAN_LIBRARY}" ]]; then
    ASAN_LIBRARY="$(gcc -print-file-name=libasan.so)"
fi
export LD_PRELOAD="${ASAN_LIBRARY}"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"
ctest --test-dir "${BUILD}" --output-on-failure -j "${JOBS}" -L integration

unset LD_PRELOAD
LSAN_WORK="${BUILD}/output/lsan"
mkdir -p "${LSAN_WORK}"
LSAN_WRAPPER="${BUILD}/gcta64-lsan"
cat > "${LSAN_WRAPPER}" <<EOF
#!/bin/sh
LD_PRELOAD="${ASAN_LIBRARY}" ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" exec "${BUILD}/gcta64" "\$@"
EOF
chmod +x "${LSAN_WRAPPER}"

echo "Running LSan on workflows that do not load genotype data; the known PgenReader leak is not suppressed."
python3 "${ROOT}/tests/regression/test_legacy_commands.py" \
    --gcta "${LSAN_WRAPPER}" --data "${ROOT}/tests/data" \
    --work "${LSAN_WORK}/ld-file" --scenario ld-file
python3 "${ROOT}/tests/regression/test_acat.py" \
    --gcta "${LSAN_WRAPPER}" --data "${ROOT}/tests/data" \
    --work "${LSAN_WORK}/acat" --scenario output
