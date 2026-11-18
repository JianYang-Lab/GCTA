#!/bin/sh
# Build and test GCTA on Linux x86-64.
set -eu

BUILD_DIR="${BUILD_DIR:-build/Release}"
MKL_DIR="${MKL_DIR:-/opt/intel/oneapi/mkl/latest/lib/cmake/mkl}"

git submodule update --init
cmake -DCMAKE_BUILD_TYPE=Release -DMKL_DIR="$MKL_DIR" -G Ninja -B "$BUILD_DIR" -S .
cmake --build "$BUILD_DIR"
ctest --test-dir "$BUILD_DIR" --output-on-failure -j 4

# To install to the configured prefix, run: cmake --install "$BUILD_DIR"
