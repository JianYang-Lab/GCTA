# Build GCTA

## Requirements

Use CMake 3.16 or later and a C++17 compiler. Initialize the plink-ng submodule before configuring:

```sh
git submodule update --init
```

CPM fetches Eigen 3.4.0, Boost 1.86.0, Spectra 1.1.0, zlib 1.3.1, zstd 1.5.6, and the SQLite 3.46.0 amalgamation. Install GSL and Intel MKL separately when required by the target platform. OpenBLAS is used on AArch64.

The regression suite is currently validated on Linux x86-64. Other platforms require separate validation.

## Linux

Configure and build:

```sh
# Example: /opt/intel/oneapi/mkl/latest/lib/cmake/mkl
cmake -DCMAKE_BUILD_TYPE=Release -DMKL_DIR="<your_mkl_cmake_path>" -G Ninja -B build/Release -S .
cmake --build build/Release
```

Run the regression tests:

```sh
ctest --test-dir build/Release --output-on-failure -j 4
```

Install to the configured prefix, or stage the installation with `DESTDIR`:

```sh
cmake --install build/Release
# Optional: installs under /tmp/gcta-install/usr/local by default
DESTDIR=/tmp/gcta-install cmake --install build/Release
```

The installed executable finds bundled libraries through `$ORIGIN/../lib`. MKL and other external libraries must remain available at their configured locations or through the system loader configuration.

## macOS

The build uses Apple Clang and OpenMP. The current configuration expects Homebrew `libomp` under `/opt/homebrew/opt/libomp`; adjust the OpenMP settings for other installations.

```sh
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja -B build/Release -S .
cmake --build build/Release
```

macOS is not covered by the current Linux regression validation.

## Windows

Use the supplied toolchain file and apply the plink-ng patch in `third_party/` when required:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="cmake/win-toolchain.cmake" -DLLVM_ROOT="<your_llvm_path>" -DMKL_DIR="<your_mkl_cmake_path>" -G Ninja -B build/Release -S .
cmake --build build/Release
```

Windows and MinGW64 are not covered by the current Linux regression validation.
