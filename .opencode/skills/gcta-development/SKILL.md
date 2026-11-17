---
name: gcta-development
description: Use when working on the GCTA (Genome-wide Complex Trait Analysis) C++ codebase — adding/modifying analysis modules, command-line options, genotype/phenotype/GRM processing, or debugging build issues. Covers the v2 module system in src/, the legacy v1 gcta class in legacy/, the CMake/CPM build, the Logger singleton, OpenMP threading, and the option-registration dispatch pattern.
---

# GCTA Development Skill

## What is GCTA

GCTA (Genome-wide Complex Trait Analysis) is a C++ command-line tool for
genome-wide association study (GWAS) analyses: estimating genetic variance
components (GREML), constructing genetic relationship matrices (GRM),
mixed-model association (fastGWA / MLMA), conditional/joint analysis (COJO),
Mendelian randomization (GSMR), and more. Output is a `.log` file plus
analysis-specific binary/text outputs.

- Upstream docs: https://yanglab.westlake.edu.cn/software/gcta/
- Entry point: `src/main.cpp` -> `main()`.

## Repository layout

```
src/main.cpp       v2 option parsing and dispatch entry point.
src/data/          v2 phenotype, marker, genotype, and covariate processing.
src/analysis/      v2 GRM, fastGWA, LD, ACAT, and statistical routines.
src/infrastructure/ v2 logging, I/O, numerical backend, threading, and utilities.
include/           v2 headers in matching data/analysis/infrastructure directories.
legacy/            v1 implementation: the legacy `gcta` class + option parser.
third_party/       plink-ng (git submodule), Pgenlib, zstd wrapper.
cmake/          CPM.cmake, toolchain files.
resources/      desktop entry + icons.
scripts/        HPC and macOS helper scripts.
docs/           (mostly empty) Sphinx-style manual.
```

The executable target is `gcta64`, built from `src/main.cpp` + per-module
libraries auto-globbed recursively from `src/`, plus the `mainV1` library from
`legacy/CMakeLists.txt`.

## Two-layer architecture (IMPORTANT)

GCTA has TWO coexisting code generations. Understand which one you are touching:

1. **v2 ("mock layer")** — `src/{data,analysis,infrastructure}/` + matching `include/` directories.
   - Modular: each module is a class (`Pheno`, `Marker`, `Geno`, `Covar`,
     `GRM`, `FastFAM`, `LD`) with static `registerOption(options)` and
     `processMain()` methods.
   - `src/main.cpp` builds an `options` map from argv, then calls each
     module's `registerOption`. A module returns 1 to claim the "main" slot;
     `processMains[mains[0]]()` then dispatches. Only ONE main module per run.
   - v2 uses bitwise genotype processing (uint64_t blocks), `AsyncBuffer`
     for I/O, and the global `Logger` singleton (`LOGGER`).
   - The `supported_flagsV2` list in `src/main.cpp` enumerates which flags
     are handled by v2. Anything else falls through to v1.

2. **v1 (legacy)** — `legacy/*.cpp` + `legacy/gcta.h`.
   - The monolithic `gcta` class (~700 lines in the header) holds all
     genotype/phenotype/GRM state and methods like `read_bimfile`,
     `make_grm`, `mlm_assoc`, `cojo`, `gsmr`, etc.
   - `legacy/option.cpp` / `option.h` parse the command line via the legacy
     `option(argc, argv)` function, called as a fallback when v2 does not
     claim the main slot.
   - v1 uses Eigen dense/sparse matrices and the `CommFunc`/`StatFunc`/
     `StrFunc` utility libraries.

**Dispatch flow:**
```
main() -> parse argv into options map
       -> for each v2 module: registerOption(options)  // returns 1 if it claims main
       -> if any module claimed main AND all flags are in supported_flagsV2:
              processMains[mains[0]]()   // v2 path
          else:
              option(argc, argv)         // v1 fallback path
```

## Module registration pattern (v2)

To add a new v2 module or understand an existing one:

1. **Header** (`include/analysis/Foo.h`): declare a class with:
   ```cpp
   static int registerOption(map<string, vector<string>>& options_in);
   static void processMain();
   ```
2. **Source** (`src/analysis/Foo.cpp`): implement `registerOption` to scan `options_in`
   for this module's flags and return `1` if a main action is requested (or
   `0` if the module only contributes data). Implement `processMain` as the
   entry point executed by `main()`.
3. **Wire into `src/main.cpp`:**
   - Add the class name to `module_names`.
   - Add `Foo::registerOption` to `registers`.
   - Add `Foo::processMain` to `processMains`.
   - Add any new flags to `supported_flagsV2`.

The `module_names` vector order matters — it defines registration order.
Modules that produce shared data (Pheno, Marker, Geno, Covar) are registered
before consumers (GRM, FastFAM, LD).

## Logger (global singleton)

`include/infrastructure/Logger.h` + `src/infrastructure/Logger.cpp`. Accessed via the `LOGGER` macro.
- `LOGGER.i(level, msg)` / `LOGGER.w(...)` / `LOGGER.e(level, msg)` —
  `e` logs an error and EXITS the program.
- `LOGGER.m(level, msg)` — terminal-only message (not written to log file).
- `LOGGER.l(level, msg)` — log-file-only.
- `LOGGER.p(level, msg)` — progress (single-line, terminal only).
- `LOGGER << value` — stream-style output.
- `LOGGER.open(filename)` opens the log file; called from `main()` after
  option parsing using `<out>.log`.
- `LOGGER.ts(marker)` / `LOGGER.tp(marker)` — timestamp start/elapsed.

## Option parsing conventions

- All CLI flags use `--flag` prefix; values are positional following the flag.
- `src/main.cpp` parses argv into `map<string, vector<string>> options`,
  where each key is a flag and the value is the list of arguments until the
  next `--flag`.
- `--out` is required (defaults to `"output"` with a warning). The log file
  is `<out>.log`.
- `--thread-num` / `--threads` / `OMP_NUM_THREADS` env var control OpenMP
  thread count.
- Modules store parsed options in static `map<string, string> options` and
  `map<string, double> options_d` members.

## Threading

- OpenMP is used throughout (`#pragma omp parallel`, `omp_set_num_threads`).
- `src/infrastructure/ThreadPool.cpp` provides a custom thread pool for task-based
  parallelism.
- Thread count is read from `OMP_NUM_THREADS` env var or `--thread-num` /
  `--threads` CLI flags.
- On macOS, OpenMP requires libomp from Homebrew (`/opt/homebrew/opt/libomp`).

## Build system

- **CMake** >= 3.16, C++17 standard.
- **CPM** (CMake Package Manager) fetches dependencies:
  - zlib 1.3.1, Eigen 3.4.0, Spectra 1.1.0, zstd 1.5.6, Boost 1.86.0,
    sqlite3 3460000.
- **External dependencies** (must be installed manually):
  - Intel MKL (Linux/Windows) or Accelerate/OpenBLAS (macOS).
  - GSL (GNU Scientific Library) — Linux only.
  - plink-ng submodule: `git submodule update --init`.
- **Build commands:**
  ```sh
  # Linux
  cmake -DCMAKE_BUILD_TYPE=Release -DMKL_DIR="<mkl_cmake_path>" -G Ninja -B build/Release -S .
  cmake --build build/Release

  # macOS (Apple Silicon)
  cmake -DCMAKE_BUILD_TYPE=Release -G Ninja -B build/Release -S .
  cmake --build build/Release
  ```
- Each v2 module `.cpp` is compiled into its own static library, then linked into
  the `gcta64` executable along with `mainV1`.
- `compile_commands.json` is generated and copied to the source root for
  IDE/clangd integration.

## Key data types

- **`GenoBuf`** (`include/data/Geno.h`): buffered genotype block with allele
  frequencies, missing-data bits, and centered/standardized genotype vectors.
- **`GenoBufItem`**: per-marker extracted genotype result (geno vector,
  missing bits, af, mean, sd, info, nValidN, nValidAllele).
- **`MarkerInfo`** / **`MarkerParam`** (`include/data/Marker.h`): SNP metadata
  (chr, name, pd, alleles, A_rev) and BGEN parsing parameters.
- **`PhenoMask`** (`include/data/Pheno.h`): bitmask + shift for phenotype flags.
- **Eigen types** (v1, `legacy/gcta.h`): `eigenMatrix`, `eigenVector`,
  `eigenSparseMat`, `eigenDiagMat` — aliased to double or float based on
  `SINGLE_PRECISION` define.

## Genotype file formats

v2 `Geno` class supports:
- **PLINK 1 binary** (`.bed`/`.bim`/`.fam`): bitwise processing, 64-bit
  block extraction.
- **BGEN** (`.bgen` + `.sample`): requires `--sample` option when using
  `--bgen` or `--mbgen`.
- **PLINK 2** (`.pgen`/`.pvar`/`.psam`): via the plink-ng pgen library
  in `third_party/`.

Multi-file input: `--bfile` (single) / `--mbfile` (multi); `--pfile` /
`--mpfile`; `--bgen` / `--mbgen`.

## ACAT module

The ACAT (Aggregated Cauchy Association Test) module has special handling in
`src/main.cpp` (lines ~107-148): if `--acat` is passed with `--gene-list`
and `--snp-list`, it calls `acat_func()` directly and exits before the normal
option dispatch. The ACAT implementation is in `src/analysis/acat.c` /
`include/analysis/acat.hpp`.

## Common development tasks

### Adding a new CLI option to an existing v2 module

1. Add the flag string to `supported_flagsV2` in `src/main.cpp`.
2. In the module's `registerOption`, check for the flag in `options_in` and
   store its value(s) in the module's static `options` map.
3. Use the parsed value in `processMain` or the module's processing methods.

### Adding a new v2 module

1. Create `include/analysis/Foo.h` and `src/analysis/Foo.cpp` following the registration
   pattern above.
2. Add the module to `module_names`, `registers`, and `processMains` in
   `src/main.cpp`.
3. Add any new flags to `supported_flagsV2`.
4. The new `.cpp` will be auto-globbed into a library by CMake.

### Modifying v1 (legacy) functionality

1. Edit `legacy/option.cpp` / `option.h` for CLI parsing changes.
2. Edit the relevant `legacy/*.cpp` file (e.g., `est_hsq.cpp` for GREML,
   `grm.cpp` for GRM construction, `gsmr.cpp` for GSMR, etc.).
3. The `gcta` class in `legacy/gcta.h` is the central data structure.

### Debugging

- Use `LOGGER.d(level, msg)` for debug messages (only visible in debug builds).
- `#ifndef NDEBUG` blocks in `include/analysis/GRM.h` show test-only file handles.
- Build in Debug mode: `cmake -DCMAKE_BUILD_TYPE=Debug ...`.
- The `--verbose` flag shows peak memory usage on Linux.

## Conventions

- **No comments** unless explicitly necessary (codebase style).
- **Error handling**: use `LOGGER.e(0, "message")` to log and exit; use
  `LOGGER.w(0, "message")` for warnings.
- **Memory**: use `posix_mem_free` for aligned allocations (see `GRM`
  destructor); use `mem.hpp` utilities for memory tracking.
- **File I/O**: use the `Logger` for log output; use standard `FILE*` or
  `std::ofstream` for data output.
- **Naming**: classes are PascalCase (`Geno`, `GRM`, `FastFAM`); methods are
  snake_case or camelCase; static members prefixed with `s` or in static
  maps; member variables often prefixed with `b` (bool) or `n` (count).
- **Includes**: v2 headers live under `include/{data,analysis,infrastructure}/`,
  v1 headers in `legacy/`. Use qualified project includes such as
  `"data/Geno.h"`, and `<Header>` for system/third-party headers.
- **Third-party headers**: PLINK2 pgen library at
  `third_party/plink-ng/2.0/include`; Eigen at `eigen-mirror/eigen#3.4.0`.
