# GCTA
GCTA (Genome-wide Complex Trait Analysis) is a software package, which was initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait but has been extensively extended for many other analyses of data from genome-wide association studies (GWASs). Please see the software website through the link below for more information.

Software website: https://yanglab.westlake.edu.cn/software/gcta/
License: GPLv3 (some parts of the code are released under LGPL as detailed in the files).


## Credits  
Jian Yang developed the original version (before v1.90) of the software (with supports from Peter Visscher, Mike Goddard and Hong Lee) and currently maintains the software.

Zhili Zheng programmed the fastGWA, fastGWA-GLMM and fastGWA-BB modules, rewrote the I/O and GRM modules, improved the GREML and bivariate GREML modules, extended the PCA module, and improved the SBLUP module.  

Zhihong Zhu programmed the mtCOJO and GSMR modules and improved the COJO module.  

Longda Jiang and Hailing Fang developed the ACAT-V module.  

Jian Zeng rewrote the GCTA-HEreg module.  

Andrew Bakshi contributed to the GCTA-fastBAT module.

Angli Xue improved the GSMR module.

Robert Maier improved the GCTA-SBLUP module.

Wujuan Zhong and Judong Shen programmed the fastGWA-GE module. 

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.


## Compilation

#### Requirements

The build requires CMake >= 3.16 and a C++17 compiler. Install the plink-ng submodule before configuring:

```sh
git submodule update --init
```

Most dependencies are fetched by CPM. The current source configuration uses Eigen 3.4.0, Boost 1.86.0, Spectra 1.1.0, zlib 1.3.1, zstd 1.5.6, and the SQLite 3.46.0 amalgamation. GSL and Intel MKL must be installed separately when required by the target platform. OpenBLAS is used on AArch64; the runtime backend can also be selected with `GCTA_NUMERIC_BACKEND`.

The build files cover Linux, macOS, and Windows. The regression suite in this checkout is validated on Linux x86-64; other platforms require separate validation.

#### Configure and build

On Linux:

```sh
# Example: /opt/intel/oneapi/mkl/latest/lib/cmake/mkl
cmake -DCMAKE_BUILD_TYPE=Release -DMKL_DIR="<your_mkl_cmake_path>" -G Ninja -B build/Release -S .
cmake --build build/Release
```

On Windows, use the provided toolchain file and apply the plink-ng patch in `third_party/` when required:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="cmake/win-toolchain.cmake" -DLLVM_ROOT="<your_llvm_path>" -DMKL_DIR="<your_mkl_cmake_path>" -G Ninja -B build/Release -S .
cmake --build build/Release
```

The executable is generated as `build/Release/gcta64`.

#### Test

Regression tests are enabled by default with `GCTA_BUILD_TESTING=ON` and require Python 3:

```sh
ctest --test-dir build/Release --output-on-failure -j 4
```

Tests include independent numerical oracles, format and error paths, and 1/4-thread comparisons. New statistical code should be covered before its legacy implementation is removed.

#### Install

Set `CMAKE_INSTALL_PREFIX` at configure time when a non-default installation location is needed:

```sh
cmake --install build/Release
```

For a staged installation without writing to the configured prefix:

```sh
DESTDIR=/tmp/gcta-install cmake --install build/Release
```

On Linux, the installed executable locates bundled libraries through `$ORIGIN/../lib` and external libraries, including MKL, through their configured link paths. A separate MKL installation is still required; the installation is not a fully self-contained package.

## Development overview

The v2 implementation is organized under `src/{data,analysis,infrastructure}` with mirrored headers in `include/`. Modules expose `registerOption` and `processMain`; unsupported options fall back to `legacy/`. See `docs/development/dev.md` for the module boundaries and migration checklist.

## MCP Server

GCTA provides an MCP (Model Context Protocol) server that exposes all GCTA functionality as tools, allowing AI assistants (Claude, Codex, opencode, etc.) to run GCTA analyses via natural language.

### Quick Start

```bash
python3 -m pip install -r .opencode/gcta-mcp/requirements.txt
python3 .opencode/gcta-mcp/server.py
```

In a source checkout, the server defaults to `build/Release/gcta64` and the
source root. Release-package layouts and custom paths are described in the
MCP README; `GCTA_BINARY_PATH` and `GCTA_WORK_DIR` override the defaults.

### Configuration

Add the following to your AI tool's MCP configuration, replacing `<GCTA_DIR>` with the full path to the GCTA source checkout:

```json
{
  "mcpServers": {
    "gcta": {
      "command": "python3",
      "args": ["<GCTA_DIR>/.opencode/gcta-mcp/server.py"]
    }
  }
}
```

### Available Tools

| Tool | Function |
|------|----------|
| `gcta_run` | Execute GCTA with arbitrary CLI arguments |
| `gcta_help` | Get GCTA documentation |
| `gcta_check` | Check GCTA configuration |
| `gcta_make_grm` | Construct GRM (Genetic Relationship Matrix) |
| `gcta_reml` | REML analysis (estimate heritability) |
| `gcta_bivariate_reml` | Bivariate REML (estimate genetic correlation) |
| `gcta_hereg` | Haseman-Elston regression |
| `gcta_pca` | Principal component analysis |
| `gcta_mlma` | Mixed linear model association |
| `gcta_cojo` | Conditional and joint analysis |
| `gcta_gsmr` | GSMR Mendelian randomization |
| `gcta_mtcojo` | Multi-trait COJO analysis |
| `gcta_fastgwa` | fastGWA genome-wide association |
| `gcta_fastbat` | fastBAT gene-based test |
| `gcta_simu_qt` | Simulate quantitative trait |
| `gcta_simu_cc` | Simulate case-control phenotype |
| `gcta_fst` | Fst population differentiation |
| `gcta_make_bed` | Data management (format conversion, filtering) |
| `gcta_ld_pruning` | LD pruning |
| `gcta_ld_score` | LD score calculation |
| `gcta_acat` | ACAT gene-based test |
| `gcta_list_files` | List files in working directory |
| `gcta_read_file` | Read file content |

See `.opencode/gcta-mcp/README.md` for detailed documentation.

