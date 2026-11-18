# Development guide

## Source layout

The v2 implementation uses mirrored source and header directories:

```text
src/data/                  include/data/
src/analysis/              include/analysis/
src/infrastructure/        include/infrastructure/
legacy/                    v1 implementation and fallback option parser
```

Current module groups include:

- `data/expression` and `analysis/expression`: expression input, ERM, and eCOJO.
- `data/summary` and `analysis/mr`: summary statistics, reference genotypes, GSMR, and mtCOJO.
- `infrastructure/io`: relationship-matrix output and other shared I/O.
- `infrastructure/numeric`: matrix operations and numerical backends.

The legacy PC-adjustment workflow still uses shared legacy summary support. `StatFunc`, string utilities, and compressed text helpers have not been fully consolidated and must not be replaced merely for architectural symmetry.

## Option dispatch

`src/main.cpp` registers modules in a fixed order. A module implements:

- `registerOption(options)`: consumes its options and returns whether it claims the main analysis slot.
- `processMain()`: constructs its data objects and runs the analysis.

Only one main module may run. A main module is selected only when every command-line flag is listed in `supported_flagsV2`; otherwise the command falls back to `legacy/option.cpp`. When two modules share an input flag, only the module that actually runs should consume it. Preserve existing CLI flags and fallback behavior when adding a module.

Use `Pheno`, `Marker`, and `Geno` for v2 genotype input. `Geno::loopDouble()` accepts extracted-marker indices, while marker-management functions may use raw BIM indices; do not mix these index spaces. For ordinary BED input, centered and standardized `GenoBufItem::geno` values use allele-frequency scaling, matching the legacy GSMR convention.

## Build rules

CMake discovers sources recursively with `CONFIGURE_DEPENDS` and creates a library name from each source basename. Keep basenames unique across directories and declare inter-library dependencies explicitly. A new directory does not by itself prevent target-name collisions.

The Linux install rule sets `gcta64`'s install RPATH directly on the target. It contains `$ORIGIN/../lib` and external link directories, including the configured MKL path.

## Regression requirements

Run tests after each change:

```sh
cmake --build build/Release -j 4
ctest --test-dir build/Release --output-on-failure -j 4
```

For a clean validation build, configure a separate Release directory and run the same tests plus `install.runtime`.

When migrating a workflow:

1. Add regression tests against the legacy implementation first. Use independent numerical oracles where practical; exit-code-only tests are insufficient for statistical code.
2. Cover primary outputs, errors, exit status, input subsets, and 1/4-thread agreement. Prioritize main workflows and high-risk branches rather than every flag combination.
3. Move the implementation without changing formulas or output formats, then compare v1/v2 results while both paths exist.
4. Confirm that tests enter the intended v2 path rather than silently falling back to v1.
5. Delete the legacy implementation in a separate commit only after compatibility is established. Do not rewrite baselines to hide differences.

Deterministic binary outputs should match byte-for-byte. Text results should be compared at their written precision, with documented handling for known nondeterministic fields.
