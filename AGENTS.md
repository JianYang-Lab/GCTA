# AGENTS.md

## Naming and include conventions

- Use `infra` for infrastructure directory names and include paths. Do not add new `infrastructure/` paths.
- Keep source files under `src/` and public headers under the corresponding `include/` module directory.
- Group includes in this order, with exactly one blank line between groups:
  1. The current module's own header.
  2. `analysis/` headers.
  3. `data/` headers.
  4. `infra/` headers.
  5. `legacy/` headers.
  6. Third-party headers such as Eigen, Boost, zlib, and MKL.
  7. C and C++ standard-library headers.
- Within each group, keep includes in lexicographic order when practical. Do not merge unrelated include blocks or move includes across conditional-compilation boundaries.

Example:

```cpp
#include "analysis/expression/ECOJO.h"

#include "analysis/StatLib.h"

#include "data/expression/ExpressionData.h"

#include "infra/Logger.h"
#include "infra/numeric/MatrixOperations.h"

#include "legacy/StatFunc.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
```

## Verification

- Apply include-only changes consistently to project-owned sources and tests; do not reformat `third_party/`.
- After renaming or reorganizing includes, run a clean Release build and the complete regression suite, including `install.runtime`.
- Do not change statistical logic, baselines, or CLI behavior for naming and formatting changes.
