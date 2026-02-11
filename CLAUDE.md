# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

EpiModel is an R package for simulating mathematical models of
infectious disease dynamics. It supports three model classes: - **DCM**
(Deterministic Compartmental Models): ODE-based, solved via deSolve -
**ICM** (Individual Contact Models): Stochastic discrete-time
microsimulation - **Network Models**: Stochastic network-based models
using ERGM/STERGM from Statnet

Each model class follows a unified API: `param.*()`, `init.*()`,
`control.*()` for inputs, then
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md),
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md), or
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) to
run simulations. All have `as.data.frame`, `plot`, `summary`, and
`print` S3 methods.

## Build & Development Commands

``` bash
# Install package locally (with compiled C/C++ code)
R CMD INSTALL .

# Generate documentation from Roxygen2 comments (updates NAMESPACE and man/)
Rscript -e 'roxygen2::roxygenise()'

# Run full test suite
Rscript -e 'testthat::test_dir("tests/testthat")'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-dcm.R")'

# Full package check (tests + CRAN compliance)
R CMD check --no-manual .

# Build tarball first, then check (the canonical way)
R CMD build . && R CMD check EpiModel_*.tar.gz

# Lint the package
Rscript -e 'lintr::lint_package()'
```

After modifying any Roxygen2 comments (`#'` blocks) or `@export` tags,
re-run `roxygen2::roxygenise()` to regenerate NAMESPACE and man/ files
before testing.

After modifying C/C++ code in `src/`, the package must be reinstalled
(`R CMD INSTALL .`) before testing.

## Code Style

- 120-character line limit (enforced by lintr)
- Roxygen2 with markdown enabled for documentation
- S3 OOP: classes are `dcm`, `icm`, `netsim`, `netest`, `netdx` with
  method dispatch
- Disabled linters: `object_name_linter`, `cyclocomp_linter`,
  `trailing_blank_lines_linter`, `return_linter`
- Lintr exclusions: `inst/`, `tests/`, `vignettes/`, `R/test.R`,
  `R/RcppExports.R`, `R/dcm.mods.R`

## Architecture

### Modular Simulation Design

Network simulations (`netsim`) use a module-based architecture. Each
simulation time step executes a sequence of modules:

1.  **Initialization** (`net.mod.init.R`) - Set up initial conditions
2.  **Network simulation** (`net.mod.simnet.R`) - Simulate network
    dynamics via ERGM/STERGM
3.  **Infection** (`net.mod.infection.R`) - Disease transmission across
    edges
4.  **Recovery** (`net.mod.recovery.R`) - Disease recovery transitions
5.  **Vital dynamics** (`net.mod.vital.R`) - Arrivals and departures
6.  **Network update** (`net.mod.nwupdate.R`) - Post-step network
    bookkeeping
7.  **Prevalence** (`net.mod.prevalence.R`) - Track disease prevalence
8.  **Summary** (`net.mod.summary.R`) - Record summary statistics

Users can replace or extend any module via
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
to build custom models. ICM follows a similar but simpler module pattern
(`icm.mod.*.R`).

### Shared Plotting Infrastructure

`plot.icm` (`R/plot.icm.R`) and `plot_netsim_epi` (`R/plot.netsim.R`)
share helper functions `draw_means()` and `draw_qnts()` defined in
`R/plot.R`. Changes to plotting behavior (e.g., axis offsets, coordinate
handling) must be coordinated across all three files. `plot.dcm`
(`R/plot.dcm.R`) is independent and uses `control$timesteps` directly
for x-coordinates.

`as.data.frame.netsim` delegates entirely to `as.data.frame.icm`, so
changes to the ICM method automatically apply to netsim. The
`as.epi.data.frame` validator in the same file checks structural
consistency of the output.

### Time Handling Differences Across Model Classes

DCM and ICM/netsim handle time differently in their output objects:

- **DCM**: `control$timesteps` stores the explicit time vector (may be
  non-integer when `dt < 1`). `control$nsteps` is the max time value.
  `control$nruns` is the number of parameter sensitivity runs. Plot and
  data frame methods use `timesteps` directly.
- **ICM/netsim**: Time is implicit via row indices. `control$nsteps` is
  both the number of rows and the max time step. `control$nsims` is the
  number of stochastic simulations. `control$start` provides a time
  offset (default 1). `as.data.frame.icm` builds the time column from
  `start` and `nsteps`.

When modifying functions that touch time or epi structure, verify
behavior across all three classes and check that `as.data.frame`,
`plot`, and `summary` methods still work correctly.

### Key Source File Groups

- **`R/net.inputs.R`** (largest file): Validates and processes all
  network model inputs (`param.net`, `init.net`, `control.net`)
- **`R/net.fn.utils.R`**, **`R/net.fn.accessor.R`**: Core utilities and
  attribute accessors for the internal `dat` simulation object
- **`R/get.R`**: Public API for extracting simulation results
  (`get_param`, `get_attr`, `get_epi`, etc.)
- **`R/edgelists.R`**: Edgelist manipulation and cumulative edgelist
  tracking
- **`R/dcm.mods.R`**: Built-in ODE systems for standard SI/SIR/SIS
  compartmental models

### The `dat` Object

Network simulations pass a central `dat` list object through all
modules. Key sub-elements: - `dat$attr` - Node-level attributes (disease
status, demographics, etc.) - `dat$epi` - Epidemic tracking variables
(counts per time step) - `dat$nwparam` - Network model parameters -
`dat$param` / `dat$init` / `dat$control` - Simulation configuration -
`dat$run` - Runtime state (current time step, simulation metadata) -
`dat$el` - Edgelists for each network layer (when using tergmLite
mode) - `dat$net` - Network objects (when not using tergmLite mode)

### Compiled Code (src/)

C/C++ source provides performance-critical functionality: -
`changestats.users.c` - Custom ERGM terms (`fuzzynodematch`,
`absdiffby`, `absdiffnodemix`) - `updatenw.cpp` - Rcpp-based network
update operations - `RcppExports.cpp` - Auto-generated Rcpp bindings (do
not edit manually)

### Multi-Network Support

Models can use multiple network layers (e.g., sexual + needle-sharing
networks). Network parameters, edgelists, and diagnostics are indexed by
network number. See `test-multinets.R` for examples.

## Dependencies

The package heavily depends on the Statnet ecosystem: `ergm`, `tergm`,
`network`, `networkDynamic`, `networkLite`, `statnet.common`.
Understanding ERGM model specification is important for working with
network model code.
