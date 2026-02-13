# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

EpiModel is an R package for simulating mathematical models of
infectious disease dynamics. Its distinguishing capability is stochastic
network epidemic modeling grounded in Exponential-family Random Graph
Models (ERGMs). It supports three model classes:

1.  **Deterministic Compartmental Models (DCMs)** – ODE-based (using
    `deSolve`). Population-level, continuous-time, no stochasticity.
2.  **Stochastic Individual Contact Models (ICMs)** – discrete-time
    agent-based with random mixing. No persistent network structure.
3.  **Stochastic Network Models** – the primary focus. Uses ERGMs/TERGMs
    to represent dynamic contact networks where partnerships form,
    persist, and dissolve according to statistically estimated rules.

Each model class follows a unified API: `param.*()`, `init.*()`,
`control.*()` for inputs, then
[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md),
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md), or
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) to
run simulations. All have `as.data.frame`, `plot`, `summary`, and
`print` S3 methods.

Built-in compartmental structures: **SI**, **SIR**, and **SIS**, in
one-group or two-group configurations, with or without demography. The
extension API allows arbitrary disease states.

## Build & Development Commands

``` bash
# Install package locally (with compiled C/C++ code)
R CMD INSTALL .

# Generate documentation from Roxygen2 comments (updates NAMESPACE and man/)
Rscript -e 'roxygen2::roxygenise()'

# Full package check -- the primary way to test (tests + examples + CRAN compliance)
# This runs all unit tests in tests/ AND all examples in man pages
R CMD build . && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual EpiModel_*.tar.gz

# Run a single test file (must load the package first for namespace access)
Rscript -e 'library(EpiModel); testthat::test_file("tests/testthat/test-dcm.R")'

# Lint the package
Rscript -e 'lintr::lint_package()'
```

After modifying any Roxygen2 comments (`#'` blocks) or `@export` tags,
re-run `roxygen2::roxygenise()` to regenerate NAMESPACE and man/ files
before testing.

After modifying C/C++ code in `src/`, the package must be reinstalled
(`R CMD INSTALL .`) before testing.

### Testing Notes

- **Use `R CMD check` as the canonical test method.** It runs all unit
  tests via
  [`testthat::test_check()`](https://testthat.r-lib.org/reference/test_package.html)
  (which properly loads the package namespace) AND runs all
  documentation examples. There is no need to separately run
  `test_check` or `test_package` – `R CMD check` already covers both.
- **Do NOT use `testthat::test_dir("tests/testthat")`** – it does not
  load the package namespace, so internal/unexported functions (e.g.,
  `tedgelist_to_toggles`, `name_saveout_elts`) and re-exported imports
  (e.g., `networkLite`) will not be found, causing spurious test
  failures.
- Use `_R_CHECK_FORCE_SUGGESTS_=false` to skip optional Suggests
  packages that may not be installed (e.g., `ergm.ego`, `ndtv`, `egor`).

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

### Shared Plotting Infrastructure

`plot.icm` (`R/plot.icm.R`) and `plot_netsim_epi` (`R/plot.netsim.R`)
share helper functions `draw_means()` and `draw_qnts()` defined in
`R/plot.R`. Changes to plotting behavior must be coordinated across all
three files. `plot.dcm` (`R/plot.dcm.R`) is independent and uses
`control$timesteps` directly for x-coordinates.

`as.data.frame.netsim` delegates entirely to `as.data.frame.icm`, so
changes to the ICM method automatically apply to netsim.

### Time Handling Differences Across Model Classes

- **DCM**: `control$timesteps` stores the explicit time vector (may be
  non-integer when `dt < 1`). `control$nsteps` is the max time value.
  `control$nruns` is the number of parameter sensitivity runs.
- **ICM/netsim**: Time is implicit via row indices. `control$nsteps` is
  both the number of rows and the max time step. `control$nsims` is the
  number of stochastic simulations. `control$start` provides a time
  offset (default 1).

When modifying functions that touch time or epi structure, verify
behavior across all three classes.

### Compiled Code (src/)

- `changestats.users.c` - Custom ERGM terms (`fuzzynodematch`,
  `absdiffby`, `absdiffnodemix`)
- `updatenw.cpp` - Rcpp-based network update operations
- `RcppExports.cpp` - Auto-generated Rcpp bindings (do not edit
  manually)

### Multi-Network Support

Models can use multiple network layers (e.g., sexual + needle-sharing
networks). Network parameters, edgelists, and diagnostics are indexed by
network number. See `test-multinets.R` for examples.

## Dependencies

The package heavily depends on the Statnet ecosystem: `ergm`, `tergm`,
`network`, `networkDynamic`, `networkLite`, `statnet.common`.
Understanding ERGM model specification is important for working with
network model code.

## Environment Setup for Testing

R is not pre-installed in the Claude Code environment but can be
installed. When asked to run a full package test:

1.  **Install R 4.5+ from the CRAN repository** (the default Ubuntu repo
    has an older version):

    ``` bash
    chmod 1777 /tmp  # fix temp dir permissions if needed
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
      gpg --dearmor -o /usr/share/keyrings/r-project.gpg
    echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] \
      https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
      > /etc/apt/sources.list.d/r-project.list
    apt-get update && apt-get install -y r-base r-base-dev
    ```

2.  **Install system libraries** needed by R package dependencies:

    ``` bash
    apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev \
      libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
      libpng-dev libtiff5-dev libjpeg-dev pandoc
    ```

3.  **Install R package dependencies**:

    ``` r
    install.packages(c("deSolve", "networkDynamic", "tergm", "statnet.common",
      "ergm", "network", "networkLite", "collections", "doParallel", "foreach",
      "RColorBrewer", "ape", "lazyeval", "ggplot2", "tibble", "rlang", "dplyr",
      "coda", "Rcpp", "testthat", "knitr", "rmarkdown", "progressr", "tidyr",
      "roxygen2", "lintr", "xml2", "DT", "shiny"),
      repos = "https://cloud.r-project.org", Ncpus = 4)
    ```

4.  **Install and check** the package:

    ``` bash
    R CMD INSTALL .
    R CMD build . && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual EpiModel_*.tar.gz
    ```

For minor changes (plotting, data frame methods, etc.) where a full
environment setup is not warranted, push and rely on CI for integration
tests.

------------------------------------------------------------------------

## Domain Context

The sections below provide broader context about EpiModel’s
epidemiological modeling concepts. This helps when working on features,
writing documentation, or understanding the purpose behind code
patterns.

### Why Network Models Matter

Traditional epidemic models assume random (“mass-action”) mixing. Real
contact patterns are structured – people have a fixed number of
partners, partnerships have durations, there is clustering, assortative
mixing, and concurrency (overlapping partnerships).

Network models capture these features: degree distributions, relational
durations, concurrency, assortative mixing by attributes, clustering,
and bidirectional feedback between network structure and epidemic
dynamics.

### The Network Model Pipeline

The core workflow for any network epidemic model:

``` r
# 1. Initialize network
nw <- network_initialize(n = 500)
nw <- set_vertex_attribute(nw, "risk", rep(0:1, each = 250))

# 2. Parameterize formation/dissolution (target edges = mean_degree * N / 2)
formation <- ~edges + concurrent + degrange(from = 4)
target.stats <- c(175, 110, 0)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)

# 3. Estimate network model
est <- netest(nw, formation, target.stats, coef.diss)

# 4. Diagnose fit
dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 4)

# 5. Simulate epidemic
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.1)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsims = 5, nsteps = 500)
sim <- netsim(est, param, init, control)
```

### ERGM Terms Reference

**Dyad-independent** (faster estimation): `edges`, `nodefactor`,
`nodematch` (homophily, use `diff=TRUE` for per-level), `nodemix`,
`nodecov`, `absdiff`.

**Dyad-dependent** (require MCMC): `degree(k)`, `concurrent` (alias for
`mindegree(2)`), `degrange(from = k)`, `triangles`, `gwesp`.

**Dissolution models** use
[`offset()`](https://rdrr.io/r/stats/offset.html) terms only:
`~offset(edges)` for homogeneous duration, add
`offset(nodematch("attr"))` for differential duration.

### Key Formulas

- **Transmission probability**:
  `finalProb = 1 - (1 - transProb)^actRate`
- **Dissolution adjustment for mortality**:
  `dissolution_coefs(~offset(edges), duration = 40, d.rate = 0.001)`
- **Estimation methods**: Approximation (`edapprox = TRUE`, default,
  best for durations \> 50) vs. direct STERGM (`edapprox = FALSE`, for
  short durations ~10-25).

### The Extension API

The extension API allows custom models with arbitrary disease states,
demographics, interventions, and feedback. This is the most used feature
for applied research.

**Module function signature:** `function(dat, at) { ...; return(dat) }`
where `dat` is the master state object and `at` is the current time
step.

**Accessor functions:** `get_attr`/`set_attr` (nodal attributes),
`get_param` (parameters), `get_epi`/`set_epi` (epidemic stats),
`append_core_attr`/`append_attr` (new nodes), `discord_edgelist` (S-I
pairs), `set_transmat` (transmission events).

**Required nodal attributes:** `active`, `status`, `infTime`,
`entrTime`, `exitTime`.

**Departure rule:** Set `active[ids] <- 0` and `exitTime[ids] <- at`.

**Arrival rule:** Call `append_core_attr(dat, at, n.new)` first, then
`append_attr` for each custom attribute.

**statusTime guard pattern:** Prevent multiple transitions per time step
with `idsElig <- which(active == 1 & status == "e" & statusTime < at)`.

**Extension model setup:** Set `type = NULL` in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
and pass custom module functions (e.g., `infection.FUN = my_infect`).

**Key control.net parameters:** `resimulate.network = TRUE` (for
feedback/demography), `tergmLite = TRUE` (20-50x faster),
`epi.by = "attr"` (stratified stats), `nwstats.formula` (track network
stats).

### Feedback Mechanisms

Network models can incorporate bidirectional feedback: demography
(births/deaths reshape network), serosorting (using `status` as ERGM
term), behavioral interventions (reduced act rates for diagnosed
individuals), and built-in interventions (`inter.eff`, `inter.start`).

When `status` is used in the ERGM formation formula, it must be
initialized directly on the network object (not via `init.net`). Network
statistics are NOT preserved at target levels – log-odds conditional on
other terms are preserved, so statistics shift as prevalence changes.

### Output Analysis

``` r
as.data.frame(sim)                            # per-simulation raw data
as.data.frame(sim, out = "mean")              # means across simulations
get_network(sim, sim = 1)                     # extract network object
get_transmat(sim, sim = 1)                    # transmission matrix
mutate_epi(sim, ir = (si.flow / s.num) * 100) # derived statistics
plot(sim)                                      # compartment sizes over time
```
