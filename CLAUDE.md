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
[`dcm()`](https://epimodel.github.io/EpiModel/reference/dcm.md),
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md), or
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md) to
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
- Long-running tests (`test-net-long.R`, `test-icm-long.R`,
  `test-dcm-long.R`) are guarded by `skip_on_cran()`.
  `R CMD check --as-cran` skips them; drop `--as-cran` locally when
  touching core simulation code and you want that coverage.

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
modules. Since the v2.5 restart refactor, all runtime state lives under
`dat$run` — `dat$run` alone is sufficient to checkpoint and resume a
simulation.

Runtime state (under `dat$run`): - `dat$run$attr` - Node-level
attributes (disease status, demographics, etc.) - `dat$run$el` -
Edgelists per network layer (tergmLite mode) - `dat$run$nw` - Network
objects per network layer (non-tergmLite mode) - `dat$run$net_attr` -
Per-network metadata (size, `lasttoggle`) - `dat$run$nwterms` - Nodal
attributes referenced by the ERGM formula - `dat$run$last_unique_id` -
Counter for newly appended nodes

Top-level (static or aggregated): - `dat$epi` - Epidemic tracking
variables (counts per time step) - `dat$stats` - Network-statistic
tracking (`dat$stats$nwstats`, `dat$stats$transmat`) - `dat$nwparam` -
Network model parameters - `dat$num.nw` - Number of network layers -
`dat$param` / `dat$init` / `dat$control` - User-provided simulation
configuration

Prefer the accessor API (`get_attr`/`set_attr`,
`get_network`/`set_network`, `get_param`, `get_epi`/`set_epi`,
`get_control`, `append_core_attr`/`append_attr`) over direct `dat$run$*`
manipulation. \#977 migrated internal code in this direction, and the
accessors handle validation, `lasttoggle` / `unique_id` bookkeeping, and
tergmLite vs. non-tergmLite dispatch.

### Shared Plotting Infrastructure

`plot.icm` (`R/plot.icm.R`) and `plot_netsim_epi` (`R/plot.netsim.R`)
are sibling plotting paths for the same stochastic-epidemic output —
different data-generating mechanisms, same summary statistics. They
share `draw_means()` / `draw_qnts()` from `R/plot.R` and should stay as
structurally symmetric as possible. Treat behavioral divergences between
them as bugs in the lagging side, not as design. Remaining known
asymmetries are tracked in \#1012.

`plot.dcm` (`R/plot.dcm.R`) is intentionally independent — different
time model, uses `control$timesteps` for x-coordinates directly.

`as.data.frame.netsim` delegates entirely to `as.data.frame.icm`, so
changes to the ICM method automatically apply to netsim.

### Time Handling Differences Across Model Classes

- **DCM**: `control$timesteps` stores the explicit time vector (may be
  non-integer when `dt < 1`). `control$nsteps` is the max time value.
  `control$nruns` is the number of sensitivity runs, derived from the
  longest vector-valued argument in
  [`param.dcm()`](https://epimodel.github.io/EpiModel/reference/param.dcm.md)
  or
  [`init.dcm()`](https://epimodel.github.io/EpiModel/reference/init.dcm.md)
  (initial-condition sensitivity was added in \#408). Default ODE solver
  is `"lsoda"` (changed from `"rk4"` in v2.6.0 for numerical stability).
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

Parallelization uses the `future` / `future.apply` framework (switched
from `foreach`/`doParallel` in v2.6.0).
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
and [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
default to `multisession`; users can override the plan via the
`future.use.plan` control argument.

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
      "ergm", "network", "networkLite", "collections", "future", "future.apply",
      "RColorBrewer", "ape", "lazyeval", "ggplot2", "tibble", "rlang", "dplyr",
      "coda", "Rcpp", "testthat", "knitr", "rmarkdown", "progressr", "tidyr",
      "roxygen2", "lintr", "DT", "shiny"),
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
for applied research. Extension support is **network-class only** —
[`control.icm()`](https://epimodel.github.io/EpiModel/reference/control.icm.md)
no longer accepts custom module functions, `skip.check`, or additional
`...` arguments (removed in \#980 for v2.6.1), so ICMs are restricted to
the built-in SI/SIR/SIS types. Users needing custom ICM-style dynamics
should migrate to
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

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
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
and pass custom module functions (e.g., `infection.FUN = my_infect`).

**Key control.net parameters:** `resimulate.network = TRUE` (for
feedback/demography), `tergmLite = TRUE` (20-50x faster),
`epi.by = "attr"` (stratified stats), `nwstats.formula` (track network
stats).

### Removed and Renamed Arguments

v2.6.1 (#992) completed removal of long-deprecated argument names — they
now hard-error with a pointer to the current name. Use the current names
when writing examples, docs, or tests:

- `b.rate` / `b.rate.g2` → `a.rate` / `a.rate.g2` (`param.dcm` /
  `param.icm` / `param.net`)
- `trans.rate` / `trans.rate.g2` → `inf.prob` / `inf.prob.g2`
- `.m2` parameter suffix → `.g2` (e.g., `i.num.m2` → `i.num.g2`)
- `births.FUN` / `deaths.FUN` → `arrivals.FUN` / `departures.FUN`
  (`control.net`)
- `depend` → `resimulate.network` (`control.net`)

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

`as.data.frame.icm` / `as.data.frame.netsim` return objects that carry
the `epi.data.frame` class alongside `data.frame`; `plot.epi.data.frame`
plots them like `plot.netsim(type = "epi")`. Build them manually with
`as.epi.data.frame(df)` or `df2epi(x)`.

Other exported accessors and helpers worth knowing:

- **Simulation slicing**:
  [`get_sims()`](https://epimodel.github.io/EpiModel/reference/get_sims.md)
  (subset by index),
  [`get_param_set()`](https://epimodel.github.io/EpiModel/reference/get_param_set.md)
  and
  [`get_attr_history()`](https://epimodel.github.io/EpiModel/reference/get_attr_history.md)
  (recorded parameter and attribute trajectories).
- **Transmission / reachability**:
  [`get_discordant_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_discordant_edgelist.md)
  (S–I pairs — generic form, supply the discordance attribute
  explicitly),
  [`get_forward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)
  /
  [`get_backward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)
  over cumulative edgelists (set
  `control.net(save.cumulative.edgelist = TRUE)` to capture them).
- **Restart & truncation**:
  [`truncate_sim()`](https://epimodel.github.io/EpiModel/reference/truncate_sim.md)
  is an S3 generic over `dcm` / `icm` / `netsim` with a `reset.time`
  argument;
  [`make_restart_point()`](https://epimodel.github.io/EpiModel/reference/make_restart_point.md)
  trims a `netsim` to the minimum state needed to resume.
- **Attribute utilities**:
  [`overwrite_attrs()`](https://epimodel.github.io/EpiModel/reference/overwrite_attrs.md)
  applies an `init_attr` data frame at initialization;
  [`get_core_attributes()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
  lists the core attributes and their types (#969).
- **Parameter tables**:
  [`param.net_from_table()`](https://epimodel.github.io/EpiModel/reference/param.net_from_table.md)
  /
  [`param.net_to_table()`](https://epimodel.github.io/EpiModel/reference/param.net_to_table.md)
  convert between list and data.frame parameter representations.
