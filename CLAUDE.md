# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EpiModel is an R package for simulating mathematical models of infectious disease dynamics. It supports three model classes:
- **DCM** (Deterministic Compartmental Models): ODE-based, solved via deSolve
- **ICM** (Individual Contact Models): Stochastic discrete-time microsimulation
- **Network Models**: Stochastic network-based models using ERGM/STERGM from Statnet

Each model class follows a unified API: `param.*()`, `init.*()`, `control.*()` for inputs, then `dcm()`, `icm()`, or `netsim()` to run simulations. All have `as.data.frame`, `plot`, `summary`, and `print` S3 methods.

## Build & Development Commands

```bash
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

After modifying any Roxygen2 comments (`#'` blocks) or `@export` tags, re-run `roxygen2::roxygenise()` to regenerate NAMESPACE and man/ files before testing.

After modifying C/C++ code in `src/`, the package must be reinstalled (`R CMD INSTALL .`) before testing.

## Code Style

- 120-character line limit (enforced by lintr)
- Roxygen2 with markdown enabled for documentation
- S3 OOP: classes are `dcm`, `icm`, `netsim`, `netest`, `netdx` with method dispatch
- Disabled linters: `object_name_linter`, `cyclocomp_linter`, `trailing_blank_lines_linter`, `return_linter`
- Lintr exclusions: `inst/`, `tests/`, `vignettes/`, `R/test.R`, `R/RcppExports.R`, `R/dcm.mods.R`

## What is EpiModel?

EpiModel is an R package for simulating mathematical models of infectious disease dynamics. Its distinguishing capability is stochastic network epidemic modeling grounded in the statistical framework of Exponential-family Random Graph Models (ERGMs). EpiModel supports three classes of epidemic models:

1. **Deterministic Compartmental Models (DCMs)** -- solved via ordinary differential equations (using the `deSolve` package). Population-level, continuous-time, no stochasticity. Good for quick explorations and closed-form analysis.

2. **Stochastic Individual Contact Models (ICMs)** -- discrete-time agent-based models with random mixing. Each individual is represented but there is no persistent network structure; contacts are drawn randomly each time step.

3. **Stochastic Network Models** -- the primary focus. These use ERGMs/TERGMs (Temporal ERGMs) to represent dynamic contact networks where partnerships form, persist, and dissolve over time according to statistically estimated rules. Epidemics then propagate across these networks. This is what makes EpiModel unique.

## Why Network Models Matter for Epidemics

Traditional epidemic models assume random ("mass-action") mixing: any individual can contact any other individual with equal probability. Real contact patterns are structured -- people have a fixed number of partners, partnerships have durations, there is clustering, assortative mixing by attributes like age or race, and concurrency (overlapping partnerships).

Network models capture these structural features because:
- **Degree distributions** control how many partners each person has (not everyone has the same number).
- **Relational durations** mean partnerships persist over time, creating sustained transmission opportunities.
- **Concurrency** (having multiple simultaneous partners) dramatically accelerates epidemic spread compared to serial monogamy, even when the total number of contacts is identical.
- **Assortative mixing** (e.g., by age, race, risk group) creates subpopulation-level dynamics where epidemics can concentrate in or bridge between groups.
- **Clustering** (triangles, shared partners) can slow epidemics by trapping transmission within already-infected clusters.
- **Feedback** between network and epidemic is possible: behavioral changes in response to infection status (serosorting), demographic turnover that reshapes the network, and interventions that alter contact patterns.

## Key Technical Concepts

### Built-in Disease Types

EpiModel has three built-in compartmental structures: **SI** (Susceptible-Infected), **SIR** (Susceptible-Infected-Recovered), and **SIS** (Susceptible-Infected-Susceptible). These can be run in one-group or two-group configurations, with or without demography (open/closed populations). The extension API allows arbitrary disease states.

### The Network Model Pipeline (5 Steps)

This is the core workflow for any network epidemic model:

**Step 1: Initialize the network**
```r
nw <- network_initialize(n = 500)
nw <- set_vertex_attribute(nw, "risk", rep(0:1, each = 250))
```

**Step 2: Parameterize formation and dissolution**
```r
formation <- ~edges + concurrent + degrange(from = 4)
target.stats <- c(175, 110, 0)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
```
Key formula: `target edges = mean_degree * N / 2`

**Step 3: Estimate the network model**
```r
est <- netest(nw, formation, target.stats, coef.diss)
```

**Step 4: Diagnose model fit**
```r
dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 4)
plot(dx)                        # formation diagnostics
plot(dx, type = "duration")     # dissolution diagnostics
```

**Step 5: Simulate the epidemic**
```r
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.1)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsims = 5, nsteps = 500)
sim <- netsim(est, param, init, control)
```

### ERGM Terms Reference

ERGM terms define network structure in the formation model. Common terms used in EpiModel:

**Dyad-independent terms** (faster estimation, edge probabilities don't depend on other edges):
- `edges` -- total edge count; controls density. Target = `mean_degree * N / 2`.
- `nodefactor("attr")` -- differential activity by attribute level. Target = `mean_degree_of_level * level_size`.
- `nodematch("attr")` -- assortative mixing (homophily). Target = `proportion_matched * total_edges`. Use `diff=TRUE` for per-level homophily.
- `nodemix("attr")` -- full mixing matrix by attribute combinations.
- `nodecov("attr")` -- sum of continuous attribute values across edges.
- `absdiff("attr")` -- sum of absolute differences in continuous attribute across edges (e.g., age homophily).

**Dyad-dependent terms** (require MCMC estimation):
- `degree(k)` -- number of nodes with exactly degree k.
- `degree(0:k, by = "group")` -- group-stratified degree distribution.
- `concurrent` (alias for `mindegree(2)`) -- number of nodes with 2+ simultaneous partners.
- `degrange(from = k)` -- constrain maximum degree (target = 0 means no nodes with degree >= k).
- `triangles` -- clustering.
- `gwesp` -- geometrically weighted edgewise shared partners (flexible clustering).

**Dissolution models** are limited to dyad-independent terms with `offset()`:
```r
~offset(edges)                                     # homogeneous duration
~offset(edges) + offset(nodematch("race"))         # differential duration by attribute match
```

### Transmission Probability Formula

Per-partnership per-timestep transmission probability:
```
finalProb = 1 - (1 - transProb)^actRate
```
Where `transProb` is the per-act transmission probability and `actRate` is the number of acts per partnership per time step.

### Dissolution Coefficient Adjustment for Mortality

When modeling open populations with demographic turnover, dissolution coefficients must be adjusted for exogenous edge loss from node departure:
```r
coef.diss <- dissolution_coefs(~offset(edges), duration = 40, d.rate = 0.001)
```

### Network Estimation Methods

- **Approximation method** (`edapprox = TRUE`, default): Uses ERGM for cross-sectional fit with analytic dissolution adjustment. Fast and stable. Best when durations > 50 time steps.
- **Direct STERGM** (`edapprox = FALSE`): Full Separable Temporal ERGM estimation. Necessary when durations are short (~10-25 time steps). Slower but more accurate for short durations.

## The EpiModel Extension API

The extension API allows building custom research-level models with arbitrary disease states, demographic processes, interventions, and feedback mechanisms. This is the most powerful and most used feature of EpiModel for applied research.

### Three Core API Rules

**Rule 1 -- Function Signature:** Every module function must have exactly this form:
```r
module_function <- function(dat, at) {
  # ... processes that update dat ...
  return(dat)
}
```
- `dat` is the master data object containing all model state.
- `at` is the current time step.

**Rule 2 -- Accessor Functions:** All data on `dat` is accessed through getter/setter functions:

| Function | Purpose |
|---|---|
| `get_attr(dat, "name")` | Read a nodal attribute vector |
| `set_attr(dat, "name", value)` | Write a nodal attribute vector |
| `get_param(dat, "name")` | Read a model parameter |
| `get_epi(dat, "name", at)` | Read an epidemic summary statistic |
| `set_epi(dat, "name", at, value)` | Write an epidemic summary statistic |
| `append_core_attr(dat, at, n.new)` | Initialize required attributes for new nodes |
| `append_attr(dat, "name", value, n.new)` | Append a custom attribute for new nodes |
| `discord_edgelist(dat, at)` | Get susceptible-infected discordant edge pairs |
| `set_transmat(dat, del, at)` | Record transmission events |

**Rule 3 -- Three-Step Design Pattern:** Each module: (a) reads inputs from `dat`, (b) performs stochastic processes on individual nodes, (c) writes updated state back to `dat`.

### Module Classification

- **Standard modules** (typically not modified): initialization, network resimulation, bookkeeping.
- **Flexible modules** (designed to be customized): infection, disease progression, departures, arrivals.

### Required Nodal Attributes

Every node must have: `active`, `status`, `infTime`, `entrTime`, `exitTime`.

### Departure Module API Rule

When nodes depart (die), two attributes MUST be set:
```r
active[idsDepts] <- 0
exitTime[idsDepts] <- at
```

### Arrival Module API Rule

When nodes arrive (are born), use `append_core_attr` first (handles `active`, `entrTime`, `exitTime`), then `append_attr` for every other attribute:
```r
dat <- append_core_attr(dat, at = at, n.new = nArrivals)
dat <- append_attr(dat, "status", "s", nArrivals)
dat <- append_attr(dat, "infTime", NA, nArrivals)
dat <- append_attr(dat, "age", 0, nArrivals)
```

### The statusTime Guard Pattern

To prevent multiple disease state transitions in a single time step, use a `statusTime` attribute that tracks when status last changed. Only nodes where `statusTime < at` are eligible for transition:
```r
idsElig <- which(active == 1 & status == "e" & statusTime < at)
```
After each transition, update: `statusTime[transitioned_ids] <- at`.

### Extension Model control.net Setup

For extension models, `type` must be `NULL` and custom module functions are passed by name:
```r
control <- control.net(
  type = NULL,
  nsteps = 500, nsims = 10, ncores = 5,
  infection.FUN = infect,
  progress.FUN = progress,
  aging.FUN = aging,
  departures.FUN = dfunc,
  arrivals.FUN = afunc,
  dx.FUN = dx_covid,
  resimulate.network = TRUE,
  tergmLite = TRUE
)
```

### Key control.net Parameters

- `resimulate.network = TRUE` -- required when network depends on changing attributes (feedback, demography).
- `tergmLite = TRUE` -- 20-50x faster simulation but does not retain full network object (no network plots at specific time points).
- `epi.by = "attr"` -- automatically generates attribute-stratified epidemic statistics (e.g., `i.num.risk0`, `i.num.risk1`).
- `nwstats.formula` -- track additional network statistics during simulation.

### The `group` Attribute (Special)

The `group` attribute (values `1` and `2` only) is treated specially by EpiModel's built-in models:
- Automatically stratifies epidemic outputs (e.g., `i.num`, `i.num.g2`).
- Allows group-specific parameters via `.g2` suffix: `inf.prob.g2`, `rec.rate.g2`.

For more flexible attribute handling (more than 2 groups, non-stratified parameters), use a regular attribute like `"risk"` with `epi.by` in `control.net`.

## Disease Model Examples in the Training Materials

### Basic Built-in Models
- **SI**: Simple epidemic growth (e.g., HIV without treatment).
- **SIR**: Epidemic with recovery and immunity (e.g., measles).
- **SIS**: Epidemic with reinfection (e.g., bacterial STIs like gonorrhea).

### Custom DCM Extensions
- **SEIR**: Adds exposed/latent stage (e.g., Ebola). Custom ODE function passed to `control.dcm(new.mod = SEIR)`.
- **Variable Mixing (Q-statistic)**: Models assortative/disassortative mixing with continuous Q parameter.

### Custom Network Model Extensions
- **SEIR on network**: Adds latent "e" state. Infection transitions S -> E (not S -> I). Separate progression module handles E -> I and I -> R.
- **COVID with demography**: SEIR plus age-structured population, age-specific mortality, births, and disease-induced excess mortality.
- **COVID with screening and interventions**: States S, E, A (asymptomatic), Ip (preclinical), Ic (clinical symptomatic), R. Includes diagnostic screening module with PCR sensitivity, symptomatic vs. asymptomatic testing rates, case isolation through act rate reduction, and age-dependent clinical pathway probability.

## Feedback Mechanisms

Network models can incorporate bidirectional feedback between epidemic and network:

1. **Demography**: Births and deaths change population size and composition, reshaping the network.
2. **Serosorting**: Using `status` as an ERGM term causes the network to dynamically respond to changing prevalence. Differential degree by disease status and preferential within-status partnering are modeled via `nodefactor("status")` and `nodematch("status")`.
3. **Behavioral interventions**: Diagnosed individuals can have reduced act rates. Interventions can be time-dependent (starting at a specific time step).
4. **Built-in vaccine/condom intervention**: `inter.eff` (efficacy as relative reduction in `inf.prob`) and `inter.start` (time step when intervention begins).

When `status` is used in the ERGM formation formula, it must be initialized directly on the network object (not via `init.net`), and `init.net()` is called empty as a placeholder.

In feedback models, network statistics are NOT preserved at their target levels over time -- what is preserved is the log-odds conditional on other terms. As prevalence changes, edge counts and degree distributions shift dynamically.

## Output Analysis

### Data Extraction
```r
print(sim)                                    # summary of inputs and outputs
summary(sim, at = 500)                        # statistics at a specific time step
as.data.frame(sim)                            # per-simulation raw data
as.data.frame(sim, out = "mean")              # means across simulations
get_network(sim, sim = 1)                     # extract networkDynamic object
get_transmat(sim, sim = 1)                    # transmission matrix (who infected whom)
as.phylo.transmat(get_transmat(sim))          # convert to phylogenetic tree
mutate_epi(sim, ir = (si.flow / s.num) * 100) # compute derived epidemic statistics
```

### Visualization
```r
plot(sim)                                      # compartment sizes over time
plot(sim, y = c("si.flow", "is.flow"))        # flow sizes
plot(sim, type = "network", col.status = TRUE) # network colored by disease status
plot(dx)                                       # formation diagnostics
plot(dx, type = "duration")                    # dissolution fit
```

## Other Broader Architecture Notes

### Shared Plotting Infrastructure

`plot.icm` (`R/plot.icm.R`) and `plot_netsim_epi` (`R/plot.netsim.R`) share helper functions `draw_means()` and `draw_qnts()` defined in `R/plot.R`. Changes to plotting behavior (e.g., axis offsets, coordinate handling) must be coordinated across all three files. `plot.dcm` (`R/plot.dcm.R`) is independent and uses `control$timesteps` directly for x-coordinates.

`as.data.frame.netsim` delegates entirely to `as.data.frame.icm`, so changes to the ICM method automatically apply to netsim. The `as.epi.data.frame` validator in the same file checks structural consistency of the output.

### Time Handling Differences Across Model Classes

DCM and ICM/netsim handle time differently in their output objects:

- **DCM**: `control$timesteps` stores the explicit time vector (may be non-integer when `dt < 1`). `control$nsteps` is the max time value. `control$nruns` is the number of parameter sensitivity runs. Plot and data frame methods use `timesteps` directly.
- **ICM/netsim**: Time is implicit via row indices. `control$nsteps` is both the number of rows and the max time step. `control$nsims` is the number of stochastic simulations. `control$start` provides a time offset (default 1). `as.data.frame.icm` builds the time column from `start` and `nsteps`.

When modifying functions that touch time or epi structure, verify behavior across all three classes and check that `as.data.frame`, `plot`, and `summary` methods still work correctly.

### Key Source File Groups

- **`R/net.inputs.R`** (largest file): Validates and processes all network model inputs (`param.net`, `init.net`, `control.net`)
- **`R/net.fn.utils.R`**, **`R/net.fn.accessor.R`**: Core utilities and attribute accessors for the internal `dat` simulation object
- **`R/get.R`**: Public API for extracting simulation results (`get_param`, `get_attr`, `get_epi`, etc.)
- **`R/edgelists.R`**: Edgelist manipulation and cumulative edgelist tracking
- **`R/dcm.mods.R`**: Built-in ODE systems for standard SI/SIR/SIS compartmental models

### The `dat` Object

Network simulations pass a central `dat` list object through all modules. Key sub-elements:
- `dat$attr` - Node-level attributes (disease status, demographics, etc.)
- `dat$epi` - Epidemic tracking variables (counts per time step)
- `dat$nwparam` - Network model parameters
- `dat$param` / `dat$init` / `dat$control` - Simulation configuration
- `dat$run` - Runtime state (current time step, simulation metadata)
- `dat$el` - Edgelists for each network layer (when using tergmLite mode)
- `dat$net` - Network objects (when not using tergmLite mode)

### Compiled Code (src/)

C/C++ source provides additional ERGM terms (like ergm's built-in nodematch or nodefactor):
- `changestats.users.c` - Custom ERGM terms (`fuzzynodematch`, `absdiffby`, `absdiffnodemix`)
- `updatenw.cpp` - Rcpp-based network update operations
- `RcppExports.cpp` - Auto-generated Rcpp bindings (do not edit manually)

### Multi-Network Support

Models can use multiple network layers (e.g., sexual + needle-sharing networks). Network parameters, edgelists, and diagnostics are indexed by network number. See `test-multinets.R` for examples.

## Dependencies

The package heavily depends on the Statnet ecosystem: `ergm`, `tergm`, `network`, `networkDynamic`, `networkLite`, `statnet.common`. Understanding ERGM model specification is important for working with network model code.

## Testing Without Full Dependencies
The Statnet dependencies (ergm, tergm, etc.) are heavy. For changes to core module code (e.g., changes to the infection or recovery module), plotting, or data frame methods, test the logic directly rather than attempting full package installation. Push and rely on CI for integration tests.
