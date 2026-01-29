# Plot Dynamic Network Model Diagnostics

Plots dynamic network model diagnostics calculated in
[`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md).

## Usage

``` r
# S3 method for class 'netdx'
plot(
  x,
  type = "formation",
  method = "l",
  sims = NULL,
  stats = NULL,
  duration.imputed = TRUE,
  sim.lines = FALSE,
  sim.col = NULL,
  sim.lwd = NULL,
  mean.line = TRUE,
  mean.smooth = TRUE,
  mean.col = NULL,
  mean.lwd = 2,
  mean.lty = 1,
  qnts = 0.5,
  qnts.col = NULL,
  qnts.alpha = 0.5,
  qnts.smooth = TRUE,
  targ.line = TRUE,
  targ.col = NULL,
  targ.lwd = 2,
  targ.lty = 2,
  plots.joined = NULL,
  legend = NULL,
  grid = FALSE,
  ...
)
```

## Arguments

- x:

  An `EpiModel` object of class `netdx`.

- type:

  Plot type, with options of `"formation"` for network model formation
  statistics, `"duration"` for dissolution model statistics for average
  edge duration, or `"dissolution"` for dissolution model statistics for
  proportion of ties dissolved per time step.

- method:

  Plot method, with options of `"l"` for line plots and `"b"` for box
  plots.

- sims:

  A vector of simulation numbers to plot.

- stats:

  Statistics to plot. For `type = "formation"`, `stats` are among those
  specified in the call to
  [`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md);
  for `type = "duration", "dissolution"`, `stats` are among those of the
  dissolution model (without
  [`offset()`](https://rdrr.io/r/stats/offset.html)). The default is to
  plot all statistics.

- duration.imputed:

  If `type = "duration"`, a logical indicating whether or not to impute
  starting times for relationships extant at the start of the
  simulation. Defaults to `TRUE` when `type = "duration"`.

- sim.lines:

  If `TRUE`, plot individual simulation lines. Default is to plot lines
  for one-group models but not for two-group models.

- sim.col:

  Vector of any standard R color format for simulation lines.

- sim.lwd:

  Line width for simulation lines.

- mean.line:

  If `TRUE`, plot mean of simulations across time.

- mean.smooth:

  If `TRUE`, use a loess smoother on the mean line.

- mean.col:

  Vector of any standard R color format for mean lines.

- mean.lwd:

  Line width for mean lines.

- mean.lty:

  Line type for mean lines.

- qnts:

  If numeric, plot polygon of simulation quantiles based on the range
  implied by the argument (see details). If `FALSE`, suppress polygon
  from plot.

- qnts.col:

  Vector of any standard R color format for polygons.

- qnts.alpha:

  Transparency level for quantile polygons, where 0 = transparent and 1
  = opaque (see `adjustcolor` function).

- qnts.smooth:

  If `TRUE`, use a loess smoother on quantile polygons.

- targ.line:

  If `TRUE`, plot target or expected value line for the statistic of
  interest.

- targ.col:

  Vector of standard R colors for target statistic lines, with default
  colors based on `RColorBrewer` color palettes.

- targ.lwd:

  Line width for the line showing the target statistic values.

- targ.lty:

  Line type for the line showing the target statistic values.

- plots.joined:

  If `TRUE`, combine all statistics in one plot, versus one plot per
  statistic if `FALSE`.

- legend:

  If `TRUE`, plot default legend.

- grid:

  If `TRUE`, a grid is added to the background of plot (see
  [`grid()`](https://rdrr.io/r/graphics/grid.html) for details), with
  default of nx by ny.

- ...:

  Additional arguments to pass.

## Details

The plot function for `netdx` objects will generate plots of two types
of model diagnostic statistics that run as part of the diagnostic tools
within that function. The `formation` plot shows the summary statistics
requested in `nwstats.formula`, where the default includes those
statistics in the network model formation formula specified in the
original call to
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md).

The `duration` plot shows the average age of existing edges at each time
step, up until the maximum time step requested. The age is used as an
estimator of the average duration of edges in the equilibrium state.
When `duration.imputed = FALSE`, edges that exist at the beginning of
the simulation are assumed to start with an age of 1, yielding a burn-in
period before the observed mean approaches its target. When
`duration.imputed = TRUE`, expected ages prior to the start of the
simulation are calculated from the dissolution model, typically
eliminating the need for a burn-in period.

The `dissolution` plot shows the proportion of the extant ties that are
dissolved at each time step, up until the maximum time step requested.
Typically, the proportion of ties that are dissolved is the reciprocal
of the mean relational duration. This plot thus contains similar
information to that in the duration plot, but should reach its expected
value more quickly, since it is not subject to censoring.

The `plots.joined` argument will control whether the statistics are
joined in one plot or plotted separately, assuming there are multiple
statistics in the model. The default is based on the number of network
statistics requested. The layout of the separate plots within the larger
plot window is also based on the number of statistics.

## See also

[`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Network initialization and model parameterization
nw <- network_initialize(n = 500)
nw <- set_vertex_attribute(nw, "sex", rbinom(500, 1, 0.5))
formation <- ~edges + nodematch("sex")
target.stats <- c(500, 300)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) +
                  offset(nodematch("sex")), duration = c(50, 40))

# Estimate the model
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Static diagnostics
dx1 <- netdx(est, nsims = 1e4, dynamic = FALSE,
             nwstats.formula = ~edges + meandeg + concurrent +
                                nodefactor("sex", levels = NULL) +
                                nodematch("sex"))
dx1

# Plot diagnostics
plot(dx1)
plot(dx1, stats = c("edges", "concurrent"), mean.col = "black",
     sim.lines = TRUE, plots.joined = FALSE)
plot(dx1, stats = "edges", method = "b",
     col = "seagreen3", grid = TRUE)

# Dynamic diagnostics
dx2 <- netdx(est, nsims = 10, nsteps = 500,
             nwstats.formula = ~edges + meandeg + concurrent +
                                nodefactor("sex", levels = NULL) +
                                nodematch("sex"))
dx2

# Formation statistics plots, joined and separate
plot(dx2, grid = TRUE)
plot(dx2, type = "formation", plots.joined = TRUE)
plot(dx2, type = "formation", sims = 1, plots.joined = TRUE,
     qnts = FALSE, sim.lines = TRUE, mean.line = FALSE)
plot(dx2, type = "formation", plots.joined = FALSE,
     stats = c("edges", "concurrent"), grid = TRUE)

plot(dx2, method = "b", col = "bisque", grid = TRUE)
plot(dx2, method = "b", stats = "meandeg", col = "dodgerblue")

# Duration statistics plot
par(mfrow = c(1, 2))
# With duration imputed
plot(dx2, type = "duration", sim.line = TRUE, sim.lwd = 0.3,
     targ.lty = 1, targ.lwd = 0.5)
# Without duration imputed
plot(dx2, type = "duration", sim.line = TRUE, sim.lwd = 0.3,
     targ.lty = 1, targ.lwd = 0.5, duration.imputed = FALSE)

# Dissolution statistics plot
plot(dx2, type = "dissolution", qnts = 0.25, grid = TRUE)
plot(dx2, type = "dissolution", method = "b", col = "pink1")
} # }
```
