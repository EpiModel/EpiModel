# Plot Epidemic Model Results From a Netsim Data.Frame

This function is a wrapper around `plot.netsim` accepting a `data.frame`
obtain with `as.data.frame(netsim_object)`.

## Usage

``` r
# S3 method for class 'epi.data.frame'
plot(
  x,
  y = NULL,
  sims = NULL,
  legend = NULL,
  mean.col = NULL,
  qnts.col = NULL,
  sim.lwd = NULL,
  sim.col = NULL,
  sim.alpha = NULL,
  popfrac = FALSE,
  qnts = 0.5,
  qnts.alpha = 0.5,
  qnts.smooth = TRUE,
  mean.line = TRUE,
  mean.smooth = TRUE,
  add = FALSE,
  mean.lwd = 2,
  mean.lty = 1,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  sim.lines = FALSE,
  grid = FALSE,
  leg.cex = 0.8,
  ...
)
```

## Arguments

- x:

  A `data.frame` obtain with `as.data.frame(netsim_object)`.

- y:

  Output compartments or flows from `netsim` object to plot.

- sims:

  If `type="epi"` or `"formation"`, a vector of simulation numbers to
  plot. If `type="network"`, a single simulation number for which to
  plot the network, or else `"min"` to plot the simulation number with
  the lowest disease prevalence, `"max"` for the simulation with the
  highest disease prevalence, or `"mean"` for the simulation with the
  prevalence closest to the mean across simulations at the specified
  time step.

- legend:

  If `TRUE`, plot default legend.

- mean.col:

  Vector of any standard R color format for mean lines.

- qnts.col:

  Vector of any standard R color format for polygons.

- sim.lwd:

  Line width for simulation lines.

- sim.col:

  Vector of any standard R color format for simulation lines.

- sim.alpha:

  Transparency level for simulation lines, where 0 = transparent and 1 =
  opaque (see `adjustcolor` function).

- popfrac:

  If `TRUE`, plot prevalence of values rather than numbers (see
  details).

- qnts:

  If numeric, plot polygon of simulation quantiles based on the range
  implied by the argument (see details). If `FALSE`, suppress polygon
  from plot.

- qnts.alpha:

  Transparency level for quantile polygons, where 0 = transparent and 1
  = opaque (see `adjustcolor` function).

- qnts.smooth:

  If `TRUE`, use a loess smoother on quantile polygons.

- mean.line:

  If `TRUE`, plot mean of simulations across time.

- mean.smooth:

  If `TRUE`, use a loess smoother on the mean line.

- add:

  If `TRUE`, new plot window is not called and lines are added to
  existing plot window.

- mean.lwd:

  Line width for mean lines.

- mean.lty:

  Line type for mean lines.

- xlim:

  the x limits (x1, x2) of the plot. Note that `x1 > x2` is allowed and
  leads to a ‘reversed axis’.

  The default value, `NULL`, indicates that the range of the
  [finite](https://rdrr.io/r/base/is.finite.html) values to be plotted
  should be used.

- ylim:

  the y limits of the plot.

- main:

  a main title for the plot, see also
  [`title`](https://rdrr.io/r/graphics/title.html).

- xlab:

  a label for the x axis, defaults to a description of `x`.

- ylab:

  a label for the y axis, defaults to a description of `y`.

- sim.lines:

  If `TRUE`, plot individual simulation lines. Default is to plot lines
  for one-group models but not for two-group models.

- grid:

  If `TRUE`, a grid is added to the background of plot (see
  [`grid()`](https://rdrr.io/r/graphics/grid.html) for details), with
  default of nx by ny.

- leg.cex:

  Legend scale size.

- ...:

  Additional arguments to pass.
