% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.netdx.R
\name{plot.netdx}
\alias{plot.netdx}
\title{Plot Dynamic Network Model Diagnostics}
\usage{
\method{plot}{netdx}(
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
}
\arguments{
\item{x}{An \code{EpiModel} object of class \code{netdx}.}

\item{type}{Plot type, with options of \code{"formation"} for network
model formation statistics, \code{"duration"} for dissolution model
statistics for average edge duration, or \code{"dissolution"} for
dissolution model statistics for proportion of ties dissolved per time
step.}

\item{method}{Plot method, with options of \code{"l"} for line plots and
\code{"b"} for box plots.}

\item{sims}{A vector of simulation numbers to plot.}

\item{stats}{Statistics to plot. For \code{type = "formation"}, \code{stats}
are among those specified in the call to \code{\link{netdx}};
for \code{type = "duration", "dissolution"}, \code{stats} are among
those of the dissolution model (without \code{offset()}). The default
is to plot all statistics.}

\item{duration.imputed}{If \code{type = "duration"}, a logical indicating
whether or not to impute starting times for relationships extant at
the start of the simulation. Defaults to \code{TRUE} when
\code{type = "duration"}.}

\item{sim.lines}{If \code{TRUE}, plot individual simulation lines. Default is
to plot lines for one-group models but not for two-group models.}

\item{sim.col}{Vector of any standard R color format for simulation lines.}

\item{sim.lwd}{Line width for simulation lines.}

\item{mean.line}{If \code{TRUE}, plot mean of simulations across time.}

\item{mean.smooth}{If \code{TRUE}, use a loess smoother on the mean line.}

\item{mean.col}{Vector of any standard R color format for mean lines.}

\item{mean.lwd}{Line width for mean lines.}

\item{mean.lty}{Line type for mean lines.}

\item{qnts}{If numeric, plot polygon of simulation quantiles based on the
range implied by the argument (see details). If \code{FALSE}, suppress
polygon from plot.}

\item{qnts.col}{Vector of any standard R color format for polygons.}

\item{qnts.alpha}{Transparency level for quantile polygons, where 0 =
transparent and 1 = opaque (see \code{adjustcolor} function).}

\item{qnts.smooth}{If \code{TRUE}, use a loess smoother on quantile polygons.}

\item{targ.line}{If \code{TRUE}, plot target or expected value line for
the statistic of interest.}

\item{targ.col}{Vector of standard R colors for target statistic lines, with
default colors based on \code{RColorBrewer} color palettes.}

\item{targ.lwd}{Line width for the line showing the target statistic values.}

\item{targ.lty}{Line type for the line showing the target statistic values.}

\item{plots.joined}{If \code{TRUE}, combine all statistics in one plot,
versus one plot per statistic if \code{FALSE}.}

\item{legend}{If \code{TRUE}, plot default legend.}

\item{grid}{If \code{TRUE}, a grid is added to the background of plot
(see \code{\link{grid}} for details), with default of nx by ny.}

\item{...}{Additional arguments to pass.}
}
\description{
Plots dynamic network model diagnostics calculated in
\code{\link{netdx}}.
}
\details{
The plot function for \code{netdx} objects will generate plots of two types
of model diagnostic statistics that run as part of the diagnostic tools
within that function. The \code{formation} plot shows the summary statistics
requested in \code{nwstats.formula}, where the default includes those
statistics in the network model formation formula specified in the original
call to \code{\link{netest}}.

The \code{duration} plot shows the average age of existing edges at each time
step, up until the maximum time step requested. The age is used as an
estimator of the average duration of edges in the equilibrium state. When
\code{duration.imputed = FALSE}, edges that exist at the beginning of the
simulation are assumed to start with an age of 1, yielding a burn-in period
before the observed mean approaches its target.  When
\code{duration.imputed = TRUE}, expected ages prior to the start of the
simulation are calculated from the dissolution model, typically eliminating
the need for a burn-in period.

The \code{dissolution} plot shows the proportion of the extant ties that are
dissolved at each time step, up until the maximum time step requested.
Typically, the proportion of ties that are dissolved is the reciprocal of the
mean relational duration. This plot thus contains similar information to that
in the duration plot, but should reach its expected value more quickly, since
it is not subject to censoring.

The \code{plots.joined} argument will control whether the statistics
are joined in one plot or plotted separately, assuming there are multiple
statistics in the model. The default is based on the number of network
statistics requested. The layout of the separate plots within the larger plot
window is also based on the number of statistics.
}
\examples{
\dontrun{
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
}

}
\seealso{
\code{\link{netdx}}
}
\keyword{plot}
