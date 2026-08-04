#' @title Plot Dynamic Network Model Diagnostics
#'
#' @description Plots dynamic network model diagnostics calculated in
#'              [netdx()].
#'
#' @param x An `EpiModel` object of class `netdx`.
#' @param type Plot type, with options of `"formation"` for network
#'        model formation statistics, `"duration"` for dissolution model
#'        statistics for average edge duration, `"dissolution"` for
#'        dissolution model statistics for proportion of ties dissolved per time
#'        step, or `"cumldeg"` for the cumulative degree distribution.
#' @param method Plot method, with options of `"l"` for line plots and
#'        `"b"` for box plots. For `type = "cumldeg"`, `"b"` gives a bar plot
#'        of the degree distribution.
#' @param sims A vector of simulation numbers to plot.
#' @param stats Statistics to plot. For `type = "formation"`, `stats`
#'        are among those specified in the call to [netdx()];
#'        for `type = "duration", "dissolution"`, `stats` are among
#'        those of the dissolution model (without `offset()`). The default
#'        is to plot all statistics.
#' @param plots.joined If `TRUE`, combine all statistics in one plot,
#'        versus one plot per statistic if `FALSE`.
#' @param momentary If `TRUE`, add the momentary degree distribution to the
#'        cumulative degree plot for comparison. Used only when
#'        `type = "cumldeg"`.
#' @param window A vector of length 2 giving the first and last time step over
#'        which degrees are counted, with the default of the full simulation.
#'        Used only when `type = "cumldeg"`.
#' @param maxdeg Highest degree to show on the x-axis, with degrees above it
#'        collapsed into a final `"maxdeg+"` category. Used only when
#'        `type = "cumldeg"`.
#' @param unique.partners If `TRUE`, a dyad that forms, dissolves, and then
#'        re-forms within the window counts as one cumulative partner rather
#'        than two. Used only when `type = "cumldeg"`.
#' @inheritParams plot.netsim
#' @inheritParams graphics::plot
#'
#' @details
#' The plot function for `netdx` objects will generate plots of two types
#' of model diagnostic statistics that run as part of the diagnostic tools
#' within that function. The `formation` plot shows the summary statistics
#' requested in `nwstats.formula`, where the default includes those
#' statistics in the network model formation formula specified in the original
#' call to [netest()].
#'
#' The `duration` plot shows the average age of existing edges at each time
#' step, up until the maximum time step requested. The age is used as an
#' estimator of the average duration of edges in the equilibrium state. When
#' `duration.imputed = FALSE`, edges that exist at the beginning of the
#' simulation are assumed to start with an age of 1, yielding a burn-in period
#' before the observed mean approaches its target.  When
#' `duration.imputed = TRUE`, expected ages prior to the start of the
#' simulation are calculated from the dissolution model, typically eliminating
#' the need for a burn-in period.
#'
#' The `dissolution` plot shows the proportion of the extant ties that are
#' dissolved at each time step, up until the maximum time step requested.
#' Typically, the proportion of ties that are dissolved is the reciprocal of the
#' mean relational duration. This plot thus contains similar information to that
#' in the duration plot, but should reach its expected value more quickly, since
#' it is not subject to censoring.
#'
#' The `plots.joined` argument will control whether the statistics
#' are joined in one plot or plotted separately, assuming there are multiple
#' statistics in the model. The default is based on the number of network
#' statistics requested. The layout of the separate plots within the larger plot
#' window is also based on the number of statistics.
#'
#' The `cumldeg` plot shows the cumulative degree distribution: the proportion
#' of nodes with each count of partners accumulated over the simulation, rather
#' than the count held at one point in time. It requires a `netdx` object run
#' with `keep.tedgelist = TRUE`, since it is calculated from the timed edgelist
#' with [get_degree_dist()]. Setting `momentary = TRUE` adds the momentary
#' degree distribution over the same window, calculated from the same edgelist,
#' which is the distribution summarized by the `degree(k)` terms of the
#' formation plot. Two models with the same mean degree and mean partnership
#' duration have the same mean cumulative degree even when their momentary
#' degree distributions differ, as the example below shows for models with high
#' and low concurrency. Legend entries report the mean of each distribution.
#'
#' Because the cumulative degree grows with the length of the observation
#' window, use `window` to restrict it to a substantively meaningful period
#' when comparing across models. The lines or bars show the mean across
#' simulations, with `qnts` giving the inter-simulation quantile band; the
#' `targ.line`, `mean.smooth`, `qnts.smooth`, and `plots.joined` arguments do
#' not apply to this plot type.
#'
#' @method plot netdx
#' @export
#'
#' @keywords plot
#' @seealso [netdx()], and [get_degree_dist()] for the degree distributions
#'   behind `type = "cumldeg"`.
#'
#' @examples
#' \dontrun{
#' # Network initialization and model parameterization
#' nw <- network_initialize(n = 500)
#' nw <- set_vertex_attribute(nw, "sex", rbinom(500, 1, 0.5))
#' formation <- ~edges + nodematch("sex")
#' target.stats <- c(500, 300)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges) +
#'                   offset(nodematch("sex")), duration = c(50, 40))
#'
#' # Estimate the model
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Static diagnostics
#' dx1 <- netdx(est, nsims = 1e4, dynamic = FALSE,
#'              nwstats.formula = ~edges + meandeg + concurrent +
#'                                 nodefactor("sex", levels = NULL) +
#'                                 nodematch("sex"))
#' dx1
#'
#' # Plot diagnostics
#' plot(dx1)
#' plot(dx1, stats = c("edges", "concurrent"), mean.col = "black",
#'      sim.lines = TRUE, plots.joined = FALSE)
#' plot(dx1, stats = "edges", method = "b",
#'      col = "seagreen3", grid = TRUE)
#'
#' # Dynamic diagnostics
#' dx2 <- netdx(est, nsims = 10, nsteps = 500,
#'              nwstats.formula = ~edges + meandeg + concurrent +
#'                                 nodefactor("sex", levels = NULL) +
#'                                 nodematch("sex"))
#' dx2
#'
#' # Formation statistics plots, joined and separate
#' plot(dx2, grid = TRUE)
#' plot(dx2, type = "formation", plots.joined = TRUE)
#' plot(dx2, type = "formation", sims = 1, plots.joined = TRUE,
#'      qnts = FALSE, sim.lines = TRUE, mean.line = FALSE)
#' plot(dx2, type = "formation", plots.joined = FALSE,
#'      stats = c("edges", "concurrent"), grid = TRUE)
#'
#' plot(dx2, method = "b", col = "bisque", grid = TRUE)
#' plot(dx2, method = "b", stats = "meandeg", col = "dodgerblue")
#'
#' # Duration statistics plot
#' par(mfrow = c(1, 2))
#' # With duration imputed
#' plot(dx2, type = "duration", sim.line = TRUE, sim.lwd = 0.3,
#'      targ.lty = 1, targ.lwd = 0.5)
#' # Without duration imputed
#' plot(dx2, type = "duration", sim.line = TRUE, sim.lwd = 0.3,
#'      targ.lty = 1, targ.lwd = 0.5, duration.imputed = FALSE)
#'
#' # Dissolution statistics plot
#' plot(dx2, type = "dissolution", qnts = 0.25, grid = TRUE)
#' plot(dx2, type = "dissolution", method = "b", col = "pink1")
#'
#' # Cumulative degree distribution: two models with the same mean degree and
#' # the same mean partnership duration, differing only in how much
#' # concurrency (overlapping partnerships) they allow
#' nw <- network_initialize(n = 500)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#'
#' est.hi <- netest(nw, formation = ~edges + concurrent,
#'                  target.stats = c(200, 150), coef.diss = coef.diss,
#'                  verbose = FALSE)
#' est.lo <- netest(nw, formation = ~edges + concurrent,
#'                  target.stats = c(200, 40), coef.diss = coef.diss,
#'                  verbose = FALSE)
#'
#' dx.hi <- netdx(est.hi, nsims = 5, nsteps = 250, keep.tedgelist = TRUE,
#'                nwstats.formula = ~edges + concurrent + degree(0:3))
#' dx.lo <- netdx(est.lo, nsims = 5, nsteps = 250, keep.tedgelist = TRUE,
#'                nwstats.formula = ~edges + concurrent + degree(0:3))
#'
#' # The momentary degree distributions are far apart: about half the nodes in
#' # the high-concurrency model have no partner at any given time and a quarter
#' # have two, while most nodes in the low-concurrency model have exactly one.
#' # Over 250 time steps both accumulate about 8.5 partners per node.
#' par(mfrow = c(1, 2))
#' plot(dx.hi, type = "cumldeg", momentary = TRUE, xlim = c(0, 25),
#'      main = "High concurrency")
#' plot(dx.lo, type = "cumldeg", momentary = TRUE, xlim = c(0, 25),
#'      main = "Low concurrency")
#'
#' # As bars, with the tail of the distribution collapsed into one category
#' # and the inter-simulation range shown
#' par(mfrow = c(1, 1))
#' plot(dx.hi, type = "cumldeg", method = "b", maxdeg = 15, qnts = 0.95)
#'
#' # Partners accumulated over the last 100 steps only
#' plot(dx.hi, type = "cumldeg", window = c(150, 250), momentary = TRUE)
#'
#' # The plotted distributions are returned invisibly
#' dd <- plot(dx.hi, type = "cumldeg")
#' tapply(dd$degree * dd$prop, dd$sim, sum)
#' }
#'
plot.netdx <- function(x, type = "formation", method = "l", sims = NULL,
                       stats = NULL, duration.imputed = TRUE, sim.lines = FALSE,
                       sim.col = NULL, sim.lwd = NULL, mean.line = TRUE,
                       mean.smooth = TRUE, mean.col = NULL, mean.lwd = 2,
                       mean.lty = 1, qnts = 0.5, qnts.col = NULL,
                       qnts.alpha = 0.5, qnts.smooth = TRUE, targ.line = TRUE,
                       targ.col = NULL, targ.lwd = 2, targ.lty = 2,
                       plots.joined = NULL, legend = NULL, grid = FALSE,
                       xlim = NULL, xlab = NULL, ylim = NULL, ylab = NULL,
                       momentary = FALSE, window = NULL, maxdeg = NULL,
                       unique.partners = TRUE, ...) {

  type <- match.arg(type, c("formation", "duration", "dissolution", "cumldeg"))
  sims <- if (is.null(sims)) seq_len(x$nsims) else sims
  if (max(sims) > x$nsims) stop("Maximum sim number is", x$nsims)

  # Cumulative Degree Plot -----------------------------------------------------
  if (type == "cumldeg") {
    dists <- list(cumulative = get_degree_dist(x, type = "cumulative",
                                               sims = sims, window = window,
                                               unique.partners = unique.partners))
    if (isTRUE(momentary)) {
      dists$momentary <- get_degree_dist(x, type = "momentary", sims = sims,
                                         window = window)
    }

    plot_degree_dist(
      dists = dists,
      method = method,
      maxdeg = maxdeg,
      sim.lines = sim.lines,
      sim.col = sim.col,
      sim.lwd = sim.lwd,
      mean.line = mean.line,
      mean.col = mean.col,
      mean.lwd = mean.lwd,
      mean.lty = mean.lty,
      qnts = qnts,
      qnts.col = qnts.col,
      qnts.alpha = qnts.alpha,
      draw_legend = legend,
      grid = grid,
      xlim = xlim, xlab = xlab,
      ylim = ylim, ylab = ylab,
      ...
    )

    out <- do.call(rbind, Map(function(nm, d) cbind(type = nm, d),
                              names(dists), dists))
    row.names(out) <- NULL
    return(invisible(out))
  }

  # Formation Plot -------------------------------------------------------------
  if (type == "formation") {
    stats_table <- x$stats.table.formation
    data <- do.call("cbind", args = x$stats)
    dim3 <- if (x$dynamic) x$nsims else 1L
    data <- array(data, dim = c(dim(data)[1], dim(data)[2] / dim3, dim3))
  } else { # duration/dissolution case -----------------------------------------
    if (!x$dynamic || is.null(x$stats.table.dissolution)) {
      stop(
        "Plots of type duration and dissolution only available if netdx ",
        "run with `dynamic = TRUE` and `skip.dissolution = FALSE`"
      )
    }

    if (x$anyNA) {
      message(
        "\nNOTE: Duration & dissolution data contains undefined values due ",
        "to zero edges of some dissolution dyad type(s) on some time step;",
        " these undefined values will be set to 0 when processing the data."
      )
    }

    if (type == "duration") {
      data <- if (duration.imputed) x$pages_imptd else x$pages
      stats_table <- x$stats.table.duration
      ylab <- NVL(ylab, "Mean Age of Active Ties")
    } else { # if type is "dissolution"
      data <- x$prop.diss
      stats_table <- x$stats.table.dissolution
      ylab <- NVL(ylab, "Proportion of Ties Dissolved")
    }
  }

  sel <- validate_stats_selection(stats_table, data, stats, "netdx object")
  data <- sel$data
  nmstats <- sel$nmstats
  targets <- sel$targets
  # sims only used to subset data in dynamic case
  if (x$dynamic) data <- data[, , sims, drop = FALSE]

  plot_stats_table(
    data = data,
    nmstats = nmstats,
    method = method,
    sim.lines = sim.lines,
    sim.col = sim.col,
    sim.lwd = sim.lwd,
    mean.line = mean.line,
    mean.smooth = mean.smooth,
    mean.col = mean.col,
    mean.lwd = mean.lwd,
    mean.lty = mean.lty,
    qnts = qnts,
    qnts.col = qnts.col,
    qnts.alpha = qnts.alpha,
    qnts.smooth = qnts.smooth,
    targ.line = targ.line,
    targ.col = targ.col,
    targ.lwd = targ.lwd,
    targ.lty = targ.lty,
    plots.joined = plots.joined,
    draw_legend = legend,
    grid = grid,
    targets = targets,
    dynamic = x$dynamic,
    xlim = xlim, xlab = xlab,
    ylim = ylim, ylab = ylab,
    ...
  )
}
