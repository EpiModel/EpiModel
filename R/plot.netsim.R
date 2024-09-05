#' @title Plot Data from a Stochastic Network Epidemic Model
#'
#' @description Plots epidemiological and network data from a stochastic network
#'              model simulated with \code{\link{netsim}}.
#'
#' @param x An \code{EpiModel} model object of class \code{netsim}.
#' @param type Type of plot: \code{"epi"} for epidemic model results,
#'        \code{"network"} for a static network plot (\code{plot.network}),
#'        or \code{"formation"}, \code{"duration"}, or \code{"dissolution"} for
#'        network formation, duration, or dissolution statistics.
#' @param y Output compartments or flows from \code{netsim} object to plot.
#' @param popfrac If \code{TRUE}, plot prevalence of values rather than numbers
#'        (see details).
#' @param sim.lines If \code{TRUE}, plot individual simulation lines. Default is
#'        to plot lines for one-group models but not for two-group models.
#' @param sims If \code{type="epi"} or \code{"formation"}, a vector of
#'        simulation numbers to plot. If \code{type="network"}, a single
#'        simulation number for which to plot the network, or else \code{"min"}
#'        to plot the simulation number with the lowest disease prevalence,
#'        \code{"max"} for the simulation with the highest disease prevalence,
#'        or \code{"mean"} for the simulation with the prevalence closest to the
#'        mean across simulations at the specified time step.
#' @param sim.col Vector of any standard R color format for simulation lines.
#' @param sim.lwd Line width for simulation lines.
#' @param sim.alpha Transparency level for simulation lines, where
#'        0 = transparent and 1 = opaque (see \code{adjustcolor} function).
#' @param mean.line If \code{TRUE}, plot mean of simulations across time.
#' @param mean.smooth If \code{TRUE}, use a loess smoother on the mean line.
#' @param mean.col Vector of any standard R color format for mean lines.
#' @param mean.lwd Line width for mean lines.
#' @param mean.lty Line type for mean lines.
#' @param qnts If numeric, plot polygon of simulation quantiles based on the
#'        range implied by the argument (see details). If \code{FALSE}, suppress
#'        polygon from plot.
#' @param qnts.col Vector of any standard R color format for polygons.
#' @param qnts.alpha Transparency level for quantile polygons, where 0 =
#'        transparent and 1 = opaque (see \code{adjustcolor} function).
#' @param qnts.smooth If \code{TRUE}, use a loess smoother on quantile polygons.
#' @param legend If \code{TRUE}, plot default legend.
#' @param leg.cex Legend scale size.
#' @param grid If \code{TRUE}, a grid is added to the background of plot
#'        (see \code{\link{grid}} for details), with default of nx by ny.
#' @param add If \code{TRUE}, new plot window is not called and lines are added
#'        to existing plot window.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#' @param at If \code{type = "network"}, time step for network graph.
#' @param col.status If \code{TRUE} and \code{type="network"}, automatic disease
#'        status colors (blue = susceptible, red = infected, green = recovered).
#' @param shp.g2 If \code{type = "network"} and \code{x} is for a two-group model,
#'        shapes for the Group 2 vertices, with acceptable inputs of "triangle"
#'        and "square". Group 1 vertices will remain circles.
#' @param vertex.cex Relative size of plotted vertices if \code{type="network"},
#'        with implicit default of 1.
#' @param stats If \code{type="formation","duration","dissolution"}, statistics
#'        to plot. For \code{type = "formation"}, \code{stats} are among those
#'        specified in \code{nwstats.formula} of \code{\link{control.net}}; for
#'        \code{type = "duration", "dissolution"}, \code{stats} are among those
#'        of the dissolution model (without \code{offset()}). The default is
#'        to plot all statistics.
#' @param targ.line If \code{TRUE}, plot target or expected value line for
#'        the statistic of interest.
#' @param targ.col Vector of standard R colors for target statistic lines, with
#'        default colors based on \code{RColorBrewer} color palettes.
#' @param targ.lwd Line width for the line showing the target statistic values.
#' @param targ.lty Line type for the line showing the target statistic values.
#' @param plots.joined If \code{TRUE} and
#'        \code{type="formation","duration","dissolution"}, combine all
#'        statistics in one plot, versus one plot per statistic if
#'        \code{FALSE}.
#' @param method Plot method for \code{type="formation", "duration", "dissolution"},
#'        with options of \code{"l"} for line plots and \code{"b"} for box plots.
#' @param duration.imputed If \code{type = "duration"}, a logical indicating
#'        whether or not to impute starting times for relationships extant at
#'        the start of the simulation. Defaults to \code{TRUE} when
#'        \code{type = "duration"}.
#' @param ... Additional arguments to pass.
#' @inheritParams graphics::plot
#'
#' @details
#' This plot function can produce three types of plots with a stochastic network
#' model simulated through \code{\link{netsim}}:
#' \enumerate{
#'  \item \strong{\code{type="epi"}}: epidemic model results (e.g., disease
#'        prevalence and incidence) may be plotted.
#'  \item \strong{\code{type="network"}}: a static network plot will be
#'        generated. A static network plot of a dynamic network is a
#'        cross-sectional extraction of that dynamic network at a specific
#'        time point. This plotting function wraps the
#'        \code{\link{plot.network}} function in the \code{network} package.
#'        Consult the help page for \code{plot.network} for all of the plotting
#'        parameters. In addition, four plotting parameters specific to
#'        \code{netsim} plots are available: \code{sim}, \code{at},
#'        \code{col.status}, and \code{shp.g2}.
#'  \item \strong{\code{type="formation"}}: summary network statistics related
#'        to the network model formation are plotted. These plots are similar
#'        to the formation plots for \code{netdx} objects. When running a
#'        \code{netsim} simulation, one must specify there that
#'        \code{save.nwstats=TRUE}; the plot here will then show the network
#'        statistics requested explicitly in \code{nwstats.formula}, or will use
#'        the formation formula set in \code{netest} otherwise.
#'  \item \strong{\code{type="duration","dissolution"}}: as in
#'        \code{\link{plot.netdx}}; supported in \code{plot.netsim} only when
#'        the dissolution model is \code{~offset(edges)}, \code{tergmLite} is
#'        \code{FALSE}, and \code{save.network} is \code{TRUE}.
#' }
#'
#' @details
#' When \code{type="epi"}, this plotting function will extract the
#' epidemiological output from a model object of class \code{netsim} and plot
#' the time series data of disease prevalence and other results. The summary
#' statistics that the function calculates and plots are individual simulation
#' lines, means of the individual simulation lines, and quantiles of those
#' individual simulation lines. The mean line, toggled on with
#' \code{mean.line=TRUE}, is calculated as the row mean across simulations at
#' each time step.
#'
#' Compartment prevalences are the size of a compartment over some denominator.
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this
#' is the default for any plots of flows. The \code{popfrac} parameter
#' calculates and plots the denominators of all specified compartments using
#' these rules: 1) for one-group models, the prevalence of any compartment is
#' the compartment size divided by the total population size; 2) for two-group
#' models, the prevalence of any compartment is the compartment size divided by
#' the group population size. For any prevalences that are not automatically
#' calculated, the \code{\link{mutate_epi}} function may be used to add new
#' variables to the \code{netsim} object to plot or analyze.
#'
#' The quantiles show the range of outcome values within a certain specified
#' quantile range. By default, the interquartile range is shown: that is the
#' middle 50\% of the data. This is specified by \code{qnts=0.5}. To show the
#' middle 95\% of the data, specify \code{qnts=0.95}. To toggle off the polygons
#' where they are plotted by default, specify \code{qnts=FALSE}.
#'
#' When \code{type="network"}, this function will plot cross sections of the
#' simulated networks at specified time steps. Because it is only possible to
#' plot one time step from one simulation at a time, it is necessary to enter
#' these in the \code{at} and \code{sims} parameters. To aid in visualizing
#' representative and extreme simulations at specific time steps, the
#' \code{sims} parameter may be set to \code{"mean"} to plot the simulation in
#' which the disease prevalence is closest to the average across all
#' simulations, \code{"min"} to plot the simulation in which the prevalence is
#' lowest, and \code{"max"} to plot the simulation in which the prevalence is
#' highest.
#'
#' @method plot netsim
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{plot.network}}, \code{\link{mutate_epi}}
#'
#' @examples
#' ## SI Model without Network Feedback
#' # Initialize network and set network model parameters
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Estimate the network model
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the epidemic model
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 20, nsims = 3,
#'                        verbose = FALSE, save.nwstats = TRUE,
#'                        nwstats.formula = ~edges + meandeg + concurrent)
#' mod <- netsim(est, param, init, control)
#'
#' # Plot epidemic trajectory
#' plot(mod)
#' plot(mod, type = "epi", grid = TRUE)
#' plot(mod, type = "epi", popfrac = TRUE)
#' plot(mod, type = "epi", y = "si.flow", qnts = 1, ylim = c(0, 4))
#'
#' # Plot static networks
#' par(mar = c(0, 0, 0, 0))
#' plot(mod, type = "network", vertex.cex = 1.5)
#'
#' # Automatic coloring of infected nodes as red
#' par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
#' plot(mod, type = "network", main = "Min Prev | Time 50",
#'      col.status = TRUE, at = 20, sims = "min", vertex.cex = 1.25)
#' plot(mod, type = "network", main = "Max Prev | Time 50",
#'      col.status = TRUE, at = 20, sims = "max", vertex.cex = 1.25)
#'
#' # Automatic shape by group number (circle = group 1)
#' par(mar = c(0, 0, 0, 0))
#' plot(mod, type = "network", at = 20, col.status = TRUE,
#'      shp.g2 = "square")
#' plot(mod, type = "network", at = 20, col.status = TRUE,
#'      shp.g2 = "triangle", vertex.cex = 2)
#'
#' # Plot formation statistics
#' par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
#' plot(mod, type = "formation", grid = TRUE)
#' plot(mod, type = "formation", plots.joined = FALSE)
#' plot(mod, type = "formation", sims = 2:3)
#' plot(mod, type = "formation", plots.joined = FALSE,
#'      stats = c("edges", "concurrent"))
#' plot(mod, type = "formation", stats = "meandeg",
#'      mean.lwd = 1, qnts.col = "seagreen", mean.col = "black")
#'
plot.netsim <- function(x, type = "epi", y = NULL, popfrac = FALSE,
                        sim.lines = FALSE, sims = NULL, sim.col = NULL,
                        sim.lwd = NULL, sim.alpha = NULL, mean.line = TRUE,
                        mean.smooth = TRUE, mean.col = NULL, mean.lwd = 2,
                        mean.lty = 1, qnts = 0.5, qnts.col = NULL,
                        qnts.alpha = 0.5, qnts.smooth = TRUE, legend = NULL,
                        leg.cex = 0.8, grid = FALSE, add = FALSE, network = 1,
                        at = 1, col.status = FALSE, shp.g2 = NULL,
                        vertex.cex = NULL, stats = NULL, targ.line = TRUE,
                        targ.col = NULL, targ.lwd = 2, targ.lty = 2,
                        plots.joined = NULL, duration.imputed = TRUE,
                        method = "l", main = NULL, xlim = NULL, xlab = NULL,
                        ylim = NULL, ylab = NULL, ...) {

  type <- match.arg(
    type,
    c("epi", "network", "formation", "duration", "dissolution")
  )

  if (type == "network") {
    plot_netsim_network(x, at, sims, network, shp.g2, col.status, vertex.cex, ...)
  } else if (type == "epi") {
    plot_netsim_epi(x, y, sims, legend, mean.col, qnts.col, sim.lwd,
                            sim.col, sim.alpha, popfrac, qnts, qnts.alpha, qnts.smooth,
                            mean.line, mean.smooth, add,
                            mean.lwd, mean.lty, xlim, ylim, main, xlab, ylab,
                            sim.lines, grid, leg.cex, ...)
  } else {
    plot_netsim_stats(
      x, type, sims, stats, network, duration.imputed,
      method, sim.lines, sim.col, sim.lwd,
      mean.line, mean.smooth, mean.col, mean.lwd,
      mean.lty, qnts, qnts.col, qnts.alpha, qnts.smooth,
      targ.line, targ.col, targ.lwd, targ.lty,
      plots.joined, legend, grid, xlim, xlab,
      ylim, ylab, ...)
  }
}

plot_netsim_network <- function(x, at, sims, network, shp.g2, col.status, vertex.cex, ...) {
  # Network plot ------------------------------------------------------------
  if (x$control$tergmLite) {
    stop("networkDyanmic object is not saved in tergmLite netsim simulation.
      Check control setting tergmLite", call. = FALSE)
  }

  nsteps <- x$control$nsteps
  if (at > nsteps) {
    stop("Specify a time step between 1 and ", nsteps, call. = FALSE)
  }

  nsims <- x$control$nsims
  sims <- if (is.null(sims)) 1 else sims
  if (length(sims) > 1 ||
    (!is.numeric(sims) && !(sims %in% c("mean", "max", "min")))) {
    stop("sims argument must be single simulation number",
      "or \"mean\", \"max\", or \"min\" ", call. = FALSE)
  }

  sims.arg <- sims
  if (sims == "mean") {
    sims <- which.min(
      abs(as.numeric(x$epi$i.num[at, ]) - mean(as.numeric(x$epi$i.num[at, ])))
    )
    sims.val <- as.numeric(x$epi$i.num[at, sims])
  } else if (sims == "max") {
    sims <- as.numeric(which.max(x$epi$i.num[at, ]))
    sims.val <- x$epi$i.num[at, sims]
  } else if (sims == "min") {
    sims <- as.numeric(which.min(x$epi$i.num[at, ]))
    sims.val <- x$epi$i.num[at, sims]
  }

  obj <- get_network(x, sims, network, collapse = TRUE, at = at)
  tergmLite <- x$control$tergmLite
  miss_vertex.cex <- is.null(vertex.cex)

  if (!is.null(shp.g2)) {
    if (all(!shp.g2 %in% c("square", "triangle"))) {
      stop("shp.g2 accepts inputs of either \"square\" or \"triangle\" ",
        call. = FALSE)
    }

    grp.flag <- length(unique(get_vertex_attribute(obj, "group")))
    if (is.numeric(grp.flag)) {
      mids <- idgroup(obj)
      if (shp.g2 == "square") {
        vertex.sides <- ifelse(mids == 1, 50, 4)
        vertex.rot <- 45
        if (miss_vertex.cex) {
          vertex.cex <- 1
        }
      }
      if (shp.g2 == "triangle") {
        vertex.sides <- ifelse(mids == 1, 50, 3)
        vertex.rot <- 90
        if (miss_vertex.cex) {
          vertex.cex <- 1
        }
      }

    } else {
      warning("shp.g2 applies to two-group networks only, so ignoring.")
      vertex.sides <- 50
      vertex.rot <- 0
      if (miss_vertex.cex) {
        vertex.cex <- 1
      }
    }
  } else {
    vertex.sides <- 50
    vertex.rot <- 0
    if (miss_vertex.cex) {
      vertex.cex <- 1
    }
  }
  if (col.status) {
    if (tergmLite) {
      stop("Plotting status colors requires tergmLite=FALSE in netsim
        control settings.", call. = FALSE)
    }
    pal <- adjustcolor(c(4, 2, 3), 0.75)
    if (!tergmLite) {
      testatus <- get.vertex.attribute.active(obj, "testatus", at = at)
      cols <- rep(pal[1], length(testatus))
      cols[testatus == "i"] <- pal[2]
      cols[testatus == "r"] <- pal[3]
    }
    plot.network(obj, vertex.col = cols, vertex.border = "grey60",
      edge.col = "grey40", vertex.sides = vertex.sides,
      vertex.rot = vertex.rot, vertex.cex = vertex.cex,
      displaylabels = FALSE, ...)
    if (sims.arg %in% c("mean", "max", "min")) {
      mtext(side = 1, text = paste("Sim =", sims, " | Prev =", sims.val))
    }
  } else {
    plot.network(obj, vertex.sides = vertex.sides, vertex.rot = vertex.rot,
      vertex.cex = vertex.cex, displaylabels = FALSE, ...)
  }
}

plot_netsim_epi <- function(x, y, sims, legend, mean.col, qnts.col, sim.lwd,
                            sim.col, sim.alpha, popfrac, qnts, qnts.alpha, qnts.smooth,
                            mean.line, mean.smooth, add,
                            mean.lwd, mean.lty, xlim, ylim, main, xlab, ylab,
                            sim.lines, grid, leg.cex, ...) {
  ## Model dimensions and class ##
  nsteps <- x$control$nsteps
  nsims <- x$control$nsims
  if (is.null(sims)) {
    sims <- seq_len(nsims)
  }
  if (max(sims) > nsims) {
    stop("Set sim to between 1 and ", nsims, call. = FALSE)
  }
  if (is.null(x$param$groups) || !is.numeric(x$param$groups)) {
    groups <- 1
    x$param$groups <- 1
  } else {
    groups <- x$param$groups
  }

  ## Compartments ##
  nocomp <- is.null(y)
  if (nocomp) {
    if (groups == 1) {
      y <- grep(".num$", names(x$epi), value = TRUE)
    }
    if (groups == 2) {
      if (inherits(x, "icm")) {
        y <- c(grep(".num$", names(x$epi), value = TRUE),
          grep(".num.g2$", names(x$epi), value = TRUE))
      }
      if (inherits(x, "netsim")) {
        y <- c(grep(".num$", names(x$epi), value = TRUE),
          grep(".num.g2$", names(x$epi), value = TRUE))
      }
    }
    if (is.null(legend)) {
      legend <- TRUE
    }
  }
  if (!nocomp) {
    if (any(y %in% names(x$epi) == FALSE)) {
      stop("Specified y is not available in object", call. = FALSE)
    }
  }
  lcomp <- length(y)


  ## Color palettes ##

  # Main color palette
  bpal <- c(4, 2, 3, 5:100)

  # Mean line
  if (is.null(mean.col)) {
    mean.col <- bpal
  }
  mean.pal <- adjustcolor(mean.col, 0.9)

  # Quantile bands
  if (is.null(qnts.col)) {
    qnts.col <- bpal
  }
  qnts.pal <- adjustcolor(qnts.col, qnts.alpha)

  # Sim lines
  if (is.null(sim.lwd)) {
    sim.lwd <- rep(0.75, lcomp)
  } else {
    if (length(sim.lwd) < lcomp) {
      sim.lwd <- rep(sim.lwd, lcomp)
    }
  }

  if (is.null(sim.col)) {
    sim.col <- bpal
  } else {
    if (length(sim.col) < lcomp) {
      sim.col <- rep(sim.col, lcomp)
    }
  }

  if (is.null(sim.alpha) && nsims == 1) {
    sim.alpha <- 0.9
  }
  if (is.null(sim.alpha) && nsims > 1) {
    sim.alpha <- max(c(0.05, 1 - log10(nsims) / 3))
  }
  sim.pal <- adjustcolor(sim.col, sim.alpha)


  ## Prevalence calculations ##
  x <- denom(x, y, popfrac)

  # Compartment max
  if (!popfrac) {
    if (lcomp == 1) {
      min.prev <- min(x$epi[[y]], na.rm = TRUE)
      max.prev <- max(x$epi[[y]], na.rm = TRUE)
    } else {
      min.prev <- min(sapply(y, function(comps) min(x$epi[[comps]], na.rm = TRUE)))
      max.prev <- max(sapply(y, function(comps) max(x$epi[[comps]], na.rm = TRUE)))
    }
  } else {
    min.prev <- 0
    max.prev <- 1
  }

  # Initialize ylim max values
  qnt.min <- 1E10
  qnt.max <- -1E10
  mean.min <- 1E10
  mean.max <- -1E10

  ## Quantiles - ylim max ##
  if (qnts == FALSE) {
    disp.qnts <- FALSE
  } else {
    disp.qnts <- TRUE
  }
  if (nsims == 1) {
    disp.qnts <- FALSE
  }

  if (disp.qnts == TRUE) {
    if (qnts > 1 || qnts < 0) {
      stop("qnts must be between 0 and 1", call. = FALSE)
    }
    qnt.max <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "max")
    qnt.min <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "min")
  }


  ## Mean lines - ylim max ##
  if (mean.line == TRUE) {

    if (!is.null(mean.lwd) && length(mean.lwd) < lcomp) {
      mean.lwd <- rep(mean.lwd, lcomp)
    }
    if (is.null(mean.lwd)) {
      mean.lwd <- rep(1.5, lcomp)
    }

    if (!is.null(mean.lty) && length(mean.lty) < lcomp) {
      mean.lty <- rep(mean.lty, lcomp)
    }
    if (is.null(mean.lty)) {
      mean.lty <- rep(1, lcomp)
    }
    mean.max <- draw_means(x, y, mean.smooth, mean.lwd, mean.pal,
      mean.lty, "epi", 0, "max")
    mean.min <- draw_means(x, y, mean.smooth, mean.lwd, mean.pal,
      mean.lty, "epi", 0, "min")
  }

  ## Missing args ##
  if (is.null(xlim)) {
    xlim <- c(0, nsteps)
  }
  if (is.null(ylim) && (popfrac == TRUE || sim.lines == TRUE)) {
    ylim <- c(min.prev, max.prev)
  } else if (is.null(ylim) && popfrac == FALSE && sim.lines == FALSE &&
    (mean.line == TRUE || qnts == TRUE)) {
    ylim <- c(min(qnt.min * 0.9, mean.min * 0.9), max(qnt.max * 1.1, mean.max * 1.1))
  }

  if (is.null(main)) {
    main <- ""
  }

  if (is.null(xlab)) {
    xlab <- "Time"
  }

  if (is.null(ylab)) {
    if (popfrac == FALSE) {
      ylab <- "Number"
    } else {
      ylab <- "Prevalence"
    }
  }

  ## Main plot window ##
  if (add == FALSE) {
    da <- list()
    da$x <- 1
    da$y <- 1
    da$type <- "n"
    da$bty <- "n"
    da$xlim <- xlim
    da$xlab <- xlab
    da$ylim <- ylim
    da$ylab <- ylab
    da$main <- main

    do.call(plot, da)
  }


  ## Quantiles ##
  ## NOTE: Why is this repeated from above?
  if (qnts == FALSE) {
    disp.qnts <- FALSE
  } else {
    disp.qnts <- TRUE
  }
  if (nsims == 1) {
    disp.qnts <- FALSE
  }

  if (disp.qnts == TRUE) {
    if (qnts > 1 || qnts < 0) {
      stop("qnts must be between 0 and 1", call. = FALSE)
    }
    y.l <- length(y)
    qnts.pal <- qnts.pal[1:y.l]
    draw_qnts(x, y, qnts, qnts.pal, qnts.smooth)
  }


  ## Simulation lines ##
  if (sim.lines == TRUE) {
    for (j in seq_len(lcomp)) {
      for (i in sims) {
        lines(x$epi[[y[j]]][, i], lwd = sim.lwd[j], col = sim.pal[j])
      }
    }
  }


  ## Mean lines ##
  if (mean.line == TRUE) {

    if (!is.null(mean.lwd) && length(mean.lwd) < lcomp) {
      mean.lwd <- rep(mean.lwd, lcomp)
    }
    if (is.null(mean.lwd)) {
      mean.lwd <- rep(2.5, lcomp)
    }

    if (!is.null(mean.lty) && length(mean.lty) < lcomp) {
      mean.lty <- rep(mean.lty, lcomp)
    }
    if (is.null(mean.lty)) {
      if (nocomp == FALSE) {
        mean.lty <- rep(1, lcomp)
      }
    }
    y.n <- length(y)
    mean.pal <- mean.pal[1:y.n]
    draw_means(x, y, mean.smooth, mean.lwd, mean.pal, mean.lty)
  }

  ## Grid
  if (grid) grid()

  ## Legends ##
  if (!is.null(legend) && legend) {
    if (groups == 2 && nocomp) {
      leg.lty <- mean.lty
    } else {
      leg.lty <- 1
    }
    legend("topright", legend = y, lty = leg.lty, lwd = 2,
      col = mean.pal, cex = leg.cex, bg = "white")
  }
}

plot_netsim_stats <- function(x, type, sims, stats, network, duration.imputed,
                              method, sim.lines, sim.col, sim.lwd,
                              mean.line, mean.smooth, mean.col, mean.lwd,
                              mean.lty, qnts, qnts.col, qnts.alpha, qnts.smooth,
                              targ.line, targ.col, targ.lwd, targ.lty,
                              plots.joined, legend, grid, xlim, xlab,
                              ylim, ylab, ...) {

  nsims <- x$control$nsims
  if (is.null(sims)) {
    sims <- seq_len(nsims)
  }
  if (max(sims) > nsims) {
    stop("Maximum sims for this object is ", nsims, call. = FALSE)
  }

  nsims <- length(sims)
  nsteps <- x$control$nsteps

  if (type == "formation") {
    # Formation plot ----------------------------------------------------------

    ## get nw stats
    data <- get_nwstats(x, sims, network, mode = "list")

    ## target stats
    nwparam <- get_nwparam(x, network)
    ts <- nwparam$target.stats
    tsn <- nwparam$target.stats.names
    names(ts) <- tsn

  } else {
    ## duration/dissolution plot
    if (isTRUE(x$control$save.diss.stats) &&
      isTRUE(x$control$save.network) &&
      isFALSE(x$control$tergmLite) &&
      isFALSE(is.null(x$diss.stats)) &&
      isTRUE(x$nwparam[[network]]$coef.diss$diss.model.type == "edgesonly")) {

      if (any(unlist(lapply(x$diss.stats, `[[`, "anyNA")))) {
        cat("\nNOTE: Duration & dissolution data contains undefined values due to zero edges of some dissolution
          dyad type(s) on some time step; these undefined values will be set to 0 when processing the data.")
      }

      if (type == "duration") {
        if (isTRUE(duration.imputed)) {
          data <- lapply(sims, function(sim) x$diss.stats[[sim]][[network]][["meanageimputed"]])
        } else {
          data <- lapply(sims, function(sim) x$diss.stats[[sim]][[network]][["meanage"]])
        }
        ts <- x$nwparam[[network]]$coef.diss$duration
      } else { # if type is "dissolution"
        data <- lapply(sims, function(sim) x$diss.stats[[sim]][[network]][["propdiss"]])
        ts <- 1 / x$nwparam[[network]]$coef.diss$duration
      }
    } else {
      stop("cannot produce duration/dissolution plot from `netsim` object ",
        "unless `save.diss.stats` is `TRUE`, `save.network` is `TRUE`, ",
        "`tergmLite` is `FALSE`, `keep.diss.stats` is `TRUE` (if ",
        "merging), and dissolution model is edges-only")
    }
  }

  stats_table <- make_stats_table(data, ts)
  data <- array(unlist(data), dim = c(dim(data[[1]]), nsims))

  ## Find available stats
  sts <- which(!is.na(stats_table[, "Sim Mean"]))
  nmstats <- rownames(stats_table)[sts]

  ## Pull and check stat argument
  if (is.null(stats)) {
    stats <- nmstats
  }
  if (any(stats %in% nmstats == FALSE)) {
    stop("One or more requested stats not contained in netsim object",
      call. = FALSE)
  }
  outsts <- which(nmstats %in% stats)
  nmstats <- nmstats[outsts]

  ## Subset data
  data <- data[, outsts, , drop = FALSE]

  ## we've already subset the data to `sims`

  ## Pull target stats
  targets <- stats_table$Target[sts][outsts]

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
    dynamic = TRUE, # always dynamic in netsim
    xlim = xlim, xlab = xlab,
    ylim = ylim, ylab = ylab,
    ...
  )
}
