#' @title Plot Dynamic Network Model Diagnostics
#'
#' @description Plots dynamic network model diagnostics calculated in
#'              \code{\link{netdx}}.
#'
#' @param x An \code{EpiModel} object of class \code{netdx}.
#' @param type Plot type, with options of \code{"formation"} for network
#'        model formation statistics, \code{"duration"} for dissolution model
#'        statistics for average edge duration, or \code{"dissolution"} for
#'        dissolution model statistics for proportion of ties dissolved per time
#'        step.
#' @param method Plot method, with options of \code{"l"} for line plots and
#'        \code{"b"} for box plots.
#' @param sims A vector of simulation numbers to plot.
#' @param stats Statistics to plot. For \code{type = "formation"}, \code{stats}
#'        are among those specified in the call to \code{\link{netdx}};
#'        for \code{type = "duration", "dissolution"}, \code{stats} are among
#'        those of the dissolution model (without \code{offset()}). The default
#'        is to plot all statistics.
#' @param plots.joined If \code{TRUE}, combine all statistics in one plot,
#'        versus one plot per statistic if \code{FALSE}.
#' @inheritParams plot.netsim
#'
#' @details
#' The plot function for \code{netdx} objects will generate plots of two types
#' of model diagnostic statistics that run as part of the diagnostic tools
#' within that function. The \code{formation} plot shows the summary statistics
#' requested in \code{nwstats.formula}, where the default includes those
#' statistics in the network model formation formula specified in the original
#' call to \code{\link{netest}}.
#'
#' The \code{duration} plot shows the average age of existing edges at each time
#' step, up until the maximum time step requested. The age is used as an
#' estimator of the average duration of edges in the equilibrium state. When
#' \code{duration.imputed = FALSE}, edges that exist at the beginning of the
#' simulation are assumed to start with an age of 1, yielding a burn-in period
#' before the observed mean approaches its target.  When
#' \code{duration.imputed = TRUE}, expected ages prior to the start of the
#' simulation are calculated from the dissolution model, typically eliminating
#' the need for a burn-in period.
#'
#' The \code{dissolution} plot shows the proportion of the extant ties that are
#' dissolved at each time step, up until the maximum time step requested.
#' Typically, the proportion of ties that are dissolved is the reciprocal of the
#' mean relational duration. This plot thus contains similar information to that
#' in the duration plot, but should reach its expected value more quickly, since
#' it is not subject to censoring.
#'
#' The \code{plots.joined} argument will control whether the statistics
#' are joined in one plot or plotted separately, assuming there are multiple
#' statistics in the model. The default is based on the number of network
#' statistics requested. The layout of the separate plots within the larger plot
#' window is also based on the number of statistics.
#'
#' @method plot netdx
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{netdx}}
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
#' }
#'
plot.netdx <- function(x, type = "formation", method = "l", sims, stats,
                       duration.imputed = TRUE, sim.lines = FALSE, sim.col, sim.lwd,
                       mean.line = TRUE, mean.smooth = TRUE, mean.col,
                       mean.lwd = 2, mean.lty = 1, qnts = 0.5, qnts.col,
                       qnts.alpha = 0.5, qnts.smooth = TRUE, targ.line = TRUE,
                       targ.col, targ.lwd = 2, targ.lty = 2,
                       plots.joined, legend, grid = FALSE, ...) {

  # Checks and Variables ----------------------------------------------------

  ## Check Object
  if (!inherits(x, "netdx")) {
    stop("x must be an object of class netdx", call. = FALSE)
  }

  if (x$dynamic == FALSE && type %in% c("duration", "dissolution")) {
    stop("Plots of type duration and dissolution only available if netdx ",
         "run with dynamic = TRUE", call. = FALSE)
  }

  if (is.null(x$stats.table.dissolution) && type %in% c("duration",
                                                        "dissolution")) {
    stop("Plots of type duration and dissolution only available if netdx ",
         "run with skip.dissolution = FALSE", call. = FALSE)
  }

  ## Check sims
  nsims <- x$nsims
  if (missing(sims)) {
    sims <- 1:nsims
  }
  if (max(sims) > nsims) {
    stop("Maximum sim number is", nsims, call. = FALSE)
  }
  dynamic <- x$dynamic

  # Get dotargs
  da <- list(...)

  type <- match.arg(type, c("formation", "duration", "dissolution"))

  # Formation Plot ----------------------------------------------------------
  if (type == "formation") {
    stats_table <- x$stats.table.formation

    data <- do.call("cbind", args = x$stats)
    dim3 <- if (isTRUE(dynamic)) nsims else 1L
    data <- array(data, dim = c(dim(data)[1], dim(data)[2] / dim3, dim3))
  } else { # duration/dissolution case
    if (x$anyNA == TRUE) {
      cat("\nNOTE: Duration & dissolution data contains undefined values due to zero edges of some dissolution
            dyad type(s) on some time step; these undefined values will be set to 0 when processing the data.")
    }

    if (type == "duration") {
      if (is.logical(duration.imputed) == FALSE) {
        stop("For plots of type duration, duration.imputed must
             be a logical value (TRUE/FALSE)", call. = FALSE)
      }

      if (isTRUE(duration.imputed)) {
        data <- x$pages_imptd
      } else {
        data <- x$pages
      }

      stats_table <- x$stats.table.duration
    } else { # if type is "dissolution"
      data <- x$prop.diss
      stats_table <- x$stats.table.dissolution
    }
  }

  ## Find available stats
  sts <- which(!is.na(stats_table[, "Sim Mean"]))
  nmstats <- rownames(stats_table)[sts]

  ## Pull and check stat argument
  if (missing(stats)) {
    stats <- nmstats
  }
  if (any(stats %in% nmstats == FALSE)) {
    stop("One or more requested stats not contained in netdx object",
         call. = FALSE)
  }
  outsts <- which(nmstats %in% stats)
  nmstats <- nmstats[outsts]

  ## Subset data
  data <- data[, outsts, , drop = FALSE]
  if (isTRUE(dynamic)) {
    # sims only used to subset data in dynamic case
    data <- data[, , sims, drop = FALSE]
  }
  ## Pull target stats
  targets <- stats_table$Target[sts][outsts]

  plot_stats_table(data = data,
                   nmstats = nmstats,
                   method = method,
                   duration.imputed = duration.imputed,
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
                   dynamic = dynamic,
                   da = da,
                   ...)
}


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
plot.netsim <- function(x, type = "epi", y, popfrac = FALSE, sim.lines = FALSE,
                        sims, sim.col, sim.lwd, sim.alpha, mean.line = TRUE,
                        mean.smooth = TRUE, mean.col, mean.lwd = 2,
                        mean.lty = 1, qnts = 0.5, qnts.col, qnts.alpha = 0.5,
                        qnts.smooth = TRUE, legend, leg.cex = 0.8,
                        grid = FALSE, add = FALSE, network = 1, at = 1,
                        col.status = FALSE, shp.g2 = NULL, vertex.cex, stats,
                        targ.line = TRUE, targ.col, targ.lwd = 2, targ.lty = 2,
                        plots.joined, duration.imputed = TRUE, method = "l", ...) {

  type <- match.arg(type, c("epi", "network", "formation", "duration", "dissolution"))

  if (type == "network") {
    # Network plot ------------------------------------------------------------

    if (x$control$tergmLite == TRUE) {
      stop("networkDyanmic object is not saved in tergmLite netsim simulation.
            Check control setting tergmLite", call. = FALSE)
    }

    nsteps <- x$control$nsteps
    if (at > x$control$nsteps) {
      stop("Specify a time step between 1 and ", nsteps, call. = FALSE)
    }

    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1
    }
    if (length(sims) > 1 || (!is.numeric(sims) &&
                               !(sims %in% c("mean", "max", "min")))) {
      stop("sims argument must be single simulation number",
           "or \"mean\", \"max\", or \"min\" ", call. = FALSE)
    }

    sims.arg <- sims
    if (sims == "mean") {
      sims <- which.min(abs(as.numeric(x$epi$i.num[at, ]) -
                              mean(as.numeric(x$epi$i.num[at, ]))))
      sims.val <- as.numeric(x$epi$i.num[at, sims])
    }
    if (sims == "max") {
      sims <- as.numeric(which.max(x$epi$i.num[at, ]))
      sims.val <- x$epi$i.num[at, sims]
    }
    if (sims == "min") {
      sims <- as.numeric(which.min(x$epi$i.num[at, ]))
      sims.val <- x$epi$i.num[at, sims]
    }

    obj <- get_network(x, sims, network, collapse = TRUE, at = at)
    tergmLite <- x$control$tergmLite

    miss_vertex.cex <- missing(vertex.cex)

    if (!is.null(shp.g2)) {
      if (all(shp.g2 != c("square", "triangle"))) {
        stop("shp.g2 accepts inputs of either \"square\" or \"triangle\" ",
             call. = FALSE)
      }

      grp.flag <- length(unique(get_vertex_attribute(obj, "group")))
      if (is.numeric(grp.flag)) {
        mids <- idgroup(obj)
        if (shp.g2 == "square") {
          vertex.sides <- ifelse(mids == 1, 50, 4)
          vertex.rot <- 45
          if (miss_vertex.cex == TRUE) {
            vertex.cex <- 1
          }
        }
        if (shp.g2 == "triangle") {
          vertex.sides <- ifelse(mids == 1, 50, 3)
          vertex.rot <- 90
          if (miss_vertex.cex == TRUE) {
            vertex.cex <- 1
          }
        }

      } else {
        warning("shp.g2 applies to two-group networks only, so ignoring.")
        vertex.sides <- 50
        vertex.rot <- 0
        if (miss_vertex.cex == TRUE) {
          vertex.cex <- 1
        }
      }
    } else {
      vertex.sides <- 50
      vertex.rot <- 0
      if (miss_vertex.cex == TRUE) {
        vertex.cex <- 1
      }
    }
    if (col.status == TRUE) {
      if (tergmLite == TRUE) {
        stop("Plotting status colors requires tergmLite=FALSE in netsim
             control settings.", call. = FALSE)
      }
      pal <- adjustcolor(c(4, 2, 3), 0.75)
      if (tergmLite == FALSE) {
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

  } else if (type == "epi") {
    # Epidemic plot -----------------------------------------------------------

    ## Model dimensions and class ##
    nsteps <- x$control$nsteps
    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1:nsims
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

    # dotargs
    da <- list(...)

    ## Compartments ##
    nocomp <- ifelse(missing(y), TRUE, FALSE)
    if (nocomp == TRUE) {
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
      if (missing(legend)) {
        legend <- TRUE
      }
    }
    if (nocomp == FALSE) {
      if (any(y %in% names(x$epi) == FALSE)) {
        stop("Specified y is not available in object", call. = FALSE)
      }
    }
    lcomp <- length(y)


    ## Color palettes ##

    # Main color palette
    bpal <- c(4, 2, 3, 5:100)

    # Mean line
    if (missing(mean.col)) {
      mean.col <- bpal
    }
    mean.pal <- adjustcolor(mean.col, 0.9)

    # Quantile bands
    if (missing(qnts.col)) {
      qnts.col <- bpal
    }
    qnts.pal <- adjustcolor(qnts.col, qnts.alpha)

    # Sim lines
    if (missing(sim.lwd)) {
      sim.lwd <- rep(0.75, lcomp)
    } else {
      if (length(sim.lwd) < lcomp) {
        sim.lwd <- rep(sim.lwd, lcomp)
      }
    }

    if (missing(sim.col)) {
      sim.col <- bpal
    } else {
      if (length(sim.col) < lcomp) {
        sim.col <- rep(sim.col, lcomp)
      }
    }

    if (missing(sim.alpha) && nsims == 1) {
      sim.alpha <- 0.9
    }
    if (missing(sim.alpha) && nsims > 1) {
      sim.alpha <- max(c(0.05, 1 - log10(nsims) / 3))
    }
    sim.pal <- adjustcolor(sim.col, sim.alpha)


    ## Prevalence calculations ##
    nopopfrac <- ifelse(missing(popfrac), TRUE, FALSE)
    if (nopopfrac == TRUE) {
      popfrac <- FALSE
    }
    if (nopopfrac == TRUE) {
      if (any(grepl(".flow", y)) ||
            (groups == 1 && all(grepl(".num$", y)) == FALSE) ||
            (groups == 2 && all(c(grepl(".num$", y), grepl(".g2$", y)) == FALSE)) ||
            any(y %in% c("num", "num.g2", "num.g2"))) {
        popfrac <- FALSE
      }
    }
    x <- denom(x, y, popfrac)

    # Compartment max
    if (popfrac == FALSE) {
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

      if (!missing(mean.lwd) && length(mean.lwd) < lcomp) {
        mean.lwd <- rep(mean.lwd, lcomp)
      }
      if (missing(mean.lwd)) {
        mean.lwd <- rep(1.5, lcomp)
      }

      if (!missing(mean.lty) && length(mean.lty) < lcomp) {
        mean.lty <- rep(mean.lty, lcomp)
      }
      if (missing(mean.lty)) {
        mean.lty <- rep(1, lcomp)
      }
      mean.max <- draw_means(x, y, mean.smooth, mean.lwd, mean.pal,
                             mean.lty, "epi", 0, "max")
      mean.min <- draw_means(x, y, mean.smooth, mean.lwd, mean.pal,
                             mean.lty, "epi", 0, "min")
    }

    ## Missing args ##
    if (is.null(da$xlim)) {
      da$xlim <- c(0, nsteps)
    }
    if (is.null(da$ylim) && (popfrac == TRUE || sim.lines == TRUE)) {
      da$ylim <- c(min.prev, max.prev)
    } else if (is.null(da$ylim) && popfrac == FALSE && sim.lines == FALSE &&
                 (mean.line == TRUE || qnts == TRUE)) {
      da$ylim <- c(min(qnt.min * 0.9, mean.min * 0.9), max(qnt.max * 1.1, mean.max * 1.1))
    }

    if (is.null(da$main)) {
      da$main <- ""
    }

    if (is.null(da$xlab)) {
      da$xlab <- "Time"
    }

    if (is.null(da$ylab)) {
      if (popfrac == FALSE) {
        da$ylab <- "Number"
      } else {
        da$ylab <- "Prevalence"
      }
    }

    ## Main plot window ##
    if (add == FALSE) {
      da$x <- 1
      da$y <- 1
      da$type <- "n"
      da$bty <- "n"

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

      if (!missing(mean.lwd) && length(mean.lwd) < lcomp) {
        mean.lwd <- rep(mean.lwd, lcomp)
      }
      if (missing(mean.lwd)) {
        mean.lwd <- rep(2.5, lcomp)
      }

      if (!missing(mean.lty) && length(mean.lty) < lcomp) {
        mean.lty <- rep(mean.lty, lcomp)
      }
      if (missing(mean.lty)) {
        if (nocomp == FALSE) {
          mean.lty <- rep(1, lcomp)
        }
      }
      y.n <- length(y)
      mean.pal <- mean.pal[1:y.n]
      draw_means(x, y, mean.smooth, mean.lwd, mean.pal, mean.lty)
    }

    ## Grid
    if (grid == TRUE) {
      grid()
    }

    ## Legends ##
    if (!missing(legend) && legend == TRUE) {
      if (groups == 2 && nocomp == TRUE) {
        leg.lty <- mean.lty
      } else {
        leg.lty <- 1
      }
      legend("topright", legend = y, lty = leg.lty, lwd = 2,
             col = mean.pal, cex = leg.cex, bg = "white")
    }
  } else {
    # stat plot

    ## Stats
    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1:nsims
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
    if (missing(stats)) {
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

    da <- list(...)

    plot_stats_table(data = data,
                     nmstats = nmstats,
                     method = method,
                     duration.imputed = duration.imputed,
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
                     da = da,
                     ...)
  }
}


#' @method comp_plot netsim
#' @rdname comp_plot
#' @export
comp_plot.netsim <- function(x, at = 1, digits = 3, ...) {

  comp_plot.icm(x = x, at = at, digits = digits, ...)

}

