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
#' @inheritParams graphics::plot
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
plot.netdx <- function(x, type = "formation", method = "l", sims = NULL,
                       stats = NULL, duration.imputed = TRUE, sim.lines = FALSE,
                       sim.col = NULL, sim.lwd = NULL, mean.line = TRUE,
                       mean.smooth = TRUE, mean.col = NULL, mean.lwd = 2,
                       mean.lty = 1, qnts = 0.5, qnts.col = NULL,
                       qnts.alpha = 0.5, qnts.smooth = TRUE, targ.line = TRUE,
                       targ.col = NULL, targ.lwd = 2, targ.lty = 2,
                       plots.joined = NULL, legend = NULL, grid = FALSE, ...) {

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
  if (is.null(sims)) {
    sims <- seq_len(nsims)
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
  if (is.null(stats)) {
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
