#' @title Plot Data from a Stochastic Individual Contact Epidemic Model
#'
#' @description Plots epidemiological data from a stochastic individual contact
#'              model simulated with \code{\link{icm}}.
#'
#' @param x An \code{EpiModel} model object of class \code{icm}.
#' @param y Output compartments or flows from \code{icm} object to plot. -------
#' @param sims A vector of simulation numbers to plot.
#' @inheritParams plot.netsim
#' @inheritParams graphics::plot
#'
#' @details
#' This plotting function will extract the epidemiological output from a model
#' object of class \code{icm} and plot the time series data of disease
#' prevalence and other results. The summary statistics that the function
#' calculates and plots are individual simulation lines, means of the individual
#' simulation lines, and quantiles of those individual simulation lines. The
#' mean line, toggled on with \code{mean.line=TRUE}, is calculated as the row
#' mean across simulations at each time step.
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
#' variables to the \code{icm} object to plot or analyze.
#'
#' The quantiles show the range of outcome values within a certain specified
#' quantile range. By default, the interquartile range is shown: that is the
#' middle 50\% of the data. This is specified by \code{qnts=0.5}. To show the
#' middle 95\% of the data, specify \code{qnts=0.95}. To toggle off the polygons
#' where they are plotted by default, specify \code{qnts=FALSE}.
#'
#' @method plot icm
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{icm}}
#'
#' @examples
#' ## Example 1: Plotting multiple compartment values from SIR model
#' param <- param.icm(inf.prob = 0.5, act.rate = 0.5, rec.rate = 0.02)
#' init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.icm(type = "SIR", nsteps = 100,
#'                        nsims = 3, verbose = FALSE)
#' mod <- icm(param, init, control)
#' plot(mod, grid = TRUE)
#'
#' ## Example 2: Plot only infected with specific output from SI model
#' param <- param.icm(inf.prob = 0.25, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 10)
#' control <- control.icm(type = "SI", nsteps = 100,
#'                        nsims = 3, verbose = FALSE)
#' mod2 <- icm(param, init, control)
#'
#' # Plot prevalence
#' plot(mod2, y = "i.num", mean.line = FALSE, sim.lines = TRUE)
#'
#' # Plot incidence
#' par(mfrow = c(1, 2))
#' plot(mod2, y = "si.flow", mean.smooth = TRUE, grid = TRUE)
#' plot(mod2, y = "si.flow", qnts.smooth = FALSE, qnts = 1)
#'
plot.icm <- function(x, y = NULL, popfrac = FALSE, sim.lines = FALSE,
                     sims = NULL, sim.col = NULL, sim.lwd = NULL,
                     sim.alpha = NULL, mean.line = TRUE, mean.smooth = TRUE,
                     mean.col = NULL, mean.lwd = 2, mean.lty = 1, qnts = 0.5,
                     qnts.col = NULL, qnts.alpha = 0.5, qnts.smooth = TRUE,
                     legend = TRUE, leg.cex = 0.8, grid = FALSE, add = FALSE,
                     xlim = NULL, ylim = NULL, main = "", xlab = "Time",
                     ylab = NULL,
                     ...) {

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
  if (nocomp == TRUE) {
    if (groups == 1) {
      y <- grep(".num$", names(x$epi), value = TRUE)
    }
    if (groups == 2) {
      if (inherits(x, "icm")) {
        y <- c(grep(".num$", names(x$epi), value = TRUE),
               grep(".num.g2$", names(x$epi), value = TRUE))
      }
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
  bpal <- c(4, 2, 3)

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

  if (is.null(sim.alpha)) {
    sim.alpha <- 1 - log10(nsims) / 3
    sim.alpha <- dplyr::between(sim.alpha, 0.05, 0.9)
  }
  sim.pal <- adjustcolor(sim.col, sim.alpha)

  ## Prevalence calculations ##
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


  ## Missing args ##
  if (is.null(xlim)) {
    xlim <- c(0, nsteps)
  }

  #Initialize ylim min max values
  qnt.min <- 1E10
  qnt.max <- -1E10
  mean.min <- 1E10
  mean.max <- -1E10

  ## Quantiles - ylim min max ##
  if (nsims > 1) {
    if (qnts > 1 || qnts < 0) {
      stop("qnts must be between 0 and 1", call. = FALSE)
    }
    qnt.min <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "min")
    qnt.max <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "max")
  }

  ## Mean lines - ylim max ##
  if (mean.line == TRUE) {

    if (length(mean.lwd) < lcomp) {
      mean.lwd <- rep(mean.lwd, lcomp)
    }

    if (length(mean.lty) < lcomp) {
      mean.lty <- rep(mean.lty, lcomp)
    }
    mean.min <- draw_means(x, y, mean.smooth, mean.lwd,
                           mean.pal, mean.lty, "epi", 0, "min")
    mean.max <- draw_means(x, y, mean.smooth, mean.lwd,
                           mean.pal, mean.lty, "epi", 0, "max")
  }

  # Dynamic scaling based on sim.lines and mean lines and quantile bands
  if (is.null(ylim)) {
    if (sim.lines == FALSE && mean.line == TRUE) {
      ylim <- c(min(qnt.min * 0.9, mean.min * 0.9),
                max(qnt.max * 1.1, mean.max * 1.1))
    } else {
      ylim <- c(min.prev, max.prev)
    }
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
    plot(1, 1, type = "n", bty = "n",
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ...)
  }


  ## Quantiles - Plotting ##
  if (nsims > 1) {
    if (qnts > 1 || qnts < 0) {
      stop("qnts must be between 0 and 1", call. = FALSE)
    }
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


  ## Mean lines - plotting ##
  if (mean.line == TRUE) {

    if (length(mean.lwd) < lcomp) {
      mean.lwd <- rep(mean.lwd, lcomp)
    }

    if (length(mean.lty) < lcomp) {
      mean.lty <- rep(mean.lty, lcomp)
    }
    draw_means(x, y, mean.smooth, mean.lwd, mean.pal, mean.lty)
  }

  ## Grid
  if (grid == TRUE) {
    grid()
  }

  ## Legends ##
  if (legend) {
    if (groups == 2 && nocomp == TRUE) {
      leg.lty <- mean.lty
    } else {
      leg.lty <- 1
    }
    legend("topright", legend = y, lty = leg.lty, lwd = 2,
           col = mean.pal, cex = leg.cex, bg = "white")
  }
}
