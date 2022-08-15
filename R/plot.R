
# Main Exported Methods ---------------------------------------------------

#' @title Plot Data from a Deterministic Compartmental Epidemic Model
#'
#' @description Plots epidemiological data from a deterministic compartment
#'              epidemic model solved with \code{\link{dcm}}.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}.
#' @param y Output compartments or flows from \code{dcm} object to plot.
#' @param popfrac If \code{TRUE}, plot prevalence of values rather than numbers
#'        (see details).
#' @param run Run number to plot, for models with multiple runs
#'        (default is run 1).
#' @param col Color for lines, either specified as a single color in a standard
#'        R color format, or alternatively as a color palette from
#'        \code{\link{RColorBrewer}} (see details).
#' @param lwd Line width for output lines.
#' @param lty Line type for output lines.
#' @param alpha Transparency level for lines, where 0 = transparent and
#'        1 = opaque (see \code{adjustcolor} function).
#' @param legend Type of legend to plot. Values are \code{"n"} for no legend,
#'        \code{"full"} for full legend, and \code{"lim"} for limited legend
#'        (see details).
#' @param leg.name Character string to use for legend, with the default
#'        determined automatically based on the \code{y} input.
#' @param leg.cex Legend scale size.
#' @param axs Plot axis type (see \code{\link{par}} for details), with default
#'        of "r".
#' @param grid If \code{TRUE}, a grid is added to the background of plot
#'        (see \code{\link{grid}} for details), with default of nx by ny.
#' @param add If \code{TRUE}, new plot window is not called and lines are added
#'        to existing plot window.
#' @param ... Additional arguments to pass to main plot window (see
#'        \code{\link{plot.default}}).
#'
#' @details
#' This function plots epidemiological outcomes from a deterministic
#' compartmental model solved with \code{\link{dcm}}. Depending on the number of
#' model runs (sensitivity analyses) and number of groups, the default plot is
#' the fractional proportion of each compartment in the model over time. The
#' specific compartments or flows to plot may be set using the \code{y}
#' parameter, and in multiple run models the specific run may also be specified.
#'
#' @section The \code{popfrac} Argument:
#' Compartment prevalence is the size of a compartment over some denominator.
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this
#' is the default. The \code{popfrac} parameter calculates
#' and plots the denominators of all specified compartments using these rules:
#' 1) for one-group models, the prevalence of any compartment is the compartment
#' size divided by the total population size; 2) for two-group models, the
#' prevalence of any compartment is the compartment size divided by the group
#' size.
#'
#' @section Color Palettes:
#' Since \code{\link{dcm}} supports multiple run sensitivity models, plotting
#' the results of such models uses a complex color scheme for distinguishing
#' runs. This is accomplished using the \code{\link{RColorBrewer}} color
#' palettes, which include a range of linked colors using named palettes. For
#' \code{plot.dcm}, one may either specify a brewer color palette listed in
#' \code{\link{brewer.pal.info}}, or, alternatively, a vector of standard R
#' colors (named, hexidecimal, or positive integers; see \code{\link{col2rgb}}).
#'
#' @section Plot Legends:
#' There are three automatic legend types available, and the legend is
#' added by default for plots. To turn off the legend, use \code{legend="n"}. To
#' plot a legend with values for every line in a sensitivity analysis, use
#' \code{legend="full"}. With models with many runs, this may be visually
#' overwhelming. In those cases, use \code{legend="lim"} to plot a legend
#' limited to the highest and lowest values of the varying parameter in the
#' model. In cases where the default legend names are not helpful, one may
#' override those names with the \code{leg.name} argument.
#'
#' @method plot dcm
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{dcm}}, \code{\link{brewer.pal.info}}
#'
#' @examples
#' # Deterministic SIR model with varying act rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = 1:10,
#'                    rec.rate = 1/3, a.rate = 0.011, ds.rate = 0.01,
#'                    di.rate = 0.03, dr.rate = 0.01)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 100, dt = 0.25)
#' mod <- dcm(param, init, control)
#'
#' # Plot disease prevalence by default
#' plot(mod)
#'
#' # Plot prevalence of susceptibles
#' plot(mod, y = "s.num", popfrac = TRUE, col = "Greys")
#'
#' # Plot number of susceptibles
#' plot(mod, y = "s.num", popfrac = FALSE, col = "Greys", grid = TRUE)
#'
#' # Plot multiple runs of multiple compartments together
#' plot(mod, y = c("s.num", "i.num"),
#'      run = 5, xlim = c(0, 50), grid = TRUE)
#' plot(mod, y = c("s.num", "i.num"),
#'      run = 10, lty = 2, legend = "n", add = TRUE)
#'
plot.dcm <- function(x, y, popfrac = FALSE, run, col, lwd, lty, alpha = 0.9,
                     legend, leg.name, leg.cex = 0.8, axs = "r", grid = FALSE,
                     add = FALSE, ...) {

  ## Set missing flags
  noy <- ifelse(missing(y), TRUE, FALSE)
  norun <- ifelse(missing(run), TRUE, FALSE)
  nocol <- ifelse(missing(col), TRUE, FALSE)
  nolwd <- ifelse(missing(lwd), TRUE, FALSE)
  nolty <- ifelse(missing(lty), TRUE, FALSE)
  noleg <- ifelse(missing(legend), TRUE, FALSE)


  ## Dot args
  da <- list(...)


  ## Model dimensions
  nsteps <- x$control$nsteps
  nruns <- x$control$nruns
  if (norun == FALSE && any(run > nruns)) {
    stop("Specify run between 1 and", nruns,
         call. = FALSE)
  }

  if (!is.null(x$control$new.mod) && noy == TRUE) {
    stop("Specify y when simulating a new model type in dcm",
         call. = FALSE)
  }

  groups <- x$param$groups
  dis.type <- x$control$type

  ## Main title default
  if (is.null(da$main)) {
    main <- ""
  } else {
    main <- da$main
  }


  ## Defaults for missing y
  if (noy == TRUE && nruns == 1) {
    y <- grep(".num", names(x$epi), value = TRUE)
  }
  if (noy == TRUE && nruns > 1) {
    y <- grep("i.num", names(x$epi), value = TRUE)
  }
  if (all(y %in% names(x$epi)) == FALSE) {
    stop("Specified y is unavailable", call. = FALSE)
  }
  lcomp <- length(y)


  ## Prevalence calculations
  x <- denom(x, y, popfrac)


  ## Compartment ymax calculations
  if (popfrac == FALSE) {
    allmax <- sapply(1:lcomp, function(i) max(x$epi[[y[i]]], na.rm = TRUE))
    ymax <- ceiling(max(allmax))
  } else {
    ymax <- 1
  }


  ## Defaults for ylim, xlim
  if (is.null(da$ylim)) {
    ylim <- c(0, ymax)
  } else {
    ylim <- da$ylim
  }
  if (is.null(da$xlim)) {
    xlim <- c(0, nsteps)
  } else {
    xlim <- da$xlim
  }


  ## Defaults for lwd
  if (nolwd == FALSE && lcomp > 1 && length(lwd) < lcomp) {
    lwd <- rep(lwd, lcomp)
  }
  if (nolwd == FALSE && lcomp == 1 && length(lwd) < nruns) {
    lwd <- rep(lwd, nruns)
  }
  if (nolwd == TRUE) {
    lwd <- rep(2.5, lcomp * nruns)
  }


  ## Defaults for lty
  if (nolty == FALSE && lcomp > 1 && length(lty) < lcomp) {
    lty <- rep(lty, lcomp)
  }
  if (nolty == FALSE && lcomp == 1 && length(lty) < nruns) {
    lty <- rep(lty, nruns)
  }
  if (nolty == TRUE) {
    lty <- rep(1, lcomp * nruns)
    if (groups == 2 && noy == TRUE) {
      lty <- rep(1:2, each = lcomp / 2)
    }
  }

  ## Defaults for xlab and ylab
  if (is.null(da$xlab)) {
    xlab <- "Time"
  } else {
    xlab <- da$xlab
  }

  if (is.null(da$ylab)) {
    if (popfrac == FALSE) {
      ylab <- "Number"
    } else {
      ylab <- "Prevalence"
    }
  } else {
    ylab <- da$ylab
  }


  ## Main plot window
  if (add == FALSE) {
    plot(1, 1, type = "n", bty = "n",
         xaxs = axs, yaxs = axs, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main)
  }


  ## Default line colors
  pal <- NULL
  # Missing col
  if (nocol == TRUE) {
    if (lcomp == 1) {
      if (nruns == 1) {
        col <- "black"
      }
      if (nruns > 1) {
        col <- "Set1"
      }
      if (nruns > 5) {
        col <- "Spectral"
      }
      if (norun == FALSE && length(run) == 1) {
        col <- "black"
      }
    }
    if (lcomp > 1) {
      col <- "Set1"
    }
  }


  # Test if using a RColorBrewer palette
  if (length(col) == 1 && col %in% row.names(brewer.pal.info)) {
    use.brewer <- TRUE
  } else {
    use.brewer <- FALSE
  }

  # Set color palette
  if (is.null(pal)) {
    if (lcomp == 1) {
      if (use.brewer == TRUE) {
        if (nruns < 6) {
          pal <- adjustcolor(brewer.pal(5, col)[1:nruns], alpha)
        } else {
          pal <- adjustcolor(brewer_ramp(nruns, col), alpha)
        }
      }
      if (use.brewer == FALSE) {
        pal <- adjustcolor(rep(col, nruns), alpha)
      }
    }
    if (lcomp > 1) {
      if (use.brewer == TRUE) {
        if (lcomp > 4) {
          pal <- adjustcolor(brewer_ramp(lcomp, col), alpha)
        } else {
          pal <- adjustcolor(brewer.pal(max(c(lcomp, 4)), col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]
          fixpal[2] <- pal[1]
          pal <- fixpal
        }
        if (groups == 2 && noy == TRUE) {
          pal <- adjustcolor(brewer.pal(3, col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]
          fixpal[2] <- pal[1]
          pal <- fixpal
          if (dis.type != "SIR") {
            pal <- pal[1:2]
          }
          pal <- rep(pal, times = lcomp / 2)
        }
      }
      if (use.brewer == FALSE) {
        pal <- adjustcolor(rep(col, lcomp), alpha)
        if (groups == 2 && noy == TRUE) {
          pal <- adjustcolor(rep(col, times = 2), alpha)
        }
      }
    }
  }


  ## Plot lines
  if (lcomp == 1) {
    if (nruns == 1) {
      lines(x$control$timesteps, x$epi[[y]][, 1],
            lwd = lwd[1], lty = lty[1], col = pal[1])
    }
    if (nruns > 1) {
      if (norun == TRUE) {
        for (i in 1:nruns) {
          lines(x$control$timesteps, x$epi[[y]][, i],
                lwd = lwd[i], lty = lty[i], col = pal[i])
        }
      } else {
        if (length(run) == 1) {
          lines(x$control$timesteps, x$epi[[y]][, run],
                lwd = lwd[1], lty = lty[1], col = pal[1])
        }
        if (length(run) > 1) {
          for (i in seq_along(run)) {
            lines(x$control$timesteps, x$epi[[y]][, run[i]],
                  lwd = lwd[i], lty = lty[i], col = pal[i])
          }
        }
      }
    }
  }
  if (lcomp > 1) {
    if (nruns == 1) {
      for (i in 1:lcomp) {
        lines(x$control$timesteps, x$epi[[y[i]]][, 1],
              lwd = lwd, lty = lty[i], col = pal[i])
      }
    }
    if (nruns > 1) {
      if (norun == TRUE) {
        for (i in 1:lcomp) {
          run <- 1
          lines(x$control$timesteps, x$epi[[y[i]]][, run],
                lwd = lwd[i], lty = lty[i], col = pal[i])
        }
      }
      if (norun == FALSE) {
        if (length(run) > 1) {
          stop("Plotting multiple runs of multiple y is not supported",
               call. = FALSE)
        }
        for (i in 1:lcomp) {
          lines(x$control$timesteps, x$epi[[y[i]]][, run],
                lwd = lwd[i], lty = lty[i], col = pal[i])
        }
      }
    }
  }

  ## Grid
  if (grid == TRUE) {
     grid()
  }

  ## Legend

  # Default legend type
  if (noleg == TRUE) {
    legend <- "n"
    if (lcomp == 1 && nruns < 3) {
      legend <- "full"
    }
    if (lcomp == 1 && nruns >= 3) {
      legend <- "lim"
    }
    if (lcomp > 1) {
      legend <- "full"
    }
    if (noy == FALSE) {
      legend <- "n"
    }
  } else {
    if (legend == "lim" && nruns < 3) {
      legend <- "full"
    }
    if (legend == "lim" && lcomp == 2) {
      legend <- "full"
    }
  }

  # Default legend names
  if (missing(leg.name)) {
    if (nruns == 1) {
      leg.names <- y
    }
    if (nruns > 1) {
      if (norun == TRUE && lcomp == 1) {
        leg.names <- names(x$epi[[y[1]]])
      }
      if (norun == FALSE && lcomp == 1) {
        if (length(run) == 1) {
          leg.names <- y
        }
        if (length(run) > 1) {
          leg.names <- names(x$epi[[y[1]]][run])
        }
      }
      if (lcomp > 1) {
        leg.names <- y
      }
    }
  } else {
    if (lcomp == 1) {
      leg.names <- paste(leg.name, 1:nruns)
    }
    if (lcomp > 1) {
      leg.names <- y
      warning("Legend names ignored for multiple y plots of multiple run
              models", call. = FALSE)
    }
  }

  # Legend
  if (norun == TRUE) {
    if (legend == "full") {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty, lwd = lwd,
             col = pal, cex = leg.cex)
    }
    if (legend == "lim") {
      legend("topright",
             legend = c(leg.names[1], "...", leg.names[nruns]),
             bg = "white",
             lty = c(lty[1], 1, lty[nruns]), lwd = lwd + 1,
             col = c(pal[1], "white", pal[nruns]), cex = leg.cex)
    }
  }
  if (norun == FALSE && legend != "n") {
    if (lcomp == 1) {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty[seq_along(run)],
             lwd = lwd[seq_along(run)],
             col = pal[seq_along(run)], cex = leg.cex)
    }
    if (lcomp > 1) {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty, lwd = lwd,
             col = pal, cex = leg.cex)
    }
  }

}


#' @title Plot Data from a Stochastic Individual Contact Epidemic Model
#'
#' @description Plots epidemiological data from a stochastic individual contact
#'              model simulated with \code{\link{icm}}.
#'
#' @param x An \code{EpiModel} model object of class \code{icm}.
#' @param y Output compartments or flows from \code{icm} object to plot.
#' @param sims A vector of simulation numbers to plot.
#' @inheritParams plot.netsim
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
plot.icm <- function(x, y, popfrac = FALSE, sim.lines = FALSE, sims, sim.col,
                     sim.lwd,  sim.alpha, mean.line = TRUE, mean.smooth = TRUE,
                     mean.col, mean.lwd = 2, mean.lty = 1, qnts = 0.5, qnts.col,
                     qnts.alpha, qnts.smooth = TRUE, legend, leg.cex = 0.8,
                     axs = "r", grid = FALSE, add = FALSE, ...) {

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
  bpal <- c(4, 2, 3)

  # Mean line
  if (missing(mean.col)) {
    mean.col <- bpal
  }
  mean.pal <- adjustcolor(mean.col, 0.9)

  # Quantile bands
  if (missing(qnts.col)) {
    qnts.col <- bpal
  }
  if (missing(qnts.alpha)) {
    qnts.alpha <- 0.5
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
  if (is.null(da$xlim)) {
    xlim <- c(0, nsteps)
  } else {
    xlim <- da$xlim
  }

  #Initialize ylim min max values
  qnt.min <- 1E10
  qnt.max <- -1E10
  mean.min <- 1E10
  mean.max <- -1E10

  ## Quantiles - ylim min max ##
  if (missing(qnts) || qnts == FALSE) {
    disp.qnts <- FALSE
  } else {
    disp.qnts <- TRUE
  }
  if (nsims == 1) {
    disp.qnts <- FALSE
  }
  if (groups == 1 && missing(qnts)) {
    disp.qnts <- TRUE
    qnts <- 0.5
  }
  if (disp.qnts == TRUE) {
    if (qnts > 1 || qnts < 0) {
      stop("qnts must be between 0 and 1", call. = FALSE)
    }
    qnt.min <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "min")
    qnt.max <- draw_qnts(x, y, qnts, qnts.pal, qnts.smooth, "epi", 0, "max")
  }

  ## Mean lines - ylim max ##
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
      mean.lty <- rep(1, lcomp)
    }
    mean.min <- draw_means(x, y, mean.smooth, mean.lwd,
                           mean.pal, mean.lty, "epi", 0, "min")
    mean.max <- draw_means(x, y, mean.smooth, mean.lwd,
                           mean.pal, mean.lty, "epi", 0, "max")
  }

  # Dynamic scaling based on sim.lines and mean lines and quantile bands
  if (!is.null(da$ylim)) {
    ylim <- da$ylim
  } else if (is.null(da$ylim) && sim.lines == FALSE &&
             (mean.line == TRUE || qnts == TRUE)) {
    ylim <- c(min(qnt.min * 0.9, mean.min * 0.9),
              max(qnt.max * 1.1, mean.max * 1.1))
  } else {
    ylim <- c(min.prev, max.prev)
  }

  if (is.null(da$main)) {
    main <- ""
  } else {
    main <- da$main
  }

  if (is.null(da$xlab)) {
    xlab <- "Time"
  } else {
    xlab <- da$xlab
  }

  if (is.null(da$ylab)) {
    if (popfrac == FALSE) {
      ylab <- "Number"
    } else {
      ylab <- "Prevalence"
    }
  } else {
    ylab <- da$ylab
  }

  ## Main plot window ##
  if (add == FALSE) {
    plot(1, 1, type = "n", bty = "n",
         xaxs = axs, yaxs = axs, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main)
  }


  ## Quantiles - Plotting ##
  if (missing(qnts) || qnts == FALSE) {
    disp.qnts <- FALSE
  } else {
    disp.qnts <- TRUE
  }
  if (nsims == 1) {
    disp.qnts <- FALSE
  }
  if (groups == 1 && missing(qnts)) {
    disp.qnts <- TRUE
    qnts <- 0.5
  }
  if (disp.qnts == TRUE) {
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
      if (nocomp == FALSE || (nocomp == TRUE && groups == 1)) {
        mean.lty <- rep(1, lcomp)
      } else {
        mean.lty <- rep(1:2, each = lcomp / 2)
      }
    }
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
}

## Helper utilities
draw_qnts <- function(x, y, qnts, qnts.pal, qnts.smooth,
                      loc = "epi", plot.qnts = 1, qnts.min_max = "max") {

  qnt.min <- 1E10
  qnt.max <- -1E10

  lcomp <- length(y)
  for (j in seq_len(lcomp)) {
    quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
    qnt.prev <- apply(x[[loc]][[y[j]]], 1,
                      function(x) {
                        quantile(x, c(quants[1], quants[2]), na.rm = TRUE)
                      })
    qnt.prev <- qnt.prev[, complete.cases(t(qnt.prev))]
    xx <- c(seq_len(ncol(qnt.prev)), rev(seq_len(ncol(qnt.prev))))
    if (qnts.smooth == FALSE) {
      yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
    } else {
      yy <- c(suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
                                      y = qnt.prev[1, ]))$y,
              rev(suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
                                          y = qnt.prev[2, ]))$y))
    }
    if (plot.qnts == 1) {
      polygon(xx, yy, col = qnts.pal[j], border = NA)
    } else {
      qnt.max[j] <-  max(yy)
      qnt.min[j] <-  min(yy)
    }
  }
  if (plot.qnts == 0 && qnts.min_max == "max") {
    return(max(qnt.max))
  } else if (plot.qnts == 0 && qnts.min_max == "min") {
    return(min(qnt.min))
  }
}


draw_means <- function(x, y, mean.smooth, mean.lwd,
                       mean.pal, mean.lty, loc = "epi",
                       plot.means = 1, mean.min_max = "max") {

  mean.min <- 1E10
  mean.max <- -1E10

  lcomp <- length(y)
  nsims <- x$control$nsims

  for (j in seq_len(lcomp)) {
    if (nsims == 1) {
      mean.prev <- x[[loc]][[y[j]]][, 1]
    } else {
      mean.prev <- rowMeans(x[[loc]][[y[j]]], na.rm = TRUE)
    }
    if (mean.smooth == TRUE) {
      mean.prev <- suppressWarnings(supsmu(x = seq_along(mean.prev),
                                           y = mean.prev))$y
    }
    if (plot.means == 1) {
      lines(mean.prev, lwd = mean.lwd[j],
            col = mean.pal[j], lty = mean.lty[j])
    } else {
      mean.max[j] <-  max(mean.prev, na.rm = TRUE)
      mean.min[j] <-  min(mean.prev, na.rm = TRUE)
    }
  }
  if (plot.means == 0 && mean.min_max == "max") {
    return(max(mean.max))
  } else if (plot.means == 0 && mean.min_max == "min") {
    return(min(mean.min))
  }
}

get_qnts <- function(data, qnts, qnts.smooth) {
  if (qnts < 0 || qnts > 1) {
    stop("qnts must be between 0 and 1", call. = FALSE)
  }
  quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
  qnt.prev <- apply(data, 1, function(x) {
    quantile(x, c(quants[1], quants[2]), na.rm = TRUE)
  })
  if (isFALSE(qnts.smooth)) {
    xx <- c(seq_len(ncol(qnt.prev)), rev(seq_len(ncol(qnt.prev))))
    yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
    xx <- xx[!is.na(yy)]
    yy <- yy[!is.na(yy)]
  } else {
    ss1 <- suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
                                   y = qnt.prev[1, ]))
    ss2 <- suppressWarnings(supsmu(x = rev(seq_len(ncol(qnt.prev))),
                                   y = qnt.prev[2, ]))

    xx <- c(ss1$x, rev(ss2$x))
    yy <- c(ss1$y, rev(ss2$y))
  }
  list(x = xx, y = yy)
}

get_means <- function(data, mean.smooth) {
  mean.prev <- rowMeans(data, na.rm = TRUE)
  if (isFALSE(mean.smooth)) {
    xx <- seq_along(mean.prev)
    yy <- mean.prev
    xx <- xx[!is.na(yy)]
    yy <- yy[!is.na(yy)]
  } else {
    ss <- suppressWarnings(supsmu(x = seq_along(mean.prev),
                                  y = mean.prev))
    xx <- ss$x
    yy <- ss$y
  }
  list(x = xx, y = yy)
}

plot_stats_table <- function(data,
                             nmstats,
                             method,
                             duration.imputed,
                             sim.lines,
                             sim.col,
                             sim.lwd,
                             mean.line,
                             mean.smooth,
                             mean.col,
                             mean.lwd,
                             mean.lty,
                             qnts,
                             qnts.col,
                             qnts.alpha,
                             qnts.smooth,
                             targ.line,
                             targ.col,
                             targ.lwd,
                             targ.lty,
                             plots.joined,
                             draw_legend,
                             grid,
                             targets,
                             dynamic,
                             da,
                             ...) {
  nstats <- length(nmstats)

  if (missing(plots.joined)) {
    plots.joined <- ifelse(nstats > 3, FALSE, TRUE)
  }

  if (nstats == 1) {
    plots.joined <- TRUE
    sim.col <- "dodgerblue3"
  }

  xlim <- NVL(da$xlim, c(1, dim(data)[1]))

  xlab <- if (isFALSE(plots.joined)) "" else NVL(da$xlab, if (isTRUE(dynamic)) "time" else "simulation number")
  ylab <- if (isFALSE(plots.joined)) "" else NVL(da$ylab, if (nstats == 1) nmstats else "Statistic")

  if (missing(sim.lwd)) {
    if (dim(data)[3] > 1) {
      sim.lwd <- max(c(1 - (dim(data)[3] * 0.05), 0.5))
    } else {
      sim.lwd <- 1
    }
  }

  ## Color Vector Validation
  # 1. Sim.col, mean.col, qnts.col, targ.col must be missing or a vector of
  #    length 1 or nstats
  # 2. If sim.col, mean.col, qnts.col, or targ.col is not missing
  #    but is a vector of length 1 and nstats is greater than 1,
  #    then replicate the color vector nstats times to achieve a vector of
  #    size nstats.
  check_len_rep <- function(object, default, name) {
    if (!missing(object)) {
      if (length(object) %in% c(1, nstats)) {
        rep(object, length.out = nstats)
      } else {
        stop(paste0(name, " must be either missing or a vector of length 1 or nstats (", nstats, ")"))
      }
    } else {
      rep(default, length.out = nstats)
    }
  }

  sim.col <- check_len_rep(sim.col,
                           if (isFALSE(plots.joined)) "dodgerblue3" else seq_len(nstats + 1L)[-1L],
                           "sim.col")

  mean.col <- check_len_rep(mean.col,
                            if (isFALSE(plots.joined)) "black" else sim.col,
                            "mean.col")

  qnts.col <- check_len_rep(qnts.col,
                            sim.col,
                            "qnts.col")
  qnts.col <- adjustcolor(qnts.col, qnts.alpha)

  targ.col <- check_len_rep(targ.col,
                            if (isFALSE(plots.joined) || nstats == 1) "black" else sim.col,
                            "targ.col")

  draw_legend <- isTRUE(plots.joined) &&
                 ((!missing(draw_legend) && isTRUE(draw_legend)) ||
                  (missing(draw_legend) && nstats == 1))

  draw_qnts <- isTRUE(dynamic) && is.numeric(qnts)

  mains <- if (isTRUE(plots.joined)) character(nstats) else nmstats

  if (method == "l") {

    qnts_list <- list()
    means <- list()
    ylims <- list()

    for (j in seq_len(nstats)) {
      dataj <- matrix(data[, j, ], nrow = dim(data)[1])

      if (isTRUE(draw_qnts)) {
        qnts_list[[j]] <- get_qnts(dataj, qnts, qnts.smooth)
      }

      if (isTRUE(mean.line)) {
        means[[j]] <- get_means(dataj, mean.smooth)
      }
    }

    for (j in seq_len(nstats)) {
      if (!is.null(da$ylim)) {
        ylims[[j]] <- da$ylim
      } else {
        limdat <- c(if (isTRUE(plots.joined) && isTRUE(sim.lines)) data,
                    if (isFALSE(plots.joined) && isTRUE(sim.lines)) data[, j, ],
                    if (isTRUE(plots.joined) && isTRUE(mean.line)) unlist(lapply(means, `[[`, "y")),
                    if (isFALSE(plots.joined) && isTRUE(mean.line)) means[[j]]$y,
                    if (isTRUE(plots.joined) && isTRUE(draw_qnts)) unlist(lapply(qnts_list, `[[`, "y")),
                    if (isFALSE(plots.joined) && isTRUE(draw_qnts)) qnts_list[[j]]$y,
                    if (isTRUE(plots.joined) && isTRUE(targ.line)) targets,
                    if (isFALSE(plots.joined) && isTRUE(targ.line)) targets[j])

        ylimsj <- suppressWarnings(c(min(limdat, na.rm = TRUE), max(limdat, na.rm = TRUE)))

        if (any(is.infinite(ylimsj))) {
          ## should both be infinite in this case, indicating no non-missing data to plot;
          ## set both limits to 0 simply to avoid errors when calling plot below
          ylimsj <- c(0, 0)
        } else {
          ## give +/- 10% buffer in a way that works for signed statistics
          ylimsj[1] <- ylimsj[1] - 0.1 * abs(ylimsj[1])
          ylimsj[2] <- ylimsj[2] + 0.1 * abs(ylimsj[2])
        }

        ylims[[j]] <- ylimsj
      }
    }

    if (isFALSE(plots.joined)) {
      if (nstats == 1) dimens <- c(1, 1)
      if (nstats == 2) dimens <- c(1, 2)
      if (nstats == 3) dimens <- c(1, 3)
      if (nstats == 4) dimens <- c(2, 2)
      if (nstats == 5) dimens <- c(2, 3)
      if (nstats == 6) dimens <- c(2, 3)
      if (nstats %in% 7:9) dimens <- c(3, 3)
      if (nstats %in% 10:12) dimens <- c(4, 3)
      if (nstats %in% 13:16) dimens <- c(4, 4)
      if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)), 2)

      # Pull graphical parameters
      ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
      par(mar = c(2.5, 2.5, 2, 1), mgp = c(2, 1, 0), mfrow = dimens)
    }

    ## do actual plotting
    for (j in seq_len(nstats)) {
      if (j == 1 || isFALSE(plots.joined)) {
        plot(NULL,
             xlim = xlim,
             ylim = ylims[[j]],
             type = "n",
             xlab = xlab,
             ylab = ylab,
             main = mains[j])
      }
      dataj <- matrix(data[, j, ], nrow = dim(data)[1])

      if (isTRUE(draw_qnts)) {
        polygon(qnts_list[[j]]$x, qnts_list[[j]]$y, col = qnts.col[j], border = NA)
      }

      if (isTRUE(sim.lines)) {
        apply(dataj,
              2,
              function(y) lines(which(!is.na(y)), y[!is.na(y)], lwd = sim.lwd, col = sim.col[j]))
      }

      if (isTRUE(mean.line)) {
        lines(means[[j]]$x,
              means[[j]]$y,
              lwd = mean.lwd,
              col = mean.col[j],
              lty = mean.lty)
      }

      if (isTRUE(targ.line)) {
        abline(h = targets[j],
               lty = targ.lty,
               lwd = targ.lwd,
               col = targ.col[j])
      }

      if (isTRUE(grid) && isFALSE(plots.joined)) {
        grid()
      }
    }

    if (isTRUE(grid) && isTRUE(plots.joined)) {
      grid()
    }

    if (isTRUE(draw_legend)) {
      legend("topleft", legend = nmstats, lwd = 2,
             col = sim.col[1:nstats], cex = 0.75, bg = "white")
    }

    if (isFALSE(plots.joined)) {
      # Reset graphical parameters
      on.exit(par(ops))
    }
  }

  if (method == "b") {

    data <- matrix(aperm(data, c(1, 3, 2)), nrow = dim(data)[1] * dim(data)[3])
    colnames(data) <- nmstats

    boxplot(data, ...)

    for (j in seq_len(nstats)) {
      points(x = j, y = targets[j],
             pch = 16, cex = 1.5, col = "blue")

      ## Grid
      if (isTRUE(grid)) {
          grid()
      }
    }
  }
}

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
#'        \code{"b"} for boxplots.
#' @param sims A vector of simulation numbers to plot.
#' @param stats Statistics to plot.  For \code{type = "formation"},
#'        \code{stats} are among those specified in the call to
#'        \code{\link{netdx}}; for \code{type = "duration","dissolution"},
#'        \code{stats} are among those of the dissolution model (without
#'        \code{offset()}).  The default is to plot all statistics.
#' @param plots.joined If \code{TRUE}, combine all statistics in one
#'        plot, versus one plot per statistic if \code{FALSE}.
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
      warning("duration/dissolution data contains undefined values due to",
              " having zero edges of some dissolution dyad type(s) on some time",
              " step(s); these undefined values will be set to 0 when",
              " processing the data; this behavior, which introduces a bias",
              " towards 0, may be changed in the future")
    }

    if (any(grepl("nodefactor", x$dissolution) == TRUE)) {
      warning("Support for dissolution models containing a nodefactor term is
              deprecated, and will be removed in a future release.",
              call. = FALSE)
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
#' @param axs Plot axis type (see \code{\link{par}} for details), with default
#'        of \code{"r"}.
#' @param grid If \code{TRUE}, a grid is added to the background of plot
#'        (see \code{\link{grid}} for details), with default of nx by ny.
#' @param add If \code{TRUE}, new plot window is not called and lines are added
#'        to existing plot window.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#' @param at If \code{type="network"}, time step for network graph.
#' @param col.status If \code{TRUE} and \code{type="network"}, automatic disease
#'        status colors (blue = susceptible, red = infected, green = recovered).
#' @param shp.g2 If \code{type="network"} and \code{x} is for a two-group model,
#'        shapes for the Group 2 vertices, with acceptable inputs of "triangle"
#'        and "square". Group 1 vertices will remain circles.
#' @param vertex.cex Relative size of plotted vertices if \code{type="network"},
#'        with implicit default of 1.
#' @param stats If \code{type="formation","duration","dissolution"}, statistics
#'        to plot.  For \code{type = "formation"}, \code{stats} are among those
#'        specified in \code{nwstats.formula} of \code{\link{control.net}}; for
#'        \code{type = "duration","dissolution"}, \code{stats} are among those
#'        of the dissolution model (without \code{offset()}).  The default is
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
#' @param method Plot method for
#'        \code{type="formation","duration","dissolution"}, with options of
#'        \code{"l"} for line plots and \code{"b"} for boxplots.
#' @param duration.imputed If \code{type="duration"}, a logical indicating
#'        whether or not to impute starting times for relationships extant at
#'        the start of the simulation. Defaults to \code{TRUE} when
#'        \code{type="duration"}.
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
                        qnts.smooth = TRUE, legend, leg.cex = 0.8, axs = "r",
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
      if (length(tsn) != length(ts)) {
        ts <- ts[which(ts > 0)]
      }
      names(ts) <- tsn
      
    } else {
      ## duration/dissolution plot
      if (isTRUE(x$control$save.diss.stats) &&
          isTRUE(x$control$save.network) &&
          isFALSE(x$control$tergmLite) &&
          isFALSE(is.null(x$diss.stats)) &&
          isTRUE(x$nwparam[[network]]$coef.diss$diss.model.type == "edgesonly")) {

        if (any(unlist(lapply(x$diss.stats, `[[`, "anyNA")))) {
          warning("duration/dissolution data contains undefined values due to",
                  " having zero edges of some dissolution dyad type(s) on some time",
                  " step(s); these undefined values will be set to 0 when",
                  " processing the data; this behavior, which introduces a bias",
                  " towards 0, may be changed in the future")
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
          ts <- 1/x$nwparam[[network]]$coef.diss$duration
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


#' @title Plot Compartment Diagram for Epidemic Models
#'
#' @description Plots a compartment flow diagram for deterministic compartmental
#'              models, stochastic individual contact models, and stochastic
#'              network models.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}, \code{icm}, or
#'        \code{netsim}.
#' @param at Time step for model statistics.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments passed to plot (not currently used).
#'
#' @details
#' The \code{comp_plot} function provides a visual summary of an epidemic model
#' at a specific time step. The information contained in \code{comp_plot} is the
#' same as in the \code{summary} functions for a model, but presented
#' graphically as a compartment flow diagram.
#'
#' For \code{dcm} class plots, specify the model run number if the model
#' contains multiple runs, as in a sensitivity analysis. For \code{icm} and
#' \code{netsim} class plots, the \code{run} argument is not used; the plots
#' show the means and standard deviations across simulations at the specified
#' time step.
#'
#' These plots are currently limited to one-group models for each of the three
#' model classes. That functionality may be expanded in future software
#' releases.
#'
#' @export
#' @keywords plot
#'
#' @examples
#' ## Example 1: DCM SIR model with varying act.rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = 5:7,
#'                    rec.rate = 1/3, a.rate = 1/90, ds.rate = 1/100,
#'                    di.rate = 1/35, dr.rate = 1/100)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 25, verbose = FALSE)
#' mod1 <- dcm(param, init, control)
#' comp_plot(mod1, at = 25, run = 3)
#'
#' ## Example 2: ICM SIR model with 3 simulations
#' param <- param.icm(inf.prob = 0.2, act.rate = 3, rec.rate = 1/50,
#'                    a.rate = 1/100, ds.rate = 1/100,
#'                    di.rate = 1/90, dr.rate = 1/100)
#' init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.icm(type = "SIR", nsteps = 25,
#'                        nsims = 3, verbose = FALSE)
#' mod2 <- icm(param, init, control)
#' comp_plot(mod2, at = 25, digits = 1)
#'
comp_plot <- function(x, at, digits, ...) {
  UseMethod("comp_plot")
}


#' @param run Model run number, for \code{dcm} class models with multiple runs
#'        (sensitivity analyses).
#' @method comp_plot dcm
#' @rdname comp_plot
#' @export
comp_plot.dcm <- function(x, at = 1, digits = 3, run = 1, ...) {


  ## Variables
  nsteps <- x$control$nsteps
  dis.type <- x$control$type
  groups <- x$param$groups
  vital <- x$param$vital

  ## Errors
  if (groups != 1) {
    stop("Only 1-group dcm models currently supported",
         call. = FALSE)
  }

  ## Time
  if (at > nsteps || at < 1) {
    stop("Specify a time step between 1 and ", nsteps)
  }
  intime <- at
  at <- which(x$control$timesteps == intime)

  ## Dataframe subsets
  df <- as.data.frame(x, run = run)
  df <- round(df[at, ], digits)

  ## Change graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(0, 0, 2, 0))
  options(scipen = 10)

  ## Main Plot
  plot(0:100, 0:100, type = "n", axes = FALSE)
  title(main = paste(dis.type, "Model Diagram"))
  mtext(paste0("time=", intime, "  |  run=", run),
        side = 3, cex = 0.8, line = -1)

  ## 1. SI Model
  if (dis.type == "SI") {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)
    harrow(22, 40, "si.flow", df$si.flow, dir = "right")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  ## 2. SIR Model
  if (dis.type == "SIR") {
    mbox(5, 40, "Susceptible", df$s.num)
    mbox(40, 40, "Infected", df$i.num)
    mbox(75, 40, "Recovered", df$r.num)
    harrow(5, 40, "si.flow", df$si.flow, dir = "right")
    harrow(40, 40, "ir.flow", df$ir.flow, dir = "right")
    if (vital == TRUE) {
      varrow(5, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(40, 40, "di.flow", df$di.flow, dir = "out")
      varrow(75, 40, "dr.flow", df$dr.flow, dir = "out")
      varrow(5, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  ## 3. SIS Model
  if (dis.type == "SIS") {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)
    harrow(22, 40, "si.flow", df$si.flow, dir = "right")
    harrow(22, 40, "is.flow", df$is.flow, dir = "left")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  # Reset graphical parameters
  on.exit(par(ops))
}


#' @method comp_plot icm
#' @rdname comp_plot
#' @export
comp_plot.icm <- function(x, at = 1, digits = 3, ...) {

  # Variables
  nsteps <- x$control$nsteps
  dis.type <- x$control$type
  vital <- x$param$vital

  # Standardize groups
  if (inherits(x, "icm")) {
    groups <- x$param$groups
  }
  if (inherits(x, "netsim")) {
    groups <- x$param$groups
  }
  if (groups != 1) {
    stop("Only 1-group models currently supported",
         call. = FALSE)
  }

  # Time
  if (at > nsteps || at < 1) {
    stop("Specify a timestep between 1 and ", nsteps,
         call. = FALSE)
  }

  ## Dataframe subsets for plots
  df.mn <- as.data.frame(x, out = "mean")
  df.mn <- round(df.mn[at == df.mn$time, ], digits)
  df.sd <- as.data.frame(x, out = "sd")
  df.sd <- round(df.sd[at == df.sd$time, ], digits)

  ## Change graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(0, 0, 2, 0))
  options(scipen = 10)

  ## Main Plot
  plot(0:100, 0:100, type = "n", axes = FALSE)
  title(main = paste(dis.type, "Model Diagram"))
  mtext(paste0("Simulation means(sd) | time=", at),
        side = 3, cex = 0.8, line = -1)

  ## 1. SI Model
  if (dis.type == "SI" && groups == 1) {
    mbox(22, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(57, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    harrow(22, 40, "si.flow", df.mn$si.flow, dir = "right")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  ## 2. SIR Model
  if (dis.type == "SIR" && groups == 1) {
    mbox(5, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(40, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    mbox(75, 40, "Recovered", paste0(df.mn$r.num, "(", df.sd$r.num, ")"))
    harrow(5, 40, "si.flow", df.mn$si.flow, dir = "right")
    harrow(40, 40, "ir.flow", df.mn$ir.flow, dir = "right")
    if (vital == TRUE) {
      varrow(5, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(40, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(75, 40, "dr.flow", df.mn$dr.flow, dir = "out")
      varrow(5, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  ## 3. SIS Model
  if (dis.type == "SIS" && groups == 1) {
    mbox(22, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(57, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    harrow(22, 40, "si.flow", df.mn$si.flow, dir = "right")
    harrow(22, 40, "is.flow", df.mn$is.flow, dir = "left")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  # Reset graphical parameters
  on.exit(par(ops))

}

#' @method comp_plot netsim
#' @rdname comp_plot
#' @export
comp_plot.netsim <- function(x, at = 1, digits = 3, ...) {

  comp_plot.icm(x = x, at = at, digits = digits, ...)

}



# ggplot ------------------------------------------------------------------

#' @title ggplot2 Geom for Quantile Bands
#'
#' @description Plots quantile bands given a data.frame with stochastic model
#'              results from \code{\link{icm}} or \code{\link{netsim}}.
#'
#' @param mapping Standard aesthetic mapping \code{aes()} input for ggplot2.
#' @param lower Lower quantile for the time series.
#' @param upper Upper quantile for the time series.
#' @param alpha Transparency of the ribbon fill.
#' @param ... Additional arguments passed to \code{stat_summary}.
#'
#' @details
#' This is a wrapper around \code{ggplot::stat_summary} with a ribbon geom as
#' aesthetic output.
#'
#' @export
#' @keywords plot
#'
#' @examples
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SI", nsteps = 250, nsims = 5)
#' mod1 <- icm(param, init, control)
#' df <- as.data.frame(mod1)
#' df.mean <- as.data.frame(mod1, out = "mean")
#'
#' library(ggplot2)
#' ggplot() +
#'    geom_line(data = df, mapping = aes(time, i.num, group = sim),
#'    alpha = 0.25, lwd = 0.25, color = "firebrick") +
#'    geom_bands(data = df, mapping = aes(time, i.num),
#'               lower = 0.1, upper = 0.9, fill = "firebrick") +
#'    geom_line(data = df.mean, mapping = aes(time, i.num)) +
#'    theme_minimal()
#'
geom_bands <- function(mapping, lower = 0.25, upper = 0.75, alpha = 0.25, ...) {
  stat_summary(mapping,
               geom = "ribbon",
               fun.min = function(x) quantile(x, lower),
               fun.max = function(x) quantile(x, upper),
               alpha = alpha, ...)
}


# Helper Functions --------------------------------------------------------


# Calculate denominators
denom <- function(x, y, popfrac) {

  cont.val <- ifelse(class(x) == "dcm", "nruns", "nsims")
  if (popfrac == TRUE) {
    for (i in seq_along(y)) {
      dname <- paste(strsplit(y[i], "[.]")[[1]][-1], collapse = ".")
      x$epi[[y[i]]] <- x$epi[[y[i]]] / x$epi[[dname]]
    }
  }
  if (popfrac == FALSE && x$control[[cont.val]] == 1) {
    for (j in seq_along(y)) {
      x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
    }
  }

  return(x)
}

## comp_plot helper utilities ##
#  Text box
mbox <- function(x, y, title, val) {
  polygon(c(x, x + 20, x + 20, x), c(y, y, y + 20, y + 20))
  text(x + 10, y + 10, paste(title, "\n n=", val, sep = ""), cex = 0.9)
}
#  Horizontal arrow
harrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "right") {
    arrows(xbox + 20, ybox + 12, xbox + 35, lwd = 2, length = 0.15)
    text(xbox + 27.5, ybox + 17, paste(title, val, sep = "="), cex = 0.8)
  }
  if (dir == "left") {
    arrows(xbox + 20 + 15, ybox + 5, xbox + 20, lwd = 2, length = 0.15)
    text(xbox + 27.5, ybox + 2, paste(title, val, sep = "="), cex = 0.8)
  }
}
#  Vertical arrow
varrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "out") {
    arrows(xbox + 10, ybox, xbox + 10, ybox - 25, lwd = 2, length = 0.15)
    text(xbox + 10, ybox - 12.5, paste(title, val, sep = "="),
         cex = 0.8, pos = 4)
  }
  if (dir == "in") {
    arrows(xbox + 10, ybox + 45, xbox + 10, ybox + 20, lwd = 2, length = 0.15)
    text(xbox + 10, ybox + 32.5, paste(title, val, sep = "="),
         cex = 0.8, pos = 4)
  }
}
