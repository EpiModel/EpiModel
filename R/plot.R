

# Main Exported Methods ---------------------------------------------------

#' @title Plot Data from a Deterministic Compartmental Epidemic Model
#'
#' @description Plots epidemiological data from a deterministic compartment
#'              epidemic model solved with \code{dcm}.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}.
#' @param y Output compartments or flows from \code{dcm} object to plot.
#' @param popfrac If \code{TRUE}, plot prevalence of values rather than numbers
#'        (see details).
#' @param run Run number to plot, for models with multiple runs (default is run 1).
#' @param col Color for lines, either specified as a single color in a standard
#'        R color format, or alternatively as a color palette from
#'        \code{\link{RColorBrewer}} (see details).
#' @param lwd Line width for output lines.
#' @param lty Line type for output lines.
#' @param alpha Transparency level for lines, where 0 = transparent and 1 = opaque
#'        (see \code{\link{transco}}).
#' @param leg Type of legend to plot. Values are "n" for no legend, "full" for
#'        full legend, and "lim" for limited legend (see details).
#' @param leg.name Character string to use for legend, with the default
#'        determined automatically based on the \code{y} input.
#' @param leg.cex Legend scale size.
#' @param axs Plot axis type (see \code{\link{par}} for details), with default
#'        of "r".
#' @param add If \code{TRUE}, new plot window is not called and lines are added to
#'        existing plot window.
#' @param ... Additional arguments to pass to main plot window (see
#'        \code{\link{plot.default}}).
#'
#' @details
#' This function plots epidemiological outcomes from a deterministic
#' compartmental model solved with \code{\link{dcm}}. Depending on the number of
#' model runs (sensitivity analyses) and number of groups, the default plot is
#' the fractional proportion of each compartment in the model over time. The
#' specific compartments or flows to plot may be set using the \code{y} parameter,
#' and in multiple run models the specific run may also be specified.
#'
#' @section The popfrac Argument:
#' Compartment prevalences are the size of a compartment over some denominator.
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this
#' is the default for any plots of flows. The \code{popfrac} parameter calculates
#' and plots the denominators of all specified compartments using these rules: 1)
#' for one-group models, the prevalence of any compartment is the compartment size
#' divided by the total population size; 2) for two-group models, the prevalence
#' of any compartment is the compartment size divided by the group size.
#'
#' @section Color Palettes:
#' Since \code{\link{dcm}} supports multiple run sensitivity models, plotting
#' the results of such models uses a complex color scheme for distinguishing runs.
#' This is accomplished using the \code{\link{RColorBrewer}} color palettes, in
#' which includes a range of linked colors using named palettes. For
#' \code{plot.dcm}, one may either specify a brewer color palette listed in
#' \code{\link{brewer.pal.info}}, or alternatively a vector of standard R colors
#' (named, hexidecimal, or positive integers; see \code{\link{col2rgb}}).
#'
#' @section Plot Legends:
#' There are three automatic legend types available, and the legend is
#' added by default for plots. To turn off the legend, use \code{leg="n"}. To
#' plot a legend with values for every line in a sensitivity analysis, use
#' \code{leg="full"}. With models with many runs, this may be visually
#' overwhelming. In those cases, use \code{leg="lim"} to plot a legend limited
#' to the highest and lowest of the varying parameter in the model. In cases
#' where the default legend names are not helpful, one may override those names
#' with the \code{leg.name} argument.
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
#'                    rec.rate = 1/3, b.rate = 0.011, ds.rate = 0.01,
#'                    di.rate = 0.03, dr.rate = 0.01)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 100, dt = 0.25)
#' mod <- dcm(param, init, control)
#'
#' # Plot disease prevalence by default
#' plot(mod)
#'
#' # Plot prevalence of susceptibles
#' plot(mod, y = "s.num", col = "Greys")
#'
#' # Plot number of susceptibles
#' plot(mod, y = "s.num", popfrac = FALSE, col = "Greys")
#'
#' # Plot multiple runs of multiple compartments together
#' plot(mod, y = c("s.num", "i.num"),
#'      run = 5, xlim = c(0, 50))
#' plot(mod, y = c("s.num", "i.num"),
#'      run = 10, lty = 2, leg = "n", add = TRUE)
#'
plot.dcm <- function(x, y, popfrac, run, col, lwd, lty, alpha = 0.9, leg,
                     leg.name, leg.cex = 0.8, axs = "r", add = FALSE, ...) {


  ## Set missing flags
  noy <- ifelse(missing(y), TRUE, FALSE)
  norun <- ifelse(missing(run), TRUE, FALSE)
  nocol <- ifelse(missing(col), TRUE, FALSE)
  nolwd <- ifelse(missing(lwd), TRUE, FALSE)
  nolty <- ifelse(missing(lty), TRUE, FALSE)
  noleg <- ifelse(missing(leg), TRUE, FALSE)


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
  if (missing(popfrac)) {
    popfrac <- TRUE
  }
  if (any(grepl(".flow", y)) | !is.null(x$control$new.mod)) {
    popfrac <- FALSE
  }
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
          pal <- transco(brewer.pal(5, col)[1:nruns], alpha)
        } else {
          pal <- transco(brewer_ramp(nruns, col), alpha)
        }
      }
      if (use.brewer == FALSE) {
        pal <- transco(rep(col, nruns), alpha)
      }
    }
    if (lcomp > 1) {
      if (use.brewer == TRUE) {
        if (lcomp > 4) {
          pal <- transco(brewer_ramp(lcomp, col), alpha)
        } else {
          pal <- transco(brewer.pal(max(c(lcomp, 4)), col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
        }
        if (groups == 2 && noy == TRUE) {
          pal <- transco(brewer.pal(3, col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
          if (dis.type != "SIR") {
            pal <- pal[1:2]
          }
          pal <- rep(pal, times = lcomp / 2)
        }
      }
      if (use.brewer == FALSE) {
        pal <- transco(rep(col, lcomp), alpha)
        if (groups == 2 && noy == TRUE) {
          pal <- transco(rep(col, times = 2), alpha)
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
          for (i in 1:length(run)) {
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


  ## Legend

  # Default legend type
  if (noleg == TRUE) {
    leg <- "n"
    if (lcomp == 1 & nruns < 3) {
      leg <- "full"
    }
    if (lcomp == 1 & nruns >= 3) {
      leg <- "lim"
    }
    if (lcomp > 1) {
      leg <- "full"
    }
    if (noy == FALSE) {
      leg <- "n"
    }
  } else {
    if (leg == "lim" & nruns < 3) {
      leg <- "full"
    }
    if (leg == "lim" & lcomp == 2) {
      leg <- "full"
    }
  }

  # Default legend names
  if (missing(leg.name)) {
    if (nruns == 1) {
      leg.names <- y
    }
    if (nruns > 1) {
      if (norun == TRUE & lcomp == 1) {
        leg.names <- names(x$epi[[y[1]]])
      }
      if (norun == FALSE & lcomp == 1) {
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
      warning("Legend names ignored for multiple y plots of multiple run models",
              call. = FALSE)
    }
  }

  # Legend
  if (norun == TRUE) {
    if (leg == "full") {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty, lwd = lwd,
             col = pal, cex = leg.cex)
    }
    if (leg == "lim") {
      legend("topright",
             legend = c(leg.names[1], "...", leg.names[nruns]),
             bg = "white",
             lty = c(lty[1], 1, lty[nruns]), lwd = lwd + 1,
             col = c(pal[1], "white", pal[nruns]), cex = leg.cex)
    }
  }
  if (norun == FALSE & leg != "n") {
    if (lcomp == 1) {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty[1:length(run)],
             lwd = lwd[1:length(run)],
             col = pal[1:length(run)], cex = leg.cex)
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
#'              model simulated with \code{icm}.
#'
#' @inheritParams plot.netsim
#'
#' @method plot icm
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{icm}}
#'
#' @examples
#' \dontrun{
#' ## Example 1: Plotting multiple compartment values from SIR model
#' param <- param.icm(inf.prob = 0.5, act.rate = 0.5, rec.rate = 0.02)
#' init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.icm(type = "SIR", nsteps = 100,
#'                        nsims = 3, verbose = FALSE)
#' mod <- icm(param, init, control)
#' plot(mod)
#'
#' ## Example 2: Plot only infected with specific output from SI model
#' param <- param.icm(inf.prob = 0.25, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 10)
#' control <- control.icm(type = "SI", nsteps = 100,
#'                        nsims = 3, verbose = FALSE)
#' mod2 <- icm(param, init, control)
#'
#' # Plot prevalence
#' plot(mod2, y = "i.num", mean.line = FALSE)
#'
#' # Plot incidence
#' plot(mod2, y = "si.flow", mean.smooth = TRUE)
#' plot(mod2, y = "si.flow", qnts.smooth = FALSE, qnts = 1)
#' }
#'
plot.icm <- function(x, y, popfrac, sim.lines = FALSE, sims, sim.col, sim.lwd,
                     sim.alpha, mean.line = TRUE, mean.smooth = TRUE,
                     mean.col, mean.lwd = 2, mean.lty = 1, qnts = 0.5, qnts.col,
                     qnts.alpha, qnts.smooth = TRUE, leg, leg.cex = 0.8,
                     axs = "r", add = FALSE, ...) {

  ## Model dimensions and class ##
  nsteps <- x$control$nsteps
  nsims <- x$control$nsims
  if (missing(sims)) {
    sims <- 1:nsims
  }
  if (max(sims) > nsims) {
    stop("Set sim to between 1 and ", nsims, call. = FALSE)
  }
  dis.type <- x$control$type
  if (is.null(x$param$groups) | !is.numeric(x$param$groups)) {
    modes <- 1
    x$param$groups <- 1
  } else {
    modes <- x$param$groups
  }

  # dotargs
  da <- list(...)


  ## Compartments ##
  nocomp <- ifelse(missing(y), TRUE, FALSE)
  if (nocomp == TRUE) {
    if (modes == 1) {
      y <- grep(".num$", names(x$epi), value = TRUE)
    }
    if (modes == 2) {
      if (class(x) == "icm") {
        y <- c(grep(".num$", names(x$epi), value = TRUE),
               grep(".num.g2$", names(x$epi), value = TRUE))
      }
    }
    if (missing(leg)) {
      leg <- TRUE
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
  bpal <- brewer.pal(3, "Set1")
  bpal <- c(bpal[2], bpal[1], bpal[3])

  # Mean line
  if (missing(mean.col)) {
    mean.col <- bpal
  }
  mean.pal <- transco(mean.col, 0.9)

  # Quantile bands
  if (missing(qnts.col)) {
    qnts.col <- bpal
  }
  if (missing(qnts.alpha)) {
    qnts.alpha <- 0.4
  }
  qnts.pal <- transco(qnts.col, qnts.alpha)

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

  if (missing(sim.alpha) & nsims == 1) {
    sim.alpha <- 0.9
  }
  if (missing(sim.alpha) & nsims > 1) {
    sim.alpha <- max(c(0.05, 1 - log10(nsims) / 3))
  }
  sim.pal <- transco(sim.col, sim.alpha)

  # Special case for 2-mode/group models
  if (modes == 2 & nocomp == TRUE) {
    pal <- brewer.pal(3, "Set1")
    if (dis.type == "SIR") {
      mean.pal <- rep(mean.pal, 2)
      qnts.pal <- rep(qnts.pal, 2)
      sim.pal <- rep(sim.pal, 2)
    } else {
      mean.pal <- rep(mean.pal[1:2], 2)
      qnts.pal <- rep(qnts.pal[1:2], 2)
      sim.pal <- rep(sim.pal[1:2], 2)
    }
  }


  ## Prevalence calculations ##
  nopopfrac <- ifelse(missing(popfrac), TRUE, FALSE)
  if (nopopfrac == TRUE) {
    popfrac <- TRUE
  }
  if (nopopfrac == TRUE) {
    if (any(grepl(".flow", y)) |
        (modes == 1 & all(grepl(".num$", y)) == FALSE) |
        (modes == 2 & all(c(grepl(".num$", y), grepl(".m2$", y)) == FALSE)) |
        any(y %in% c("num", "num.m2", "num.g2"))) {
      popfrac <- FALSE
    }
  }
  x <- denom(x, y, popfrac)

  # Compartment max
  if (popfrac == FALSE) {
    if (lcomp == 1) {
      max.prev <- max(x$epi[[y]], na.rm = TRUE)
    } else {
      max.prev <- max(sapply(y, function(comps) max(x$epi[[comps]], na.rm = TRUE)))
    }
  } else {
    max.prev <- 1
  }


  ## Missing args ##
  if (is.null(da$xlim)) {
    xlim <- c(0, nsteps)
  } else {
    xlim <- da$xlim
  }
  if (is.null(da$ylim)) {
    ylim <- c(0, max.prev)
  } else {
    ylim <- da$ylim
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


  ## Quantiles ##
  if (missing(qnts) || qnts == FALSE) {
    disp.qnts <- FALSE
  } else {
    disp.qnts <- TRUE
  }
  if (nsims == 1) {
    disp.qnts <- FALSE
  }
  if (modes == 1 & missing(qnts)) {
    disp.qnts <- TRUE
    qnts <- 0.5
  }
  if (disp.qnts == TRUE) {
    if (qnts > 1 | qnts < 0) {
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
      if (nocomp == FALSE || (nocomp == TRUE && modes == 1)) {
        mean.lty <- rep(1, lcomp)
      } else {
        mean.lty <- rep(1:2, each = lcomp / 2)
      }
    }
    draw_means(x, y, mean.smooth, mean.lwd, mean.pal, mean.lty)
  }


  ## Legends ##
  if (!missing(leg) && leg == TRUE) {
    if (modes == 2 & nocomp == TRUE) {
      leg.lty <- mean.lty
    } else {
      leg.lty <- 1
    }
    legend("topright", legend = y, lty = leg.lty, lwd = 3,
           col = mean.pal, cex = leg.cex, bg = "white")
  }


}


## Helper utilities
draw_qnts <- function(x, y, qnts, qnts.pal, qnts.smooth, loc = "epi") {

  lcomp <- length(y)
  for (j in seq_len(lcomp)) {
    quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
    qnt.prev <- apply(x[[loc]][[y[j]]], 1,
                      function(x) {
                        quantile(x, c(quants[1], quants[2]), na.rm = TRUE)
                      })
    qnt.prev <- qnt.prev[, complete.cases(t(qnt.prev))]
    xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
    if (qnts.smooth == FALSE) {
      yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
    } else {
      yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
              rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
    }
    polygon(xx, yy, col = qnts.pal[j], border = NA)
  }

}


draw_means <- function(x, y, mean.smooth, mean.lwd,
                       mean.pal, mean.lty, loc = "epi") {

  lcomp <- length(y)
  nsims <- x$control$nsims

  for (j in seq_len(lcomp)) {
    if (nsims == 1) {
      mean.prev <- x[[loc]][[y[j]]][, 1]
    } else {
      mean.prev <- rowMeans(x[[loc]][[y[j]]])
    }
    if (mean.smooth == TRUE) {
      mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
    }
    lines(mean.prev, lwd = mean.lwd[j],
          col = mean.pal[j], lty = mean.lty[j])
  }

}


#' @title Plot Dynamic Network Model Diagnostics
#'
#' @description Plots dynamic network model diagnostics calculated in
#'              \code{netdx}.
#'
#' @param x An \code{EpiModel} object of class \code{netdx}.
#' @param type Plot type, with options of \code{"formation"} for network
#'        model formation statistics, \code{"duration"} for dissolution model
#'        statistics for average edge duration, or \code{"dissolution"} for
#'        dissolution model statistics for proportion of ties dissolved per time
#'        step.
#' @param method Plot method, with options of \code{"l"} for line plots and
#'        \code{"b"} for boxplots.
#' @param stats Network statistics to plot, among those specified in the call
#'        to \code{\link{netdx}}, with the default to plot all statistics
#'        contained in the object.
#' @inheritParams plot.netsim
#'
#' @details
#' The plot function for \code{netdx} objects will generate plots of two types of
#' model diagnostic statistics that run as part of the diagnostic tools within
#' that function. The \code{formation} plot shows the summary statistics
#' requested in \code{nwstats.formula}, where the default includes those
#' statistics in the network model formation formula specified in the original
#' call to \code{\link{netest}}.
#'
#' The \code{duration} plot shows the average age of existing edges at each time
#' step, up until the maximum time step requested. This is calculated with the
#' \code{\link{edgelist_meanage}} function. The age is used as an estimator of
#' the average duration of edges in the equilibrium state.
#'
#' The \code{dissolution} plot shows the proportion of the extant ties that are
#' dissolved at each time step, up until the maximum time step requested.
#' Typically the proportion of ties that are dissolved is the reciprocal of the
#' mean relational duration. This plot thus contains similar information to that
#' in the duration plot, but should reach its expected value more quickly, since
#' it is not subject to censoring.
#'
#' The \code{plots.joined} argument will control whether the statistics
#' in the \code{formation} plot are joined in one plot or plotted separately.
#' The default is based on the number of network statistics requested. The
#' layout of the separate plots within the larger plot window is also based on
#' the number of statistics.
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
#' nw <- network.initialize(100, directed = FALSE)
#' nw <- set.vertex.attribute(nw, "sex", rbinom(100, 1, 0.5))
#' formation <- ~edges + nodematch("sex")
#' target.stats <- c(50, 40)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
#'
#' # Estimate the model
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Static diagnostics
#' dx1 <- netdx(est, nsims = 1e4, dynamic = FALSE,
#'              nwstats.formula = ~edges + meandeg + concurrent +
#'                                 nodefactor("sex", base = 0) +
#'                                 nodematch("sex"))
#' dx1
#'
#' # Only formation diagnostics are available to plot
#' plot(dx1, stats = "edges")
#' plot(dx1, stats = c("edges", "concurrent"))
#' plot(dx1, stats = "edges", method = "b", col = "seagreen3")
#' plot(dx1, stats = c("nodefactor.sex.0", "nodefactor.sex.1"),
#'      method = "b", col = transco(2:3, 0.5))
#'
#' # Dynamic diagnostics
#' dx2 <- netdx(est, nsims = 10, nsteps = 500,
#'              nwstats.formula = ~edges + meandeg + concurrent +
#'                                 nodefactor("sex", base = 0) +
#'                                 nodematch("sex"))
#' dx2
#'
#' # Formation statistics plots, joined and separate
#' plot(dx2)
#' plot(dx2, type = "formation", plots.joined = TRUE)
#' plot(dx2, type = "formation", sim = 1, plots.joined = TRUE,
#'      qnts = FALSE, sim.lines = TRUE, mean.line = FALSE)
#' plot(dx2, type = "formation", plots.joined = FALSE,
#'      stats = c("edges", "concurrent"))
#' plot(dx2, type = "formation", stats = "nodefactor.sex.0",
#'      sim = 1, sim.lwd = 5, sim.col = "darkmagenta")
#'
#' plot(dx2, method = "b", col = "bisque")
#' plot(dx2, method = "b", stats = "meandeg", col = "dodgerblue")
#'
#' # Duration statistics plot
#' plot(dx2, type = "duration", mean.col = "black")
#' plot(dx2, type = "duration", sim = 10, mean.line = FALSE, sim.line = TRUE,
#'      sim.col = "steelblue", sim.lwd = 3, targ.lty = 1, targ.lwd = 0.5)
#'
#' # Dissolution statistics plot
#' plot(dx2, type = "dissolution", mean.col = "black")
#' plot(dx2, type = "dissolution", method = "b", col = "pink1")
#' }
#'
plot.netdx <- function(x, type = "formation", method = "l", sims, stats,
                       sim.lines, sim.col, sim.lwd, mean.line = TRUE,
                       mean.smooth = TRUE, mean.col, mean.lwd = 2, mean.lty = 1,
                       qnts = 0.5, qnts.col, qnts.alpha, qnts.smooth = TRUE,
                       targ.line = TRUE, targ.col, targ.lwd = 2, targ.lty = 2,
                       plots.joined, leg, ...) {

  # Checks and Variables ----------------------------------------------------

  ## Check Object
  if (class(x) != "netdx") {
    stop("x must be an object of class netdx", call. = FALSE)
  }

  if (x$dynamic == FALSE && type %in% c("duration", "dissolution")) {
    stop("Plots of type duration and dissolution only available if netdx ",
         "run with dynamic = TRUE", call. = FALSE)
  }

  ## Check sims
  nsims <- x$nsims
  if (missing(sims)) {
    sims <- 1:nsims
  }
  if (max(sims) > nsims) {
    stop("Maximum sim number is", nsims, call. = FALSE)
  }
  nsteps <- x$nsteps
  dynamic <- x$dynamic

  # Get dotargs
  da <- list(...)


  # Formation Plot ----------------------------------------------------------
  if (type == "formation") {

    ## Stats
    nwstats <- x$stats
    nwstats.table <- x$stats.table.formation

    ## Find available stats
    sts <- which(!is.na(nwstats.table[, "Sim Mean"]))
    nmstats <- rownames(nwstats.table)[sts]

    ## Pull and check stat argument
    if (missing(stats)) {
      stats <- nmstats
    }
    if (any(stats %in% nmstats == FALSE)) {
      stop("One or more requested stats not contained in netdx object",
           call. = FALSE)
    }
    outsts <- which(nmstats %in% stats)

    ## Subset data
    nstats <- length(outsts)
    data <- do.call("cbind", args = nwstats)
    data <- data[, colnames(data) %in% nmstats[outsts], drop = FALSE]

    ## Pull target stats
    targs <- which(!is.na(nwstats.table$Target))


    ## Plotting
    if (missing(plots.joined)) {
      plots.joined <- ifelse(nstats > 3, FALSE, TRUE)
    }

    if (nstats == 1) {
      plots.joined <- TRUE
    }

    if (dynamic == TRUE) {
      xlim <- c(1, nsteps)
      if (length(da) > 0 && !is.null(da$xlim)) {
        xlim <- da$xlim
      }
    } else {
      xlim <- c(1, nsims)
      if (length(da) > 0 && !is.null(da$xlim)) {
        xlim <- da$xlim
      }
    }

    if (missing(sim.lwd)) {
      if (nsims == 1 | dynamic == FALSE) {
        sim.lwd <- 1
      } else {
        sim.lwd <- max(c(1 - (nsims * 0.05), 0.5))
      }
    }

    # Default colors
    if (missing(sim.col)) {
      if (nstats > 8) {
        sim.col <- brewer_ramp(nstats, "Set1")
      } else {
        sim.col <- brewer.pal(9, "Set1")[1:(nstats + 1)]
        if (nstats >= 6) {
          sim.col <- sim.col[-which(sim.col == "#FFFF33")]
        }
      }
    }


    ## Joined Plots
    if (method == "l") {

      if (missing(sim.lines)) {
        if (dynamic == TRUE) {
          sim.lines <- FALSE
        } else {
          sim.lines <- TRUE
        }
      }

      if (plots.joined == TRUE) {

        ## Default legend
        if (nstats == 1) {
          if (missing(leg)) {
            leg <- FALSE
          }
        } else {
          leg <- TRUE
        }

        ## Default ylim
        ylim <- c(min(data) * 0.8, max(data) * 1.2)
        if (length(da) > 0 && !is.null(da$ylim)) {
          ylim <- da$ylim
        }

        ## Default ylab
        if (length(da) > 0 && !is.null(da$ylab)) {
          ylab <- da$ylab
        } else {
          if (nstats == 1) {
            ylab <- nmstats[outsts]
          } else {
            ylab <- "Statistic"
          }
        }

        ## Default xlab
        if (dynamic == TRUE) {
          if (length(da) > 0 && !is.null(da$xlab)) {
            xlab <- da$xlab
          } else {
            xlab <- "time"
          }
        } else {
          if (length(da) > 0 && !is.null(da$xlab)) {
            xlab <- da$xlab
          } else {
            xlab <- "simulation number"
          }
        }


        ## Default target line color
        if (missing(targ.col)) {
          if (nstats == 1) {
            targ.col <- "black"
          } else {
            targ.col <- sim.col
          }
        }


        ## Main plot window
        plot(1, 1, xlim = xlim, ylim = ylim,
             type = "n", xlab = xlab, ylab = ylab)
        for (j in outsts) {
          dataj <- data[, colnames(data) %in% nmstats[j], drop = FALSE]

          if (is.numeric(qnts)) {
            if (qnts < 0 | qnts > 1) {
              stop("qnts must be between 0 and 1", call. = FALSE)
            }
            if (missing(qnts.col)) {
              qnts.col <- sim.col
            }
            if (missing(qnts.alpha)) {
              qnts.alpha <- 0.35
            }
            qnts.col <- transco(qnts.col, qnts.alpha)
            quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
            qnt.prev <- apply(dataj, 1, function(x) {
                 quantile(x, c(quants[1], quants[2]))
              })
            xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
            if (qnts.smooth == FALSE) {
              yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
            } else {
              yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                      rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
            }
            polygon(xx, yy, col = qnts.col[which(j == outsts)], border = NA)
          }

          if (sim.lines == TRUE) {
            if (dynamic == TRUE) {
              for (i in sims) {
                lines(dataj[,i], lwd = sim.lwd, col = sim.col[which(j == outsts)])
              }
            } else {
              lines(dataj, lwd = sim.lwd, col = sim.col[which(j == outsts)])
            }
          }
          if (mean.line == TRUE) {
            if (missing(mean.col)) {
              mean.col <- sim.col
            }
            mean.prev <- rowMeans(dataj)
            if (mean.smooth == TRUE) {
              mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
            }
            lines(mean.prev, lwd = mean.lwd,
                  col = mean.col[which(j == outsts)], lty = mean.lty)
          }

          if (targ.line == TRUE) {
            if (j %in% targs) {
              abline(h = nwstats.table$Target[j],
                     lty = targ.lty, lwd = targ.lwd,
                     col = targ.col[which(j == outsts)])
            }
          }

        }
        if (leg == TRUE) {
          legend("topleft", legend = nmstats[outsts], lwd = 3,
                 col = sim.col[1:nstats], cex = 0.75, bg = "white")
        }

      }

      ## Split plots
      if (plots.joined == FALSE) {

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

        if (missing(targ.col)) {
          targ.col <- rep("black", nstats)
        }

        for (j in outsts) {
          dataj <- data[, colnames(data) %in% nmstats[j], drop = FALSE]
          plot(x = 1, y = 1,
               xlim = xlim,
               ylim = c(min(dataj) * 0.8, max(dataj) * 1.2),
               type = "n", main = nmstats[j],
               xlab = "", ylab = "")

          if (is.numeric(qnts)) {
            if (qnts < 0 | qnts > 1) {
              stop("qnts must be between 0 and 1", call. = FALSE)
            }
            if (missing(qnts.col)) {
              qnts.col <- sim.col
            }
            if (missing(qnts.alpha)) {
              qnts.alpha <- 0.35
            }
            qnts.col <- transco(qnts.col, qnts.alpha)
            quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
            qnt.prev <- apply(dataj, 1, function(x) {
                  quantile(x, c(quants[1], quants[2]))
              })
            xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
            if (qnts.smooth == FALSE) {
              yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
            } else {
              yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                      rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
            }
            polygon(xx, yy, col = qnts.col[which(j == outsts)], border = NA)
          }

          if (sim.lines == TRUE) {
            if (dynamic == TRUE) {
              for (i in sims) {
                lines(dataj[, i], lwd = sim.lwd, col = sim.col[which(j == outsts)])
              }
            } else {
              lines(dataj, lwd = sim.lwd, col = sim.col[which(j == outsts)])
            }
          }

          if (mean.line == TRUE) {
            if (missing(mean.col)) {
              mean.col <- sim.col
            }
            mean.prev <- rowMeans(dataj)
            if (mean.smooth == TRUE) {
              mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
            }
            lines(mean.prev, lwd = mean.lwd,
                  col = mean.col[which(j == outsts)], lty = mean.lty)
          }

          if (targ.line == TRUE) {
            if (j %in% targs) {
              abline(h = nwstats.table$Target[j],
                     lty = targ.lty, lwd = targ.lwd,
                     col = targ.col[which(j == outsts)])
            }
          }
        }

        # Reset graphical parameters
        on.exit(par(ops))
      }
    }
    if (method == "b") {

      data <- list()
      for (j in seq_along(stats)) {
        data[[j]] <- unlist(lapply(nwstats, function(x) x[,stats[j]]))
      }
      data <- do.call("cbind", data)
      colnames(data) <- stats

      boxplot(data, ...)

      for (j in outsts) {
        if (j %in% targs) {
          points(x = outsts[j], y = nwstats.table$Target[j],
                 pch = 16, cex = 1.5, col = "blue")
        }
      }

    }

  }

  # Duration plot -----------------------------------------------------------

  if (type == "duration") {

    if (x$coef.diss$model.type == "hetero") {
      stop("Duration plots for heterogeneous dissolution models not currently available",
           call. = FALSE)
    }

    pages <- x$pages

    xlim <- c(1, nsteps)
    if (length(da) > 0 & !is.null(da$xlim)) {
     xlim <- da$xlim
    }

    ylim <- c(0, max(sapply(pages, max, na.rm = TRUE)) * 1.2)
    if (length(da) > 0 & !is.null(da$ylim)) {
      ylim <- da$ylim
    }

    if (missing(sim.lwd)) {
      sim.lwd <- max(c(1 - (nsims * 0.05), 0.5))
    }

    if (missing(sim.col)) {
      sim.col <- rep("grey20", nsims)
    }
    if (missing(targ.col)) {
      targ.col <- "black"
    }

    # Default ylab
    if (length(da) > 0 && !is.null(da$ylab)) {
      ylab <- da$ylab
    } else {
      ylab <- "Edge Age"
    }

    # Default xlab
    if (length(da) > 0 && !is.null(da$xlab)) {
      xlab <- da$xlab
    } else {
      xlab <- "time"
    }

    if (method == "l") {

      if (missing(sim.lines)) {
        sim.lines <- FALSE
      }

      plot(x = 1, y = 1, type = "n",
           xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab)

      if (is.numeric(qnts) & nsims > 1) {
        if (qnts < 0 | qnts > 1) {
          stop("qnts must be between 0 and 1", call. = FALSE)
        }
        if (missing(qnts.col)) {
          qnts.col <- sim.col
        }
        if (missing(qnts.alpha)) {
          qnts.alpha <- 0.35
        }
        qnts.col <- transco(qnts.col, qnts.alpha)
        quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
        dataj <- as.data.frame(pages)
        dataj <- dataj[complete.cases(dataj), , drop = FALSE]
        qnt.prev <- apply(dataj, 1, function(x)
                          quantile(x, c(quants[1], quants[2]),
                                   na.rm = TRUE))
        xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
        if (qnts.smooth == FALSE) {
          yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
        } else {
          yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                  rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
        }
        polygon(xx, yy, col = qnts.col, border = NA)
      }

      ## Sim lines
      if (sim.lines == TRUE) {
        for (i in sims) {
          lines(pages[[i]], lwd = sim.lwd, col = sim.col)
        }
      }

      if (mean.line == TRUE) {
        if (missing(mean.col)) {
          mean.col <- sim.col
        }
        dataj <- as.data.frame(pages)
        dataj <- dataj[complete.cases(dataj), , drop = FALSE]
        mean.prev <- rowMeans(dataj)
        if (mean.smooth == TRUE) {
          mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
        }
        lines(mean.prev, lwd = mean.lwd,
              col = mean.col, lty = mean.lty)
      }

      if (targ.line == TRUE) {
        abline(h = as.numeric(x$coef.diss[2]),
               lty = targ.lty, lwd = targ.lwd,
               col = targ.col)
      }
    }

    if (method == "b") {
      data <- do.call("c", x$pages)
      boxplot(data, ...)
      points(x = 1, y = as.numeric(x$coef.diss[2]),
             pch = 16, cex = 1.5, col = "blue")
    }

  }

  # Dissolution plot -----------------------------------------------------------
  if (type == "dissolution") {

    if (x$coef.diss$model.type == "hetero") {
      stop("Dissolution plots for heterogeneous dissolution models not currently available",
           call. = FALSE)
    }

    prop.diss <- x$prop.diss

    xlim <- c(1, nsteps)
    if (length(da) > 0 & !is.null(da$xlim)) {
      xlim <- da$xlim
    }

    ylim <- c(0, max(sapply(prop.diss, max, na.rm = TRUE)) * 1.1)
    if (length(da) > 0 & !is.null(da$ylim)) {
      ylim <- da$ylim
    }

    if (missing(sim.lwd)) {
      sim.lwd <- max(c(1 - (nsims * 0.05), 0.5))
    }

    if (missing(sim.col)) {
      sim.col <- rep("grey20", nsims)
    }
    if (missing(targ.col)) {
      targ.col <- "black"
    }

    # Default ylab
    if (length(da) > 0 && !is.null(da$ylab)) {
      ylab <- da$ylab
    } else {
      ylab <- "Pct Edges Diss"
    }

    # Default xlab
    if (length(da) > 0 && !is.null(da$xlab)) {
      xlab <- da$xlab
    } else {
      xlab <- "time"
    }

    if (method == "l") {

      if (missing(sim.lines)) {
        sim.lines <- FALSE
      }

      plot(x = 1, y = 1, type = "n",
           xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab)

      if (is.numeric(qnts)) {
        if (qnts < 0 | qnts > 1) {
          stop("qnts must be between 0 and 1", call. = FALSE)
        }
        if (missing(qnts.col)) {
          qnts.col <- sim.col
        }
        if (missing(qnts.alpha)) {
          qnts.alpha <- 0.35
        }
        qnts.col <- transco(qnts.col, qnts.alpha)
        quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
        dataj <- as.data.frame(prop.diss)
        qnt.prev <- apply(dataj, 1, function(x)
          quantile(x, c(quants[1], quants[2]),
                   na.rm = TRUE))
        xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
        if (qnts.smooth == FALSE) {
          yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
        } else {
          yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                  rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
        }
        polygon(xx, yy, col = qnts.col, border = NA)
      }

      # Sim lines
      if (sim.lines == TRUE) {
        for (i in sims) {
          lines(prop.diss[[i]], lwd = sim.lwd, col = sim.col)
        }
      }

      if (mean.line == TRUE) {
        if (missing(mean.col)) {
          mean.col <- sim.col
        }
        dataj <- as.data.frame(prop.diss)
        mean.prev <- rowMeans(dataj)
        if (mean.smooth == TRUE) {
          mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
        }
        lines(mean.prev, lwd = mean.lwd,
              col = mean.col, lty = mean.lty)
      }

      if (targ.line == TRUE) {
        abline(h = as.numeric(1 / (x$coef.diss[2]$duration)),
               lty = targ.lty, lwd = targ.lwd, col = targ.col)
      }
    }

    if (method == "b") {
      data <- do.call("c", x$prop.diss)
      boxplot(data, ...)
      points(x = 1, y = as.numeric(1 / (x$coef.diss[2]$duration)),
             pch = 16, cex = 1.5, col = "blue")
    }

  }

}


#' @title Plot Data from a Stochastic Network Epidemic Model
#'
#' @description Plots epidemiological and network data from a stochastic network
#'              model simulated with \code{netsim}.
#'
#' @param x An \code{EpiModel} model object of class \code{netsim}.
#' @param type Type of plot: \code{"epi"} for epidemic model results,
#'        \code{"network"} for a static network plot (\code{plot.network}),
#'        or \code{"formation"} for network formation statistics.
#' @param y Output compartments or flows from \code{icm} object to plot.
#' @param popfrac If \code{TRUE}, plot prevalence of values rather than numbers
#'        (see details).
#' @param sim.lines If \code{TRUE}, plot individual simulation lines. Default is
#'        to plot lines for one-group models but not for two-group models.
#' @param sims If \code{type="epi"} or \code{"formation"}, a vector of
#'        simulation numbers to plot. If \code{type="network"}, a single
#'        simulation number for network plot, or else \code{"min"} to plot the
#'        simulation number with the lowest disease prevalence, \code{"max"} for
#'        the simulation with the highest disease prevalence, or \code{"mean"}
#'        for the simulation with the prevalance closest to the mean across
#'        simulations at the specified time step.
#' @param sim.col Vector of any standard R color format for simulation lines.
#' @param sim.lwd Line width for simulation lines.
#' @param sim.alpha Transparency level for simulation lines, where 0 = transparent
#'        and 1 = opaque (see \code{\link{transco}}).
#' @param mean.line If \code{TRUE}, plot mean of simulations across time.
#' @param mean.smooth If \code{TRUE}, use a lowess smoother on the mean line.
#' @param mean.col Vector of any standard R color format for mean lines.
#' @param mean.lwd Line width for mean lines.
#' @param mean.lty Line type for mean lines.
#' @param qnts If numeric, plot polygon of simulation quantiles based on the
#'        range implied by the argument (see details). If \code{FALSE}, suppress
#'        polygon from plot.
#' @param qnts.col Vector of any standard R color format for polygons.
#' @param qnts.alpha Transparency level for quantile polygons, where 0 =
#'        transparent and 1 = opaque (see \code{\link{transco}}).
#' @param qnts.smooth If \code{TRUE}, use a lowess smoother on quantile polygons.
#' @param leg If \code{TRUE}, plot default legend.
#' @param leg.cex Legend scale size.
#' @param axs Plot axis type (see \code{\link{par}} for details), with default
#'        to \code{"r"}.
#' @param add If \code{TRUE}, new plot window is not called and lines are added to
#'        existing plot window.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#' @param at If \code{type="network"}, time step for network graph.
#' @param col.status If \code{TRUE} and \code{type="network"}, automatic disease
#'        status colors (blue = susceptible, red = infected, , green = recovered).
#' @param shp.bip If \code{type="network"} and a bipartite simulation, shapes
#'        for the mode 2 vertices, with acceptable inputs of "triangle" and
#'        "square". Mode 1 vertices will be circles.
#' @param stats If \code{type="formation"}, network statistics to plot, among
#'        those specified in \code{nwstats.formula} of \code{\link{control.net}},
#'        with the default to plot all statistics.
#' @param targ.line If \code{TRUE}, plot target or expected value line for
#'        the statistic of interest.
#' @param targ.col Vector of standard R colors for target statistic lines, with
#'        default colors based on \code{RColorBrewer} color palettes.
#' @param targ.lwd Line width for the line showing the target statistic values.
#' @param targ.lty Line type for the line showing the target statistic values.
#' @param plots.joined If \code{TRUE} and \code{type="formation"}, combine all
#'        target statistics in one plot, versus one plot per target statistic if
#'        \code{FALSE}.
#' @param ... additional arguments to pass.
#'
#' @details
#' This plot function can produce three types of plots with a stochastic network
#' model simulated through \code{\link{netsim}}:
#' \enumerate{
#'  \item \strong{\code{type="epi"}}: epidemic model results (e.g., disease
#'        prevalence and incidence) may be plotted.
#'  \item \strong{\code{type="network"}}: a static network plot will be generated.
#'        A static network plot of a dynamic network is a cross-sectional
#'        extraction of that dynamic network at a specific time point. This
#'        plotting function wraps the \code{\link{plot.network}} function in the
#'        \code{network} package. Consult the help page for \code{plot.network}
#'        for all the plotting parameters. In addition, four plotting parameters
#'        specific to \code{netsim} plots are available: \code{sim}, \code{at},
#'        \code{col.status}, and \code{shp.bip}.
#'  \item \strong{\code{type="formation"}}: summary network statistics related to
#'        the network model formation are plotted. These plots are similar to the
#'        formation plots for \code{netdx} objects. When running a \code{netsim}
#'        simulation, one must specify there that \code{save.nwstats=TRUE}; the
#'        plot here will then show the network statistics requested explicitly in
#'        \code{nwstats.formula}, or will use the formation formula set in
#'        \code{netest} otherwise.
#' }
#'
#' @details
#' When \code{type="epi"}, this plotting function will extract the epidemiological
#' output from a model object of class \code{netsim} and plot the time series
#' data of disease prevalence and other results. The summary statistics that the
#' function calculates and plots are individual simulation lines, means of the
#' individual simulation lines, and quantiles of those individual simulation lines.
#' The mean line, toggled on with \code{mean.line=TRUE} is calculated as the row mean
#' across simulations at each time step.
#'
#' Compartment prevalences are the size of a compartment over some denominator.
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this
#' is the default for any plots of flows. The \code{popfrac} parameter calculates
#' and plots the denominators of all specified compartments using these rules: 1)
#' for one-group models, the prevalence of any compartment is the compartment size
#' divided by the total population size; 2) for two-group models, the prevalence
#' of any compartment is the compartment size divided by the group population size.
#'
#' The quantiles show the range of outcome values within a certain specified
#' quantile range. By default, the interquartile range is shown: that is the
#' middle 50\% of the data. This is specified by \code{qnts=0.5}. To show the
#' middle 95\% of the data, specify \code{qnts=0.95}. To toggle off the polygons
#' where they are plotted by default, specify \code{qnts=FALSE}.
#'
#' @method plot netsim
#' @export
#'
#' @keywords plot
#' @seealso \code{\link{plot.network}}
#'
#' @examples
#' \dontrun{
#' ## Independent SI Model
#' # Initialize network and set network model parameters
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Estimate the network model
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the epidemic model
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        verbose = FALSE, save.nwstats = TRUE,
#'                        nwstats.formula = ~edges + meandeg + concurrent)
#' mod <- netsim(est, param, init, control)
#'
#' # Plot epidemic trajectory (default type)
#' plot(mod, type = "epi")
#' plot(mod, type = "epi", popfrac = FALSE)
#' plot(mod, type = "epi", y = "si.flow", qnts = 1)
#'
#' # Plot static networks
#' par(mar = c(0,0,0,0))
#' plot(mod, type = "network")
#'
#' # Automatic coloring of infected nodes as red
#' par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
#' plot(mod, type = "network", main = "Min Prev | Time 50",
#'      col.status = TRUE, at = 50, sims = "min")
#' plot(mod, type = "network", main = "Max Prev | Time 50",
#'      col.status = TRUE, at = 50, sims = "max")
#'
#' # Automatic shape by mode number (circle = mode 1)
#' par(mar = c(0,0,0,0))
#' plot(mod, type = "network", at = 50, col.status = TRUE, shp.bip = "square")
#' plot(mod, type = "network", at = 50, col.status = TRUE, shp.bip = "triangle")
#'
#' # Plot formation statistics
#' par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
#' plot(mod, type = "formation")
#' plot(mod, type = "formation", plots.joined = FALSE)
#' plot(mod, type = "formation", sim = 2:4)
#' plot(mod, type = "formation", plots.joined = FALSE,
#'      stats = c("edges", "concurrent"))
#' plot(mod, type = "formation", stats = "meandeg",
#'      sim.lwd = 2, sim.col = "seagreen")
#' }
#'
plot.netsim <- function(x, type = "epi", y, popfrac, sim.lines = FALSE, sims, sim.col,
                        sim.lwd, sim.alpha, mean.line = TRUE, mean.smooth = TRUE,
                        mean.col, mean.lwd = 2, mean.lty = 1, qnts = 0.5, qnts.col,
                        qnts.alpha, qnts.smooth = TRUE, leg, leg.cex = 0.8, axs = "r",
                        add = FALSE, network = 1, at = 1, col.status = FALSE,
                        shp.bip = NULL, stats, targ.line = TRUE, targ.col,
                        targ.lwd = 2, targ.lty = 2, plots.joined, ...) {

  # Network plot ------------------------------------------------------------
  if (type == "network") {

    if (is.null(x$network)) {
      stop("networkDynamic object not saved in netsim simulation",
           call. = FALSE)
    }

    nsteps <- x$control$nsteps
    if (at > x$control$nsteps) {
      stop("Specify a time step between 1 and", nsteps)
    }

    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1
    }
    if (length(sims) > 1 || (!is.numeric(sims) && !(sims %in% c("mean", "max", "min")))) {
      stop("sims argument must be single simulation number",
           "or \"mean\", \"max\", or \"min\" ", call. = FALSE)
    }

    if (sims == "mean") {
      sims <- which.min(abs(as.numeric(x$epi$i.num[at, ]) -
                           mean(as.numeric(x$epi$i.num[at, ]))))
    }
    if (sims == "max") {
      sims <- as.numeric(which.max(x$epi$i.num[at, ]))
    }
    if (sims == "min") {
      sims <- as.numeric(which.min(x$epi$i.num[at, ]))
    }

    obj <- get_network(x, sims, network, collapse = TRUE, at = at)
    tea.status <- x$control$tea.status

    if (!is.null(shp.bip)) {
      if (all(shp.bip != c("square", "triangle"))) {
        stop("shp.bip accepts inputs of either \"square\" or \"triangle\" ",
             call. = FALSE)
      }

      if (is.numeric(obj$gal$bipartite)) {
        mids <- idmode(obj)
        if (shp.bip == "square") {
          vertex.sides <- ifelse(mids == 1, 50, 4)
          vertex.rot <- 45
          vertex.cex <- ifelse(mids == 1, 1, 1.3)
        }
        if (shp.bip == "triangle") {
          vertex.sides <- ifelse(mids == 1, 50, 3)
          vertex.rot <- 90
          vertex.cex <- ifelse(mids == 1, 1, 1.4)
        }

      } else {
        warning("shp.bip applies to bipartite networks only, so ignoring argument")
        vertex.sides <- 50
        vertex.rot <- 0
        vertex.cex <- 1
      }
    } else {
      vertex.sides <- 50
      vertex.rot <- 0
      vertex.cex <- 1
    }
    if (col.status == TRUE) {
      if (is.null(tea.status) || tea.status == FALSE) {
        stop("Plotting status colors requires tea.status=TRUE in netsim control settings",
             call. = FALSE)
      }
      pal <- transco(c("firebrick", "steelblue", "seagreen"), 0.75)
      if (tea.status == TRUE) {
        cols <- ifelse(get.vertex.attribute.active(obj, "testatus", at = at) == "i",
                       pal[1], pal[2])
        cols <- ifelse(get.vertex.attribute.active(obj, "testatus", at = at) == "r",
                       pal[3], cols)
      }
      plot.network(obj, vertex.col = cols, vertex.border = "grey60",
                   edge.col = "grey40", vertex.sides = vertex.sides,
                   vertex.rot = vertex.rot, vertex.cex = vertex.cex,
                   displaylabels = FALSE, ...)
    } else {
      plot.network(obj, vertex.sides = vertex.sides, vertex.rot = vertex.rot,
                   vertex.cex = vertex.cex, displaylabels = FALSE, ...)
    }

  }


  # Epidemic plot -----------------------------------------------------------
  if (type == "epi") {

    ## Model dimensions and class ##
    nsteps <- x$control$nsteps
    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1:nsims
    }
    if (max(sims) > nsims) {
      stop("Set sim to between 1 and ", nsims, call. = FALSE)
    }
    dis.type <- x$control$type
    if (is.null(x$param$modes) | !is.numeric(x$param$modes)) {
      modes <- 1
      x$param$modes <- 1
    } else {
      modes <- x$param$modes
    }


    # dotargs
    da <- list(...)


    ## Compartments ##
    nocomp <- ifelse(missing(y), TRUE, FALSE)
    if (nocomp == TRUE) {
      if (modes == 1) {
        y <- grep(".num$", names(x$epi), value = TRUE)
      }
      if (modes == 2) {
        if (class(x) == "icm") {
          y <- c(grep(".num$", names(x$epi), value = TRUE),
                 grep(".num.g2$", names(x$epi), value = TRUE))
        }
        if (class(x) == "netsim") {
          y <- c(grep(".num$", names(x$epi), value = TRUE),
                 grep(".num.m2$", names(x$epi), value = TRUE))
        }
      }
      if (missing(leg)) {
        leg <- TRUE
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
    bpal <- brewer.pal(3, "Set1")
    bpal <- c(bpal[2], bpal[1], bpal[3])

    # Mean line
    if (missing(mean.col)) {
      mean.col <- bpal
    }
    mean.pal <- transco(mean.col, 0.9)

    # Quantile bands
    if (missing(qnts.col)) {
      qnts.col <- bpal
    }
    if (missing(qnts.alpha)) {
      qnts.alpha <- 0.4
    }
    qnts.pal <- transco(qnts.col, qnts.alpha)

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

    if (missing(sim.alpha) & nsims == 1) {
      sim.alpha <- 0.9
    }
    if (missing(sim.alpha) & nsims > 1) {
      sim.alpha <- max(c(0.05, 1 - log10(nsims) / 3))
    }
    sim.pal <- transco(sim.col, sim.alpha)

    # Special case for 2-mode/group models
    if (modes == 2 & nocomp == TRUE) {
      pal <- brewer.pal(3, "Set1")
      if (dis.type == "SIR") {
        mean.pal <- rep(mean.pal, 2)
        qnts.pal <- rep(qnts.pal, 2)
        sim.pal <- rep(sim.pal, 2)
      } else {
        mean.pal <- rep(mean.pal[1:2], 2)
        qnts.pal <- rep(qnts.pal[1:2], 2)
        sim.pal <- rep(sim.pal[1:2], 2)
      }
    }


    ## Prevalence calculations ##
    nopopfrac <- ifelse(missing(popfrac), TRUE, FALSE)
    if (nopopfrac == TRUE) {
      popfrac <- TRUE
    }
    if (nopopfrac == TRUE) {
      if (any(grepl(".flow", y)) |
          (modes == 1 & all(grepl(".num$", y)) == FALSE) |
          (modes == 2 & all(c(grepl(".num$", y), grepl(".m2$", y)) == FALSE)) |
          any(y %in% c("num", "num.m2", "num.g2"))) {
        popfrac <- FALSE
      }
    }
    x <- denom(x, y, popfrac)

    # Compartment max
    if (popfrac == FALSE) {
      if (lcomp == 1) {
        max.prev <- max(x$epi[[y]], na.rm = TRUE)
      } else {
        max.prev <- max(sapply(y, function(comps) max(x$epi[[comps]], na.rm = TRUE)))
      }
    } else {
      max.prev <- 1
    }


    ## Missing args ##
    if (is.null(da$xlim)) {
      xlim <- c(0, nsteps)
    } else {
      xlim <- da$xlim
    }
    if (is.null(da$ylim)) {
      ylim <- c(0, max.prev)
    } else {
      ylim <- da$ylim
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


    ## Quantiles ##
    if (missing(qnts) || qnts == FALSE) {
      disp.qnts <- FALSE
    } else {
      disp.qnts <- TRUE
    }
    if (nsims == 1) {
      disp.qnts <- FALSE
    }
    if (modes == 1 & missing(qnts)) {
      disp.qnts <- TRUE
      qnts <- 0.5
    }
    if (disp.qnts == TRUE) {
      if (qnts > 1 | qnts < 0) {
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
        if (nocomp == FALSE || (nocomp == TRUE && modes == 1)) {
          mean.lty <- rep(1, lcomp)
        } else {
          mean.lty <- rep(1:2, each = lcomp / 2)
        }
      }
      draw_means(x, y, mean.smooth, mean.lwd, mean.pal, mean.lty)
    }


    ## Legends ##
    if (!missing(leg) && leg == TRUE) {
      if (modes == 2 & nocomp == TRUE) {
        leg.lty <- mean.lty
      } else {
        leg.lty <- 1
      }
      legend("topright", legend = y, lty = leg.lty, lwd = 3,
             col = mean.pal, cex = leg.cex, bg = "white")
    }
  }


  # Formation plot ----------------------------------------------------------
  if (type == "formation") {

    ## Stats
    nsims <- x$control$nsims
    if (missing(sims)) {
      sims <- 1:nsims
    }
    if (max(sims) > nsims) {
      stop("Maximum sims for this object is ", nsims, call. = FALSE)
    }

    nwstats <- get_nwstats(x, sims, network)
    nsims <- length(sims)
    nsteps <- x$control$nsteps

    ## Find available stats
    nmstats <- colnames(nwstats[[1]])

    ## Pull and check stat argument
    if (missing(stats)) {
      stats <- nmstats
    }
    if (any(stats %in% nmstats == FALSE)) {
      stop("One or more requested stats not contained in netsim object",
           call. = FALSE)
    }
    outsts <- which(nmstats %in% stats)
    nstats <- length(outsts)

    names(nwstats) <- rep("", length(nwstats))
    data <- as.matrix(do.call("cbind", args = nwstats))
    data <- data[, colnames(data) %in% nmstats[outsts], drop = FALSE]

    ## target stats
    nwparam <- get_nwparam(x, network)
    formation.terms <- nwparam$target.stats.names
    target.stats <- nwparam$target.stats

    st <- data.frame(sorder = 1:length(nmstats), names = nmstats)
    ts <- data.frame(names = formation.terms, targets = target.stats)

    m <- merge(st, ts, all = TRUE)
    m <- m[do.call("order", m[, "sorder", drop = FALSE]), , drop = FALSE]
    targs <- which(!is.na(m$targets))

    ## Get dotargs
    da <- list(...)

    ## Plotting
    if (missing(plots.joined)) {
      plots.joined <- ifelse(nstats > 3, FALSE, TRUE)
    }
    if (nstats == 1) {
      plots.joined <- TRUE
    }
    xlim <- c(1, nsteps)
    if (!is.null(da$xlim)) {
      xlim <- da$xlim
    }
    if (missing(sim.lwd)) {
      if (nsims == 1) {
        sim.lwd <- 1
      } else {
        sim.lwd <- max(c(1 - (nsims * 0.05), 0.5))
      }
    }

    # Default colors
    if (missing(sim.col)) {
      if (nstats > 8) {
        sim.col <- brewer_ramp(nstats, "Set1")
      } else {
        sim.col <- brewer.pal(9, "Set1")[1:(nstats + 1)]
        if (nstats >= 6) {
          sim.col <- sim.col[-which(sim.col == "#FFFF33")]
        }
      }
    }


    ## Joined Plots
    if (plots.joined == TRUE) {

      ## Default legend
      if (nstats == 1) {
        if (missing(leg)) {
          leg <- FALSE
        }
      } else {
        leg <- TRUE
      }

      ## Default ylim
      if (ncol(nwstats[[1]]) == 1) {
        mins <- sapply(nwstats, function(x) apply(x, 2, min))
        maxs <- sapply(nwstats, function(x) apply(x, 2, max))
      } else {
        mins <- sapply(nwstats, function(x) apply(x, 2, min))[outsts, ]
        maxs <- sapply(nwstats, function(x) apply(x, 2, max))[outsts, ]
      }
      ymin <- min(mins)
      ymax <- max(maxs)
      if (!is.null(da$ylim)) {
        ylim <- da$ylim
      } else {
        ylim <- c(ymin * 0.8, ymax * 1.2)
      }

      ## Default ylab
      if (!is.null(da$ylab)) {
        ylab <- da$ylab
      } else {
        if (nstats == 1) {
          ylab <- nmstats[outsts]
        } else {
          ylab <- "Statistic"
        }
      }

      ## Default xlab
      if (!is.null(da$xlab)) {
        xlab <- da$xlab
      } else {
        xlab <- "time"
      }

      ## Default target line color
      if (missing(targ.col)) {
        if (nstats == 1) {
          targ.col <- "black"
        } else {
          targ.col <- sim.col
        }
      }


      ## Main plot window
      plot(1, 1, xlim = xlim, ylim = ylim,
           type = "n", xlab = xlab, ylab = ylab)
      for (j in outsts) {

        dataj <- data[, colnames(data) %in% nmstats[j], drop = FALSE]

        if (is.numeric(qnts)) {
          if (qnts < 0 | qnts > 1) {
            stop("qnts must be between 0 and 1", call. = FALSE)
          }
          if (missing(qnts.col)) {
            qnts.col <- sim.col
          }
          if (missing(qnts.alpha)) {
            qnts.alpha <- 0.35
          }
          qnts.col <- transco(qnts.col, qnts.alpha)
          quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
          qnt.prev <- apply(dataj, 1, function(x) {
            quantile(x, c(quants[1], quants[2]))
          })
          xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
          if (qnts.smooth == FALSE) {
            yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
          } else {
            yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                    rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
          }
          polygon(xx, yy, col = qnts.col[which(j == outsts)], border = NA)
        }

        if (sim.lines == TRUE) {
          for (i in sims) {
            lines(dataj[,i], lwd = sim.lwd,
                  col = sim.col[which(j == outsts)])
          }
        }
        if (mean.line == TRUE) {
          if (missing(mean.col)) {
            mean.col <- sim.col
          }
          mean.prev <- rowMeans(dataj)
          if (mean.smooth == TRUE) {
            mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
          }
          lines(mean.prev, lwd = mean.lwd,
                col = mean.col[which(j == outsts)], lty = mean.lty)
        }

        if (targ.line == TRUE) {
          if (j %in% targs) {
            abline(h = m$targets[j], lty = targ.lty, lwd = targ.lwd,
                   col = targ.col[which(j == outsts)])
          }
        }

      }
      if (leg == TRUE) {
        legend("topleft", legend = nmstats[outsts], lwd = 3,
               col = sim.col[1:nstats], cex = 0.75, bg = "white")
      }
    }

    ## Split plots
    if (plots.joined == FALSE) {

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

      mins <- sapply(nwstats, function(x) apply(x, 2, min))
      maxs <- sapply(nwstats, function(x) apply(x, 2, max))

      if (missing(targ.col)) {
        targ.col <- rep("black", nstats)
      }

      for (j in outsts) {
        dataj <- data[, colnames(data) %in% nmstats[j], drop = FALSE]
        plot(x = 1, y = 1,
             xlim = xlim,
             ylim = c(min(dataj) * 0.8, max(dataj) * 1.2),
             type = "n", main = nmstats[j],
             xlab = "", ylab = "")

        if (is.numeric(qnts)) {
          if (qnts < 0 | qnts > 1) {
            stop("qnts must be between 0 and 1", call. = FALSE)
          }
          if (missing(qnts.col)) {
            qnts.col <- sim.col
          }
          if (missing(qnts.alpha)) {
            qnts.alpha <- 0.35
          }
          qnts.col <- transco(qnts.col, qnts.alpha)
          quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
          qnt.prev <- apply(dataj, 1, function(x) {
            quantile(x, c(quants[1], quants[2]))
          })
          xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
          if (qnts.smooth == FALSE) {
            yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
          } else {
            yy <- c(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[1, ]))$y,
                    rev(suppressWarnings(supsmu(x = 1:(ncol(qnt.prev)), y = qnt.prev[2, ]))$y))
          }
          polygon(xx, yy, col = qnts.col[which(j == outsts)], border = NA)
        }

        if (sim.lines == TRUE) {
          for (i in sims) {
            lines(dataj[, i], lwd = sim.lwd, col = sim.col[which(j == outsts)])
          }
        }

        if (mean.line == TRUE) {
          if (missing(mean.col)) {
            mean.col <- sim.col
          }
          mean.prev <- rowMeans(dataj)
          if (mean.smooth == TRUE) {
            mean.prev <- suppressWarnings(supsmu(x = 1:length(mean.prev), y = mean.prev))$y
          }
          lines(mean.prev, lwd = mean.lwd,
                col = mean.col[which(j == outsts)], lty = mean.lty)
        }

        if (targ.line == TRUE) {
          if (j %in% targs) {
            abline(h = m$targets[j], lty = targ.lty, lwd = targ.lwd,
                   col = targ.col[which(j == outsts)])
          }
        }
      }

      # Reset graphical parameters
      on.exit(par(ops))
    }
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
#' same as in the \code{summary} functions for a model, but presented graphically
#' as a compartment flow diagram.
#'
#' For \code{dcm} class plots, specify the model run number if the model contains
#' multiple runs, as in a sensitivity analysis. For \code{icm} and \code{netsim}
#' class plots, the \code{run} argument is not used; the plots show the means and
#' standard deviations across simulations at the specified time step.
#'
#' These plots are currently limited to one-group and one-mode models for each of
#' the three model classes. That functionality may be expanded in future software
#' releases.
#'
#' @export
#' @keywords plot
#'
#' @examples
#' ## Example 1: DCM SIR model with varying act.rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = 5:7,
#'                    rec.rate = 1/3, b.rate = 1/90, ds.rate = 1/100,
#'                    di.rate = 1/35, dr.rate = 1/100)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 25, verbose = FALSE)
#' mod1 <- dcm(param, init, control)
#' comp_plot(mod1, at = 25, run = 3)
#'
#' ## Example 2: ICM SIR model with 3 simulations
#' param <- param.icm(inf.prob = 0.2, act.rate = 3, rec.rate = 1/50,
#'                    b.rate = 1/100, ds.rate = 1/100,
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
  if (at > nsteps | at < 1) {
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
      varrow(22, 40, "b.flow", df$b.flow, dir = "in")
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
      varrow(5, 40, "b.flow", df$b.flow, dir = "in")
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
      varrow(22, 40, "b.flow", df$b.flow, dir = "in")
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
  if (class(x) == "icm") {
    groups <- x$param$groups
  }
  if (class(x) == "netsim") {
    groups <- x$param$modes
  }
  if (groups != 1) {
    stop("Only 1-group/mode models currently supported",
         call. = FALSE)
  }

  # Time
  if (at > nsteps | at < 1) {
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
      varrow(22, 40, "b.flow", df.mn$b.flow, dir = "in")
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
      varrow(5, 40, "b.flow", df.mn$b.flow, dir = "in")
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
      varrow(22, 40, "b.flow", df.mn$b.flow, dir = "in")
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



# Helper Functions --------------------------------------------------------


# Calculate denominators
denom <- function(x, y, popfrac) {

  if (class(x) == "dcm") {
    if (popfrac == TRUE) {
      den <- data.frame(den = rep(NA, length(x$control$timesteps)))
      if (x$param$groups == 1) {
        for (i in 1:x$control$nruns) {
          den[, i] <- as.data.frame(x, run = i)$num
        }
        for (i in 1:length(y)) {
          x$epi[[y[i]]] <- x$epi[[y[i]]] / den
        }
      }
      if (x$param$groups == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$control$nruns) {
          den[[1]][, i] <- as.data.frame(x, run = i)$num
          den[[2]][, i] <- as.data.frame(x, run = i)$num.g2
        }
        y.group.num <- ifelse(grepl("num.g2$", y), 2, 1)
        for (j in 1:length(y)) {
          x$epi[[y[j]]] <- x$epi[[y[j]]] / den[[y.group.num[j]]]
        }
      }
    }
    if (popfrac == FALSE && x$control$nruns == 1) {
      for (j in 1:length(y)) {
        x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
      }
    }
  }

  if (class(x) == "icm") {
    if (popfrac == TRUE) {
      den <- data.frame(den = rep(NA, x$control$nsteps))
      if (x$param$groups == 1) {
        for (i in 1:x$control$nsims) {
          den[, i] <- as.data.frame(x, sim = i, out = "vals")$num
        }
        for (i in 1:length(y)) {
          x$epi[[y[i]]] <- x$epi[[y[i]]] / den
        }
      }
      if (x$param$groups == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$control$nsims) {
          den[[1]][, i] <- as.data.frame(x, sim = i, out = "vals")$num
          den[[2]][, i] <- as.data.frame(x, sim = i, out = "vals")$num.g2
        }
        y.group.num <- ifelse(grepl("num.g2$", y), 2, 1)
        for (j in 1:length(y)) {
          x$epi[[y[j]]] <- x$epi[[y[j]]] / den[[y.group.num[j]]]
        }
      }
    }
    if (popfrac == FALSE && x$control$nsims == 1) {
      for (j in 1:length(y)) {
        x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
      }
    }
  }

  if (class(x) == "netsim") {
    if (popfrac == TRUE) {
      den <- data.frame(den = rep(NA, x$control$nsteps))
      if (x$param$modes == 1) {
        for (i in 1:x$control$nsims) {
          den[, i] <- as.data.frame(x, sim = i, out = "vals")$num
        }
        for (i in 1:length(y)) {
          x$epi[[y[i]]] <- x$epi[[y[i]]] / den
        }
      }
      if (x$param$modes == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$control$nsims) {
          den[[1]][, i] <- as.data.frame(x, sim = i, out = "vals")$num
          den[[2]][, i] <- as.data.frame(x, sim = i, out = "vals")$num.m2
        }
        y.group.num <- ifelse(grepl("num.m2$", y), 2, 1)
        for (j in 1:length(y)) {
          x$epi[[y[j]]] <- x$epi[[y[j]]] / den[[y.group.num[j]]]
        }
      }
    }
    if (popfrac == FALSE && x$control$nsims == 1) {
      for (j in 1:length(y)) {
        x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
      }
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
    text(xbox + 10, ybox - 12.5, paste(title, val, sep = "="), cex = 0.8, pos = 4)
  }
  if (dir == "in") {
    arrows(xbox + 10, ybox + 45, xbox + 10, ybox + 20, lwd = 2, length = 0.15)
    text(xbox + 10, ybox + 32.5, paste(title, val, sep = "="), cex = 0.8, pos = 4)
  }
}
