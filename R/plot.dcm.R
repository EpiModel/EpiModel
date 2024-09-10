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
#' @param grid If \code{TRUE}, a grid is added to the background of plot
#'        (see \code{\link{grid}} for details), with default of nx by ny.
#' @param add If \code{TRUE}, new plot window is not called and lines are added
#'        to existing plot window.
#' @param ... Additional arguments to pass to main plot window (see
#'        \code{\link{plot.default}}).
#'
#' @inheritParams plot.netsim
#' @inheritParams graphics::plot
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
plot.dcm <- function(x, y = NULL, popfrac = FALSE, run = NULL, col = NULL,
                     lwd = NULL, lty = NULL, alpha = 0.9, legend = NULL,
                     leg.name = NULL, leg.cex = 0.8, grid = FALSE, add = FALSE,
                     main = "", xlim = NULL, ylim = NULL, xlab = "Time", ylab = NULL,
                     ...) {

  ## Set missing flags
  noy <- is.null(y)
  norun <- is.null(run)
  nocol <- is.null(col)
  nolwd <- is.null(lwd)
  nolty <- is.null(lty)
  noleg <- is.null(legend)

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
  if (is.null(ylim)) {
    ylim <- c(0, ymax)
  }

  if (is.null(xlim)) {
    xlim <- c(0, nsteps)
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

  ## Defaults for ylab
  if (is.null(ylab)) {
    if (popfrac == FALSE) {
      ylab <- "Number"
    } else {
      ylab <- "Prevalence"
    }
  }

  ## Main plot window
  if (add == FALSE) {
    plot(1, 1, type = "n", bty = "n",
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ...)
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
  if (is.null(leg.name)) {
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
