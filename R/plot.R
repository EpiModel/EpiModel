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



## Helper utilities ------------------------------------------------------------
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
    mean.prev <- rowMeans(x[[loc]][[y[j]]], na.rm = TRUE)
    xs <- seq_len(length(mean.prev))
    if (mean.smooth == TRUE) {
      smoother <- suppressWarnings(supsmu(x = xs, y = mean.prev))
      mean.prev <- smoother$y
      xs <- smoother$x
    }
    if (plot.means == 1) {
      lines(x = xs, y = mean.prev, lwd = mean.lwd[j],
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
