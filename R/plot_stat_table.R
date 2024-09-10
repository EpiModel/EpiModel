plot_stats_table <- function(data, nmstats, method,
                             sim.lines,
                             sim.col = NULL, sim.lwd = NULL, mean.line,
                             mean.smooth, mean.col = NULL,
                             mean.lwd, mean.lty, qnts, qnts.col = NULL,
                             qnts.alpha,
                             qnts.smooth, targ.line, targ.col = NULL, targ.lwd,
                             targ.lty, plots.joined = NULL, draw_legend = NULL,
                             grid, targets,
                             dynamic,
                             xlim = NULL, ylim = NULL,
                             xlab = NULL, ylab = NULL,
                             ...) {
  nstats <- length(nmstats)

  if (is.null(plots.joined)) {
    plots.joined <- nstats <= 3
  }

  if (nstats == 1) {
    plots.joined <- TRUE
    sim.col <- "dodgerblue3"
  }

  xlim <- NVL(xlim, c(1, dim(data)[1]))

  xlab <- if (isFALSE(plots.joined)) "" else NVL(xlab, if (isTRUE(dynamic)) "time" else "simulation number")
  ylab <- if (isFALSE(plots.joined)) "" else NVL(ylab, if (nstats == 1) nmstats else "Statistic")

  if (is.null(sim.lwd)) {
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
    if (!is.null(object)) {
      if (length(object) %in% c(1, nstats)) {
        rep(object, length.out = nstats)
      } else {
        stop(paste0(name, " must be either missing or a vector of length 1 or nstats (", nstats, ")"))
      }
    } else {
      rep(default, length.out = nstats)
    }
  }

  sim.col <- check_len_rep(
    sim.col,
    if (isFALSE(plots.joined)) "dodgerblue3" else seq_len(nstats + 1L)[-1L],
    "sim.col"
  )

  mean.col <- check_len_rep(
    mean.col,
    if (isFALSE(plots.joined)) "black" else sim.col,
    "mean.col"
  )

  qnts.col <- check_len_rep(qnts.col, sim.col, "qnts.col")
  qnts.col <- adjustcolor(qnts.col, qnts.alpha)

  targ.col <- check_len_rep(
    targ.col,
    if (isFALSE(plots.joined) || nstats == 1) "black" else sim.col,
    "targ.col"
  )

  draw_legend <- isTRUE(plots.joined) &&
    ((!is.null(draw_legend) && isTRUE(draw_legend)) ||
       (is.null(draw_legend) && nstats == 1))

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
      if (!is.null(ylim)) {
        ylims[[j]] <- ylim
      } else {
        limdat <- c(
          if (isTRUE(plots.joined) && isTRUE(sim.lines)) data,
          if (isFALSE(plots.joined) && isTRUE(sim.lines)) data[, j, ],
          if (isTRUE(plots.joined) && isTRUE(mean.line)) sapply(means, `[[`, "y"),
          if (isFALSE(plots.joined) && isTRUE(mean.line)) means[[j]]$y,
          if (isTRUE(plots.joined) && isTRUE(draw_qnts)) sapply(qnts_list, `[[`, "y"),
          if (isFALSE(plots.joined) && isTRUE(draw_qnts)) qnts_list[[j]]$y,
          if (isTRUE(plots.joined) && isTRUE(targ.line)) targets,
          if (isFALSE(plots.joined) && isTRUE(targ.line)) targets[j]
        )

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
