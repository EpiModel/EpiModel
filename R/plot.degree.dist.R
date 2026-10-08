## Plotting helper for the degree distributions returned by get_degree_dist.
##
## `dists` is a named list of data.frames, one per distribution to draw, in the
## long format returned by `get_degree_dist`. All are drawn in a single panel
## on a shared degree axis, so that the momentary and cumulative distributions
## of one model can be compared directly.
plot_degree_dist <- function(dists, method = "l", maxdeg = NULL,
                             sim.lines = FALSE, sim.col = NULL, sim.lwd = NULL,
                             mean.line = TRUE, mean.col = NULL, mean.lwd = 2,
                             mean.lty = 1, qnts = 0.5, qnts.col = NULL,
                             qnts.alpha = 0.5, draw_legend = NULL, grid = FALSE,
                             xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
                             ...) {

  method <- match.arg(method, c("l", "b"))
  nseries <- length(dists)
  nsims <- length(unique(dists[[1]]$sim))

  ## shared degree axis across the distributions being drawn
  ndeg <- max(vapply(dists, function(d) max(d$degree), numeric(1))) + 1

  ## proportion of nodes at each degree, as a (degree x sim) matrix
  mats <- lapply(dists, function(d) {
    m <- matrix(0, nrow = ndeg, ncol = nsims)
    obs <- seq_len(length(d$degree) / nsims)
    m[obs, ] <- matrix(d$prop, ncol = nsims)
    m
  })

  ## mean degree of each distribution, before any truncation of the axis
  means <- vapply(mats, function(m) mean(colSums(m * (seq_len(ndeg) - 1))),
                  numeric(1))

  ## collapse the tail of the axis into a final "maxdeg+" category
  lumped <- !is.null(maxdeg) && maxdeg < ndeg - 1
  if (lumped) {
    mats <- lapply(mats, function(m) {
      out <- m[seq_len(maxdeg + 1), , drop = FALSE]
      out[maxdeg + 1, ] <- colSums(m[seq(maxdeg + 1, ndeg), , drop = FALSE])
      out
    })
    ndeg <- maxdeg + 1
  }
  degs <- seq_len(ndeg) - 1

  ## Color palettes, following plot.netsim: blue then red
  bpal <- c(4, 2, 3, 5:100)
  sim.col <- rep(NVL(sim.col, bpal), length.out = nseries)
  mean.col <- rep(NVL(mean.col, bpal), length.out = nseries)
  qnts.col <- adjustcolor(rep(NVL(qnts.col, bpal), length.out = nseries),
                          qnts.alpha)
  sim.lwd <- NVL(sim.lwd, max(c(1 - nsims * 0.05, 0.5)))

  draw_qnts <- is.numeric(qnts) && nsims > 1
  if (draw_qnts) {
    if (qnts < 0 || qnts > 1) {
      stop("qnts must be between 0 and 1")
    }
    qnt.probs <- c((1 - qnts) / 2, 1 - (1 - qnts) / 2)
    qnt.bands <- lapply(mats, function(m) {
      apply(m, 1, quantile, probs = qnt.probs, na.rm = TRUE)
    })
  }

  xlab <- NVL(xlab, "Degree")
  ylab <- NVL(ylab, "Proportion of Nodes")
  ylim <- NVL(ylim, c(0, 1.05 * max(unlist(c(
    lapply(mats, rowMeans),
    if (isTRUE(sim.lines)) mats,
    if (isTRUE(draw_qnts)) qnt.bands
  )))))

  draw_legend <- NVL(draw_legend, nseries > 1)
  legend.txt <- paste0(names(dists), " (mean = ", round(means, 2), ")")

  ## axis labels: every degree if the axis is short, otherwise pretty ticks,
  ## with the collapsed final category flagged with a "+"
  ticks <- if (ndeg <= 15) degs else degs[degs %in% pretty(degs)]
  if (lumped && !maxdeg %in% ticks) {
    ticks <- c(ticks, maxdeg)
  }
  tick.labs <- as.character(ticks)
  if (lumped) {
    tick.labs[ticks == maxdeg] <- paste0(maxdeg, "+")
  }

  if (method == "l") {
    plot(NULL, xlim = NVL(xlim, range(degs)), ylim = ylim, type = "n",
         xlab = xlab, ylab = ylab, xaxt = "n", ...)
    axis(1, at = ticks, labels = tick.labs)

    for (j in seq_len(nseries)) {
      if (isTRUE(draw_qnts)) {
        polygon(c(degs, rev(degs)),
                c(qnt.bands[[j]][1, ], rev(qnt.bands[[j]][2, ])),
                col = qnts.col[j], border = NA)
      }
      if (isTRUE(sim.lines)) {
        apply(mats[[j]], 2, function(y) {
          lines(degs, y, lwd = sim.lwd, col = sim.col[j])
        })
      }
      if (isTRUE(mean.line)) {
        lines(degs, rowMeans(mats[[j]]), lwd = mean.lwd, col = mean.col[j],
              lty = mean.lty)
        points(degs, rowMeans(mats[[j]]), pch = 16, col = mean.col[j])
      }
    }

    if (isTRUE(grid)) {
      grid()
    }
    if (isTRUE(draw_legend)) {
      legend("topright", legend = legend.txt, lwd = 2, col = mean.col,
             cex = 0.75, bg = "white")
    }
  }

  if (method == "b") {
    heights <- do.call(rbind, lapply(mats, rowMeans))
    dimnames(heights) <- list(names(dists), degs)

    ## bar labels are drawn from the tick set rather than by barplot itself,
    ## which thins them unpredictably when the degree axis is long
    bar.args <- list(height = heights, beside = TRUE, ylim = ylim,
                     xlab = xlab, ylab = ylab, axisnames = FALSE, ...)
    bar.args$col <- NVL(bar.args$col, mean.col)
    mids <- do.call(barplot, bar.args)
    axis(1, at = colMeans(mids)[match(ticks, degs)], labels = tick.labs,
         lty = 0)

    if (isTRUE(draw_qnts)) {
      for (j in seq_len(nseries)) {
        segments(mids[j, ], qnt.bands[[j]][1, ], mids[j, ],
                 qnt.bands[[j]][2, ])
      }
    }

    if (isTRUE(grid)) {
      grid(nx = NA, ny = NULL)
    }
    if (isTRUE(draw_legend)) {
      legend("topright", legend = legend.txt, pch = 15, col = bar.args$col,
             cex = 0.75, bg = "white")
    }
  }
}
