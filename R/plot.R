## Helper utilities ------------------------------------------------------------
draw_qnts <- function(x, y, qnts, qnts.pal, qnts.smooth,
                      loc = "epi", plot.qnts = 1, qnts.min_max = "max") {

  qnt.min <- 1E10
  qnt.max <- -1E10

  lcomp <- length(y)
  for (j in seq_len(lcomp)) {
    quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
    qnt.prev <- apply(
      x[[loc]][[y[j]]], 1,
      function(x) quantile(x, c(quants[1], quants[2]), na.rm = TRUE)
    )
    complete_rows <- complete.cases(t(qnt.prev))
    x_cords <- seq_len(ncol(qnt.prev))[complete_rows]
    xx <- c(x_cords, rev(x_cords))
    qnt.prev <- qnt.prev[, complete_rows]
    if (qnts.smooth == FALSE) {
      yy <- c(qnt.prev[1, ], rev(qnt.prev[2, ]))
    } else {
      yy <- c(
        suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
                                y = qnt.prev[1, ]))$y,
        rev(suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
                                    y = qnt.prev[2, ]))$y)
      )
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
    ss2 <- suppressWarnings(supsmu(x = seq_len(ncol(qnt.prev)),
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
    ss <- suppressWarnings(
      supsmu(x = seq_along(mean.prev), y = mean.prev)
    )
    xx <- ss$x
    yy <- ss$y
  }
  list(x = xx, y = yy)
}

# Helper Functions --------------------------------------------------------


# Calculate denominators
denom <- function(x, y, popfrac) {

  cont.val <- if (inherits(x, "dcm")) "nruns" else "nsims"
  if (popfrac) {
    for (i in seq_along(y)) {
      dname <- paste(strsplit(y[i], "[.]")[[1]][-1], collapse = ".")
      x$epi[[y[i]]] <- x$epi[[y[i]]] / x$epi[[dname]]
    }
  }
  if (!popfrac && x$control[[cont.val]] == 1) {
    for (j in seq_along(y)) {
      x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
    }
  }

  return(x)
}
