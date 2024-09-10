# Helper function to convert a `data.frame` obtained with
# `as.data.frame(netsim)` back into an `epi` list format.
df2epi <- function(x) {
  n_steps <- max(x$time)
  sims <- unique(x$sim)
  n_sims <- length(sims)
  colnames <- paste0("sim", sims)
  epi_names <- setdiff(names(x), c("time", "sim"))
  sapply(
    epi_names,
    \(nme) {
      t <- as.data.frame(matrix(x[[nme]], nrow = n_steps, ncol = n_sims))
      names(t) <- colnames
      t
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
}

#' @title Plot Epidemic Model Results From a Netsim Data.Frame
#'
#' @description This function is a wrapper around `plot.netsim` accepting a
#'   `data.frame` obtain with `as.data.frame(netsim_object)`.
#'
#' @param df A `data.frame` obtain with `as.data.frame(netsim_object)`.
#' @inheritParams plot.netsim
#'
#' @method plot epi.data.frame
#' @export
plot.epi.data.frame <- function(df, y = NULL, sims = NULL, legend = NULL,
                                mean.col = NULL, qnts.col = NULL, sim.lwd = NULL,
                                sim.col = NULL, sim.alpha = NULL, popfrac = FALSE,
                                qnts = 0.5, qnts.alpha = 0.5, qnts.smooth = TRUE,
                                mean.line = TRUE, mean.smooth = TRUE, add = FALSE,
                                mean.lwd = 2, mean.lty = 1, xlim = NULL,
                                ylim = NULL, main = NULL, xlab = NULL, ylab = NULL,
                                sim.lines = FALSE, grid = FALSE, leg.cex = 0.8,
                                ...) {
  ntemp <- list(
    epi = df2epi(df),
    control = list(
      nsteps = max(df$time),
      nsims = max(df$sim)
    )
  )
  EpiModel:::plot_netsim_epi(
    ntemp, y, sims, legend, mean.col, qnts.col, sim.lwd, sim.col, sim.alpha,
    popfrac, qnts, qnts.alpha, qnts.smooth, mean.line, mean.smooth, add,
    mean.lwd, mean.lty, xlim, ylim, main, xlab, ylab, sim.lines, grid,
    leg.cex, ...
  )
}
