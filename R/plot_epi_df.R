# Helper function to convert a `data.frame` obtained with
# `as.data.frame(netsim)` back into an `epi` list format.
epi2df <- function(epi) {
  n_sims <- ncol(epi[[1]])
  n_steps <- nrow(epi[[1]])
  out <- vector(mode = "list", length = length(epi) + 2)
  names(out) <- c("sim", "time", names(epi))
  out$sim <- rep(seq_len(n_sims), each = n_steps)
  out$time <- rep(seq_len(n_steps), n_sims)
  for (nme in names(epi)) {
    out[[nme]] <- as.vector(as.matrix(epi[[nme]]))
  }
  return(as.data.frame(out))
}

#' @title Plot Epidemic Model Results From a Netsim Data.Frame
#'
#' @description This function is a wrapper around `plot.netsim` accepting a
#'   `data.frame` obtain with `as.data.frame(netsim_object)`.
#'
#' @param df A `data.frame` obtain with `as.data.frame(netsim_object)`.
#' @inheritParams plot.netsim
#'
#' @export
plot_epi_df <- function(df, y = NULL, sims = NULL, legend = NULL,
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
