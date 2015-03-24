
#' @title Save dcm Data to Output List Format
#'
#' @description This function transfers the data from the master \code{df}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{dcm}}.
#'
#' @param df Master object in \code{\link{dcm}} simulations.
#' @param s Current run number.
#' @param param Param list set in \code{\link{param.dcm}}
#' @param control Control list set in \code{\link{control.dcm}}.
#' @param out Out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.dcm <- function(df, s, param, control, out) {

  if (s == 1) {
    out <- list()
    out$param <- param
    out$control <- control
    out$control$timesteps <- seq(1, control$nsteps, control$dt)
    out$epi <- list()
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]] <- data.frame(df[, j])
    }
  } else {
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]][, s] <- data.frame(df[, j])
    }
  }

  ## Processing for final run
  if (s == control$nruns) {

    # Set column names for varying list elements
    for (s in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[s]]) <- paste0("run", 1:control$nruns)
    }

    # Reset partnership parameters if used
    if (!is.null(param$act.inf.prob)) {
      out$param$inf.prob <- out$param$act.inf.prob
      out$param$act.inf.prob <- NULL
      out$param$act.rate <- NULL
    }

  }

  return(out)
}


#' @title Save icm Data to Output List Format
#'
#' @description This function transfers the data from the master \code{all}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{icm}}.
#'
#' @param dat Master object in \code{icm} simulations.
#' @param s Current run number.
#' @param out Out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.icm <- function(dat, s, out) {

  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$epi <- list()
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }
  } else {
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }
  }

  ## Processing for final run
  if (s == dat$control$nsims) {

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL

    # Set column names for varying list elements
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- paste0("sim", 1:dat$control$nsims)
    }

  }

  return(out)
}


#' @title Save netsim Data to Output List Format
#'
#' @description This function transfers the data from the master \code{dat}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{netsim}}.
#'
#' @param dat Master object in \code{netsim} simulations.
#' @param s Current simulation number.
#' @param out Out list passed back in for updating at simulations 2+.
#'
#' @keywords internal
#' @export
#'
saveout.net <- function(dat, s, out) {

  # Counts number of simulated networks
  num.nw <- ifelse(any(class(dat$nw) == "network"), 1, length(dat$nw))

  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$nwparam <- dat$nwparam
    out$control$num.nw <- num.nw

    out$epi <- list()
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }

    out$stats <- list()
    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats <- list(dat$stats$nwstats)
    }
    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat <- list(dat$stats$transmat)
      } else {
        out$stats$transmat <- list(data.frame())
      }
    }

    if (dat$control$save.network == TRUE) {
      out$network <- list(dat$nw)
    }

    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]] <- list(dat[[el.name]])
      }
    }
  }

  if (s > 1) {
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }

    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats[[s]] <- dat$stats$nwstats
    }
    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat[[s]] <- dat$stats$transmat
      } else {
        out$stats$transmat[[s]] <- data.frame()
      }
    }
    if (dat$control$save.network == TRUE) {
      out$network[[s]] <- dat$nw
    }

    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]][[s]] <- dat[[el.name]]
      }
    }
  }

  ## Final processing
  if (s == dat$control$nsims) {

    # Set names for out
    simnames <- paste0("sim", 1:dat$control$nsims)
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- simnames
    }
    if (dat$control$save.nwstats == TRUE) {
      names(out$stats$nwstats) <- simnames
    }
    if (dat$control$save.transmat == TRUE) {
      names(out$stats$transmat) <- simnames[1:length(out$stats$transmat)]
    }
    if (dat$control$save.network == TRUE) {
      names(out$network) <- simnames
    }

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL
    out$control$currsim <- NULL

    out$temp <- NULL
  }

  return(out)
}


